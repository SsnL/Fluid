#include "pathtracer.h"
#include "bsdf.h"
#include "ray.h"

#include <stack>
#include <random>
#include <algorithm>
#include <sstream>

#include "CGL/CGL.h"
#include "CGL/vector3D.h"
#include "CGL/matrix3x3.h"
#include "CGL/lodepng.h"

#include "GL/glew.h"

#include "static_scene/sphere.h"
#include "static_scene/triangle.h"
#include "static_scene/light.h"

using namespace CGL::StaticScene;

using std::min;
using std::max;

namespace CGL {

PathTracer::PathTracer(size_t ns_aa,
                       size_t max_ray_depth, size_t ns_area_light,
                       size_t ns_diff, size_t ns_glsy, size_t ns_refr,
                       size_t num_threads, HDRImageBuffer* envmap) {
  state = INIT,
  this->ns_aa = ns_aa;
  this->max_ray_depth = max_ray_depth;
  this->ns_area_light = ns_area_light;
  this->ns_diff = ns_diff;
  this->ns_glsy = ns_diff;
  this->ns_refr = ns_refr;
  // this->fluid_particles = particles;


  if (envmap) {
    this->envLight = new EnvironmentLight(envmap);
  } else {
    this->envLight = NULL;
  }

  bvh = NULL;
  scene = NULL;
  camera = NULL;

  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  show_rays = true;

  imageTileSize = 32;
  numWorkerThreads = num_threads;
  workerThreads.resize(numWorkerThreads);

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;

}

PathTracer::~PathTracer() {

  delete bvh;
  delete gridSampler;
  delete hemisphereSampler;

}

void PathTracer::set_scene(Scene *scene) {

  if (state != INIT) {
    return;
  }

  if (this->scene != nullptr) {
    delete scene;
    delete bvh;
    selectionHistory.pop();
  }

  if (this->envLight != nullptr) {
    scene->lights.push_back(this->envLight);
  }

  this->scene = scene;
  build_accel();

  if (has_valid_configuration()) {
    state = READY;
  }
}

void PathTracer::set_camera(Camera *camera) {

  if (state != INIT) {
    return;
  }

  this->camera = camera;
  if (has_valid_configuration()) {
    state = READY;
  }

}

void PathTracer::set_frame_size(size_t width, size_t height) {
  if (state != INIT && state != READY) {
    stop();
  }
  sampleBuffer.resize(width, height);
  frameBuffer.resize(width, height);
  if (has_valid_configuration()) {
    state = READY;
  }
}

bool PathTracer::has_valid_configuration() {
  return scene && camera && gridSampler && hemisphereSampler &&
         (!sampleBuffer.is_empty());
}

void PathTracer::update_screen() {
  switch (state) {
    case INIT:
    case READY:
      break;
    case VISUALIZE:
      visualize_accel();
      break;
    case RENDERING:
      glDrawPixels(frameBuffer.w, frameBuffer.h, GL_RGBA,
                   GL_UNSIGNED_BYTE, &frameBuffer.data[0]);
      break;
    case DONE:
        //sampleBuffer.tonemap(frameBuffer, tm_gamma, tm_level, tm_key, tm_wht);
      glDrawPixels(frameBuffer.w, frameBuffer.h, GL_RGBA,
                   GL_UNSIGNED_BYTE, &frameBuffer.data[0]);
      break;
  }
}

void PathTracer::stop() {
  switch (state) {
    case INIT:
    case READY:
      break;
    case VISUALIZE:
      while (selectionHistory.size() > 1) {
        selectionHistory.pop();
      }
      state = READY;
      break;
    case RENDERING:
      continueRaytracing = false;
    case DONE:
      for (int i=0; i<numWorkerThreads; i++) {
            workerThreads[i]->join();
            delete workerThreads[i];
        }
      state = READY;
      break;
  }
}

void PathTracer::clear() {
  if (state != READY) return;
  delete bvh;
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  selectionHistory.pop();
  sampleBuffer.resize(0, 0);
  frameBuffer.resize(0, 0);
  state = INIT;
}

void PathTracer::start_visualizing() {
  if (state != READY) {
    return;
  }
  state = VISUALIZE;
  // fluid sim params //
  fprintf(stdout, "[Fluid Simulation] %s", fluid_particles->paramsString().c_str()); fflush(stdout);
}

void PathTracer::start_raytracing() {
  if (state != READY) return;

  build_accel(true);

  rayLog.clear();
  workQueue.clear();

  state = RENDERING;
  continueRaytracing = true;
  workerDoneCount = 0;

  sampleBuffer.clear();
  frameBuffer.clear();
  num_tiles_w = sampleBuffer.w / imageTileSize + 1;
  num_tiles_h = sampleBuffer.h / imageTileSize + 1;
  tilesTotal = num_tiles_w * num_tiles_h;
  tilesDone = 0;
  tile_samples.resize(num_tiles_w * num_tiles_h);
  memset(&tile_samples[0], 0, num_tiles_w * num_tiles_h * sizeof(int));

  // populate the tile work queue
  for (size_t y = 0; y < sampleBuffer.h; y += imageTileSize) {
      for (size_t x = 0; x < sampleBuffer.w; x += imageTileSize) {
          workQueue.put_work(WorkItem(x, y, imageTileSize, imageTileSize));
      }
  }

  bvh->total_isects = 0; bvh->total_rays = 0;
  // launch threads
  fprintf(stdout, "[PathTracer] Rendering... "); fflush(stdout);
  for (int i=0; i<numWorkerThreads; i++) {
      workerThreads[i] = new std::thread(&PathTracer::worker_thread, this);
  }
}

void PathTracer::render_to_file(string filename) {
  cout << "RENDER" << endl;
  unique_lock<std::mutex> lk(m_done);
  start_raytracing();
  cv_done.wait(lk, [this]{ return state == DONE; });
  lk.unlock();
  save_image(filename);
  fprintf(stdout, "[PathTracer] Job completed.\n");
}


void PathTracer::build_accel(bool includeSurface) {

  // collect primitives //
  fprintf(stdout, "[PathTracer] Collecting primitives... "); fflush(stdout);
  timer.start();
  vector<Primitive *> primitives;
  for (SceneObject *obj : scene->objects) {
    const vector<Primitive *> &obj_prims = obj->get_primitives();
    primitives.reserve(primitives.size() + obj_prims.size());
    primitives.insert(primitives.end(), obj_prims.begin(), obj_prims.end());
  }
  if (includeSurface) {
    fluid_particles->updateSurface();
    fprintf(stdout, "including %lu marching cube surfacing triangles... ", fluid_particles->surface.size()); fflush(stdout);
    primitives.reserve(primitives.size() + fluid_particles->surface.size());
    primitives.insert(primitives.end(), fluid_particles->surface.begin(), fluid_particles->surface.end());
    bvh_has_surface = true;
  } else {
    bvh_has_surface = false;
  }
  cout << bvh_has_surface << ' ';
  timer.stop();
  fprintf(stdout, "Done! (%.4f sec)\n", timer.duration());

  // build BVH //
  fprintf(stdout, "[PathTracer] Building BVH from %lu primitives... ", primitives.size());
  fflush(stdout);
  timer.start();
  bvh = new BVHAccel(primitives);
  fluid_particles->bvh = bvh;
  timer.stop();
  fprintf(stdout, "Done! (%.4f sec)\n", timer.duration());

  // initial visualization //
  selectionHistory.push(bvh->get_root());
}

void PathTracer::visualize_accel() const {

  glPushAttrib(GL_ENABLE_BIT);
  glDisable(GL_LIGHTING);
  glLineWidth(1);
  glEnable(GL_DEPTH_TEST);

  // hardcoded color settings
  Color cnode = Color(.5, .5, .5, .25);
  Color cnode_hl = Color(1., .25, .0, .6);
  Color cnode_hl_child = Color(1., 1., 1., .6);

  Color cprim_hl_left = Color(.6, .6, 1., 1);
  Color cprim_hl_right = Color(.8, .8, 1., 1);
  Color cprim_hl_edges = Color(0., 0., 0., 0.5);

  BVHNode *selected = selectionHistory.top();

  // render solid geometry (with depth offset)
  glPolygonOffset(1.0, 1.0);
  glEnable(GL_POLYGON_OFFSET_FILL);

  fluid_particles->redraw(Color(1, 1, 1, 1));
  if (selected->isLeaf()) {
    bvh->draw(selected, cprim_hl_left);
  } else {
    bvh->draw(selected->l, cprim_hl_left);
    bvh->draw(selected->r, cprim_hl_right);
  }

  glDisable(GL_POLYGON_OFFSET_FILL);

  // draw geometry outline
  bvh->drawOutline(selected, cprim_hl_edges);

  // keep depth buffer check enabled so that mesh occluded bboxes, but
  // disable depth write so that bboxes don't occlude each other.
  glDepthMask(GL_FALSE);

  // create traversal stack
  stack<BVHNode *> tstack;

  // push initial traversal data
  tstack.push(bvh->get_root());

  // draw all BVH bboxes with non-highlighted color
  // while (!tstack.empty()) {

  //   BVHNode *current = tstack.top();
  //   tstack.pop();

  //   current->bb.draw(cnode);
  //   if (current->l) tstack.push(current->l);
  //   if (current->r) tstack.push(current->r);
  // }

  // draw selected node bbox and primitives
  // if (selected->l) selected->l->bb.draw(cnode_hl_child);
  // if (selected->r) selected->r->bb.draw(cnode_hl_child);

  glLineWidth(3.f);
  // selected->bb.draw(cnode_hl);

  // now perform visualization of the rays
  if (show_rays) {
      glLineWidth(1.f);
      glBegin(GL_LINES);

      for (size_t i=0; i<rayLog.size(); i+=500) {

          const static double VERY_LONG = 10e4;
          double ray_t = VERY_LONG;

          // color rays that are hits yellow
          // and rays this miss all geometry red
          if (rayLog[i].hit_t >= 0.0) {
              ray_t = rayLog[i].hit_t;
              glColor4f(1.f, 1.f, 0.f, 0.1f);
          } else {
              glColor4f(1.f, 0.f, 0.f, 0.1f);
          }

          Vector3D end = rayLog[i].o + ray_t * rayLog[i].d;

          glVertex3f(rayLog[i].o[0], rayLog[i].o[1], rayLog[i].o[2]);
          glVertex3f(end[0], end[1], end[2]);
      }
      glEnd();
  }

  glDepthMask(GL_TRUE);
  glPopAttrib();
}

void PathTracer::key_press(int key) {

  BVHNode *current = selectionHistory.top();
  switch (key) {
  case ']':
      ns_aa *=2;
      fprintf(stdout, "[PathTracer] Samples per pixel changed to %lu\n", ns_aa);
      //tm_key = clamp(tm_key + 0.02f, 0.0f, 1.0f);
      break;
  case '[':
      //tm_key = clamp(tm_key - 0.02f, 0.0f, 1.0f);
      ns_aa /=2;
      if (ns_aa < 1) ns_aa = 1;
      fprintf(stdout, "[PathTracer] Samples per pixel changed to %lu\n", ns_aa);
      break;
  case '=': case '+':
      ns_area_light *= 2;
      fprintf(stdout, "[PathTracer] Area light sample count increased to %zu.\n", ns_area_light);
      break;
  case '-': case '_':
      if (ns_area_light > 1) ns_area_light /= 2;
      fprintf(stdout, "[PathTracer] Area light sample count decreased to %zu.\n", ns_area_light);
      break;
  case '.': case '>':
      max_ray_depth++;
      fprintf(stdout, "[PathTracer] Max ray depth increased to %zu.\n", max_ray_depth);
      break;
  case ',': case '<':
      if (max_ray_depth) max_ray_depth--;
      fprintf(stdout, "[PathTracer] Max ray depth decreased to %zu.\n", max_ray_depth);
      break;
  case KEYBOARD_UP:
      if (current != bvh->get_root()) {
          selectionHistory.pop();
      }
      break;
  case KEYBOARD_LEFT:
      if (current->l) {
          selectionHistory.push(current->l);
      }
      break;
  case KEYBOARD_RIGHT:
      if (current->l) {
          selectionHistory.push(current->r);
      }
      break;
  case 'm': case 'M':
      fluid_simulate_time_step();
      visualize_accel();
      fprintf(stdout, "[Fluid Simulation] Fluid particles updated.\n");
      break;
  case 'g': case 'G':
      fluid_simulate_time(1, true);
      visualize_accel();
      fprintf(stdout, "[Fluid Simulation] Fluid particles updated, screenshots saved.\n");
      break;
  case 'c': case 'C':
      fluid_simulate_time(1);
      visualize_accel();
      fprintf(stdout, "[Fluid Simulation] Fluid particles updated.\n");
      break;
  case 'h': case 'H': // wonder happens here
      fluid_simulate_render_time(1);
      visualize_accel();
      fprintf(stdout, "[Fluid Simulation] Fluid particles updated, screenshots saved.\n");
      break;
  case 's': case 'S':
      save_glimage();
      break;
  case 'a': case 'A':
      show_rays = !show_rays;
  default:
      return;
  }
}

void PathTracer::fluid_simulate_time(double delta_t, bool save_png) {
  if (bvh_has_surface) {
    build_accel();
  }
  double start_t = fluid_particles->simulate_time;
  if (save_png) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    visualize_accel();
    save_glimage();
  }
  while (fluid_particles->simulate_time < start_t + delta_t) {
    fluid_particles->timeStep();
    if (save_png) {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      visualize_accel();
      save_glimage();
    }
  }
}

void PathTracer::fluid_simulate_time_step() {
  if (bvh_has_surface) {
    build_accel();
  }
  fluid_particles->timeStep();
}

void PathTracer::fluid_simulate_render_time(double delta_t) {
  double start_t = fluid_particles->simulate_time;
  stop();
  render_to_file(timestamp_based_png_file_name());
  while (fluid_particles->simulate_time < start_t + delta_t) {
    fluid_simulate_time_step();
    stop();
    render_to_file(timestamp_based_png_file_name());
  }
  start_visualizing();
}


Spectrum PathTracer::estimate_direct_lighting(const Ray& r, const Intersection& isect) {
  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  Vector3D n = isect.n;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d) /* object coordinates */;

  Spectrum L_out, L_light_sample;
  size_t num_samples, i;
  Vector3D w_i /* world coordinates */, w_in /* object coordinates */;
  float dist_to_light, pdf;
  BSDF *bsdf = isect.bsdf;

  for (SceneLight *light : scene->lights) {
    if (light->is_delta_light()) {
      num_samples = 1;
    } else {
      num_samples = ns_area_light;
    }
    Spectrum L_light;
    for (i = 0; i < num_samples; i++) {
      L_light_sample = light->sample_L(hit_p, &w_i, &dist_to_light, &pdf);
      w_in = w2o * w_i;
      if (w_in.z < 0.0) {
        continue;
      }
      Ray shadow_ray(hit_p, w_i, EPS_D, dist_to_light);
      if (bvh->intersect(shadow_ray)) {
        continue;
      }
      L_light += \
        L_light_sample * bsdf->f(w_out, w_in) * abs_cos_theta(w_in) / pdf;
    }
    L_out += L_light / ((double) num_samples);
  }

  return L_out;
}

Spectrum PathTracer::estimate_indirect_lighting(
  const Ray& r, const Intersection& isect
) {

  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d), w_in;
  float pdf;
  BSDF *bsdf = isect.bsdf;

  Spectrum bsdf_spectrum = bsdf->sample_f(w_out, &w_in, &pdf);

  // Russian roulette continue probability
  double p_cont = fmin(bsdf_spectrum.illum() * 15.0, 1.0);
  if (coin_flip(p_cont)) {
    return bsdf_spectrum * trace_ray(
      Ray(hit_p, o2w * w_in, EPS_D, INF_D, r.depth - 1),
      bsdf->is_delta()
    ) * abs_cos_theta(w_in) / p_cont / pdf;
  }

  return Spectrum();

}

Spectrum PathTracer::trace_ray(const Ray &r, bool includeLe) {

  Intersection isect;
  Spectrum L_out;

  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.
  if (!bvh->intersect(r, &isect)) {
    return L_out;
  }

  // return normal_shading(isect.n);

  // We only include the emitted light if the previous BSDF was a delta
  // distribution or if the previous ray came from the camera.
  if (includeLe)
    L_out += isect.bsdf->get_emission();

  // Delta BSDFs have no direct lighting since they are zero with probability 1
  // -- their values get accumulated through indirect lighting, where the BSDF
  // gets to sample itself.
  if (!isect.bsdf->is_delta())
    L_out += estimate_direct_lighting(r, isect);

  // If the ray's depth is zero, then the path must terminate
  // and no further indirect lighting is calculated.
  if (r.depth > 0)
    L_out += estimate_indirect_lighting(r, isect);

  return L_out;

}

Spectrum PathTracer::raytrace_pixel(size_t x, size_t y) {

  // Make a loop that generates num_samples camera rays and traces them
  // through the scene. Return the average Spectrum.

  int num_samples = ns_aa; // total samples to evaluate
  double w = sampleBuffer.w, h = sampleBuffer.h;
  if (num_samples == 1) {
    Ray r = camera->generate_ray((x + 0.5) / w, (y + 0.5) / h);
    r.depth = max_ray_depth;
    return trace_ray(r, true);
  }
  Spectrum s = Spectrum();
  Vector2D d;
  for (size_t i = 0; i < num_samples; i++) {
    d = gridSampler->get_sample();
    Ray r = camera->generate_ray((x + d.x) / w, (y + d.y) / h);
    r.depth = max_ray_depth;
    s += trace_ray(r, true);
  }
  return s / ((double) ns_aa);
}

void PathTracer::raytrace_tile(int tile_x, int tile_y,
                               int tile_w, int tile_h) {

  size_t w = sampleBuffer.w;
  size_t h = sampleBuffer.h;

  size_t tile_start_x = tile_x;
  size_t tile_start_y = tile_y;

  size_t tile_end_x = std::min(tile_start_x + tile_w, w);
  size_t tile_end_y = std::min(tile_start_y + tile_h, h);

  size_t tile_idx_x = tile_x / imageTileSize;
  size_t tile_idx_y = tile_y / imageTileSize;
  size_t num_samples_tile = tile_samples[tile_idx_x + tile_idx_y * num_tiles_w];

  for (size_t y = tile_start_y; y < tile_end_y; y++) {
    if (!continueRaytracing) return;
    for (size_t x = tile_start_x; x < tile_end_x; x++) {
        Spectrum s = raytrace_pixel(x, y);
        sampleBuffer.update_pixel(s, x, y);
    }
  }

  tile_samples[tile_idx_x + tile_idx_y * num_tiles_w] += 1;
  sampleBuffer.toColor(frameBuffer, tile_start_x, tile_start_y, tile_end_x, tile_end_y);
}

void PathTracer::worker_thread() {

  Timer timer;
  timer.start();

  WorkItem work;
  while (continueRaytracing && workQueue.try_get_work(&work)) {
    raytrace_tile(work.tile_x, work.tile_y, work.tile_w, work.tile_h);
    {
      lock_guard<std::mutex> lk(m_done);
      ++tilesDone;
      cout << "\r[PathTracer] Rendering... " << int((double)tilesDone/tilesTotal * 100) << '%';
      cout.flush();
    }
  }

  workerDoneCount++;
  if (!continueRaytracing && workerDoneCount == numWorkerThreads) {
    timer.stop();
    fprintf(stdout, "\n[PathTracer] Rendering canceled!\n");
    state = READY;
  }

  if (continueRaytracing && workerDoneCount == numWorkerThreads) {
    timer.stop();
    fprintf(stdout, "\r[PathTracer] Rendering... 100%%! (%.4fs)\n", timer.duration());
    fprintf(stdout, "[PathTracer] BVH traced %llu rays.\n", bvh->total_rays);
    fprintf(stdout, "[PathTracer] Averaged %f intersection tests per ray.\n", (((double)bvh->total_isects)/bvh->total_rays));

    lock_guard<std::mutex> lk(m_done);
    state = DONE;
    cv_done.notify_one();
  }
}

void PathTracer::save_image(string filename) {

  if (state != DONE) return;

  if (filename == "") {
    filename = timestamp_based_png_file_name();
  }


  uint32_t* frame = &frameBuffer.data[0];
  size_t w = frameBuffer.w;
  size_t h = frameBuffer.h;
  uint32_t* frame_out = new uint32_t[w * h];
  for(size_t i = 0; i < h; ++i) {
    memcpy(frame_out + i * w, frame + (h - i - 1) * w, 4 * w);
  }

  fprintf(stderr, "[PathTracer] Saving to file: %s... ", filename.c_str());
  lodepng::encode(filename, (unsigned char*) frame_out, w, h);
  fprintf(stderr, "Done!\n");
}

void PathTracer::save_glimage(string filename) {

  if (state != READY && state != VISUALIZE) return;

  if (filename == "") {
    filename = timestamp_based_png_file_name();
  }

  int width = sampleBuffer.w, height = sampleBuffer.h;
  vector<unsigned char> windowPixels( 4*width*height );
  glReadPixels(0, 0,
              width,
              height,
              GL_RGBA,
              GL_UNSIGNED_BYTE,
              &windowPixels[0] );

  vector<unsigned char> flippedPixels( 4*width*height );
  for (int row = 0; row < height; ++row)
    memcpy(&flippedPixels[row * width * 4], &windowPixels[(height - row - 1) * width * 4], 4*width);

  stringstream ss;
  fprintf(stderr, "[PathTracer] Saving to file: %s... ", filename.c_str());
  if (lodepng::encode(filename, flippedPixels, width, height))
    cerr << "Could not be written" << endl;
  else
    cout << "Success!" << endl;
}

string PathTracer::timestamp_based_png_file_name() {
  time_t rawtime;
  time (&rawtime);

  time_t t = time(nullptr);
  tm *lt = localtime(&t);
  stringstream ss;
  ss << "screenshot_" << lt->tm_mon+1 << "-" << lt->tm_mday << "_"
    << lt->tm_hour << "-" << lt->tm_min << "-" << lt->tm_sec << ".png";
  return ss.str();
}

}  // namespace CGL
