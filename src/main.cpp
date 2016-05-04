#include "CGL/CGL.h"
#include "CGL/viewer.h"

#define TINYEXR_IMPLEMENTATION
#include "CGL/tinyexr.h"

#include "application.h"
#include "image.h"

#include <iostream>
typedef uint32_t gid_t;
#include <unistd.h>

using namespace std;
using namespace CGL;

#define msg(s) cerr << "[PathTracer] " << s << endl;

void usage(const char* binaryName) {
  printf("Usage: %s [options] <scenefile>\n", binaryName);
  printf("Program Options:\n");
  printf("  -p  <PATH>       Path to particle definition\n");
  printf("  -d  <INT>        Simulation time of fluid rendered and saved as png files.\n");
  printf("  -s  <INT>        Number of camera rays per pixel\n");
  printf("  -l  <INT>        Number of samples per area light\n");
  printf("  -t  <INT>        Number of render threads\n");
  printf("  -m  <INT>        Maximum ray depth\n");
  printf("  -e  <PATH>       Path to environment map\n");
  printf("  -f  <FILENAME>   Image (.png) file to save output to in windowless mode\n");
  printf("  -r  <INT> <INT>  Width and height of output image (if windowless)\n");
  printf("  -h               Print this help message\n");
  printf("\n");
}

HDRImageBuffer* load_exr(const char* file_path) {

  const char* err;

  EXRImage exr;
  InitEXRImage(&exr);

  int ret = ParseMultiChannelEXRHeaderFromFile(&exr, file_path, &err);
  if (ret != 0) {
    msg("Error parsing OpenEXR file: " << err);
    return NULL;
  }

  for (int i = 0; i < exr.num_channels; i++) {
    if (exr.pixel_types[i] == TINYEXR_PIXELTYPE_HALF) {
      exr.requested_pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT;
    }
  }

  ret = LoadMultiChannelEXRFromFile(&exr, file_path, &err);
  if (ret != 0) {
    msg("Error loading OpenEXR file: " << err);
    exit(EXIT_FAILURE);
  }

  HDRImageBuffer* envmap = new HDRImageBuffer();
  envmap->resize(exr.width, exr.height);
  float* channel_r = (float*) exr.images[2];
  float* channel_g = (float*) exr.images[1];
  float* channel_b = (float*) exr.images[0];
  for (size_t i = 0; i < exr.width * exr.height; i++) {
    envmap->data[i] = Spectrum(channel_r[i],
                               channel_g[i],
                               channel_b[i]);
  }

  return envmap;
}

int main( int argc, char** argv ) {

  // get the options
  AppConfig config; int opt;
  bool write_to_file = false;
  bool simulate_render_to_file = false;
  double simulate_time = 0.0;
  size_t w = 0, h = 0;
  string filename;
  string particle_path;
  while ( (opt = getopt(argc, argv, "s:l:t:m:e:h:f:r:p:d:")) != -1 ) {  // for each option...
    switch ( opt ) {
    case 'f':
        write_to_file = true;
        filename  = string(optarg);
        break;
    case 'r':
        w = atoi(argv[optind-1]);
        h = atoi(argv[optind]);
        optind++;
        break;
    case 's':
        config.pathtracer_ns_aa = atoi(optarg);
        break;
    case 'l':
        config.pathtracer_ns_area_light = atoi(optarg);
        break;
    case 't':
        config.pathtracer_num_threads = atoi(optarg);
        break;
    case 'm':
        config.pathtracer_max_ray_depth = atoi(optarg);
        break;
    case 'e':
        config.pathtracer_envmap = load_exr(optarg);
        break;
    case 'p':
        particle_path = string(optarg);
        msg(particle_path);
        break;
    case 'd':
        simulate_render_to_file = true;
        simulate_time = (double) atoi(optarg);
        break;
    default:
        usage(argv[0]);
        return 1;
    }
  }

  // print usage if no argument given
  if (optind >= argc) {
    usage(argv[0]);
    return 1;
  }

  string sceneFilePath = argv[optind];
  msg("Input scene file: " << sceneFilePath);

  // parse scene
  Collada::SceneInfo *sceneInfo = new Collada::SceneInfo();
  if (Collada::ColladaParser::load(sceneFilePath.c_str(), sceneInfo) < 0) {
    delete sceneInfo;
    exit(0);
  }


  // create application
  Application *app  = new Application(config, !write_to_file && !simulate_render_to_file);

  // write straight to file without opening a window if -f option provided
  if (write_to_file) {
    app->init();
    app->load(sceneInfo);
    // load particles
    app->load_particles(particle_path.c_str());
    delete sceneInfo;

    if (w && h)
      app->resize(w, h);

    app->render_to_file(filename); return 0;
    return 0;
  }

  if (simulate_render_to_file) {
    app->init();
    app->load(sceneInfo);
    // load particles
    app->load_particles(particle_path.c_str());
    delete sceneInfo;

    if (w && h)
      app->resize(w, h);

    app->fluid_simulate_render_to_file(simulate_time); return 0;
    return 0;
  }

  // create viewer
  Viewer viewer = Viewer();

  // set renderer
  viewer.set_renderer(app);

  // init viewer
  viewer.init();

  // load scene
  app->load(sceneInfo);

  // load particles
  app->load_particles(particle_path.c_str());
  // msg(app->particles->ps[0]->origin());
  // msg(app->particles->ps[0]->velocity);
  // msg(app->particles->ps[0]->radius());
  delete sceneInfo;

  // start viewer
  viewer.start();

  return 0;

}
