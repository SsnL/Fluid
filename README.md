# Fluid renderer with subsurface structure

## Summary
In this project, we want to support fluid simulation and display the simulated result with rendering algorithms. In order to achieve more realistic rendering, we also plan to implement subsurface scattering of fluid.

## Team members
+ Chuqian Li
+ Tongzhou Wang
+ Huirong Zhu

## Background
### Challenges
+ Fluid simulation is widely used by artists and producers in providing amazing visual effects. It's also used in engineering area in studying fluid dynamics.
+ However, fluid movement follows complicated physical model and thus requires a huge amount of computation power to simulate. Furthermore, we want to make the rendering of fluid more realistic by applying the subsurface scattering model. Thus, efficient simulation and rendering become critical in our project.

### Papers
+ [Position Based Fluids](http://mmacklin.com/pbf_sig_preprint.pdf)
+ [Position Based Fluids [slides]](http://mmacklin.com/pbf_slides.pdf)
+ [Fluid Surface Reconstruction
from Particles](https://www.cs.ubc.ca/~rbridson/docs/brentw_msc.pdf)

## Resources
We decide to implement the project in C++, and use assignment 3 as starting code.

## Goals
### Baseline plan
+ Display the geometry of results from fluid simulation in edit mode.
+ Render an image of fluid moving in the Cornell box.

### Aspirational plan
+ High quality ray tracing with fluid will take a long time. If time permits, we want to implement advanced algorithms to speed up rendering, or investigate photon mapping to get better rendering result.

## Schedule
+ 0.75 week:
	+ Implement fluid simulation, investigate ways to specify a fluid model to renderer.
	+ Achieve baseline plan #1.
+ 0.75 week:
	+ Implement fluid surfacing.
	+ Achieve baseline plan #2.
+ 0.5 week: prepare presentation.
+ 0.5 week: flexible time in case any of the above parts takes more time than expected.
