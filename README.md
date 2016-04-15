# Fluid renderer with subsurface structure

## Summary
In this project, we want to support fluid simulation and display the simulated result with rendering algorithms. In order to achieve more realistic rendering, we also plan to implement subsurface scattering of fluid.

## Team members
+ Chuqian Li
+ Tongzhou Wang
+ Huirong Zhu

## Background

+ [Position Based Fluids](http://mmacklin.com/pbf_sig_preprint.pdf)
+ [A Practical Model for Subsurface Light Transport](https://graphics.stanford.edu/papers/bssrdf/bssrdf.pdf)
+ [Fluid Surface Reconstruction
from Particles](https://www.cs.ubc.ca/~rbridson/docs/brentw_msc.pdf)

## Resources
We decide to implement the project in C++ basing on code from assignment 3.

## Goals
### Baseline plan
+ Display the geometry of result from fluid simulation in edit mode.
+ Render an image of fluid moving in the Cornell box.
+ Render an image of fluid with reflected light showing the subsurface structure.

### Aspirational plan
+ High quality ray tracing with fluid will take a long time. If time permits, we want to use advanced algorithm to speed up rendering, or investigate photon mapping to get better result.

## Schedule
+ 1 week: 
	+ Implement fluid simulation, investigate ways to specify a fluid model to renderer. 
	+ Achieve baseline plan #1.
+ 0.5 week: 
	+ Implement fluid surfacing. 
	+ Achieve baseline plan #2.
+ 0.5 week: 
	+ Implement subsurface structure.
	+ Achieve baseline plan #3.
+ 0.5 week: prepare presentation.
+ 0.5 week: flexible time in case any of the above parts takes more time than expected.