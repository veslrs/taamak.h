# taamak.h - Semi-analytical Truncated Elastic Half-space Model
`taamak.h` is a header only C library for calculating responses in a truncated multi-layer linear-elastic pavement system. Its features include:

- Evaluation of stresses and displacements in a trimmed homogeneous linear elastic pavement model, subject to a square patch stress.
- Evaluation of surface deflection on a truncated multi-layered linear elastic pavement model, subject to a square patch stress.
- Back-calculation of Effective Slope Angle (ESA) of a homogeneous/multi-layered linear elastic pavement model, utilising FWD deflection measurements tested near the edge/discontinuity of the pavement system.

## Documentation
This [documentation](./docs/documentation.pdf) consists a tutorial section with a few trivial examples, and `taamak.h` API reference. 

## Project in action
[BCRRA-26](https://github.com/veslrs/bcrra-26): Backcalculation of deflections obtained near asphalt pavement edges.
