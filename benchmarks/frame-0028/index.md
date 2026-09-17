---
title: Distributed moments
thumbnail: img/isometric.png
---

`FrameLoad` is a simple class offering a clean API that implements generalized body loading  in frame elements.
This includes:

- Distributed moments (distr=`Heaviside`)
- Distributed forces  (distr=`Heaviside`)
- Point loads (distr=`Dirac`)
- Point moments (distr=`Dirac`)

All loads can:
- offset arbitrarily from the frame's centroid, 
- Defined in the local or global coordinate system
- Defined as follower or conservative loads.

The implementation is feature complete for displacement formulations, and partially complete for mixed (force) formulations.


For example, the problem of a distributed couple over an element can be implemented in two different ways:

1. Define a uniform distributed axial force, `force=[w,0,0]` with `offset=[0,0,h/2]`:
   ```python
   ops.eleLoad("Frame",
               "Heaviside",
               basis = "local",
               offset=[0,h/2,0],
               force = [w, 0, 0],
               pattern=1,
               elements=1
   )
   ```

2. Define a uniform force `force=[w,0,0]` and a uniform moment `couple=[0,-m,0]`

   ```python
   ops.eleLoad("Frame",
               "Heaviside",
               basis = "local",
               force = [w, 0, 0],
               couple= [0,-m,0],
               pattern=1,
               elements=1
   )
   ```

This fixes some existing bugs and defects in the current implementation:

- Conservative distributed loading works correctly with the corotational transformation
- Elements dont require privileged knowlege of the layout of an arbitrary "data" vector


## To-Dos

- Implement offsets for force formulation
