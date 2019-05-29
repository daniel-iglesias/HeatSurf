# HeatSurf

HeatSurf is a particle tracer code developed originally for modelling the IFMIF/EVEDA beam stopping surface devices. It has been extended for its application within tokamaks, for modelling surfaces loads in divertors and neutral beam dumps 

# New Features!

  - Revolute topologies: Cylinder, cone and truncated cone surfaces
  - Open topologies: Plate and double plate surfaces

Output:
  - VTK surface structured grid
  - VTK 3D structured grid
  - Average density profiles in 2D and 1D
  - ANSYS heat flux density maps

### Tech

HeatFlux uses a number of open source projects to work properly:

* [VTK] - Visualization Toolkit from Kitware Inc.
* [LMX] - Linked Matrix Methods library
* [GSL] - GNU Scientific Library

And of course HeatFlux itself is open source with a [public repository][HeatFlux] on GitHub.

### Documentation

Reference documentation is in doc/html/index.html and in doc/latex/refman.pdf

### Development

Want to contribute? Great! Feel free to branch and merge, and get in touch in case you need help or are looking for joining up efforts.

#### Building for source
HeatSurf uses the CMake system for easy of cross-platform build. For example, in a linux system:
```sh
$ cmake -B<build directory> -H. -DCMAKE_BUILD_TYPE=<Debug|Release>
```
Generating pre-built zip archives for distribution:
```sh
$ cmake -Bbuild -H. -DCMAKE_BUILD_TYPE=Debug
```

License
----
Copyright (C) 2015 by Daniel Iglesias
Copyright (C) 2009-2012 by CIEMAT

HeatSurf is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2.1 of the License, or (at your option) any later version.

HEatSurf is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with HEatSurf. If not, see http://www.gnu.org/licenses/.

This software includes (in src/LMX/) the [LMX] library, which is distruted under the GNU LGPL v2.1 license. See COPYING and AUTHORS files in its respective directory.


[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


   [VTK]: <https://vtk.org/>
   [LMX]: <http://daniel-iglesias.github.io/lmx/>
   [GSL]: <https://www.gnu.org/software/gsl/>
   [HeatFlux]: <https://github.com/daniel-iglesias/HeatSurf.git>

