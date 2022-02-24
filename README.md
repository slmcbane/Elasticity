This repository contains the basic primitives needed to solve an isotropic
linear elasticity problem. You probably want to include "elasticity.hpp",
"read\_mesh.hpp", and perhaps "mesh\_tools.hpp" in order to use the code.
The validation problem in "validation/manufactured.cpp" provides guidance
on how to use the primitives provided here as it solves a problem with
the Lame parameters defined as scalar fields (not constants), non-homogeneous
Dirichlet boundary conditions, and non-zero traction boundary conditions.

License
-------
The contents of this repository are Copyright 2022 Sean McBane, under the
terms of the MIT license:

Permission is hereby granted, free of charge, to any person obtaining a copy of this
software and associated documentation files (the "Software"), to deal in the
Software without restriction, including without limitation the rights to use, copy,
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so, subject to the
following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
