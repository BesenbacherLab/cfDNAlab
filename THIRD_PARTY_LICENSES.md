# Third-party licenses

## SciPy 1.13.0

The SciPy-derived implementation is in
`src/commands/overlapping_lengths_correction/scipy_bfgs.rs`, with adapted regression tests in
`src/commands/overlapping_lengths_correction/scipy_bfgs_tests.rs`. It translates parts of SciPy
1.13.0's `scipy/optimize/_optimize.py`, `_linesearch.py`, `_dcsrch.py`, and the
numerical-differentiation path.

Any defects or behavioral differences introduced by this Rust translation are the responsibility
of the cfDNAlab developers, not the SciPy or MINPACK contributors.

Copyright (c) 2001-2002 Enthought, Inc. 2003-2024, SciPy Developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following
   disclaimer in the documentation and/or other materials provided
   with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

### Additional `_optimize.py` notice

```text
******NOTICE***************
SciPy's scipy/optimize/_optimize.py module was originally written by Travis E. Oliphant

You may copy and use that SciPy module as you see fit with no
guarantee implied provided you keep this notice in all copies.
*****END NOTICE************
```

### DCSRCH and DCSTEP provenance

SciPy 1.13.0's `_dcsrch.py` records that its Python implementation was ported from the
MINPACK-1 and MINPACK-2 Fortran routines:

```text
MINPACK-1 Project. June 1983.
Argonne National Laboratory.
Jorge J. More' and David J. Thuente.

MINPACK-2 Project. November 1993.
Argonne National Laboratory and University of Minnesota.
Brett M. Averick, Richard G. Carter, and Jorge J. More'.
```
