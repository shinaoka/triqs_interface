# ALPS/CT-HYB solver for TRIQS applications
This is a Python wapper of the efficient and general ALPSCore/CT-HYB solver for TRIQS applications.

## Installation
Please follow the following steps.
You must use the same C++ compiler and the same C++ standard (i.e., C++14) to build TRIQS, ALPSCore, ALPSCore/CT-HYB, and this wapper. We strongly recommend to use GCC or clang compilers because TRIQS uses very new features of C++. As of now, we did not successfully build TRIQS using Intel C++ compiler.

### TRIQS
TRIQS and ALPSCore are independent sets of libraries. But, we strongly recommed you to install TRIQS first, because it requires the newer C++ standard. As of now (12/14/2016), TRIQS requires C++14.

TRIQS requires TRIQS-dependent applications to be compiled with exactly the same C++ compiler.
We strongly recommend you to explicitly pass the full path of your C++ compiler to CMake. TRIQS even cares about differencies between the full path and a relative path.<br><br>
Example: cmake -DCMAKE\_CXX\_COMPILER=/opt/local/bin/mpicxx-openmpi-clang38 ...

### ALPSCore
Then, please install ALPSCore. Although ALPSCore compiles with C++03, you must use the same C++ standard (and the compiler) as TRIQS for compatibility.

### ALPSCore/CT-HYB
ALPSCore/CT-HYB is an application depending on ALPSCore. You must use the same C++ compiler and C++ standard as TRIQS and ALPSCore. A C++ interface will be installed into the installation directory as a shared library as well.

### Python wrapper
Now, you are ready to install the wrapper. Please clone the repository and move into a build directory.
```
export ALPSCoreCTHYB_DIR=/installation/directory/of/alpscthyb
cmake  -DTRIQS_PATH=/installation/directory/of/triqs -DCMAKE_CXX_COMPILER==/opt/local/bin/mpicxx-openmpi-clang38 /path/to/source/directory
```

Then, you can make and make install.
The wapper will be installed into the TRIQS installation directory, and will be visible from the TRIQS Python interface "pytriqs" immediately.


## How to use
Now, you can call ALPSCore/CT-HYB from pytriqs in a very similar manner to TRIQS/cthyb.
You just need to import the solver instead of TRIQS/cthyb as follows.

```
from pytriqs.applications.impurity_solvers.alps_cthyb import Solver
```

You can create a solver object "S" in exactly the same manner to TRIQS/cthyb.
Then, you run the solver. Please pass the maximum simulation time (in units of second) as an integer to solver.

```
S.solve(h_int = U * n('up',0) * n('down',0), max_time = 60, perform_post_proc = True)
```

The wrapper does not support the full functionalities of the ALPSCore/CT-HYB.
The single-particle Green's function is measured using Legendre polynomials.
After the measurement process,
the imaginary-time and Matsubara-frequency data are computed from the Legendre basis data.
If you set perform\_post\_proc = True, the self-energy will be computed in addition.

An example Python script for a single-site impurity model is available [here](https://github.com/shinaoka/triqs_interface/blob/master/samples/aim_alps.py).

## License
This application is licensed under GPLv3.

## Acknowlegements
We acknowlege Emanuel Gull, Olivier Parcollet, Michel Ferrero and all other contributors to ALPS and TRIQS.
We also acknowledge Junya Otsuki for testing the code.
We took some source code from TRIQS/cthyb.
