Armadillo C++ Linear Algebra Library
http://arma.sourceforge.net



Contents
========

 1: Introduction
 2: Citation Details
 3: License

 4: Requirements

 5: Linux and Mac OS X: Installation
 6: Linux and Mac OS X: Compiling & Linking

 7: Windows: Installation
 8: Windows: Compiling & Linking

 9: Support for OpenBLAS, Intel MKL and AMD ACML
10: Support for ATLAS

11: API Documentation
12: MEX Interface to Octave/Matlab

13: Bug Reports and Frequently Asked Questions
14: Related Software



1: Introduction
===============

Armadillo is a high quality C++ linear algebra library,
aiming towards a good balance between speed and ease of use.

It's useful for algorithm development directly in C++,
and/or quick conversion of research code into production environments.
The syntax (API) is deliberately similar to Matlab.

The library provides efficient classes for vectors, matrices and cubes,
as well as 200+ associated functions (eg. contiguous and non-contiguous
submatrix views).  Various matrix decompositions are provided through
integration with LAPACK, or one of its high performance drop-in replacements
(eg. OpenBLAS, Intel MKL, AMD ACML, Apple Accelerate framework, etc).

A sophisticated expression evaluator (via C++ template meta-programming)
automatically combines several operations (at compile time) to increase speed
and efficiency.

The library can be used for machine learning, pattern recognition, computer vision,
signal processing, bioinformatics, statistics, econometrics, etc.

Authors:
  Conrad Sanderson - http://conradsanderson.id.au
  Ryan Curtin      - http://ratml.org



2: Citation Details
===================

Please cite the following article if you use Armadillo in your
research and/or software. Citations are useful for the continued
development and maintenance of the library.

  Conrad Sanderson and Ryan Curtin.
  Armadillo: a template-based C++ library for linear algebra.
  Journal of Open Source Software, Vol. 1, pp. 26, 2016.



3: License
==========

Armadillo is available as open source, under the restrictions of the
Mozilla Public License (MPL) 2.0. See the "LICENSE.txt" file for details.



4: Requirements
===============

Armadillo makes extensive use of template meta-programming, recursive templates
and template based function overloading.  As such, C++ compilers which do not
fully implement the C++ standard may not work correctly.

The functionality of Armadillo is partly dependent on other libraries:
LAPACK, BLAS, ARPACK and SuperLU.  The LAPACK and BLAS libraries are
used for dense matrices, while the ARPACK and SuperLU libraries are
used for sparse matrices.  Armadillo can work without these libraries,
but its functionality will be reduced. In particular, basic functionality
will be available (eg. matrix addition and multiplication), but things
like eigen decomposition or matrix inversion will not be.
Matrix multiplication (mainly for big matrices) may not be as fast.

As Armadillo is a template library, we recommended that optimisation
is enabled during compilation of programs that use Armadillo.
For example, for GCC and Clang compilers use -O2 or -O3



5: Linux and Mac OS X: Installation
===================================

* Step 1:
  Ensure a C++ compiler is installed on your system.
  
  Caveat: on Mac OS X you will need to install Xcode
  and then type the following command in a terminal window:
  xcode-select --install
  
* Step 2:
  Ensure the CMake tool is installed on your system.
  You can download it from http://www.cmake.org
  or (preferably) install it using your package manager.
  
  On Linux-based systems, you can get CMake using yum, dnf, apt, aptitude, ...
  
  On Mac OS X systems, you can get CMake through MacPorts or Homebrew.
  
* Step 3:
  Ensure LAPACK and BLAS are installed on your system.
  On Mac OS X this is not necessary.
  
  For better performance, we recommend installing the OpenBLAS library.
  See http://xianyi.github.com/OpenBLAS/
  
  If you are using sparse matrices, also install ARPACK and SuperLU.
  Caveat: only SuperLU version 5.2 can be used!
  
  On Linux-based systems, the following libraries are recommended
  to be present: OpenBLAS, LAPACK, SuperLU and ARPACK.
  It is also necessary to install the corresponding development
  files for each library.  For example, when installing the "lapack"
  package, also install the "lapack-devel" or "lapack-dev" package.
  
* Step 4:
  Open a terminal window, change into the directory that was created
  by unpacking the armadillo archive, and type the following commands:
  
  cmake .
  make 
  
  The full stop separated from "cmake" by a space is important.
  CMake will detect which relevant libraries are installed on your system
  (eg. OpenBLAS, LAPACK, SuperLU, ARPACK, etc)
  and will modify Armadillo's configuration correspondingly.
  CMake will also generate a run-time armadillo library,
  which is a wrapper for all the detected libraries.
  
  If you need to re-run cmake, it's a good idea to first delete the
  "CMakeCache.txt" file (not "CMakeLists.txt").
  
  Caveat: out-of-tree builds are currently not fully supported;
  eg, creating a sub-directory called "build" and running cmake ..
  from within "build" is currently not supported.
  
* Step 5:
  If you have access to root/administrator/superuser privileges
  (ie. able to use "sudo"), type the following command:
  
  sudo make install
  
  If you don't have root/administrator/superuser privileges, 
  type the following command:
  
  make install DESTDIR=my_usr_dir
  
  where "my_usr_dir" is for storing C++ headers and library files.
  Caveat: make sure your C++ compiler is configured to use the
  "lib" and "include" sub-directories present within this directory.



6: Linux and Mac OS X: Compiling & Linking
==========================================

The "examples" directory contains several quick example programs
that use the Armadillo library.

In general, programs which use Armadillo are compiled along these lines:
  
  g++ example1.cpp -o example1 -O2 -larmadillo
  
If you want to use Armadillo without installation (not recommended),
compile along these lines:
  
  g++ example1.cpp -o example1 -O2 -I /home/blah/armadillo-7.200.3/include -DARMA_DONT_USE_WRAPPER -lblas -llapack
  
The above command line assumes that you have unpacked the armadillo archive into /home/blah/
You will need to adjust this for later versions of Armadillo (ie. change the 7.200.3 part)
and/or if you have unpacked the armadillo archive into a different directory.

Replace -lblas with -lopenblas if you have OpenBLAS.
On Mac OS X, replace -lblas -llapack with -framework Accelerate



7: Windows: Installation
========================

The installation is comprised of 3 steps:

* Step 1:
  Copy the entire "include" folder to a convenient location
  and tell your compiler to use that location for header files
  (in addition to the locations it uses already).
  Alternatively, you can use the "include" folder directly.
  
* Step 2:
  Modify "include/armadillo_bits/config.hpp" to indicate which
  libraries are currently available on your system. For example,
  if you have LAPACK, BLAS (or OpenBLAS), ARPACK and SuperLU present,
  uncomment the following lines:
  
  #define ARMA_USE_LAPACK
  #define ARMA_USE_BLAS
  #define ARMA_USE_ARPACK
  #define ARMA_USE_SUPERLU
  
  If you don't need sparse matrices, don't worry about ARPACK or SuperLU.
  
* Step 3:
  Configure your compiler to link with LAPACK and BLAS
  (and optionally ARPACK and SuperLU).



8: Windows: Compiling & Linking
===============================

Within the "examples" folder, we have included an MSVC project named "example1_win64"
which can be used to compile "example1.cpp".  The project needs to be compiled as a
64 bit program: the active solution platform must be set to x64, instead of win32.

The MSCV project was tested on 64 bit Windows 7 with Visual C++ 2012.
You may need to make adaptations for 32 bit systems, later versions of Windows
and/or the compiler.  For example, you may have to enable or disable
ARMA_BLAS_LONG and ARMA_BLAS_UNDERSCORE macros in "armadillo_bits/config.hpp".

The folder "examples/lib_win64" contains standard LAPACK and BLAS libraries compiled
for 64 bit Windows.  The compilation was done by a third party.  USE AT YOUR OWN RISK.
The compiled versions of LAPACK and BLAS were obtained from:
  http://ylzhao.blogspot.com.au/2013/10/blas-lapack-precompiled-binaries-for.html

You can find the original sources for standard BLAS and LAPACK at:
  http://www.netlib.org/blas/
  http://www.netlib.org/lapack/
  
Faster and/or alternative implementations of BLAS and LAPACK are available:
  http://xianyi.github.com/OpenBLAS/
  http://software.intel.com/en-us/intel-mkl/
  http://developer.amd.com/tools-and-sdks/cpu-development/amd-core-math-library-acml/
  http://icl.cs.utk.edu/lapack-for-windows/lapack/

The OpenBLAS, MKL and ACML libraries are generally the fastest.

For better performance, we recommend the following high-quality C++ compilers:
  GCC from MinGW:     http://www.mingw.org/
  GCC from CygWin:    http://www.cygwin.com/
  Intel C++ compiler: http://software.intel.com/en-us/intel-compilers/

For the GCC compiler, use version 4.6 or later.
For the Intel compiler, use version 13.0 or later.

For best results we also recommend using an operating system
that's more reliable and more suitable for heavy duty work,
such as Mac OS X, or various Linux-based systems:
  Ubuntu                    http://www.ubuntu.com/
  Debian                    http://www.debian.org/
  OpenSUSE                  http://www.opensuse.org/
  Fedora                    http://fedoraproject.org/
  Scientific Linux          http://www.scientificlinux.org/
  CentOS                    http://centos.org/
  Red Hat Enterprise Linux  http://www.redhat.com/



9: Support for OpenBLAS, Intel MKL and AMD ACML
===============================================

Armadillo can use OpenBLAS, or Intel Math Kernel Library (MKL),
or the AMD Core Math Library (ACML) as high-speed replacements
for BLAS and LAPACK.  In essence this involves linking with the
replacement libraries instead of BLAS and LAPACK.

You may need to make minor modifications to include/armadillo_bits/config.hpp
to make sure Armadillo uses the same integer sizes and style of function names
as used by MKL or ACML. Specifically, you may need comment or uncomment the
following defines:

ARMA_USE_WRAPPER
ARMA_BLAS_CAPITALS
ARMA_BLAS_UNDERSCORE
ARMA_BLAS_LONG
ARMA_BLAS_LONG_LONG

See the associated descriptions in include/armadillo_bits/config.hpp
for more information about the above defines.  

On Linux systems, MKL and ACML might be installed in a non-standard
location, such as /opt, which can cause problems during linking.
Before installing Armadillo, the system should know where the MKL or ACML
libraries are located. For example, "/opt/intel/mkl/lib/intel64/".
This can be achieved by setting the LD_LIBRARY_PATH environment variable,
or for a more permanent solution, adding the directory locations
to "/etc/ld.so.conf". It may also be possible to store a text file 
with the locations in the "/etc/ld.so.conf.d" directory.
For example, "/etc/ld.so.conf.d/mkl.conf".
If you modify "/etc/ld.so.conf" or create "/etc/ld.so.conf.d/mkl.conf",
you will need to run "/sbin/ldconfig" afterwards.

Example of the contents of "/etc/ld.so.conf.d/mkl.conf" on a RHEL 6 system,
where Intel MKL version 11.0.3 is installed in "/opt/intel":

/opt/intel/lib/intel64
/opt/intel/mkl/lib/intel64

The default installations of ACML 4.4.0 and MKL 10.2.2.025 are known 
to have issues with SELinux, which is turned on by default in Fedora
(and possibly RHEL). The problem may manifest itself during run-time,
where the run-time linker reports permission problems.
It is possible to work around the problem by applying an appropriate
SELinux type to all ACML and MKL libraries.

If you have ACML or MKL installed and they are persistently giving
you problems during linking, you can disable the support for them
by editing the "CMakeLists.txt" file, deleting "CMakeCache.txt" and
re-running the CMake based installation. Specifically, comment out
the lines containing:
  INCLUDE(ARMA_FindMKL)
  INCLUDE(ARMA_FindACMLMP)
  INCLUDE(ARMA_FindACML)



10: Support for ATLAS
=====================

Armadillo can use the ATLAS library for faster versions of
certain LAPACK and BLAS functions. Not all ATLAS functions are
currently used, and as such LAPACK should still be installed.

The minimum recommended version of ATLAS is 3.8.
Old versions (eg. 3.6) can produce incorrect results
as well as corrupting memory, leading to random crashes.

Users of older Ubuntu and Debian based systems should explicitly
check that ATLAS 3.6 is not installed.  It's better to
remove the old version and use the standard LAPACK library.



11: API Documentation
=====================

Documentation of functions and classes is available at:
  
  http://arma.sourceforge.net/docs.html
  
The documentation is also in the "docs.html" file in this folder,
which can be viewed with a web browser.



12: MEX Interface to Octave/Matlab
==================================

The "mex_interface" folder contains examples of how to interface
Octave/Matlab with C++ code that uses Armadillo matrices.



13: Bug Reports and Frequently Asked Questions
==============================================

Answers to frequently asked questions can be found at:

  http://arma.sourceforge.net/faq.html

This library has gone through extensive testing and
has been successfully used in production environments.
However, as with almost all software, it's impossible
to guarantee 100% correct functionality.

If you find a bug in the library (or the documentation),
we are interested in hearing about it. Please make a
_small_ and _self-contained_ program which exposes the bug,
and then send the program source (as well as the bug description)
to the developers.  The developers' contact details are at:

  http://arma.sourceforge.net/contact.html



14: Related Software
====================

* MLPACK: C++ library for machine learning and pattern recognition, built on top of Armadillo.
  http://mlpack.org
  
* libpca: C++ library for principal component analysis
  http://sourceforge.net/projects/libpca/
  
* ArmaNpy: interfaces Armadillo matrices with Python
  http://sourceforge.net/projects/armanpy/
  
* matlab2cpp: conversion of Matlab code to Armadillo based C++ code
  https://github.com/jonathf/matlab2cpp
