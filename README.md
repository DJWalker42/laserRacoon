# laserRacoon
A mostly object orientated library providing C++ code for physics computations

The library can be constructed from the header and source files found in the 'inc' and 'src' directories respectively. As provided the Makefile will build the library using the g++ compiler the output wil be placed into a directory called 'lib' that it will also create (tested on Ubuntu 14.04). For Windows and other operating systems using an IDE add the source files to your project and point the IDE to the location of the header files.

For the Visualisation module to work you will need a copy of OpenCV. Ubuntu 14.04 comes with a version of OpenCV2. For other OSs you may have to download a version from the offical website. This code will work with either OpenCV2 or OpenCV3.  Remember to point your IDE to the location of the OpenCV headers. If you so wish you may exclude the visulisation module from the library; you will just have to plot any data elsewhere.

In the 'progs' directory is a list of example programs using the library; the majority of these use the Visulisation module. For Linux users there is a MakeFile for your convenience; this will automatically link to the required OpenCV library files for the Visulisation module (at least it will on Ubuntu 14.04). For others using an IDE you will have to provide the relevant information to your IDE to link with the correct OpenCV libraries namely, opencv_core2XX and opencv_highgui2XX where the Xs refer to the minor version you're using. For OpenCV3 individual libraries have been replaced with an over-arching single file called opencv_world3XX. 


Note that this code was written with education and learning in mind so do not expect a fully optimised library. For example the linear algebra solvers (Cholesky, LU decomposition, etc) are much better done elsewhere (Google Eigen3 for an example) but the code as presented is intended to be as clear and straightforward as possible.

Why laserRacoon? Well lasers are to do with physics and racoons are cute. Plus now you've all got an image racoons with lasers strapped to their furry little backs. It's memorable. 