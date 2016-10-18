# laserRacoon
## A mostly object orientated library providing C++ code for physics computations.

### Notice: This code was written with education and learning in mind so do not expect a fully optimised library. For example, the linear algebra solvers (Cholesky, LU decomposition, etc) are much better done elsewhere ([Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page) for an example) but the code as presented is intended to be as clear and straightforward as possible.

The library can be constructed from the header and source files found in the _inc_ and _src_ directories respectively. As provided the Makefile will build the library using the g++ compiler the output will be placed into a directory called _lib_ that will be created on the first call to make (tested on Ubuntu 14.04 and Mac OS X El Capitan). The intermediate object files are saved in directory called 'obj' that also will be created on the first call to make. For Windows and other operating systems were you might be using an IDE, add the source files to your project and point the IDE to the location of the header files. The IDE will manage intermediate files and you can specify where the output it built.

For the Visualisation module to work you will need a copy of OpenCV. Ubuntu 14.04 comes with a version of OpenCV2. For other OSs you may have to download a version from the official [OpenCV website](http://opencv.org/). This code will work with either OpenCV2 or OpenCV3.  Remember to point your IDE to the location of the OpenCV headers. If you so wish you may exclude the Visulisation module from the library; you will just have to plot any data elsewhere.

In the 'progs' directory is a list of example programs using the library; the majority of these use the Visulisation module. For Linux users there is a MakeFile for your convenience; this will automatically link to the required OpenCV library files for the Visulisation module (tested on Ubuntu 14.04 and Mac OS X El Capitan). For others using an IDE you will have to provide the relevant information to your IDE to link with the correct OpenCV libraries namely, opencv_core, opencv_highgui, opencv_imgproc, and opencv_imgcodecs (when saving the viewer plot as an image). On Windows these library names append the version number and end in a .lib extension. For example if you have OpenCV2.4.9 on Windows then the "core" library is called "opencv_core249.lib".

I wrote this library when trying to write a blog/tutorial on the C++ programming language and how to use it to write physics algorithms. That blog can be found [here](https://cppcomputationalphysics.wordpress.com/).


Why laserRacoon? Well lasers are to do with physics and racoons are cute. Plus now you've all got an image in your head of racoons with lasers strapped to their furry little backs. It's memorable.
