# laserRacoon Programs

If creating new programs in this directory that make use of the Viewer class then append the filename 
with '_ocv' so that the Makefile can include and link the required OpenCV headers and libraries. For example,
the 'doublePendulum' program uses the Viewer class so we give the source file the name 'doublePendulum_ocv.cpp'.

**NOTE: The header and library default locations for OpenCV are different for different platforms:**

- Mac OSX: 	/opt/local/include, /opt/local/lib
- Unix/Linux: /usr/local/include, /usr/local/lib (or /usr/include, /usr/lib)
- Windows: 	depends on CMAKE options when buidling/installing OpenCV

The Makefile as provided has the OpenCV location set to the Mac OSX default.


## Build and Install

To build the programs in this directory invoke 'make' on the command line. 
This will compile the source code into binaries, create './bin', then move those binaries into './bin'.

To compile individual programs invoke 'make binary_name' 

where the binary_name is the source filename without the .cpp extension, and in the case of OpenCV programs also
without the '_ocv'. 

