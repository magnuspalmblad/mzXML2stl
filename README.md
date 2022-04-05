# mzXML2stl
Coverter from mzXML to STL format, e.g. for 3D printing LC-MS data. This is an experimental tool, but may still be useful for 3D printing LC-MS data for teaching and artistic purposes. It can also be used to visualize LC-MS data using generic 3D visualization/CAD software.

# building mzXML2stl

mzXML2stl can be compiled using GCC on Windows/Cygwin and Linux using:

```
gcc -o mzXML2stl mzXML2stl.c base64.c ramp.c -I. -lgd -lm -lz
```
