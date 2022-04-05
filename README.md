# mzXML2stl
Coverter from mzXML to STL format, e.g. for 3D printing LC-MS data. This is an experimental tool, but may still be useful for 3D printing LC-MS data for teaching and artistic purposes. It can also be used to visualize LC-MS data using generic 3D visualization/CAD software.

## building mzXML2stl

mzXML2stl can be compiled using GCC on Windows/Cygwin and Linux using:

```
gcc -o mzXML2stl mzXML2stl.c base64.c ramp.c -I. -lgd -lm -lz
```


## running mzXML2stl

mzXML2stl is run on the command line using:

```
mzXML2stl -i<mzXML file> -o<STL filename> -s<first scan>,<last scan> -m<lowest m/z>,<highest m/z>
```

mzXML2stl stretches each spectrum (scan) in the time dimension for a discrete, stepwise, look true to the real data rather than interpolating between spectra. mzXML2stl also creates a 5 mm high platform for the LC-MS data, the height of which can be adjusted in most CAD/3D printing software (or in the code).
