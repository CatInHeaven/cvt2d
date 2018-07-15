## Build

To build the matrix generator and the native C++ version of the test, any decent C++ compiler should be fine. My makefile assumes that g++ is installed, though it can be changed easily.

```
make gen
make cpp
./generator 1000 5000 > A.txt
./matrix A.txt L.txt
```

To generate the Javascript version, install the Emscripten framework. If you want to also minify the resulting JS file, Node.JS must be installed to run UglifyJS. The resulting HTML/JS files should work on any decent browser.

```
make emc
make min
```
