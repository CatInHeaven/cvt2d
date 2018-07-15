
LIBPATH = ../../../applications/webgl/lib/

all: normal minified

normal:
	emcc -I $(LIBPATH) -O2 -o webgl.js webgl.cpp -s EXPORTED_FUNCTIONS="['_init', '_loadBegin', '_loadVertex', '_loadEnd', '_mouseMoved']"

minified: normal
	node $(LIBPATH)UglifyJS/bin/uglifyjs webgl.js > webgl.min.js
	sed 's/webgl\.js/webgl.min.js/g' webgl.html > webgl.min.html
