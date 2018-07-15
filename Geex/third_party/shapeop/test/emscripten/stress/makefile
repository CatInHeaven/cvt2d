
LIBPATH = ../../../applications/webgl/lib/

all: normal minified

normal:
	emcc -I $(LIBPATH) -O2 -o stress.js stress.cpp -s EXPORTED_FUNCTIONS="['_init']"

minified: normal
	node $(LIBPATH)UglifyJS/bin/uglifyjs stress.js > stress.min.js
	sed 's/stress\.js/stress.min.js/g' stress.html > stress.min.html
