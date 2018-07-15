
LIBPATH = ../../../applications/webgl/lib/

all: cpp emc gen min

cpp:
	g++ -I $(LIBPATH) -o matrix matrix.cpp -std=c++11

emc:
	emcc -I $(LIBPATH) -O2 -o matrix.js matrix.cpp -s EXPORTED_FUNCTIONS="['_process']"

gen:
	g++ -I $(LIBPATH) -o generator generator.cpp

min: emc
	node $(LIBPATH)UglifyJS/bin/uglifyjs matrix.js > matrix.min.js
	sed 's/matrix\.js/matrix.min.js/g' matrix.html > matrix.min.html
