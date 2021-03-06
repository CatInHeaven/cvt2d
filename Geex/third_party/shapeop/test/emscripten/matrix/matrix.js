// The Module object: Our interface to the outside world. We import
// and export values on it, and do the work to get that through
// closure compiler if necessary. There are various ways Module can be used:
// 1. Not defined. We create it here
// 2. A function parameter, function(Module) { ..generated code.. }
// 3. pre-run appended it, var Module = {}; ..generated code..
// 4. External script tag defines var Module.
// We need to do an eval in order to handle the closure compiler
// case, where this code here is minified but Module was defined
// elsewhere (e.g. case 4 above). We also need to check if Module
// already exists (e.g. case 3 above).
// Note that if you want to run closure, and also to use Module
// after the generated code, you will need to define   var Module = {};
// before the code. Then that object will be used in the code, and you
// can continue to use Module afterwards as well.
var Module;
if (!Module) Module = (typeof Module !== 'undefined' ? Module : null) || {};

// Sometimes an existing Module object exists with properties
// meant to overwrite the default module functionality. Here
// we collect those properties and reapply _after_ we configure
// the current environment's defaults to avoid having to be so
// defensive during initialization.
var moduleOverrides = {};
for (var key in Module) {
  if (Module.hasOwnProperty(key)) {
    moduleOverrides[key] = Module[key];
  }
}

// The environment setup code below is customized to use Module.
// *** Environment setup code ***
var ENVIRONMENT_IS_NODE = typeof process === 'object' && typeof require === 'function';
var ENVIRONMENT_IS_WEB = typeof window === 'object';
var ENVIRONMENT_IS_WORKER = typeof importScripts === 'function';
var ENVIRONMENT_IS_SHELL = !ENVIRONMENT_IS_WEB && !ENVIRONMENT_IS_NODE && !ENVIRONMENT_IS_WORKER;

if (ENVIRONMENT_IS_NODE) {
  // Expose functionality in the same simple way that the shells work
  // Note that we pollute the global namespace here, otherwise we break in node
  if (!Module['print']) Module['print'] = function print(x) {
    process['stdout'].write(x + '\n');
  };
  if (!Module['printErr']) Module['printErr'] = function printErr(x) {
    process['stderr'].write(x + '\n');
  };

  var nodeFS = require('fs');
  var nodePath = require('path');

  Module['read'] = function read(filename, binary) {
    filename = nodePath['normalize'](filename);
    var ret = nodeFS['readFileSync'](filename);
    // The path is absolute if the normalized version is the same as the resolved.
    if (!ret && filename != nodePath['resolve'](filename)) {
      filename = path.join(__dirname, '..', 'src', filename);
      ret = nodeFS['readFileSync'](filename);
    }
    if (ret && !binary) ret = ret.toString();
    return ret;
  };

  Module['readBinary'] = function readBinary(filename) { return Module['read'](filename, true) };

  Module['load'] = function load(f) {
    globalEval(read(f));
  };

  Module['arguments'] = process['argv'].slice(2);

  module['exports'] = Module;
}
else if (ENVIRONMENT_IS_SHELL) {
  if (!Module['print']) Module['print'] = print;
  if (typeof printErr != 'undefined') Module['printErr'] = printErr; // not present in v8 or older sm

  if (typeof read != 'undefined') {
    Module['read'] = read;
  } else {
    Module['read'] = function read() { throw 'no read() available (jsc?)' };
  }

  Module['readBinary'] = function readBinary(f) {
    return read(f, 'binary');
  };

  if (typeof scriptArgs != 'undefined') {
    Module['arguments'] = scriptArgs;
  } else if (typeof arguments != 'undefined') {
    Module['arguments'] = arguments;
  }

  this['Module'] = Module;

  eval("if (typeof gc === 'function' && gc.toString().indexOf('[native code]') > 0) var gc = undefined"); // wipe out the SpiderMonkey shell 'gc' function, which can confuse closure (uses it as a minified name, and it is then initted to a non-falsey value unexpectedly)
}
else if (ENVIRONMENT_IS_WEB || ENVIRONMENT_IS_WORKER) {
  Module['read'] = function read(url) {
    var xhr = new XMLHttpRequest();
    xhr.open('GET', url, false);
    xhr.send(null);
    return xhr.responseText;
  };

  if (typeof arguments != 'undefined') {
    Module['arguments'] = arguments;
  }

  if (typeof console !== 'undefined') {
    if (!Module['print']) Module['print'] = function print(x) {
      console.log(x);
    };
    if (!Module['printErr']) Module['printErr'] = function printErr(x) {
      console.log(x);
    };
  } else {
    // Probably a worker, and without console.log. We can do very little here...
    var TRY_USE_DUMP = false;
    if (!Module['print']) Module['print'] = (TRY_USE_DUMP && (typeof(dump) !== "undefined") ? (function(x) {
      dump(x);
    }) : (function(x) {
      // self.postMessage(x); // enable this if you want stdout to be sent as messages
    }));
  }

  if (ENVIRONMENT_IS_WEB) {
    window['Module'] = Module;
  } else {
    Module['load'] = importScripts;
  }
}
else {
  // Unreachable because SHELL is dependant on the others
  throw 'Unknown runtime environment. Where are we?';
}

function globalEval(x) {
  eval.call(null, x);
}
if (!Module['load'] == 'undefined' && Module['read']) {
  Module['load'] = function load(f) {
    globalEval(Module['read'](f));
  };
}
if (!Module['print']) {
  Module['print'] = function(){};
}
if (!Module['printErr']) {
  Module['printErr'] = Module['print'];
}
if (!Module['arguments']) {
  Module['arguments'] = [];
}
// *** Environment setup code ***

// Closure helpers
Module.print = Module['print'];
Module.printErr = Module['printErr'];

// Callbacks
Module['preRun'] = [];
Module['postRun'] = [];

// Merge back in the overrides
for (var key in moduleOverrides) {
  if (moduleOverrides.hasOwnProperty(key)) {
    Module[key] = moduleOverrides[key];
  }
}



// === Auto-generated preamble library stuff ===

//========================================
// Runtime code shared with compiler
//========================================

var Runtime = {
  stackSave: function () {
    return STACKTOP;
  },
  stackRestore: function (stackTop) {
    STACKTOP = stackTop;
  },
  forceAlign: function (target, quantum) {
    quantum = quantum || 4;
    if (quantum == 1) return target;
    if (isNumber(target) && isNumber(quantum)) {
      return Math.ceil(target/quantum)*quantum;
    } else if (isNumber(quantum) && isPowerOfTwo(quantum)) {
      return '(((' +target + ')+' + (quantum-1) + ')&' + -quantum + ')';
    }
    return 'Math.ceil((' + target + ')/' + quantum + ')*' + quantum;
  },
  isNumberType: function (type) {
    return type in Runtime.INT_TYPES || type in Runtime.FLOAT_TYPES;
  },
  isPointerType: function isPointerType(type) {
  return type[type.length-1] == '*';
},
  isStructType: function isStructType(type) {
  if (isPointerType(type)) return false;
  if (isArrayType(type)) return true;
  if (/<?\{ ?[^}]* ?\}>?/.test(type)) return true; // { i32, i8 } etc. - anonymous struct types
  // See comment in isStructPointerType()
  return type[0] == '%';
},
  INT_TYPES: {"i1":0,"i8":0,"i16":0,"i32":0,"i64":0},
  FLOAT_TYPES: {"float":0,"double":0},
  or64: function (x, y) {
    var l = (x | 0) | (y | 0);
    var h = (Math.round(x / 4294967296) | Math.round(y / 4294967296)) * 4294967296;
    return l + h;
  },
  and64: function (x, y) {
    var l = (x | 0) & (y | 0);
    var h = (Math.round(x / 4294967296) & Math.round(y / 4294967296)) * 4294967296;
    return l + h;
  },
  xor64: function (x, y) {
    var l = (x | 0) ^ (y | 0);
    var h = (Math.round(x / 4294967296) ^ Math.round(y / 4294967296)) * 4294967296;
    return l + h;
  },
  getNativeTypeSize: function (type) {
    switch (type) {
      case 'i1': case 'i8': return 1;
      case 'i16': return 2;
      case 'i32': return 4;
      case 'i64': return 8;
      case 'float': return 4;
      case 'double': return 8;
      default: {
        if (type[type.length-1] === '*') {
          return Runtime.QUANTUM_SIZE; // A pointer
        } else if (type[0] === 'i') {
          var bits = parseInt(type.substr(1));
          assert(bits % 8 === 0);
          return bits/8;
        } else {
          return 0;
        }
      }
    }
  },
  getNativeFieldSize: function (type) {
    return Math.max(Runtime.getNativeTypeSize(type), Runtime.QUANTUM_SIZE);
  },
  dedup: function dedup(items, ident) {
  var seen = {};
  if (ident) {
    return items.filter(function(item) {
      if (seen[item[ident]]) return false;
      seen[item[ident]] = true;
      return true;
    });
  } else {
    return items.filter(function(item) {
      if (seen[item]) return false;
      seen[item] = true;
      return true;
    });
  }
},
  set: function set() {
  var args = typeof arguments[0] === 'object' ? arguments[0] : arguments;
  var ret = {};
  for (var i = 0; i < args.length; i++) {
    ret[args[i]] = 0;
  }
  return ret;
},
  STACK_ALIGN: 8,
  getAlignSize: function (type, size, vararg) {
    // we align i64s and doubles on 64-bit boundaries, unlike x86
    if (!vararg && (type == 'i64' || type == 'double')) return 8;
    if (!type) return Math.min(size, 8); // align structures internally to 64 bits
    return Math.min(size || (type ? Runtime.getNativeFieldSize(type) : 0), Runtime.QUANTUM_SIZE);
  },
  calculateStructAlignment: function calculateStructAlignment(type) {
    type.flatSize = 0;
    type.alignSize = 0;
    var diffs = [];
    var prev = -1;
    var index = 0;
    type.flatIndexes = type.fields.map(function(field) {
      index++;
      var size, alignSize;
      if (Runtime.isNumberType(field) || Runtime.isPointerType(field)) {
        size = Runtime.getNativeTypeSize(field); // pack char; char; in structs, also char[X]s.
        alignSize = Runtime.getAlignSize(field, size);
      } else if (Runtime.isStructType(field)) {
        if (field[1] === '0') {
          // this is [0 x something]. When inside another structure like here, it must be at the end,
          // and it adds no size
          // XXX this happens in java-nbody for example... assert(index === type.fields.length, 'zero-length in the middle!');
          size = 0;
          if (Types.types[field]) {
            alignSize = Runtime.getAlignSize(null, Types.types[field].alignSize);
          } else {
            alignSize = type.alignSize || QUANTUM_SIZE;
          }
        } else {
          size = Types.types[field].flatSize;
          alignSize = Runtime.getAlignSize(null, Types.types[field].alignSize);
        }
      } else if (field[0] == 'b') {
        // bN, large number field, like a [N x i8]
        size = field.substr(1)|0;
        alignSize = 1;
      } else if (field[0] === '<') {
        // vector type
        size = alignSize = Types.types[field].flatSize; // fully aligned
      } else if (field[0] === 'i') {
        // illegal integer field, that could not be legalized because it is an internal structure field
        // it is ok to have such fields, if we just use them as markers of field size and nothing more complex
        size = alignSize = parseInt(field.substr(1))/8;
        assert(size % 1 === 0, 'cannot handle non-byte-size field ' + field);
      } else {
        assert(false, 'invalid type for calculateStructAlignment');
      }
      if (type.packed) alignSize = 1;
      type.alignSize = Math.max(type.alignSize, alignSize);
      var curr = Runtime.alignMemory(type.flatSize, alignSize); // if necessary, place this on aligned memory
      type.flatSize = curr + size;
      if (prev >= 0) {
        diffs.push(curr-prev);
      }
      prev = curr;
      return curr;
    });
    if (type.name_ && type.name_[0] === '[') {
      // arrays have 2 elements, so we get the proper difference. then we scale here. that way we avoid
      // allocating a potentially huge array for [999999 x i8] etc.
      type.flatSize = parseInt(type.name_.substr(1))*type.flatSize/2;
    }
    type.flatSize = Runtime.alignMemory(type.flatSize, type.alignSize);
    if (diffs.length == 0) {
      type.flatFactor = type.flatSize;
    } else if (Runtime.dedup(diffs).length == 1) {
      type.flatFactor = diffs[0];
    }
    type.needsFlattening = (type.flatFactor != 1);
    return type.flatIndexes;
  },
  generateStructInfo: function (struct, typeName, offset) {
    var type, alignment;
    if (typeName) {
      offset = offset || 0;
      type = (typeof Types === 'undefined' ? Runtime.typeInfo : Types.types)[typeName];
      if (!type) return null;
      if (type.fields.length != struct.length) {
        printErr('Number of named fields must match the type for ' + typeName + ': possibly duplicate struct names. Cannot return structInfo');
        return null;
      }
      alignment = type.flatIndexes;
    } else {
      var type = { fields: struct.map(function(item) { return item[0] }) };
      alignment = Runtime.calculateStructAlignment(type);
    }
    var ret = {
      __size__: type.flatSize
    };
    if (typeName) {
      struct.forEach(function(item, i) {
        if (typeof item === 'string') {
          ret[item] = alignment[i] + offset;
        } else {
          // embedded struct
          var key;
          for (var k in item) key = k;
          ret[key] = Runtime.generateStructInfo(item[key], type.fields[i], alignment[i]);
        }
      });
    } else {
      struct.forEach(function(item, i) {
        ret[item[1]] = alignment[i];
      });
    }
    return ret;
  },
  dynCall: function (sig, ptr, args) {
    if (args && args.length) {
      if (!args.splice) args = Array.prototype.slice.call(args);
      args.splice(0, 0, ptr);
      return Module['dynCall_' + sig].apply(null, args);
    } else {
      return Module['dynCall_' + sig].call(null, ptr);
    }
  },
  functionPointers: [],
  addFunction: function (func) {
    for (var i = 0; i < Runtime.functionPointers.length; i++) {
      if (!Runtime.functionPointers[i]) {
        Runtime.functionPointers[i] = func;
        return 2*(1 + i);
      }
    }
    throw 'Finished up all reserved function pointers. Use a higher value for RESERVED_FUNCTION_POINTERS.';
  },
  removeFunction: function (index) {
    Runtime.functionPointers[(index-2)/2] = null;
  },
  getAsmConst: function (code, numArgs) {
    // code is a constant string on the heap, so we can cache these
    if (!Runtime.asmConstCache) Runtime.asmConstCache = {};
    var func = Runtime.asmConstCache[code];
    if (func) return func;
    var args = [];
    for (var i = 0; i < numArgs; i++) {
      args.push(String.fromCharCode(36) + i); // $0, $1 etc
    }
    var source = Pointer_stringify(code);
    if (source[0] === '"') {
      // tolerate EM_ASM("..code..") even though EM_ASM(..code..) is correct
      if (source.indexOf('"', 1) === source.length-1) {
        source = source.substr(1, source.length-2);
      } else {
        // something invalid happened, e.g. EM_ASM("..code($0)..", input)
        abort('invalid EM_ASM input |' + source + '|. Please use EM_ASM(..code..) (no quotes) or EM_ASM({ ..code($0).. }, input) (to input values)');
      }
    }
    try {
      var evalled = eval('(function(' + args.join(',') + '){ ' + source + ' })'); // new Function does not allow upvars in node
    } catch(e) {
      Module.printErr('error in executing inline EM_ASM code: ' + e + ' on: \n\n' + source + '\n\nwith args |' + args + '| (make sure to use the right one out of EM_ASM, EM_ASM_ARGS, etc.)');
      throw e;
    }
    return Runtime.asmConstCache[code] = evalled;
  },
  warnOnce: function (text) {
    if (!Runtime.warnOnce.shown) Runtime.warnOnce.shown = {};
    if (!Runtime.warnOnce.shown[text]) {
      Runtime.warnOnce.shown[text] = 1;
      Module.printErr(text);
    }
  },
  funcWrappers: {},
  getFuncWrapper: function (func, sig) {
    assert(sig);
    if (!Runtime.funcWrappers[func]) {
      Runtime.funcWrappers[func] = function dynCall_wrapper() {
        return Runtime.dynCall(sig, func, arguments);
      };
    }
    return Runtime.funcWrappers[func];
  },
  UTF8Processor: function () {
    var buffer = [];
    var needed = 0;
    this.processCChar = function (code) {
      code = code & 0xFF;

      if (buffer.length == 0) {
        if ((code & 0x80) == 0x00) {        // 0xxxxxxx
          return String.fromCharCode(code);
        }
        buffer.push(code);
        if ((code & 0xE0) == 0xC0) {        // 110xxxxx
          needed = 1;
        } else if ((code & 0xF0) == 0xE0) { // 1110xxxx
          needed = 2;
        } else {                            // 11110xxx
          needed = 3;
        }
        return '';
      }

      if (needed) {
        buffer.push(code);
        needed--;
        if (needed > 0) return '';
      }

      var c1 = buffer[0];
      var c2 = buffer[1];
      var c3 = buffer[2];
      var c4 = buffer[3];
      var ret;
      if (buffer.length == 2) {
        ret = String.fromCharCode(((c1 & 0x1F) << 6)  | (c2 & 0x3F));
      } else if (buffer.length == 3) {
        ret = String.fromCharCode(((c1 & 0x0F) << 12) | ((c2 & 0x3F) << 6)  | (c3 & 0x3F));
      } else {
        // http://mathiasbynens.be/notes/javascript-encoding#surrogate-formulae
        var codePoint = ((c1 & 0x07) << 18) | ((c2 & 0x3F) << 12) |
                        ((c3 & 0x3F) << 6)  | (c4 & 0x3F);
        ret = String.fromCharCode(
          Math.floor((codePoint - 0x10000) / 0x400) + 0xD800,
          (codePoint - 0x10000) % 0x400 + 0xDC00);
      }
      buffer.length = 0;
      return ret;
    }
    this.processJSString = function processJSString(string) {
      /* TODO: use TextEncoder when present,
        var encoder = new TextEncoder();
        encoder['encoding'] = "utf-8";
        var utf8Array = encoder['encode'](aMsg.data);
      */
      string = unescape(encodeURIComponent(string));
      var ret = [];
      for (var i = 0; i < string.length; i++) {
        ret.push(string.charCodeAt(i));
      }
      return ret;
    }
  },
  getCompilerSetting: function (name) {
    throw 'You must build with -s RETAIN_COMPILER_SETTINGS=1 for Runtime.getCompilerSetting or emscripten_get_compiler_setting to work';
  },
  stackAlloc: function (size) { var ret = STACKTOP;STACKTOP = (STACKTOP + size)|0;STACKTOP = (((STACKTOP)+7)&-8); return ret; },
  staticAlloc: function (size) { var ret = STATICTOP;STATICTOP = (STATICTOP + size)|0;STATICTOP = (((STATICTOP)+7)&-8); return ret; },
  dynamicAlloc: function (size) { var ret = DYNAMICTOP;DYNAMICTOP = (DYNAMICTOP + size)|0;DYNAMICTOP = (((DYNAMICTOP)+7)&-8); if (DYNAMICTOP >= TOTAL_MEMORY) enlargeMemory();; return ret; },
  alignMemory: function (size,quantum) { var ret = size = Math.ceil((size)/(quantum ? quantum : 8))*(quantum ? quantum : 8); return ret; },
  makeBigInt: function (low,high,unsigned) { var ret = (unsigned ? ((+((low>>>0)))+((+((high>>>0)))*(+4294967296))) : ((+((low>>>0)))+((+((high|0)))*(+4294967296)))); return ret; },
  GLOBAL_BASE: 8,
  QUANTUM_SIZE: 4,
  __dummy__: 0
}


Module['Runtime'] = Runtime;









//========================================
// Runtime essentials
//========================================

var __THREW__ = 0; // Used in checking for thrown exceptions.

var ABORT = false; // whether we are quitting the application. no code should run after this. set in exit() and abort()
var EXITSTATUS = 0;

var undef = 0;
// tempInt is used for 32-bit signed values or smaller. tempBigInt is used
// for 32-bit unsigned values or more than 32 bits. TODO: audit all uses of tempInt
var tempValue, tempInt, tempBigInt, tempInt2, tempBigInt2, tempPair, tempBigIntI, tempBigIntR, tempBigIntS, tempBigIntP, tempBigIntD, tempDouble, tempFloat;
var tempI64, tempI64b;
var tempRet0, tempRet1, tempRet2, tempRet3, tempRet4, tempRet5, tempRet6, tempRet7, tempRet8, tempRet9;

function assert(condition, text) {
  if (!condition) {
    abort('Assertion failed: ' + text);
  }
}

var globalScope = this;

// C calling interface. A convenient way to call C functions (in C files, or
// defined with extern "C").
//
// Note: LLVM optimizations can inline and remove functions, after which you will not be
//       able to call them. Closure can also do so. To avoid that, add your function to
//       the exports using something like
//
//         -s EXPORTED_FUNCTIONS='["_main", "_myfunc"]'
//
// @param ident      The name of the C function (note that C++ functions will be name-mangled - use extern "C")
// @param returnType The return type of the function, one of the JS types 'number', 'string' or 'array' (use 'number' for any C pointer, and
//                   'array' for JavaScript arrays and typed arrays; note that arrays are 8-bit).
// @param argTypes   An array of the types of arguments for the function (if there are no arguments, this can be ommitted). Types are as in returnType,
//                   except that 'array' is not possible (there is no way for us to know the length of the array)
// @param args       An array of the arguments to the function, as native JS values (as in returnType)
//                   Note that string arguments will be stored on the stack (the JS string will become a C string on the stack).
// @return           The return value, as a native JS value (as in returnType)
function ccall(ident, returnType, argTypes, args) {
  return ccallFunc(getCFunc(ident), returnType, argTypes, args);
}
Module["ccall"] = ccall;

// Returns the C function with a specified identifier (for C++, you need to do manual name mangling)
function getCFunc(ident) {
  try {
    var func = Module['_' + ident]; // closure exported function
    if (!func) func = eval('_' + ident); // explicit lookup
  } catch(e) {
  }
  assert(func, 'Cannot call unknown function ' + ident + ' (perhaps LLVM optimizations or closure removed it?)');
  return func;
}

// Internal function that does a C call using a function, not an identifier
function ccallFunc(func, returnType, argTypes, args) {
  var stack = 0;
  function toC(value, type) {
    if (type == 'string') {
      if (value === null || value === undefined || value === 0) return 0; // null string
      value = intArrayFromString(value);
      type = 'array';
    }
    if (type == 'array') {
      if (!stack) stack = Runtime.stackSave();
      var ret = Runtime.stackAlloc(value.length);
      writeArrayToMemory(value, ret);
      return ret;
    }
    return value;
  }
  function fromC(value, type) {
    if (type == 'string') {
      return Pointer_stringify(value);
    }
    assert(type != 'array');
    return value;
  }
  var i = 0;
  var cArgs = args ? args.map(function(arg) {
    return toC(arg, argTypes[i++]);
  }) : [];
  var ret = fromC(func.apply(null, cArgs), returnType);
  if (stack) Runtime.stackRestore(stack);
  return ret;
}

// Returns a native JS wrapper for a C function. This is similar to ccall, but
// returns a function you can call repeatedly in a normal way. For example:
//
//   var my_function = cwrap('my_c_function', 'number', ['number', 'number']);
//   alert(my_function(5, 22));
//   alert(my_function(99, 12));
//
function cwrap(ident, returnType, argTypes) {
  var func = getCFunc(ident);
  return function() {
    return ccallFunc(func, returnType, argTypes, Array.prototype.slice.call(arguments));
  }
}
Module["cwrap"] = cwrap;

// Sets a value in memory in a dynamic way at run-time. Uses the
// type data. This is the same as makeSetValue, except that
// makeSetValue is done at compile-time and generates the needed
// code then, whereas this function picks the right code at
// run-time.
// Note that setValue and getValue only do *aligned* writes and reads!
// Note that ccall uses JS types as for defining types, while setValue and
// getValue need LLVM types ('i8', 'i32') - this is a lower-level operation
function setValue(ptr, value, type, noSafe) {
  type = type || 'i8';
  if (type.charAt(type.length-1) === '*') type = 'i32'; // pointers are 32-bit
    switch(type) {
      case 'i1': HEAP8[(ptr)]=value; break;
      case 'i8': HEAP8[(ptr)]=value; break;
      case 'i16': HEAP16[((ptr)>>1)]=value; break;
      case 'i32': HEAP32[((ptr)>>2)]=value; break;
      case 'i64': (tempI64 = [value>>>0,(tempDouble=value,(+(Math_abs(tempDouble))) >= (+1) ? (tempDouble > (+0) ? ((Math_min((+(Math_floor((tempDouble)/(+4294967296)))), (+4294967295)))|0)>>>0 : (~~((+(Math_ceil((tempDouble - +(((~~(tempDouble)))>>>0))/(+4294967296))))))>>>0) : 0)],HEAP32[((ptr)>>2)]=tempI64[0],HEAP32[(((ptr)+(4))>>2)]=tempI64[1]); break;
      case 'float': HEAPF32[((ptr)>>2)]=value; break;
      case 'double': HEAPF64[((ptr)>>3)]=value; break;
      default: abort('invalid type for setValue: ' + type);
    }
}
Module['setValue'] = setValue;

// Parallel to setValue.
function getValue(ptr, type, noSafe) {
  type = type || 'i8';
  if (type.charAt(type.length-1) === '*') type = 'i32'; // pointers are 32-bit
    switch(type) {
      case 'i1': return HEAP8[(ptr)];
      case 'i8': return HEAP8[(ptr)];
      case 'i16': return HEAP16[((ptr)>>1)];
      case 'i32': return HEAP32[((ptr)>>2)];
      case 'i64': return HEAP32[((ptr)>>2)];
      case 'float': return HEAPF32[((ptr)>>2)];
      case 'double': return HEAPF64[((ptr)>>3)];
      default: abort('invalid type for setValue: ' + type);
    }
  return null;
}
Module['getValue'] = getValue;

var ALLOC_NORMAL = 0; // Tries to use _malloc()
var ALLOC_STACK = 1; // Lives for the duration of the current function call
var ALLOC_STATIC = 2; // Cannot be freed
var ALLOC_DYNAMIC = 3; // Cannot be freed except through sbrk
var ALLOC_NONE = 4; // Do not allocate
Module['ALLOC_NORMAL'] = ALLOC_NORMAL;
Module['ALLOC_STACK'] = ALLOC_STACK;
Module['ALLOC_STATIC'] = ALLOC_STATIC;
Module['ALLOC_DYNAMIC'] = ALLOC_DYNAMIC;
Module['ALLOC_NONE'] = ALLOC_NONE;

// allocate(): This is for internal use. You can use it yourself as well, but the interface
//             is a little tricky (see docs right below). The reason is that it is optimized
//             for multiple syntaxes to save space in generated code. So you should
//             normally not use allocate(), and instead allocate memory using _malloc(),
//             initialize it with setValue(), and so forth.
// @slab: An array of data, or a number. If a number, then the size of the block to allocate,
//        in *bytes* (note that this is sometimes confusing: the next parameter does not
//        affect this!)
// @types: Either an array of types, one for each byte (or 0 if no type at that position),
//         or a single type which is used for the entire block. This only matters if there
//         is initial data - if @slab is a number, then this does not matter at all and is
//         ignored.
// @allocator: How to allocate memory, see ALLOC_*
function allocate(slab, types, allocator, ptr) {
  var zeroinit, size;
  if (typeof slab === 'number') {
    zeroinit = true;
    size = slab;
  } else {
    zeroinit = false;
    size = slab.length;
  }

  var singleType = typeof types === 'string' ? types : null;

  var ret;
  if (allocator == ALLOC_NONE) {
    ret = ptr;
  } else {
    ret = [_malloc, Runtime.stackAlloc, Runtime.staticAlloc, Runtime.dynamicAlloc][allocator === undefined ? ALLOC_STATIC : allocator](Math.max(size, singleType ? 1 : types.length));
  }

  if (zeroinit) {
    var ptr = ret, stop;
    assert((ret & 3) == 0);
    stop = ret + (size & ~3);
    for (; ptr < stop; ptr += 4) {
      HEAP32[((ptr)>>2)]=0;
    }
    stop = ret + size;
    while (ptr < stop) {
      HEAP8[((ptr++)|0)]=0;
    }
    return ret;
  }

  if (singleType === 'i8') {
    if (slab.subarray || slab.slice) {
      HEAPU8.set(slab, ret);
    } else {
      HEAPU8.set(new Uint8Array(slab), ret);
    }
    return ret;
  }

  var i = 0, type, typeSize, previousType;
  while (i < size) {
    var curr = slab[i];

    if (typeof curr === 'function') {
      curr = Runtime.getFunctionIndex(curr);
    }

    type = singleType || types[i];
    if (type === 0) {
      i++;
      continue;
    }

    if (type == 'i64') type = 'i32'; // special case: we have one i32 here, and one i32 later

    setValue(ret+i, curr, type);

    // no need to look up size unless type changes, so cache it
    if (previousType !== type) {
      typeSize = Runtime.getNativeTypeSize(type);
      previousType = type;
    }
    i += typeSize;
  }

  return ret;
}
Module['allocate'] = allocate;

function Pointer_stringify(ptr, /* optional */ length) {
  // TODO: use TextDecoder
  // Find the length, and check for UTF while doing so
  var hasUtf = false;
  var t;
  var i = 0;
  while (1) {
    t = HEAPU8[(((ptr)+(i))|0)];
    if (t >= 128) hasUtf = true;
    else if (t == 0 && !length) break;
    i++;
    if (length && i == length) break;
  }
  if (!length) length = i;

  var ret = '';

  if (!hasUtf) {
    var MAX_CHUNK = 1024; // split up into chunks, because .apply on a huge string can overflow the stack
    var curr;
    while (length > 0) {
      curr = String.fromCharCode.apply(String, HEAPU8.subarray(ptr, ptr + Math.min(length, MAX_CHUNK)));
      ret = ret ? ret + curr : curr;
      ptr += MAX_CHUNK;
      length -= MAX_CHUNK;
    }
    return ret;
  }

  var utf8 = new Runtime.UTF8Processor();
  for (i = 0; i < length; i++) {
    t = HEAPU8[(((ptr)+(i))|0)];
    ret += utf8.processCChar(t);
  }
  return ret;
}
Module['Pointer_stringify'] = Pointer_stringify;

// Given a pointer 'ptr' to a null-terminated UTF16LE-encoded string in the emscripten HEAP, returns
// a copy of that string as a Javascript String object.
function UTF16ToString(ptr) {
  var i = 0;

  var str = '';
  while (1) {
    var codeUnit = HEAP16[(((ptr)+(i*2))>>1)];
    if (codeUnit == 0)
      return str;
    ++i;
    // fromCharCode constructs a character from a UTF-16 code unit, so we can pass the UTF16 string right through.
    str += String.fromCharCode(codeUnit);
  }
}
Module['UTF16ToString'] = UTF16ToString;

// Copies the given Javascript String object 'str' to the emscripten HEAP at address 'outPtr',
// null-terminated and encoded in UTF16LE form. The copy will require at most (str.length*2+1)*2 bytes of space in the HEAP.
function stringToUTF16(str, outPtr) {
  for(var i = 0; i < str.length; ++i) {
    // charCodeAt returns a UTF-16 encoded code unit, so it can be directly written to the HEAP.
    var codeUnit = str.charCodeAt(i); // possibly a lead surrogate
    HEAP16[(((outPtr)+(i*2))>>1)]=codeUnit;
  }
  // Null-terminate the pointer to the HEAP.
  HEAP16[(((outPtr)+(str.length*2))>>1)]=0;
}
Module['stringToUTF16'] = stringToUTF16;

// Given a pointer 'ptr' to a null-terminated UTF32LE-encoded string in the emscripten HEAP, returns
// a copy of that string as a Javascript String object.
function UTF32ToString(ptr) {
  var i = 0;

  var str = '';
  while (1) {
    var utf32 = HEAP32[(((ptr)+(i*4))>>2)];
    if (utf32 == 0)
      return str;
    ++i;
    // Gotcha: fromCharCode constructs a character from a UTF-16 encoded code (pair), not from a Unicode code point! So encode the code point to UTF-16 for constructing.
    if (utf32 >= 0x10000) {
      var ch = utf32 - 0x10000;
      str += String.fromCharCode(0xD800 | (ch >> 10), 0xDC00 | (ch & 0x3FF));
    } else {
      str += String.fromCharCode(utf32);
    }
  }
}
Module['UTF32ToString'] = UTF32ToString;

// Copies the given Javascript String object 'str' to the emscripten HEAP at address 'outPtr',
// null-terminated and encoded in UTF32LE form. The copy will require at most (str.length+1)*4 bytes of space in the HEAP,
// but can use less, since str.length does not return the number of characters in the string, but the number of UTF-16 code units in the string.
function stringToUTF32(str, outPtr) {
  var iChar = 0;
  for(var iCodeUnit = 0; iCodeUnit < str.length; ++iCodeUnit) {
    // Gotcha: charCodeAt returns a 16-bit word that is a UTF-16 encoded code unit, not a Unicode code point of the character! We must decode the string to UTF-32 to the heap.
    var codeUnit = str.charCodeAt(iCodeUnit); // possibly a lead surrogate
    if (codeUnit >= 0xD800 && codeUnit <= 0xDFFF) {
      var trailSurrogate = str.charCodeAt(++iCodeUnit);
      codeUnit = 0x10000 + ((codeUnit & 0x3FF) << 10) | (trailSurrogate & 0x3FF);
    }
    HEAP32[(((outPtr)+(iChar*4))>>2)]=codeUnit;
    ++iChar;
  }
  // Null-terminate the pointer to the HEAP.
  HEAP32[(((outPtr)+(iChar*4))>>2)]=0;
}
Module['stringToUTF32'] = stringToUTF32;

function demangle(func) {
  var i = 3;
  // params, etc.
  var basicTypes = {
    'v': 'void',
    'b': 'bool',
    'c': 'char',
    's': 'short',
    'i': 'int',
    'l': 'long',
    'f': 'float',
    'd': 'double',
    'w': 'wchar_t',
    'a': 'signed char',
    'h': 'unsigned char',
    't': 'unsigned short',
    'j': 'unsigned int',
    'm': 'unsigned long',
    'x': 'long long',
    'y': 'unsigned long long',
    'z': '...'
  };
  var subs = [];
  var first = true;
  function dump(x) {
    //return;
    if (x) Module.print(x);
    Module.print(func);
    var pre = '';
    for (var a = 0; a < i; a++) pre += ' ';
    Module.print (pre + '^');
  }
  function parseNested() {
    i++;
    if (func[i] === 'K') i++; // ignore const
    var parts = [];
    while (func[i] !== 'E') {
      if (func[i] === 'S') { // substitution
        i++;
        var next = func.indexOf('_', i);
        var num = func.substring(i, next) || 0;
        parts.push(subs[num] || '?');
        i = next+1;
        continue;
      }
      if (func[i] === 'C') { // constructor
        parts.push(parts[parts.length-1]);
        i += 2;
        continue;
      }
      var size = parseInt(func.substr(i));
      var pre = size.toString().length;
      if (!size || !pre) { i--; break; } // counter i++ below us
      var curr = func.substr(i + pre, size);
      parts.push(curr);
      subs.push(curr);
      i += pre + size;
    }
    i++; // skip E
    return parts;
  }
  function parse(rawList, limit, allowVoid) { // main parser
    limit = limit || Infinity;
    var ret = '', list = [];
    function flushList() {
      return '(' + list.join(', ') + ')';
    }
    var name;
    if (func[i] === 'N') {
      // namespaced N-E
      name = parseNested().join('::');
      limit--;
      if (limit === 0) return rawList ? [name] : name;
    } else {
      // not namespaced
      if (func[i] === 'K' || (first && func[i] === 'L')) i++; // ignore const and first 'L'
      var size = parseInt(func.substr(i));
      if (size) {
        var pre = size.toString().length;
        name = func.substr(i + pre, size);
        i += pre + size;
      }
    }
    first = false;
    if (func[i] === 'I') {
      i++;
      var iList = parse(true);
      var iRet = parse(true, 1, true);
      ret += iRet[0] + ' ' + name + '<' + iList.join(', ') + '>';
    } else {
      ret = name;
    }
    paramLoop: while (i < func.length && limit-- > 0) {
      //dump('paramLoop');
      var c = func[i++];
      if (c in basicTypes) {
        list.push(basicTypes[c]);
      } else {
        switch (c) {
          case 'P': list.push(parse(true, 1, true)[0] + '*'); break; // pointer
          case 'R': list.push(parse(true, 1, true)[0] + '&'); break; // reference
          case 'L': { // literal
            i++; // skip basic type
            var end = func.indexOf('E', i);
            var size = end - i;
            list.push(func.substr(i, size));
            i += size + 2; // size + 'EE'
            break;
          }
          case 'A': { // array
            var size = parseInt(func.substr(i));
            i += size.toString().length;
            if (func[i] !== '_') throw '?';
            i++; // skip _
            list.push(parse(true, 1, true)[0] + ' [' + size + ']');
            break;
          }
          case 'E': break paramLoop;
          default: ret += '?' + c; break paramLoop;
        }
      }
    }
    if (!allowVoid && list.length === 1 && list[0] === 'void') list = []; // avoid (void)
    if (rawList) {
      if (ret) {
        list.push(ret + '?');
      }
      return list;
    } else {
      return ret + flushList();
    }
  }
  try {
    // Special-case the entry point, since its name differs from other name mangling.
    if (func == 'Object._main' || func == '_main') {
      return 'main()';
    }
    if (typeof func === 'number') func = Pointer_stringify(func);
    if (func[0] !== '_') return func;
    if (func[1] !== '_') return func; // C function
    if (func[2] !== 'Z') return func;
    switch (func[3]) {
      case 'n': return 'operator new()';
      case 'd': return 'operator delete()';
    }
    return parse();
  } catch(e) {
    return func;
  }
}

function demangleAll(text) {
  return text.replace(/__Z[\w\d_]+/g, function(x) { var y = demangle(x); return x === y ? x : (x + ' [' + y + ']') });
}

function stackTrace() {
  var stack = new Error().stack;
  return stack ? demangleAll(stack) : '(no stack trace available)'; // Stack trace is not available at least on IE10 and Safari 6.
}

// Memory management

var PAGE_SIZE = 4096;
function alignMemoryPage(x) {
  return (x+4095)&-4096;
}

var HEAP;
var HEAP8, HEAPU8, HEAP16, HEAPU16, HEAP32, HEAPU32, HEAPF32, HEAPF64;

var STATIC_BASE = 0, STATICTOP = 0, staticSealed = false; // static area
var STACK_BASE = 0, STACKTOP = 0, STACK_MAX = 0; // stack area
var DYNAMIC_BASE = 0, DYNAMICTOP = 0; // dynamic area handled by sbrk

function enlargeMemory() {
  abort('Cannot enlarge memory arrays. Either (1) compile with -s TOTAL_MEMORY=X with X higher than the current value ' + TOTAL_MEMORY + ', (2) compile with ALLOW_MEMORY_GROWTH which adjusts the size at runtime but prevents some optimizations, or (3) set Module.TOTAL_MEMORY before the program runs.');
}

var TOTAL_STACK = Module['TOTAL_STACK'] || 5242880;
var TOTAL_MEMORY = Module['TOTAL_MEMORY'] || 16777216;
var FAST_MEMORY = Module['FAST_MEMORY'] || 2097152;

var totalMemory = 4096;
while (totalMemory < TOTAL_MEMORY || totalMemory < 2*TOTAL_STACK) {
  if (totalMemory < 16*1024*1024) {
    totalMemory *= 2;
  } else {
    totalMemory += 16*1024*1024
  }
}
if (totalMemory !== TOTAL_MEMORY) {
  Module.printErr('increasing TOTAL_MEMORY to ' + totalMemory + ' to be more reasonable');
  TOTAL_MEMORY = totalMemory;
}

// Initialize the runtime's memory
// check for full engine support (use string 'subarray' to avoid closure compiler confusion)
assert(typeof Int32Array !== 'undefined' && typeof Float64Array !== 'undefined' && !!(new Int32Array(1)['subarray']) && !!(new Int32Array(1)['set']),
       'JS engine does not provide full typed array support');

var buffer = new ArrayBuffer(TOTAL_MEMORY);
HEAP8 = new Int8Array(buffer);
HEAP16 = new Int16Array(buffer);
HEAP32 = new Int32Array(buffer);
HEAPU8 = new Uint8Array(buffer);
HEAPU16 = new Uint16Array(buffer);
HEAPU32 = new Uint32Array(buffer);
HEAPF32 = new Float32Array(buffer);
HEAPF64 = new Float64Array(buffer);

// Endianness check (note: assumes compiler arch was little-endian)
HEAP32[0] = 255;
assert(HEAPU8[0] === 255 && HEAPU8[3] === 0, 'Typed arrays 2 must be run on a little-endian system');

Module['HEAP'] = HEAP;
Module['HEAP8'] = HEAP8;
Module['HEAP16'] = HEAP16;
Module['HEAP32'] = HEAP32;
Module['HEAPU8'] = HEAPU8;
Module['HEAPU16'] = HEAPU16;
Module['HEAPU32'] = HEAPU32;
Module['HEAPF32'] = HEAPF32;
Module['HEAPF64'] = HEAPF64;

function callRuntimeCallbacks(callbacks) {
  while(callbacks.length > 0) {
    var callback = callbacks.shift();
    if (typeof callback == 'function') {
      callback();
      continue;
    }
    var func = callback.func;
    if (typeof func === 'number') {
      if (callback.arg === undefined) {
        Runtime.dynCall('v', func);
      } else {
        Runtime.dynCall('vi', func, [callback.arg]);
      }
    } else {
      func(callback.arg === undefined ? null : callback.arg);
    }
  }
}

var __ATPRERUN__  = []; // functions called before the runtime is initialized
var __ATINIT__    = []; // functions called during startup
var __ATMAIN__    = []; // functions called when main() is to be run
var __ATEXIT__    = []; // functions called during shutdown
var __ATPOSTRUN__ = []; // functions called after the runtime has exited

var runtimeInitialized = false;

function preRun() {
  // compatibility - merge in anything from Module['preRun'] at this time
  if (Module['preRun']) {
    if (typeof Module['preRun'] == 'function') Module['preRun'] = [Module['preRun']];
    while (Module['preRun'].length) {
      addOnPreRun(Module['preRun'].shift());
    }
  }
  callRuntimeCallbacks(__ATPRERUN__);
}

function ensureInitRuntime() {
  if (runtimeInitialized) return;
  runtimeInitialized = true;
  callRuntimeCallbacks(__ATINIT__);
}

function preMain() {
  callRuntimeCallbacks(__ATMAIN__);
}

function exitRuntime() {
  callRuntimeCallbacks(__ATEXIT__);
}

function postRun() {
  // compatibility - merge in anything from Module['postRun'] at this time
  if (Module['postRun']) {
    if (typeof Module['postRun'] == 'function') Module['postRun'] = [Module['postRun']];
    while (Module['postRun'].length) {
      addOnPostRun(Module['postRun'].shift());
    }
  }
  callRuntimeCallbacks(__ATPOSTRUN__);
}

function addOnPreRun(cb) {
  __ATPRERUN__.unshift(cb);
}
Module['addOnPreRun'] = Module.addOnPreRun = addOnPreRun;

function addOnInit(cb) {
  __ATINIT__.unshift(cb);
}
Module['addOnInit'] = Module.addOnInit = addOnInit;

function addOnPreMain(cb) {
  __ATMAIN__.unshift(cb);
}
Module['addOnPreMain'] = Module.addOnPreMain = addOnPreMain;

function addOnExit(cb) {
  __ATEXIT__.unshift(cb);
}
Module['addOnExit'] = Module.addOnExit = addOnExit;

function addOnPostRun(cb) {
  __ATPOSTRUN__.unshift(cb);
}
Module['addOnPostRun'] = Module.addOnPostRun = addOnPostRun;

// Tools

// This processes a JS string into a C-line array of numbers, 0-terminated.
// For LLVM-originating strings, see parser.js:parseLLVMString function
function intArrayFromString(stringy, dontAddNull, length /* optional */) {
  var ret = (new Runtime.UTF8Processor()).processJSString(stringy);
  if (length) {
    ret.length = length;
  }
  if (!dontAddNull) {
    ret.push(0);
  }
  return ret;
}
Module['intArrayFromString'] = intArrayFromString;

function intArrayToString(array) {
  var ret = [];
  for (var i = 0; i < array.length; i++) {
    var chr = array[i];
    if (chr > 0xFF) {
      chr &= 0xFF;
    }
    ret.push(String.fromCharCode(chr));
  }
  return ret.join('');
}
Module['intArrayToString'] = intArrayToString;

// Write a Javascript array to somewhere in the heap
function writeStringToMemory(string, buffer, dontAddNull) {
  var array = intArrayFromString(string, dontAddNull);
  var i = 0;
  while (i < array.length) {
    var chr = array[i];
    HEAP8[(((buffer)+(i))|0)]=chr;
    i = i + 1;
  }
}
Module['writeStringToMemory'] = writeStringToMemory;

function writeArrayToMemory(array, buffer) {
  for (var i = 0; i < array.length; i++) {
    HEAP8[(((buffer)+(i))|0)]=array[i];
  }
}
Module['writeArrayToMemory'] = writeArrayToMemory;

function writeAsciiToMemory(str, buffer, dontAddNull) {
  for (var i = 0; i < str.length; i++) {
    HEAP8[(((buffer)+(i))|0)]=str.charCodeAt(i);
  }
  if (!dontAddNull) HEAP8[(((buffer)+(str.length))|0)]=0;
}
Module['writeAsciiToMemory'] = writeAsciiToMemory;

function unSign(value, bits, ignore) {
  if (value >= 0) {
    return value;
  }
  return bits <= 32 ? 2*Math.abs(1 << (bits-1)) + value // Need some trickery, since if bits == 32, we are right at the limit of the bits JS uses in bitshifts
                    : Math.pow(2, bits)         + value;
}
function reSign(value, bits, ignore) {
  if (value <= 0) {
    return value;
  }
  var half = bits <= 32 ? Math.abs(1 << (bits-1)) // abs is needed if bits == 32
                        : Math.pow(2, bits-1);
  if (value >= half && (bits <= 32 || value > half)) { // for huge values, we can hit the precision limit and always get true here. so don't do that
                                                       // but, in general there is no perfect solution here. With 64-bit ints, we get rounding and errors
                                                       // TODO: In i64 mode 1, resign the two parts separately and safely
    value = -2*half + value; // Cannot bitshift half, as it may be at the limit of the bits JS uses in bitshifts
  }
  return value;
}

// check for imul support, and also for correctness ( https://bugs.webkit.org/show_bug.cgi?id=126345 )
if (!Math['imul'] || Math['imul'](0xffffffff, 5) !== -5) Math['imul'] = function imul(a, b) {
  var ah  = a >>> 16;
  var al = a & 0xffff;
  var bh  = b >>> 16;
  var bl = b & 0xffff;
  return (al*bl + ((ah*bl + al*bh) << 16))|0;
};
Math.imul = Math['imul'];


var Math_abs = Math.abs;
var Math_cos = Math.cos;
var Math_sin = Math.sin;
var Math_tan = Math.tan;
var Math_acos = Math.acos;
var Math_asin = Math.asin;
var Math_atan = Math.atan;
var Math_atan2 = Math.atan2;
var Math_exp = Math.exp;
var Math_log = Math.log;
var Math_sqrt = Math.sqrt;
var Math_ceil = Math.ceil;
var Math_floor = Math.floor;
var Math_pow = Math.pow;
var Math_imul = Math.imul;
var Math_fround = Math.fround;
var Math_min = Math.min;

// A counter of dependencies for calling run(). If we need to
// do asynchronous work before running, increment this and
// decrement it. Incrementing must happen in a place like
// PRE_RUN_ADDITIONS (used by emcc to add file preloading).
// Note that you can add dependencies in preRun, even though
// it happens right before run - run will be postponed until
// the dependencies are met.
var runDependencies = 0;
var runDependencyWatcher = null;
var dependenciesFulfilled = null; // overridden to take different actions when all run dependencies are fulfilled

function addRunDependency(id) {
  runDependencies++;
  if (Module['monitorRunDependencies']) {
    Module['monitorRunDependencies'](runDependencies);
  }
}
Module['addRunDependency'] = addRunDependency;
function removeRunDependency(id) {
  runDependencies--;
  if (Module['monitorRunDependencies']) {
    Module['monitorRunDependencies'](runDependencies);
  }
  if (runDependencies == 0) {
    if (runDependencyWatcher !== null) {
      clearInterval(runDependencyWatcher);
      runDependencyWatcher = null;
    }
    if (dependenciesFulfilled) {
      var callback = dependenciesFulfilled;
      dependenciesFulfilled = null;
      callback(); // can add another dependenciesFulfilled
    }
  }
}
Module['removeRunDependency'] = removeRunDependency;

Module["preloadedImages"] = {}; // maps url to image data
Module["preloadedAudios"] = {}; // maps url to audio data


var memoryInitializer = null;

// === Body ===
var __ZTVN10__cxxabiv117__class_type_infoE = 16344;
var __ZTVN10__cxxabiv120__si_class_type_infoE = 16384;




STATIC_BASE = 8;

STATICTOP = STATIC_BASE + Runtime.alignMemory(17211);
/* global initializers */ __ATINIT__.push({ func: function() { __GLOBAL__I_a() } });


/* memory initializer */ allocate([123,114,101,116,117,114,110,32,105,110,112,117,116,83,105,122,101,59,125,0,0,0,0,0,123,114,101,116,117,114,110,32,105,110,112,117,116,73,46,108,101,110,103,116,104,59,125,0,123,114,101,116,117,114,110,32,105,110,112,117,116,73,91,36,48,93,59,125,0,0,0,0,123,114,101,116,117,114,110,32,105,110,112,117,116,74,91,36,48,93,59,125,0,0,0,0,123,114,101,116,117,114,110,32,105,110,112,117,116,86,91,36,48,93,59,125,0,0,0,0,123,32,111,117,116,112,117,116,73,91,36,48,93,32,61,32,36,49,59,32,111,117,116,112,117,116,74,91,36,48,93,32,61,32,36,50,59,32,111,117,116,112,117,116,86,91,36,48,93,32,61,32,36,51,59,32,125,0,0,0,0,0,0,0,123,110,111,114,109,32,61,32,36,48,59,125,0,0,0,0,114,111,119,115,40,41,62,48,32,38,38,32,99,111,108,115,40,41,62,48,32,38,38,32,34,121,111,117,32,97,114,101,32,117,115,105,110,103,32,97,32,110,111,110,32,105,110,105,116,105,97,108,105,122,101,100,32,109,97,116,114,105,120,34,0,0,0,0,0,0,0,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,83,112,97,114,115,101,67,111,114,101,47,83,112,97,114,115,101,82,101,100,117,120,46,104,0,0,0,0,0,115,117,109,0,0,0,0,0,97,76,104,115,46,114,111,119,115,40,41,32,61,61,32,97,82,104,115,46,114,111,119,115,40,41,32,38,38,32,97,76,104,115,46,99,111,108,115,40,41,32,61,61,32,97,82,104,115,46,99,111,108,115,40,41,0,0,0,0,0,0,0,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,67,111,114,101,47,67,119,105,115,101,66,105,110,97,114,121,79,112,46,104,0,67,119,105,115,101,66,105,110,97,114,121,79,112,0,0,0,108,104,115,46,111,117,116,101,114,83,105,122,101,40,41,32,61,61,32,114,104,115,46,105,110,110,101,114,83,105,122,101,40,41,0,0,0,0,0,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,83,112,97,114,115,101,67,111,114,101,47,83,112,97,114,115,101,83,112,97,114,115,101,80,114,111,100,117,99,116,87,105,116,104,80,114,117,110,105,110,103,46,104,0,0,115,112,97,114,115,101,95,115,112,97,114,115,101,95,112,114,111,100,117,99,116,95,119,105,116,104,95,112,114,117,110,105,110,103,95,105,109,112,108,0,115,105,122,101,95,116,40,109,95,111,117,116,101,114,73,110,100,101,120,91,111,117,116,101,114,43,49,93,41,32,61,61,32,109,95,100,97,116,97,46,115,105,122,101,40,41,32,38,38,32,34,73,110,118,97,108,105,100,32,111,114,100,101,114,101,100,32,105,110,115,101,114,116,105,111,110,32,40,105,110,118,97,108,105,100,32,111,117,116,101,114,32,105,110,100,101,120,41,34,0,0,0,0,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,83,112,97,114,115,101,67,111,114,101,47,83,112,97,114,115,101,77,97,116,114,105,120,46,104,0,0,0,0,105,110,115,101,114,116,66,97,99,107,66,121,79,117,116,101,114,73,110,110,101,114,0,0,40,109,95,111,117,116,101,114,73,110,100,101,120,91,111,117,116,101,114,43,49,93,45,109,95,111,117,116,101,114,73,110,100,101,120,91,111,117,116,101,114,93,61,61,48,32,124,124,32,109,95,100,97,116,97,46,105,110,100,101,120,40,109,95,100,97,116,97,46,115,105,122,101,40,41,45,49,41,60,105,110,110,101,114,41,32,38,38,32,34,73,110,118,97,108,105,100,32,111,114,100,101,114,101,100,32,105,110,115,101,114,116,105,111,110,32,40,105,110,118,97,108,105,100,32,105,110,110,101,114,32,105,110,100,101,120,41,34,0,0,0,0,0,0,109,95,111,117,116,101,114,73,110,100,101,120,91,111,117,116,101,114,93,61,61,73,110,100,101,120,40,109,95,100,97,116,97,46,115,105,122,101,40,41,41,32,38,38,32,34,89,111,117,32,109,117,115,116,32,99,97,108,108,32,115,116,97,114,116,86,101,99,32,102,111,114,32,101,97,99,104,32,105,110,110,101,114,32,118,101,99,116,111,114,32,115,101,113,117,101,110,116,105,97,108,108,121,34,0,0,0,0,0,0,0,0,115,116,97,114,116,86,101,99,0,0,0,0,0,0,0,0,109,95,111,117,116,101,114,73,110,100,101,120,91,111,117,116,101,114,43,49,93,61,61,48,32,38,38,32,34,89,111,117,32,109,117,115,116,32,99,97,108,108,32,115,116,97,114,116,86,101,99,32,102,111,114,32,101,97,99,104,32,105,110,110,101,114,32,118,101,99,116,111,114,32,115,101,113,117,101,110,116,105,97,108,108,121,34,0,109,95,109,111,100,101,61,61,73,115,83,112,97,114,115,101,0,0,0,0,0,0,0,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,83,112,97,114,115,101,67,111,114,101,47,65,109,98,105,86,101,99,116,111,114,46,104,0,0,0,0,0,0,99,111,101,102,102,82,101,102,0,0,0,0,0,0,0,0,105,62,61,108,108,69,108,101,109,101,110,116,115,91,109,95,108,108,67,117,114,114,101,110,116,93,46,105,110,100,101,120,32,38,38,32,34,121,111,117,32,109,117,115,116,32,99,97,108,108,32,114,101,115,116,97,114,116,40,41,32,98,101,102,111,114,101,32,105,110,115,101,114,116,105,110,103,32,97,110,32,101,108,101,109,101,110,116,32,119,105,116,104,32,108,111,119,101,114,32,111,114,32,101,113,117,97,108,32,105,110,100,101,120,34,0,0,0,0,0,105,115,67,111,109,112,114,101,115,115,101,100,40,41,32,38,38,32,34,84,104,105,115,32,102,117,110,99,116,105,111,110,32,100,111,101,115,32,110,111,116,32,109,97,107,101,32,115,101,110,115,101,32,105,110,32,110,111,110,32,99,111,109,112,114,101,115,115,101,100,32,109,111,100,101,46,34,0,0,0,114,101,115,101,114,118,101,0,109,97,116,46,114,111,119,115,40,41,62,48,32,38,38,32,109,97,116,46,99,111,108,115,40,41,62,48,32,38,38,32,34,121,111,117,32,97,114,101,32,117,115,105,110,103,32,97,110,32,101,109,112,116,121,32,109,97,116,114,105,120,34,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,67,111,114,101,47,82,101,100,117,120,46,104,0,114,117,110,0,0,0,0,0,118,101,99,83,105,122,101,32,62,61,32,48,0,0,0,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,67,111,114,101,47,77,97,112,66,97,115,101,46,104,0,0,0,0,0,0,0,77,97,112,66,97,115,101,0,105,110,100,101,120,32,62,61,32,48,32,38,38,32,105,110,100,101,120,32,60,32,115,105,122,101,40,41,0,0,0,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,67,111,114,101,47,68,101,110,115,101,67,111,101,102,102,115,66,97,115,101,46,104,0,0,0,0,0,0,0,111,112,101,114,97,116,111,114,91,93,0,0,0,0,0,0,100,105,109,32,62,61,32,48,0,0,0,0,0,0,0,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,67,111,114,101,47,77,97,116,114,105,120,46,104,0,0,0,0,0,0,0,0,77,97,116,114,105,120,0,0,114,111,119,115,40,41,32,61,61,32,111,116,104,101,114,46,114,111,119,115,40,41,32,38,38,32,99,111,108,115,40,41,32,61,61,32,111,116,104,101,114,46,99,111,108,115,40,41,0,0,0,0,0,0,0,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,67,111,114,101,47,65,115,115,105,103,110,46,104,0,0,0,0,0,0,0,0,108,97,122,121,65,115,115,105,103,110,0,0,0,0,0,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,83,112,97,114,115,101,67,111,114,101,47,67,111,110,115,101,114,118,97,116,105,118,101,83,112,97,114,115,101,83,112,97,114,115,101,80,114,111,100,117,99,116,46,104,0,99,111,110,115,101,114,118,97,116,105,118,101,95,115,112,97,114,115,101,95,115,112,97,114,115,101,95,112,114,111,100,117,99,116,95,105,109,112,108,0,109,95,108,104,115,46,99,111,108,115,40,41,32,61,61,32,109,95,114,104,115,46,114,111,119,115,40,41,0,0,0,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,83,112,97,114,115,101,67,111,114,101,47,83,112,97,114,115,101,80,114,111,100,117,99,116,46,104,0,0,0,105,110,105,116,0,0,0,0,66,97,115,101,58,58,109,95,102,97,99,116,111,114,105,122,97,116,105,111,110,73,115,79,107,32,38,38,32,34,83,105,109,112,108,105,99,105,97,108,32,76,76,84,32,110,111,116,32,102,97,99,116,111,114,105,122,101,100,34,0,0,0,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,83,112,97,114,115,101,67,104,111,108,101,115,107,121,47,83,105,109,112,108,105,99,105,97,108,67,104,111,108,101,115,107,121,46,104,0,0,109,97,116,114,105,120,76,0,109,97,116,114,105,120,46,114,111,119,115,40,41,61,61,109,97,116,114,105,120,46,99,111,108,115,40,41,0,0,0,0,99,111,109,112,117,116,101,0,109,95,97,110,97,108,121,115,105,115,73,115,79,107,32,38,38,32,34,89,111,117,32,109,117,115,116,32,102,105,114,115,116,32,99,97,108,108,32,97,110,97,108,121,122,101,80,97,116,116,101,114,110,40,41,34,0,0,0,0,0,0,0,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,83,112,97,114,115,101,67,104,111,108,101,115,107,121,47,83,105,109,112,108,105,99,105,97,108,67,104,111,108,101,115,107,121,95,105,109,112,108,46,104,0,0,0,0,0,102,97,99,116,111,114,105,122,101,95,112,114,101,111,114,100,101,114,101,100,0,0,0,0,97,112,46,114,111,119,115,40,41,61,61,97,112,46,99,111,108,115,40,41,0,0,0,0,109,95,112,97,114,101,110,116,46,115,105,122,101,40,41,61,61,115,105,122,101,0,0,0,109,95,110,111,110,90,101,114,111,115,80,101,114,67,111,108,46,115,105,122,101,40,41,61,61,115,105,122,101,0,0,0,40,40,83,105,122,101,65,116,67,111,109,112,105,108,101,84,105,109,101,32,61,61,32,68,121,110,97,109,105,99,32,38,38,32,40,77,97,120,83,105,122,101,65,116,67,111,109,112,105,108,101,84,105,109,101,61,61,68,121,110,97,109,105,99,32,124,124,32,115,105,122,101,60,61,77,97,120,83,105,122,101,65,116,67,111,109,112,105,108,101,84,105,109,101,41,41,32,124,124,32,83,105,122,101,65,116,67,111,109,112,105,108,101,84,105,109,101,32,61,61,32,115,105,122,101,41,32,38,38,32,115,105,122,101,62,61,48,0,0,0,0,0,0,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,67,111,114,101,47,80,108,97,105,110,79,98,106,101,99,116,66,97,115,101,46,104,0,0,0,0,0,0,0,114,101,115,105,122,101,0,0,97,46,114,111,119,115,40,41,61,61,97,46,99,111,108,115,40,41,0,0,0,0,0,0,111,114,100,101,114,105,110,103,0,0,0,0,0,0,0,0,40,33,40,82,111,119,115,65,116,67,111,109,112,105,108,101,84,105,109,101,33,61,68,121,110,97,109,105,99,41,32,124,124,32,40,110,98,82,111,119,115,61,61,82,111,119,115,65,116,67,111,109,112,105,108,101,84,105,109,101,41,41,32,38,38,32,40,33,40,67,111,108,115,65,116,67,111,109,112,105,108,101,84,105,109,101,33,61,68,121,110,97,109,105,99,41,32,124,124,32,40,110,98,67,111,108,115,61,61,67,111,108,115,65,116,67,111,109,112,105,108,101,84,105,109,101,41,41,32,38,38,32,40,33,40,82,111,119,115,65,116,67,111,109,112,105,108,101,84,105,109,101,61,61,68,121,110,97,109,105,99,32,38,38,32,77,97,120,82,111,119,115,65,116,67,111,109,112,105,108,101,84,105,109,101,33,61,68,121,110,97,109,105,99,41,32,124,124,32,40,110,98,82,111,119,115,60,61,77,97,120,82,111,119,115,65,116,67,111,109,112,105,108,101,84,105,109,101,41,41,32,38,38,32,40,33,40,67,111,108,115,65,116,67,111,109,112,105,108,101,84,105,109,101,61,61,68,121,110,97,109,105,99,32,38,38,32,77,97,120,67,111,108,115,65,116,67,111,109,112,105,108,101,84,105,109,101,33,61,68,121,110,97,109,105,99,41,32,124,124,32,40,110,98,67,111,108,115,60,61,77,97,120,67,111,108,115,65,116,67,111,109,112,105,108,101,84,105,109,101,41,41,32,38,38,32,110,98,82,111,119,115,62,61,48,32,38,38,32,110,98,67,111,108,115,62,61,48,32,38,38,32,34,73,110,118,97,108,105,100,32,115,105,122,101,115,32,119,104,101,110,32,114,101,115,105,122,105,110,103,32,97,32,109,97,116,114,105,120,32,111,114,32,97,114,114,97,121,46,34,0,0,0,0,0,0,114,111,119,115,40,41,61,61,99,111,108,115,40,41,32,38,38,32,34,83,101,108,102,65,100,106,111,105,110,116,86,105,101,119,32,105,115,32,111,110,108,121,32,102,111,114,32,115,113,117,97,114,101,100,32,109,97,116,114,105,99,101,115,34,0,0,0,0,0,0,0,0,46,47,46,46,47,108,105,98,47,69,105,103,101,110,47,115,114,99,47,83,112,97,114,115,101,67,111,114,101,47,83,112,97,114,115,101,83,101,108,102,65,100,106,111,105,110,116,86,105,101,119,46,104,0,0,0,83,112,97,114,115,101,83,101,108,102,65,100,106,111,105,110,116,86,105,101,119,0,0,0,105,116,45,62,114,111,119,40,41,62,61,48,32,38,38,32,105,116,45,62,114,111,119,40,41,60,109,97,116,46,114,111,119,115,40,41,32,38,38,32,105,116,45,62,99,111,108,40,41,62,61,48,32,38,38,32,105,116,45,62,99,111,108,40,41,60,109,97,116,46,99,111,108,115,40,41,0,0,0,0,115,101,116,95,102,114,111,109,95,116,114,105,112,108,101,116,115,0,0,0,0,0,0,0,33,105,115,67,111,109,112,114,101,115,115,101,100,40,41,0,115,117,109,117,112,68,117,112,108,105,99,97,116,101,115,0,105,110,115,101,114,116,66,97,99,107,85,110,99,111,109,112,114,101,115,115,101,100,0,0,109,95,105,110,110,101,114,78,111,110,90,101,114,111,115,91,111,117,116,101,114,93,60,61,40,109,95,111,117,116,101,114,73,110,100,101,120,91,111,117,116,101,114,43,49,93,32,45,32,109,95,111,117,116,101,114,73,110,100,101,120,91,111,117,116,101,114,93,41,0,0,0,111,112,101,114,97,116,111,114,40,41], "i8", ALLOC_NONE, Runtime.GLOBAL_BASE);
/* memory initializer */ allocate([72,18,0,0,1,0,0,0,2,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,2,0,0,0,2,0,0,0,3,0,0,0,4,0,0,0,1,0,0,0,3,0,0,0,2,0,0,0,78,83,116,51,95,95,49,49,49,95,95,115,116,100,111,117,116,98,117,102,73,119,69,69,0,0,0,0,0,0,0,0,8,64,0,0,40,18,0,0,120,24,0,0,0,0,0,0,0,0,0,0,176,18,0,0,3,0,0,0,4,0,0,0,2,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,5,0,0,0,2,0,0,0,2,0,0,0,6,0,0,0,7,0,0,0,3,0,0,0,4,0,0,0,4,0,0,0,78,83,116,51,95,95,49,49,48,95,95,115,116,100,105,110,98,117,102,73,119,69,69,0,8,64,0,0,152,18,0,0,120,24,0,0,0,0,0,0,117,110,115,117,112,112,111,114,116,101,100,32,108,111,99,97,108,101,32,102,111,114,32,115,116,97,110,100,97,114,100,32,105,110,112,117,116,0,0,0,0,0,0,0,72,19,0,0,5,0,0,0,6,0,0,0,3,0,0,0,5,0,0,0,2,0,0,0,2,0,0,0,8,0,0,0,9,0,0,0,6,0,0,0,10,0,0,0,11,0,0,0,5,0,0,0,7,0,0,0,6,0,0,0,78,83,116,51,95,95,49,49,49,95,95,115,116,100,111,117,116,98,117,102,73,99,69,69,0,0,0,0,0,0,0,0,8,64,0,0,40,19,0,0,56,24,0,0,0,0,0,0,0,0,0,0,176,19,0,0,7,0,0,0,8,0,0,0,4,0,0,0,5,0,0,0,2,0,0,0,2,0,0,0,12,0,0,0,9,0,0,0,6,0,0,0,13,0,0,0,14,0,0,0,7,0,0,0,8,0,0,0,8,0,0,0,78,83,116,51,95,95,49,49,48,95,95,115,116,100,105,110,98,117,102,73,99,69,69,0,8,64,0,0,152,19,0,0,56,24,0,0,0,0,0,0,78,83,116,51,95,95,49,49,52,95,95,115,104,97,114,101,100,95,99,111,117,110,116,69,0,0,0,0,0,0,0,0,224,63,0,0,192,19,0,0,0,0,0,0,40,20,0,0,9,0,0,0,10,0,0,0,15,0,0,0,0,0,0,0,0,0,0,0,144,20,0,0,11,0,0,0,12,0,0,0,16,0,0,0,0,0,0,0,83,116,49,49,108,111,103,105,99,95,101,114,114,111,114,0,8,64,0,0,24,20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,104,20,0,0,9,0,0,0,13,0,0,0,15,0,0,0,0,0,0,0,83,116,49,50,108,101,110,103,116,104,95,101,114,114,111,114,0,0,0,0,0,0,0,0,8,64,0,0,80,20,0,0,40,20,0,0,0,0,0,0,83,116,49,51,114,117,110,116,105,109,101,95,101,114,114,111,114,0,0,0,0,0,0,0,8,64,0,0,120,20,0,0,0,0,0,0,0,0,0,0,58,32,0,0,0,0,0,0,0,0,0,0,216,20,0,0,14,0,0,0,15,0,0,0,16,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,50,115,121,115,116,101,109,95,101,114,114,111,114,69,0,0,8,64,0,0,192,20,0,0,144,20,0,0,0,0,0,0,78,83,116,51,95,95,49,49,52,101,114,114,111,114,95,99,97,116,101,103,111,114,121,69,0,0,0,0,0,0,0,0,224,63,0,0,232,20,0,0,78,83,116,51,95,95,49,49,50,95,95,100,111,95,109,101,115,115,97,103,101,69,0,0,8,64,0,0,16,21,0,0,8,21,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,98,97,115,105,99,95,115,116,114,105,110,103,0,0,0,0,0,0,0,0,56,24,0,0,16,0,0,0,17,0,0,0,5,0,0,0,5,0,0,0,2,0,0,0,2,0,0,0,12,0,0,0,9,0,0,0,6,0,0,0,10,0,0,0,11,0,0,0,5,0,0,0,8,0,0,0,8,0,0,0,0,0,0,0,120,24,0,0,18,0,0,0,19,0,0,0,6,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,5,0,0,0,2,0,0,0,2,0,0,0,3,0,0,0,4,0,0,0,1,0,0,0,4,0,0,0,4,0,0,0,8,0,0,0,0,0,0,0,176,24,0,0,20,0,0,0,21,0,0,0,248,255,255,255,248,255,255,255,176,24,0,0,22,0,0,0,23,0,0,0,8,0,0,0,0,0,0,0,248,24,0,0,24,0,0,0,25,0,0,0,248,255,255,255,248,255,255,255,248,24,0,0,26,0,0,0,27,0,0,0,4,0,0,0,0,0,0,0,64,25,0,0,28,0,0,0,29,0,0,0,252,255,255,255,252,255,255,255,64,25,0,0,30,0,0,0,31,0,0,0,4,0,0,0,0,0,0,0,136,25,0,0,32,0,0,0,33,0,0,0,252,255,255,255,252,255,255,255,136,25,0,0,34,0,0,0,35,0,0,0,105,111,115,116,114,101,97,109,0,0,0,0,0,0,0,0,117,110,115,112,101,99,105,102,105,101,100,32,105,111,115,116,114,101,97,109,95,99,97,116,101,103,111,114,121,32,101,114,114,111,114,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,23,0,0,36,0,0,0,37,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,120,23,0,0,38,0,0,0,39,0,0,0,105,111,115,95,98,97,115,101,58,58,99,108,101,97,114,0,78,83,116,51,95,95,49,56,105,111,115,95,98,97,115,101,55,102,97,105,108,117,114,101,69,0,0,0,0,0,0,0,8,64,0,0,48,23,0,0,216,20,0,0,0,0,0,0,78,83,116,51,95,95,49,56,105,111,115,95,98,97,115,101,69,0,0,0,0,0,0,0,224,63,0,0,96,23,0,0,78,83,116,51,95,95,49,57,98,97,115,105,99,95,105,111,115,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,0,0,0,0,0,0,0,8,64,0,0,128,23,0,0,120,23,0,0,0,0,0,0,78,83,116,51,95,95,49,57,98,97,115,105,99,95,105,111,115,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,0,0,0,0,0,0,0,8,64,0,0,192,23,0,0,120,23,0,0,0,0,0,0,78,83,116,51,95,95,49,49,53,98,97,115,105,99,95,115,116,114,101,97,109,98,117,102,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,0,0,0,0,0,0,0,0,224,63,0,0,0,24,0,0,78,83,116,51,95,95,49,49,53,98,97,115,105,99,95,115,116,114,101,97,109,98,117,102,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,0,0,0,0,0,0,0,0,224,63,0,0,64,24,0,0,78,83,116,51,95,95,49,49,51,98,97,115,105,99,95,105,115,116,114,101,97,109,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,0,0,104,64,0,0,128,24,0,0,0,0,0,0,1,0,0,0,176,23,0,0,3,244,255,255,78,83,116,51,95,95,49,49,51,98,97,115,105,99,95,105,115,116,114,101,97,109,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,0,0,104,64,0,0,200,24,0,0,0,0,0,0,1,0,0,0,240,23,0,0,3,244,255,255,78,83,116,51,95,95,49,49,51,98,97,115,105,99,95,111,115,116,114,101,97,109,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,0,0,104,64,0,0,16,25,0,0,0,0,0,0,1,0,0,0,176,23,0,0,3,244,255,255,78,83,116,51,95,95,49,49,51,98,97,115,105,99,95,111,115,116,114,101,97,109,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,0,0,104,64,0,0,88,25,0,0,0,0,0,0,1,0,0,0,240,23,0,0,3,244,255,255,0,0,0,0,232,25,0,0,40,0,0,0,41,0,0,0,17,0,0,0,1,0,0,0,9,0,0,0,10,0,0,0,2,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,57,95,95,105,111,115,116,114,101,97,109,95,99,97,116,101,103,111,114,121,69,0,0,0,8,64,0,0,200,25,0,0,40,21,0,0,0,0,0,0,0,0,0,0,16,40,0,0,42,0,0,0,43,0,0,0,44,0,0,0,1,0,0,0,3,0,0,0,11,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,56,40,0,0,45,0,0,0,46,0,0,0,44,0,0,0,2,0,0,0,4,0,0,0,12,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,232,44,0,0,47,0,0,0,48,0,0,0,44,0,0,0,1,0,0,0,2,0,0,0,3,0,0,0,4,0,0,0,5,0,0,0,6,0,0,0,7,0,0,0,8,0,0,0,9,0,0,0,10,0,0,0,11,0,0,0,0,0,0,0,0,0,0,0,48,49,50,51,52,53,54,55,56,57,97,98,99,100,101,102,65,66,67,68,69,70,120,88,43,45,112,80,105,73,110,78,0,0,0,0,0,0,0,0,37,112,0,0,0,0,0,0,0,0,0,0,128,45,0,0,49,0,0,0,50,0,0,0,44,0,0,0,12,0,0,0,13,0,0,0,14,0,0,0,15,0,0,0,16,0,0,0,17,0,0,0,18,0,0,0,19,0,0,0,20,0,0,0,21,0,0,0,22,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,64,46,0,0,51,0,0,0,52,0,0,0,44,0,0,0,3,0,0,0,4,0,0,0,23,0,0,0,5,0,0,0,24,0,0,0,1,0,0,0,2,0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,216,46,0,0,53,0,0,0,54,0,0,0,44,0,0,0,7,0,0,0,8,0,0,0,25,0,0,0,9,0,0,0,26,0,0,0,3,0,0,0,4,0,0,0,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,37,0,0,0,0,0,0,0,37,112,0,0,0,0,0,0,0,0,0,0,0,42,0,0,55,0,0,0,56,0,0,0,44,0,0,0,18,0,0,0,27,0,0,0,28,0,0,0,29,0,0,0,30,0,0,0,31,0,0,0,1,0,0,0,248,255,255,255,0,42,0,0,19,0,0,0,20,0,0,0,21,0,0,0,22,0,0,0,23,0,0,0,24,0,0,0,25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,37,72,58,37,77,58,37,83,37,109,47,37,100,47,37,121,37,89,45,37,109,45,37,100,37,73,58,37,77,58,37,83,32,37,112,0,0,0,0,0,37,72,58,37,77,0,0,0,37,72,58,37,77,58,37,83,0,0,0,0,160,42,0,0,57,0,0,0,58,0,0,0,44,0,0,0,26,0,0,0,32,0,0,0,33,0,0,0,34,0,0,0,35,0,0,0,36,0,0,0,2,0,0,0,248,255,255,255,160,42,0,0,27,0,0,0,28,0,0,0,29,0,0,0,30,0,0,0,31,0,0,0,32,0,0,0,33,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,37,0,0,0,72,0,0,0,58,0,0,0,37,0,0,0,77,0,0,0,58,0,0,0,37,0,0,0,83,0,0,0,37,0,0,0,109,0,0,0,47,0,0,0,37,0,0,0,100,0,0,0,47,0,0,0,37,0,0,0,121,0,0,0,37,0,0,0,89,0,0,0,45,0,0,0,37,0,0,0,109,0,0,0,45,0,0,0,37,0,0,0,100,0,0,0,37,0,0,0,73,0,0,0,58,0,0,0,37,0,0,0,77,0,0,0,58,0,0,0,37,0,0,0,83,0,0,0,32,0,0,0,37,0,0,0,112,0,0,0,0,0,0,0,37,0,0,0,72,0,0,0,58,0,0,0,37,0,0,0,77,0,0,0,0,0,0,0,37,0,0,0,72,0,0,0,58,0,0,0,37,0,0,0,77,0,0,0,58,0,0,0,37,0,0,0,83,0,0,0,0,0,0,0,48,43,0,0,59,0,0,0,60,0,0,0,44,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,152,43,0,0,61,0,0,0,62,0,0,0,44,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,136,40,0,0,63,0,0,0,64,0,0,0,44,0,0,0,34,0,0,0,35,0,0,0,7,0,0,0,8,0,0,0,9,0,0,0,10,0,0,0,36,0,0,0,11,0,0,0,12,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,200,40,0,0,65,0,0,0,66,0,0,0,44,0,0,0,37,0,0,0,38,0,0,0,13,0,0,0,14,0,0,0,15,0,0,0,16,0,0,0,39,0,0,0,17,0,0,0,18,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,41,0,0,67,0,0,0,68,0,0,0,44,0,0,0,40,0,0,0,41,0,0,0,19,0,0,0,20,0,0,0,21,0,0,0,22,0,0,0,42,0,0,0,23,0,0,0,24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,72,41,0,0,69,0,0,0,70,0,0,0,44,0,0,0,43,0,0,0,44,0,0,0,25,0,0,0,26,0,0,0,27,0,0,0,28,0,0,0,45,0,0,0,29,0,0,0,30,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,104,47,0,0,71,0,0,0,72,0,0,0,44,0,0,0,3,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,48,49,50,51,52,53,54,55,56,57,0,0,0,0,0,0,37,76,102,0,0,0,0,0,109,111,110,101,121,95,103,101,116,32,101,114,114,111,114,0,0,0,0,0,248,47,0,0,73,0,0,0,74,0,0,0,44,0,0,0,5,0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,48,49,50,51,52,53,54,55,56,57,0,0,0,0,0,0,0,0,0,0,136,48,0,0,75,0,0,0,76,0,0,0,44,0,0,0,1,0,0,0,37,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,37,46,48,76,102,0,0,0,0,0,0,0,24,49,0,0,77,0,0,0,78,0,0,0,44,0,0,0,2,0,0,0,38,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,240,43,0,0,79,0,0,0,80,0,0,0,44,0,0,0,13,0,0,0,11,0,0,0,31,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,40,44,0,0,81,0,0,0,82,0,0,0,44,0,0,0,14,0,0,0,12,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,118,101,99,116,111,114,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,67,0,0,0,0,0,0,0,0,0,0,0,232,39,0,0,83,0,0,0,84,0,0,0,44,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,192,36,0,0,85,0,0,0,86,0,0,0,44,0,0,0,9,0,0,0,15,0,0,0,10,0,0,0,16,0,0,0,11,0,0,0,1,0,0,0,17,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,160,37,0,0,87,0,0,0,88,0,0,0,44,0,0,0,1,0,0,0,2,0,0,0,4,0,0,0,46,0,0,0,47,0,0,0,5,0,0,0,48,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,152,39,0,0,89,0,0,0,90,0,0,0,44,0,0,0,49,0,0,0,50,0,0,0,33,0,0,0,34,0,0,0,35,0,0,0,0,0,0,0,192,39,0,0,91,0,0,0,92,0,0,0,44,0,0,0,51,0,0,0,52,0,0,0,36,0,0,0,37,0,0,0,38,0,0,0,116,114,117,101,0,0,0,0,116,0,0,0,114,0,0,0,117,0,0,0,101,0,0,0,0,0,0,0,0,0,0,0,102,97,108,115,101,0,0,0,102,0,0,0,97,0,0,0,108,0,0,0,115,0,0,0,101,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,37,109,47,37,100,47,37,121,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,37,0,0,0,109,0,0,0,47,0,0,0,37,0,0,0,100,0,0,0,47,0,0,0,37,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,37,72,58,37,77,58,37,83,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,37,0,0,0,72,0,0,0,58,0,0,0,37,0,0,0,77,0,0,0,58,0,0,0,37,0,0,0,83,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,37,97,32,37,98,32,37,100,32,37,72,58,37,77,58,37,83,32,37,89,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,37,0,0,0,97,0,0,0,32,0,0,0,37,0,0,0,98,0,0,0,32,0,0,0,37,0,0,0,100,0,0,0,32,0,0,0,37,0,0,0,72,0,0,0,58,0,0,0,37,0,0,0,77,0,0,0,58,0,0,0,37,0,0,0,83,0,0,0,32,0,0,0,37,0,0,0,89,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,37,73,58,37,77,58,37,83,32,37,112,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,37,0,0,0,73,0,0,0,58,0,0,0,37,0,0,0,77,0,0,0,58,0,0,0,37,0,0,0,83,0,0,0,32,0,0,0,37,0,0,0,112,0,0,0,0,0,0,0,108,111,99,97,108,101,32,110,111,116,32,115,117,112,112,111,114,116,101,100,0,0,0,0,0,0,0,0,248,35,0,0,93,0,0,0,94,0,0,0,44,0,0,0,0,0,0,0,78,83,116,51,95,95,49,54,108,111,99,97,108,101,53,102,97,99,101,116,69,0,0,0,8,64,0,0,224,35,0,0,224,19,0,0,0,0,0,0,0,0,0,0,136,36,0,0,93,0,0,0,95,0,0,0,44,0,0,0,18,0,0,0,2,0,0,0,3,0,0,0,4,0,0,0,12,0,0,0,19,0,0,0,13,0,0,0,20,0,0,0,14,0,0,0,5,0,0,0,21,0,0,0,6,0,0,0,0,0,0,0,78,83,116,51,95,95,49,53,99,116,121,112,101,73,119,69,69,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,48,99,116,121,112,101,95,98,97,115,101,69,0,0,0,0,224,63,0,0,104,36,0,0,104,64,0,0,80,36,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,128,36,0,0,2,0,0,0,78,83,116,51,95,95,49,53,99,116,121,112,101,73,99,69,69,0,0,0,0,0,0,0,104,64,0,0,168,36,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,128,36,0,0,2,0,0,0,0,0,0,0,88,37,0,0,93,0,0,0,96,0,0,0,44,0,0,0,3,0,0,0,4,0,0,0,7,0,0,0,53,0,0,0,54,0,0,0,8,0,0,0,55,0,0,0,78,83,116,51,95,95,49,55,99,111,100,101,99,118,116,73,99,99,49,49,95,95,109,98,115,116,97,116,101,95,116,69,69,0,0,0,0,0,0,0,78,83,116,51,95,95,49,49,50,99,111,100,101,99,118,116,95,98,97,115,101,69,0,0,224,63,0,0,56,37,0,0,104,64,0,0,16,37,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,80,37,0,0,2,0,0,0,78,83,116,51,95,95,49,55,99,111,100,101,99,118,116,73,119,99,49,49,95,95,109,98,115,116,97,116,101,95,116,69,69,0,0,0,0,0,0,0,104,64,0,0,120,37,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,80,37,0,0,2,0,0,0,0,0,0,0,24,38,0,0,93,0,0,0,97,0,0,0,44,0,0,0,5,0,0,0,6,0,0,0,9,0,0,0,56,0,0,0,57,0,0,0,10,0,0,0,58,0,0,0,78,83,116,51,95,95,49,55,99,111,100,101,99,118,116,73,68,115,99,49,49,95,95,109,98,115,116,97,116,101,95,116,69,69,0,0,0,0,0,0,104,64,0,0,240,37,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,80,37,0,0,2,0,0,0,0,0,0,0,144,38,0,0,93,0,0,0,98,0,0,0,44,0,0,0,7,0,0,0,8,0,0,0,11,0,0,0,59,0,0,0,60,0,0,0,12,0,0,0,61,0,0,0,78,83,116,51,95,95,49,55,99,111,100,101,99,118,116,73,68,105,99,49,49,95,95,109,98,115,116,97,116,101,95,116,69,69,0,0,0,0,0,0,104,64,0,0,104,38,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,80,37,0,0,2,0,0,0,0,0,0,0,8,39,0,0,93,0,0,0,99,0,0,0,44,0,0,0,7,0,0,0,8,0,0,0,11,0,0,0,59,0,0,0,60,0,0,0,12,0,0,0,61,0,0,0,78,83,116,51,95,95,49,49,54,95,95,110,97,114,114,111,119,95,116,111,95,117,116,102,56,73,76,106,51,50,69,69,69,0,0,0,0,0,0,0,8,64,0,0,224,38,0,0,144,38,0,0,0,0,0,0,0,0,0,0,112,39,0,0,93,0,0,0,100,0,0,0,44,0,0,0,7,0,0,0,8,0,0,0,11,0,0,0,59,0,0,0,60,0,0,0,12,0,0,0,61,0,0,0,78,83,116,51,95,95,49,49,55,95,95,119,105,100,101,110,95,102,114,111,109,95,117,116,102,56,73,76,106,51,50,69,69,69,0,0,0,0,0,0,8,64,0,0,72,39,0,0,144,38,0,0,0,0,0,0,78,83,116,51,95,95,49,56,110,117,109,112,117,110,99,116,73,99,69,69,0,0,0,0,8,64,0,0,128,39,0,0,248,35,0,0,0,0,0,0,78,83,116,51,95,95,49,56,110,117,109,112,117,110,99,116,73,119,69,69,0,0,0,0,8,64,0,0,168,39,0,0,248,35,0,0,0,0,0,0,78,83,116,51,95,95,49,54,108,111,99,97,108,101,53,95,95,105,109,112,69,0,0,0,8,64,0,0,208,39,0,0,248,35,0,0,0,0,0,0,78,83,116,51,95,95,49,55,99,111,108,108,97,116,101,73,99,69,69,0,0,0,0,0,8,64,0,0,248,39,0,0,248,35,0,0,0,0,0,0,78,83,116,51,95,95,49,55,99,111,108,108,97,116,101,73,119,69,69,0,0,0,0,0,8,64,0,0,32,40,0,0,248,35,0,0,0,0,0,0,78,83,116,51,95,95,49,49,48,109,111,110,101,121,112,117,110,99,116,73,99,76,98,48,69,69,69,0,0,0,0,0,78,83,116,51,95,95,49,49,48,109,111,110,101,121,95,98,97,115,101,69,0,0,0,0,224,63,0,0,104,40,0,0,104,64,0,0,72,40,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,128,40,0,0,2,0,0,0,78,83,116,51,95,95,49,49,48,109,111,110,101,121,112,117,110,99,116,73,99,76,98,49,69,69,69,0,0,0,0,0,104,64,0,0,168,40,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,128,40,0,0,2,0,0,0,78,83,116,51,95,95,49,49,48,109,111,110,101,121,112,117,110,99,116,73,119,76,98,48,69,69,69,0,0,0,0,0,104,64,0,0,232,40,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,128,40,0,0,2,0,0,0,78,83,116,51,95,95,49,49,48,109,111,110,101,121,112,117,110,99,116,73,119,76,98,49,69,69,69,0,0,0,0,0,104,64,0,0,40,41,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,128,40,0,0,2,0,0,0,78,83,116,51,95,95,49,56,116,105,109,101,95,103,101,116,73,99,78,83,95,49,57,105,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,69,69,0,0,0,0,78,83,116,51,95,95,49,57,116,105,109,101,95,98,97,115,101,69,0,0,0,0,0,0,224,63,0,0,176,41,0,0,78,83,116,51,95,95,49,50,48,95,95,116,105,109,101,95,103,101,116,95,99,95,115,116,111,114,97,103,101,73,99,69,69,0,0,0,0,0,0,0,224,63,0,0,208,41,0,0,104,64,0,0,104,41,0,0,0,0,0,0,3,0,0,0,248,35,0,0,2,0,0,0,200,41,0,0,2,0,0,0,248,41,0,0,0,8,0,0,78,83,116,51,95,95,49,56,116,105,109,101,95,103,101,116,73,119,78,83,95,49,57,105,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,69,69,0,0,0,0,78,83,116,51,95,95,49,50,48,95,95,116,105,109,101,95,103,101,116,95,99,95,115,116,111,114,97,103,101,73,119,69,69,0,0,0,0,0,0,0,224,63,0,0,112,42,0,0,104,64,0,0,40,42,0,0,0,0,0,0,3,0,0,0,248,35,0,0,2,0,0,0,200,41,0,0,2,0,0,0,152,42,0,0,0,8,0,0,78,83,116,51,95,95,49,56,116,105,109,101,95,112,117,116,73,99,78,83,95,49,57,111,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,69,69,0,0,0,0,78,83,116,51,95,95,49,49,48,95,95,116,105,109,101,95,112,117,116,69,0,0,0,0,224,63,0,0,16,43,0,0,104,64,0,0,200,42,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,40,43,0,0,0,8,0,0,78,83,116,51,95,95,49,56,116,105,109,101,95,112,117,116,73,119,78,83,95,49,57,111,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,69,69,0,0,0,0,104,64,0,0,80,43,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,40,43,0,0,0,8,0,0,78,83,116,51,95,95,49,56,109,101,115,115,97,103,101,115,73,99,69,69,0,0,0,0,78,83,116,51,95,95,49,49,51,109,101,115,115,97,103,101,115,95,98,97,115,101,69,0,224,63,0,0,208,43,0,0,104,64,0,0,184,43,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,232,43,0,0,2,0,0,0,78,83,116,51,95,95,49,56,109,101,115,115,97,103,101,115,73,119,69,69,0,0,0,0,104,64,0,0,16,44,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,232,43,0,0,2,0,0,0,78,83,116,51,95,95,49,55,110,117,109,95,103,101,116,73,99,78,83,95,49,57,105,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,69,69,0,0,0,0,0,78,83,116,51,95,95,49,57,95,95,110,117,109,95,103,101,116,73,99,69,69,0,0,0,78,83,116,51,95,95,49,49,52,95,95,110,117,109,95,103,101,116,95,98,97,115,101,69,0,0,0,0,0,0,0,0,224,63,0,0,168,44,0,0,104,64,0,0,144,44,0,0,0,0,0,0,1,0,0,0,200,44,0,0,0,0,0,0,104,64,0,0,72,44,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,208,44,0,0,0,0,0,0,78,83,116,51,95,95,49,55,110,117,109,95,103,101,116,73,119,78,83,95,49,57,105,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,69,69,0,0,0,0,0,78,83,116,51,95,95,49,57,95,95,110,117,109,95,103,101,116,73,119,69,69,0,0,0,104,64,0,0,80,45,0,0,0,0,0,0,1,0,0,0,200,44,0,0,0,0,0,0,104,64,0,0,8,45,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,104,45,0,0,0,0,0,0,78,83,116,51,95,95,49,55,110,117,109,95,112,117,116,73,99,78,83,95,49,57,111,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,69,69,0,0,0,0,0,78,83,116,51,95,95,49,57,95,95,110,117,109,95,112,117,116,73,99,69,69,0,0,0,78,83,116,51,95,95,49,49,52,95,95,110,117,109,95,112,117,116,95,98,97,115,101,69,0,0,0,0,0,0,0,0,224,63,0,0,0,46,0,0,104,64,0,0,232,45,0,0,0,0,0,0,1,0,0,0,32,46,0,0,0,0,0,0,104,64,0,0,160,45,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,40,46,0,0,0,0,0,0,78,83,116,51,95,95,49,55,110,117,109,95,112,117,116,73,119,78,83,95,49,57,111,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,69,69,0,0,0,0,0,78,83,116,51,95,95,49,57,95,95,110,117,109,95,112,117,116,73,119,69,69,0,0,0,104,64,0,0,168,46,0,0,0,0,0,0,1,0,0,0,32,46,0,0,0,0,0,0,104,64,0,0,96,46,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,192,46,0,0,0,0,0,0,78,83,116,51,95,95,49,57,109,111,110,101,121,95,103,101,116,73,99,78,83,95,49,57,105,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,69,69,0,0,0,78,83,116,51,95,95,49,49,49,95,95,109,111,110,101,121,95,103,101,116,73,99,69,69,0,0,0,0,0,0,0,0,224,63,0,0,64,47,0,0,104,64,0,0,248,46,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,96,47,0,0,0,0,0,0,78,83,116,51,95,95,49,57,109,111,110,101,121,95,103,101,116,73,119,78,83,95,49,57,105,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,69,69,0,0,0,78,83,116,51,95,95,49,49,49,95,95,109,111,110,101,121,95,103,101,116,73,119,69,69,0,0,0,0,0,0,0,0,224,63,0,0,208,47,0,0,104,64,0,0,136,47,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,240,47,0,0,0,0,0,0,78,83,116,51,95,95,49,57,109,111,110,101,121,95,112,117,116,73,99,78,83,95,49,57,111,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,99,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,99,69,69,69,69,69,69,0,0,0,78,83,116,51,95,95,49,49,49,95,95,109,111,110,101,121,95,112,117,116,73,99,69,69,0,0,0,0,0,0,0,0,224,63,0,0,96,48,0,0,104,64,0,0,24,48,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,128,48,0,0,0,0,0,0,78,83,116,51,95,95,49,57,109,111,110,101,121,95,112,117,116,73,119,78,83,95,49,57,111,115,116,114,101,97,109,98,117,102,95,105,116,101,114,97,116,111,114,73,119,78,83,95,49,49,99,104,97,114,95,116,114,97,105,116,115,73,119,69,69,69,69,69,69,0,0,0,78,83,116,51,95,95,49,49,49,95,95,109,111,110,101,121,95,112,117,116,73,119,69,69,0,0,0,0,0,0,0,0,224,63,0,0,240,48,0,0,104,64,0,0,168,48,0,0,0,0,0,0,2,0,0,0,248,35,0,0,2,0,0,0,16,49,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,65,0,0,0,77,0,0,0,0,0,0,0,0,0,0,0,80,0,0,0,77,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,65,77,0,0,0,0,0,0,80,77,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,74,0,0,0,97,0,0,0,110,0,0,0,117,0,0,0,97,0,0,0,114,0,0,0,121,0,0,0,0,0,0,0,70,0,0,0,101,0,0,0,98,0,0,0,114,0,0,0,117,0,0,0,97,0,0,0,114,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,77,0,0,0,97,0,0,0,114,0,0,0,99,0,0,0,104,0,0,0,0,0,0,0,65,0,0,0,112,0,0,0,114,0,0,0,105,0,0,0,108,0,0,0,0,0,0,0,74,0,0,0,117,0,0,0,110,0,0,0,101,0,0,0,0,0,0,0,0,0,0,0,74,0,0,0,117,0,0,0,108,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,65,0,0,0,117,0,0,0,103,0,0,0,117,0,0,0,115,0,0,0,116,0,0,0,0,0,0,0,0,0,0,0,83,0,0,0,101,0,0,0,112,0,0,0,116,0,0,0,101,0,0,0,109,0,0,0,98,0,0,0,101,0,0,0,114,0,0,0,0,0,0,0,79,0,0,0,99,0,0,0,116,0,0,0,111,0,0,0,98,0,0,0,101,0,0,0,114,0,0,0,0,0,0,0,78,0,0,0,111,0,0,0,118,0,0,0,101,0,0,0,109,0,0,0,98,0,0,0,101,0,0,0,114,0,0,0,0,0,0,0,0,0,0,0,68,0,0,0,101,0,0,0,99,0,0,0,101,0,0,0,109,0,0,0,98,0,0,0,101,0,0,0,114,0,0,0,0,0,0,0,0,0,0,0,74,0,0,0,97,0,0,0,110,0,0,0,0,0,0,0,70,0,0,0,101,0,0,0,98,0,0,0,0,0,0,0,77,0,0,0,97,0,0,0,114,0,0,0,0,0,0,0,65,0,0,0,112,0,0,0,114,0,0,0,0,0,0,0,77,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,74,0,0,0,117,0,0,0,110,0,0,0,0,0,0,0,74,0,0,0,117,0,0,0,108,0,0,0,0,0,0,0,65,0,0,0,117,0,0,0,103,0,0,0,0,0,0,0,83,0,0,0,101,0,0,0,112,0,0,0,0,0,0,0,79,0,0,0,99,0,0,0,116,0,0,0,0,0,0,0,78,0,0,0,111,0,0,0,118,0,0,0,0,0,0,0,68,0,0,0,101,0,0,0,99], "i8", ALLOC_NONE, Runtime.GLOBAL_BASE+4580);
/* memory initializer */ allocate([74,97,110,117,97,114,121,0,70,101,98,114,117,97,114,121,0,0,0,0,0,0,0,0,77,97,114,99,104,0,0,0,65,112,114,105,108,0,0,0,77,97,121,0,0,0,0,0,74,117,110,101,0,0,0,0,74,117,108,121,0,0,0,0,65,117,103,117,115,116,0,0,83,101,112,116,101,109,98,101,114,0,0,0,0,0,0,0,79,99,116,111,98,101,114,0,78,111,118,101,109,98,101,114,0,0,0,0,0,0,0,0,68,101,99,101,109,98,101,114,0,0,0,0,0,0,0,0,74,97,110,0,0,0,0,0,70,101,98,0,0,0,0,0,77,97,114,0,0,0,0,0,65,112,114,0,0,0,0,0,74,117,110,0,0,0,0,0,74,117,108,0,0,0,0,0,65,117,103,0,0,0,0,0,83,101,112,0,0,0,0,0,79,99,116,0,0,0,0,0,78,111,118,0,0,0,0,0,68,101,99,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,83,0,0,0,117,0,0,0,110,0,0,0,100,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,77,0,0,0,111,0,0,0,110,0,0,0,100,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,84,0,0,0,117,0,0,0,101,0,0,0,115,0,0,0,100,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,87,0,0,0,101,0,0,0,100,0,0,0,110,0,0,0,101,0,0,0,115,0,0,0,100,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,84,0,0,0,104,0,0,0,117,0,0,0,114,0,0,0,115,0,0,0,100,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,70,0,0,0,114,0,0,0,105,0,0,0,100,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,83,0,0,0,97,0,0,0,116,0,0,0,117,0,0,0,114,0,0,0,100,0,0,0,97,0,0,0,121,0,0,0,0,0,0,0,0,0,0,0,83,0,0,0,117,0,0,0,110,0,0,0,0,0,0,0,77,0,0,0,111,0,0,0,110,0,0,0,0,0,0,0,84,0,0,0,117,0,0,0,101,0,0,0,0,0,0,0,87,0,0,0,101,0,0,0,100,0,0,0,0,0,0,0,84,0,0,0,104,0,0,0,117,0,0,0,0,0,0,0,70,0,0,0,114,0,0,0,105,0,0,0,0,0,0,0,83,0,0,0,97,0,0,0,116,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,83,117,110,100,97,121,0,0,77,111,110,100,97,121,0,0,84,117,101,115,100,97,121,0,87,101,100,110,101,115,100,97,121,0,0,0,0,0,0,0,84,104,117,114,115,100,97,121,0,0,0,0,0,0,0,0,70,114,105,100,97,121,0,0,83,97,116,117,114,100,97,121,0,0,0,0,0,0,0,0,83,117,110,0,0,0,0,0,77,111,110,0,0,0,0,0,84,117,101,0,0,0,0,0,87,101,100,0,0,0,0,0,84,104,117,0,0,0,0,0,70,114,105,0,0,0,0,0,83,97,116,0,0,0,0,0,2,0,0,192,3,0,0,192,4,0,0,192,5,0,0,192,6,0,0,192,7,0,0,192,8,0,0,192,9,0,0,192,10,0,0,192,11,0,0,192,12,0,0,192,13,0,0,192,14,0,0,192,15,0,0,192,16,0,0,192,17,0,0,192,18,0,0,192,19,0,0,192,20,0,0,192,21,0,0,192,22,0,0,192,23,0,0,192,24,0,0,192,25,0,0,192,26,0,0,192,27,0,0,192,28,0,0,192,29,0,0,192,30,0,0,192,31,0,0,192,0,0,0,179,1,0,0,195,2,0,0,195,3,0,0,195,4,0,0,195,5,0,0,195,6,0,0,195,7,0,0,195,8,0,0,195,9,0,0,195,10,0,0,195,11,0,0,195,12,0,0,195,13,0,0,211,14,0,0,195,15,0,0,195,0,0,12,187,1,0,12,195,2,0,12,195,3,0,12,195,4,0,12,211,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,88,63,0,0,101,0,0,0,102,0,0,0,62,0,0,0,0,0,0,0,115,116,100,58,58,98,97,100,95,99,97,115,116,0,0,0,83,116,57,116,121,112,101,95,105,110,102,111,0,0,0,0,224,63,0,0,48,63,0,0,83,116,56,98,97,100,95,99,97,115,116,0,0,0,0,0,8,64,0,0,72,63,0,0,0,0,0,0,0,0,0,0,78,49,48,95,95,99,120,120,97,98,105,118,49,49,54,95,95,115,104,105,109,95,116,121,112,101,95,105,110,102,111,69,0,0,0,0,0,0,0,0,8,64,0,0,104,63,0,0,64,63,0,0,0,0,0,0,78,49,48,95,95,99,120,120,97,98,105,118,49,49,55,95,95,99,108,97,115,115,95,116,121,112,101,95,105,110,102,111,69,0,0,0,0,0,0,0,8,64,0,0,160,63,0,0,144,63,0,0,0,0,0,0,0,0,0,0,200,63,0,0,103,0,0,0,104,0,0,0,105,0,0,0,106,0,0,0,22,0,0,0,13,0,0,0,1,0,0,0,5,0,0,0,0,0,0,0,80,64,0,0,103,0,0,0,107,0,0,0,105,0,0,0,106,0,0,0,22,0,0,0,14,0,0,0,2,0,0,0,6,0,0,0,78,49,48,95,95,99,120,120,97,98,105,118,49,50,48,95,95,115,105,95,99,108,97,115,115,95,116,121,112,101,95,105,110,102,111,69,0,0,0,0,8,64,0,0,40,64,0,0,200,63,0,0,0,0,0,0,0,0,0,0,176,64,0,0,103,0,0,0,108,0,0,0,105,0,0,0,106,0,0,0,22,0,0,0,15,0,0,0,3,0,0,0,7,0,0,0,78,49,48,95,95,99,120,120,97,98,105,118,49,50,49,95,95,118,109,105,95,99,108,97,115,115,95,116,121,112,101,95,105,110,102,111,69,0,0,0,8,64,0,0,136,64,0,0,200,63,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,240,66,0,0,109,0,0,0,110,0,0,0,63,0,0,0,0,0,0,0,115,116,100,58,58,98,97,100,95,97,108,108,111,99,0,0,83,116,57,98,97,100,95,97,108,108,111,99,0,0,0,0,8,64,0,0,224,66,0,0,0,0,0,0,0,0,0,0,105,110,102,105,110,105,116,121,0,0,0,0,0,0,0,0,110,97,110,0,0,0,0,0,95,112,137,0,255,9,47,15,10,0,0,0,100,0,0,0,232,3,0,0,16,39,0,0,160,134,1,0,64,66,15,0,128,150,152,0,0,225,245,5], "i8", ALLOC_NONE, Runtime.GLOBAL_BASE+14840);




var tempDoublePtr = Runtime.alignMemory(allocate(12, "i8", ALLOC_STATIC), 8);

assert(tempDoublePtr % 8 == 0);

function copyTempFloat(ptr) { // functions, because inlining this code increases code size too much

  HEAP8[tempDoublePtr] = HEAP8[ptr];

  HEAP8[tempDoublePtr+1] = HEAP8[ptr+1];

  HEAP8[tempDoublePtr+2] = HEAP8[ptr+2];

  HEAP8[tempDoublePtr+3] = HEAP8[ptr+3];

}

function copyTempDouble(ptr) {

  HEAP8[tempDoublePtr] = HEAP8[ptr];

  HEAP8[tempDoublePtr+1] = HEAP8[ptr+1];

  HEAP8[tempDoublePtr+2] = HEAP8[ptr+2];

  HEAP8[tempDoublePtr+3] = HEAP8[ptr+3];

  HEAP8[tempDoublePtr+4] = HEAP8[ptr+4];

  HEAP8[tempDoublePtr+5] = HEAP8[ptr+5];

  HEAP8[tempDoublePtr+6] = HEAP8[ptr+6];

  HEAP8[tempDoublePtr+7] = HEAP8[ptr+7];

}


  
  function _atexit(func, arg) {
      __ATEXIT__.unshift({ func: func, arg: arg });
    }var ___cxa_atexit=_atexit;

   
  Module["_i64Subtract"] = _i64Subtract;

  var _fabsf=Math_abs;

   
  Module["_i64Add"] = _i64Add;

  
  function __ZSt18uncaught_exceptionv() { // std::uncaught_exception()
      return !!__ZSt18uncaught_exceptionv.uncaught_exception;
    }
  
  
  
  function ___cxa_is_number_type(type) {
      var isNumber = false;
      try { if (type == __ZTIi) isNumber = true } catch(e){}
      try { if (type == __ZTIj) isNumber = true } catch(e){}
      try { if (type == __ZTIl) isNumber = true } catch(e){}
      try { if (type == __ZTIm) isNumber = true } catch(e){}
      try { if (type == __ZTIx) isNumber = true } catch(e){}
      try { if (type == __ZTIy) isNumber = true } catch(e){}
      try { if (type == __ZTIf) isNumber = true } catch(e){}
      try { if (type == __ZTId) isNumber = true } catch(e){}
      try { if (type == __ZTIe) isNumber = true } catch(e){}
      try { if (type == __ZTIc) isNumber = true } catch(e){}
      try { if (type == __ZTIa) isNumber = true } catch(e){}
      try { if (type == __ZTIh) isNumber = true } catch(e){}
      try { if (type == __ZTIs) isNumber = true } catch(e){}
      try { if (type == __ZTIt) isNumber = true } catch(e){}
      return isNumber;
    }function ___cxa_does_inherit(definiteType, possibilityType, possibility) {
      if (possibility == 0) return false;
      if (possibilityType == 0 || possibilityType == definiteType)
        return true;
      var possibility_type_info;
      if (___cxa_is_number_type(possibilityType)) {
        possibility_type_info = possibilityType;
      } else {
        var possibility_type_infoAddr = HEAP32[((possibilityType)>>2)] - 8;
        possibility_type_info = HEAP32[((possibility_type_infoAddr)>>2)];
      }
      switch (possibility_type_info) {
      case 0: // possibility is a pointer
        // See if definite type is a pointer
        var definite_type_infoAddr = HEAP32[((definiteType)>>2)] - 8;
        var definite_type_info = HEAP32[((definite_type_infoAddr)>>2)];
        if (definite_type_info == 0) {
          // Also a pointer; compare base types of pointers
          var defPointerBaseAddr = definiteType+8;
          var defPointerBaseType = HEAP32[((defPointerBaseAddr)>>2)];
          var possPointerBaseAddr = possibilityType+8;
          var possPointerBaseType = HEAP32[((possPointerBaseAddr)>>2)];
          return ___cxa_does_inherit(defPointerBaseType, possPointerBaseType, possibility);
        } else
          return false; // one pointer and one non-pointer
      case 1: // class with no base class
        return false;
      case 2: // class with base class
        var parentTypeAddr = possibilityType + 8;
        var parentType = HEAP32[((parentTypeAddr)>>2)];
        return ___cxa_does_inherit(definiteType, parentType, possibility);
      default:
        return false; // some unencountered type
      }
    }
  
  
  
  var ___cxa_last_thrown_exception=0;function ___resumeException(ptr) {
      if (!___cxa_last_thrown_exception) { ___cxa_last_thrown_exception = ptr; }
      throw ptr + " - Exception catching is disabled, this exception cannot be caught. Compile with -s DISABLE_EXCEPTION_CATCHING=0 or DISABLE_EXCEPTION_CATCHING=2 to catch.";
    }
  
  var ___cxa_exception_header_size=8;function ___cxa_find_matching_catch(thrown, throwntype) {
      if (thrown == -1) thrown = ___cxa_last_thrown_exception;
      header = thrown - ___cxa_exception_header_size;
      if (throwntype == -1) throwntype = HEAP32[((header)>>2)];
      var typeArray = Array.prototype.slice.call(arguments, 2);
  
      // If throwntype is a pointer, this means a pointer has been
      // thrown. When a pointer is thrown, actually what's thrown
      // is a pointer to the pointer. We'll dereference it.
      if (throwntype != 0 && !___cxa_is_number_type(throwntype)) {
        var throwntypeInfoAddr= HEAP32[((throwntype)>>2)] - 8;
        var throwntypeInfo= HEAP32[((throwntypeInfoAddr)>>2)];
        if (throwntypeInfo == 0)
          thrown = HEAP32[((thrown)>>2)];
      }
      // The different catch blocks are denoted by different types.
      // Due to inheritance, those types may not precisely match the
      // type of the thrown object. Find one which matches, and
      // return the type of the catch block which should be called.
      for (var i = 0; i < typeArray.length; i++) {
        if (___cxa_does_inherit(typeArray[i], throwntype, thrown))
          return ((asm["setTempRet0"](typeArray[i]),thrown)|0);
      }
      // Shouldn't happen unless we have bogus data in typeArray
      // or encounter a type for which emscripten doesn't have suitable
      // typeinfo defined. Best-efforts match just in case.
      return ((asm["setTempRet0"](throwntype),thrown)|0);
    }function ___cxa_throw(ptr, type, destructor) {
      if (!___cxa_throw.initialized) {
        try {
          HEAP32[((__ZTVN10__cxxabiv119__pointer_type_infoE)>>2)]=0; // Workaround for libcxxabi integration bug
        } catch(e){}
        try {
          HEAP32[((__ZTVN10__cxxabiv117__class_type_infoE)>>2)]=1; // Workaround for libcxxabi integration bug
        } catch(e){}
        try {
          HEAP32[((__ZTVN10__cxxabiv120__si_class_type_infoE)>>2)]=2; // Workaround for libcxxabi integration bug
        } catch(e){}
        ___cxa_throw.initialized = true;
      }
      var header = ptr - ___cxa_exception_header_size;
      HEAP32[((header)>>2)]=type;
      HEAP32[(((header)+(4))>>2)]=destructor;
      ___cxa_last_thrown_exception = ptr;
      if (!("uncaught_exception" in __ZSt18uncaught_exceptionv)) {
        __ZSt18uncaught_exceptionv.uncaught_exception = 1;
      } else {
        __ZSt18uncaught_exceptionv.uncaught_exception++;
      }
      throw ptr + " - Exception catching is disabled, this exception cannot be caught. Compile with -s DISABLE_EXCEPTION_CATCHING=0 or DISABLE_EXCEPTION_CATCHING=2 to catch.";
    }

  function _pthread_mutex_lock() {}

  
  
  var ERRNO_CODES={EPERM:1,ENOENT:2,ESRCH:3,EINTR:4,EIO:5,ENXIO:6,E2BIG:7,ENOEXEC:8,EBADF:9,ECHILD:10,EAGAIN:11,EWOULDBLOCK:11,ENOMEM:12,EACCES:13,EFAULT:14,ENOTBLK:15,EBUSY:16,EEXIST:17,EXDEV:18,ENODEV:19,ENOTDIR:20,EISDIR:21,EINVAL:22,ENFILE:23,EMFILE:24,ENOTTY:25,ETXTBSY:26,EFBIG:27,ENOSPC:28,ESPIPE:29,EROFS:30,EMLINK:31,EPIPE:32,EDOM:33,ERANGE:34,ENOMSG:42,EIDRM:43,ECHRNG:44,EL2NSYNC:45,EL3HLT:46,EL3RST:47,ELNRNG:48,EUNATCH:49,ENOCSI:50,EL2HLT:51,EDEADLK:35,ENOLCK:37,EBADE:52,EBADR:53,EXFULL:54,ENOANO:55,EBADRQC:56,EBADSLT:57,EDEADLOCK:35,EBFONT:59,ENOSTR:60,ENODATA:61,ETIME:62,ENOSR:63,ENONET:64,ENOPKG:65,EREMOTE:66,ENOLINK:67,EADV:68,ESRMNT:69,ECOMM:70,EPROTO:71,EMULTIHOP:72,EDOTDOT:73,EBADMSG:74,ENOTUNIQ:76,EBADFD:77,EREMCHG:78,ELIBACC:79,ELIBBAD:80,ELIBSCN:81,ELIBMAX:82,ELIBEXEC:83,ENOSYS:38,ENOTEMPTY:39,ENAMETOOLONG:36,ELOOP:40,EOPNOTSUPP:95,EPFNOSUPPORT:96,ECONNRESET:104,ENOBUFS:105,EAFNOSUPPORT:97,EPROTOTYPE:91,ENOTSOCK:88,ENOPROTOOPT:92,ESHUTDOWN:108,ECONNREFUSED:111,EADDRINUSE:98,ECONNABORTED:103,ENETUNREACH:101,ENETDOWN:100,ETIMEDOUT:110,EHOSTDOWN:112,EHOSTUNREACH:113,EINPROGRESS:115,EALREADY:114,EDESTADDRREQ:89,EMSGSIZE:90,EPROTONOSUPPORT:93,ESOCKTNOSUPPORT:94,EADDRNOTAVAIL:99,ENETRESET:102,EISCONN:106,ENOTCONN:107,ETOOMANYREFS:109,EUSERS:87,EDQUOT:122,ESTALE:116,ENOTSUP:95,ENOMEDIUM:123,EILSEQ:84,EOVERFLOW:75,ECANCELED:125,ENOTRECOVERABLE:131,EOWNERDEAD:130,ESTRPIPE:86};
  
  var ERRNO_MESSAGES={0:"Success",1:"Not super-user",2:"No such file or directory",3:"No such process",4:"Interrupted system call",5:"I/O error",6:"No such device or address",7:"Arg list too long",8:"Exec format error",9:"Bad file number",10:"No children",11:"No more processes",12:"Not enough core",13:"Permission denied",14:"Bad address",15:"Block device required",16:"Mount device busy",17:"File exists",18:"Cross-device link",19:"No such device",20:"Not a directory",21:"Is a directory",22:"Invalid argument",23:"Too many open files in system",24:"Too many open files",25:"Not a typewriter",26:"Text file busy",27:"File too large",28:"No space left on device",29:"Illegal seek",30:"Read only file system",31:"Too many links",32:"Broken pipe",33:"Math arg out of domain of func",34:"Math result not representable",35:"File locking deadlock error",36:"File or path name too long",37:"No record locks available",38:"Function not implemented",39:"Directory not empty",40:"Too many symbolic links",42:"No message of desired type",43:"Identifier removed",44:"Channel number out of range",45:"Level 2 not synchronized",46:"Level 3 halted",47:"Level 3 reset",48:"Link number out of range",49:"Protocol driver not attached",50:"No CSI structure available",51:"Level 2 halted",52:"Invalid exchange",53:"Invalid request descriptor",54:"Exchange full",55:"No anode",56:"Invalid request code",57:"Invalid slot",59:"Bad font file fmt",60:"Device not a stream",61:"No data (for no delay io)",62:"Timer expired",63:"Out of streams resources",64:"Machine is not on the network",65:"Package not installed",66:"The object is remote",67:"The link has been severed",68:"Advertise error",69:"Srmount error",70:"Communication error on send",71:"Protocol error",72:"Multihop attempted",73:"Cross mount point (not really error)",74:"Trying to read unreadable message",75:"Value too large for defined data type",76:"Given log. name not unique",77:"f.d. invalid for this operation",78:"Remote address changed",79:"Can   access a needed shared lib",80:"Accessing a corrupted shared lib",81:".lib section in a.out corrupted",82:"Attempting to link in too many libs",83:"Attempting to exec a shared library",84:"Illegal byte sequence",86:"Streams pipe error",87:"Too many users",88:"Socket operation on non-socket",89:"Destination address required",90:"Message too long",91:"Protocol wrong type for socket",92:"Protocol not available",93:"Unknown protocol",94:"Socket type not supported",95:"Not supported",96:"Protocol family not supported",97:"Address family not supported by protocol family",98:"Address already in use",99:"Address not available",100:"Network interface is not configured",101:"Network is unreachable",102:"Connection reset by network",103:"Connection aborted",104:"Connection reset by peer",105:"No buffer space available",106:"Socket is already connected",107:"Socket is not connected",108:"Can't send after socket shutdown",109:"Too many references",110:"Connection timed out",111:"Connection refused",112:"Host is down",113:"Host is unreachable",114:"Socket already connected",115:"Connection already in progress",116:"Stale file handle",122:"Quota exceeded",123:"No medium (in tape drive)",125:"Operation canceled",130:"Previous owner died",131:"State not recoverable"};
  
  
  var ___errno_state=0;function ___setErrNo(value) {
      // For convenient setting and returning of errno.
      HEAP32[((___errno_state)>>2)]=value;
      return value;
    }
  
  var PATH={splitPath:function (filename) {
        var splitPathRe = /^(\/?|)([\s\S]*?)((?:\.{1,2}|[^\/]+?|)(\.[^.\/]*|))(?:[\/]*)$/;
        return splitPathRe.exec(filename).slice(1);
      },normalizeArray:function (parts, allowAboveRoot) {
        // if the path tries to go above the root, `up` ends up > 0
        var up = 0;
        for (var i = parts.length - 1; i >= 0; i--) {
          var last = parts[i];
          if (last === '.') {
            parts.splice(i, 1);
          } else if (last === '..') {
            parts.splice(i, 1);
            up++;
          } else if (up) {
            parts.splice(i, 1);
            up--;
          }
        }
        // if the path is allowed to go above the root, restore leading ..s
        if (allowAboveRoot) {
          for (; up--; up) {
            parts.unshift('..');
          }
        }
        return parts;
      },normalize:function (path) {
        var isAbsolute = path.charAt(0) === '/',
            trailingSlash = path.substr(-1) === '/';
        // Normalize the path
        path = PATH.normalizeArray(path.split('/').filter(function(p) {
          return !!p;
        }), !isAbsolute).join('/');
        if (!path && !isAbsolute) {
          path = '.';
        }
        if (path && trailingSlash) {
          path += '/';
        }
        return (isAbsolute ? '/' : '') + path;
      },dirname:function (path) {
        var result = PATH.splitPath(path),
            root = result[0],
            dir = result[1];
        if (!root && !dir) {
          // No dirname whatsoever
          return '.';
        }
        if (dir) {
          // It has a dirname, strip trailing slash
          dir = dir.substr(0, dir.length - 1);
        }
        return root + dir;
      },basename:function (path) {
        // EMSCRIPTEN return '/'' for '/', not an empty string
        if (path === '/') return '/';
        var lastSlash = path.lastIndexOf('/');
        if (lastSlash === -1) return path;
        return path.substr(lastSlash+1);
      },extname:function (path) {
        return PATH.splitPath(path)[3];
      },join:function () {
        var paths = Array.prototype.slice.call(arguments, 0);
        return PATH.normalize(paths.join('/'));
      },join2:function (l, r) {
        return PATH.normalize(l + '/' + r);
      },resolve:function () {
        var resolvedPath = '',
          resolvedAbsolute = false;
        for (var i = arguments.length - 1; i >= -1 && !resolvedAbsolute; i--) {
          var path = (i >= 0) ? arguments[i] : FS.cwd();
          // Skip empty and invalid entries
          if (typeof path !== 'string') {
            throw new TypeError('Arguments to path.resolve must be strings');
          } else if (!path) {
            continue;
          }
          resolvedPath = path + '/' + resolvedPath;
          resolvedAbsolute = path.charAt(0) === '/';
        }
        // At this point the path should be resolved to a full absolute path, but
        // handle relative paths to be safe (might happen when process.cwd() fails)
        resolvedPath = PATH.normalizeArray(resolvedPath.split('/').filter(function(p) {
          return !!p;
        }), !resolvedAbsolute).join('/');
        return ((resolvedAbsolute ? '/' : '') + resolvedPath) || '.';
      },relative:function (from, to) {
        from = PATH.resolve(from).substr(1);
        to = PATH.resolve(to).substr(1);
        function trim(arr) {
          var start = 0;
          for (; start < arr.length; start++) {
            if (arr[start] !== '') break;
          }
          var end = arr.length - 1;
          for (; end >= 0; end--) {
            if (arr[end] !== '') break;
          }
          if (start > end) return [];
          return arr.slice(start, end - start + 1);
        }
        var fromParts = trim(from.split('/'));
        var toParts = trim(to.split('/'));
        var length = Math.min(fromParts.length, toParts.length);
        var samePartsLength = length;
        for (var i = 0; i < length; i++) {
          if (fromParts[i] !== toParts[i]) {
            samePartsLength = i;
            break;
          }
        }
        var outputParts = [];
        for (var i = samePartsLength; i < fromParts.length; i++) {
          outputParts.push('..');
        }
        outputParts = outputParts.concat(toParts.slice(samePartsLength));
        return outputParts.join('/');
      }};
  
  var TTY={ttys:[],init:function () {
        // https://github.com/kripken/emscripten/pull/1555
        // if (ENVIRONMENT_IS_NODE) {
        //   // currently, FS.init does not distinguish if process.stdin is a file or TTY
        //   // device, it always assumes it's a TTY device. because of this, we're forcing
        //   // process.stdin to UTF8 encoding to at least make stdin reading compatible
        //   // with text files until FS.init can be refactored.
        //   process['stdin']['setEncoding']('utf8');
        // }
      },shutdown:function () {
        // https://github.com/kripken/emscripten/pull/1555
        // if (ENVIRONMENT_IS_NODE) {
        //   // inolen: any idea as to why node -e 'process.stdin.read()' wouldn't exit immediately (with process.stdin being a tty)?
        //   // isaacs: because now it's reading from the stream, you've expressed interest in it, so that read() kicks off a _read() which creates a ReadReq operation
        //   // inolen: I thought read() in that case was a synchronous operation that just grabbed some amount of buffered data if it exists?
        //   // isaacs: it is. but it also triggers a _read() call, which calls readStart() on the handle
        //   // isaacs: do process.stdin.pause() and i'd think it'd probably close the pending call
        //   process['stdin']['pause']();
        // }
      },register:function (dev, ops) {
        TTY.ttys[dev] = { input: [], output: [], ops: ops };
        FS.registerDevice(dev, TTY.stream_ops);
      },stream_ops:{open:function (stream) {
          var tty = TTY.ttys[stream.node.rdev];
          if (!tty) {
            throw new FS.ErrnoError(ERRNO_CODES.ENODEV);
          }
          stream.tty = tty;
          stream.seekable = false;
        },close:function (stream) {
          // flush any pending line data
          if (stream.tty.output.length) {
            stream.tty.ops.put_char(stream.tty, 10);
          }
        },read:function (stream, buffer, offset, length, pos /* ignored */) {
          if (!stream.tty || !stream.tty.ops.get_char) {
            throw new FS.ErrnoError(ERRNO_CODES.ENXIO);
          }
          var bytesRead = 0;
          for (var i = 0; i < length; i++) {
            var result;
            try {
              result = stream.tty.ops.get_char(stream.tty);
            } catch (e) {
              throw new FS.ErrnoError(ERRNO_CODES.EIO);
            }
            if (result === undefined && bytesRead === 0) {
              throw new FS.ErrnoError(ERRNO_CODES.EAGAIN);
            }
            if (result === null || result === undefined) break;
            bytesRead++;
            buffer[offset+i] = result;
          }
          if (bytesRead) {
            stream.node.timestamp = Date.now();
          }
          return bytesRead;
        },write:function (stream, buffer, offset, length, pos) {
          if (!stream.tty || !stream.tty.ops.put_char) {
            throw new FS.ErrnoError(ERRNO_CODES.ENXIO);
          }
          for (var i = 0; i < length; i++) {
            try {
              stream.tty.ops.put_char(stream.tty, buffer[offset+i]);
            } catch (e) {
              throw new FS.ErrnoError(ERRNO_CODES.EIO);
            }
          }
          if (length) {
            stream.node.timestamp = Date.now();
          }
          return i;
        }},default_tty_ops:{get_char:function (tty) {
          if (!tty.input.length) {
            var result = null;
            if (ENVIRONMENT_IS_NODE) {
              result = process['stdin']['read']();
              if (!result) {
                if (process['stdin']['_readableState'] && process['stdin']['_readableState']['ended']) {
                  return null;  // EOF
                }
                return undefined;  // no data available
              }
            } else if (typeof window != 'undefined' &&
              typeof window.prompt == 'function') {
              // Browser.
              result = window.prompt('Input: ');  // returns null on cancel
              if (result !== null) {
                result += '\n';
              }
            } else if (typeof readline == 'function') {
              // Command line.
              result = readline();
              if (result !== null) {
                result += '\n';
              }
            }
            if (!result) {
              return null;
            }
            tty.input = intArrayFromString(result, true);
          }
          return tty.input.shift();
        },put_char:function (tty, val) {
          if (val === null || val === 10) {
            Module['print'](tty.output.join(''));
            tty.output = [];
          } else {
            tty.output.push(TTY.utf8.processCChar(val));
          }
        }},default_tty1_ops:{put_char:function (tty, val) {
          if (val === null || val === 10) {
            Module['printErr'](tty.output.join(''));
            tty.output = [];
          } else {
            tty.output.push(TTY.utf8.processCChar(val));
          }
        }}};
  
  var MEMFS={ops_table:null,CONTENT_OWNING:1,CONTENT_FLEXIBLE:2,CONTENT_FIXED:3,mount:function (mount) {
        return MEMFS.createNode(null, '/', 16384 | 511 /* 0777 */, 0);
      },createNode:function (parent, name, mode, dev) {
        if (FS.isBlkdev(mode) || FS.isFIFO(mode)) {
          // no supported
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        if (!MEMFS.ops_table) {
          MEMFS.ops_table = {
            dir: {
              node: {
                getattr: MEMFS.node_ops.getattr,
                setattr: MEMFS.node_ops.setattr,
                lookup: MEMFS.node_ops.lookup,
                mknod: MEMFS.node_ops.mknod,
                rename: MEMFS.node_ops.rename,
                unlink: MEMFS.node_ops.unlink,
                rmdir: MEMFS.node_ops.rmdir,
                readdir: MEMFS.node_ops.readdir,
                symlink: MEMFS.node_ops.symlink
              },
              stream: {
                llseek: MEMFS.stream_ops.llseek
              }
            },
            file: {
              node: {
                getattr: MEMFS.node_ops.getattr,
                setattr: MEMFS.node_ops.setattr
              },
              stream: {
                llseek: MEMFS.stream_ops.llseek,
                read: MEMFS.stream_ops.read,
                write: MEMFS.stream_ops.write,
                allocate: MEMFS.stream_ops.allocate,
                mmap: MEMFS.stream_ops.mmap
              }
            },
            link: {
              node: {
                getattr: MEMFS.node_ops.getattr,
                setattr: MEMFS.node_ops.setattr,
                readlink: MEMFS.node_ops.readlink
              },
              stream: {}
            },
            chrdev: {
              node: {
                getattr: MEMFS.node_ops.getattr,
                setattr: MEMFS.node_ops.setattr
              },
              stream: FS.chrdev_stream_ops
            },
          };
        }
        var node = FS.createNode(parent, name, mode, dev);
        if (FS.isDir(node.mode)) {
          node.node_ops = MEMFS.ops_table.dir.node;
          node.stream_ops = MEMFS.ops_table.dir.stream;
          node.contents = {};
        } else if (FS.isFile(node.mode)) {
          node.node_ops = MEMFS.ops_table.file.node;
          node.stream_ops = MEMFS.ops_table.file.stream;
          node.contents = [];
          node.contentMode = MEMFS.CONTENT_FLEXIBLE;
        } else if (FS.isLink(node.mode)) {
          node.node_ops = MEMFS.ops_table.link.node;
          node.stream_ops = MEMFS.ops_table.link.stream;
        } else if (FS.isChrdev(node.mode)) {
          node.node_ops = MEMFS.ops_table.chrdev.node;
          node.stream_ops = MEMFS.ops_table.chrdev.stream;
        }
        node.timestamp = Date.now();
        // add the new node to the parent
        if (parent) {
          parent.contents[name] = node;
        }
        return node;
      },ensureFlexible:function (node) {
        if (node.contentMode !== MEMFS.CONTENT_FLEXIBLE) {
          var contents = node.contents;
          node.contents = Array.prototype.slice.call(contents);
          node.contentMode = MEMFS.CONTENT_FLEXIBLE;
        }
      },node_ops:{getattr:function (node) {
          var attr = {};
          // device numbers reuse inode numbers.
          attr.dev = FS.isChrdev(node.mode) ? node.id : 1;
          attr.ino = node.id;
          attr.mode = node.mode;
          attr.nlink = 1;
          attr.uid = 0;
          attr.gid = 0;
          attr.rdev = node.rdev;
          if (FS.isDir(node.mode)) {
            attr.size = 4096;
          } else if (FS.isFile(node.mode)) {
            attr.size = node.contents.length;
          } else if (FS.isLink(node.mode)) {
            attr.size = node.link.length;
          } else {
            attr.size = 0;
          }
          attr.atime = new Date(node.timestamp);
          attr.mtime = new Date(node.timestamp);
          attr.ctime = new Date(node.timestamp);
          // NOTE: In our implementation, st_blocks = Math.ceil(st_size/st_blksize),
          //       but this is not required by the standard.
          attr.blksize = 4096;
          attr.blocks = Math.ceil(attr.size / attr.blksize);
          return attr;
        },setattr:function (node, attr) {
          if (attr.mode !== undefined) {
            node.mode = attr.mode;
          }
          if (attr.timestamp !== undefined) {
            node.timestamp = attr.timestamp;
          }
          if (attr.size !== undefined) {
            MEMFS.ensureFlexible(node);
            var contents = node.contents;
            if (attr.size < contents.length) contents.length = attr.size;
            else while (attr.size > contents.length) contents.push(0);
          }
        },lookup:function (parent, name) {
          throw FS.genericErrors[ERRNO_CODES.ENOENT];
        },mknod:function (parent, name, mode, dev) {
          return MEMFS.createNode(parent, name, mode, dev);
        },rename:function (old_node, new_dir, new_name) {
          // if we're overwriting a directory at new_name, make sure it's empty.
          if (FS.isDir(old_node.mode)) {
            var new_node;
            try {
              new_node = FS.lookupNode(new_dir, new_name);
            } catch (e) {
            }
            if (new_node) {
              for (var i in new_node.contents) {
                throw new FS.ErrnoError(ERRNO_CODES.ENOTEMPTY);
              }
            }
          }
          // do the internal rewiring
          delete old_node.parent.contents[old_node.name];
          old_node.name = new_name;
          new_dir.contents[new_name] = old_node;
          old_node.parent = new_dir;
        },unlink:function (parent, name) {
          delete parent.contents[name];
        },rmdir:function (parent, name) {
          var node = FS.lookupNode(parent, name);
          for (var i in node.contents) {
            throw new FS.ErrnoError(ERRNO_CODES.ENOTEMPTY);
          }
          delete parent.contents[name];
        },readdir:function (node) {
          var entries = ['.', '..']
          for (var key in node.contents) {
            if (!node.contents.hasOwnProperty(key)) {
              continue;
            }
            entries.push(key);
          }
          return entries;
        },symlink:function (parent, newname, oldpath) {
          var node = MEMFS.createNode(parent, newname, 511 /* 0777 */ | 40960, 0);
          node.link = oldpath;
          return node;
        },readlink:function (node) {
          if (!FS.isLink(node.mode)) {
            throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
          }
          return node.link;
        }},stream_ops:{read:function (stream, buffer, offset, length, position) {
          var contents = stream.node.contents;
          if (position >= contents.length)
            return 0;
          var size = Math.min(contents.length - position, length);
          assert(size >= 0);
          if (size > 8 && contents.subarray) { // non-trivial, and typed array
            buffer.set(contents.subarray(position, position + size), offset);
          } else
          {
            for (var i = 0; i < size; i++) {
              buffer[offset + i] = contents[position + i];
            }
          }
          return size;
        },write:function (stream, buffer, offset, length, position, canOwn) {
          var node = stream.node;
          node.timestamp = Date.now();
          var contents = node.contents;
          if (length && contents.length === 0 && position === 0 && buffer.subarray) {
            // just replace it with the new data
            if (canOwn && offset === 0) {
              node.contents = buffer; // this could be a subarray of Emscripten HEAP, or allocated from some other source.
              node.contentMode = (buffer.buffer === HEAP8.buffer) ? MEMFS.CONTENT_OWNING : MEMFS.CONTENT_FIXED;
            } else {
              node.contents = new Uint8Array(buffer.subarray(offset, offset+length));
              node.contentMode = MEMFS.CONTENT_FIXED;
            }
            return length;
          }
          MEMFS.ensureFlexible(node);
          var contents = node.contents;
          while (contents.length < position) contents.push(0);
          for (var i = 0; i < length; i++) {
            contents[position + i] = buffer[offset + i];
          }
          return length;
        },llseek:function (stream, offset, whence) {
          var position = offset;
          if (whence === 1) {  // SEEK_CUR.
            position += stream.position;
          } else if (whence === 2) {  // SEEK_END.
            if (FS.isFile(stream.node.mode)) {
              position += stream.node.contents.length;
            }
          }
          if (position < 0) {
            throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
          }
          stream.ungotten = [];
          stream.position = position;
          return position;
        },allocate:function (stream, offset, length) {
          MEMFS.ensureFlexible(stream.node);
          var contents = stream.node.contents;
          var limit = offset + length;
          while (limit > contents.length) contents.push(0);
        },mmap:function (stream, buffer, offset, length, position, prot, flags) {
          if (!FS.isFile(stream.node.mode)) {
            throw new FS.ErrnoError(ERRNO_CODES.ENODEV);
          }
          var ptr;
          var allocated;
          var contents = stream.node.contents;
          // Only make a new copy when MAP_PRIVATE is specified.
          if ( !(flags & 2) &&
                (contents.buffer === buffer || contents.buffer === buffer.buffer) ) {
            // We can't emulate MAP_SHARED when the file is not backed by the buffer
            // we're mapping to (e.g. the HEAP buffer).
            allocated = false;
            ptr = contents.byteOffset;
          } else {
            // Try to avoid unnecessary slices.
            if (position > 0 || position + length < contents.length) {
              if (contents.subarray) {
                contents = contents.subarray(position, position + length);
              } else {
                contents = Array.prototype.slice.call(contents, position, position + length);
              }
            }
            allocated = true;
            ptr = _malloc(length);
            if (!ptr) {
              throw new FS.ErrnoError(ERRNO_CODES.ENOMEM);
            }
            buffer.set(contents, ptr);
          }
          return { ptr: ptr, allocated: allocated };
        }}};
  
  var IDBFS={dbs:{},indexedDB:function () {
        return window.indexedDB || window.mozIndexedDB || window.webkitIndexedDB || window.msIndexedDB;
      },DB_VERSION:21,DB_STORE_NAME:"FILE_DATA",mount:function (mount) {
        // reuse all of the core MEMFS functionality
        return MEMFS.mount.apply(null, arguments);
      },syncfs:function (mount, populate, callback) {
        IDBFS.getLocalSet(mount, function(err, local) {
          if (err) return callback(err);
  
          IDBFS.getRemoteSet(mount, function(err, remote) {
            if (err) return callback(err);
  
            var src = populate ? remote : local;
            var dst = populate ? local : remote;
  
            IDBFS.reconcile(src, dst, callback);
          });
        });
      },getDB:function (name, callback) {
        // check the cache first
        var db = IDBFS.dbs[name];
        if (db) {
          return callback(null, db);
        }
  
        var req;
        try {
          req = IDBFS.indexedDB().open(name, IDBFS.DB_VERSION);
        } catch (e) {
          return callback(e);
        }
        req.onupgradeneeded = function(e) {
          var db = e.target.result;
          var transaction = e.target.transaction;
  
          var fileStore;
  
          if (db.objectStoreNames.contains(IDBFS.DB_STORE_NAME)) {
            fileStore = transaction.objectStore(IDBFS.DB_STORE_NAME);
          } else {
            fileStore = db.createObjectStore(IDBFS.DB_STORE_NAME);
          }
  
          fileStore.createIndex('timestamp', 'timestamp', { unique: false });
        };
        req.onsuccess = function() {
          db = req.result;
  
          // add to the cache
          IDBFS.dbs[name] = db;
          callback(null, db);
        };
        req.onerror = function() {
          callback(this.error);
        };
      },getLocalSet:function (mount, callback) {
        var entries = {};
  
        function isRealDir(p) {
          return p !== '.' && p !== '..';
        };
        function toAbsolute(root) {
          return function(p) {
            return PATH.join2(root, p);
          }
        };
  
        var check = FS.readdir(mount.mountpoint).filter(isRealDir).map(toAbsolute(mount.mountpoint));
  
        while (check.length) {
          var path = check.pop();
          var stat;
  
          try {
            stat = FS.stat(path);
          } catch (e) {
            return callback(e);
          }
  
          if (FS.isDir(stat.mode)) {
            check.push.apply(check, FS.readdir(path).filter(isRealDir).map(toAbsolute(path)));
          }
  
          entries[path] = { timestamp: stat.mtime };
        }
  
        return callback(null, { type: 'local', entries: entries });
      },getRemoteSet:function (mount, callback) {
        var entries = {};
  
        IDBFS.getDB(mount.mountpoint, function(err, db) {
          if (err) return callback(err);
  
          var transaction = db.transaction([IDBFS.DB_STORE_NAME], 'readonly');
          transaction.onerror = function() { callback(this.error); };
  
          var store = transaction.objectStore(IDBFS.DB_STORE_NAME);
          var index = store.index('timestamp');
  
          index.openKeyCursor().onsuccess = function(event) {
            var cursor = event.target.result;
  
            if (!cursor) {
              return callback(null, { type: 'remote', db: db, entries: entries });
            }
  
            entries[cursor.primaryKey] = { timestamp: cursor.key };
  
            cursor.continue();
          };
        });
      },loadLocalEntry:function (path, callback) {
        var stat, node;
  
        try {
          var lookup = FS.lookupPath(path);
          node = lookup.node;
          stat = FS.stat(path);
        } catch (e) {
          return callback(e);
        }
  
        if (FS.isDir(stat.mode)) {
          return callback(null, { timestamp: stat.mtime, mode: stat.mode });
        } else if (FS.isFile(stat.mode)) {
          return callback(null, { timestamp: stat.mtime, mode: stat.mode, contents: node.contents });
        } else {
          return callback(new Error('node type not supported'));
        }
      },storeLocalEntry:function (path, entry, callback) {
        try {
          if (FS.isDir(entry.mode)) {
            FS.mkdir(path, entry.mode);
          } else if (FS.isFile(entry.mode)) {
            FS.writeFile(path, entry.contents, { encoding: 'binary', canOwn: true });
          } else {
            return callback(new Error('node type not supported'));
          }
  
          FS.utime(path, entry.timestamp, entry.timestamp);
        } catch (e) {
          return callback(e);
        }
  
        callback(null);
      },removeLocalEntry:function (path, callback) {
        try {
          var lookup = FS.lookupPath(path);
          var stat = FS.stat(path);
  
          if (FS.isDir(stat.mode)) {
            FS.rmdir(path);
          } else if (FS.isFile(stat.mode)) {
            FS.unlink(path);
          }
        } catch (e) {
          return callback(e);
        }
  
        callback(null);
      },loadRemoteEntry:function (store, path, callback) {
        var req = store.get(path);
        req.onsuccess = function(event) { callback(null, event.target.result); };
        req.onerror = function() { callback(this.error); };
      },storeRemoteEntry:function (store, path, entry, callback) {
        var req = store.put(entry, path);
        req.onsuccess = function() { callback(null); };
        req.onerror = function() { callback(this.error); };
      },removeRemoteEntry:function (store, path, callback) {
        var req = store.delete(path);
        req.onsuccess = function() { callback(null); };
        req.onerror = function() { callback(this.error); };
      },reconcile:function (src, dst, callback) {
        var total = 0;
  
        var create = [];
        Object.keys(src.entries).forEach(function (key) {
          var e = src.entries[key];
          var e2 = dst.entries[key];
          if (!e2 || e.timestamp > e2.timestamp) {
            create.push(key);
            total++;
          }
        });
  
        var remove = [];
        Object.keys(dst.entries).forEach(function (key) {
          var e = dst.entries[key];
          var e2 = src.entries[key];
          if (!e2) {
            remove.push(key);
            total++;
          }
        });
  
        if (!total) {
          return callback(null);
        }
  
        var errored = false;
        var completed = 0;
        var db = src.type === 'remote' ? src.db : dst.db;
        var transaction = db.transaction([IDBFS.DB_STORE_NAME], 'readwrite');
        var store = transaction.objectStore(IDBFS.DB_STORE_NAME);
  
        function done(err) {
          if (err) {
            if (!done.errored) {
              done.errored = true;
              return callback(err);
            }
            return;
          }
          if (++completed >= total) {
            return callback(null);
          }
        };
  
        transaction.onerror = function() { done(this.error); };
  
        // sort paths in ascending order so directory entries are created
        // before the files inside them
        create.sort().forEach(function (path) {
          if (dst.type === 'local') {
            IDBFS.loadRemoteEntry(store, path, function (err, entry) {
              if (err) return done(err);
              IDBFS.storeLocalEntry(path, entry, done);
            });
          } else {
            IDBFS.loadLocalEntry(path, function (err, entry) {
              if (err) return done(err);
              IDBFS.storeRemoteEntry(store, path, entry, done);
            });
          }
        });
  
        // sort paths in descending order so files are deleted before their
        // parent directories
        remove.sort().reverse().forEach(function(path) {
          if (dst.type === 'local') {
            IDBFS.removeLocalEntry(path, done);
          } else {
            IDBFS.removeRemoteEntry(store, path, done);
          }
        });
      }};
  
  var NODEFS={isWindows:false,staticInit:function () {
        NODEFS.isWindows = !!process.platform.match(/^win/);
      },mount:function (mount) {
        assert(ENVIRONMENT_IS_NODE);
        return NODEFS.createNode(null, '/', NODEFS.getMode(mount.opts.root), 0);
      },createNode:function (parent, name, mode, dev) {
        if (!FS.isDir(mode) && !FS.isFile(mode) && !FS.isLink(mode)) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        var node = FS.createNode(parent, name, mode);
        node.node_ops = NODEFS.node_ops;
        node.stream_ops = NODEFS.stream_ops;
        return node;
      },getMode:function (path) {
        var stat;
        try {
          stat = fs.lstatSync(path);
          if (NODEFS.isWindows) {
            // On Windows, directories return permission bits 'rw-rw-rw-', even though they have 'rwxrwxrwx', so 
            // propagate write bits to execute bits.
            stat.mode = stat.mode | ((stat.mode & 146) >> 1);
          }
        } catch (e) {
          if (!e.code) throw e;
          throw new FS.ErrnoError(ERRNO_CODES[e.code]);
        }
        return stat.mode;
      },realPath:function (node) {
        var parts = [];
        while (node.parent !== node) {
          parts.push(node.name);
          node = node.parent;
        }
        parts.push(node.mount.opts.root);
        parts.reverse();
        return PATH.join.apply(null, parts);
      },flagsToPermissionStringMap:{0:"r",1:"r+",2:"r+",64:"r",65:"r+",66:"r+",129:"rx+",193:"rx+",514:"w+",577:"w",578:"w+",705:"wx",706:"wx+",1024:"a",1025:"a",1026:"a+",1089:"a",1090:"a+",1153:"ax",1154:"ax+",1217:"ax",1218:"ax+",4096:"rs",4098:"rs+"},flagsToPermissionString:function (flags) {
        if (flags in NODEFS.flagsToPermissionStringMap) {
          return NODEFS.flagsToPermissionStringMap[flags];
        } else {
          return flags;
        }
      },node_ops:{getattr:function (node) {
          var path = NODEFS.realPath(node);
          var stat;
          try {
            stat = fs.lstatSync(path);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
          // node.js v0.10.20 doesn't report blksize and blocks on Windows. Fake them with default blksize of 4096.
          // See http://support.microsoft.com/kb/140365
          if (NODEFS.isWindows && !stat.blksize) {
            stat.blksize = 4096;
          }
          if (NODEFS.isWindows && !stat.blocks) {
            stat.blocks = (stat.size+stat.blksize-1)/stat.blksize|0;
          }
          return {
            dev: stat.dev,
            ino: stat.ino,
            mode: stat.mode,
            nlink: stat.nlink,
            uid: stat.uid,
            gid: stat.gid,
            rdev: stat.rdev,
            size: stat.size,
            atime: stat.atime,
            mtime: stat.mtime,
            ctime: stat.ctime,
            blksize: stat.blksize,
            blocks: stat.blocks
          };
        },setattr:function (node, attr) {
          var path = NODEFS.realPath(node);
          try {
            if (attr.mode !== undefined) {
              fs.chmodSync(path, attr.mode);
              // update the common node structure mode as well
              node.mode = attr.mode;
            }
            if (attr.timestamp !== undefined) {
              var date = new Date(attr.timestamp);
              fs.utimesSync(path, date, date);
            }
            if (attr.size !== undefined) {
              fs.truncateSync(path, attr.size);
            }
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },lookup:function (parent, name) {
          var path = PATH.join2(NODEFS.realPath(parent), name);
          var mode = NODEFS.getMode(path);
          return NODEFS.createNode(parent, name, mode);
        },mknod:function (parent, name, mode, dev) {
          var node = NODEFS.createNode(parent, name, mode, dev);
          // create the backing node for this in the fs root as well
          var path = NODEFS.realPath(node);
          try {
            if (FS.isDir(node.mode)) {
              fs.mkdirSync(path, node.mode);
            } else {
              fs.writeFileSync(path, '', { mode: node.mode });
            }
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
          return node;
        },rename:function (oldNode, newDir, newName) {
          var oldPath = NODEFS.realPath(oldNode);
          var newPath = PATH.join2(NODEFS.realPath(newDir), newName);
          try {
            fs.renameSync(oldPath, newPath);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },unlink:function (parent, name) {
          var path = PATH.join2(NODEFS.realPath(parent), name);
          try {
            fs.unlinkSync(path);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },rmdir:function (parent, name) {
          var path = PATH.join2(NODEFS.realPath(parent), name);
          try {
            fs.rmdirSync(path);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },readdir:function (node) {
          var path = NODEFS.realPath(node);
          try {
            return fs.readdirSync(path);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },symlink:function (parent, newName, oldPath) {
          var newPath = PATH.join2(NODEFS.realPath(parent), newName);
          try {
            fs.symlinkSync(oldPath, newPath);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },readlink:function (node) {
          var path = NODEFS.realPath(node);
          try {
            return fs.readlinkSync(path);
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        }},stream_ops:{open:function (stream) {
          var path = NODEFS.realPath(stream.node);
          try {
            if (FS.isFile(stream.node.mode)) {
              stream.nfd = fs.openSync(path, NODEFS.flagsToPermissionString(stream.flags));
            }
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },close:function (stream) {
          try {
            if (FS.isFile(stream.node.mode) && stream.nfd) {
              fs.closeSync(stream.nfd);
            }
          } catch (e) {
            if (!e.code) throw e;
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
        },read:function (stream, buffer, offset, length, position) {
          // FIXME this is terrible.
          var nbuffer = new Buffer(length);
          var res;
          try {
            res = fs.readSync(stream.nfd, nbuffer, 0, length, position);
          } catch (e) {
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
          if (res > 0) {
            for (var i = 0; i < res; i++) {
              buffer[offset + i] = nbuffer[i];
            }
          }
          return res;
        },write:function (stream, buffer, offset, length, position) {
          // FIXME this is terrible.
          var nbuffer = new Buffer(buffer.subarray(offset, offset + length));
          var res;
          try {
            res = fs.writeSync(stream.nfd, nbuffer, 0, length, position);
          } catch (e) {
            throw new FS.ErrnoError(ERRNO_CODES[e.code]);
          }
          return res;
        },llseek:function (stream, offset, whence) {
          var position = offset;
          if (whence === 1) {  // SEEK_CUR.
            position += stream.position;
          } else if (whence === 2) {  // SEEK_END.
            if (FS.isFile(stream.node.mode)) {
              try {
                var stat = fs.fstatSync(stream.nfd);
                position += stat.size;
              } catch (e) {
                throw new FS.ErrnoError(ERRNO_CODES[e.code]);
              }
            }
          }
  
          if (position < 0) {
            throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
          }
  
          stream.position = position;
          return position;
        }}};
  
  var _stdin=allocate(1, "i32*", ALLOC_STATIC);
  
  var _stdout=allocate(1, "i32*", ALLOC_STATIC);
  
  var _stderr=allocate(1, "i32*", ALLOC_STATIC);var FS={root:null,mounts:[],devices:[null],streams:[],nextInode:1,nameTable:null,currentPath:"/",initialized:false,ignorePermissions:true,ErrnoError:null,genericErrors:{},handleFSError:function (e) {
        if (!(e instanceof FS.ErrnoError)) throw e + ' : ' + stackTrace();
        return ___setErrNo(e.errno);
      },lookupPath:function (path, opts) {
        path = PATH.resolve(FS.cwd(), path);
        opts = opts || {};
  
        var defaults = {
          follow_mount: true,
          recurse_count: 0
        };
        for (var key in defaults) {
          if (opts[key] === undefined) {
            opts[key] = defaults[key];
          }
        }
  
        if (opts.recurse_count > 8) {  // max recursive lookup of 8
          throw new FS.ErrnoError(ERRNO_CODES.ELOOP);
        }
  
        // split the path
        var parts = PATH.normalizeArray(path.split('/').filter(function(p) {
          return !!p;
        }), false);
  
        // start at the root
        var current = FS.root;
        var current_path = '/';
  
        for (var i = 0; i < parts.length; i++) {
          var islast = (i === parts.length-1);
          if (islast && opts.parent) {
            // stop resolving
            break;
          }
  
          current = FS.lookupNode(current, parts[i]);
          current_path = PATH.join2(current_path, parts[i]);
  
          // jump to the mount's root node if this is a mountpoint
          if (FS.isMountpoint(current)) {
            if (!islast || (islast && opts.follow_mount)) {
              current = current.mounted.root;
            }
          }
  
          // by default, lookupPath will not follow a symlink if it is the final path component.
          // setting opts.follow = true will override this behavior.
          if (!islast || opts.follow) {
            var count = 0;
            while (FS.isLink(current.mode)) {
              var link = FS.readlink(current_path);
              current_path = PATH.resolve(PATH.dirname(current_path), link);
              
              var lookup = FS.lookupPath(current_path, { recurse_count: opts.recurse_count });
              current = lookup.node;
  
              if (count++ > 40) {  // limit max consecutive symlinks to 40 (SYMLOOP_MAX).
                throw new FS.ErrnoError(ERRNO_CODES.ELOOP);
              }
            }
          }
        }
  
        return { path: current_path, node: current };
      },getPath:function (node) {
        var path;
        while (true) {
          if (FS.isRoot(node)) {
            var mount = node.mount.mountpoint;
            if (!path) return mount;
            return mount[mount.length-1] !== '/' ? mount + '/' + path : mount + path;
          }
          path = path ? node.name + '/' + path : node.name;
          node = node.parent;
        }
      },hashName:function (parentid, name) {
        var hash = 0;
  
  
        for (var i = 0; i < name.length; i++) {
          hash = ((hash << 5) - hash + name.charCodeAt(i)) | 0;
        }
        return ((parentid + hash) >>> 0) % FS.nameTable.length;
      },hashAddNode:function (node) {
        var hash = FS.hashName(node.parent.id, node.name);
        node.name_next = FS.nameTable[hash];
        FS.nameTable[hash] = node;
      },hashRemoveNode:function (node) {
        var hash = FS.hashName(node.parent.id, node.name);
        if (FS.nameTable[hash] === node) {
          FS.nameTable[hash] = node.name_next;
        } else {
          var current = FS.nameTable[hash];
          while (current) {
            if (current.name_next === node) {
              current.name_next = node.name_next;
              break;
            }
            current = current.name_next;
          }
        }
      },lookupNode:function (parent, name) {
        var err = FS.mayLookup(parent);
        if (err) {
          throw new FS.ErrnoError(err);
        }
        var hash = FS.hashName(parent.id, name);
        for (var node = FS.nameTable[hash]; node; node = node.name_next) {
          var nodeName = node.name;
          if (node.parent.id === parent.id && nodeName === name) {
            return node;
          }
        }
        // if we failed to find it in the cache, call into the VFS
        return FS.lookup(parent, name);
      },createNode:function (parent, name, mode, rdev) {
        if (!FS.FSNode) {
          FS.FSNode = function(parent, name, mode, rdev) {
            if (!parent) {
              parent = this;  // root node sets parent to itself
            }
            this.parent = parent;
            this.mount = parent.mount;
            this.mounted = null;
            this.id = FS.nextInode++;
            this.name = name;
            this.mode = mode;
            this.node_ops = {};
            this.stream_ops = {};
            this.rdev = rdev;
          };
  
          FS.FSNode.prototype = {};
  
          // compatibility
          var readMode = 292 | 73;
          var writeMode = 146;
  
          // NOTE we must use Object.defineProperties instead of individual calls to
          // Object.defineProperty in order to make closure compiler happy
          Object.defineProperties(FS.FSNode.prototype, {
            read: {
              get: function() { return (this.mode & readMode) === readMode; },
              set: function(val) { val ? this.mode |= readMode : this.mode &= ~readMode; }
            },
            write: {
              get: function() { return (this.mode & writeMode) === writeMode; },
              set: function(val) { val ? this.mode |= writeMode : this.mode &= ~writeMode; }
            },
            isFolder: {
              get: function() { return FS.isDir(this.mode); },
            },
            isDevice: {
              get: function() { return FS.isChrdev(this.mode); },
            },
          });
        }
  
        var node = new FS.FSNode(parent, name, mode, rdev);
  
        FS.hashAddNode(node);
  
        return node;
      },destroyNode:function (node) {
        FS.hashRemoveNode(node);
      },isRoot:function (node) {
        return node === node.parent;
      },isMountpoint:function (node) {
        return !!node.mounted;
      },isFile:function (mode) {
        return (mode & 61440) === 32768;
      },isDir:function (mode) {
        return (mode & 61440) === 16384;
      },isLink:function (mode) {
        return (mode & 61440) === 40960;
      },isChrdev:function (mode) {
        return (mode & 61440) === 8192;
      },isBlkdev:function (mode) {
        return (mode & 61440) === 24576;
      },isFIFO:function (mode) {
        return (mode & 61440) === 4096;
      },isSocket:function (mode) {
        return (mode & 49152) === 49152;
      },flagModes:{"r":0,"rs":1052672,"r+":2,"w":577,"wx":705,"xw":705,"w+":578,"wx+":706,"xw+":706,"a":1089,"ax":1217,"xa":1217,"a+":1090,"ax+":1218,"xa+":1218},modeStringToFlags:function (str) {
        var flags = FS.flagModes[str];
        if (typeof flags === 'undefined') {
          throw new Error('Unknown file open mode: ' + str);
        }
        return flags;
      },flagsToPermissionString:function (flag) {
        var accmode = flag & 2097155;
        var perms = ['r', 'w', 'rw'][accmode];
        if ((flag & 512)) {
          perms += 'w';
        }
        return perms;
      },nodePermissions:function (node, perms) {
        if (FS.ignorePermissions) {
          return 0;
        }
        // return 0 if any user, group or owner bits are set.
        if (perms.indexOf('r') !== -1 && !(node.mode & 292)) {
          return ERRNO_CODES.EACCES;
        } else if (perms.indexOf('w') !== -1 && !(node.mode & 146)) {
          return ERRNO_CODES.EACCES;
        } else if (perms.indexOf('x') !== -1 && !(node.mode & 73)) {
          return ERRNO_CODES.EACCES;
        }
        return 0;
      },mayLookup:function (dir) {
        return FS.nodePermissions(dir, 'x');
      },mayCreate:function (dir, name) {
        try {
          var node = FS.lookupNode(dir, name);
          return ERRNO_CODES.EEXIST;
        } catch (e) {
        }
        return FS.nodePermissions(dir, 'wx');
      },mayDelete:function (dir, name, isdir) {
        var node;
        try {
          node = FS.lookupNode(dir, name);
        } catch (e) {
          return e.errno;
        }
        var err = FS.nodePermissions(dir, 'wx');
        if (err) {
          return err;
        }
        if (isdir) {
          if (!FS.isDir(node.mode)) {
            return ERRNO_CODES.ENOTDIR;
          }
          if (FS.isRoot(node) || FS.getPath(node) === FS.cwd()) {
            return ERRNO_CODES.EBUSY;
          }
        } else {
          if (FS.isDir(node.mode)) {
            return ERRNO_CODES.EISDIR;
          }
        }
        return 0;
      },mayOpen:function (node, flags) {
        if (!node) {
          return ERRNO_CODES.ENOENT;
        }
        if (FS.isLink(node.mode)) {
          return ERRNO_CODES.ELOOP;
        } else if (FS.isDir(node.mode)) {
          if ((flags & 2097155) !== 0 ||  // opening for write
              (flags & 512)) {
            return ERRNO_CODES.EISDIR;
          }
        }
        return FS.nodePermissions(node, FS.flagsToPermissionString(flags));
      },MAX_OPEN_FDS:4096,nextfd:function (fd_start, fd_end) {
        fd_start = fd_start || 0;
        fd_end = fd_end || FS.MAX_OPEN_FDS;
        for (var fd = fd_start; fd <= fd_end; fd++) {
          if (!FS.streams[fd]) {
            return fd;
          }
        }
        throw new FS.ErrnoError(ERRNO_CODES.EMFILE);
      },getStream:function (fd) {
        return FS.streams[fd];
      },createStream:function (stream, fd_start, fd_end) {
        if (!FS.FSStream) {
          FS.FSStream = function(){};
          FS.FSStream.prototype = {};
          // compatibility
          Object.defineProperties(FS.FSStream.prototype, {
            object: {
              get: function() { return this.node; },
              set: function(val) { this.node = val; }
            },
            isRead: {
              get: function() { return (this.flags & 2097155) !== 1; }
            },
            isWrite: {
              get: function() { return (this.flags & 2097155) !== 0; }
            },
            isAppend: {
              get: function() { return (this.flags & 1024); }
            }
          });
        }
        // clone it, so we can return an instance of FSStream
        var newStream = new FS.FSStream();
        for (var p in stream) {
          newStream[p] = stream[p];
        }
        stream = newStream;
        var fd = FS.nextfd(fd_start, fd_end);
        stream.fd = fd;
        FS.streams[fd] = stream;
        return stream;
      },closeStream:function (fd) {
        FS.streams[fd] = null;
      },getStreamFromPtr:function (ptr) {
        return FS.streams[ptr - 1];
      },getPtrForStream:function (stream) {
        return stream ? stream.fd + 1 : 0;
      },chrdev_stream_ops:{open:function (stream) {
          var device = FS.getDevice(stream.node.rdev);
          // override node's stream ops with the device's
          stream.stream_ops = device.stream_ops;
          // forward the open call
          if (stream.stream_ops.open) {
            stream.stream_ops.open(stream);
          }
        },llseek:function () {
          throw new FS.ErrnoError(ERRNO_CODES.ESPIPE);
        }},major:function (dev) {
        return ((dev) >> 8);
      },minor:function (dev) {
        return ((dev) & 0xff);
      },makedev:function (ma, mi) {
        return ((ma) << 8 | (mi));
      },registerDevice:function (dev, ops) {
        FS.devices[dev] = { stream_ops: ops };
      },getDevice:function (dev) {
        return FS.devices[dev];
      },getMounts:function (mount) {
        var mounts = [];
        var check = [mount];
  
        while (check.length) {
          var m = check.pop();
  
          mounts.push(m);
  
          check.push.apply(check, m.mounts);
        }
  
        return mounts;
      },syncfs:function (populate, callback) {
        if (typeof(populate) === 'function') {
          callback = populate;
          populate = false;
        }
  
        var mounts = FS.getMounts(FS.root.mount);
        var completed = 0;
  
        function done(err) {
          if (err) {
            if (!done.errored) {
              done.errored = true;
              return callback(err);
            }
            return;
          }
          if (++completed >= mounts.length) {
            callback(null);
          }
        };
  
        // sync all mounts
        mounts.forEach(function (mount) {
          if (!mount.type.syncfs) {
            return done(null);
          }
          mount.type.syncfs(mount, populate, done);
        });
      },mount:function (type, opts, mountpoint) {
        var root = mountpoint === '/';
        var pseudo = !mountpoint;
        var node;
  
        if (root && FS.root) {
          throw new FS.ErrnoError(ERRNO_CODES.EBUSY);
        } else if (!root && !pseudo) {
          var lookup = FS.lookupPath(mountpoint, { follow_mount: false });
  
          mountpoint = lookup.path;  // use the absolute path
          node = lookup.node;
  
          if (FS.isMountpoint(node)) {
            throw new FS.ErrnoError(ERRNO_CODES.EBUSY);
          }
  
          if (!FS.isDir(node.mode)) {
            throw new FS.ErrnoError(ERRNO_CODES.ENOTDIR);
          }
        }
  
        var mount = {
          type: type,
          opts: opts,
          mountpoint: mountpoint,
          mounts: []
        };
  
        // create a root node for the fs
        var mountRoot = type.mount(mount);
        mountRoot.mount = mount;
        mount.root = mountRoot;
  
        if (root) {
          FS.root = mountRoot;
        } else if (node) {
          // set as a mountpoint
          node.mounted = mount;
  
          // add the new mount to the current mount's children
          if (node.mount) {
            node.mount.mounts.push(mount);
          }
        }
  
        return mountRoot;
      },unmount:function (mountpoint) {
        var lookup = FS.lookupPath(mountpoint, { follow_mount: false });
  
        if (!FS.isMountpoint(lookup.node)) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
  
        // destroy the nodes for this mount, and all its child mounts
        var node = lookup.node;
        var mount = node.mounted;
        var mounts = FS.getMounts(mount);
  
        Object.keys(FS.nameTable).forEach(function (hash) {
          var current = FS.nameTable[hash];
  
          while (current) {
            var next = current.name_next;
  
            if (mounts.indexOf(current.mount) !== -1) {
              FS.destroyNode(current);
            }
  
            current = next;
          }
        });
  
        // no longer a mountpoint
        node.mounted = null;
  
        // remove this mount from the child mounts
        var idx = node.mount.mounts.indexOf(mount);
        assert(idx !== -1);
        node.mount.mounts.splice(idx, 1);
      },lookup:function (parent, name) {
        return parent.node_ops.lookup(parent, name);
      },mknod:function (path, mode, dev) {
        var lookup = FS.lookupPath(path, { parent: true });
        var parent = lookup.node;
        var name = PATH.basename(path);
        var err = FS.mayCreate(parent, name);
        if (err) {
          throw new FS.ErrnoError(err);
        }
        if (!parent.node_ops.mknod) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        return parent.node_ops.mknod(parent, name, mode, dev);
      },create:function (path, mode) {
        mode = mode !== undefined ? mode : 438 /* 0666 */;
        mode &= 4095;
        mode |= 32768;
        return FS.mknod(path, mode, 0);
      },mkdir:function (path, mode) {
        mode = mode !== undefined ? mode : 511 /* 0777 */;
        mode &= 511 | 512;
        mode |= 16384;
        return FS.mknod(path, mode, 0);
      },mkdev:function (path, mode, dev) {
        if (typeof(dev) === 'undefined') {
          dev = mode;
          mode = 438 /* 0666 */;
        }
        mode |= 8192;
        return FS.mknod(path, mode, dev);
      },symlink:function (oldpath, newpath) {
        var lookup = FS.lookupPath(newpath, { parent: true });
        var parent = lookup.node;
        var newname = PATH.basename(newpath);
        var err = FS.mayCreate(parent, newname);
        if (err) {
          throw new FS.ErrnoError(err);
        }
        if (!parent.node_ops.symlink) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        return parent.node_ops.symlink(parent, newname, oldpath);
      },rename:function (old_path, new_path) {
        var old_dirname = PATH.dirname(old_path);
        var new_dirname = PATH.dirname(new_path);
        var old_name = PATH.basename(old_path);
        var new_name = PATH.basename(new_path);
        // parents must exist
        var lookup, old_dir, new_dir;
        try {
          lookup = FS.lookupPath(old_path, { parent: true });
          old_dir = lookup.node;
          lookup = FS.lookupPath(new_path, { parent: true });
          new_dir = lookup.node;
        } catch (e) {
          throw new FS.ErrnoError(ERRNO_CODES.EBUSY);
        }
        // need to be part of the same mount
        if (old_dir.mount !== new_dir.mount) {
          throw new FS.ErrnoError(ERRNO_CODES.EXDEV);
        }
        // source must exist
        var old_node = FS.lookupNode(old_dir, old_name);
        // old path should not be an ancestor of the new path
        var relative = PATH.relative(old_path, new_dirname);
        if (relative.charAt(0) !== '.') {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        // new path should not be an ancestor of the old path
        relative = PATH.relative(new_path, old_dirname);
        if (relative.charAt(0) !== '.') {
          throw new FS.ErrnoError(ERRNO_CODES.ENOTEMPTY);
        }
        // see if the new path already exists
        var new_node;
        try {
          new_node = FS.lookupNode(new_dir, new_name);
        } catch (e) {
          // not fatal
        }
        // early out if nothing needs to change
        if (old_node === new_node) {
          return;
        }
        // we'll need to delete the old entry
        var isdir = FS.isDir(old_node.mode);
        var err = FS.mayDelete(old_dir, old_name, isdir);
        if (err) {
          throw new FS.ErrnoError(err);
        }
        // need delete permissions if we'll be overwriting.
        // need create permissions if new doesn't already exist.
        err = new_node ?
          FS.mayDelete(new_dir, new_name, isdir) :
          FS.mayCreate(new_dir, new_name);
        if (err) {
          throw new FS.ErrnoError(err);
        }
        if (!old_dir.node_ops.rename) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        if (FS.isMountpoint(old_node) || (new_node && FS.isMountpoint(new_node))) {
          throw new FS.ErrnoError(ERRNO_CODES.EBUSY);
        }
        // if we are going to change the parent, check write permissions
        if (new_dir !== old_dir) {
          err = FS.nodePermissions(old_dir, 'w');
          if (err) {
            throw new FS.ErrnoError(err);
          }
        }
        // remove the node from the lookup hash
        FS.hashRemoveNode(old_node);
        // do the underlying fs rename
        try {
          old_dir.node_ops.rename(old_node, new_dir, new_name);
        } catch (e) {
          throw e;
        } finally {
          // add the node back to the hash (in case node_ops.rename
          // changed its name)
          FS.hashAddNode(old_node);
        }
      },rmdir:function (path) {
        var lookup = FS.lookupPath(path, { parent: true });
        var parent = lookup.node;
        var name = PATH.basename(path);
        var node = FS.lookupNode(parent, name);
        var err = FS.mayDelete(parent, name, true);
        if (err) {
          throw new FS.ErrnoError(err);
        }
        if (!parent.node_ops.rmdir) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        if (FS.isMountpoint(node)) {
          throw new FS.ErrnoError(ERRNO_CODES.EBUSY);
        }
        parent.node_ops.rmdir(parent, name);
        FS.destroyNode(node);
      },readdir:function (path) {
        var lookup = FS.lookupPath(path, { follow: true });
        var node = lookup.node;
        if (!node.node_ops.readdir) {
          throw new FS.ErrnoError(ERRNO_CODES.ENOTDIR);
        }
        return node.node_ops.readdir(node);
      },unlink:function (path) {
        var lookup = FS.lookupPath(path, { parent: true });
        var parent = lookup.node;
        var name = PATH.basename(path);
        var node = FS.lookupNode(parent, name);
        var err = FS.mayDelete(parent, name, false);
        if (err) {
          // POSIX says unlink should set EPERM, not EISDIR
          if (err === ERRNO_CODES.EISDIR) err = ERRNO_CODES.EPERM;
          throw new FS.ErrnoError(err);
        }
        if (!parent.node_ops.unlink) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        if (FS.isMountpoint(node)) {
          throw new FS.ErrnoError(ERRNO_CODES.EBUSY);
        }
        parent.node_ops.unlink(parent, name);
        FS.destroyNode(node);
      },readlink:function (path) {
        var lookup = FS.lookupPath(path);
        var link = lookup.node;
        if (!link.node_ops.readlink) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        return link.node_ops.readlink(link);
      },stat:function (path, dontFollow) {
        var lookup = FS.lookupPath(path, { follow: !dontFollow });
        var node = lookup.node;
        if (!node.node_ops.getattr) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        return node.node_ops.getattr(node);
      },lstat:function (path) {
        return FS.stat(path, true);
      },chmod:function (path, mode, dontFollow) {
        var node;
        if (typeof path === 'string') {
          var lookup = FS.lookupPath(path, { follow: !dontFollow });
          node = lookup.node;
        } else {
          node = path;
        }
        if (!node.node_ops.setattr) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        node.node_ops.setattr(node, {
          mode: (mode & 4095) | (node.mode & ~4095),
          timestamp: Date.now()
        });
      },lchmod:function (path, mode) {
        FS.chmod(path, mode, true);
      },fchmod:function (fd, mode) {
        var stream = FS.getStream(fd);
        if (!stream) {
          throw new FS.ErrnoError(ERRNO_CODES.EBADF);
        }
        FS.chmod(stream.node, mode);
      },chown:function (path, uid, gid, dontFollow) {
        var node;
        if (typeof path === 'string') {
          var lookup = FS.lookupPath(path, { follow: !dontFollow });
          node = lookup.node;
        } else {
          node = path;
        }
        if (!node.node_ops.setattr) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        node.node_ops.setattr(node, {
          timestamp: Date.now()
          // we ignore the uid / gid for now
        });
      },lchown:function (path, uid, gid) {
        FS.chown(path, uid, gid, true);
      },fchown:function (fd, uid, gid) {
        var stream = FS.getStream(fd);
        if (!stream) {
          throw new FS.ErrnoError(ERRNO_CODES.EBADF);
        }
        FS.chown(stream.node, uid, gid);
      },truncate:function (path, len) {
        if (len < 0) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        var node;
        if (typeof path === 'string') {
          var lookup = FS.lookupPath(path, { follow: true });
          node = lookup.node;
        } else {
          node = path;
        }
        if (!node.node_ops.setattr) {
          throw new FS.ErrnoError(ERRNO_CODES.EPERM);
        }
        if (FS.isDir(node.mode)) {
          throw new FS.ErrnoError(ERRNO_CODES.EISDIR);
        }
        if (!FS.isFile(node.mode)) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        var err = FS.nodePermissions(node, 'w');
        if (err) {
          throw new FS.ErrnoError(err);
        }
        node.node_ops.setattr(node, {
          size: len,
          timestamp: Date.now()
        });
      },ftruncate:function (fd, len) {
        var stream = FS.getStream(fd);
        if (!stream) {
          throw new FS.ErrnoError(ERRNO_CODES.EBADF);
        }
        if ((stream.flags & 2097155) === 0) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        FS.truncate(stream.node, len);
      },utime:function (path, atime, mtime) {
        var lookup = FS.lookupPath(path, { follow: true });
        var node = lookup.node;
        node.node_ops.setattr(node, {
          timestamp: Math.max(atime, mtime)
        });
      },open:function (path, flags, mode, fd_start, fd_end) {
        flags = typeof flags === 'string' ? FS.modeStringToFlags(flags) : flags;
        mode = typeof mode === 'undefined' ? 438 /* 0666 */ : mode;
        if ((flags & 64)) {
          mode = (mode & 4095) | 32768;
        } else {
          mode = 0;
        }
        var node;
        if (typeof path === 'object') {
          node = path;
        } else {
          path = PATH.normalize(path);
          try {
            var lookup = FS.lookupPath(path, {
              follow: !(flags & 131072)
            });
            node = lookup.node;
          } catch (e) {
            // ignore
          }
        }
        // perhaps we need to create the node
        if ((flags & 64)) {
          if (node) {
            // if O_CREAT and O_EXCL are set, error out if the node already exists
            if ((flags & 128)) {
              throw new FS.ErrnoError(ERRNO_CODES.EEXIST);
            }
          } else {
            // node doesn't exist, try to create it
            node = FS.mknod(path, mode, 0);
          }
        }
        if (!node) {
          throw new FS.ErrnoError(ERRNO_CODES.ENOENT);
        }
        // can't truncate a device
        if (FS.isChrdev(node.mode)) {
          flags &= ~512;
        }
        // check permissions
        var err = FS.mayOpen(node, flags);
        if (err) {
          throw new FS.ErrnoError(err);
        }
        // do truncation if necessary
        if ((flags & 512)) {
          FS.truncate(node, 0);
        }
        // we've already handled these, don't pass down to the underlying vfs
        flags &= ~(128 | 512);
  
        // register the stream with the filesystem
        var stream = FS.createStream({
          node: node,
          path: FS.getPath(node),  // we want the absolute path to the node
          flags: flags,
          seekable: true,
          position: 0,
          stream_ops: node.stream_ops,
          // used by the file family libc calls (fopen, fwrite, ferror, etc.)
          ungotten: [],
          error: false
        }, fd_start, fd_end);
        // call the new stream's open function
        if (stream.stream_ops.open) {
          stream.stream_ops.open(stream);
        }
        if (Module['logReadFiles'] && !(flags & 1)) {
          if (!FS.readFiles) FS.readFiles = {};
          if (!(path in FS.readFiles)) {
            FS.readFiles[path] = 1;
            Module['printErr']('read file: ' + path);
          }
        }
        return stream;
      },close:function (stream) {
        try {
          if (stream.stream_ops.close) {
            stream.stream_ops.close(stream);
          }
        } catch (e) {
          throw e;
        } finally {
          FS.closeStream(stream.fd);
        }
      },llseek:function (stream, offset, whence) {
        if (!stream.seekable || !stream.stream_ops.llseek) {
          throw new FS.ErrnoError(ERRNO_CODES.ESPIPE);
        }
        return stream.stream_ops.llseek(stream, offset, whence);
      },read:function (stream, buffer, offset, length, position) {
        if (length < 0 || position < 0) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        if ((stream.flags & 2097155) === 1) {
          throw new FS.ErrnoError(ERRNO_CODES.EBADF);
        }
        if (FS.isDir(stream.node.mode)) {
          throw new FS.ErrnoError(ERRNO_CODES.EISDIR);
        }
        if (!stream.stream_ops.read) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        var seeking = true;
        if (typeof position === 'undefined') {
          position = stream.position;
          seeking = false;
        } else if (!stream.seekable) {
          throw new FS.ErrnoError(ERRNO_CODES.ESPIPE);
        }
        var bytesRead = stream.stream_ops.read(stream, buffer, offset, length, position);
        if (!seeking) stream.position += bytesRead;
        return bytesRead;
      },write:function (stream, buffer, offset, length, position, canOwn) {
        if (length < 0 || position < 0) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        if ((stream.flags & 2097155) === 0) {
          throw new FS.ErrnoError(ERRNO_CODES.EBADF);
        }
        if (FS.isDir(stream.node.mode)) {
          throw new FS.ErrnoError(ERRNO_CODES.EISDIR);
        }
        if (!stream.stream_ops.write) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        var seeking = true;
        if (typeof position === 'undefined') {
          position = stream.position;
          seeking = false;
        } else if (!stream.seekable) {
          throw new FS.ErrnoError(ERRNO_CODES.ESPIPE);
        }
        if (stream.flags & 1024) {
          // seek to the end before writing in append mode
          FS.llseek(stream, 0, 2);
        }
        var bytesWritten = stream.stream_ops.write(stream, buffer, offset, length, position, canOwn);
        if (!seeking) stream.position += bytesWritten;
        return bytesWritten;
      },allocate:function (stream, offset, length) {
        if (offset < 0 || length <= 0) {
          throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
        }
        if ((stream.flags & 2097155) === 0) {
          throw new FS.ErrnoError(ERRNO_CODES.EBADF);
        }
        if (!FS.isFile(stream.node.mode) && !FS.isDir(node.mode)) {
          throw new FS.ErrnoError(ERRNO_CODES.ENODEV);
        }
        if (!stream.stream_ops.allocate) {
          throw new FS.ErrnoError(ERRNO_CODES.EOPNOTSUPP);
        }
        stream.stream_ops.allocate(stream, offset, length);
      },mmap:function (stream, buffer, offset, length, position, prot, flags) {
        // TODO if PROT is PROT_WRITE, make sure we have write access
        if ((stream.flags & 2097155) === 1) {
          throw new FS.ErrnoError(ERRNO_CODES.EACCES);
        }
        if (!stream.stream_ops.mmap) {
          throw new FS.ErrnoError(ERRNO_CODES.ENODEV);
        }
        return stream.stream_ops.mmap(stream, buffer, offset, length, position, prot, flags);
      },ioctl:function (stream, cmd, arg) {
        if (!stream.stream_ops.ioctl) {
          throw new FS.ErrnoError(ERRNO_CODES.ENOTTY);
        }
        return stream.stream_ops.ioctl(stream, cmd, arg);
      },readFile:function (path, opts) {
        opts = opts || {};
        opts.flags = opts.flags || 'r';
        opts.encoding = opts.encoding || 'binary';
        if (opts.encoding !== 'utf8' && opts.encoding !== 'binary') {
          throw new Error('Invalid encoding type "' + opts.encoding + '"');
        }
        var ret;
        var stream = FS.open(path, opts.flags);
        var stat = FS.stat(path);
        var length = stat.size;
        var buf = new Uint8Array(length);
        FS.read(stream, buf, 0, length, 0);
        if (opts.encoding === 'utf8') {
          ret = '';
          var utf8 = new Runtime.UTF8Processor();
          for (var i = 0; i < length; i++) {
            ret += utf8.processCChar(buf[i]);
          }
        } else if (opts.encoding === 'binary') {
          ret = buf;
        }
        FS.close(stream);
        return ret;
      },writeFile:function (path, data, opts) {
        opts = opts || {};
        opts.flags = opts.flags || 'w';
        opts.encoding = opts.encoding || 'utf8';
        if (opts.encoding !== 'utf8' && opts.encoding !== 'binary') {
          throw new Error('Invalid encoding type "' + opts.encoding + '"');
        }
        var stream = FS.open(path, opts.flags, opts.mode);
        if (opts.encoding === 'utf8') {
          var utf8 = new Runtime.UTF8Processor();
          var buf = new Uint8Array(utf8.processJSString(data));
          FS.write(stream, buf, 0, buf.length, 0, opts.canOwn);
        } else if (opts.encoding === 'binary') {
          FS.write(stream, data, 0, data.length, 0, opts.canOwn);
        }
        FS.close(stream);
      },cwd:function () {
        return FS.currentPath;
      },chdir:function (path) {
        var lookup = FS.lookupPath(path, { follow: true });
        if (!FS.isDir(lookup.node.mode)) {
          throw new FS.ErrnoError(ERRNO_CODES.ENOTDIR);
        }
        var err = FS.nodePermissions(lookup.node, 'x');
        if (err) {
          throw new FS.ErrnoError(err);
        }
        FS.currentPath = lookup.path;
      },createDefaultDirectories:function () {
        FS.mkdir('/tmp');
      },createDefaultDevices:function () {
        // create /dev
        FS.mkdir('/dev');
        // setup /dev/null
        FS.registerDevice(FS.makedev(1, 3), {
          read: function() { return 0; },
          write: function() { return 0; }
        });
        FS.mkdev('/dev/null', FS.makedev(1, 3));
        // setup /dev/tty and /dev/tty1
        // stderr needs to print output using Module['printErr']
        // so we register a second tty just for it.
        TTY.register(FS.makedev(5, 0), TTY.default_tty_ops);
        TTY.register(FS.makedev(6, 0), TTY.default_tty1_ops);
        FS.mkdev('/dev/tty', FS.makedev(5, 0));
        FS.mkdev('/dev/tty1', FS.makedev(6, 0));
        // we're not going to emulate the actual shm device,
        // just create the tmp dirs that reside in it commonly
        FS.mkdir('/dev/shm');
        FS.mkdir('/dev/shm/tmp');
      },createStandardStreams:function () {
        // TODO deprecate the old functionality of a single
        // input / output callback and that utilizes FS.createDevice
        // and instead require a unique set of stream ops
  
        // by default, we symlink the standard streams to the
        // default tty devices. however, if the standard streams
        // have been overwritten we create a unique device for
        // them instead.
        if (Module['stdin']) {
          FS.createDevice('/dev', 'stdin', Module['stdin']);
        } else {
          FS.symlink('/dev/tty', '/dev/stdin');
        }
        if (Module['stdout']) {
          FS.createDevice('/dev', 'stdout', null, Module['stdout']);
        } else {
          FS.symlink('/dev/tty', '/dev/stdout');
        }
        if (Module['stderr']) {
          FS.createDevice('/dev', 'stderr', null, Module['stderr']);
        } else {
          FS.symlink('/dev/tty1', '/dev/stderr');
        }
  
        // open default streams for the stdin, stdout and stderr devices
        var stdin = FS.open('/dev/stdin', 'r');
        HEAP32[((_stdin)>>2)]=FS.getPtrForStream(stdin);
        assert(stdin.fd === 0, 'invalid handle for stdin (' + stdin.fd + ')');
  
        var stdout = FS.open('/dev/stdout', 'w');
        HEAP32[((_stdout)>>2)]=FS.getPtrForStream(stdout);
        assert(stdout.fd === 1, 'invalid handle for stdout (' + stdout.fd + ')');
  
        var stderr = FS.open('/dev/stderr', 'w');
        HEAP32[((_stderr)>>2)]=FS.getPtrForStream(stderr);
        assert(stderr.fd === 2, 'invalid handle for stderr (' + stderr.fd + ')');
      },ensureErrnoError:function () {
        if (FS.ErrnoError) return;
        FS.ErrnoError = function ErrnoError(errno) {
          this.errno = errno;
          for (var key in ERRNO_CODES) {
            if (ERRNO_CODES[key] === errno) {
              this.code = key;
              break;
            }
          }
          this.message = ERRNO_MESSAGES[errno];
        };
        FS.ErrnoError.prototype = new Error();
        FS.ErrnoError.prototype.constructor = FS.ErrnoError;
        // Some errors may happen quite a bit, to avoid overhead we reuse them (and suffer a lack of stack info)
        [ERRNO_CODES.ENOENT].forEach(function(code) {
          FS.genericErrors[code] = new FS.ErrnoError(code);
          FS.genericErrors[code].stack = '<generic error, no stack>';
        });
      },staticInit:function () {
        FS.ensureErrnoError();
  
        FS.nameTable = new Array(4096);
  
        FS.mount(MEMFS, {}, '/');
  
        FS.createDefaultDirectories();
        FS.createDefaultDevices();
      },init:function (input, output, error) {
        assert(!FS.init.initialized, 'FS.init was previously called. If you want to initialize later with custom parameters, remove any earlier calls (note that one is automatically added to the generated code)');
        FS.init.initialized = true;
  
        FS.ensureErrnoError();
  
        // Allow Module.stdin etc. to provide defaults, if none explicitly passed to us here
        Module['stdin'] = input || Module['stdin'];
        Module['stdout'] = output || Module['stdout'];
        Module['stderr'] = error || Module['stderr'];
  
        FS.createStandardStreams();
      },quit:function () {
        FS.init.initialized = false;
        for (var i = 0; i < FS.streams.length; i++) {
          var stream = FS.streams[i];
          if (!stream) {
            continue;
          }
          FS.close(stream);
        }
      },getMode:function (canRead, canWrite) {
        var mode = 0;
        if (canRead) mode |= 292 | 73;
        if (canWrite) mode |= 146;
        return mode;
      },joinPath:function (parts, forceRelative) {
        var path = PATH.join.apply(null, parts);
        if (forceRelative && path[0] == '/') path = path.substr(1);
        return path;
      },absolutePath:function (relative, base) {
        return PATH.resolve(base, relative);
      },standardizePath:function (path) {
        return PATH.normalize(path);
      },findObject:function (path, dontResolveLastLink) {
        var ret = FS.analyzePath(path, dontResolveLastLink);
        if (ret.exists) {
          return ret.object;
        } else {
          ___setErrNo(ret.error);
          return null;
        }
      },analyzePath:function (path, dontResolveLastLink) {
        // operate from within the context of the symlink's target
        try {
          var lookup = FS.lookupPath(path, { follow: !dontResolveLastLink });
          path = lookup.path;
        } catch (e) {
        }
        var ret = {
          isRoot: false, exists: false, error: 0, name: null, path: null, object: null,
          parentExists: false, parentPath: null, parentObject: null
        };
        try {
          var lookup = FS.lookupPath(path, { parent: true });
          ret.parentExists = true;
          ret.parentPath = lookup.path;
          ret.parentObject = lookup.node;
          ret.name = PATH.basename(path);
          lookup = FS.lookupPath(path, { follow: !dontResolveLastLink });
          ret.exists = true;
          ret.path = lookup.path;
          ret.object = lookup.node;
          ret.name = lookup.node.name;
          ret.isRoot = lookup.path === '/';
        } catch (e) {
          ret.error = e.errno;
        };
        return ret;
      },createFolder:function (parent, name, canRead, canWrite) {
        var path = PATH.join2(typeof parent === 'string' ? parent : FS.getPath(parent), name);
        var mode = FS.getMode(canRead, canWrite);
        return FS.mkdir(path, mode);
      },createPath:function (parent, path, canRead, canWrite) {
        parent = typeof parent === 'string' ? parent : FS.getPath(parent);
        var parts = path.split('/').reverse();
        while (parts.length) {
          var part = parts.pop();
          if (!part) continue;
          var current = PATH.join2(parent, part);
          try {
            FS.mkdir(current);
          } catch (e) {
            // ignore EEXIST
          }
          parent = current;
        }
        return current;
      },createFile:function (parent, name, properties, canRead, canWrite) {
        var path = PATH.join2(typeof parent === 'string' ? parent : FS.getPath(parent), name);
        var mode = FS.getMode(canRead, canWrite);
        return FS.create(path, mode);
      },createDataFile:function (parent, name, data, canRead, canWrite, canOwn) {
        var path = name ? PATH.join2(typeof parent === 'string' ? parent : FS.getPath(parent), name) : parent;
        var mode = FS.getMode(canRead, canWrite);
        var node = FS.create(path, mode);
        if (data) {
          if (typeof data === 'string') {
            var arr = new Array(data.length);
            for (var i = 0, len = data.length; i < len; ++i) arr[i] = data.charCodeAt(i);
            data = arr;
          }
          // make sure we can write to the file
          FS.chmod(node, mode | 146);
          var stream = FS.open(node, 'w');
          FS.write(stream, data, 0, data.length, 0, canOwn);
          FS.close(stream);
          FS.chmod(node, mode);
        }
        return node;
      },createDevice:function (parent, name, input, output) {
        var path = PATH.join2(typeof parent === 'string' ? parent : FS.getPath(parent), name);
        var mode = FS.getMode(!!input, !!output);
        if (!FS.createDevice.major) FS.createDevice.major = 64;
        var dev = FS.makedev(FS.createDevice.major++, 0);
        // Create a fake device that a set of stream ops to emulate
        // the old behavior.
        FS.registerDevice(dev, {
          open: function(stream) {
            stream.seekable = false;
          },
          close: function(stream) {
            // flush any pending line data
            if (output && output.buffer && output.buffer.length) {
              output(10);
            }
          },
          read: function(stream, buffer, offset, length, pos /* ignored */) {
            var bytesRead = 0;
            for (var i = 0; i < length; i++) {
              var result;
              try {
                result = input();
              } catch (e) {
                throw new FS.ErrnoError(ERRNO_CODES.EIO);
              }
              if (result === undefined && bytesRead === 0) {
                throw new FS.ErrnoError(ERRNO_CODES.EAGAIN);
              }
              if (result === null || result === undefined) break;
              bytesRead++;
              buffer[offset+i] = result;
            }
            if (bytesRead) {
              stream.node.timestamp = Date.now();
            }
            return bytesRead;
          },
          write: function(stream, buffer, offset, length, pos) {
            for (var i = 0; i < length; i++) {
              try {
                output(buffer[offset+i]);
              } catch (e) {
                throw new FS.ErrnoError(ERRNO_CODES.EIO);
              }
            }
            if (length) {
              stream.node.timestamp = Date.now();
            }
            return i;
          }
        });
        return FS.mkdev(path, mode, dev);
      },createLink:function (parent, name, target, canRead, canWrite) {
        var path = PATH.join2(typeof parent === 'string' ? parent : FS.getPath(parent), name);
        return FS.symlink(target, path);
      },forceLoadFile:function (obj) {
        if (obj.isDevice || obj.isFolder || obj.link || obj.contents) return true;
        var success = true;
        if (typeof XMLHttpRequest !== 'undefined') {
          throw new Error("Lazy loading should have been performed (contents set) in createLazyFile, but it was not. Lazy loading only works in web workers. Use --embed-file or --preload-file in emcc on the main thread.");
        } else if (Module['read']) {
          // Command-line.
          try {
            // WARNING: Can't read binary files in V8's d8 or tracemonkey's js, as
            //          read() will try to parse UTF8.
            obj.contents = intArrayFromString(Module['read'](obj.url), true);
          } catch (e) {
            success = false;
          }
        } else {
          throw new Error('Cannot load without read() or XMLHttpRequest.');
        }
        if (!success) ___setErrNo(ERRNO_CODES.EIO);
        return success;
      },createLazyFile:function (parent, name, url, canRead, canWrite) {
        // Lazy chunked Uint8Array (implements get and length from Uint8Array). Actual getting is abstracted away for eventual reuse.
        function LazyUint8Array() {
          this.lengthKnown = false;
          this.chunks = []; // Loaded chunks. Index is the chunk number
        }
        LazyUint8Array.prototype.get = function LazyUint8Array_get(idx) {
          if (idx > this.length-1 || idx < 0) {
            return undefined;
          }
          var chunkOffset = idx % this.chunkSize;
          var chunkNum = Math.floor(idx / this.chunkSize);
          return this.getter(chunkNum)[chunkOffset];
        }
        LazyUint8Array.prototype.setDataGetter = function LazyUint8Array_setDataGetter(getter) {
          this.getter = getter;
        }
        LazyUint8Array.prototype.cacheLength = function LazyUint8Array_cacheLength() {
            // Find length
            var xhr = new XMLHttpRequest();
            xhr.open('HEAD', url, false);
            xhr.send(null);
            if (!(xhr.status >= 200 && xhr.status < 300 || xhr.status === 304)) throw new Error("Couldn't load " + url + ". Status: " + xhr.status);
            var datalength = Number(xhr.getResponseHeader("Content-length"));
            var header;
            var hasByteServing = (header = xhr.getResponseHeader("Accept-Ranges")) && header === "bytes";
            var chunkSize = 1024*1024; // Chunk size in bytes
  
            if (!hasByteServing) chunkSize = datalength;
  
            // Function to get a range from the remote URL.
            var doXHR = (function(from, to) {
              if (from > to) throw new Error("invalid range (" + from + ", " + to + ") or no bytes requested!");
              if (to > datalength-1) throw new Error("only " + datalength + " bytes available! programmer error!");
  
              // TODO: Use mozResponseArrayBuffer, responseStream, etc. if available.
              var xhr = new XMLHttpRequest();
              xhr.open('GET', url, false);
              if (datalength !== chunkSize) xhr.setRequestHeader("Range", "bytes=" + from + "-" + to);
  
              // Some hints to the browser that we want binary data.
              if (typeof Uint8Array != 'undefined') xhr.responseType = 'arraybuffer';
              if (xhr.overrideMimeType) {
                xhr.overrideMimeType('text/plain; charset=x-user-defined');
              }
  
              xhr.send(null);
              if (!(xhr.status >= 200 && xhr.status < 300 || xhr.status === 304)) throw new Error("Couldn't load " + url + ". Status: " + xhr.status);
              if (xhr.response !== undefined) {
                return new Uint8Array(xhr.response || []);
              } else {
                return intArrayFromString(xhr.responseText || '', true);
              }
            });
            var lazyArray = this;
            lazyArray.setDataGetter(function(chunkNum) {
              var start = chunkNum * chunkSize;
              var end = (chunkNum+1) * chunkSize - 1; // including this byte
              end = Math.min(end, datalength-1); // if datalength-1 is selected, this is the last block
              if (typeof(lazyArray.chunks[chunkNum]) === "undefined") {
                lazyArray.chunks[chunkNum] = doXHR(start, end);
              }
              if (typeof(lazyArray.chunks[chunkNum]) === "undefined") throw new Error("doXHR failed!");
              return lazyArray.chunks[chunkNum];
            });
  
            this._length = datalength;
            this._chunkSize = chunkSize;
            this.lengthKnown = true;
        }
        if (typeof XMLHttpRequest !== 'undefined') {
          if (!ENVIRONMENT_IS_WORKER) throw 'Cannot do synchronous binary XHRs outside webworkers in modern browsers. Use --embed-file or --preload-file in emcc';
          var lazyArray = new LazyUint8Array();
          Object.defineProperty(lazyArray, "length", {
              get: function() {
                  if(!this.lengthKnown) {
                      this.cacheLength();
                  }
                  return this._length;
              }
          });
          Object.defineProperty(lazyArray, "chunkSize", {
              get: function() {
                  if(!this.lengthKnown) {
                      this.cacheLength();
                  }
                  return this._chunkSize;
              }
          });
  
          var properties = { isDevice: false, contents: lazyArray };
        } else {
          var properties = { isDevice: false, url: url };
        }
  
        var node = FS.createFile(parent, name, properties, canRead, canWrite);
        // This is a total hack, but I want to get this lazy file code out of the
        // core of MEMFS. If we want to keep this lazy file concept I feel it should
        // be its own thin LAZYFS proxying calls to MEMFS.
        if (properties.contents) {
          node.contents = properties.contents;
        } else if (properties.url) {
          node.contents = null;
          node.url = properties.url;
        }
        // override each stream op with one that tries to force load the lazy file first
        var stream_ops = {};
        var keys = Object.keys(node.stream_ops);
        keys.forEach(function(key) {
          var fn = node.stream_ops[key];
          stream_ops[key] = function forceLoadLazyFile() {
            if (!FS.forceLoadFile(node)) {
              throw new FS.ErrnoError(ERRNO_CODES.EIO);
            }
            return fn.apply(null, arguments);
          };
        });
        // use a custom read function
        stream_ops.read = function stream_ops_read(stream, buffer, offset, length, position) {
          if (!FS.forceLoadFile(node)) {
            throw new FS.ErrnoError(ERRNO_CODES.EIO);
          }
          var contents = stream.node.contents;
          if (position >= contents.length)
            return 0;
          var size = Math.min(contents.length - position, length);
          assert(size >= 0);
          if (contents.slice) { // normal array
            for (var i = 0; i < size; i++) {
              buffer[offset + i] = contents[position + i];
            }
          } else {
            for (var i = 0; i < size; i++) { // LazyUint8Array from sync binary XHR
              buffer[offset + i] = contents.get(position + i);
            }
          }
          return size;
        };
        node.stream_ops = stream_ops;
        return node;
      },createPreloadedFile:function (parent, name, url, canRead, canWrite, onload, onerror, dontCreateFile, canOwn) {
        Browser.init();
        // TODO we should allow people to just pass in a complete filename instead
        // of parent and name being that we just join them anyways
        var fullname = name ? PATH.resolve(PATH.join2(parent, name)) : parent;
        function processData(byteArray) {
          function finish(byteArray) {
            if (!dontCreateFile) {
              FS.createDataFile(parent, name, byteArray, canRead, canWrite, canOwn);
            }
            if (onload) onload();
            removeRunDependency('cp ' + fullname);
          }
          var handled = false;
          Module['preloadPlugins'].forEach(function(plugin) {
            if (handled) return;
            if (plugin['canHandle'](fullname)) {
              plugin['handle'](byteArray, fullname, finish, function() {
                if (onerror) onerror();
                removeRunDependency('cp ' + fullname);
              });
              handled = true;
            }
          });
          if (!handled) finish(byteArray);
        }
        addRunDependency('cp ' + fullname);
        if (typeof url == 'string') {
          Browser.asyncLoad(url, function(byteArray) {
            processData(byteArray);
          }, onerror);
        } else {
          processData(url);
        }
      },indexedDB:function () {
        return window.indexedDB || window.mozIndexedDB || window.webkitIndexedDB || window.msIndexedDB;
      },DB_NAME:function () {
        return 'EM_FS_' + window.location.pathname;
      },DB_VERSION:20,DB_STORE_NAME:"FILE_DATA",saveFilesToDB:function (paths, onload, onerror) {
        onload = onload || function(){};
        onerror = onerror || function(){};
        var indexedDB = FS.indexedDB();
        try {
          var openRequest = indexedDB.open(FS.DB_NAME(), FS.DB_VERSION);
        } catch (e) {
          return onerror(e);
        }
        openRequest.onupgradeneeded = function openRequest_onupgradeneeded() {
          console.log('creating db');
          var db = openRequest.result;
          db.createObjectStore(FS.DB_STORE_NAME);
        };
        openRequest.onsuccess = function openRequest_onsuccess() {
          var db = openRequest.result;
          var transaction = db.transaction([FS.DB_STORE_NAME], 'readwrite');
          var files = transaction.objectStore(FS.DB_STORE_NAME);
          var ok = 0, fail = 0, total = paths.length;
          function finish() {
            if (fail == 0) onload(); else onerror();
          }
          paths.forEach(function(path) {
            var putRequest = files.put(FS.analyzePath(path).object.contents, path);
            putRequest.onsuccess = function putRequest_onsuccess() { ok++; if (ok + fail == total) finish() };
            putRequest.onerror = function putRequest_onerror() { fail++; if (ok + fail == total) finish() };
          });
          transaction.onerror = onerror;
        };
        openRequest.onerror = onerror;
      },loadFilesFromDB:function (paths, onload, onerror) {
        onload = onload || function(){};
        onerror = onerror || function(){};
        var indexedDB = FS.indexedDB();
        try {
          var openRequest = indexedDB.open(FS.DB_NAME(), FS.DB_VERSION);
        } catch (e) {
          return onerror(e);
        }
        openRequest.onupgradeneeded = onerror; // no database to load from
        openRequest.onsuccess = function openRequest_onsuccess() {
          var db = openRequest.result;
          try {
            var transaction = db.transaction([FS.DB_STORE_NAME], 'readonly');
          } catch(e) {
            onerror(e);
            return;
          }
          var files = transaction.objectStore(FS.DB_STORE_NAME);
          var ok = 0, fail = 0, total = paths.length;
          function finish() {
            if (fail == 0) onload(); else onerror();
          }
          paths.forEach(function(path) {
            var getRequest = files.get(path);
            getRequest.onsuccess = function getRequest_onsuccess() {
              if (FS.analyzePath(path).exists) {
                FS.unlink(path);
              }
              FS.createDataFile(PATH.dirname(path), PATH.basename(path), getRequest.result, true, true, true);
              ok++;
              if (ok + fail == total) finish();
            };
            getRequest.onerror = function getRequest_onerror() { fail++; if (ok + fail == total) finish() };
          });
          transaction.onerror = onerror;
        };
        openRequest.onerror = onerror;
      }};function _fflush(stream) {
      // int fflush(FILE *stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fflush.html
      // we don't currently perform any user-space buffering of data
    }

  var _sqrtf=Math_sqrt;

  
  
  
  function _isspace(chr) {
      return (chr == 32) || (chr >= 9 && chr <= 13);
    }
  function __parseInt64(str, endptr, base, min, max, unsign) {
      var isNegative = false;
      // Skip space.
      while (_isspace(HEAP8[(str)])) str++;
  
      // Check for a plus/minus sign.
      if (HEAP8[(str)] == 45) {
        str++;
        isNegative = true;
      } else if (HEAP8[(str)] == 43) {
        str++;
      }
  
      // Find base.
      var ok = false;
      var finalBase = base;
      if (!finalBase) {
        if (HEAP8[(str)] == 48) {
          if (HEAP8[((str+1)|0)] == 120 ||
              HEAP8[((str+1)|0)] == 88) {
            finalBase = 16;
            str += 2;
          } else {
            finalBase = 8;
            ok = true; // we saw an initial zero, perhaps the entire thing is just "0"
          }
        }
      } else if (finalBase==16) {
        if (HEAP8[(str)] == 48) {
          if (HEAP8[((str+1)|0)] == 120 ||
              HEAP8[((str+1)|0)] == 88) {
            str += 2;
          }
        }
      }
      if (!finalBase) finalBase = 10;
      var start = str;
  
      // Get digits.
      var chr;
      while ((chr = HEAP8[(str)]) != 0) {
        var digit = parseInt(String.fromCharCode(chr), finalBase);
        if (isNaN(digit)) {
          break;
        } else {
          str++;
          ok = true;
        }
      }
  
      if (!ok) {
        ___setErrNo(ERRNO_CODES.EINVAL);
        return ((asm["setTempRet0"](0),0)|0);
      }
  
      // Set end pointer.
      if (endptr) {
        HEAP32[((endptr)>>2)]=str;
      }
  
      try {
        var numberString = isNegative ? '-'+Pointer_stringify(start, str - start) : Pointer_stringify(start, str - start);
        i64Math.fromString(numberString, finalBase, min, max, unsign);
      } catch(e) {
        ___setErrNo(ERRNO_CODES.ERANGE); // not quite correct
      }
  
      return ((asm["setTempRet0"](((HEAP32[(((tempDoublePtr)+(4))>>2)])|0)),((HEAP32[((tempDoublePtr)>>2)])|0))|0);
    }function _strtoull(str, endptr, base) {
      return __parseInt64(str, endptr, base, 0, '18446744073709551615', true);  // ULONG_MAX.
    }function _strtoull_l(str, endptr, base) {
      return _strtoull(str, endptr, base); // no locale support yet
    }

  
  
   
  Module["_strlen"] = _strlen;
  
  function __reallyNegative(x) {
      return x < 0 || (x === 0 && (1/x) === -Infinity);
    }function __formatString(format, varargs) {
      var textIndex = format;
      var argIndex = 0;
      function getNextArg(type) {
        // NOTE: Explicitly ignoring type safety. Otherwise this fails:
        //       int x = 4; printf("%c\n", (char)x);
        var ret;
        if (type === 'double') {
          ret = (HEAP32[((tempDoublePtr)>>2)]=HEAP32[(((varargs)+(argIndex))>>2)],HEAP32[(((tempDoublePtr)+(4))>>2)]=HEAP32[(((varargs)+((argIndex)+(4)))>>2)],(+(HEAPF64[(tempDoublePtr)>>3])));
        } else if (type == 'i64') {
          ret = [HEAP32[(((varargs)+(argIndex))>>2)],
                 HEAP32[(((varargs)+(argIndex+4))>>2)]];
  
        } else {
          type = 'i32'; // varargs are always i32, i64, or double
          ret = HEAP32[(((varargs)+(argIndex))>>2)];
        }
        argIndex += Runtime.getNativeFieldSize(type);
        return ret;
      }
  
      var ret = [];
      var curr, next, currArg;
      while(1) {
        var startTextIndex = textIndex;
        curr = HEAP8[(textIndex)];
        if (curr === 0) break;
        next = HEAP8[((textIndex+1)|0)];
        if (curr == 37) {
          // Handle flags.
          var flagAlwaysSigned = false;
          var flagLeftAlign = false;
          var flagAlternative = false;
          var flagZeroPad = false;
          var flagPadSign = false;
          flagsLoop: while (1) {
            switch (next) {
              case 43:
                flagAlwaysSigned = true;
                break;
              case 45:
                flagLeftAlign = true;
                break;
              case 35:
                flagAlternative = true;
                break;
              case 48:
                if (flagZeroPad) {
                  break flagsLoop;
                } else {
                  flagZeroPad = true;
                  break;
                }
              case 32:
                flagPadSign = true;
                break;
              default:
                break flagsLoop;
            }
            textIndex++;
            next = HEAP8[((textIndex+1)|0)];
          }
  
          // Handle width.
          var width = 0;
          if (next == 42) {
            width = getNextArg('i32');
            textIndex++;
            next = HEAP8[((textIndex+1)|0)];
          } else {
            while (next >= 48 && next <= 57) {
              width = width * 10 + (next - 48);
              textIndex++;
              next = HEAP8[((textIndex+1)|0)];
            }
          }
  
          // Handle precision.
          var precisionSet = false, precision = -1;
          if (next == 46) {
            precision = 0;
            precisionSet = true;
            textIndex++;
            next = HEAP8[((textIndex+1)|0)];
            if (next == 42) {
              precision = getNextArg('i32');
              textIndex++;
            } else {
              while(1) {
                var precisionChr = HEAP8[((textIndex+1)|0)];
                if (precisionChr < 48 ||
                    precisionChr > 57) break;
                precision = precision * 10 + (precisionChr - 48);
                textIndex++;
              }
            }
            next = HEAP8[((textIndex+1)|0)];
          }
          if (precision < 0) {
            precision = 6; // Standard default.
            precisionSet = false;
          }
  
          // Handle integer sizes. WARNING: These assume a 32-bit architecture!
          var argSize;
          switch (String.fromCharCode(next)) {
            case 'h':
              var nextNext = HEAP8[((textIndex+2)|0)];
              if (nextNext == 104) {
                textIndex++;
                argSize = 1; // char (actually i32 in varargs)
              } else {
                argSize = 2; // short (actually i32 in varargs)
              }
              break;
            case 'l':
              var nextNext = HEAP8[((textIndex+2)|0)];
              if (nextNext == 108) {
                textIndex++;
                argSize = 8; // long long
              } else {
                argSize = 4; // long
              }
              break;
            case 'L': // long long
            case 'q': // int64_t
            case 'j': // intmax_t
              argSize = 8;
              break;
            case 'z': // size_t
            case 't': // ptrdiff_t
            case 'I': // signed ptrdiff_t or unsigned size_t
              argSize = 4;
              break;
            default:
              argSize = null;
          }
          if (argSize) textIndex++;
          next = HEAP8[((textIndex+1)|0)];
  
          // Handle type specifier.
          switch (String.fromCharCode(next)) {
            case 'd': case 'i': case 'u': case 'o': case 'x': case 'X': case 'p': {
              // Integer.
              var signed = next == 100 || next == 105;
              argSize = argSize || 4;
              var currArg = getNextArg('i' + (argSize * 8));
              var origArg = currArg;
              var argText;
              // Flatten i64-1 [low, high] into a (slightly rounded) double
              if (argSize == 8) {
                currArg = Runtime.makeBigInt(currArg[0], currArg[1], next == 117);
              }
              // Truncate to requested size.
              if (argSize <= 4) {
                var limit = Math.pow(256, argSize) - 1;
                currArg = (signed ? reSign : unSign)(currArg & limit, argSize * 8);
              }
              // Format the number.
              var currAbsArg = Math.abs(currArg);
              var prefix = '';
              if (next == 100 || next == 105) {
                if (argSize == 8 && i64Math) argText = i64Math.stringify(origArg[0], origArg[1], null); else
                argText = reSign(currArg, 8 * argSize, 1).toString(10);
              } else if (next == 117) {
                if (argSize == 8 && i64Math) argText = i64Math.stringify(origArg[0], origArg[1], true); else
                argText = unSign(currArg, 8 * argSize, 1).toString(10);
                currArg = Math.abs(currArg);
              } else if (next == 111) {
                argText = (flagAlternative ? '0' : '') + currAbsArg.toString(8);
              } else if (next == 120 || next == 88) {
                prefix = (flagAlternative && currArg != 0) ? '0x' : '';
                if (argSize == 8 && i64Math) {
                  if (origArg[1]) {
                    argText = (origArg[1]>>>0).toString(16);
                    var lower = (origArg[0]>>>0).toString(16);
                    while (lower.length < 8) lower = '0' + lower;
                    argText += lower;
                  } else {
                    argText = (origArg[0]>>>0).toString(16);
                  }
                } else
                if (currArg < 0) {
                  // Represent negative numbers in hex as 2's complement.
                  currArg = -currArg;
                  argText = (currAbsArg - 1).toString(16);
                  var buffer = [];
                  for (var i = 0; i < argText.length; i++) {
                    buffer.push((0xF - parseInt(argText[i], 16)).toString(16));
                  }
                  argText = buffer.join('');
                  while (argText.length < argSize * 2) argText = 'f' + argText;
                } else {
                  argText = currAbsArg.toString(16);
                }
                if (next == 88) {
                  prefix = prefix.toUpperCase();
                  argText = argText.toUpperCase();
                }
              } else if (next == 112) {
                if (currAbsArg === 0) {
                  argText = '(nil)';
                } else {
                  prefix = '0x';
                  argText = currAbsArg.toString(16);
                }
              }
              if (precisionSet) {
                while (argText.length < precision) {
                  argText = '0' + argText;
                }
              }
  
              // Add sign if needed
              if (currArg >= 0) {
                if (flagAlwaysSigned) {
                  prefix = '+' + prefix;
                } else if (flagPadSign) {
                  prefix = ' ' + prefix;
                }
              }
  
              // Move sign to prefix so we zero-pad after the sign
              if (argText.charAt(0) == '-') {
                prefix = '-' + prefix;
                argText = argText.substr(1);
              }
  
              // Add padding.
              while (prefix.length + argText.length < width) {
                if (flagLeftAlign) {
                  argText += ' ';
                } else {
                  if (flagZeroPad) {
                    argText = '0' + argText;
                  } else {
                    prefix = ' ' + prefix;
                  }
                }
              }
  
              // Insert the result into the buffer.
              argText = prefix + argText;
              argText.split('').forEach(function(chr) {
                ret.push(chr.charCodeAt(0));
              });
              break;
            }
            case 'f': case 'F': case 'e': case 'E': case 'g': case 'G': {
              // Float.
              var currArg = getNextArg('double');
              var argText;
              if (isNaN(currArg)) {
                argText = 'nan';
                flagZeroPad = false;
              } else if (!isFinite(currArg)) {
                argText = (currArg < 0 ? '-' : '') + 'inf';
                flagZeroPad = false;
              } else {
                var isGeneral = false;
                var effectivePrecision = Math.min(precision, 20);
  
                // Convert g/G to f/F or e/E, as per:
                // http://pubs.opengroup.org/onlinepubs/9699919799/functions/printf.html
                if (next == 103 || next == 71) {
                  isGeneral = true;
                  precision = precision || 1;
                  var exponent = parseInt(currArg.toExponential(effectivePrecision).split('e')[1], 10);
                  if (precision > exponent && exponent >= -4) {
                    next = ((next == 103) ? 'f' : 'F').charCodeAt(0);
                    precision -= exponent + 1;
                  } else {
                    next = ((next == 103) ? 'e' : 'E').charCodeAt(0);
                    precision--;
                  }
                  effectivePrecision = Math.min(precision, 20);
                }
  
                if (next == 101 || next == 69) {
                  argText = currArg.toExponential(effectivePrecision);
                  // Make sure the exponent has at least 2 digits.
                  if (/[eE][-+]\d$/.test(argText)) {
                    argText = argText.slice(0, -1) + '0' + argText.slice(-1);
                  }
                } else if (next == 102 || next == 70) {
                  argText = currArg.toFixed(effectivePrecision);
                  if (currArg === 0 && __reallyNegative(currArg)) {
                    argText = '-' + argText;
                  }
                }
  
                var parts = argText.split('e');
                if (isGeneral && !flagAlternative) {
                  // Discard trailing zeros and periods.
                  while (parts[0].length > 1 && parts[0].indexOf('.') != -1 &&
                         (parts[0].slice(-1) == '0' || parts[0].slice(-1) == '.')) {
                    parts[0] = parts[0].slice(0, -1);
                  }
                } else {
                  // Make sure we have a period in alternative mode.
                  if (flagAlternative && argText.indexOf('.') == -1) parts[0] += '.';
                  // Zero pad until required precision.
                  while (precision > effectivePrecision++) parts[0] += '0';
                }
                argText = parts[0] + (parts.length > 1 ? 'e' + parts[1] : '');
  
                // Capitalize 'E' if needed.
                if (next == 69) argText = argText.toUpperCase();
  
                // Add sign.
                if (currArg >= 0) {
                  if (flagAlwaysSigned) {
                    argText = '+' + argText;
                  } else if (flagPadSign) {
                    argText = ' ' + argText;
                  }
                }
              }
  
              // Add padding.
              while (argText.length < width) {
                if (flagLeftAlign) {
                  argText += ' ';
                } else {
                  if (flagZeroPad && (argText[0] == '-' || argText[0] == '+')) {
                    argText = argText[0] + '0' + argText.slice(1);
                  } else {
                    argText = (flagZeroPad ? '0' : ' ') + argText;
                  }
                }
              }
  
              // Adjust case.
              if (next < 97) argText = argText.toUpperCase();
  
              // Insert the result into the buffer.
              argText.split('').forEach(function(chr) {
                ret.push(chr.charCodeAt(0));
              });
              break;
            }
            case 's': {
              // String.
              var arg = getNextArg('i8*');
              var argLength = arg ? _strlen(arg) : '(null)'.length;
              if (precisionSet) argLength = Math.min(argLength, precision);
              if (!flagLeftAlign) {
                while (argLength < width--) {
                  ret.push(32);
                }
              }
              if (arg) {
                for (var i = 0; i < argLength; i++) {
                  ret.push(HEAPU8[((arg++)|0)]);
                }
              } else {
                ret = ret.concat(intArrayFromString('(null)'.substr(0, argLength), true));
              }
              if (flagLeftAlign) {
                while (argLength < width--) {
                  ret.push(32);
                }
              }
              break;
            }
            case 'c': {
              // Character.
              if (flagLeftAlign) ret.push(getNextArg('i8'));
              while (--width > 0) {
                ret.push(32);
              }
              if (!flagLeftAlign) ret.push(getNextArg('i8'));
              break;
            }
            case 'n': {
              // Write the length written so far to the next parameter.
              var ptr = getNextArg('i32*');
              HEAP32[((ptr)>>2)]=ret.length;
              break;
            }
            case '%': {
              // Literal percent sign.
              ret.push(curr);
              break;
            }
            default: {
              // Unknown specifiers remain untouched.
              for (var i = startTextIndex; i < textIndex + 2; i++) {
                ret.push(HEAP8[(i)]);
              }
            }
          }
          textIndex += 2;
          // TODO: Support a/A (hex float) and m (last error) specifiers.
          // TODO: Support %1${specifier} for arg selection.
        } else {
          ret.push(curr);
          textIndex += 1;
        }
      }
      return ret;
    }
  
  function _malloc(bytes) {
      /* Over-allocate to make sure it is byte-aligned by 8.
       * This will leak memory, but this is only the dummy
       * implementation (replaced by dlmalloc normally) so
       * not an issue.
       */
      var ptr = Runtime.dynamicAlloc(bytes + 8);
      return (ptr+8) & 0xFFFFFFF8;
    }
  Module["_malloc"] = _malloc;function _snprintf(s, n, format, varargs) {
      // int snprintf(char *restrict s, size_t n, const char *restrict format, ...);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/printf.html
      var result = __formatString(format, varargs);
      var limit = (n === undefined) ? result.length
                                    : Math.min(result.length, Math.max(n - 1, 0));
      if (s < 0) {
        s = -s;
        var buf = _malloc(limit+1);
        HEAP32[((s)>>2)]=buf;
        s = buf;
      }
      for (var i = 0; i < limit; i++) {
        HEAP8[(((s)+(i))|0)]=result[i];
      }
      if (limit < n || (n === undefined)) HEAP8[(((s)+(i))|0)]=0;
      return result.length;
    }

  
  
  
  
  function _mkport() { throw 'TODO' }var SOCKFS={mount:function (mount) {
        return FS.createNode(null, '/', 16384 | 511 /* 0777 */, 0);
      },createSocket:function (family, type, protocol) {
        var streaming = type == 1;
        if (protocol) {
          assert(streaming == (protocol == 6)); // if SOCK_STREAM, must be tcp
        }
  
        // create our internal socket structure
        var sock = {
          family: family,
          type: type,
          protocol: protocol,
          server: null,
          peers: {},
          pending: [],
          recv_queue: [],
          sock_ops: SOCKFS.websocket_sock_ops
        };
  
        // create the filesystem node to store the socket structure
        var name = SOCKFS.nextname();
        var node = FS.createNode(SOCKFS.root, name, 49152, 0);
        node.sock = sock;
  
        // and the wrapping stream that enables library functions such
        // as read and write to indirectly interact with the socket
        var stream = FS.createStream({
          path: name,
          node: node,
          flags: FS.modeStringToFlags('r+'),
          seekable: false,
          stream_ops: SOCKFS.stream_ops
        });
  
        // map the new stream to the socket structure (sockets have a 1:1
        // relationship with a stream)
        sock.stream = stream;
  
        return sock;
      },getSocket:function (fd) {
        var stream = FS.getStream(fd);
        if (!stream || !FS.isSocket(stream.node.mode)) {
          return null;
        }
        return stream.node.sock;
      },stream_ops:{poll:function (stream) {
          var sock = stream.node.sock;
          return sock.sock_ops.poll(sock);
        },ioctl:function (stream, request, varargs) {
          var sock = stream.node.sock;
          return sock.sock_ops.ioctl(sock, request, varargs);
        },read:function (stream, buffer, offset, length, position /* ignored */) {
          var sock = stream.node.sock;
          var msg = sock.sock_ops.recvmsg(sock, length);
          if (!msg) {
            // socket is closed
            return 0;
          }
          buffer.set(msg.buffer, offset);
          return msg.buffer.length;
        },write:function (stream, buffer, offset, length, position /* ignored */) {
          var sock = stream.node.sock;
          return sock.sock_ops.sendmsg(sock, buffer, offset, length);
        },close:function (stream) {
          var sock = stream.node.sock;
          sock.sock_ops.close(sock);
        }},nextname:function () {
        if (!SOCKFS.nextname.current) {
          SOCKFS.nextname.current = 0;
        }
        return 'socket[' + (SOCKFS.nextname.current++) + ']';
      },websocket_sock_ops:{createPeer:function (sock, addr, port) {
          var ws;
  
          if (typeof addr === 'object') {
            ws = addr;
            addr = null;
            port = null;
          }
  
          if (ws) {
            // for sockets that've already connected (e.g. we're the server)
            // we can inspect the _socket property for the address
            if (ws._socket) {
              addr = ws._socket.remoteAddress;
              port = ws._socket.remotePort;
            }
            // if we're just now initializing a connection to the remote,
            // inspect the url property
            else {
              var result = /ws[s]?:\/\/([^:]+):(\d+)/.exec(ws.url);
              if (!result) {
                throw new Error('WebSocket URL must be in the format ws(s)://address:port');
              }
              addr = result[1];
              port = parseInt(result[2], 10);
            }
          } else {
            // create the actual websocket object and connect
            try {
              // runtimeConfig gets set to true if WebSocket runtime configuration is available.
              var runtimeConfig = (Module['websocket'] && ('object' === typeof Module['websocket']));
  
              // The default value is 'ws://' the replace is needed because the compiler replaces "//" comments with '#'
              // comments without checking context, so we'd end up with ws:#, the replace swaps the "#" for "//" again.
              var url = 'ws:#'.replace('#', '//');
  
              if (runtimeConfig) {
                if ('string' === typeof Module['websocket']['url']) {
                  url = Module['websocket']['url']; // Fetch runtime WebSocket URL config.
                }
              }
  
              if (url === 'ws://' || url === 'wss://') { // Is the supplied URL config just a prefix, if so complete it.
                url = url + addr + ':' + port;
              }
  
              // Make the WebSocket subprotocol (Sec-WebSocket-Protocol) default to binary if no configuration is set.
              var subProtocols = 'binary'; // The default value is 'binary'
  
              if (runtimeConfig) {
                if ('string' === typeof Module['websocket']['subprotocol']) {
                  subProtocols = Module['websocket']['subprotocol']; // Fetch runtime WebSocket subprotocol config.
                }
              }
  
              // The regex trims the string (removes spaces at the beginning and end, then splits the string by
              // <any space>,<any space> into an Array. Whitespace removal is important for Websockify and ws.
              subProtocols = subProtocols.replace(/^ +| +$/g,"").split(/ *, */);
  
              // The node ws library API for specifying optional subprotocol is slightly different than the browser's.
              var opts = ENVIRONMENT_IS_NODE ? {'protocol': subProtocols.toString()} : subProtocols;
  
              // If node we use the ws library.
              var WebSocket = ENVIRONMENT_IS_NODE ? require('ws') : window['WebSocket'];
              ws = new WebSocket(url, opts);
              ws.binaryType = 'arraybuffer';
            } catch (e) {
              throw new FS.ErrnoError(ERRNO_CODES.EHOSTUNREACH);
            }
          }
  
  
          var peer = {
            addr: addr,
            port: port,
            socket: ws,
            dgram_send_queue: []
          };
  
          SOCKFS.websocket_sock_ops.addPeer(sock, peer);
          SOCKFS.websocket_sock_ops.handlePeerEvents(sock, peer);
  
          // if this is a bound dgram socket, send the port number first to allow
          // us to override the ephemeral port reported to us by remotePort on the
          // remote end.
          if (sock.type === 2 && typeof sock.sport !== 'undefined') {
            peer.dgram_send_queue.push(new Uint8Array([
                255, 255, 255, 255,
                'p'.charCodeAt(0), 'o'.charCodeAt(0), 'r'.charCodeAt(0), 't'.charCodeAt(0),
                ((sock.sport & 0xff00) >> 8) , (sock.sport & 0xff)
            ]));
          }
  
          return peer;
        },getPeer:function (sock, addr, port) {
          return sock.peers[addr + ':' + port];
        },addPeer:function (sock, peer) {
          sock.peers[peer.addr + ':' + peer.port] = peer;
        },removePeer:function (sock, peer) {
          delete sock.peers[peer.addr + ':' + peer.port];
        },handlePeerEvents:function (sock, peer) {
          var first = true;
  
          var handleOpen = function () {
            try {
              var queued = peer.dgram_send_queue.shift();
              while (queued) {
                peer.socket.send(queued);
                queued = peer.dgram_send_queue.shift();
              }
            } catch (e) {
              // not much we can do here in the way of proper error handling as we've already
              // lied and said this data was sent. shut it down.
              peer.socket.close();
            }
          };
  
          function handleMessage(data) {
            assert(typeof data !== 'string' && data.byteLength !== undefined);  // must receive an ArrayBuffer
            data = new Uint8Array(data);  // make a typed array view on the array buffer
  
  
            // if this is the port message, override the peer's port with it
            var wasfirst = first;
            first = false;
            if (wasfirst &&
                data.length === 10 &&
                data[0] === 255 && data[1] === 255 && data[2] === 255 && data[3] === 255 &&
                data[4] === 'p'.charCodeAt(0) && data[5] === 'o'.charCodeAt(0) && data[6] === 'r'.charCodeAt(0) && data[7] === 't'.charCodeAt(0)) {
              // update the peer's port and it's key in the peer map
              var newport = ((data[8] << 8) | data[9]);
              SOCKFS.websocket_sock_ops.removePeer(sock, peer);
              peer.port = newport;
              SOCKFS.websocket_sock_ops.addPeer(sock, peer);
              return;
            }
  
            sock.recv_queue.push({ addr: peer.addr, port: peer.port, data: data });
          };
  
          if (ENVIRONMENT_IS_NODE) {
            peer.socket.on('open', handleOpen);
            peer.socket.on('message', function(data, flags) {
              if (!flags.binary) {
                return;
              }
              handleMessage((new Uint8Array(data)).buffer);  // copy from node Buffer -> ArrayBuffer
            });
            peer.socket.on('error', function() {
              // don't throw
            });
          } else {
            peer.socket.onopen = handleOpen;
            peer.socket.onmessage = function peer_socket_onmessage(event) {
              handleMessage(event.data);
            };
          }
        },poll:function (sock) {
          if (sock.type === 1 && sock.server) {
            // listen sockets should only say they're available for reading
            // if there are pending clients.
            return sock.pending.length ? (64 | 1) : 0;
          }
  
          var mask = 0;
          var dest = sock.type === 1 ?  // we only care about the socket state for connection-based sockets
            SOCKFS.websocket_sock_ops.getPeer(sock, sock.daddr, sock.dport) :
            null;
  
          if (sock.recv_queue.length ||
              !dest ||  // connection-less sockets are always ready to read
              (dest && dest.socket.readyState === dest.socket.CLOSING) ||
              (dest && dest.socket.readyState === dest.socket.CLOSED)) {  // let recv return 0 once closed
            mask |= (64 | 1);
          }
  
          if (!dest ||  // connection-less sockets are always ready to write
              (dest && dest.socket.readyState === dest.socket.OPEN)) {
            mask |= 4;
          }
  
          if ((dest && dest.socket.readyState === dest.socket.CLOSING) ||
              (dest && dest.socket.readyState === dest.socket.CLOSED)) {
            mask |= 16;
          }
  
          return mask;
        },ioctl:function (sock, request, arg) {
          switch (request) {
            case 21531:
              var bytes = 0;
              if (sock.recv_queue.length) {
                bytes = sock.recv_queue[0].data.length;
              }
              HEAP32[((arg)>>2)]=bytes;
              return 0;
            default:
              return ERRNO_CODES.EINVAL;
          }
        },close:function (sock) {
          // if we've spawned a listen server, close it
          if (sock.server) {
            try {
              sock.server.close();
            } catch (e) {
            }
            sock.server = null;
          }
          // close any peer connections
          var peers = Object.keys(sock.peers);
          for (var i = 0; i < peers.length; i++) {
            var peer = sock.peers[peers[i]];
            try {
              peer.socket.close();
            } catch (e) {
            }
            SOCKFS.websocket_sock_ops.removePeer(sock, peer);
          }
          return 0;
        },bind:function (sock, addr, port) {
          if (typeof sock.saddr !== 'undefined' || typeof sock.sport !== 'undefined') {
            throw new FS.ErrnoError(ERRNO_CODES.EINVAL);  // already bound
          }
          sock.saddr = addr;
          sock.sport = port || _mkport();
          // in order to emulate dgram sockets, we need to launch a listen server when
          // binding on a connection-less socket
          // note: this is only required on the server side
          if (sock.type === 2) {
            // close the existing server if it exists
            if (sock.server) {
              sock.server.close();
              sock.server = null;
            }
            // swallow error operation not supported error that occurs when binding in the
            // browser where this isn't supported
            try {
              sock.sock_ops.listen(sock, 0);
            } catch (e) {
              if (!(e instanceof FS.ErrnoError)) throw e;
              if (e.errno !== ERRNO_CODES.EOPNOTSUPP) throw e;
            }
          }
        },connect:function (sock, addr, port) {
          if (sock.server) {
            throw new FS.ErrnoError(ERRNO_CODS.EOPNOTSUPP);
          }
  
          // TODO autobind
          // if (!sock.addr && sock.type == 2) {
          // }
  
          // early out if we're already connected / in the middle of connecting
          if (typeof sock.daddr !== 'undefined' && typeof sock.dport !== 'undefined') {
            var dest = SOCKFS.websocket_sock_ops.getPeer(sock, sock.daddr, sock.dport);
            if (dest) {
              if (dest.socket.readyState === dest.socket.CONNECTING) {
                throw new FS.ErrnoError(ERRNO_CODES.EALREADY);
              } else {
                throw new FS.ErrnoError(ERRNO_CODES.EISCONN);
              }
            }
          }
  
          // add the socket to our peer list and set our
          // destination address / port to match
          var peer = SOCKFS.websocket_sock_ops.createPeer(sock, addr, port);
          sock.daddr = peer.addr;
          sock.dport = peer.port;
  
          // always "fail" in non-blocking mode
          throw new FS.ErrnoError(ERRNO_CODES.EINPROGRESS);
        },listen:function (sock, backlog) {
          if (!ENVIRONMENT_IS_NODE) {
            throw new FS.ErrnoError(ERRNO_CODES.EOPNOTSUPP);
          }
          if (sock.server) {
             throw new FS.ErrnoError(ERRNO_CODES.EINVAL);  // already listening
          }
          var WebSocketServer = require('ws').Server;
          var host = sock.saddr;
          sock.server = new WebSocketServer({
            host: host,
            port: sock.sport
            // TODO support backlog
          });
  
          sock.server.on('connection', function(ws) {
            if (sock.type === 1) {
              var newsock = SOCKFS.createSocket(sock.family, sock.type, sock.protocol);
  
              // create a peer on the new socket
              var peer = SOCKFS.websocket_sock_ops.createPeer(newsock, ws);
              newsock.daddr = peer.addr;
              newsock.dport = peer.port;
  
              // push to queue for accept to pick up
              sock.pending.push(newsock);
            } else {
              // create a peer on the listen socket so calling sendto
              // with the listen socket and an address will resolve
              // to the correct client
              SOCKFS.websocket_sock_ops.createPeer(sock, ws);
            }
          });
          sock.server.on('closed', function() {
            sock.server = null;
          });
          sock.server.on('error', function() {
            // don't throw
          });
        },accept:function (listensock) {
          if (!listensock.server) {
            throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
          }
          var newsock = listensock.pending.shift();
          newsock.stream.flags = listensock.stream.flags;
          return newsock;
        },getname:function (sock, peer) {
          var addr, port;
          if (peer) {
            if (sock.daddr === undefined || sock.dport === undefined) {
              throw new FS.ErrnoError(ERRNO_CODES.ENOTCONN);
            }
            addr = sock.daddr;
            port = sock.dport;
          } else {
            // TODO saddr and sport will be set for bind()'d UDP sockets, but what
            // should we be returning for TCP sockets that've been connect()'d?
            addr = sock.saddr || 0;
            port = sock.sport || 0;
          }
          return { addr: addr, port: port };
        },sendmsg:function (sock, buffer, offset, length, addr, port) {
          if (sock.type === 2) {
            // connection-less sockets will honor the message address,
            // and otherwise fall back to the bound destination address
            if (addr === undefined || port === undefined) {
              addr = sock.daddr;
              port = sock.dport;
            }
            // if there was no address to fall back to, error out
            if (addr === undefined || port === undefined) {
              throw new FS.ErrnoError(ERRNO_CODES.EDESTADDRREQ);
            }
          } else {
            // connection-based sockets will only use the bound
            addr = sock.daddr;
            port = sock.dport;
          }
  
          // find the peer for the destination address
          var dest = SOCKFS.websocket_sock_ops.getPeer(sock, addr, port);
  
          // early out if not connected with a connection-based socket
          if (sock.type === 1) {
            if (!dest || dest.socket.readyState === dest.socket.CLOSING || dest.socket.readyState === dest.socket.CLOSED) {
              throw new FS.ErrnoError(ERRNO_CODES.ENOTCONN);
            } else if (dest.socket.readyState === dest.socket.CONNECTING) {
              throw new FS.ErrnoError(ERRNO_CODES.EAGAIN);
            }
          }
  
          // create a copy of the incoming data to send, as the WebSocket API
          // doesn't work entirely with an ArrayBufferView, it'll just send
          // the entire underlying buffer
          var data;
          if (buffer instanceof Array || buffer instanceof ArrayBuffer) {
            data = buffer.slice(offset, offset + length);
          } else {  // ArrayBufferView
            data = buffer.buffer.slice(buffer.byteOffset + offset, buffer.byteOffset + offset + length);
          }
  
          // if we're emulating a connection-less dgram socket and don't have
          // a cached connection, queue the buffer to send upon connect and
          // lie, saying the data was sent now.
          if (sock.type === 2) {
            if (!dest || dest.socket.readyState !== dest.socket.OPEN) {
              // if we're not connected, open a new connection
              if (!dest || dest.socket.readyState === dest.socket.CLOSING || dest.socket.readyState === dest.socket.CLOSED) {
                dest = SOCKFS.websocket_sock_ops.createPeer(sock, addr, port);
              }
              dest.dgram_send_queue.push(data);
              return length;
            }
          }
  
          try {
            // send the actual data
            dest.socket.send(data);
            return length;
          } catch (e) {
            throw new FS.ErrnoError(ERRNO_CODES.EINVAL);
          }
        },recvmsg:function (sock, length) {
          // http://pubs.opengroup.org/onlinepubs/7908799/xns/recvmsg.html
          if (sock.type === 1 && sock.server) {
            // tcp servers should not be recv()'ing on the listen socket
            throw new FS.ErrnoError(ERRNO_CODES.ENOTCONN);
          }
  
          var queued = sock.recv_queue.shift();
          if (!queued) {
            if (sock.type === 1) {
              var dest = SOCKFS.websocket_sock_ops.getPeer(sock, sock.daddr, sock.dport);
  
              if (!dest) {
                // if we have a destination address but are not connected, error out
                throw new FS.ErrnoError(ERRNO_CODES.ENOTCONN);
              }
              else if (dest.socket.readyState === dest.socket.CLOSING || dest.socket.readyState === dest.socket.CLOSED) {
                // return null if the socket has closed
                return null;
              }
              else {
                // else, our socket is in a valid state but truly has nothing available
                throw new FS.ErrnoError(ERRNO_CODES.EAGAIN);
              }
            } else {
              throw new FS.ErrnoError(ERRNO_CODES.EAGAIN);
            }
          }
  
          // queued.data will be an ArrayBuffer if it's unadulterated, but if it's
          // requeued TCP data it'll be an ArrayBufferView
          var queuedLength = queued.data.byteLength || queued.data.length;
          var queuedOffset = queued.data.byteOffset || 0;
          var queuedBuffer = queued.data.buffer || queued.data;
          var bytesRead = Math.min(length, queuedLength);
          var res = {
            buffer: new Uint8Array(queuedBuffer, queuedOffset, bytesRead),
            addr: queued.addr,
            port: queued.port
          };
  
  
          // push back any unread data for TCP connections
          if (sock.type === 1 && bytesRead < queuedLength) {
            var bytesRemaining = queuedLength - bytesRead;
            queued.data = new Uint8Array(queuedBuffer, queuedOffset + bytesRead, bytesRemaining);
            sock.recv_queue.unshift(queued);
          }
  
          return res;
        }}};function _send(fd, buf, len, flags) {
      var sock = SOCKFS.getSocket(fd);
      if (!sock) {
        ___setErrNo(ERRNO_CODES.EBADF);
        return -1;
      }
      // TODO honor flags
      return _write(fd, buf, len);
    }
  
  function _pwrite(fildes, buf, nbyte, offset) {
      // ssize_t pwrite(int fildes, const void *buf, size_t nbyte, off_t offset);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/write.html
      var stream = FS.getStream(fildes);
      if (!stream) {
        ___setErrNo(ERRNO_CODES.EBADF);
        return -1;
      }
      try {
        var slab = HEAP8;
        return FS.write(stream, slab, buf, nbyte, offset);
      } catch (e) {
        FS.handleFSError(e);
        return -1;
      }
    }function _write(fildes, buf, nbyte) {
      // ssize_t write(int fildes, const void *buf, size_t nbyte);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/write.html
      var stream = FS.getStream(fildes);
      if (!stream) {
        ___setErrNo(ERRNO_CODES.EBADF);
        return -1;
      }
  
  
      try {
        var slab = HEAP8;
        return FS.write(stream, slab, buf, nbyte);
      } catch (e) {
        FS.handleFSError(e);
        return -1;
      }
    }
  
  function _fileno(stream) {
      // int fileno(FILE *stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fileno.html
      stream = FS.getStreamFromPtr(stream);
      if (!stream) return -1;
      return stream.fd;
    }function _fwrite(ptr, size, nitems, stream) {
      // size_t fwrite(const void *restrict ptr, size_t size, size_t nitems, FILE *restrict stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fwrite.html
      var bytesToWrite = nitems * size;
      if (bytesToWrite == 0) return 0;
      var fd = _fileno(stream);
      var bytesWritten = _write(fd, ptr, bytesToWrite);
      if (bytesWritten == -1) {
        var streamObj = FS.getStreamFromPtr(stream);
        if (streamObj) streamObj.error = true;
        return 0;
      } else {
        return Math.floor(bytesWritten / size);
      }
    }

  
  function _strtoll(str, endptr, base) {
      return __parseInt64(str, endptr, base, '-9223372036854775808', '9223372036854775807');  // LLONG_MIN, LLONG_MAX.
    }function _strtoll_l(str, endptr, base) {
      return _strtoll(str, endptr, base); // no locale support yet
    }


  function _sysconf(name) {
      // long sysconf(int name);
      // http://pubs.opengroup.org/onlinepubs/009695399/functions/sysconf.html
      switch(name) {
        case 30: return PAGE_SIZE;
        case 132:
        case 133:
        case 12:
        case 137:
        case 138:
        case 15:
        case 235:
        case 16:
        case 17:
        case 18:
        case 19:
        case 20:
        case 149:
        case 13:
        case 10:
        case 236:
        case 153:
        case 9:
        case 21:
        case 22:
        case 159:
        case 154:
        case 14:
        case 77:
        case 78:
        case 139:
        case 80:
        case 81:
        case 79:
        case 82:
        case 68:
        case 67:
        case 164:
        case 11:
        case 29:
        case 47:
        case 48:
        case 95:
        case 52:
        case 51:
        case 46:
          return 200809;
        case 27:
        case 246:
        case 127:
        case 128:
        case 23:
        case 24:
        case 160:
        case 161:
        case 181:
        case 182:
        case 242:
        case 183:
        case 184:
        case 243:
        case 244:
        case 245:
        case 165:
        case 178:
        case 179:
        case 49:
        case 50:
        case 168:
        case 169:
        case 175:
        case 170:
        case 171:
        case 172:
        case 97:
        case 76:
        case 32:
        case 173:
        case 35:
          return -1;
        case 176:
        case 177:
        case 7:
        case 155:
        case 8:
        case 157:
        case 125:
        case 126:
        case 92:
        case 93:
        case 129:
        case 130:
        case 131:
        case 94:
        case 91:
          return 1;
        case 74:
        case 60:
        case 69:
        case 70:
        case 4:
          return 1024;
        case 31:
        case 42:
        case 72:
          return 32;
        case 87:
        case 26:
        case 33:
          return 2147483647;
        case 34:
        case 1:
          return 47839;
        case 38:
        case 36:
          return 99;
        case 43:
        case 37:
          return 2048;
        case 0: return 2097152;
        case 3: return 65536;
        case 28: return 32768;
        case 44: return 32767;
        case 75: return 16384;
        case 39: return 1000;
        case 89: return 700;
        case 71: return 256;
        case 40: return 255;
        case 2: return 100;
        case 180: return 64;
        case 25: return 20;
        case 5: return 16;
        case 6: return 6;
        case 73: return 4;
        case 84: return 1;
      }
      ___setErrNo(ERRNO_CODES.EINVAL);
      return -1;
    }

  function _pthread_cond_broadcast() {
      return 0;
    }

  
  
  function __exit(status) {
      // void _exit(int status);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/exit.html
      Module['exit'](status);
    }function _exit(status) {
      __exit(status);
    }function __ZSt9terminatev() {
      _exit(-1234);
    }

  function _emscripten_asm_const_double(code) {
      var args = Array.prototype.slice.call(arguments, 1);
      return +Runtime.getAsmConst(code, args.length).apply(null, args);
    }

  function _pthread_mutex_unlock() {}

  
  function _isxdigit(chr) {
      return (chr >= 48 && chr <= 57) ||
             (chr >= 97 && chr <= 102) ||
             (chr >= 65 && chr <= 70);
    }function _isxdigit_l(chr) {
      return _isxdigit(chr); // no locale support yet
    }

  
  function _emscripten_memcpy_big(dest, src, num) {
      HEAPU8.set(HEAPU8.subarray(src, src+num), dest);
      return dest;
    } 
  Module["_memcpy"] = _memcpy;

  function _sbrk(bytes) {
      // Implement a Linux-like 'memory area' for our 'process'.
      // Changes the size of the memory area by |bytes|; returns the
      // address of the previous top ('break') of the memory area
      // We control the "dynamic" memory - DYNAMIC_BASE to DYNAMICTOP
      var self = _sbrk;
      if (!self.called) {
        DYNAMICTOP = alignMemoryPage(DYNAMICTOP); // make sure we start out aligned
        self.called = true;
        assert(Runtime.dynamicAlloc);
        self.alloc = Runtime.dynamicAlloc;
        Runtime.dynamicAlloc = function() { abort('cannot dynamically allocate, sbrk now has control') };
      }
      var ret = DYNAMICTOP;
      if (bytes != 0) self.alloc(bytes);
      return ret;  // Previous break location.
    }


  function _newlocale(mask, locale, base) {
      return _malloc(4);
    }

   
  Module["_memmove"] = _memmove;

  function ___errno_location() {
      return ___errno_state;
    }

  var _BItoD=true;

  function _catclose(catd) {
      // int catclose (nl_catd catd)
      return 0;
    }

  function _fmod(x, y) {
      return x % y;
    }

  function ___cxa_guard_release() {}

  function _ungetc(c, stream) {
      // int ungetc(int c, FILE *stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/ungetc.html
      stream = FS.getStreamFromPtr(stream);
      if (!stream) {
        return -1;
      }
      if (c === -1) {
        // do nothing for EOF character
        return c;
      }
      c = unSign(c & 0xFF);
      stream.ungotten.push(c);
      stream.eof = false;
      return c;
    }

  function _uselocale(locale) {
      return 0;
    }

  function _vsnprintf(s, n, format, va_arg) {
      return _snprintf(s, n, format, HEAP32[((va_arg)>>2)]);
    }

  
  
  
  function __getFloat(text) {
      return /^[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?/.exec(text);
    }function __scanString(format, get, unget, varargs) {
      if (!__scanString.whiteSpace) {
        __scanString.whiteSpace = {};
        __scanString.whiteSpace[32] = 1;
        __scanString.whiteSpace[9] = 1;
        __scanString.whiteSpace[10] = 1;
        __scanString.whiteSpace[11] = 1;
        __scanString.whiteSpace[12] = 1;
        __scanString.whiteSpace[13] = 1;
      }
      // Supports %x, %4x, %d.%d, %lld, %s, %f, %lf.
      // TODO: Support all format specifiers.
      format = Pointer_stringify(format);
      var soFar = 0;
      if (format.indexOf('%n') >= 0) {
        // need to track soFar
        var _get = get;
        get = function get() {
          soFar++;
          return _get();
        }
        var _unget = unget;
        unget = function unget() {
          soFar--;
          return _unget();
        }
      }
      var formatIndex = 0;
      var argsi = 0;
      var fields = 0;
      var argIndex = 0;
      var next;
  
      mainLoop:
      for (var formatIndex = 0; formatIndex < format.length;) {
        if (format[formatIndex] === '%' && format[formatIndex+1] == 'n') {
          var argPtr = HEAP32[(((varargs)+(argIndex))>>2)];
          argIndex += Runtime.getAlignSize('void*', null, true);
          HEAP32[((argPtr)>>2)]=soFar;
          formatIndex += 2;
          continue;
        }
  
        if (format[formatIndex] === '%') {
          var nextC = format.indexOf('c', formatIndex+1);
          if (nextC > 0) {
            var maxx = 1;
            if (nextC > formatIndex+1) {
              var sub = format.substring(formatIndex+1, nextC);
              maxx = parseInt(sub);
              if (maxx != sub) maxx = 0;
            }
            if (maxx) {
              var argPtr = HEAP32[(((varargs)+(argIndex))>>2)];
              argIndex += Runtime.getAlignSize('void*', null, true);
              fields++;
              for (var i = 0; i < maxx; i++) {
                next = get();
                HEAP8[((argPtr++)|0)]=next;
                if (next === 0) return i > 0 ? fields : fields-1; // we failed to read the full length of this field
              }
              formatIndex += nextC - formatIndex + 1;
              continue;
            }
          }
        }
  
        // handle %[...]
        if (format[formatIndex] === '%' && format.indexOf('[', formatIndex+1) > 0) {
          var match = /\%([0-9]*)\[(\^)?(\]?[^\]]*)\]/.exec(format.substring(formatIndex));
          if (match) {
            var maxNumCharacters = parseInt(match[1]) || Infinity;
            var negateScanList = (match[2] === '^');
            var scanList = match[3];
  
            // expand "middle" dashs into character sets
            var middleDashMatch;
            while ((middleDashMatch = /([^\-])\-([^\-])/.exec(scanList))) {
              var rangeStartCharCode = middleDashMatch[1].charCodeAt(0);
              var rangeEndCharCode = middleDashMatch[2].charCodeAt(0);
              for (var expanded = ''; rangeStartCharCode <= rangeEndCharCode; expanded += String.fromCharCode(rangeStartCharCode++));
              scanList = scanList.replace(middleDashMatch[1] + '-' + middleDashMatch[2], expanded);
            }
  
            var argPtr = HEAP32[(((varargs)+(argIndex))>>2)];
            argIndex += Runtime.getAlignSize('void*', null, true);
            fields++;
  
            for (var i = 0; i < maxNumCharacters; i++) {
              next = get();
              if (negateScanList) {
                if (scanList.indexOf(String.fromCharCode(next)) < 0) {
                  HEAP8[((argPtr++)|0)]=next;
                } else {
                  unget();
                  break;
                }
              } else {
                if (scanList.indexOf(String.fromCharCode(next)) >= 0) {
                  HEAP8[((argPtr++)|0)]=next;
                } else {
                  unget();
                  break;
                }
              }
            }
  
            // write out null-terminating character
            HEAP8[((argPtr++)|0)]=0;
            formatIndex += match[0].length;
            
            continue;
          }
        }      
        // remove whitespace
        while (1) {
          next = get();
          if (next == 0) return fields;
          if (!(next in __scanString.whiteSpace)) break;
        }
        unget();
  
        if (format[formatIndex] === '%') {
          formatIndex++;
          var suppressAssignment = false;
          if (format[formatIndex] == '*') {
            suppressAssignment = true;
            formatIndex++;
          }
          var maxSpecifierStart = formatIndex;
          while (format[formatIndex].charCodeAt(0) >= 48 &&
                 format[formatIndex].charCodeAt(0) <= 57) {
            formatIndex++;
          }
          var max_;
          if (formatIndex != maxSpecifierStart) {
            max_ = parseInt(format.slice(maxSpecifierStart, formatIndex), 10);
          }
          var long_ = false;
          var half = false;
          var longLong = false;
          if (format[formatIndex] == 'l') {
            long_ = true;
            formatIndex++;
            if (format[formatIndex] == 'l') {
              longLong = true;
              formatIndex++;
            }
          } else if (format[formatIndex] == 'h') {
            half = true;
            formatIndex++;
          }
          var type = format[formatIndex];
          formatIndex++;
          var curr = 0;
          var buffer = [];
          // Read characters according to the format. floats are trickier, they may be in an unfloat state in the middle, then be a valid float later
          if (type == 'f' || type == 'e' || type == 'g' ||
              type == 'F' || type == 'E' || type == 'G') {
            next = get();
            while (next > 0 && (!(next in __scanString.whiteSpace)))  {
              buffer.push(String.fromCharCode(next));
              next = get();
            }
            var m = __getFloat(buffer.join(''));
            var last = m ? m[0].length : 0;
            for (var i = 0; i < buffer.length - last + 1; i++) {
              unget();
            }
            buffer.length = last;
          } else {
            next = get();
            var first = true;
            
            // Strip the optional 0x prefix for %x.
            if ((type == 'x' || type == 'X') && (next == 48)) {
              var peek = get();
              if (peek == 120 || peek == 88) {
                next = get();
              } else {
                unget();
              }
            }
            
            while ((curr < max_ || isNaN(max_)) && next > 0) {
              if (!(next in __scanString.whiteSpace) && // stop on whitespace
                  (type == 's' ||
                   ((type === 'd' || type == 'u' || type == 'i') && ((next >= 48 && next <= 57) ||
                                                                     (first && next == 45))) ||
                   ((type === 'x' || type === 'X') && (next >= 48 && next <= 57 ||
                                     next >= 97 && next <= 102 ||
                                     next >= 65 && next <= 70))) &&
                  (formatIndex >= format.length || next !== format[formatIndex].charCodeAt(0))) { // Stop when we read something that is coming up
                buffer.push(String.fromCharCode(next));
                next = get();
                curr++;
                first = false;
              } else {
                break;
              }
            }
            unget();
          }
          if (buffer.length === 0) return 0;  // Failure.
          if (suppressAssignment) continue;
  
          var text = buffer.join('');
          var argPtr = HEAP32[(((varargs)+(argIndex))>>2)];
          argIndex += Runtime.getAlignSize('void*', null, true);
          switch (type) {
            case 'd': case 'u': case 'i':
              if (half) {
                HEAP16[((argPtr)>>1)]=parseInt(text, 10);
              } else if (longLong) {
                (tempI64 = [parseInt(text, 10)>>>0,(tempDouble=parseInt(text, 10),(+(Math_abs(tempDouble))) >= (+1) ? (tempDouble > (+0) ? ((Math_min((+(Math_floor((tempDouble)/(+4294967296)))), (+4294967295)))|0)>>>0 : (~~((+(Math_ceil((tempDouble - +(((~~(tempDouble)))>>>0))/(+4294967296))))))>>>0) : 0)],HEAP32[((argPtr)>>2)]=tempI64[0],HEAP32[(((argPtr)+(4))>>2)]=tempI64[1]);
              } else {
                HEAP32[((argPtr)>>2)]=parseInt(text, 10);
              }
              break;
            case 'X':
            case 'x':
              HEAP32[((argPtr)>>2)]=parseInt(text, 16);
              break;
            case 'F':
            case 'f':
            case 'E':
            case 'e':
            case 'G':
            case 'g':
            case 'E':
              // fallthrough intended
              if (long_) {
                HEAPF64[((argPtr)>>3)]=parseFloat(text);
              } else {
                HEAPF32[((argPtr)>>2)]=parseFloat(text);
              }
              break;
            case 's':
              var array = intArrayFromString(text);
              for (var j = 0; j < array.length; j++) {
                HEAP8[(((argPtr)+(j))|0)]=array[j];
              }
              break;
          }
          fields++;
        } else if (format[formatIndex].charCodeAt(0) in __scanString.whiteSpace) {
          next = get();
          while (next in __scanString.whiteSpace) {
            if (next <= 0) break mainLoop;  // End of input.
            next = get();
          }
          unget(next);
          formatIndex++;
        } else {
          // Not a specifier.
          next = get();
          if (format[formatIndex].charCodeAt(0) !== next) {
            unget(next);
            break mainLoop;
          }
          formatIndex++;
        }
      }
      return fields;
    }function _sscanf(s, format, varargs) {
      // int sscanf(const char *restrict s, const char *restrict format, ... );
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/scanf.html
      var index = 0;
      function get() { return HEAP8[(((s)+(index++))|0)]; };
      function unget() { index--; };
      return __scanString(format, get, unget, varargs);
    }function _vsscanf(s, format, va_arg) {
      return _sscanf(s, format, HEAP32[((va_arg)>>2)]);
    }


  function ___assert_fail(condition, filename, line, func) {
      ABORT = true;
      throw 'Assertion failed: ' + Pointer_stringify(condition) + ', at: ' + [filename ? Pointer_stringify(filename) : 'unknown filename', line, func ? Pointer_stringify(func) : 'unknown function'] + ' at ' + stackTrace();
    }

  
  
  
  function _recv(fd, buf, len, flags) {
      var sock = SOCKFS.getSocket(fd);
      if (!sock) {
        ___setErrNo(ERRNO_CODES.EBADF);
        return -1;
      }
      // TODO honor flags
      return _read(fd, buf, len);
    }
  
  function _pread(fildes, buf, nbyte, offset) {
      // ssize_t pread(int fildes, void *buf, size_t nbyte, off_t offset);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/read.html
      var stream = FS.getStream(fildes);
      if (!stream) {
        ___setErrNo(ERRNO_CODES.EBADF);
        return -1;
      }
      try {
        var slab = HEAP8;
        return FS.read(stream, slab, buf, nbyte, offset);
      } catch (e) {
        FS.handleFSError(e);
        return -1;
      }
    }function _read(fildes, buf, nbyte) {
      // ssize_t read(int fildes, void *buf, size_t nbyte);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/read.html
      var stream = FS.getStream(fildes);
      if (!stream) {
        ___setErrNo(ERRNO_CODES.EBADF);
        return -1;
      }
  
  
      try {
        var slab = HEAP8;
        return FS.read(stream, slab, buf, nbyte);
      } catch (e) {
        FS.handleFSError(e);
        return -1;
      }
    }function _fread(ptr, size, nitems, stream) {
      // size_t fread(void *restrict ptr, size_t size, size_t nitems, FILE *restrict stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fread.html
      var bytesToRead = nitems * size;
      if (bytesToRead == 0) {
        return 0;
      }
      var bytesRead = 0;
      var streamObj = FS.getStreamFromPtr(stream);
      if (!streamObj) {
        ___setErrNo(ERRNO_CODES.EBADF);
        return 0;
      }
      while (streamObj.ungotten.length && bytesToRead > 0) {
        HEAP8[((ptr++)|0)]=streamObj.ungotten.pop();
        bytesToRead--;
        bytesRead++;
      }
      var err = _read(streamObj.fd, ptr, bytesToRead);
      if (err == -1) {
        if (streamObj) streamObj.error = true;
        return 0;
      }
      bytesRead += err;
      if (bytesRead < bytesToRead) streamObj.eof = true;
      return Math.floor(bytesRead / size);
    }function _fgetc(stream) {
      // int fgetc(FILE *stream);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/fgetc.html
      var streamObj = FS.getStreamFromPtr(stream);
      if (!streamObj) return -1;
      if (streamObj.eof || streamObj.error) return -1;
      var ret = _fread(_fgetc.ret, 1, 1, stream);
      if (ret == 0) {
        return -1;
      } else if (ret == -1) {
        streamObj.error = true;
        return -1;
      } else {
        return HEAPU8[((_fgetc.ret)|0)];
      }
    }

   
  Module["_memset"] = _memset;

  
  function _strerror_r(errnum, strerrbuf, buflen) {
      if (errnum in ERRNO_MESSAGES) {
        if (ERRNO_MESSAGES[errnum].length > buflen - 1) {
          return ___setErrNo(ERRNO_CODES.ERANGE);
        } else {
          var msg = ERRNO_MESSAGES[errnum];
          writeAsciiToMemory(msg, strerrbuf);
          return 0;
        }
      } else {
        return ___setErrNo(ERRNO_CODES.EINVAL);
      }
    }function _strerror(errnum) {
      if (!_strerror.buffer) _strerror.buffer = _malloc(256);
      _strerror_r(errnum, _strerror.buffer, 256);
      return _strerror.buffer;
    }

  
  
  function __isLeapYear(year) {
        return year%4 === 0 && (year%100 !== 0 || year%400 === 0);
    }
  
  function __arraySum(array, index) {
      var sum = 0;
      for (var i = 0; i <= index; sum += array[i++]);
      return sum;
    }
  
  
  var __MONTH_DAYS_LEAP=[31,29,31,30,31,30,31,31,30,31,30,31];
  
  var __MONTH_DAYS_REGULAR=[31,28,31,30,31,30,31,31,30,31,30,31];function __addDays(date, days) {
      var newDate = new Date(date.getTime());
      while(days > 0) {
        var leap = __isLeapYear(newDate.getFullYear());
        var currentMonth = newDate.getMonth();
        var daysInCurrentMonth = (leap ? __MONTH_DAYS_LEAP : __MONTH_DAYS_REGULAR)[currentMonth];
  
        if (days > daysInCurrentMonth-newDate.getDate()) {
          // we spill over to next month
          days -= (daysInCurrentMonth-newDate.getDate()+1);
          newDate.setDate(1);
          if (currentMonth < 11) {
            newDate.setMonth(currentMonth+1)
          } else {
            newDate.setMonth(0);
            newDate.setFullYear(newDate.getFullYear()+1);
          }
        } else {
          // we stay in current month 
          newDate.setDate(newDate.getDate()+days);
          return newDate;
        }
      }
  
      return newDate;
    }function _strftime(s, maxsize, format, tm) {
      // size_t strftime(char *restrict s, size_t maxsize, const char *restrict format, const struct tm *restrict timeptr);
      // http://pubs.opengroup.org/onlinepubs/009695399/functions/strftime.html
      
      var date = {
        tm_sec: HEAP32[((tm)>>2)],
        tm_min: HEAP32[(((tm)+(4))>>2)],
        tm_hour: HEAP32[(((tm)+(8))>>2)],
        tm_mday: HEAP32[(((tm)+(12))>>2)],
        tm_mon: HEAP32[(((tm)+(16))>>2)],
        tm_year: HEAP32[(((tm)+(20))>>2)],
        tm_wday: HEAP32[(((tm)+(24))>>2)],
        tm_yday: HEAP32[(((tm)+(28))>>2)],
        tm_isdst: HEAP32[(((tm)+(32))>>2)]
      };
  
      var pattern = Pointer_stringify(format);
  
      // expand format
      var EXPANSION_RULES_1 = {
        '%c': '%a %b %d %H:%M:%S %Y',     // Replaced by the locale's appropriate date and time representation - e.g., Mon Aug  3 14:02:01 2013
        '%D': '%m/%d/%y',                 // Equivalent to %m / %d / %y
        '%F': '%Y-%m-%d',                 // Equivalent to %Y - %m - %d
        '%h': '%b',                       // Equivalent to %b
        '%r': '%I:%M:%S %p',              // Replaced by the time in a.m. and p.m. notation
        '%R': '%H:%M',                    // Replaced by the time in 24-hour notation
        '%T': '%H:%M:%S',                 // Replaced by the time
        '%x': '%m/%d/%y',                 // Replaced by the locale's appropriate date representation
        '%X': '%H:%M:%S',                 // Replaced by the locale's appropriate date representation
      };
      for (var rule in EXPANSION_RULES_1) {
        pattern = pattern.replace(new RegExp(rule, 'g'), EXPANSION_RULES_1[rule]);
      }
  
      var WEEKDAYS = ['Sunday', 'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday'];
      var MONTHS = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'];
  
      function leadingSomething(value, digits, character) {
        var str = typeof value === 'number' ? value.toString() : (value || '');
        while (str.length < digits) {
          str = character[0]+str;
        }
        return str;
      };
  
      function leadingNulls(value, digits) {
        return leadingSomething(value, digits, '0');
      };
  
      function compareByDay(date1, date2) {
        function sgn(value) {
          return value < 0 ? -1 : (value > 0 ? 1 : 0);
        };
  
        var compare;
        if ((compare = sgn(date1.getFullYear()-date2.getFullYear())) === 0) {
          if ((compare = sgn(date1.getMonth()-date2.getMonth())) === 0) {
            compare = sgn(date1.getDate()-date2.getDate());
          }
        }
        return compare;
      };
  
      function getFirstWeekStartDate(janFourth) {
          switch (janFourth.getDay()) {
            case 0: // Sunday
              return new Date(janFourth.getFullYear()-1, 11, 29);
            case 1: // Monday
              return janFourth;
            case 2: // Tuesday
              return new Date(janFourth.getFullYear(), 0, 3);
            case 3: // Wednesday
              return new Date(janFourth.getFullYear(), 0, 2);
            case 4: // Thursday
              return new Date(janFourth.getFullYear(), 0, 1);
            case 5: // Friday
              return new Date(janFourth.getFullYear()-1, 11, 31);
            case 6: // Saturday
              return new Date(janFourth.getFullYear()-1, 11, 30);
          }
      };
  
      function getWeekBasedYear(date) {
          var thisDate = __addDays(new Date(date.tm_year+1900, 0, 1), date.tm_yday);
  
          var janFourthThisYear = new Date(thisDate.getFullYear(), 0, 4);
          var janFourthNextYear = new Date(thisDate.getFullYear()+1, 0, 4);
  
          var firstWeekStartThisYear = getFirstWeekStartDate(janFourthThisYear);
          var firstWeekStartNextYear = getFirstWeekStartDate(janFourthNextYear);
  
          if (compareByDay(firstWeekStartThisYear, thisDate) <= 0) {
            // this date is after the start of the first week of this year
            if (compareByDay(firstWeekStartNextYear, thisDate) <= 0) {
              return thisDate.getFullYear()+1;
            } else {
              return thisDate.getFullYear();
            }
          } else { 
            return thisDate.getFullYear()-1;
          }
      };
  
      var EXPANSION_RULES_2 = {
        '%a': function(date) {
          return WEEKDAYS[date.tm_wday].substring(0,3);
        },
        '%A': function(date) {
          return WEEKDAYS[date.tm_wday];
        },
        '%b': function(date) {
          return MONTHS[date.tm_mon].substring(0,3);
        },
        '%B': function(date) {
          return MONTHS[date.tm_mon];
        },
        '%C': function(date) {
          var year = date.tm_year+1900;
          return leadingNulls(Math.floor(year/100),2);
        },
        '%d': function(date) {
          return leadingNulls(date.tm_mday, 2);
        },
        '%e': function(date) {
          return leadingSomething(date.tm_mday, 2, ' ');
        },
        '%g': function(date) {
          // %g, %G, and %V give values according to the ISO 8601:2000 standard week-based year. 
          // In this system, weeks begin on a Monday and week 1 of the year is the week that includes 
          // January 4th, which is also the week that includes the first Thursday of the year, and 
          // is also the first week that contains at least four days in the year. 
          // If the first Monday of January is the 2nd, 3rd, or 4th, the preceding days are part of 
          // the last week of the preceding year; thus, for Saturday 2nd January 1999, 
          // %G is replaced by 1998 and %V is replaced by 53. If December 29th, 30th, 
          // or 31st is a Monday, it and any following days are part of week 1 of the following year. 
          // Thus, for Tuesday 30th December 1997, %G is replaced by 1998 and %V is replaced by 01.
          
          return getWeekBasedYear(date).toString().substring(2);
        },
        '%G': function(date) {
          return getWeekBasedYear(date);
        },
        '%H': function(date) {
          return leadingNulls(date.tm_hour, 2);
        },
        '%I': function(date) {
          return leadingNulls(date.tm_hour < 13 ? date.tm_hour : date.tm_hour-12, 2);
        },
        '%j': function(date) {
          // Day of the year (001-366)
          return leadingNulls(date.tm_mday+__arraySum(__isLeapYear(date.tm_year+1900) ? __MONTH_DAYS_LEAP : __MONTH_DAYS_REGULAR, date.tm_mon-1), 3);
        },
        '%m': function(date) {
          return leadingNulls(date.tm_mon+1, 2);
        },
        '%M': function(date) {
          return leadingNulls(date.tm_min, 2);
        },
        '%n': function() {
          return '\n';
        },
        '%p': function(date) {
          if (date.tm_hour > 0 && date.tm_hour < 13) {
            return 'AM';
          } else {
            return 'PM';
          }
        },
        '%S': function(date) {
          return leadingNulls(date.tm_sec, 2);
        },
        '%t': function() {
          return '\t';
        },
        '%u': function(date) {
          var day = new Date(date.tm_year+1900, date.tm_mon+1, date.tm_mday, 0, 0, 0, 0);
          return day.getDay() || 7;
        },
        '%U': function(date) {
          // Replaced by the week number of the year as a decimal number [00,53]. 
          // The first Sunday of January is the first day of week 1; 
          // days in the new year before this are in week 0. [ tm_year, tm_wday, tm_yday]
          var janFirst = new Date(date.tm_year+1900, 0, 1);
          var firstSunday = janFirst.getDay() === 0 ? janFirst : __addDays(janFirst, 7-janFirst.getDay());
          var endDate = new Date(date.tm_year+1900, date.tm_mon, date.tm_mday);
          
          // is target date after the first Sunday?
          if (compareByDay(firstSunday, endDate) < 0) {
            // calculate difference in days between first Sunday and endDate
            var februaryFirstUntilEndMonth = __arraySum(__isLeapYear(endDate.getFullYear()) ? __MONTH_DAYS_LEAP : __MONTH_DAYS_REGULAR, endDate.getMonth()-1)-31;
            var firstSundayUntilEndJanuary = 31-firstSunday.getDate();
            var days = firstSundayUntilEndJanuary+februaryFirstUntilEndMonth+endDate.getDate();
            return leadingNulls(Math.ceil(days/7), 2);
          }
  
          return compareByDay(firstSunday, janFirst) === 0 ? '01': '00';
        },
        '%V': function(date) {
          // Replaced by the week number of the year (Monday as the first day of the week) 
          // as a decimal number [01,53]. If the week containing 1 January has four 
          // or more days in the new year, then it is considered week 1. 
          // Otherwise, it is the last week of the previous year, and the next week is week 1. 
          // Both January 4th and the first Thursday of January are always in week 1. [ tm_year, tm_wday, tm_yday]
          var janFourthThisYear = new Date(date.tm_year+1900, 0, 4);
          var janFourthNextYear = new Date(date.tm_year+1901, 0, 4);
  
          var firstWeekStartThisYear = getFirstWeekStartDate(janFourthThisYear);
          var firstWeekStartNextYear = getFirstWeekStartDate(janFourthNextYear);
  
          var endDate = __addDays(new Date(date.tm_year+1900, 0, 1), date.tm_yday);
  
          if (compareByDay(endDate, firstWeekStartThisYear) < 0) {
            // if given date is before this years first week, then it belongs to the 53rd week of last year
            return '53';
          } 
  
          if (compareByDay(firstWeekStartNextYear, endDate) <= 0) {
            // if given date is after next years first week, then it belongs to the 01th week of next year
            return '01';
          }
  
          // given date is in between CW 01..53 of this calendar year
          var daysDifference;
          if (firstWeekStartThisYear.getFullYear() < date.tm_year+1900) {
            // first CW of this year starts last year
            daysDifference = date.tm_yday+32-firstWeekStartThisYear.getDate()
          } else {
            // first CW of this year starts this year
            daysDifference = date.tm_yday+1-firstWeekStartThisYear.getDate();
          }
          return leadingNulls(Math.ceil(daysDifference/7), 2);
        },
        '%w': function(date) {
          var day = new Date(date.tm_year+1900, date.tm_mon+1, date.tm_mday, 0, 0, 0, 0);
          return day.getDay();
        },
        '%W': function(date) {
          // Replaced by the week number of the year as a decimal number [00,53]. 
          // The first Monday of January is the first day of week 1; 
          // days in the new year before this are in week 0. [ tm_year, tm_wday, tm_yday]
          var janFirst = new Date(date.tm_year, 0, 1);
          var firstMonday = janFirst.getDay() === 1 ? janFirst : __addDays(janFirst, janFirst.getDay() === 0 ? 1 : 7-janFirst.getDay()+1);
          var endDate = new Date(date.tm_year+1900, date.tm_mon, date.tm_mday);
  
          // is target date after the first Monday?
          if (compareByDay(firstMonday, endDate) < 0) {
            var februaryFirstUntilEndMonth = __arraySum(__isLeapYear(endDate.getFullYear()) ? __MONTH_DAYS_LEAP : __MONTH_DAYS_REGULAR, endDate.getMonth()-1)-31;
            var firstMondayUntilEndJanuary = 31-firstMonday.getDate();
            var days = firstMondayUntilEndJanuary+februaryFirstUntilEndMonth+endDate.getDate();
            return leadingNulls(Math.ceil(days/7), 2);
          }
          return compareByDay(firstMonday, janFirst) === 0 ? '01': '00';
        },
        '%y': function(date) {
          // Replaced by the last two digits of the year as a decimal number [00,99]. [ tm_year]
          return (date.tm_year+1900).toString().substring(2);
        },
        '%Y': function(date) {
          // Replaced by the year as a decimal number (for example, 1997). [ tm_year]
          return date.tm_year+1900;
        },
        '%z': function(date) {
          // Replaced by the offset from UTC in the ISO 8601:2000 standard format ( +hhmm or -hhmm ),
          // or by no characters if no timezone is determinable. 
          // For example, "-0430" means 4 hours 30 minutes behind UTC (west of Greenwich). 
          // If tm_isdst is zero, the standard time offset is used. 
          // If tm_isdst is greater than zero, the daylight savings time offset is used. 
          // If tm_isdst is negative, no characters are returned. 
          // FIXME: we cannot determine time zone (or can we?)
          return '';
        },
        '%Z': function(date) {
          // Replaced by the timezone name or abbreviation, or by no bytes if no timezone information exists. [ tm_isdst]
          // FIXME: we cannot determine time zone (or can we?)
          return '';
        },
        '%%': function() {
          return '%';
        }
      };
      for (var rule in EXPANSION_RULES_2) {
        if (pattern.indexOf(rule) >= 0) {
          pattern = pattern.replace(new RegExp(rule, 'g'), EXPANSION_RULES_2[rule](date));
        }
      }
  
      var bytes = intArrayFromString(pattern, false);
      if (bytes.length > maxsize) {
        return 0;
      } 
  
      writeArrayToMemory(bytes, s);
      return bytes.length-1;
    }function _strftime_l(s, maxsize, format, tm) {
      return _strftime(s, maxsize, format, tm); // no locale support yet
    }

  function _abort() {
      Module['abort']();
    }


  function _pthread_cond_wait() {
      return 0;
    }

  
  function _isdigit(chr) {
      return chr >= 48 && chr <= 57;
    }function _isdigit_l(chr) {
      return _isdigit(chr); // no locale support yet
    }

  var _fabs=Math_abs;

  var _getc=_fgetc;

  var _sqrt=Math_sqrt;

  function _emscripten_asm_const_int(code) {
      var args = Array.prototype.slice.call(arguments, 1);
      return Runtime.getAsmConst(code, args.length).apply(null, args) | 0;
    }

  var Browser={mainLoop:{scheduler:null,method:"",shouldPause:false,paused:false,queue:[],pause:function () {
          Browser.mainLoop.shouldPause = true;
        },resume:function () {
          if (Browser.mainLoop.paused) {
            Browser.mainLoop.paused = false;
            Browser.mainLoop.scheduler();
          }
          Browser.mainLoop.shouldPause = false;
        },updateStatus:function () {
          if (Module['setStatus']) {
            var message = Module['statusMessage'] || 'Please wait...';
            var remaining = Browser.mainLoop.remainingBlockers;
            var expected = Browser.mainLoop.expectedBlockers;
            if (remaining) {
              if (remaining < expected) {
                Module['setStatus'](message + ' (' + (expected - remaining) + '/' + expected + ')');
              } else {
                Module['setStatus'](message);
              }
            } else {
              Module['setStatus']('');
            }
          }
        }},isFullScreen:false,pointerLock:false,moduleContextCreatedCallbacks:[],workers:[],init:function () {
        if (!Module["preloadPlugins"]) Module["preloadPlugins"] = []; // needs to exist even in workers
  
        if (Browser.initted || ENVIRONMENT_IS_WORKER) return;
        Browser.initted = true;
  
        try {
          new Blob();
          Browser.hasBlobConstructor = true;
        } catch(e) {
          Browser.hasBlobConstructor = false;
          console.log("warning: no blob constructor, cannot create blobs with mimetypes");
        }
        Browser.BlobBuilder = typeof MozBlobBuilder != "undefined" ? MozBlobBuilder : (typeof WebKitBlobBuilder != "undefined" ? WebKitBlobBuilder : (!Browser.hasBlobConstructor ? console.log("warning: no BlobBuilder") : null));
        Browser.URLObject = typeof window != "undefined" ? (window.URL ? window.URL : window.webkitURL) : undefined;
        if (!Module.noImageDecoding && typeof Browser.URLObject === 'undefined') {
          console.log("warning: Browser does not support creating object URLs. Built-in browser image decoding will not be available.");
          Module.noImageDecoding = true;
        }
  
        // Support for plugins that can process preloaded files. You can add more of these to
        // your app by creating and appending to Module.preloadPlugins.
        //
        // Each plugin is asked if it can handle a file based on the file's name. If it can,
        // it is given the file's raw data. When it is done, it calls a callback with the file's
        // (possibly modified) data. For example, a plugin might decompress a file, or it
        // might create some side data structure for use later (like an Image element, etc.).
  
        var imagePlugin = {};
        imagePlugin['canHandle'] = function imagePlugin_canHandle(name) {
          return !Module.noImageDecoding && /\.(jpg|jpeg|png|bmp)$/i.test(name);
        };
        imagePlugin['handle'] = function imagePlugin_handle(byteArray, name, onload, onerror) {
          var b = null;
          if (Browser.hasBlobConstructor) {
            try {
              b = new Blob([byteArray], { type: Browser.getMimetype(name) });
              if (b.size !== byteArray.length) { // Safari bug #118630
                // Safari's Blob can only take an ArrayBuffer
                b = new Blob([(new Uint8Array(byteArray)).buffer], { type: Browser.getMimetype(name) });
              }
            } catch(e) {
              Runtime.warnOnce('Blob constructor present but fails: ' + e + '; falling back to blob builder');
            }
          }
          if (!b) {
            var bb = new Browser.BlobBuilder();
            bb.append((new Uint8Array(byteArray)).buffer); // we need to pass a buffer, and must copy the array to get the right data range
            b = bb.getBlob();
          }
          var url = Browser.URLObject.createObjectURL(b);
          var img = new Image();
          img.onload = function img_onload() {
            assert(img.complete, 'Image ' + name + ' could not be decoded');
            var canvas = document.createElement('canvas');
            canvas.width = img.width;
            canvas.height = img.height;
            var ctx = canvas.getContext('2d');
            ctx.drawImage(img, 0, 0);
            Module["preloadedImages"][name] = canvas;
            Browser.URLObject.revokeObjectURL(url);
            if (onload) onload(byteArray);
          };
          img.onerror = function img_onerror(event) {
            console.log('Image ' + url + ' could not be decoded');
            if (onerror) onerror();
          };
          img.src = url;
        };
        Module['preloadPlugins'].push(imagePlugin);
  
        var audioPlugin = {};
        audioPlugin['canHandle'] = function audioPlugin_canHandle(name) {
          return !Module.noAudioDecoding && name.substr(-4) in { '.ogg': 1, '.wav': 1, '.mp3': 1 };
        };
        audioPlugin['handle'] = function audioPlugin_handle(byteArray, name, onload, onerror) {
          var done = false;
          function finish(audio) {
            if (done) return;
            done = true;
            Module["preloadedAudios"][name] = audio;
            if (onload) onload(byteArray);
          }
          function fail() {
            if (done) return;
            done = true;
            Module["preloadedAudios"][name] = new Audio(); // empty shim
            if (onerror) onerror();
          }
          if (Browser.hasBlobConstructor) {
            try {
              var b = new Blob([byteArray], { type: Browser.getMimetype(name) });
            } catch(e) {
              return fail();
            }
            var url = Browser.URLObject.createObjectURL(b); // XXX we never revoke this!
            var audio = new Audio();
            audio.addEventListener('canplaythrough', function() { finish(audio) }, false); // use addEventListener due to chromium bug 124926
            audio.onerror = function audio_onerror(event) {
              if (done) return;
              console.log('warning: browser could not fully decode audio ' + name + ', trying slower base64 approach');
              function encode64(data) {
                var BASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';
                var PAD = '=';
                var ret = '';
                var leftchar = 0;
                var leftbits = 0;
                for (var i = 0; i < data.length; i++) {
                  leftchar = (leftchar << 8) | data[i];
                  leftbits += 8;
                  while (leftbits >= 6) {
                    var curr = (leftchar >> (leftbits-6)) & 0x3f;
                    leftbits -= 6;
                    ret += BASE[curr];
                  }
                }
                if (leftbits == 2) {
                  ret += BASE[(leftchar&3) << 4];
                  ret += PAD + PAD;
                } else if (leftbits == 4) {
                  ret += BASE[(leftchar&0xf) << 2];
                  ret += PAD;
                }
                return ret;
              }
              audio.src = 'data:audio/x-' + name.substr(-3) + ';base64,' + encode64(byteArray);
              finish(audio); // we don't wait for confirmation this worked - but it's worth trying
            };
            audio.src = url;
            // workaround for chrome bug 124926 - we do not always get oncanplaythrough or onerror
            Browser.safeSetTimeout(function() {
              finish(audio); // try to use it even though it is not necessarily ready to play
            }, 10000);
          } else {
            return fail();
          }
        };
        Module['preloadPlugins'].push(audioPlugin);
  
        // Canvas event setup
  
        var canvas = Module['canvas'];
        
        // forced aspect ratio can be enabled by defining 'forcedAspectRatio' on Module
        // Module['forcedAspectRatio'] = 4 / 3;
        
        canvas.requestPointerLock = canvas['requestPointerLock'] ||
                                    canvas['mozRequestPointerLock'] ||
                                    canvas['webkitRequestPointerLock'] ||
                                    canvas['msRequestPointerLock'] ||
                                    function(){};
        canvas.exitPointerLock = document['exitPointerLock'] ||
                                 document['mozExitPointerLock'] ||
                                 document['webkitExitPointerLock'] ||
                                 document['msExitPointerLock'] ||
                                 function(){}; // no-op if function does not exist
        canvas.exitPointerLock = canvas.exitPointerLock.bind(document);
  
        function pointerLockChange() {
          Browser.pointerLock = document['pointerLockElement'] === canvas ||
                                document['mozPointerLockElement'] === canvas ||
                                document['webkitPointerLockElement'] === canvas ||
                                document['msPointerLockElement'] === canvas;
        }
  
        document.addEventListener('pointerlockchange', pointerLockChange, false);
        document.addEventListener('mozpointerlockchange', pointerLockChange, false);
        document.addEventListener('webkitpointerlockchange', pointerLockChange, false);
        document.addEventListener('mspointerlockchange', pointerLockChange, false);
  
        if (Module['elementPointerLock']) {
          canvas.addEventListener("click", function(ev) {
            if (!Browser.pointerLock && canvas.requestPointerLock) {
              canvas.requestPointerLock();
              ev.preventDefault();
            }
          }, false);
        }
      },createContext:function (canvas, useWebGL, setInModule, webGLContextAttributes) {
        var ctx;
        var errorInfo = '?';
        function onContextCreationError(event) {
          errorInfo = event.statusMessage || errorInfo;
        }
        try {
          if (useWebGL) {
            var contextAttributes = {
              antialias: false,
              alpha: false
            };
  
            if (webGLContextAttributes) {
              for (var attribute in webGLContextAttributes) {
                contextAttributes[attribute] = webGLContextAttributes[attribute];
              }
            }
  
  
            canvas.addEventListener('webglcontextcreationerror', onContextCreationError, false);
            try {
              ['experimental-webgl', 'webgl'].some(function(webglId) {
                return ctx = canvas.getContext(webglId, contextAttributes);
              });
            } finally {
              canvas.removeEventListener('webglcontextcreationerror', onContextCreationError, false);
            }
          } else {
            ctx = canvas.getContext('2d');
          }
          if (!ctx) throw ':(';
        } catch (e) {
          Module.print('Could not create canvas: ' + [errorInfo, e]);
          return null;
        }
        if (useWebGL) {
          // Set the background of the WebGL canvas to black
          canvas.style.backgroundColor = "black";
  
          // Warn on context loss
          canvas.addEventListener('webglcontextlost', function(event) {
            alert('WebGL context lost. You will need to reload the page.');
          }, false);
        }
        if (setInModule) {
          GLctx = Module.ctx = ctx;
          Module.useWebGL = useWebGL;
          Browser.moduleContextCreatedCallbacks.forEach(function(callback) { callback() });
          Browser.init();
        }
        return ctx;
      },destroyContext:function (canvas, useWebGL, setInModule) {},fullScreenHandlersInstalled:false,lockPointer:undefined,resizeCanvas:undefined,requestFullScreen:function (lockPointer, resizeCanvas) {
        Browser.lockPointer = lockPointer;
        Browser.resizeCanvas = resizeCanvas;
        if (typeof Browser.lockPointer === 'undefined') Browser.lockPointer = true;
        if (typeof Browser.resizeCanvas === 'undefined') Browser.resizeCanvas = false;
  
        var canvas = Module['canvas'];
        function fullScreenChange() {
          Browser.isFullScreen = false;
          var canvasContainer = canvas.parentNode;
          if ((document['webkitFullScreenElement'] || document['webkitFullscreenElement'] ||
               document['mozFullScreenElement'] || document['mozFullscreenElement'] ||
               document['fullScreenElement'] || document['fullscreenElement'] ||
               document['msFullScreenElement'] || document['msFullscreenElement'] ||
               document['webkitCurrentFullScreenElement']) === canvasContainer) {
            canvas.cancelFullScreen = document['cancelFullScreen'] ||
                                      document['mozCancelFullScreen'] ||
                                      document['webkitCancelFullScreen'] ||
                                      document['msExitFullscreen'] ||
                                      document['exitFullscreen'] ||
                                      function() {};
            canvas.cancelFullScreen = canvas.cancelFullScreen.bind(document);
            if (Browser.lockPointer) canvas.requestPointerLock();
            Browser.isFullScreen = true;
            if (Browser.resizeCanvas) Browser.setFullScreenCanvasSize();
          } else {
            
            // remove the full screen specific parent of the canvas again to restore the HTML structure from before going full screen
            canvasContainer.parentNode.insertBefore(canvas, canvasContainer);
            canvasContainer.parentNode.removeChild(canvasContainer);
            
            if (Browser.resizeCanvas) Browser.setWindowedCanvasSize();
          }
          if (Module['onFullScreen']) Module['onFullScreen'](Browser.isFullScreen);
          Browser.updateCanvasDimensions(canvas);
        }
  
        if (!Browser.fullScreenHandlersInstalled) {
          Browser.fullScreenHandlersInstalled = true;
          document.addEventListener('fullscreenchange', fullScreenChange, false);
          document.addEventListener('mozfullscreenchange', fullScreenChange, false);
          document.addEventListener('webkitfullscreenchange', fullScreenChange, false);
          document.addEventListener('MSFullscreenChange', fullScreenChange, false);
        }
  
        // create a new parent to ensure the canvas has no siblings. this allows browsers to optimize full screen performance when its parent is the full screen root
        var canvasContainer = document.createElement("div");
        canvas.parentNode.insertBefore(canvasContainer, canvas);
        canvasContainer.appendChild(canvas);
        
        // use parent of canvas as full screen root to allow aspect ratio correction (Firefox stretches the root to screen size)
        canvasContainer.requestFullScreen = canvasContainer['requestFullScreen'] ||
                                            canvasContainer['mozRequestFullScreen'] ||
                                            canvasContainer['msRequestFullscreen'] ||
                                           (canvasContainer['webkitRequestFullScreen'] ? function() { canvasContainer['webkitRequestFullScreen'](Element['ALLOW_KEYBOARD_INPUT']) } : null);
        canvasContainer.requestFullScreen();
      },requestAnimationFrame:function requestAnimationFrame(func) {
        if (typeof window === 'undefined') { // Provide fallback to setTimeout if window is undefined (e.g. in Node.js)
          setTimeout(func, 1000/60);
        } else {
          if (!window.requestAnimationFrame) {
            window.requestAnimationFrame = window['requestAnimationFrame'] ||
                                           window['mozRequestAnimationFrame'] ||
                                           window['webkitRequestAnimationFrame'] ||
                                           window['msRequestAnimationFrame'] ||
                                           window['oRequestAnimationFrame'] ||
                                           window['setTimeout'];
          }
          window.requestAnimationFrame(func);
        }
      },safeCallback:function (func) {
        return function() {
          if (!ABORT) return func.apply(null, arguments);
        };
      },safeRequestAnimationFrame:function (func) {
        return Browser.requestAnimationFrame(function() {
          if (!ABORT) func();
        });
      },safeSetTimeout:function (func, timeout) {
        return setTimeout(function() {
          if (!ABORT) func();
        }, timeout);
      },safeSetInterval:function (func, timeout) {
        return setInterval(function() {
          if (!ABORT) func();
        }, timeout);
      },getMimetype:function (name) {
        return {
          'jpg': 'image/jpeg',
          'jpeg': 'image/jpeg',
          'png': 'image/png',
          'bmp': 'image/bmp',
          'ogg': 'audio/ogg',
          'wav': 'audio/wav',
          'mp3': 'audio/mpeg'
        }[name.substr(name.lastIndexOf('.')+1)];
      },getUserMedia:function (func) {
        if(!window.getUserMedia) {
          window.getUserMedia = navigator['getUserMedia'] ||
                                navigator['mozGetUserMedia'];
        }
        window.getUserMedia(func);
      },getMovementX:function (event) {
        return event['movementX'] ||
               event['mozMovementX'] ||
               event['webkitMovementX'] ||
               0;
      },getMovementY:function (event) {
        return event['movementY'] ||
               event['mozMovementY'] ||
               event['webkitMovementY'] ||
               0;
      },getMouseWheelDelta:function (event) {
        return Math.max(-1, Math.min(1, event.type === 'DOMMouseScroll' ? event.detail : -event.wheelDelta));
      },mouseX:0,mouseY:0,mouseMovementX:0,mouseMovementY:0,touches:{},lastTouches:{},calculateMouseEvent:function (event) { // event should be mousemove, mousedown or mouseup
        if (Browser.pointerLock) {
          // When the pointer is locked, calculate the coordinates
          // based on the movement of the mouse.
          // Workaround for Firefox bug 764498
          if (event.type != 'mousemove' &&
              ('mozMovementX' in event)) {
            Browser.mouseMovementX = Browser.mouseMovementY = 0;
          } else {
            Browser.mouseMovementX = Browser.getMovementX(event);
            Browser.mouseMovementY = Browser.getMovementY(event);
          }
          
          // check if SDL is available
          if (typeof SDL != "undefined") {
          	Browser.mouseX = SDL.mouseX + Browser.mouseMovementX;
          	Browser.mouseY = SDL.mouseY + Browser.mouseMovementY;
          } else {
          	// just add the mouse delta to the current absolut mouse position
          	// FIXME: ideally this should be clamped against the canvas size and zero
          	Browser.mouseX += Browser.mouseMovementX;
          	Browser.mouseY += Browser.mouseMovementY;
          }        
        } else {
          // Otherwise, calculate the movement based on the changes
          // in the coordinates.
          var rect = Module["canvas"].getBoundingClientRect();
          var cw = Module["canvas"].width;
          var ch = Module["canvas"].height;
  
          // Neither .scrollX or .pageXOffset are defined in a spec, but
          // we prefer .scrollX because it is currently in a spec draft.
          // (see: http://www.w3.org/TR/2013/WD-cssom-view-20131217/)
          var scrollX = ((typeof window.scrollX !== 'undefined') ? window.scrollX : window.pageXOffset);
          var scrollY = ((typeof window.scrollY !== 'undefined') ? window.scrollY : window.pageYOffset);
  
          if (event.type === 'touchstart' || event.type === 'touchend' || event.type === 'touchmove') {
            var touch = event.touch;
            if (touch === undefined) {
              return; // the "touch" property is only defined in SDL
  
            }
            var adjustedX = touch.pageX - (scrollX + rect.left);
            var adjustedY = touch.pageY - (scrollY + rect.top);
  
            adjustedX = adjustedX * (cw / rect.width);
            adjustedY = adjustedY * (ch / rect.height);
  
            var coords = { x: adjustedX, y: adjustedY };
            
            if (event.type === 'touchstart') {
              Browser.lastTouches[touch.identifier] = coords;
              Browser.touches[touch.identifier] = coords;
            } else if (event.type === 'touchend' || event.type === 'touchmove') {
              Browser.lastTouches[touch.identifier] = Browser.touches[touch.identifier];
              Browser.touches[touch.identifier] = { x: adjustedX, y: adjustedY };
            } 
            return;
          }
  
          var x = event.pageX - (scrollX + rect.left);
          var y = event.pageY - (scrollY + rect.top);
  
          // the canvas might be CSS-scaled compared to its backbuffer;
          // SDL-using content will want mouse coordinates in terms
          // of backbuffer units.
          x = x * (cw / rect.width);
          y = y * (ch / rect.height);
  
          Browser.mouseMovementX = x - Browser.mouseX;
          Browser.mouseMovementY = y - Browser.mouseY;
          Browser.mouseX = x;
          Browser.mouseY = y;
        }
      },xhrLoad:function (url, onload, onerror) {
        var xhr = new XMLHttpRequest();
        xhr.open('GET', url, true);
        xhr.responseType = 'arraybuffer';
        xhr.onload = function xhr_onload() {
          if (xhr.status == 200 || (xhr.status == 0 && xhr.response)) { // file URLs can return 0
            onload(xhr.response);
          } else {
            onerror();
          }
        };
        xhr.onerror = onerror;
        xhr.send(null);
      },asyncLoad:function (url, onload, onerror, noRunDep) {
        Browser.xhrLoad(url, function(arrayBuffer) {
          assert(arrayBuffer, 'Loading data file "' + url + '" failed (no arrayBuffer).');
          onload(new Uint8Array(arrayBuffer));
          if (!noRunDep) removeRunDependency('al ' + url);
        }, function(event) {
          if (onerror) {
            onerror();
          } else {
            throw 'Loading data file "' + url + '" failed.';
          }
        });
        if (!noRunDep) addRunDependency('al ' + url);
      },resizeListeners:[],updateResizeListeners:function () {
        var canvas = Module['canvas'];
        Browser.resizeListeners.forEach(function(listener) {
          listener(canvas.width, canvas.height);
        });
      },setCanvasSize:function (width, height, noUpdates) {
        var canvas = Module['canvas'];
        Browser.updateCanvasDimensions(canvas, width, height);
        if (!noUpdates) Browser.updateResizeListeners();
      },windowedWidth:0,windowedHeight:0,setFullScreenCanvasSize:function () {
        // check if SDL is available   
        if (typeof SDL != "undefined") {
        	var flags = HEAPU32[((SDL.screen+Runtime.QUANTUM_SIZE*0)>>2)];
        	flags = flags | 0x00800000; // set SDL_FULLSCREEN flag
        	HEAP32[((SDL.screen+Runtime.QUANTUM_SIZE*0)>>2)]=flags
        }
        Browser.updateResizeListeners();
      },setWindowedCanvasSize:function () {
        // check if SDL is available       
        if (typeof SDL != "undefined") {
        	var flags = HEAPU32[((SDL.screen+Runtime.QUANTUM_SIZE*0)>>2)];
        	flags = flags & ~0x00800000; // clear SDL_FULLSCREEN flag
        	HEAP32[((SDL.screen+Runtime.QUANTUM_SIZE*0)>>2)]=flags
        }
        Browser.updateResizeListeners();
      },updateCanvasDimensions:function (canvas, wNative, hNative) {
        if (wNative && hNative) {
          canvas.widthNative = wNative;
          canvas.heightNative = hNative;
        } else {
          wNative = canvas.widthNative;
          hNative = canvas.heightNative;
        }
        var w = wNative;
        var h = hNative;
        if (Module['forcedAspectRatio'] && Module['forcedAspectRatio'] > 0) {
          if (w/h < Module['forcedAspectRatio']) {
            w = Math.round(h * Module['forcedAspectRatio']);
          } else {
            h = Math.round(w / Module['forcedAspectRatio']);
          }
        }
        if (((document['webkitFullScreenElement'] || document['webkitFullscreenElement'] ||
             document['mozFullScreenElement'] || document['mozFullscreenElement'] ||
             document['fullScreenElement'] || document['fullscreenElement'] ||
             document['msFullScreenElement'] || document['msFullscreenElement'] ||
             document['webkitCurrentFullScreenElement']) === canvas.parentNode) && (typeof screen != 'undefined')) {
           var factor = Math.min(screen.width / w, screen.height / h);
           w = Math.round(w * factor);
           h = Math.round(h * factor);
        }
        if (Browser.resizeCanvas) {
          if (canvas.width  != w) canvas.width  = w;
          if (canvas.height != h) canvas.height = h;
          if (typeof canvas.style != 'undefined') {
            canvas.style.removeProperty( "width");
            canvas.style.removeProperty("height");
          }
        } else {
          if (canvas.width  != wNative) canvas.width  = wNative;
          if (canvas.height != hNative) canvas.height = hNative;
          if (typeof canvas.style != 'undefined') {
            if (w != wNative || h != hNative) {
              canvas.style.setProperty( "width", w + "px", "important");
              canvas.style.setProperty("height", h + "px", "important");
            } else {
              canvas.style.removeProperty( "width");
              canvas.style.removeProperty("height");
            }
          }
        }
      }};

  function ___ctype_b_loc() {
      // http://refspecs.freestandards.org/LSB_3.0.0/LSB-Core-generic/LSB-Core-generic/baselib---ctype-b-loc.html
      var me = ___ctype_b_loc;
      if (!me.ret) {
        var values = [
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,8195,8194,8194,8194,8194,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,24577,49156,49156,49156,
          49156,49156,49156,49156,49156,49156,49156,49156,49156,49156,49156,49156,55304,55304,55304,55304,55304,55304,55304,55304,
          55304,55304,49156,49156,49156,49156,49156,49156,49156,54536,54536,54536,54536,54536,54536,50440,50440,50440,50440,50440,
          50440,50440,50440,50440,50440,50440,50440,50440,50440,50440,50440,50440,50440,50440,50440,49156,49156,49156,49156,49156,
          49156,54792,54792,54792,54792,54792,54792,50696,50696,50696,50696,50696,50696,50696,50696,50696,50696,50696,50696,50696,
          50696,50696,50696,50696,50696,50696,50696,49156,49156,49156,49156,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
          0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
        ];
        var i16size = 2;
        var arr = _malloc(values.length * i16size);
        for (var i = 0; i < values.length; i++) {
          HEAP16[(((arr)+(i * i16size))>>1)]=values[i];
        }
        me.ret = allocate([arr + 128 * i16size], 'i16*', ALLOC_NORMAL);
      }
      return me.ret;
    }

  
  function _free() {
  }
  Module["_free"] = _free;function _freelocale(locale) {
      _free(locale);
    }

  function ___cxa_allocate_exception(size) {
      var ptr = _malloc(size + ___cxa_exception_header_size);
      return ptr + ___cxa_exception_header_size;
    }

  var _fmodl=_fmod;

  function _catopen(name, oflag) {
      // nl_catd catopen (const char *name, int oflag)
      return -1;
    }

  function _catgets(catd, set_id, msg_id, s) {
      // char *catgets (nl_catd catd, int set_id, int msg_id, const char *s)
      return s;
    }

  
  
  function _sprintf(s, format, varargs) {
      // int sprintf(char *restrict s, const char *restrict format, ...);
      // http://pubs.opengroup.org/onlinepubs/000095399/functions/printf.html
      return _snprintf(s, undefined, format, varargs);
    }function _asprintf(s, format, varargs) {
      return _sprintf(-s, format, varargs);
    }function _vasprintf(s, format, va_arg) {
      return _asprintf(s, format, HEAP32[((va_arg)>>2)]);
    }

  function _time(ptr) {
      var ret = Math.floor(Date.now()/1000);
      if (ptr) {
        HEAP32[((ptr)>>2)]=ret;
      }
      return ret;
    }

  function _copysign(a, b) {
      return __reallyNegative(a) === __reallyNegative(b) ? a : -a;
    }

  function ___ctype_toupper_loc() {
      // http://refspecs.freestandards.org/LSB_3.1.1/LSB-Core-generic/LSB-Core-generic/libutil---ctype-toupper-loc.html
      var me = ___ctype_toupper_loc;
      if (!me.ret) {
        var values = [
          128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,
          158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,
          188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,
          218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,
          248,249,250,251,252,253,254,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
          33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,
          73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,
          81,82,83,84,85,86,87,88,89,90,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,
          145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,
          175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,
          205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,
          235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255
        ];
        var i32size = 4;
        var arr = _malloc(values.length * i32size);
        for (var i = 0; i < values.length; i++) {
          HEAP32[(((arr)+(i * i32size))>>2)]=values[i];
        }
        me.ret = allocate([arr + 128 * i32size], 'i32*', ALLOC_NORMAL);
      }
      return me.ret;
    }

  function ___cxa_guard_acquire(variable) {
      if (!HEAP8[(variable)]) { // ignore SAFE_HEAP stuff because llvm mixes i64 and i8 here
        HEAP8[(variable)]=1;
        return 1;
      }
      return 0;
    }

  function ___ctype_tolower_loc() {
      // http://refspecs.freestandards.org/LSB_3.1.1/LSB-Core-generic/LSB-Core-generic/libutil---ctype-tolower-loc.html
      var me = ___ctype_tolower_loc;
      if (!me.ret) {
        var values = [
          128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,
          158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,
          188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,
          218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,
          248,249,250,251,252,253,254,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
          33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,97,98,99,100,101,102,103,
          104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,91,92,93,94,95,96,97,98,99,100,101,102,103,
          104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,
          134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,
          164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,
          194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,
          224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,
          254,255
        ];
        var i32size = 4;
        var arr = _malloc(values.length * i32size);
        for (var i = 0; i < values.length; i++) {
          HEAP32[(((arr)+(i * i32size))>>2)]=values[i];
        }
        me.ret = allocate([arr + 128 * i32size], 'i32*', ALLOC_NORMAL);
      }
      return me.ret;
    }

  
  var ___cxa_caught_exceptions=[];function ___cxa_begin_catch(ptr) {
      __ZSt18uncaught_exceptionv.uncaught_exception--;
      ___cxa_caught_exceptions.push(___cxa_last_thrown_exception);
      return ptr;
    }

   
  Module["_bitshift64Shl"] = _bitshift64Shl;


  function __ZNSt9exceptionD2Ev() {}

  var _copysignl=_copysign;

  var __ZTISt9exception=allocate([allocate([1,0,0,0,0,0,0], "i8", ALLOC_STATIC)+8, 0], "i32", ALLOC_STATIC);

  var ___dso_handle=allocate(1, "i32*", ALLOC_STATIC);



FS.staticInit();__ATINIT__.unshift({ func: function() { if (!Module["noFSInit"] && !FS.init.initialized) FS.init() } });__ATMAIN__.push({ func: function() { FS.ignorePermissions = false } });__ATEXIT__.push({ func: function() { FS.quit() } });Module["FS_createFolder"] = FS.createFolder;Module["FS_createPath"] = FS.createPath;Module["FS_createDataFile"] = FS.createDataFile;Module["FS_createPreloadedFile"] = FS.createPreloadedFile;Module["FS_createLazyFile"] = FS.createLazyFile;Module["FS_createLink"] = FS.createLink;Module["FS_createDevice"] = FS.createDevice;
___errno_state = Runtime.staticAlloc(4); HEAP32[((___errno_state)>>2)]=0;
__ATINIT__.unshift({ func: function() { TTY.init() } });__ATEXIT__.push({ func: function() { TTY.shutdown() } });TTY.utf8 = new Runtime.UTF8Processor();
if (ENVIRONMENT_IS_NODE) { var fs = require("fs"); NODEFS.staticInit(); }
__ATINIT__.push({ func: function() { SOCKFS.root = FS.mount(SOCKFS, {}, null); } });
_fgetc.ret = allocate([0], "i8", ALLOC_STATIC);
Module["requestFullScreen"] = function Module_requestFullScreen(lockPointer, resizeCanvas) { Browser.requestFullScreen(lockPointer, resizeCanvas) };
  Module["requestAnimationFrame"] = function Module_requestAnimationFrame(func) { Browser.requestAnimationFrame(func) };
  Module["setCanvasSize"] = function Module_setCanvasSize(width, height, noUpdates) { Browser.setCanvasSize(width, height, noUpdates) };
  Module["pauseMainLoop"] = function Module_pauseMainLoop() { Browser.mainLoop.pause() };
  Module["resumeMainLoop"] = function Module_resumeMainLoop() { Browser.mainLoop.resume() };
  Module["getUserMedia"] = function Module_getUserMedia() { Browser.getUserMedia() }
STACK_BASE = STACKTOP = Runtime.alignMemory(STATICTOP);

staticSealed = true; // seal the static portion of memory

STACK_MAX = STACK_BASE + 5242880;

DYNAMIC_BASE = DYNAMICTOP = Runtime.alignMemory(STACK_MAX);

assert(DYNAMIC_BASE < TOTAL_MEMORY, "TOTAL_MEMORY not big enough for stack");

 var ctlz_i8 = allocate([8,7,6,6,5,5,5,5,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], "i8", ALLOC_DYNAMIC);
 var cttz_i8 = allocate([8,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0], "i8", ALLOC_DYNAMIC);

var Math_min = Math.min;
function invoke_iiii(index,a1,a2,a3) {
  try {
    return Module["dynCall_iiii"](index,a1,a2,a3);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_viiiiiii(index,a1,a2,a3,a4,a5,a6,a7) {
  try {
    Module["dynCall_viiiiiii"](index,a1,a2,a3,a4,a5,a6,a7);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_viiiii(index,a1,a2,a3,a4,a5) {
  try {
    Module["dynCall_viiiii"](index,a1,a2,a3,a4,a5);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_vi(index,a1) {
  try {
    Module["dynCall_vi"](index,a1);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_vii(index,a1,a2) {
  try {
    Module["dynCall_vii"](index,a1,a2);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_viiiiiiiii(index,a1,a2,a3,a4,a5,a6,a7,a8,a9) {
  try {
    Module["dynCall_viiiiiiiii"](index,a1,a2,a3,a4,a5,a6,a7,a8,a9);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_ii(index,a1) {
  try {
    return Module["dynCall_ii"](index,a1);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_viiiiiid(index,a1,a2,a3,a4,a5,a6,a7) {
  try {
    Module["dynCall_viiiiiid"](index,a1,a2,a3,a4,a5,a6,a7);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_viii(index,a1,a2,a3) {
  try {
    Module["dynCall_viii"](index,a1,a2,a3);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_viiiiid(index,a1,a2,a3,a4,a5,a6) {
  try {
    Module["dynCall_viiiiid"](index,a1,a2,a3,a4,a5,a6);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_v(index) {
  try {
    Module["dynCall_v"](index);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_iiiiiiiii(index,a1,a2,a3,a4,a5,a6,a7,a8) {
  try {
    return Module["dynCall_iiiiiiiii"](index,a1,a2,a3,a4,a5,a6,a7,a8);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_iiiii(index,a1,a2,a3,a4) {
  try {
    return Module["dynCall_iiiii"](index,a1,a2,a3,a4);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_viiiiiiii(index,a1,a2,a3,a4,a5,a6,a7,a8) {
  try {
    Module["dynCall_viiiiiiii"](index,a1,a2,a3,a4,a5,a6,a7,a8);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_viiiiii(index,a1,a2,a3,a4,a5,a6) {
  try {
    Module["dynCall_viiiiii"](index,a1,a2,a3,a4,a5,a6);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_iii(index,a1,a2) {
  try {
    return Module["dynCall_iii"](index,a1,a2);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_iiiiii(index,a1,a2,a3,a4,a5) {
  try {
    return Module["dynCall_iiiiii"](index,a1,a2,a3,a4,a5);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function invoke_viiii(index,a1,a2,a3,a4) {
  try {
    Module["dynCall_viiii"](index,a1,a2,a3,a4);
  } catch(e) {
    if (typeof e !== 'number' && e !== 'longjmp') throw e;
    asm["setThrew"](1, 0);
  }
}

function asmPrintInt(x, y) {
  Module.print('int ' + x + ',' + y);// + ' ' + new Error().stack);
}
function asmPrintFloat(x, y) {
  Module.print('float ' + x + ',' + y);// + ' ' + new Error().stack);
}
// EMSCRIPTEN_START_ASM
var asm=(function(global,env,buffer){"use asm";var a=new global.Int8Array(buffer);var b=new global.Int16Array(buffer);var c=new global.Int32Array(buffer);var d=new global.Uint8Array(buffer);var e=new global.Uint16Array(buffer);var f=new global.Uint32Array(buffer);var g=new global.Float32Array(buffer);var h=new global.Float64Array(buffer);var i=env.STACKTOP|0;var j=env.STACK_MAX|0;var k=env.tempDoublePtr|0;var l=env.ABORT|0;var m=env.cttz_i8|0;var n=env.ctlz_i8|0;var o=env.__ZTISt9exception|0;var p=env.___dso_handle|0;var q=env._stderr|0;var r=env._stdin|0;var s=env._stdout|0;var t=0;var u=0;var v=0;var w=0;var x=+env.NaN,y=+env.Infinity;var z=0,A=0,B=0,C=0,D=0.0,E=0,F=0,G=0,H=0.0;var I=0;var J=0;var K=0;var L=0;var M=0;var N=0;var O=0;var P=0;var Q=0;var R=0;var S=global.Math.floor;var T=global.Math.abs;var U=global.Math.sqrt;var V=global.Math.pow;var W=global.Math.cos;var X=global.Math.sin;var Y=global.Math.tan;var Z=global.Math.acos;var _=global.Math.asin;var $=global.Math.atan;var aa=global.Math.atan2;var ba=global.Math.exp;var ca=global.Math.log;var da=global.Math.ceil;var ea=global.Math.imul;var fa=env.abort;var ga=env.assert;var ha=env.asmPrintInt;var ia=env.asmPrintFloat;var ja=env.min;var ka=env.invoke_iiii;var la=env.invoke_viiiiiii;var ma=env.invoke_viiiii;var na=env.invoke_vi;var oa=env.invoke_vii;var pa=env.invoke_viiiiiiiii;var qa=env.invoke_ii;var ra=env.invoke_viiiiiid;var sa=env.invoke_viii;var ta=env.invoke_viiiiid;var ua=env.invoke_v;var va=env.invoke_iiiiiiiii;var wa=env.invoke_iiiii;var xa=env.invoke_viiiiiiii;var ya=env.invoke_viiiiii;var za=env.invoke_iii;var Aa=env.invoke_iiiiii;var Ba=env.invoke_viiii;var Ca=env._fabs;var Da=env._pthread_cond_wait;var Ea=env._freelocale;var Fa=env.__formatString;var Ga=env._asprintf;var Ha=env.__reallyNegative;var Ia=env._send;var Ja=env._strtoll_l;var Ka=env._vsscanf;var La=env.___ctype_b_loc;var Ma=env.__ZSt9terminatev;var Na=env._fmod;var Oa=env.___cxa_guard_acquire;var Pa=env.___setErrNo;var Qa=env.___cxa_is_number_type;var Ra=env._atexit;var Sa=env._copysign;var Ta=env._emscripten_asm_const_int;var Ua=env.___assert_fail;var Va=env._ungetc;var Wa=env.___cxa_allocate_exception;var Xa=env.__ZSt18uncaught_exceptionv;var Ya=env._isxdigit_l;var Za=env.___ctype_toupper_loc;var _a=env._fflush;var $a=env.___cxa_guard_release;var ab=env.__addDays;var bb=env.___errno_location;var cb=env._strtoll;var db=env._strerror_r;var eb=env._emscripten_asm_const_double;var fb=env._strftime_l;var gb=env._sscanf;var hb=env._fabsf;var ib=env._sbrk;var jb=env._uselocale;var kb=env._catgets;var lb=env._newlocale;var mb=env._snprintf;var nb=env.___cxa_begin_catch;var ob=env._emscripten_memcpy_big;var pb=env._fileno;var qb=env.___resumeException;var rb=env.___cxa_find_matching_catch;var sb=env._sysconf;var tb=env._strtoull;var ub=env._isdigit_l;var vb=env._strftime;var wb=env.__arraySum;var xb=env.___cxa_throw;var yb=env.___ctype_tolower_loc;var zb=env._isspace;var Ab=env._pthread_mutex_unlock;var Bb=env._fread;var Cb=env._pthread_cond_broadcast;var Db=env._isxdigit;var Eb=env._sprintf;var Fb=env._sqrtf;var Gb=env._vasprintf;var Hb=env._write;var Ib=env.__isLeapYear;var Jb=env.__scanString;var Kb=env._recv;var Lb=env._vsnprintf;var Mb=env.__ZNSt9exceptionD2Ev;var Nb=env._fgetc;var Ob=env._strtoull_l;var Pb=env._mkport;var Qb=env.___cxa_does_inherit;var Rb=env.__exit;var Sb=env._read;var Tb=env.__parseInt64;var Ub=env._abort;var Vb=env._catclose;var Wb=env._fwrite;var Xb=env._time;var Yb=env._pthread_mutex_lock;var Zb=env._strerror;var _b=env._pread;var $b=env._pwrite;var ac=env._sqrt;var bc=env._catopen;var cc=env._exit;var dc=env._isdigit;var ec=env.__getFloat;var fc=0.0;
// EMSCRIPTEN_START_FUNCS
function yc(a){a=a|0;var b=0;b=i;i=i+a|0;i=i+7&-8;return b|0}function zc(){return i|0}function Ac(a){a=a|0;i=a}function Bc(a,b){a=a|0;b=b|0;if((t|0)==0){t=a;u=b}}function Cc(b){b=b|0;a[k]=a[b];a[k+1|0]=a[b+1|0];a[k+2|0]=a[b+2|0];a[k+3|0]=a[b+3|0]}function Dc(b){b=b|0;a[k]=a[b];a[k+1|0]=a[b+1|0];a[k+2|0]=a[b+2|0];a[k+3|0]=a[b+3|0];a[k+4|0]=a[b+4|0];a[k+5|0]=a[b+5|0];a[k+6|0]=a[b+6|0];a[k+7|0]=a[b+7|0]}function Ec(a){a=a|0;I=a}function Fc(a){a=a|0;J=a}function Gc(a){a=a|0;K=a}function Hc(a){a=a|0;L=a}function Ic(a){a=a|0;M=a}function Jc(a){a=a|0;N=a}function Kc(a){a=a|0;O=a}function Lc(a){a=a|0;P=a}function Mc(a){a=a|0;Q=a}function Nc(a){a=a|0;R=a}function Oc(){var b=0,d=0,e=0,f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0.0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0;b=i;i=i+304|0;d=b+280|0;e=b+240|0;f=b+236|0;h=b+232|0;j=b+140|0;k=b+104|0;l=b+88|0;m=b+52|0;n=b;o=b+96|0;p=Ta(8,0)|0;q=Ta(32,0)|0;a:do{if((q|0)>0){r=0;s=0;t=0;u=0;while(1){v=Ta(56,t|0)|0;w=Ta(80,t|0)|0;x=+eb(104,t|0);if((u|0)==(r|0)){y=u-s|0;z=(y|0)/12|0;A=z+1|0;if(A>>>0>357913941){break}if(z>>>0<178956970){B=z<<1;C=B>>>0<A>>>0?A:B;if((C|0)==0){D=0;E=0}else{F=C;G=10}}else{F=357913941;G=10}if((G|0)==10){G=0;D=F;E=Sm(F*12|0)|0}C=E+(z*12|0)|0;B=E+(D*12|0)|0;if((C|0)!=0){c[C>>2]=v;c[E+(z*12|0)+4>>2]=w;g[E+(z*12|0)+8>>2]=x}C=E+(A*12|0)|0;A=E+((((y|0)/-12|0)+z|0)*12|0)|0;nn(A|0,s|0,y|0)|0;if((s|0)==0){H=B;I=A;J=C}else{Um(s);H=B;I=A;J=C}}else{if((u|0)==0){K=0}else{c[u>>2]=v;c[u+4>>2]=w;g[u+8>>2]=x;K=u}H=r;I=s;J=K+12|0}w=t+1|0;if((w|0)<(q|0)){r=H;s=I;t=w;u=J}else{L=I;M=J;break a}}Sj(0)}else{L=0;M=0}}while(0);a[e]=0;J=e+4|0;c[J+0>>2]=0;c[J+4>>2]=0;c[J+8>>2]=0;c[J+12>>2]=0;c[J+16>>2]=0;c[J+20>>2]=0;c[J+24>>2]=0;c[J+28>>2]=0;c[e+8>>2]=p;c[e+28>>2]=0;I=e+12|0;H=Lm((p<<2)+4|0)|0;c[I>>2]=H;if((H|0)==0){Wc()}c[J>>2]=p;J=e+16|0;pn(H|0,0,(p<<2)+4|0)|0;c[f>>2]=L;c[h>>2]=M;rd(f,h,e,0);gd(j,e);if((a[j+5|0]|0)==0){Ua(2064,2128,324,2184)}a[l]=0;c[l+4>>2]=j+8;a[k]=0;h=k+4|0;c[h+0>>2]=0;c[h+4>>2]=0;c[h+8>>2]=0;c[h+12>>2]=0;c[h+16>>2]=0;c[h+20>>2]=0;c[h+24>>2]=0;c[h+28>>2]=0;ad(k,l)|0;a[o]=0;c[o+4>>2]=k;_c(n,k,o);a[m]=0;o=m+4|0;l=n+48|0;c[o+0>>2]=0;c[o+4>>2]=0;c[o+8>>2]=0;c[o+12>>2]=0;c[o+16>>2]=0;c[o+20>>2]=0;c[o+24>>2]=0;c[o+28>>2]=0;o=c[n+4>>2]|0;f=n+8|0;if((a[l]|0)==0){Sc(o,f,m,n+44|0)}else{Rc(o,f,m)}Mm(c[n+20>>2]|0);Mm(c[n+24>>2]|0);f=c[n+28>>2]|0;if((f|0)!=0){Vm(f)}f=c[n+32>>2]|0;if((f|0)!=0){Vm(f)}if((c[m+8>>2]|0)!=(c[e+8>>2]|0)){Ua(336,400,131,440)}if((c[m+4>>2]|0)!=(c[e+4>>2]|0)){Ua(336,400,131,440)}a[d]=0;a[d+4|0]=0;f=d+8|0;c[f>>2]=m;c[f+4>>2]=e;x=+U(+(+Qc(d)));d=c[h>>2]|0;if((d|0)>0){f=k+20|0;n=k+24|0;o=k+12|0;l=k+16|0;M=d;d=0;p=0;while(1){H=c[f>>2]|0;q=c[n>>2]|0;K=c[o>>2]|0;E=c[K+(d<<2)>>2]|0;D=c[l>>2]|0;if((D|0)==0){N=c[K+(d+1<<2)>>2]|0}else{N=(c[D+(d<<2)>>2]|0)+E|0}if((E|0)<(N|0)){D=E;E=p;do{Ta(128,E|0,c[q+(D<<2)>>2]|0,d|0,+(+g[H+(D<<2)>>2]))|0;E=E+1|0;D=D+1|0}while((D|0)<(N|0));O=c[h>>2]|0;P=E}else{O=M;P=p}D=d+1|0;if((D|0)<(O|0)){M=O;d=D;p=P}else{Q=P;break}}}else{Q=0}Ta(192,+x)|0;Mm(c[m+12>>2]|0);Mm(c[m+16>>2]|0);P=c[m+20>>2]|0;if((P|0)!=0){Vm(P)}P=c[m+24>>2]|0;if((P|0)!=0){Vm(P)}Mm(c[k+12>>2]|0);Mm(c[k+16>>2]|0);P=c[k+20>>2]|0;if((P|0)!=0){Vm(P)}P=c[k+24>>2]|0;if((P|0)!=0){Vm(P)}ud(j);Mm(c[I>>2]|0);Mm(c[J>>2]|0);J=c[e+20>>2]|0;if((J|0)!=0){Vm(J)}J=c[e+24>>2]|0;if((J|0)!=0){Vm(J)}if((L|0)==0){i=b;return Q|0}Um(L);i=b;return Q|0}function Pc(a){a=a|0;nb(a|0)|0;Ma()}function Qc(a){a=a|0;var b=0,d=0,e=0,f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0.0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0.0,B=0,C=0,D=0,E=0.0,F=0.0,G=0.0,H=0.0,I=0,J=0.0,K=0,L=0,M=0.0,N=0,O=0,P=0,Q=0.0,R=0,S=0;b=i;d=c[a+12>>2]|0;if((c[d+8>>2]|0)<=0){Ua(208,280,19,328)}e=c[d+4>>2]|0;if((e|0)<=0){Ua(208,280,19,328)}f=c[a+8>>2]|0;a=c[f+20>>2]|0;h=c[f+24>>2]|0;j=c[f+12>>2]|0;k=c[f+16>>2]|0;f=(k|0)==0;l=c[d+20>>2]|0;m=c[d+24>>2]|0;n=c[d+12>>2]|0;o=c[d+16>>2]|0;d=(o|0)==0;p=0;q=0.0;while(1){r=c[j+(p<<2)>>2]|0;if(f){s=c[j+(p+1<<2)>>2]|0}else{s=(c[k+(p<<2)>>2]|0)+r|0}t=c[n+(p<<2)>>2]|0;if(d){u=c[n+(p+1<<2)>>2]|0}else{u=(c[o+(p<<2)>>2]|0)+t|0}v=(r|0)<(s|0);a:do{if(v){w=c[h+(r<<2)>>2]|0;do{if((t|0)<(u|0)){x=c[m+(t<<2)>>2]|0;if((w|0)!=(x|0)){if((w|0)<(x|0)){break}else{y=17;break a}}else{z=w;A=+g[a+(r<<2)>>2]- +g[l+(t<<2)>>2];B=t+1|0;C=r+1|0;y=22;break a}}}while(0);z=w;A=+g[a+(r<<2)>>2];B=t;C=r+1|0;y=22}else{y=17}}while(0);do{if((y|0)==17){y=0;if((t|0)<(u|0)){if(v){x=c[m+(t<<2)>>2]|0;if((c[h+(r<<2)>>2]|0)>(x|0)){D=x}else{E=q;break}}else{D=c[m+(t<<2)>>2]|0}z=D;A=0.0- +g[l+(t<<2)>>2];B=t+1|0;C=r;y=22}else{E=q}}}while(0);b:do{if((y|0)==22){y=0;if((z|0)>-1){F=A;r=B;t=C;G=q;while(1){v=(r|0)<(u|0);x=m+(r<<2)|0;H=F;I=t;J=G;while(1){K=(I|0)<(s|0);L=h+(I<<2)|0;M=H*H+J;if(!K){y=32;break}N=c[L>>2]|0;if(v){O=c[x>>2]|0;if((N|0)==(O|0)){y=28;break}if((N|0)>=(O|0)){y=32;break}}if((N|0)>-1){H=+g[a+(I<<2)>>2];I=I+1|0;J=M}else{E=M;break b}}if((y|0)==28){y=0;P=N;Q=+g[a+(I<<2)>>2]- +g[l+(r<<2)>>2];R=I+1|0}else if((y|0)==32){y=0;if(!v){E=M;break b}if(K){O=c[x>>2]|0;if((c[L>>2]|0)>(O|0)){S=O}else{E=M;break b}}else{S=c[x>>2]|0}P=S;Q=0.0- +g[l+(r<<2)>>2];R=I}if((P|0)>-1){F=Q;r=r+1|0;t=R;G=M}else{E=M;break}}}else{E=q}}}while(0);p=p+1|0;if((p|0)>=(e|0)){break}else{q=E}}i=b;return+E}function Rc(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0,X=0,Y=0,Z=0,_=0,$=0,aa=0.0,ba=0,ca=0,da=0,ea=0,fa=0.0,ga=0,ha=0,ia=0,ja=0,ka=0,la=0,ma=0,na=0,oa=0,pa=0,qa=0,ra=0,sa=0,ta=0;f=i;i=i+96|0;h=f+48|0;j=f;k=f+84|0;l=f+12|0;m=b+8|0;n=c[m>>2]|0;o=d+4|0;p=c[o>>2]|0;a[l]=0;q=l+4|0;c[q+0>>2]=0;c[q+4>>2]=0;c[q+8>>2]=0;c[q+12>>2]=0;c[q+16>>2]=0;c[q+20>>2]=0;c[q+24>>2]=0;c[q+28>>2]=0;c[l+8>>2]=n;n=l+28|0;c[n>>2]=0;r=l+12|0;s=(p<<2)+4|0;t=Lm(s)|0;c[r>>2]=t;if((t|0)==0){Wc()}c[q>>2]=p;p=l+16|0;pn(t|0,0,s|0)|0;s=c[m>>2]|0;m=c[o>>2]|0;u=b+4|0;if((c[u>>2]|0)!=(c[d+8>>2]|0)){Ua(456,1872,26,1936)}a[k]=0;Zc(j,s,k);if(s>>>0>1073741823){Wc()}k=s<<2;if((Qm(h,16,k)|0)==0){v=c[h>>2]|0}else{c[h>>2]=0;v=0}w=(k|0)==0;if(!((v|0)!=0|w)){Wc()}if(!((s|0)>-1)){Ua(1696,1712,225,1752)}if((Qm(h,16,k)|0)==0){x=c[h>>2]|0}else{c[h>>2]=0;x=0}if(!((x|0)!=0|w)){Wc()}w=b+16|0;k=c[w>>2]|0;if((k|0)!=0){y=c[u>>2]|0;if(!((y|0)>-1)){Ua(1536,1552,138,1592)}if((y|0)!=0){u=c[k>>2]|0;if((y|0)>1){z=u;A=1;while(1){B=(c[k+(A<<2)>>2]|0)+z|0;C=A+1|0;if((C|0)==(y|0)){D=B;break}else{z=B;A=C}}}else{D=u}}else{D=0}}else{D=c[b+28>>2]|0}u=d+16|0;A=c[u>>2]|0;if((A|0)!=0){z=c[o>>2]|0;if(!((z|0)>-1)){Ua(1536,1552,138,1592)}if((z|0)!=0){o=c[A>>2]|0;if((z|0)>1){y=o;k=1;while(1){C=(c[A+(k<<2)>>2]|0)+y|0;B=k+1|0;if((B|0)==(z|0)){E=C;break}else{y=C;k=B}}}else{E=o}}else{E=0}}else{E=c[d+28>>2]|0}o=E+D|0;c[n>>2]=0;D=c[q>>2]<<2;pn(t|0,0,D+4|0)|0;E=c[p>>2]|0;if((E|0)!=0){pn(E|0,0,D|0)|0;Ua(1344,704,258,1424)}D=l+32|0;if((o|0)==0){F=0;G=0}else{E=o>>>0>1073741823?-1:o<<2;k=Tm(E)|0;y=Tm(E)|0;E=l+20|0;z=c[E>>2]|0;A=l+24|0;if((z|0)==0){H=c[A>>2]|0}else{Vm(z);H=c[A>>2]|0}if((H|0)!=0){Vm(H)}c[E>>2]=k;c[A>>2]=y;c[D>>2]=o;F=y;G=k}a:do{if((m|0)>0){k=d+20|0;y=d+24|0;o=d+12|0;A=l+20|0;E=l+24|0;H=c[j>>2]|0;z=b+20|0;B=b+24|0;C=b+12|0;I=F;J=F;K=G;L=F;M=G;N=0;O=0;b:while(1){if((c[t+(O<<2)>>2]|0)!=(N|0)){P=46;break}Q=O+1|0;R=t+(Q<<2)|0;if((c[R>>2]|0)!=0){P=48;break}c[R>>2]=N;S=c[k>>2]|0;T=c[y>>2]|0;U=c[o>>2]|0;V=c[U+(O<<2)>>2]|0;W=c[u>>2]|0;if((W|0)==0){X=c[U+(Q<<2)>>2]|0}else{X=(c[W+(O<<2)>>2]|0)+V|0}if((V|0)<(X|0)){W=c[z>>2]|0;U=c[B>>2]|0;Y=c[C>>2]|0;Z=c[w>>2]|0;_=(Z|0)==0;$=V;V=0;while(1){aa=+g[S+($<<2)>>2];ba=c[T+($<<2)>>2]|0;ca=c[Y+(ba<<2)>>2]|0;if(_){da=c[Y+(ba+1<<2)>>2]|0}else{da=(c[Z+(ba<<2)>>2]|0)+ca|0}if((ca|0)<(da|0)){ba=ca;ca=V;while(1){ea=c[U+(ba<<2)>>2]|0;fa=+g[W+(ba<<2)>>2];ga=H+(ea>>>5<<2)|0;ha=1<<(ea&31);ia=c[ga>>2]|0;if((ia&ha|0)==0){c[ga>>2]=ia|ha;if(!((ea|0)>-1&(s|0)>(ea|0))){P=60;break b}g[v+(ea<<2)>>2]=aa*fa;if(!((ca|0)>-1&(s|0)>(ca|0))){P=62;break b}c[x+(ca<<2)>>2]=ea;ja=ca+1|0}else{if(!((ea|0)>-1&(s|0)>(ea|0))){P=65;break b}ha=v+(ea<<2)|0;g[ha>>2]=aa*fa+ +g[ha>>2];ja=ca}ha=ba+1|0;if((ha|0)<(da|0)){ba=ha;ca=ja}else{ka=ja;break}}}else{ka=V}$=$+1|0;if(($|0)>=(X|0)){break}else{V=ka}}if((ka|0)>0){V=I;$=M;W=L;U=J;Z=K;Y=0;while(1){if((s|0)<=(Y|0)){P=71;break b}_=c[x+(Y<<2)>>2]|0;if(!((_|0)>-1&(s|0)>(_|0))){P=73;break b}aa=+g[v+(_<<2)>>2];T=c[R>>2]|0;c[R>>2]=T+1;S=c[n>>2]|0;ca=S+1|0;if((c[D>>2]|0)>>>0<ca>>>0){ba=(~~+(ca>>>0)>>>0)+ca|0;ha=ba>>>0>1073741823?-1:ba<<2;ea=Tm(ha)|0;ia=Tm(ha)|0;ha=(S>>>0<ba>>>0?S:ba)<<2;nn(ea|0,$|0,ha|0)|0;nn(ia|0,W|0,ha|0)|0;if(($|0)==0){la=W}else{Vm($);la=V}if((la|0)!=0){Vm(la)}c[A>>2]=ea;c[E>>2]=ia;c[D>>2]=ba;ma=ea;na=ia}else{ma=Z;na=U}c[n>>2]=ca;g[ma+(S<<2)>>2]=0.0;c[na+(S<<2)>>2]=_;g[ma+(T<<2)>>2]=aa;T=H+(_>>>5<<2)|0;c[T>>2]=c[T>>2]&~(1<<(_&31));_=Y+1|0;if((_|0)<(ka|0)){V=na;$=ma;W=na;U=na;Z=ma;Y=_}else{oa=na;pa=na;qa=ma;ra=na;sa=ma;break}}}else{oa=I;pa=J;qa=K;ra=L;sa=M}}else{oa=I;pa=J;qa=K;ra=L;sa=M}if((Q|0)>=(m|0)){break a}I=oa;J=pa;K=qa;L=ra;M=sa;N=c[n>>2]|0;O=Q}if((P|0)==46){Ua(920,704,405,1032)}else if((P|0)==48){Ua(1048,704,406,1032)}else if((P|0)==60){Ua(1600,1632,378,1680)}else if((P|0)==62){Ua(1600,1632,378,1680)}else if((P|0)==65){Ua(1600,1632,378,1680)}else if((P|0)==71){Ua(1600,1632,378,1680)}else if((P|0)==73){Ua(1600,1632,378,1680)}}}while(0);P=c[n>>2]|0;n=c[q>>2]|0;if((n|0)>-1){q=n;while(1){sa=q+ -1|0;if((c[t+(q<<2)>>2]|0)!=0){ta=q;break}if((q|0)>0){q=sa}else{ta=sa;break}}if((ta|0)<(n|0)){q=ta;do{q=q+1|0;c[t+(q<<2)>>2]=P}while((q|0)<(n|0))}}Mm(x);Mm(v);v=c[j>>2]|0;if((v|0)!=0){Um(v)}a[h]=0;v=h+4|0;c[v+0>>2]=0;c[v+4>>2]=0;c[v+8>>2]=0;c[v+12>>2]=0;c[v+16>>2]=0;c[v+20>>2]=0;c[v+24>>2]=0;c[v+28>>2]=0;Yc(h,l)|0;Xc(e,h)|0;Mm(c[h+12>>2]|0);Mm(c[h+16>>2]|0);e=c[h+20>>2]|0;if((e|0)!=0){Vm(e)}e=c[h+24>>2]|0;if((e|0)!=0){Vm(e)}Mm(c[r>>2]|0);Mm(c[p>>2]|0);p=c[l+20>>2]|0;if((p|0)!=0){Vm(p)}p=c[l+24>>2]|0;if((p|0)==0){i=f;return}Vm(p);i=f;return}function Sc(a,b,d,e){a=a|0;b=b|0;d=d|0;e=e|0;var f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,U=0,V=0,W=0,X=0,Y=0.0,Z=0,_=0,$=0,aa=0,ba=0,ca=0.0,da=0,fa=0.0,ga=0,ha=0,ia=0,ja=0,ka=0.0,la=0,ma=0,na=0,oa=0,pa=0,qa=0,ra=0,sa=0,ta=0,ua=0,va=0,wa=0,xa=0,ya=0,za=0,Aa=0.0,Ba=0.0,Ca=0,Da=0,Ea=0,Fa=0.0,Ga=0,Ha=0,Ia=0,Ja=0,Ka=0,La=0;f=i;i=i+48|0;h=f+4|0;j=f;k=d+8|0;l=d+4|0;m=c[l>>2]|0;c[j>>2]=0;n=(m<<2)+4|0;o=Lm(n)|0;if((o|0)==0){Wc()}pn(o|0,0,n|0)|0;n=a+8|0;p=c[n>>2]|0;q=b+4|0;r=c[q>>2]|0;s=c[a+4>>2]|0;if((s|0)!=(c[b+8>>2]|0)){Ua(456,496,31,560)}c[h>>2]=0;g[h+4>>2]=0.0;t=h+8|0;c[t>>2]=0;c[h+20>>2]=0;u=h+24|0;c[u>>2]=0;v=h+28|0;c[v>>2]=-1;if((p|0)>0){if((p|0)<1e3){w=p*12|0;c[u>>2]=(w>>>0)/12|0;x=Tm(w)|0}else{c[u>>2]=(p<<2>>>0)/12|0;x=Tm(p>>>0>1073741823?-1:p<<2)|0}c[h>>2]=x;c[t>>2]=p;c[h+12>>2]=0;c[h+16>>2]=p;y=x}else{y=0}c[t>>2]=p;t=a+16|0;x=c[t>>2]|0;if((x|0)!=0){if(!((s|0)>-1)){Ua(1536,1552,138,1592)}if((s|0)!=0){u=c[x>>2]|0;if((s|0)>1){w=u;z=1;while(1){A=(c[x+(z<<2)>>2]|0)+w|0;B=z+1|0;if((B|0)==(s|0)){C=A;break}else{w=A;z=B}}}else{C=u}}else{C=0}}else{C=c[a+28>>2]|0}u=b+16|0;z=c[u>>2]|0;if((z|0)!=0){if(!((r|0)>-1)){Ua(1536,1552,138,1592)}if((r|0)!=0){w=c[z>>2]|0;if((r|0)>1){s=w;x=1;while(1){B=(c[z+(x<<2)>>2]|0)+s|0;A=x+1|0;if((A|0)==(r|0)){D=B;break}else{s=B;x=A}}}else{D=w}}else{D=0}}else{D=c[b+28>>2]|0}w=D+C|0;c[j>>2]=0;if((m|0)!=(r|0)|(m|0)==0){Mm(o);C=(r<<2)+4|0;D=Lm(C)|0;if((D|0)==0){Wc()}else{E=C;F=r;G=D}}else{E=(r<<2)+4|0;F=m;G=o}pn(G|0,0,E|0)|0;if((w|0)==0){H=0;I=0;J=0}else{E=w>>>0>1073741823?-1:w<<2;o=Tm(E)|0;H=w;I=Tm(E)|0;J=o}do{if((r|0)>0){o=+(w|0)/+(ea(c[n>>2]|0,c[q>>2]|0)|0)>.1;E=h+12|0;m=h+16|0;D=b+20|0;C=b+24|0;x=b+12|0;s=h+32|0;z=h+36|0;A=a+20|0;B=a+24|0;K=a+12|0;L=h+40|0;M=H;N=I;O=J;P=y;Q=0;a:while(1){if(o){c[v>>2]=0;R=c[E>>2]|0;S=c[m>>2]|0;if((R|0)<(S|0)){U=R+1|0;pn(P+(R<<2)|0,0,((S|0)>(U|0)?S:U)-R<<2|0)|0}}else{c[v>>2]=1;c[L>>2]=0;c[s>>2]=-1}R=c[D>>2]|0;U=c[C>>2]|0;S=c[x>>2]|0;V=c[S+(Q<<2)>>2]|0;W=c[u>>2]|0;if((W|0)==0){X=c[S+(Q+1<<2)>>2]|0}else{X=(c[W+(Q<<2)>>2]|0)+V|0}if((V|0)<(X|0)){W=V;do{c[z>>2]=c[s>>2];Y=+g[R+(W<<2)>>2];V=c[U+(W<<2)>>2]|0;S=c[A>>2]|0;Z=c[B>>2]|0;_=c[K>>2]|0;$=c[_+(V<<2)>>2]|0;aa=c[t>>2]|0;if((aa|0)==0){ba=c[_+(V+1<<2)>>2]|0}else{ba=(c[aa+(V<<2)>>2]|0)+$|0}if(($|0)<(ba|0)){V=$;do{ca=+g[S+(V<<2)>>2];$=Uc(h,c[Z+(V<<2)>>2]|0)|0;g[$>>2]=Y*ca+ +g[$>>2];V=V+1|0}while((V|0)<(ba|0))}W=W+1|0}while((W|0)<(X|0))}W=G+(Q<<2)|0;U=c[W>>2]|0;if((U|0)!=(c[j>>2]|0)){da=48;break}R=Q+1|0;V=G+(R<<2)|0;if((c[V>>2]|0)!=0){da=50;break}c[V>>2]=U;Y=+g[e>>2];Z=(c[v>>2]|0)==0;b:do{if(!Z){S=c[h>>2]|0;$=c[s>>2]|0;if(($|0)>-1){aa=$;while(1){fa=+g[S+(aa*12|0)+8>>2];if(!(+T(+fa)<=Y)){break}$=c[S+(aa*12|0)>>2]|0;if(($|0)>-1){aa=$}else{ga=S;ha=M;ia=N;ja=O;break b}}ka=fa;la=c[S+(aa*12|0)+4>>2]|0;ma=aa;na=S;da=59}else{ga=S;ha=M;ia=N;ja=O}}else{$=c[m>>2]|0;_=c[h>>2]|0;oa=(c[E>>2]|0)+ -1|0;while(1){pa=oa+1|0;if((pa|0)>=($|0)){ga=_;ha=M;ia=N;ja=O;break b}ca=+g[_+(pa<<2)>>2];if(+T(+ca)<Y){oa=pa}else{ka=ca;la=pa;ma=0;na=_;da=59;break}}}}while(0);c:do{if((da|0)==59){da=0;if((la|0)>-1){ca=ka;_=la;oa=ma;$=U;S=U;aa=N;pa=M;qa=O;ra=na;while(1){if(($|0)!=(S|0)){da=61;break a}if((S|0)!=(c[W>>2]|0)?(c[aa+(S+ -1<<2)>>2]|0)>=(_|0):0){da=64;break a}sa=S+1|0;c[V>>2]=sa;if(pa>>>0<sa>>>0){ta=(~~+(sa>>>0)>>>0)+sa|0;ua=ta>>>0>1073741823?-1:ta<<2;va=Tm(ua)|0;wa=Tm(ua)|0;ua=(S>>>0<ta>>>0?S:ta)<<2;nn(va|0,qa|0,ua|0)|0;nn(wa|0,aa|0,ua|0)|0;if((qa|0)!=0){Vm(qa)}if((aa|0)==0){xa=va;ya=wa;za=ta}else{Vm(aa);xa=va;ya=wa;za=ta}}else{xa=qa;ya=aa;za=pa}c[j>>2]=sa;c[ya+(S<<2)>>2]=_;g[xa+(S<<2)>>2]=ca;if(Z){sa=c[m>>2]|0;ta=_;while(1){wa=ta+1|0;if((wa|0)>=(sa|0)){ga=na;ha=za;ia=ya;ja=xa;break c}Aa=+g[na+(wa<<2)>>2];if(+T(+Aa)<Y){ta=wa}else{Ba=Aa;Ca=wa;Da=oa;Ea=na;break}}}else{ta=oa;do{ta=c[ra+(ta*12|0)>>2]|0;if(!((ta|0)>-1)){ga=na;ha=za;ia=ya;ja=xa;break c}Fa=+g[ra+(ta*12|0)+8>>2]}while(+T(+Fa)<Y);if((ta|0)<0){ga=na;ha=za;ia=ya;ja=xa;break c}Ba=Fa;Ca=c[ra+(ta*12|0)+4>>2]|0;Da=ta;Ea=ra}if(!((Ca|0)>-1)){ga=na;ha=za;ia=ya;ja=xa;break c}ca=Ba;_=Ca;oa=Da;$=c[V>>2]|0;S=c[j>>2]|0;aa=ya;pa=za;qa=xa;ra=Ea}}else{ga=na;ha=M;ia=N;ja=O}}}while(0);if((R|0)<(r|0)){M=ha;N=ia;O=ja;P=ga;Q=R}else{da=81;break}}if((da|0)==48){Ua(920,704,405,1032)}else if((da|0)==50){Ua(1048,704,406,1032)}else if((da|0)==61){Ua(600,704,383,752)}else if((da|0)==64){Ua(776,704,384,752)}else if((da|0)==81){Ga=c[j>>2]|0;Ha=ga;Ia=ja;Ja=ia;Ka=ha;break}}else{Ga=0;Ha=y;Ia=J;Ja=I;Ka=H}}while(0);if((F|0)>-1){H=F;while(1){I=H+ -1|0;if((c[G+(H<<2)>>2]|0)!=0){La=H;break}if((H|0)>0){H=I}else{La=I;break}}if((La|0)<(F|0)){H=La;do{H=H+1|0;c[G+(H<<2)>>2]=Ga}while((H|0)!=(F|0))}}if((Ha|0)!=0){Vm(Ha)}Ha=d+12|0;H=c[Ha>>2]|0;c[Ha>>2]=G;c[k>>2]=p;c[l>>2]=F;F=d+16|0;l=c[F>>2]|0;c[F>>2]=0;F=d+20|0;p=c[F>>2]|0;c[F>>2]=Ia;Ia=d+24|0;F=c[Ia>>2]|0;c[Ia>>2]=Ja;Ja=d+28|0;Ia=c[Ja>>2]|0;c[Ja>>2]=Ga;c[j>>2]=Ia;c[d+32>>2]=Ka;Mm(H);Mm(l);if((p|0)!=0){Vm(p)}if((F|0)==0){i=f;return}Vm(F);i=f;return}function Tc(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0;d=i;if((c[a+16>>2]|0)!=0){Ua(1344,704,258,1424)}e=c[a+28>>2]|0;f=e+b|0;b=a+32|0;if(!(f>>>0>(c[b>>2]|0)>>>0)){i=d;return}g=f>>>0>1073741823?-1:f<<2;h=Tm(g)|0;j=Tm(g)|0;g=a+20|0;k=c[g>>2]|0;l=(e>>>0<f>>>0?e:f)<<2;nn(h|0,k|0,l|0)|0;e=a+24|0;a=c[e>>2]|0;nn(j|0,a|0,l|0)|0;if((k|0)==0){m=a}else{Vm(k);m=c[e>>2]|0}if((m|0)!=0){Vm(m)}c[g>>2]=h;c[e>>2]=j;c[b>>2]=f;i=d;return}function Uc(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0;d=i;e=c[a+28>>2]|0;f=c[a>>2]|0;if((e|0)==0){h=f+(b<<2)|0;i=d;return h|0}if((e|0)!=1){Ua(1136,1160,190,1208)}e=a+40|0;j=c[e>>2]|0;k=a+32|0;if((j|0)==0){c[k>>2]=0;c[a+36>>2]=0;c[e>>2]=1;l=f+8|0;g[l>>2]=0.0;c[f+4>>2]=b;c[f>>2]=-1;h=l;i=d;return h|0}if((c[f+((c[k>>2]|0)*12|0)+4>>2]|0)>(b|0)){l=f+(j*12|0)+8|0;g[l>>2]=0.0;c[f+(j*12|0)+4>>2]=b;c[f+(j*12|0)>>2]=c[k>>2];m=c[e>>2]|0;c[k>>2]=m;c[e>>2]=m+1;c[a+36>>2]=m;h=l;i=d;return h|0}l=a+36|0;m=c[l>>2]|0;if((c[f+(m*12|0)+4>>2]|0)>(b|0)){Ua(1224,1160,217,1208)}k=c[f+(m*12|0)>>2]|0;a:do{if((k|0)>-1){n=m;o=k;while(1){if((c[f+(o*12|0)+4>>2]|0)>(b|0)){p=n;break a}c[l>>2]=o;q=c[f+(o*12|0)>>2]|0;if((q|0)>-1){r=o;o=q;n=r}else{p=o;break}}}else{p=m}}while(0);if((c[f+(p*12|0)+4>>2]|0)==(b|0)){h=f+(p*12|0)+8|0;i=d;return h|0}p=a+24|0;m=c[p>>2]|0;if((j|0)<(m|0)){s=j;t=f}else{k=~~(+(m|0)*1.5);o=c[a+8>>2]|0;n=(o|0)<(k|0)?o:k;c[p>>2]=n;p=Tm(n*12|0)|0;nn(p|0,f|0,m*12|0)|0;if((f|0)==0){u=j}else{Vm(f);u=c[e>>2]|0}c[a>>2]=p;s=u;t=p}p=t+(s*12|0)+8|0;g[p>>2]=0.0;c[t+(s*12|0)+4>>2]=b;c[t+(s*12|0)>>2]=c[t+((c[l>>2]|0)*12|0)>>2];c[t+((c[l>>2]|0)*12|0)>>2]=c[e>>2];c[e>>2]=(c[e>>2]|0)+1;h=p;i=d;return h|0}function Vc(a,b,d){a=a|0;b=b|0;d=+d;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;e=i;f=a+12|0;if(!((c[f>>2]|0)>>>0<b>>>0)){g=a+8|0;c[g>>2]=b;i=e;return}h=(~~(+(b>>>0)*d)>>>0)+b|0;j=h>>>0>1073741823?-1:h<<2;k=Tm(j)|0;l=Tm(j)|0;j=a+8|0;m=c[j>>2]|0;n=c[a>>2]|0;o=(m>>>0<h>>>0?m:h)<<2;nn(k|0,n|0,o|0)|0;m=a+4|0;p=c[m>>2]|0;nn(l|0,p|0,o|0)|0;if((n|0)==0){q=p}else{Vm(n);q=c[m>>2]|0}if((q|0)!=0){Vm(q)}c[a>>2]=k;c[m>>2]=l;c[f>>2]=h;g=j;c[g>>2]=b;i=e;return}function Wc(){var a=0;a=Wa(4)|0;Wm(a);xb(a|0,17136,109)}function Xc(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0;d=i;i=i+16|0;e=d;f=b+4|0;h=c[f>>2]|0;j=c[b+8>>2]|0;k=j<<2;l=k+4|0;m=Lm(l)|0;if((m|0)==0){Wc()}pn(m|0,0,l|0)|0;if(!((j|0)>-1)){Ua(1536,1552,138,1592)}l=(j|0)>0;if(l){pn(m|0,0,k|0)|0}if((h|0)>0){n=c[b+24>>2]|0;o=c[b+12>>2]|0;p=c[b+16>>2]|0;if((p|0)==0){q=0;do{r=c[o+(q<<2)>>2]|0;s=c[o+(q+1<<2)>>2]|0;if((r|0)<(s|0)){t=r;do{r=m+(c[n+(t<<2)>>2]<<2)|0;c[r>>2]=(c[r>>2]|0)+1;t=t+1|0}while((t|0)!=(s|0))}q=q+1|0}while((q|0)<(h|0))}else{q=0;do{s=c[o+(q<<2)>>2]|0;t=c[p+(q<<2)>>2]|0;r=t+s|0;if((t|0)>0){t=s;do{s=m+(c[n+(t<<2)>>2]<<2)|0;c[s>>2]=(c[s>>2]|0)+1;t=t+1|0}while((t|0)!=(r|0))}q=q+1|0}while((q|0)<(h|0))}}if(j>>>0>1073741823){Wc()}if((Qm(e,16,k)|0)==0){u=c[e>>2]|0}else{c[e>>2]=0;u=0}if(!((u|0)!=0|(k|0)==0)){Wc()}if(l){l=0;k=0;do{e=m+(k<<2)|0;q=c[e>>2]|0;c[e>>2]=l;c[u+(k<<2)>>2]=l;l=q+l|0;k=k+1|0}while((k|0)<(j|0));c[m+(j<<2)>>2]=l;if((l|0)==0){v=0;w=0;x=0;y=0}else{k=(~~(+(l>>>0)*0.0)>>>0)+l|0;q=k>>>0>1073741823?-1:k<<2;e=Tm(q)|0;v=e;w=l;x=k;y=Tm(q)|0}}else{c[m+(j<<2)>>2]=0;v=0;w=0;x=0;y=0}q=c[f>>2]|0;a:do{if((q|0)>0){k=c[b+20>>2]|0;l=c[b+24>>2]|0;e=c[b+12>>2]|0;n=c[b+16>>2]|0;if((n|0)==0){p=q;o=0;b:while(1){r=c[e+(o<<2)>>2]|0;t=c[e+(o+1<<2)>>2]|0;if((r|0)<(t|0)){s=r;do{r=c[l+(s<<2)>>2]|0;if(!((r|0)>-1&(j|0)>(r|0))){break b}z=u+(r<<2)|0;r=c[z>>2]|0;c[z>>2]=r+1;c[y+(r<<2)>>2]=o;g[v+(r<<2)>>2]=+g[k+(s<<2)>>2];s=s+1|0}while((s|0)<(t|0));A=c[f>>2]|0}else{A=p}o=o+1|0;if((o|0)>=(A|0)){break a}else{p=A}}Ua(1600,1632,378,1680)}else{p=q;o=0;c:while(1){t=c[e+(o<<2)>>2]|0;s=c[n+(o<<2)>>2]|0;r=s+t|0;if((s|0)>0){s=t;do{t=c[l+(s<<2)>>2]|0;if(!((t|0)>-1&(j|0)>(t|0))){break c}z=u+(t<<2)|0;t=c[z>>2]|0;c[z>>2]=t+1;c[y+(t<<2)>>2]=o;g[v+(t<<2)>>2]=+g[k+(s<<2)>>2];s=s+1|0}while((s|0)<(r|0));B=c[f>>2]|0}else{B=p}o=o+1|0;if((o|0)>=(B|0)){break a}else{p=B}}Ua(1600,1632,378,1680)}}}while(0);B=a+12|0;f=c[B>>2]|0;c[B>>2]=m;c[a+8>>2]=h;c[a+4>>2]=j;j=a+16|0;h=c[j>>2]|0;c[j>>2]=0;j=a+20|0;m=c[j>>2]|0;c[j>>2]=v;v=a+24|0;j=c[v>>2]|0;c[v>>2]=y;c[a+28>>2]=w;c[a+32>>2]=x;Mm(u);Mm(f);Mm(h);if((m|0)!=0){Vm(m)}if((j|0)==0){i=d;return a|0}Vm(j);i=d;return a|0}function Yc(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0;d=i;i=i+16|0;e=d;f=c[b+8>>2]|0;h=b+4|0;j=c[h>>2]|0;k=f<<2;l=k+4|0;m=Lm(l)|0;if((m|0)==0){Wc()}pn(m|0,0,l|0)|0;if(!((f|0)>-1)){Ua(1536,1552,138,1592)}l=(f|0)>0;if(l){pn(m|0,0,k|0)|0}if((j|0)>0){n=c[b+24>>2]|0;o=c[b+12>>2]|0;p=c[b+16>>2]|0;q=(p|0)==0;r=0;do{s=c[o+(r<<2)>>2]|0;if(q){t=c[o+(r+1<<2)>>2]|0}else{t=(c[p+(r<<2)>>2]|0)+s|0}if((s|0)<(t|0)){u=s;do{s=m+(c[n+(u<<2)>>2]<<2)|0;c[s>>2]=(c[s>>2]|0)+1;u=u+1|0}while((u|0)!=(t|0))}r=r+1|0}while((r|0)<(j|0))}if(f>>>0>1073741823){Wc()}if((Qm(e,16,k)|0)==0){v=c[e>>2]|0}else{c[e>>2]=0;v=0}if(!((v|0)!=0|(k|0)==0)){Wc()}if(l){l=0;k=0;do{e=m+(k<<2)|0;r=c[e>>2]|0;c[e>>2]=l;c[v+(k<<2)>>2]=l;l=r+l|0;k=k+1|0}while((k|0)<(f|0));c[m+(f<<2)>>2]=l;if((l|0)==0){w=0;x=0;y=0;z=0}else{k=(~~(+(l>>>0)*0.0)>>>0)+l|0;r=k>>>0>1073741823?-1:k<<2;e=Tm(r)|0;w=l;x=k;y=Tm(r)|0;z=e}}else{c[m+(f<<2)>>2]=0;w=0;x=0;y=0;z=0}e=c[h>>2]|0;a:do{if((e|0)>0){r=c[b+20>>2]|0;k=c[b+24>>2]|0;l=c[b+12>>2]|0;t=c[b+16>>2]|0;n=(t|0)==0;p=e;o=0;b:while(1){q=c[l+(o<<2)>>2]|0;if(n){A=c[l+(o+1<<2)>>2]|0}else{A=(c[t+(o<<2)>>2]|0)+q|0}if((q|0)<(A|0)){u=q;do{q=c[k+(u<<2)>>2]|0;if(!((q|0)>-1&(f|0)>(q|0))){break b}s=v+(q<<2)|0;q=c[s>>2]|0;c[s>>2]=q+1;c[y+(q<<2)>>2]=o;g[z+(q<<2)>>2]=+g[r+(u<<2)>>2];u=u+1|0}while((u|0)<(A|0));B=c[h>>2]|0}else{B=p}o=o+1|0;if((o|0)>=(B|0)){break a}else{p=B}}Ua(1600,1632,378,1680)}}while(0);B=a+12|0;h=c[B>>2]|0;c[B>>2]=m;c[a+8>>2]=j;c[a+4>>2]=f;f=a+16|0;j=c[f>>2]|0;c[f>>2]=0;f=a+20|0;m=c[f>>2]|0;c[f>>2]=z;z=a+24|0;f=c[z>>2]|0;c[z>>2]=y;c[a+28>>2]=w;c[a+32>>2]=x;Mm(v);Mm(h);Mm(j);if((m|0)!=0){Vm(m)}if((f|0)==0){i=d;return a|0}Vm(f);i=d;return a|0}function Zc(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0,k=0;f=i;c[b>>2]=0;g=b+4|0;c[g>>2]=0;h=b+8|0;c[h>>2]=0;if((d|0)==0){i=f;return}if((d|0)<0){Sj(0)}j=((d+ -1|0)>>>5)+1|0;k=Sm(j<<2)|0;c[b>>2]=k;c[h>>2]=j;j=(a[e]|0)==0;c[g>>2]=d;g=d>>>5;e=g<<2;if(j){pn(k|0,0,e|0)|0;j=g<<5;if((j|0)==(d|0)){i=f;return}h=k+(g<<2)|0;c[h>>2]=c[h>>2]&~(-1>>>(32-d+j|0));i=f;return}else{pn(k|0,-1,e|0)|0;e=g<<5;if((e|0)==(d|0)){i=f;return}j=k+(g<<2)|0;c[j>>2]=c[j>>2]|-1>>>(32-d+e|0);i=f;return}}function _c(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,h=0,j=0;f=i;a[b]=0;h=b+4|0;c[h>>2]=d;d=b+8|0;a[d]=0;j=b+12|0;c[j+0>>2]=0;c[j+4>>2]=0;c[j+8>>2]=0;c[j+12>>2]=0;c[j+16>>2]=0;c[j+20>>2]=0;c[j+24>>2]=0;c[j+28>>2]=0;$c(d,e)|0;g[b+44>>2]=0.0;a[b+48|0]=1;if((c[(c[h>>2]|0)+4>>2]|0)==(c[b+16>>2]|0)){i=f;return}else{Ua(1976,2008,133,2056)}}function $c(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0;d=i;i=i+16|0;e=d;f=c[b+4>>2]|0;h=c[b+4>>2]|0;b=c[h+4>>2]|0;j=c[h+8>>2]|0;h=j<<2;k=h+4|0;l=Lm(k)|0;if((l|0)==0){Wc()}pn(l|0,0,k|0)|0;if(!((j|0)>-1)){Ua(1536,1552,138,1592)}k=(j|0)>0;if(k){pn(l|0,0,h|0)|0}m=f+4|0;n=c[m>>2]|0;if((n|0)>0){o=c[f+24>>2]|0;p=c[f+12>>2]|0;q=c[f+16>>2]|0;r=(q|0)==0;s=n;n=0;while(1){t=c[p+(n<<2)>>2]|0;if(r){u=c[p+(n+1<<2)>>2]|0}else{u=(c[q+(n<<2)>>2]|0)+t|0}if((t|0)<(u|0)){v=t;do{t=l+(c[o+(v<<2)>>2]<<2)|0;c[t>>2]=(c[t>>2]|0)+1;v=v+1|0}while((v|0)!=(u|0));w=c[m>>2]|0}else{w=s}n=n+1|0;if((n|0)>=(w|0)){break}else{s=w}}}if(j>>>0>1073741823){Wc()}if((Qm(e,16,h)|0)==0){x=c[e>>2]|0}else{c[e>>2]=0;x=0}if(!((x|0)!=0|(h|0)==0)){Wc()}if(k){k=0;h=0;do{e=l+(h<<2)|0;w=c[e>>2]|0;c[e>>2]=k;c[x+(h<<2)>>2]=k;k=w+k|0;h=h+1|0}while((h|0)<(j|0));c[l+(j<<2)>>2]=k;if((k|0)==0){y=0;z=0;A=0;B=0}else{h=(~~(+(k>>>0)*0.0)>>>0)+k|0;w=h>>>0>1073741823?-1:h<<2;e=Tm(w)|0;y=Tm(w)|0;z=e;A=k;B=h}}else{c[l+(j<<2)>>2]=0;y=0;z=0;A=0;B=0}h=c[m>>2]|0;a:do{if((h|0)>0){k=c[f+20>>2]|0;e=c[f+24>>2]|0;w=c[f+12>>2]|0;s=c[f+16>>2]|0;n=(s|0)==0;u=h;o=0;b:while(1){q=c[w+(o<<2)>>2]|0;if(n){C=c[w+(o+1<<2)>>2]|0}else{C=(c[s+(o<<2)>>2]|0)+q|0}if((q|0)<(C|0)){p=q;do{q=c[e+(p<<2)>>2]|0;if(!((q|0)>-1&(j|0)>(q|0))){break b}r=x+(q<<2)|0;q=c[r>>2]|0;c[r>>2]=q+1;c[y+(q<<2)>>2]=o;g[z+(q<<2)>>2]=+g[k+(p<<2)>>2];p=p+1|0}while((p|0)<(C|0));D=c[m>>2]|0}else{D=u}o=o+1|0;if((o|0)>=(D|0)){break a}else{u=D}}Ua(1600,1632,378,1680)}}while(0);D=a+12|0;m=c[D>>2]|0;c[D>>2]=l;c[a+8>>2]=b;c[a+4>>2]=j;j=a+16|0;b=c[j>>2]|0;c[j>>2]=0;j=a+20|0;l=c[j>>2]|0;c[j>>2]=z;z=a+24|0;j=c[z>>2]|0;c[z>>2]=y;c[a+28>>2]=A;c[a+32>>2]=B;Mm(x);Mm(m);Mm(b);if((l|0)!=0){Vm(l)}if((j|0)==0){i=d;return a|0}Vm(j);i=d;return a|0}function ad(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0;e=i;if((a[d]|0)!=0){f=c[d+4>>2]|0;g=c[f+4>>2]|0;c[b+8>>2]=c[f+8>>2];c[b+28>>2]=0;f=b+4|0;h=c[f>>2]|0;do{if((h|0)!=(g|0)|(h|0)==0){j=b+12|0;Mm(c[j>>2]|0);k=Lm((g<<2)+4|0)|0;c[j>>2]=k;if((k|0)==0){Wc()}else{c[f>>2]=g;break}}}while(0);h=b+16|0;k=c[h>>2]|0;if((k|0)==0){l=g}else{Mm(k);c[h>>2]=0;l=c[f>>2]|0}pn(c[b+12>>2]|0,0,(l<<2)+4|0)|0;l=c[h>>2]|0;if((l|0)!=0){Mm(l);c[h>>2]=0}}h=bd(b,d)|0;i=e;return h|0}function bd(b,d){b=b|0;d=d|0;var e=0,f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0.0,D=0,E=0,F=0,G=0;e=i;f=d+4|0;h=c[f>>2]|0;j=c[h+4>>2]|0;if((a[d]|0)==0){cd(b,d);k=b;i=e;return k|0}d=b+8|0;c[d>>2]=c[h+8>>2];h=b+28|0;c[h>>2]=0;l=b+4|0;m=c[l>>2]|0;do{if((m|0)!=(j|0)|(m|0)==0){n=b+12|0;Mm(c[n>>2]|0);o=Lm((j<<2)+4|0)|0;c[n>>2]=o;if((o|0)==0){Wc()}else{c[l>>2]=j;break}}}while(0);m=b+16|0;o=c[m>>2]|0;if((o|0)==0){p=j}else{Mm(o);c[m>>2]=0;p=c[l>>2]|0}o=b+12|0;pn(c[o>>2]|0,0,(p<<2)+4|0)|0;c[h>>2]=0;pn(c[o>>2]|0,0,(c[l>>2]<<2)+4|0)|0;p=c[m>>2]|0;if((p|0)!=0){pn(p|0,0,c[l>>2]<<2|0)|0}p=c[d>>2]|0;d=c[l>>2]|0;Tc(b,((p|0)<(d|0)?d:p)<<1);a:do{if((j|0)>0){p=b+20|0;d=b+24|0;n=0;b:while(1){q=c[o>>2]|0;r=c[q+(n<<2)>>2]|0;if((r|0)!=(c[h>>2]|0)){s=14;break}t=n+1|0;u=q+(t<<2)|0;if((c[u>>2]|0)!=0){s=16;break}c[u>>2]=r;u=c[f>>2]|0;v=c[u+20>>2]|0;w=c[u+24>>2]|0;x=c[u+12>>2]|0;y=c[x+(n<<2)>>2]|0;z=c[u+16>>2]|0;if((z|0)==0){A=c[x+(t<<2)>>2]|0}else{A=(c[z+(n<<2)>>2]|0)+y|0}c:do{if((y|0)<(A|0)){z=y;while(1){if((c[w+(z<<2)>>2]|0)>=(n|0)){B=z;break c}x=z+1|0;if((x|0)<(A|0)){z=x}else{B=x;break}}}else{B=y}}while(0);d:do{if((B|0)<(A|0)){y=B;z=q;x=r;while(1){C=+g[v+(y<<2)>>2];u=c[w+(y<<2)>>2]|0;D=z+(t<<2)|0;if((x|0)!=(c[h>>2]|0)){s=25;break b}if((x|0)!=(c[z+(n<<2)>>2]|0)?(c[(c[d>>2]|0)+(x+ -1<<2)>>2]|0)>=(u|0):0){s=28;break b}c[D>>2]=x+1;D=c[h>>2]|0;Vc(p,D+1|0,1.0);E=c[p>>2]|0;g[E+(D<<2)>>2]=0.0;c[(c[d>>2]|0)+(D<<2)>>2]=u;g[E+(x<<2)>>2]=C;E=y+1|0;if((E|0)>=(A|0)){break d}u=c[o>>2]|0;y=E;z=u;x=c[u+(t<<2)>>2]|0}}}while(0);if((t|0)<(j|0)){n=t}else{break a}}if((s|0)==14){Ua(920,704,405,1032)}else if((s|0)==16){Ua(1048,704,406,1032)}else if((s|0)==25){Ua(600,704,383,752)}else if((s|0)==28){Ua(776,704,384,752)}}}while(0);if((c[m>>2]|0)!=0){k=b;i=e;return k|0}m=c[h>>2]|0;h=c[l>>2]|0;if(!((h|0)>-1)){k=b;i=e;return k|0}s=c[o>>2]|0;o=h;while(1){j=o+ -1|0;if((c[s+(o<<2)>>2]|0)!=0){F=o;break}if((o|0)>0){o=j}else{F=j;break}}if((F|0)<(h|0)){G=F}else{k=b;i=e;return k|0}while(1){F=G+1|0;c[s+(F<<2)>>2]=m;if((F|0)<(c[l>>2]|0)){G=F}else{k=b;break}}i=e;return k|0}function cd(b,d){b=b|0;d=d|0;var e=0,f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0.0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0;e=i;i=i+48|0;f=e;h=c[d+4>>2]|0;j=c[h+4>>2]|0;k=d+4|0;d=c[h+8>>2]|0;a[f]=0;h=f+4|0;c[h+0>>2]=0;c[h+4>>2]=0;c[h+8>>2]=0;c[h+12>>2]=0;c[h+16>>2]=0;c[h+20>>2]=0;c[h+24>>2]=0;c[h+28>>2]=0;c[f+8>>2]=d;d=f+28|0;c[d>>2]=0;l=f+12|0;m=Lm((j<<2)+4|0)|0;c[l>>2]=m;if((m|0)==0){Wc()}c[h>>2]=j;n=f+16|0;pn(m|0,0,(j<<2)+4|0)|0;m=c[b+8>>2]|0;o=c[b+4>>2]|0;if((c[n>>2]|0)!=0){Ua(1344,704,258,1424)}p=c[d>>2]|0;q=p+(((m|0)<(o|0)?o:m)<<1)|0;m=f+32|0;if(q>>>0>(c[m>>2]|0)>>>0){o=q>>>0>1073741823?-1:q<<2;r=Tm(o)|0;s=Tm(o)|0;o=f+20|0;t=c[o>>2]|0;u=(p>>>0<q>>>0?p:q)<<2;nn(r|0,t|0,u|0)|0;p=f+24|0;v=c[p>>2]|0;nn(s|0,v|0,u|0)|0;if((t|0)==0){w=v}else{Vm(t);w=c[p>>2]|0}if((w|0)!=0){Vm(w)}c[o>>2]=r;c[p>>2]=s;c[m>>2]=q}a:do{if((j|0)>0){q=f+24|0;s=f+20|0;p=0;b:while(1){r=c[l>>2]|0;o=c[r+(p<<2)>>2]|0;if((o|0)!=(c[d>>2]|0)){x=15;break}w=p+1|0;t=r+(w<<2)|0;if((c[t>>2]|0)!=0){x=17;break}c[t>>2]=o;t=c[k>>2]|0;v=c[t+20>>2]|0;u=c[t+24>>2]|0;y=c[t+12>>2]|0;z=c[y+(p<<2)>>2]|0;A=c[t+16>>2]|0;if((A|0)==0){B=c[y+(w<<2)>>2]|0}else{B=(c[A+(p<<2)>>2]|0)+z|0}c:do{if((z|0)<(B|0)){A=z;while(1){if((c[u+(A<<2)>>2]|0)>=(p|0)){C=A;break c}y=A+1|0;if((y|0)<(B|0)){A=y}else{C=y;break}}}else{C=z}}while(0);d:do{if((C|0)<(B|0)){z=C;A=r;y=o;while(1){D=+g[v+(z<<2)>>2];t=c[u+(z<<2)>>2]|0;E=A+(w<<2)|0;if((y|0)!=(c[d>>2]|0)){x=26;break b}if((y|0)!=(c[A+(p<<2)>>2]|0)?(c[(c[q>>2]|0)+(y+ -1<<2)>>2]|0)>=(t|0):0){x=29;break b}c[E>>2]=y+1;E=c[d>>2]|0;F=E+1|0;if((c[m>>2]|0)>>>0<F>>>0){G=(~~+(F>>>0)>>>0)+F|0;H=G>>>0>1073741823?-1:G<<2;I=Tm(H)|0;J=Tm(H)|0;H=c[s>>2]|0;K=(E>>>0<G>>>0?E:G)<<2;nn(I|0,H|0,K|0)|0;L=c[q>>2]|0;nn(J|0,L|0,K|0)|0;if((H|0)==0){M=L}else{Vm(H);M=c[q>>2]|0}if((M|0)!=0){Vm(M)}c[s>>2]=I;c[q>>2]=J;c[m>>2]=G;N=I;O=J}else{N=c[s>>2]|0;O=c[q>>2]|0}c[d>>2]=F;g[N+(E<<2)>>2]=0.0;c[O+(E<<2)>>2]=t;g[N+(y<<2)>>2]=D;t=z+1|0;if((t|0)>=(B|0)){break d}E=c[l>>2]|0;z=t;A=E;y=c[E+(w<<2)>>2]|0}}}while(0);if((w|0)<(j|0)){p=w}else{break a}}if((x|0)==15){Ua(920,704,405,1032)}else if((x|0)==17){Ua(1048,704,406,1032)}else if((x|0)==26){Ua(600,704,383,752)}else if((x|0)==29){Ua(776,704,384,752)}}}while(0);if((c[n>>2]|0)==0?(x=c[d>>2]|0,d=c[h>>2]|0,(d|0)>-1):0){j=c[l>>2]|0;B=d;while(1){N=B+ -1|0;if((c[j+(B<<2)>>2]|0)!=0){P=B;break}if((B|0)>0){B=N}else{P=N;break}}if((P|0)<(d|0)){d=P;do{d=d+1|0;c[j+(d<<2)>>2]=x}while((d|0)<(c[h>>2]|0))}}a[f]=1;dd(b,f)|0;Mm(c[l>>2]|0);Mm(c[n>>2]|0);n=c[f+20>>2]|0;if((n|0)!=0){Vm(n)}n=c[f+24>>2]|0;if((n|0)==0){i=e;return}Vm(n);i=e;return}function dd(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0;e=i;if((a[d]|0)!=0){f=b+12|0;g=d+12|0;h=c[f>>2]|0;c[f>>2]=c[g>>2];c[g>>2]=h;h=b+8|0;g=d+8|0;f=c[h>>2]|0;c[h>>2]=c[g>>2];c[g>>2]=f;f=b+4|0;g=d+4|0;h=c[f>>2]|0;c[f>>2]=c[g>>2];c[g>>2]=h;h=b+16|0;g=d+16|0;f=c[h>>2]|0;c[h>>2]=c[g>>2];c[g>>2]=f;f=b+20|0;g=d+20|0;h=c[f>>2]|0;c[f>>2]=c[g>>2];c[g>>2]=h;h=b+24|0;g=d+24|0;f=c[h>>2]|0;c[h>>2]=c[g>>2];c[g>>2]=f;f=b+28|0;g=d+28|0;h=c[f>>2]|0;c[f>>2]=c[g>>2];c[g>>2]=h;h=b+32|0;g=d+32|0;f=c[h>>2]|0;c[h>>2]=c[g>>2];c[g>>2]=f;i=e;return b|0}if((b|0)==(d|0)){i=e;return b|0}f=c[d+4>>2]|0;c[b+8>>2]=c[d+8>>2];g=b+28|0;c[g>>2]=0;h=b+4|0;j=c[h>>2]|0;do{if((j|0)!=(f|0)|(j|0)==0){k=b+12|0;Mm(c[k>>2]|0);l=Lm((f<<2)+4|0)|0;c[k>>2]=l;if((l|0)==0){Wc()}else{c[h>>2]=f;break}}}while(0);j=b+16|0;l=c[j>>2]|0;if((l|0)==0){m=f}else{Mm(l);c[j>>2]=0;m=c[h>>2]|0}l=b+12|0;pn(c[l>>2]|0,0,(m<<2)+4|0)|0;m=c[j>>2]|0;if((m|0)!=0){Mm(m);c[j>>2]=0}if((c[d+16>>2]|0)==0){nn(c[l>>2]|0,c[d+12>>2]|0,(c[h>>2]<<2)+4|0)|0;h=b+20|0;Vc(h,c[d+28>>2]|0,0.0);nn(c[h>>2]|0,c[d+20>>2]|0,c[g>>2]<<2|0)|0;nn(c[b+24>>2]|0,c[d+24>>2]|0,c[g>>2]<<2|0)|0;i=e;return b|0}else{ed(b,d)|0;i=e;return b|0}return 0}function ed(b,d){b=b|0;d=d|0;var e=0,f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0.0,G=0,H=0,I=0,J=0,K=0;e=i;f=c[d+4>>2]|0;if((a[d]|0)==0){fd(b,d);h=b;i=e;return h|0}j=b+8|0;c[j>>2]=c[d+8>>2];k=b+28|0;c[k>>2]=0;l=b+4|0;m=c[l>>2]|0;do{if((m|0)!=(f|0)|(m|0)==0){n=b+12|0;Mm(c[n>>2]|0);o=Lm((f<<2)+4|0)|0;c[n>>2]=o;if((o|0)==0){Wc()}else{c[l>>2]=f;break}}}while(0);m=b+16|0;o=c[m>>2]|0;if((o|0)==0){p=f}else{Mm(o);c[m>>2]=0;p=c[l>>2]|0}o=b+12|0;pn(c[o>>2]|0,0,(p<<2)+4|0)|0;c[k>>2]=0;pn(c[o>>2]|0,0,(c[l>>2]<<2)+4|0)|0;p=c[m>>2]|0;if((p|0)!=0){pn(p|0,0,c[l>>2]<<2|0)|0}p=c[j>>2]|0;j=c[l>>2]|0;Tc(b,((p|0)<(j|0)?j:p)<<1);a:do{if((f|0)>0){p=d+20|0;j=d+24|0;n=d+12|0;q=d+16|0;r=b+20|0;s=b+24|0;t=0;b:while(1){u=c[o>>2]|0;v=c[u+(t<<2)>>2]|0;if((v|0)!=(c[k>>2]|0)){w=14;break}x=t+1|0;y=u+(x<<2)|0;if((c[y>>2]|0)!=0){w=16;break}c[y>>2]=v;y=c[p>>2]|0;z=c[j>>2]|0;A=c[n>>2]|0;B=c[A+(t<<2)>>2]|0;C=c[q>>2]|0;if((C|0)==0){D=c[A+(x<<2)>>2]|0}else{D=(c[C+(t<<2)>>2]|0)+B|0}c:do{if((B|0)<(D|0)){C=B;A=u;E=v;while(1){F=+g[y+(C<<2)>>2];G=c[z+(C<<2)>>2]|0;H=A+(x<<2)|0;if((E|0)!=(c[k>>2]|0)){w=22;break b}if((E|0)!=(c[A+(t<<2)>>2]|0)?(c[(c[s>>2]|0)+(E+ -1<<2)>>2]|0)>=(G|0):0){w=25;break b}c[H>>2]=E+1;H=c[k>>2]|0;Vc(r,H+1|0,1.0);I=c[r>>2]|0;g[I+(H<<2)>>2]=0.0;c[(c[s>>2]|0)+(H<<2)>>2]=G;g[I+(E<<2)>>2]=F;I=C+1|0;if((I|0)>=(D|0)){break c}G=c[o>>2]|0;C=I;A=G;E=c[G+(x<<2)>>2]|0}}}while(0);if((x|0)<(f|0)){t=x}else{break a}}if((w|0)==14){Ua(920,704,405,1032)}else if((w|0)==16){Ua(1048,704,406,1032)}else if((w|0)==22){Ua(600,704,383,752)}else if((w|0)==25){Ua(776,704,384,752)}}}while(0);if((c[m>>2]|0)!=0){h=b;i=e;return h|0}m=c[k>>2]|0;k=c[l>>2]|0;if(!((k|0)>-1)){h=b;i=e;return h|0}w=c[o>>2]|0;o=k;while(1){f=o+ -1|0;if((c[w+(o<<2)>>2]|0)!=0){J=o;break}if((o|0)>0){o=f}else{J=f;break}}if((J|0)<(k|0)){K=J}else{h=b;i=e;return h|0}while(1){J=K+1|0;c[w+(J<<2)>>2]=m;if((J|0)<(c[l>>2]|0)){K=J}else{h=b;break}}i=e;return h|0}function fd(b,d){b=b|0;d=d|0;var e=0,f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0.0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0;e=i;i=i+48|0;f=e;h=c[d+4>>2]|0;j=c[d+8>>2]|0;a[f]=0;k=f+4|0;c[k+0>>2]=0;c[k+4>>2]=0;c[k+8>>2]=0;c[k+12>>2]=0;c[k+16>>2]=0;c[k+20>>2]=0;c[k+24>>2]=0;c[k+28>>2]=0;c[f+8>>2]=j;j=f+28|0;c[j>>2]=0;l=f+12|0;m=Lm((h<<2)+4|0)|0;c[l>>2]=m;if((m|0)==0){Wc()}c[k>>2]=h;n=f+16|0;pn(m|0,0,(h<<2)+4|0)|0;m=c[b+8>>2]|0;o=c[b+4>>2]|0;if((c[n>>2]|0)!=0){Ua(1344,704,258,1424)}p=c[j>>2]|0;q=p+(((m|0)<(o|0)?o:m)<<1)|0;m=f+32|0;if(q>>>0>(c[m>>2]|0)>>>0){o=q>>>0>1073741823?-1:q<<2;r=Tm(o)|0;s=Tm(o)|0;o=f+20|0;t=c[o>>2]|0;u=(p>>>0<q>>>0?p:q)<<2;nn(r|0,t|0,u|0)|0;p=f+24|0;v=c[p>>2]|0;nn(s|0,v|0,u|0)|0;if((t|0)==0){w=v}else{Vm(t);w=c[p>>2]|0}if((w|0)!=0){Vm(w)}c[o>>2]=r;c[p>>2]=s;c[m>>2]=q}a:do{if((h|0)>0){q=d+20|0;s=d+24|0;p=d+12|0;r=d+16|0;o=f+24|0;w=f+20|0;t=0;b:while(1){v=c[l>>2]|0;u=c[v+(t<<2)>>2]|0;if((u|0)!=(c[j>>2]|0)){x=15;break}y=t+1|0;z=v+(y<<2)|0;if((c[z>>2]|0)!=0){x=17;break}c[z>>2]=u;z=c[q>>2]|0;A=c[s>>2]|0;B=c[p>>2]|0;C=c[B+(t<<2)>>2]|0;D=c[r>>2]|0;if((D|0)==0){E=c[B+(y<<2)>>2]|0}else{E=(c[D+(t<<2)>>2]|0)+C|0}c:do{if((C|0)<(E|0)){D=C;B=v;F=u;while(1){G=+g[z+(D<<2)>>2];H=c[A+(D<<2)>>2]|0;I=B+(y<<2)|0;if((F|0)!=(c[j>>2]|0)){x=23;break b}if((F|0)!=(c[B+(t<<2)>>2]|0)?(c[(c[o>>2]|0)+(F+ -1<<2)>>2]|0)>=(H|0):0){x=26;break b}c[I>>2]=F+1;I=c[j>>2]|0;J=I+1|0;if((c[m>>2]|0)>>>0<J>>>0){K=(~~+(J>>>0)>>>0)+J|0;L=K>>>0>1073741823?-1:K<<2;M=Tm(L)|0;N=Tm(L)|0;L=c[w>>2]|0;O=(I>>>0<K>>>0?I:K)<<2;nn(M|0,L|0,O|0)|0;P=c[o>>2]|0;nn(N|0,P|0,O|0)|0;if((L|0)==0){Q=P}else{Vm(L);Q=c[o>>2]|0}if((Q|0)!=0){Vm(Q)}c[w>>2]=M;c[o>>2]=N;c[m>>2]=K;R=M;S=N}else{R=c[w>>2]|0;S=c[o>>2]|0}c[j>>2]=J;g[R+(I<<2)>>2]=0.0;c[S+(I<<2)>>2]=H;g[R+(F<<2)>>2]=G;H=D+1|0;if((H|0)>=(E|0)){break c}I=c[l>>2]|0;D=H;B=I;F=c[I+(y<<2)>>2]|0}}}while(0);if((y|0)<(h|0)){t=y}else{break a}}if((x|0)==15){Ua(920,704,405,1032)}else if((x|0)==17){Ua(1048,704,406,1032)}else if((x|0)==23){Ua(600,704,383,752)}else if((x|0)==26){Ua(776,704,384,752)}}}while(0);if((c[n>>2]|0)==0?(x=c[j>>2]|0,j=c[k>>2]|0,(j|0)>-1):0){h=c[l>>2]|0;E=j;while(1){R=E+ -1|0;if((c[h+(E<<2)>>2]|0)!=0){T=E;break}if((E|0)>0){E=R}else{T=R;break}}if((T|0)<(j|0)){j=T;do{j=j+1|0;c[h+(j<<2)>>2]=x}while((j|0)<(c[k>>2]|0))}}a[f]=1;dd(b,f)|0;Mm(c[l>>2]|0);Mm(c[n>>2]|0);n=c[f+20>>2]|0;if((n|0)!=0){Vm(n)}n=c[f+24>>2]|0;if((n|0)==0){i=e;return}Vm(n);i=e;return}function gd(b,d){b=b|0;d=d|0;var e=0,f=0,h=0,j=0;e=i;c[b>>2]=0;a[b+4|0]=0;a[b+8|0]=0;f=b+12|0;c[f>>2]=-1;h=b+16|0;j=b+20|0;c[h+0>>2]=0;c[h+4>>2]=0;c[h+8>>2]=0;c[h+12>>2]=0;c[h+16>>2]=0;c[h+20>>2]=0;c[h+24>>2]=0;h=Lm(4)|0;c[j>>2]=h;if((h|0)==0){Wc()}else{c[f>>2]=0;c[h>>2]=0;h=b+88|0;f=b+44|0;j=f+44|0;do{c[f>>2]=0;f=f+4|0}while((f|0)<(j|0));g[h>>2]=1.0;hd(b,d);i=e;return}}function hd(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0;e=i;i=i+48|0;f=e;g=c[d+8>>2]|0;if((g|0)!=(c[d+4>>2]|0)){Ua(2192,2128,186,2224)}a[f]=0;h=f+4|0;c[h+0>>2]=0;c[h+4>>2]=0;c[h+8>>2]=0;c[h+12>>2]=0;c[h+16>>2]=0;c[h+20>>2]=0;c[h+24>>2]=0;c[h+28>>2]=0;c[f+8>>2]=g;c[f+28>>2]=0;j=f+12|0;k=Lm((g<<2)+4|0)|0;c[j>>2]=k;if((k|0)==0){Wc()}c[h>>2]=g;pn(k|0,0,(g<<2)+4|0)|0;id(b,d,f);jd(b,f,0);kd(b,f);Mm(c[j>>2]|0);Mm(c[f+16>>2]|0);j=c[f+20>>2]|0;if((j|0)!=0){Vm(j)}j=c[f+24>>2]|0;if((j|0)==0){i=e;return}Vm(j);i=e;return}function id(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;f=i;i=i+48|0;g=f;h=c[d+8>>2]|0;if((h|0)!=(c[d+4>>2]|0)){Ua(2664,2128,613,2688)}a[g]=0;j=g+4|0;c[j>>2]=-1;k=g+8|0;l=g+12|0;c[k+0>>2]=0;c[k+4>>2]=0;c[k+8>>2]=0;c[k+12>>2]=0;c[k+16>>2]=0;c[k+20>>2]=0;c[k+24>>2]=0;k=Lm(4)|0;c[l>>2]=k;if((k|0)==0){Wc()}c[j>>2]=0;c[k>>2]=0;k=d+8|0;j=d+4|0;if((c[k>>2]|0)!=(c[j>>2]|0)){Ua(3120,3192,62,3248)}qd(d,g,0);m=b+76|0;md(g,m);Mm(c[l>>2]|0);Mm(c[g+16>>2]|0);l=c[g+20>>2]|0;if((l|0)!=0){Vm(l)}l=c[g+24>>2]|0;if((l|0)!=0){Vm(l)}l=c[b+80>>2]|0;n=b+68|0;if((l|0)>0){if(l>>>0>1073741823){Wc()}o=l<<2;if((Qm(g,16,o)|0)==0){p=c[g>>2]|0}else{c[g>>2]=0;p=0}if(!((p|0)!=0|(o|0)==0)){Wc()}if(!((l|0)>-1)){Ua(1696,1712,225,1752)}o=c[m>>2]|0;m=0;do{c[p+(c[o+(m<<2)>>2]<<2)>>2]=m;m=m+1|0}while((m|0)!=(l|0));od(n,l,1);if((c[b+72>>2]|0)!=(l|0)){Ua(1760,1816,498,1856)}m=c[n>>2]|0;o=0;do{c[m+(o<<2)>>2]=c[p+(o<<2)>>2];o=o+1|0}while((o|0)!=(l|0));Mm(p)}else{p=b+72|0;if((c[p>>2]|0)!=0){Mm(c[n>>2]|0);c[n>>2]=0}c[p>>2]=0}c[e+8>>2]=h;c[e+28>>2]=0;p=e+4|0;b=c[p>>2]|0;do{if((b|0)!=(h|0)|(b|0)==0){l=e+12|0;Mm(c[l>>2]|0);o=Lm((h<<2)+4|0)|0;c[l>>2]=o;if((o|0)==0){Wc()}else{c[p>>2]=h;break}}}while(0);b=e+16|0;o=c[b>>2]|0;if((o|0)==0){q=h}else{Mm(o);c[b>>2]=0;q=c[p>>2]|0}pn(c[e+12>>2]|0,0,(q<<2)+4|0)|0;if((c[e+8>>2]|0)!=(c[e+4>>2]|0)){Ua(3120,3192,62,3248)}if((c[k>>2]|0)==(c[j>>2]|0)){nd(d,e,c[n>>2]|0);i=f;return}else{Ua(3120,3192,62,3248)}}function jd(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0;f=i;i=i+16|0;g=f;h=c[d+8>>2]|0;c[b+16>>2]=h;c[b+36>>2]=0;j=b+12|0;k=c[j>>2]|0;do{if((k|0)!=(h|0)|(k|0)==0){l=b+20|0;Mm(c[l>>2]|0);m=Lm((h<<2)+4|0)|0;c[l>>2]=m;if((m|0)==0){Wc()}else{c[j>>2]=h;break}}}while(0);k=b+24|0;m=c[k>>2]|0;if((m|0)==0){n=h}else{Mm(m);c[k>>2]=0;n=c[j>>2]|0}j=b+20|0;pn(c[j>>2]|0,0,(n<<2)+4|0)|0;n=b+52|0;ld(n,h);k=b+60|0;ld(k,h);if(h>>>0>1073741823){Wc()}m=h<<2;if((Qm(g,16,m)|0)==0){o=c[g>>2]|0}else{c[g>>2]=0;o=0}if(!((o|0)!=0|(m|0)==0)){Wc()}m=(h|0)>0;if(!m){g=c[j>>2]|0;c[g>>2]=0;p=g;q=p+(h<<2)|0;r=c[q>>2]|0;s=b+28|0;Vc(s,r,0.0);t=b+4|0;a[t]=1;c[b>>2]=0;u=b+6|0;a[u]=1;v=b+5|0;a[v]=0;Mm(o);i=f;return}g=b+56|0;l=b+64|0;w=d+24|0;x=d+12|0;y=d+16|0;d=0;a:while(1){if((c[g>>2]|0)<=(d|0)){z=18;break}A=c[n>>2]|0;c[A+(d<<2)>>2]=-1;c[o+(d<<2)>>2]=d;if((c[l>>2]|0)<=(d|0)){z=20;break}B=c[k>>2]|0;c[B+(d<<2)>>2]=0;C=c[w>>2]|0;D=c[x>>2]|0;E=c[D+(d<<2)>>2]|0;F=c[y>>2]|0;if((F|0)==0){G=c[D+(d+1<<2)>>2]|0}else{G=(c[F+(d<<2)>>2]|0)+E|0}if((E|0)<(G|0)){F=E;do{E=c[C+(F<<2)>>2]|0;if((E|0)<(d|0)?(D=o+(E<<2)|0,(c[D>>2]|0)!=(d|0)):0){H=D;D=E;do{if(!((D|0)>-1)){z=29;break a}if((c[g>>2]|0)<=(D|0)){z=29;break a}E=A+(D<<2)|0;if((c[E>>2]|0)==-1){c[E>>2]=d}if((c[l>>2]|0)<=(D|0)){z=33;break a}I=B+(D<<2)|0;c[I>>2]=(c[I>>2]|0)+1;c[H>>2]=d;if((c[g>>2]|0)<=(D|0)){z=35;break a}D=c[E>>2]|0;H=o+(D<<2)|0}while((c[H>>2]|0)!=(d|0))}F=F+1|0}while((F|0)<(G|0))}F=d+1|0;if((F|0)<(h|0)){d=F}else{z=39;break}}if((z|0)==18){Ua(1600,1632,378,1680)}else if((z|0)==20){Ua(1600,1632,378,1680)}else if((z|0)==29){Ua(1600,1632,378,1680)}else if((z|0)==33){Ua(1600,1632,378,1680)}else if((z|0)==35){Ua(1600,1632,378,1680)}else if((z|0)==39){d=c[j>>2]|0;c[d>>2]=0;if(!m){p=d;q=p+(h<<2)|0;r=c[q>>2]|0;s=b+28|0;Vc(s,r,0.0);t=b+4|0;a[t]=1;c[b>>2]=0;u=b+6|0;a[u]=1;v=b+5|0;a[v]=0;Mm(o);i=f;return}m=b+64|0;j=e&1^1;e=0;G=0;while(1){if((c[m>>2]|0)<=(G|0)){z=42;break}g=e+j+(c[(c[k>>2]|0)+(G<<2)>>2]|0)|0;l=G+1|0;c[d+(l<<2)>>2]=g;if((l|0)<(h|0)){e=g;G=l}else{p=d;z=44;break}}if((z|0)==42){Ua(1600,1632,378,1680)}else if((z|0)==44){q=p+(h<<2)|0;r=c[q>>2]|0;s=b+28|0;Vc(s,r,0.0);t=b+4|0;a[t]=1;c[b>>2]=0;u=b+6|0;a[u]=1;v=b+5|0;a[v]=0;Mm(o);i=f;return}}}function kd(b,d){b=b|0;d=d|0;var e=0,f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0.0,P=0.0,Q=0.0,R=0.0,S=0,T=0.0;e=i;i=i+16|0;f=e;if((a[b+6|0]|0)==0){Ua(2232,2296,105,2360)}h=c[d+8>>2]|0;if((h|0)!=(c[d+4>>2]|0)){Ua(2384,2296,106,2360)}j=b+56|0;if((c[j>>2]|0)!=(h|0)){Ua(2408,2296,108,2360)}k=b+64|0;if((c[k>>2]|0)!=(h|0)){Ua(2432,2296,109,2360)}l=c[b+20>>2]|0;m=c[b+32>>2]|0;n=c[b+28>>2]|0;if(h>>>0>1073741823){Wc()}o=h<<2;if((Qm(f,16,o)|0)==0){p=c[f>>2]|0}else{c[f>>2]=0;p=0}q=(o|0)==0;if(!((p|0)!=0|q)){Wc()}if((Qm(f,16,o)|0)==0){r=c[f>>2]|0}else{c[f>>2]=0;r=0}if(!((r|0)!=0|q)){Wc()}if((Qm(f,16,o)|0)==0){s=c[f>>2]|0}else{c[f>>2]=0;s=0}if(!((s|0)!=0|q)){Wc()}q=b+48|0;if((c[q>>2]|0)!=0){f=b+44|0;Mm(c[f>>2]|0);c[f>>2]=0}c[q>>2]=0;if((h|0)<=0){t=0;c[b>>2]=t;u=b+5|0;a[u]=1;Mm(s);Mm(r);Mm(p);i=e;return}q=b+60|0;f=d+20|0;o=d+24|0;v=d+12|0;w=d+16|0;d=b+88|0;x=b+84|0;y=b+52|0;z=0;a:while(1){A=p+(z<<2)|0;g[A>>2]=0.0;c[s+(z<<2)>>2]=z;if((c[k>>2]|0)<=(z|0)){B=31;break}C=c[q>>2]|0;D=C+(z<<2)|0;c[D>>2]=0;E=c[f>>2]|0;F=c[o>>2]|0;G=c[v>>2]|0;H=c[G+(z<<2)>>2]|0;I=c[w>>2]|0;if((I|0)==0){J=c[G+(z+1<<2)>>2]|0}else{J=(c[I+(z<<2)>>2]|0)+H|0}if((H|0)<(J|0)){I=H;H=h;while(1){G=c[F+(I<<2)>>2]|0;if((G|0)<=(z|0)?(K=p+(G<<2)|0,g[K>>2]=+g[E+(I<<2)>>2]+ +g[K>>2],K=s+(G<<2)|0,(c[K>>2]|0)!=(z|0)):0){L=K;K=G;G=0;while(1){c[r+(G<<2)>>2]=K;c[L>>2]=z;if(!((K|0)>-1)){B=41;break a}if((c[j>>2]|0)<=(K|0)){B=41;break a}M=G+1|0;K=c[(c[y>>2]|0)+(K<<2)>>2]|0;L=s+(K<<2)|0;if((c[L>>2]|0)==(z|0)){break}else{G=M}}if((G|0)>-1){L=M;K=H;do{L=L+ -1|0;K=K+ -1|0;c[r+(K<<2)>>2]=c[r+(L<<2)>>2]}while((L|0)>0);N=H-M|0}else{N=H}}else{N=H}I=I+1|0;if((I|0)>=(J|0)){break}else{H=N}}O=+g[A>>2]*+g[d>>2]+ +g[x>>2];g[A>>2]=0.0;if((N|0)<(h|0)){P=O;H=N;while(1){I=c[r+(H<<2)>>2]|0;E=p+(I<<2)|0;Q=+g[E>>2];g[E>>2]=0.0;E=c[l+(I<<2)>>2]|0;R=Q/+g[n+(E<<2)>>2];if(!((I|0)>-1)){B=50;break a}if((c[k>>2]|0)<=(I|0)){B=50;break a}F=C+(I<<2)|0;L=(c[F>>2]|0)+E|0;K=E+1|0;if((K|0)<(L|0)){E=K;while(1){G=p+(c[m+(E<<2)>>2]<<2)|0;g[G>>2]=+g[G>>2]-R*+g[n+(E<<2)>>2];G=E+1|0;if((G|0)==(L|0)){S=L;break}else{E=G}}}else{S=K}c[m+(S<<2)>>2]=z;g[n+(S<<2)>>2]=R;if((c[k>>2]|0)<=(I|0)){B=54;break a}Q=P-R*R;c[F>>2]=(c[F>>2]|0)+1;E=H+1|0;if((E|0)<(h|0)){P=Q;H=E}else{T=Q;break}}}else{T=O}}else{P=+g[d>>2]*0.0+ +g[x>>2];g[A>>2]=0.0;T=P}if((c[k>>2]|0)<=(z|0)){B=57;break}H=c[l+(z<<2)>>2]|0;C=c[D>>2]|0;c[D>>2]=C+1;E=C+H|0;c[m+(E<<2)>>2]=z;if(T<=0.0){t=1;B=60;break}g[n+(E<<2)>>2]=+U(+T);E=z+1|0;if((E|0)<(h|0)){z=E}else{t=0;B=60;break}}if((B|0)==31){Ua(1600,1632,378,1680)}else if((B|0)==41){Ua(1600,1632,378,1680)}else if((B|0)==50){Ua(1600,1632,378,1680)}else if((B|0)==54){Ua(1600,1632,378,1680)}else if((B|0)==57){Ua(1600,1632,378,1680)}else if((B|0)==60){c[b>>2]=t;u=b+5|0;a[u]=1;Mm(s);Mm(r);Mm(p);i=e;return}}function ld(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,g=0,h=0;d=i;i=i+16|0;e=d;if(!((b|0)>-1)){Ua(2464,2608,268,2656)}f=a+4|0;if((c[f>>2]|0)==(b|0)){c[f>>2]=b;i=d;return}Mm(c[a>>2]|0);if((b|0)==0){c[a>>2]=0;c[f>>2]=b;i=d;return}if(b>>>0>1073741823){Wc()}g=b<<2;if((Qm(e,16,g)|0)==0){h=c[e>>2]|0}else{c[e>>2]=0;h=0}if(!((h|0)!=0|(g|0)==0)){Wc()}c[a>>2]=h;c[f>>2]=b;i=d;return}function md(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,V=0,W=0,X=0,Y=0,Z=0,_=0,$=0,aa=0,ba=0,ca=0,da=0,ea=0,fa=0,ga=0,ha=0,ia=0,ja=0,ka=0,la=0,ma=0,na=0,oa=0,pa=0,qa=0,ra=0,sa=0,ta=0,ua=0,va=0,wa=0,xa=0,ya=0,za=0,Aa=0,Ba=0,Ca=0,Da=0,Ea=0,Fa=0,Ga=0,Ha=0,Ia=0,Ja=0,Ka=0,La=0,Ma=0,Na=0,Oa=0,Pa=0,Qa=0,Ra=0,Sa=0,Ta=0,Va=0,Wa=0,Xa=0,Ya=0,Za=0,_a=0,$a=0,ab=0,bb=0,cb=0,db=0,eb=0,fb=0,gb=0,hb=0,ib=0,jb=0,kb=0,lb=0;d=i;e=c[a+4>>2]|0;f=~~(+U(+(+(e|0)))*10.0);g=(f|0)>16?f:16;f=e+ -2|0;h=(g|0)<(f|0)?g:f;f=c[a+16>>2]|0;if((f|0)!=0){if(!((e|0)>-1)){Ua(1536,1552,138,1592)}if((e|0)!=0){g=c[f>>2]|0;if((e|0)>1){j=g;k=1;while(1){l=(c[f+(k<<2)>>2]|0)+j|0;m=k+1|0;if((m|0)==(e|0)){n=l;break}else{j=l;k=m}}}else{n=g}}else{n=0}}else{n=c[a+28>>2]|0}g=e+1|0;ld(b,g);k=n+(e<<1)+((n|0)/5|0)|0;Vc(a+20|0,k,0.0);j=(e<<3)+8|0;f=Tm(j>>>0>1073741823?-1:j<<2)|0;j=g<<1;m=f+(j<<2)|0;l=g*3|0;o=f+(l<<2)|0;p=g<<2;q=g*5|0;r=g*6|0;s=f+(r<<2)|0;t=g*7|0;u=c[b>>2]|0;v=c[a+12>>2]|0;w=c[a+24>>2]|0;a=(e|0)>0;if(a){x=0;while(1){y=x+1|0;c[f+(x<<2)>>2]=(c[v+(y<<2)>>2]|0)-(c[v+(x<<2)>>2]|0);if((y|0)==(e|0)){break}else{x=y}}}c[f+(e<<2)>>2]=0;x=(e|0)<0;if(!x){y=0;do{c[f+(y+l<<2)>>2]=-1;c[u+(y<<2)>>2]=-1;c[f+(y+j<<2)>>2]=-1;c[f+(y+t<<2)>>2]=-1;c[f+(y+g<<2)>>2]=1;c[f+(y+r<<2)>>2]=1;c[f+(y+p<<2)>>2]=0;c[f+(y+q<<2)>>2]=c[f+(y<<2)>>2];y=y+1|0}while((y|0)!=(g|0))}if(a){y=0;do{z=f+(y+r<<2)|0;if((c[z>>2]|0)!=0){c[z>>2]=1}y=y+1|0}while((y|0)!=(e|0));c[f+(p+e<<2)>>2]=-2;c[v+(e<<2)>>2]=-1;c[f+(r+e<<2)>>2]=0;if(a){y=-2-e|0;z=f+(g+e<<2)|0;A=0;B=0;while(1){C=c[f+(B+q<<2)>>2]|0;do{if((C|0)==0){c[f+(B+p<<2)>>2]=-2;c[v+(B<<2)>>2]=-1;c[f+(B+r<<2)>>2]=0;D=A+1|0}else{if((C|0)>(h|0)){c[f+(B+g<<2)>>2]=0;c[f+(B+p<<2)>>2]=-1;c[v+(B<<2)>>2]=y;c[z>>2]=(c[z>>2]|0)+1;D=A+1|0;break}E=f+(C+l<<2)|0;F=c[E>>2]|0;if((F|0)==-1){G=-1}else{c[u+(F<<2)>>2]=B;G=c[E>>2]|0}c[f+(B+j<<2)>>2]=G;c[E>>2]=B;D=A}}while(0);C=B+1|0;if((C|0)==(e|0)){H=D;break}else{A=D;B=C}}}else{H=0}}else{c[f+(p+e<<2)>>2]=-2;c[v+(e<<2)>>2]=-1;c[f+(r+e<<2)>>2]=0;H=0}if((H|0)<(e|0)){B=a^1;D=0;A=0;G=n;n=2;z=H;while(1){a:do{if((D|0)<(e|0)){H=D;while(1){y=c[f+(H+l<<2)>>2]|0;h=H+1|0;if(!((y|0)==-1)){I=H;J=y;break a}if((h|0)<(e|0)){H=h}else{I=h;J=-1;break}}}else{I=D;J=-1}}while(0);H=f+(J+j<<2)|0;h=c[H>>2]|0;if((h|0)==-1){K=-1}else{c[u+(h<<2)>>2]=-1;K=c[H>>2]|0}c[f+(I+l<<2)>>2]=K;H=f+(J+p<<2)|0;h=c[H>>2]|0;y=f+(J+g<<2)|0;C=c[y>>2]|0;E=C+z|0;b:do{if((h|0)>0?(I+G|0)>=(k|0):0){if(a){F=0;while(1){L=v+(F<<2)|0;M=c[L>>2]|0;if((M|0)>-1){N=w+(M<<2)|0;c[L>>2]=c[N>>2];c[N>>2]=-2-F}N=F+1|0;if((N|0)==(e|0)){O=0;P=0;break}else{F=N}}}else{O=0;P=0}c:while(1){F=O;while(1){if((F|0)>=(G|0)){Q=P;break b}R=F+1|0;S=-2-(c[w+(F<<2)>>2]|0)|0;if((S|0)>-1){break}else{F=R}}F=v+(S<<2)|0;c[w+(P<<2)>>2]=c[F>>2];c[F>>2]=P;F=P+1|0;N=f+(S<<2)|0;if(((c[N>>2]|0)+ -1|0)>0){T=0;V=R;W=F}else{O=R;P=F;continue}while(1){F=V+1|0;c[w+(W<<2)>>2]=c[w+(V<<2)>>2];L=T+1|0;M=W+1|0;if((L|0)<((c[N>>2]|0)+ -1|0)){T=L;V=F;W=M}else{O=F;P=M;continue c}}}}else{Q=G}}while(0);c[y>>2]=0-C;N=v+(J<<2)|0;M=c[N>>2]|0;F=(h|0)==0?M:Q;L=h+1|0;X=f+(J<<2)|0;if((h|0)<0){Y=0;Z=F}else{_=-2-J|0;$=0;aa=1;ba=M;M=F;while(1){if((aa|0)>(h|0)){ca=J;da=(c[X>>2]|0)-h|0;ea=ba;fa=ba}else{ga=c[w+(ba<<2)>>2]|0;ca=ga;da=c[f+(ga<<2)>>2]|0;ea=ba+1|0;fa=c[v+(ga<<2)>>2]|0}if((da|0)<1){ha=$;ia=M}else{ga=$;ja=1;ka=fa;la=M;while(1){ma=ka+1|0;na=c[w+(ka<<2)>>2]|0;oa=f+(na+g<<2)|0;pa=c[oa>>2]|0;do{if((pa|0)>=1){qa=pa+ga|0;c[oa>>2]=0-pa;ra=la+1|0;c[w+(la<<2)>>2]=na;sa=f+(na+j<<2)|0;ta=c[sa>>2]|0;ua=u+(na<<2)|0;if((ta|0)==-1){va=-1}else{c[u+(ta<<2)>>2]=c[ua>>2];va=c[sa>>2]|0}sa=c[ua>>2]|0;if((sa|0)==-1){c[f+((c[f+(na+q<<2)>>2]|0)+l<<2)>>2]=va;wa=qa;xa=ra;break}else{c[f+(sa+j<<2)>>2]=va;wa=qa;xa=ra;break}}else{wa=ga;xa=la}}while(0);if((ja|0)==(da|0)){ha=wa;ia=xa;break}else{ga=wa;ja=ja+1|0;ka=ma;la=xa}}}if((ca|0)!=(J|0)){c[v+(ca<<2)>>2]=_;c[f+(ca+r<<2)>>2]=0}if((aa|0)<(L|0)){$=ha;aa=aa+1|0;ba=ea;M=ia}else{Y=ha;Z=ia;break}}}M=(h|0)!=0;ba=f+(J+q<<2)|0;c[ba>>2]=Y;c[N>>2]=F;c[X>>2]=Z-F;c[H>>2]=-2;if((n|0)<2){if(a){ya=0;za=65}else{Aa=2}}else{aa=(n+A|0)>-1;if(aa|B){Aa=aa?n:2}else{ya=0;za=65}}if((za|0)==65){while(1){za=0;aa=f+(ya+r<<2)|0;if((c[aa>>2]|0)!=0){c[aa>>2]=1}aa=ya+1|0;if((aa|0)==(e|0)){Aa=2;break}else{ya=aa;za=65}}}H=(F|0)<(Z|0);if(H){h=F;do{aa=c[w+(h<<2)>>2]|0;$=c[f+(aa+p<<2)>>2]|0;if(($|0)>=1?(L=c[f+(aa+g<<2)>>2]|0,_=L+Aa|0,la=v+(aa<<2)|0,aa=c[la>>2]|0,ka=$+ -1|0,(aa|0)<=(ka+aa|0)):0){$=aa;while(1){aa=c[w+($<<2)>>2]|0;ja=f+(aa+r<<2)|0;ga=c[ja>>2]|0;if((ga|0)<(Aa|0)){if((ga|0)!=0){c[ja>>2]=_+(c[f+(aa+q<<2)>>2]|0)}}else{c[ja>>2]=ga+L}if(($|0)<(ka+(c[la>>2]|0)|0)){$=$+1|0}else{break}}}h=h+1|0}while((h|0)!=(Z|0));if(H){h=-2-J|0;$=Y;la=E;ka=C;L=F;while(1){_=c[w+(L<<2)>>2]|0;ga=v+(_<<2)|0;ja=c[ga>>2]|0;aa=f+(_+p<<2)|0;na=(c[aa>>2]|0)+ja|0;pa=na+ -1|0;if((ja|0)>(pa|0)){Ba=0;Ca=0;Da=ja}else{oa=0;ra=0;qa=ja;sa=ja;while(1){ua=c[w+(qa<<2)>>2]|0;ta=f+(ua+r<<2)|0;Ea=c[ta>>2]|0;do{if((Ea|0)!=0){Fa=Ea-Aa|0;if((Fa|0)>0){c[w+(sa<<2)>>2]=ua;Ga=Fa+oa|0;Ha=ua+ra|0;Ia=sa+1|0;break}else{c[v+(ua<<2)>>2]=h;c[ta>>2]=0;Ga=oa;Ha=ra;Ia=sa;break}}else{Ga=oa;Ha=ra;Ia=sa}}while(0);if((qa|0)<(pa|0)){oa=Ga;ra=Ha;qa=qa+1|0;sa=Ia}else{Ba=Ga;Ca=Ha;Da=Ia;break}}}sa=1-ja|0;c[aa>>2]=sa+Da;qa=f+(_<<2)|0;ra=(c[qa>>2]|0)+ja|0;if((na|0)<(ra|0)){oa=Ba;pa=Ca;ta=na;ua=Da;while(1){Ea=c[w+(ta<<2)>>2]|0;ma=c[f+(Ea+g<<2)>>2]|0;if((ma|0)<1){Ja=oa;Ka=pa;La=ua}else{c[w+(ua<<2)>>2]=Ea;Ja=ma+oa|0;Ka=Ea+pa|0;La=ua+1|0}Ea=ta+1|0;if((Ea|0)==(ra|0)){Ma=Ja;Na=Ka;Oa=La;break}else{oa=Ja;pa=Ka;ta=Ea;ua=La}}}else{Ma=Ba;Na=Ca;Oa=Da}if((Ma|0)==0){c[ga>>2]=h;ua=f+(_+g<<2)|0;ta=c[ua>>2]|0;c[ua>>2]=0;c[aa>>2]=-1;Pa=ta+$|0;Qa=la-ta|0;Ra=ka-ta|0}else{ta=f+(_+q<<2)|0;ua=c[ta>>2]|0;c[ta>>2]=(Ma|0)<(ua|0)?Ma:ua;ua=w+(Da<<2)|0;c[w+(Oa<<2)>>2]=c[ua>>2];ta=w+(ja<<2)|0;c[ua>>2]=c[ta>>2];c[ta>>2]=J;c[qa>>2]=sa+Oa;ta=(Na>>>0)%(e>>>0)|0;ua=f+(ta+t<<2)|0;c[f+(_+j<<2)>>2]=c[ua>>2];c[ua>>2]=_;c[u+(_<<2)>>2]=ta;Pa=$;Qa=la;Ra=ka}ta=L+1|0;if((ta|0)==(Z|0)){Sa=Pa;Ta=Qa;Va=Ra;break}else{$=Pa;la=Qa;ka=Ra;L=ta}}}else{Sa=Y;Ta=E;Va=C}}else{Sa=Y;Ta=E;Va=C}c[ba>>2]=Sa;A=(A|0)<(Sa|0)?Sa:A;L=A+Aa|0;ka=(L+A|0)>-1&(L|0)>1;if(ka|B){Wa=ka?L:2}else{L=0;while(1){ka=f+(L+r<<2)|0;if((c[ka>>2]|0)!=0){c[ka>>2]=1}ka=L+1|0;if((ka|0)==(e|0)){Wa=2;break}else{L=ka}}}if(H){L=Wa;ba=F;while(1){C=c[w+(ba<<2)>>2]|0;d:do{if(!((c[f+(C+g<<2)>>2]|0)>-1)?(E=f+((c[u+(C<<2)>>2]|0)+t<<2)|0,ka=c[E>>2]|0,c[E>>2]=-1,!((ka|0)==-1)):0){E=ka;ka=L;while(1){la=f+(E+j<<2)|0;$=c[la>>2]|0;if(($|0)==-1){Xa=ka;break d}h=c[f+(E<<2)>>2]|0;ta=c[f+(E+p<<2)>>2]|0;ua=v+(E<<2)|0;pa=c[ua>>2]|0;oa=h+ -1|0;if((pa|0)<(oa+pa|0)){ra=pa;do{ra=ra+1|0;c[f+((c[w+(ra<<2)>>2]|0)+r<<2)>>2]=ka}while((ra|0)<(oa+(c[ua>>2]|0)|0));Ya=c[la>>2]|0}else{Ya=$}e:do{if(!((Ya|0)==-1)){ua=-2-E|0;ra=f+(E+g<<2)|0;pa=Ya;na=E;while(1){Ea=f+(na+j<<2)|0;ma=pa;f:while(1){if((c[f+(ma<<2)>>2]|0)!=(h|0)){break}Fa=v+(ma<<2)|0;Za=c[Fa>>2]|0;if((c[f+(ma+p<<2)>>2]|0)!=(ta|0)){break}_a=oa+Za|0;$a=Za;while(1){Za=$a+1|0;if(($a|0)>=(_a|0)){break}if((c[f+((c[w+(Za<<2)>>2]|0)+r<<2)>>2]|0)==(ka|0)){$a=Za}else{break f}}c[Fa>>2]=ua;$a=f+(ma+g<<2)|0;c[ra>>2]=(c[ra>>2]|0)+(c[$a>>2]|0);c[$a>>2]=0;c[f+(ma+p<<2)>>2]=-1;$a=c[f+(ma+j<<2)>>2]|0;c[Ea>>2]=$a;if(($a|0)==-1){break e}else{ma=$a}}pa=c[f+(ma+j<<2)>>2]|0;if((pa|0)==-1){break}else{na=ma}}}}while(0);oa=c[la>>2]|0;ta=ka+1|0;if((oa|0)==-1){Xa=ta;break}else{E=oa;ka=ta}}}else{Xa=L}}while(0);ba=ba+1|0;if((ba|0)==(Z|0)){break}else{L=Xa}}if(H){L=e-Ta|0;ba=I;C=F;_=F;while(1){sa=c[w+(_<<2)>>2]|0;qa=f+(sa+g<<2)|0;ja=c[qa>>2]|0;aa=0-ja|0;if((aa|0)<1){ab=ba;bb=C}else{c[qa>>2]=aa;aa=f+(sa+q<<2)|0;qa=ja+Sa+(c[aa>>2]|0)|0;ga=ja+L|0;ja=(ga|0)<(qa|0)?ga:qa;qa=f+(ja+l<<2)|0;ga=c[qa>>2]|0;if((ga|0)==-1){cb=-1}else{c[u+(ga<<2)>>2]=sa;cb=c[qa>>2]|0}c[f+(sa+j<<2)>>2]=cb;c[u+(sa<<2)>>2]=-1;c[qa>>2]=sa;c[aa>>2]=ja;c[w+(C<<2)>>2]=sa;ab=(ja|0)<(ba|0)?ja:ba;bb=C+1|0}_=_+1|0;if((_|0)==(Z|0)){break}else{ba=ab;C=bb}}c[y>>2]=Va;c[X>>2]=bb-F;if((bb|0)==(F|0)){db=ab;eb=Xa;fb=bb;za=125}else{gb=ab;hb=Xa;ib=bb}}else{jb=Xa;za=123}}else{jb=Wa;za=123}if((za|0)==123){za=0;c[y>>2]=Va;c[X>>2]=0;db=I;eb=jb;fb=F;za=125}if((za|0)==125){za=0;c[N>>2]=-1;c[f+(J+r<<2)>>2]=0;gb=db;hb=eb;ib=fb}if((Ta|0)>=(e|0)){break}else{D=gb;G=M?ib:Q;n=hb;z=Ta}}}if(a){a=0;do{Ta=v+(a<<2)|0;c[Ta>>2]=-2-(c[Ta>>2]|0);a=a+1|0}while((a|0)!=(e|0))}if(x){Wc()}pn(f+((e*12|0)+12)|0,-1,(e<<2)+4|0)|0;a=e;while(1){if((c[f+(a+g<<2)>>2]|0)<=0){Ta=v+(a<<2)|0;c[f+(a+j<<2)>>2]=c[f+((c[Ta>>2]|0)+l<<2)>>2];c[f+((c[Ta>>2]|0)+l<<2)>>2]=a}if((a|0)>0){a=a+ -1|0}else{kb=e;break}}while(1){if((c[f+(kb+g<<2)>>2]|0)>=1?(a=v+(kb<<2)|0,Ta=c[a>>2]|0,!((Ta|0)==-1)):0){c[f+(kb+j<<2)>>2]=c[f+(Ta+l<<2)>>2];c[f+((c[a>>2]|0)+l<<2)>>2]=kb}if((kb|0)>0){kb=kb+ -1|0}else{break}}if(x){Wc()}x=(o|0)==0|(m|0)==0;m=(s|0)==0;o=0;kb=0;while(1){g:do{if((c[v+(kb<<2)>>2]|0)==-1){if(x|(u|0)==0|m){lb=-1}else{c[s>>2]=kb;a=o;Ta=kb;z=0;while(1){hb=Ta;n=z;while(1){Q=f+(hb+l<<2)|0;ib=c[Q>>2]|0;if((ib|0)==-1){break}c[Q>>2]=c[f+(ib+j<<2)>>2];Q=n+1|0;c[f+(Q+r<<2)>>2]=ib;if((n|0)>-2){hb=ib;n=Q}else{lb=a;break g}}Q=n+ -1|0;ib=a+1|0;c[u+(a<<2)>>2]=hb;if((n|0)<=0){lb=ib;break g}a=ib;Ta=c[f+(Q+r<<2)>>2]|0;z=Q}}}else{lb=o}}while(0);kb=kb+1|0;if((kb|0)==(g|0)){break}else{o=lb}}lb=b+4|0;o=c[lb>>2]|0;if(e>>>0>1073741823){Wc()}if(o>>>0>1073741823){Wc()}g=e<<2;kb=pd(u,g,o<<2)|0;if(!((kb|0)!=0|(g|0)==0)){Wc()}c[b>>2]=kb;c[lb>>2]=e;if((f|0)==0){i=d;return}Vm(f);i=d;return}function nd(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0;e=i;i=i+16|0;f=e+8|0;h=e;j=c[a+8>>2]|0;if(j>>>0>1073741823){Wc()}k=j<<2;if((Qm(f,16,k)|0)==0){l=c[f>>2]|0}else{c[f>>2]=0;l=0}if(!((l|0)!=0|(k|0)==0)){Wc()}c[h>>2]=l;l=h+4|0;c[l>>2]=j;if(!((j|0)>-1)){Ua(1696,1712,225,1752)}od(h,j,1);if((c[l>>2]|0)!=(j|0)){Ua(1760,1816,498,1856)}l=(j|0)>0;if(l){pn(c[h>>2]|0,0,k|0)|0}c[b+8>>2]=j;c[b+28>>2]=0;f=b+4|0;m=c[f>>2]|0;do{if((m|0)!=(j|0)|(m|0)==0){n=b+12|0;Mm(c[n>>2]|0);o=Lm(k+4|0)|0;c[n>>2]=o;if((o|0)==0){Wc()}else{c[f>>2]=j;break}}}while(0);k=b+16|0;m=c[k>>2]|0;if((m|0)==0){p=j}else{Mm(m);c[k>>2]=0;p=c[f>>2]|0}f=b+12|0;pn(c[f>>2]|0,0,(p<<2)+4|0)|0;a:do{if(l){p=(d|0)!=0;k=c[a+24>>2]|0;m=c[a+12>>2]|0;o=c[a+16>>2]|0;n=(o|0)==0;q=c[h>>2]|0;r=0;b:do{if(p){s=c[d+(r<<2)>>2]|0}else{s=r}t=c[m+(r<<2)>>2]|0;if(n){u=c[m+(r+1<<2)>>2]|0}else{u=(c[o+(r<<2)>>2]|0)+t|0}if((t|0)<(u|0)){v=t;do{t=c[k+(v<<2)>>2]|0;if((t|0)>=(r|0)){if(p){w=c[d+(t<<2)>>2]|0}else{w=t}t=(w|0)<(s|0)?s:w;if(!((t|0)>-1&(j|0)>(t|0))){x=33;break b}y=q+(t<<2)|0;c[y>>2]=(c[y>>2]|0)+1}v=v+1|0}while((v|0)<(u|0))}r=r+1|0}while((r|0)<(j|0));if((x|0)==33){Ua(1600,1632,378,1680)}r=c[f>>2]|0;c[r>>2]=0;if(!l){Vc(b+20|0,c[r+(j<<2)>>2]|0,0.0);break}q=c[h>>2]|0;p=0;k=0;while(1){if((j|0)<=(k|0)){break}o=(c[q+(k<<2)>>2]|0)+p|0;m=k+1|0;c[r+(m<<2)>>2]=o;if((m|0)<(j|0)){p=o;k=m}else{z=r;x=42;break a}}Ua(1600,1632,378,1680)}else{r=c[f>>2]|0;c[r>>2]=0;z=r;x=42}}while(0);if((x|0)==42?(Vc(b+20|0,c[z+(j<<2)>>2]|0,0.0),l):0){z=c[f>>2]|0;f=c[h>>2]|0;u=0;do{if((j|0)<=(u|0)){x=49;break}c[f+(u<<2)>>2]=c[z+(u<<2)>>2];u=u+1|0}while((u|0)<(j|0));if((x|0)==49){Ua(1600,1632,378,1680)}if(l){l=c[a+20>>2]|0;u=c[a+24>>2]|0;z=c[a+12>>2]|0;f=c[a+16>>2]|0;a=(f|0)==0;w=(d|0)==0;s=c[h>>2]|0;r=b+24|0;k=b+20|0;b=0;c:while(1){p=c[z+(b<<2)>>2]|0;if(a){A=c[z+(b+1<<2)>>2]|0}else{A=(c[f+(b<<2)>>2]|0)+p|0}if((p|0)<(A|0)){q=d+(b<<2)|0;if(w){m=p;do{o=c[u+(m<<2)>>2]|0;if((o|0)>=(b|0)){if(!((o|0)>-1&(j|0)>(o|0))){x=62;break c}n=s+(o<<2)|0;v=c[n>>2]|0;c[n>>2]=v+1;c[(c[r>>2]|0)+(v<<2)>>2]=(b|0)<(o|0)?b:o;g[(c[k>>2]|0)+(v<<2)>>2]=+g[l+(m<<2)>>2]}m=m+1|0}while((m|0)<(A|0))}else{m=p;do{v=c[u+(m<<2)>>2]|0;if((v|0)>=(b|0)){o=c[q>>2]|0;n=c[d+(v<<2)>>2]|0;v=(n|0)<(o|0)?o:n;if(!((v|0)>-1&(j|0)>(v|0))){x=62;break c}y=s+(v<<2)|0;v=c[y>>2]|0;c[y>>2]=v+1;c[(c[r>>2]|0)+(v<<2)>>2]=(o|0)<(n|0)?o:n;g[(c[k>>2]|0)+(v<<2)>>2]=+g[l+(m<<2)>>2]}m=m+1|0}while((m|0)<(A|0))}}m=b+1|0;if((m|0)<(j|0)){b=m}else{B=s;x=66;break}}if((x|0)==62){Ua(1600,1632,378,1680)}else if((x|0)==66){Mm(B);i=e;return}}}B=c[h>>2]|0;Mm(B);i=e;return}function od(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,g=0,h=0;e=i;i=i+16|0;f=e;if(!((d|0)==1&(b|0)>-1&(d|0)>-1)){Ua(2704,2608,241,2656)}d=a+4|0;if((c[d>>2]|0)==(b|0)){c[d>>2]=b;i=e;return}Mm(c[a>>2]|0);if((b|0)==0){c[a>>2]=0;c[d>>2]=b;i=e;return}if(b>>>0>1073741823){Wc()}g=b<<2;if((Qm(f,16,g)|0)==0){h=c[f>>2]|0}else{c[f>>2]=0;h=0}if(!((h|0)!=0|(g|0)==0)){Wc()}c[a>>2]=h;c[d>>2]=b;i=e;return}function pd(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,g=0,h=0;e=i;i=i+16|0;f=e;if((a|0)==0){if((Qm(f,16,b)|0)==0){g=c[f>>2]|0}else{c[f>>2]=0;g=0}if((g|0)!=0|(b|0)==0){h=g;i=e;return h|0}else{Wc()}}if((b|0)==0){Mm(a);h=0;i=e;return h|0}if((Qm(f,16,b)|0)!=0){c[f>>2]=0;Wc()}g=c[f>>2]|0;if((g|0)==0){Wc()}nn(g|0,a|0,(d>>>0<b>>>0?d:b)|0)|0;Mm(a);h=g;i=e;return h|0}function qd(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0;e=i;i=i+16|0;f=e+8|0;h=e;j=c[a+8>>2]|0;c[h>>2]=0;k=h+4|0;c[k>>2]=0;if(!((j|0)>-1)){Ua(2464,2608,268,2656)}l=(j|0)==0;do{if(!l){if(j>>>0>1073741823){Wc()}m=j<<2;if((Qm(f,16,m)|0)==0){n=c[f>>2]|0}else{c[f>>2]=0;n=0}if((n|0)!=0|(m|0)==0){c[h>>2]=n;break}else{Wc()}}}while(0);c[k>>2]=j;od(h,j,1);if((c[k>>2]|0)!=(j|0)){Ua(1760,1816,498,1856)}k=(j|0)>0;if(k){pn(c[h>>2]|0,0,j<<2|0)|0}c[b+8>>2]=j;c[b+28>>2]=0;n=b+4|0;f=c[n>>2]|0;do{if((f|0)!=(j|0)|(f|0)==0){m=b+12|0;Mm(c[m>>2]|0);o=Lm((j<<2)+4|0)|0;c[m>>2]=o;if((o|0)==0){Wc()}else{c[n>>2]=j;break}}}while(0);f=b+16|0;o=c[f>>2]|0;if((o|0)==0){p=j}else{Mm(o);c[f>>2]=0;p=c[n>>2]|0}n=b+12|0;pn(c[n>>2]|0,0,(p<<2)+4|0)|0;a:do{if(k){p=(d|0)!=0;f=c[a+24>>2]|0;o=c[a+12>>2]|0;m=c[a+16>>2]|0;q=(m|0)==0;r=c[h>>2]|0;s=0;b:while(1){if(p){t=c[d+(s<<2)>>2]|0}else{t=s}u=c[o+(s<<2)>>2]|0;if(q){v=c[o+(s+1<<2)>>2]|0}else{v=(c[m+(s<<2)>>2]|0)+u|0}if((u|0)<(v|0)){w=r+(t<<2)|0;x=u;do{u=c[f+(x<<2)>>2]|0;if(p){y=c[d+(u<<2)>>2]|0}else{y=u}if((u|0)!=(s|0)){if((u|0)>(s|0)){if(!((y|0)>-1&(j|0)>(y|0))){z=39;break b}u=r+(y<<2)|0;c[u>>2]=(c[u>>2]|0)+1;if((t|0)<0|(j|0)<=(t|0)){z=41;break b}c[w>>2]=(c[w>>2]|0)+1}}else{if(!((y|0)>-1&(j|0)>(y|0))){z=35;break b}u=r+(y<<2)|0;c[u>>2]=(c[u>>2]|0)+1}x=x+1|0}while((x|0)<(v|0))}s=s+1|0;if((s|0)>=(j|0)){break a}}if((z|0)==35){Ua(1600,1632,378,1680)}else if((z|0)==39){Ua(1600,1632,378,1680)}else if((z|0)==41){Ua(1600,1632,378,1680)}}}while(0);if(!l){if(!k){Ua(1432,1496,175,1528)}l=c[h>>2]|0;v=c[l>>2]|0;if((j|0)>1){y=v;t=1;while(1){s=(c[l+(t<<2)>>2]|0)+y|0;r=t+1|0;if((r|0)==(j|0)){A=s;break}else{y=s;t=r}}}else{A=v}}else{A=0}Vc(b+20|0,A,0.0);A=c[n>>2]|0;c[A>>2]=0;if(k){n=c[h>>2]|0;v=0;t=0;do{if((j|0)<=(t|0)){z=55;break}v=(c[n+(t<<2)>>2]|0)+v|0;t=t+1|0;c[A+(t<<2)>>2]=v}while((t|0)<(j|0));if((z|0)==55){Ua(1600,1632,378,1680)}if(k){t=c[h>>2]|0;v=0;do{if((j|0)<=(v|0)){z=61;break}c[t+(v<<2)>>2]=c[A+(v<<2)>>2];v=v+1|0}while((v|0)<(j|0));if((z|0)==61){Ua(1600,1632,378,1680)}if(k){k=c[a+20>>2]|0;v=c[a+24>>2]|0;A=c[a+12>>2]|0;t=c[a+16>>2]|0;a=(t|0)==0;n=(d|0)==0;y=c[h>>2]|0;l=b+24|0;r=b+20|0;b=0;c:while(1){s=c[A+(b<<2)>>2]|0;if(a){B=c[A+(b+1<<2)>>2]|0}else{B=(c[t+(b<<2)>>2]|0)+s|0}if((s|0)<(B|0)){p=d+(b<<2)|0;f=s;do{s=c[v+(f<<2)>>2]|0;if(n){C=s;D=b}else{C=c[d+(s<<2)>>2]|0;D=c[p>>2]|0}if((s|0)!=(b|0)){if((s|0)>(b|0)){if(!((D|0)>-1&(j|0)>(D|0))){z=76;break c}s=y+(D<<2)|0;m=c[s>>2]|0;c[s>>2]=m+1;s=c[l>>2]|0;c[s+(m<<2)>>2]=C;o=k+(f<<2)|0;q=c[r>>2]|0;g[q+(m<<2)>>2]=+g[o>>2];if(!((C|0)>-1&(j|0)>(C|0))){z=78;break c}m=y+(C<<2)|0;x=c[m>>2]|0;c[m>>2]=x+1;c[s+(x<<2)>>2]=D;g[q+(x<<2)>>2]=+g[o>>2]}}else{if(!((C|0)>-1&(j|0)>(C|0))){z=72;break c}o=y+(C<<2)|0;x=c[o>>2]|0;c[o>>2]=x+1;c[(c[l>>2]|0)+(x<<2)>>2]=C;g[(c[r>>2]|0)+(x<<2)>>2]=+g[k+(f<<2)>>2]}f=f+1|0}while((f|0)<(B|0))}f=b+1|0;if((f|0)<(j|0)){b=f}else{E=y;z=82;break}}if((z|0)==72){Ua(1600,1632,378,1680)}else if((z|0)==76){Ua(1600,1632,378,1680)}else if((z|0)==78){Ua(1600,1632,378,1680)}else if((z|0)==82){Mm(E);i=e;return}}}}E=c[h>>2]|0;Mm(E);i=e;return}function rd(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0.0;f=i;i=i+64|0;h=f+48|0;j=f;k=f+40|0;l=e+8|0;m=c[l>>2]|0;n=e+4|0;o=c[n>>2]|0;a[j]=0;p=j+4|0;c[p+0>>2]=0;c[p+4>>2]=0;c[p+8>>2]=0;c[p+12>>2]=0;c[p+16>>2]=0;c[p+20>>2]=0;c[p+24>>2]=0;c[p+28>>2]=0;c[j+8>>2]=o;c[j+28>>2]=0;o=Lm((m<<2)+4|0)|0;c[j+12>>2]=o;if((o|0)==0){Wc()}c[p>>2]=m;p=j+16|0;q=j+12|0;pn(o|0,0,(m<<2)+4|0)|0;if((c[b>>2]|0)!=(c[d>>2]|0)){if(m>>>0>1073741823){Wc()}o=m<<2;if((Qm(h,16,o)|0)==0){r=c[h>>2]|0}else{c[h>>2]=0;r=0}if(!((r|0)!=0|(o|0)==0)){Wc()}c[k>>2]=r;r=k+4|0;c[r>>2]=m;if(!((m|0)>-1)){Ua(1696,1712,225,1752)}od(k,m,1);if((c[r>>2]|0)!=(m|0)){Ua(1760,1816,498,1856)}if((m|0)>0){pn(c[k>>2]|0,0,o|0)|0}o=c[b>>2]|0;r=c[d>>2]|0;a:do{if((o|0)!=(r|0)){h=c[k>>2]|0;s=o;while(1){t=c[s>>2]|0;if(!((t|0)>-1)){u=23;break}if((t|0)>=(c[l>>2]|0)){u=23;break}v=c[s+4>>2]|0;if(!((v|0)>-1)){u=23;break}if((v|0)>=(c[n>>2]|0)){u=23;break}if((m|0)<=(t|0)){u=25;break}v=h+(t<<2)|0;c[v>>2]=(c[v>>2]|0)+1;s=s+12|0;if((s|0)==(r|0)){break a}}if((u|0)==23){Ua(3272,704,952,3352)}else if((u|0)==25){Ua(1600,1632,394,3504)}}}while(0);td(j,k);u=c[b>>2]|0;b=c[d>>2]|0;b:do{if((u|0)!=(b|0)){d=c[p>>2]|0;r=c[q>>2]|0;m=c[j+24>>2]|0;n=c[j+20>>2]|0;if((d|0)==0){Ua(3376,704,842,3408)}else{w=u}while(1){l=c[w>>2]|0;o=d+(l<<2)|0;s=c[o>>2]|0;h=c[r+(l<<2)>>2]|0;if((s|0)>((c[r+(l+1<<2)>>2]|0)-h|0)){break}l=c[w+4>>2]|0;x=+g[w+8>>2];c[o>>2]=s+1;o=h+s|0;c[m+(o<<2)>>2]=l;g[n+(o<<2)>>2]=x;w=w+12|0;if((w|0)==(b|0)){break b}}Ua(3432,704,843,3408)}}while(0);sd(j);Mm(c[k>>2]|0)}Xc(e,j)|0;Mm(c[q>>2]|0);Mm(c[p>>2]|0);p=c[j+20>>2]|0;if((p|0)!=0){Vm(p)}p=c[j+24>>2]|0;if((p|0)==0){i=f;return}Vm(p);i=f;return}function sd(a){a=a|0;var b=0,d=0,e=0,f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0.0,L=0,M=0;b=i;i=i+16|0;d=b+8|0;e=b;f=a+16|0;if((c[f>>2]|0)==0){Ua(3376,704,1020,3392)}h=c[a+8>>2]|0;if(h>>>0>1073741823){Wc()}j=h<<2;if((Qm(d,16,j)|0)==0){k=c[d>>2]|0}else{c[d>>2]=0;k=0}if(!((k|0)!=0|(j|0)==0)){Wc()}c[e>>2]=k;k=e+4|0;c[k>>2]=h;if(!((h|0)>-1)){Ua(1696,1712,225,1752)}od(e,h,1);if((c[k>>2]|0)!=(h|0)){Ua(1760,1816,498,1856)}if((h|0)>0){pn(c[e>>2]|0,-1,h<<2|0)|0}k=a+4|0;j=c[k>>2]|0;d=a+12|0;l=c[d>>2]|0;m=c[f>>2]|0;if((j|0)<=0){n=j;o=0;p=l+(n<<2)|0;c[p>>2]=o;Mm(m);c[f>>2]=0;q=a+20|0;r=c[k>>2]|0;s=c[d>>2]|0;t=s+(r<<2)|0;u=c[t>>2]|0;Vc(q,u,0.0);v=c[e>>2]|0;Mm(v);i=b;return}j=a+24|0;w=c[e>>2]|0;x=a+20|0;y=0;z=0;a:while(1){A=l+(z<<2)|0;B=c[A>>2]|0;C=c[m+(z<<2)>>2]|0;D=C+B|0;if((C|0)>0){C=c[j>>2]|0;E=y;F=B;while(1){B=c[C+(F<<2)>>2]|0;if(!((B|0)>-1&(h|0)>(B|0))){G=21;break a}H=w+(B<<2)|0;I=c[H>>2]|0;J=c[x>>2]|0;K=+g[J+(F<<2)>>2];if((I|0)<(y|0)){g[J+(E<<2)>>2]=K;c[C+(E<<2)>>2]=B;c[H>>2]=E;L=E+1|0}else{H=J+(I<<2)|0;g[H>>2]=K+ +g[H>>2];L=E}H=F+1|0;if((H|0)<(D|0)){E=L;F=H}else{M=L;break}}}else{M=y}c[A>>2]=y;F=z+1|0;E=c[k>>2]|0;if((F|0)<(E|0)){y=M;z=F}else{n=E;o=M;G=27;break}}if((G|0)==21){Ua(1600,1632,394,3504)}else if((G|0)==27){p=l+(n<<2)|0;c[p>>2]=o;Mm(m);c[f>>2]=0;q=a+20|0;r=c[k>>2]|0;s=c[d>>2]|0;t=s+(r<<2)|0;u=c[t>>2]|0;Vc(q,u,0.0);v=c[e>>2]|0;Mm(v);i=b;return}}function td(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0;d=i;e=a+16|0;f=c[e>>2]|0;h=a+4|0;j=c[h>>2]|0;k=j<<2;if((f|0)!=0){l=Lm(k+4|0)|0;if((l|0)==0){Wc()}a:do{if((j|0)>0){m=c[a+12>>2]|0;n=b+4|0;o=0;p=0;while(1){c[l+(p<<2)>>2]=o;q=p+1|0;r=c[f+(p<<2)>>2]|0;s=(c[m+(q<<2)>>2]|0)-(c[m+(p<<2)>>2]|0)-r|0;if((c[n>>2]|0)<=(p|0)){break}t=c[(c[b>>2]|0)+(p<<2)>>2]|0;u=r+o+((t|0)<(s|0)?s:t)|0;if((q|0)<(j|0)){o=u;p=q}else{v=u;break a}}Ua(1600,1632,156,1680)}else{v=0}}while(0);c[l+(j<<2)>>2]=v;f=a+20|0;Vc(f,v,0.0);v=c[h>>2]|0;p=c[a+12>>2]|0;if((v|0)>0){o=a+24|0;n=v;do{n=n+ -1|0;v=l+(n<<2)|0;m=c[v>>2]|0;u=p+(n<<2)|0;q=c[u>>2]|0;if((m-q|0)>0?(t=c[(c[e>>2]|0)+(n<<2)>>2]|0,(t|0)>0):0){s=c[o>>2]|0;r=c[f>>2]|0;w=q;q=m;m=t;do{m=m+ -1|0;c[s+(q+m<<2)>>2]=c[s+(w+m<<2)>>2];w=c[u>>2]|0;q=c[v>>2]|0;g[r+(q+m<<2)>>2]=+g[r+(w+m<<2)>>2]}while((m|0)>0)}}while((n|0)>0)}c[a+12>>2]=l;Mm(p);i=d;return}p=Lm(k)|0;c[e>>2]=p;if((p|0)==0){Wc()}b:do{if((j|0)>0){k=b+4|0;l=a+12|0;n=0;f=0;o=0;while(1){c[p+(f<<2)>>2]=n;if((c[k>>2]|0)<=(f|0)){break}m=c[(c[b>>2]|0)+(f<<2)>>2]|0;w=f+1|0;r=c[l>>2]|0;q=m+o|0;if((w|0)<(j|0)){n=m+n-(c[r+(f<<2)>>2]|0)+(c[r+(w<<2)>>2]|0)|0;f=w;o=q}else{x=q;break b}}Ua(1600,1632,156,1680)}else{x=0}}while(0);o=a+20|0;f=c[a+28>>2]|0;n=f+x|0;x=a+32|0;if(n>>>0>(c[x>>2]|0)>>>0){l=n>>>0>1073741823?-1:n<<2;k=Tm(l)|0;q=Tm(l)|0;l=c[o>>2]|0;w=(f>>>0<n>>>0?f:n)<<2;nn(k|0,l|0,w|0)|0;f=a+24|0;r=c[f>>2]|0;nn(q|0,r|0,w|0)|0;if((l|0)==0){y=r}else{Vm(l);y=c[f>>2]|0}if((y|0)!=0){Vm(y)}c[o>>2]=k;c[f>>2]=q;c[x>>2]=n;z=c[h>>2]|0}else{z=j}j=c[a+12>>2]|0;if((z|0)>0){n=c[e>>2]|0;x=a+24|0;a=z;q=c[j+(z<<2)>>2]|0;while(1){a=a+ -1|0;f=j+(a<<2)|0;k=c[f>>2]|0;y=q-k|0;if((y|0)>0){l=c[x>>2]|0;r=p+(a<<2)|0;w=c[o>>2]|0;m=k;v=c[r>>2]|0;u=y;while(1){s=u+ -1|0;c[l+(v+s<<2)>>2]=c[l+(m+s<<2)>>2];t=c[f>>2]|0;A=c[r>>2]|0;g[w+(A+s<<2)>>2]=+g[w+(t+s<<2)>>2];if((s|0)>0){m=t;v=A;u=s}else{B=t;C=A;break}}}else{B=k;C=c[p+(a<<2)>>2]|0}c[f>>2]=C;c[n+(a<<2)>>2]=y;if((a|0)<=0){break}else{q=B}}D=c[h>>2]|0;E=n}else{D=z;E=c[e>>2]|0}e=D+ -1|0;z=(c[E+(e<<2)>>2]|0)+(c[j+(e<<2)>>2]|0)|0;if((D|0)<=0){Ua(1600,1632,156,1680)}if((c[b+4>>2]|0)<=(e|0)){Ua(1600,1632,156,1680)}c[j+(D<<2)>>2]=z+(c[(c[b>>2]|0)+(e<<2)>>2]|0);Vc(o,c[j+(c[h>>2]<<2)>>2]|0,0.0);i=d;return}function ud(a){a=a|0;var b=0,d=0;b=i;Mm(c[a+76>>2]|0);Mm(c[a+68>>2]|0);Mm(c[a+60>>2]|0);Mm(c[a+52>>2]|0);Mm(c[a+44>>2]|0);Mm(c[a+20>>2]|0);Mm(c[a+24>>2]|0);d=c[a+28>>2]|0;if((d|0)!=0){Vm(d)}d=c[a+32>>2]|0;if((d|0)==0){i=b;return}Vm(d);i=b;return}function vd(b){b=b|0;var d=0,e=0,f=0,g=0,h=0;b=i;i=i+16|0;d=b;e=c[r>>2]|0;Rd(4232,e,4288);c[880]=5660;c[3528>>2]=5680;c[3524>>2]=0;Ue(3528|0,4232);c[3600>>2]=0;c[3604>>2]=-1;f=c[s>>2]|0;c[1084]=5528;Yj(4340|0);c[4344>>2]=0;c[4348>>2]=0;c[4352>>2]=0;c[4356>>2]=0;c[4360>>2]=0;c[4364>>2]=0;c[1084]=4848;c[4368>>2]=f;Zj(d,4340|0);g=$j(d,8320)|0;_j(d);c[4372>>2]=g;c[4376>>2]=4296;a[4380|0]=(mc[c[(c[g>>2]|0)+28>>2]&63](g)|0)&1;c[902]=5740;c[3612>>2]=5760;Ue(3612|0,4336);c[3684>>2]=0;c[3688>>2]=-1;g=c[q>>2]|0;c[1096]=5528;Yj(4388|0);c[4392>>2]=0;c[4396>>2]=0;c[4400>>2]=0;c[4404>>2]=0;c[4408>>2]=0;c[4412>>2]=0;c[1096]=4848;c[4416>>2]=g;Zj(d,4388|0);h=$j(d,8320)|0;_j(d);c[4420>>2]=h;c[4424>>2]=4304;a[4428|0]=(mc[c[(c[h>>2]|0)+28>>2]&63](h)|0)&1;c[924]=5740;c[3700>>2]=5760;Ue(3700|0,4384);c[3772>>2]=0;c[3776>>2]=-1;h=c[(c[(c[924]|0)+ -12>>2]|0)+3720>>2]|0;c[946]=5740;c[3788>>2]=5760;Ue(3788|0,h);c[3860>>2]=0;c[3864>>2]=-1;c[(c[(c[880]|0)+ -12>>2]|0)+3592>>2]=3608;h=(c[(c[924]|0)+ -12>>2]|0)+3700|0;c[h>>2]=c[h>>2]|8192;c[(c[(c[924]|0)+ -12>>2]|0)+3768>>2]=3608;Dd(4432,e,4312|0);c[968]=5700;c[3880>>2]=5720;c[3876>>2]=0;Ue(3880|0,4432);c[3952>>2]=0;c[3956>>2]=-1;c[1122]=5592;Yj(4492|0);c[4496>>2]=0;c[4500>>2]=0;c[4504>>2]=0;c[4508>>2]=0;c[4512>>2]=0;c[4516>>2]=0;c[1122]=4592;c[4520>>2]=f;Zj(d,4492|0);f=$j(d,8328)|0;_j(d);c[4524>>2]=f;c[4528>>2]=4320;a[4532|0]=(mc[c[(c[f>>2]|0)+28>>2]&63](f)|0)&1;c[990]=5780;c[3964>>2]=5800;Ue(3964|0,4488);c[4036>>2]=0;c[4040>>2]=-1;c[1134]=5592;Yj(4540|0);c[4544>>2]=0;c[4548>>2]=0;c[4552>>2]=0;c[4556>>2]=0;c[4560>>2]=0;c[4564>>2]=0;c[1134]=4592;c[4568>>2]=g;Zj(d,4540|0);g=$j(d,8328)|0;_j(d);c[4572>>2]=g;c[4576>>2]=4328;a[4580|0]=(mc[c[(c[g>>2]|0)+28>>2]&63](g)|0)&1;c[1012]=5780;c[4052>>2]=5800;Ue(4052|0,4536);c[4124>>2]=0;c[4128>>2]=-1;g=c[(c[(c[1012]|0)+ -12>>2]|0)+4072>>2]|0;c[1034]=5780;c[4140>>2]=5800;Ue(4140|0,g);c[4212>>2]=0;c[4216>>2]=-1;c[(c[(c[968]|0)+ -12>>2]|0)+3944>>2]=3960;g=(c[(c[1012]|0)+ -12>>2]|0)+4052|0;c[g>>2]=c[g>>2]|8192;c[(c[(c[1012]|0)+ -12>>2]|0)+4120>>2]=3960;i=b;return}function wd(a){a=a|0;a=i;Af(3608)|0;Af(3784)|0;Ff(3960)|0;Ff(4136)|0;i=a;return}function xd(a){a=a|0;var b=0;b=i;c[a>>2]=5592;_j(a+4|0);i=b;return}function yd(a){a=a|0;var b=0;b=i;c[a>>2]=5592;_j(a+4|0);Um(a);i=b;return}function zd(b,d){b=b|0;d=d|0;var e=0,f=0;e=i;mc[c[(c[b>>2]|0)+24>>2]&63](b)|0;f=$j(d,8328)|0;c[b+36>>2]=f;a[b+44|0]=(mc[c[(c[f>>2]|0)+28>>2]&63](f)|0)&1;i=e;return}function Ad(a){a=a|0;var b=0,d=0,e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0;b=i;i=i+16|0;d=b+8|0;e=b;f=a+36|0;g=a+40|0;h=d+8|0;j=d;k=a+32|0;while(1){a=c[f>>2]|0;l=wc[c[(c[a>>2]|0)+20>>2]&15](a,c[g>>2]|0,d,h,e)|0;a=(c[e>>2]|0)-j|0;if((Wb(d|0,1,a|0,c[k>>2]|0)|0)!=(a|0)){m=-1;n=5;break}if((l|0)==2){m=-1;n=5;break}else if((l|0)!=1){n=4;break}}if((n|0)==4){m=((_a(c[k>>2]|0)|0)!=0)<<31>>31;i=b;return m|0}else if((n|0)==5){i=b;return m|0}return 0}function Bd(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0,k=0;f=i;a:do{if((a[b+44|0]|0)==0){if((e|0)>0){g=d;h=0;while(1){if((vc[c[(c[b>>2]|0)+52>>2]&15](b,c[g>>2]|0)|0)==-1){j=h;break a}k=h+1|0;if((k|0)<(e|0)){g=g+4|0;h=k}else{j=k;break}}}else{j=0}}else{j=Wb(d|0,4,e|0,c[b+32>>2]|0)|0}}while(0);i=f;return j|0}function Cd(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0;e=i;i=i+32|0;f=e+16|0;g=e;h=e+4|0;j=e+8|0;k=(d|0)==-1;a:do{if(!k){c[g>>2]=d;if((a[b+44|0]|0)!=0){if((Wb(g|0,4,1,c[b+32>>2]|0)|0)==1){break}else{l=-1}i=e;return l|0}c[h>>2]=f;m=g+4|0;n=b+36|0;o=b+40|0;p=f+8|0;q=f;r=b+32|0;s=g;while(1){t=c[n>>2]|0;u=rc[c[(c[t>>2]|0)+12>>2]&15](t,c[o>>2]|0,s,m,j,f,p,h)|0;if((c[j>>2]|0)==(s|0)){l=-1;v=12;break}if((u|0)==3){v=7;break}t=(u|0)==1;if(!(u>>>0<2)){l=-1;v=12;break}u=(c[h>>2]|0)-q|0;if((Wb(f|0,1,u|0,c[r>>2]|0)|0)!=(u|0)){l=-1;v=12;break}if(t){s=t?c[j>>2]|0:s}else{break a}}if((v|0)==7){if((Wb(s|0,1,1,c[r>>2]|0)|0)==1){break}else{l=-1}i=e;return l|0}else if((v|0)==12){i=e;return l|0}}}while(0);l=k?0:d;i=e;return l|0}function Dd(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0;f=i;i=i+16|0;g=f;c[b>>2]=5592;h=b+4|0;Yj(h);j=b+8|0;c[j+0>>2]=0;c[j+4>>2]=0;c[j+8>>2]=0;c[j+12>>2]=0;c[j+16>>2]=0;c[j+20>>2]=0;c[b>>2]=4704;c[b+32>>2]=d;c[b+40>>2]=e;c[b+48>>2]=-1;a[b+52|0]=0;Zj(g,h);h=$j(g,8328)|0;e=b+36|0;c[e>>2]=h;d=b+44|0;c[d>>2]=mc[c[(c[h>>2]|0)+24>>2]&63](h)|0;h=c[e>>2]|0;a[b+53|0]=(mc[c[(c[h>>2]|0)+28>>2]&63](h)|0)&1;if((c[d>>2]|0)>8){jj(4800)}else{_j(g);i=f;return}}function Ed(a){a=a|0;var b=0;b=i;c[a>>2]=5592;_j(a+4|0);i=b;return}function Fd(a){a=a|0;var b=0;b=i;c[a>>2]=5592;_j(a+4|0);Um(a);i=b;return}function Gd(b,d){b=b|0;d=d|0;var e=0,f=0,g=0;e=i;f=$j(d,8328)|0;d=b+36|0;c[d>>2]=f;g=b+44|0;c[g>>2]=mc[c[(c[f>>2]|0)+24>>2]&63](f)|0;f=c[d>>2]|0;a[b+53|0]=(mc[c[(c[f>>2]|0)+28>>2]&63](f)|0)&1;if((c[g>>2]|0)>8){jj(4800)}else{i=e;return}}function Hd(a){a=a|0;var b=0,c=0;b=i;c=Kd(a,0)|0;i=b;return c|0}function Id(a){a=a|0;var b=0,c=0;b=i;c=Kd(a,1)|0;i=b;return c|0}function Jd(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;e=i;i=i+32|0;f=e+16|0;g=e;h=e+4|0;j=e+8|0;k=b+52|0;l=(a[k]|0)!=0;if((d|0)==-1){if(l){m=-1;i=e;return m|0}n=c[b+48>>2]|0;a[k]=(n|0)!=-1|0;m=n;i=e;return m|0}n=b+48|0;a:do{if(l){c[h>>2]=c[n>>2];o=c[b+36>>2]|0;p=rc[c[(c[o>>2]|0)+12>>2]&15](o,c[b+40>>2]|0,h,h+4|0,j,f,f+8|0,g)|0;if((p|0)==1|(p|0)==2){m=-1;i=e;return m|0}else if((p|0)==3){a[f]=c[n>>2];c[g>>2]=f+1}p=b+32|0;while(1){o=c[g>>2]|0;if(!(o>>>0>f>>>0)){break a}q=o+ -1|0;c[g>>2]=q;if((Va(a[q]|0,c[p>>2]|0)|0)==-1){m=-1;break}}i=e;return m|0}}while(0);c[n>>2]=d;a[k]=1;m=d;i=e;return m|0}function Kd(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0;e=i;i=i+32|0;f=e+16|0;g=e;h=e+4|0;j=e+8|0;k=b+52|0;if((a[k]|0)!=0){l=b+48|0;m=c[l>>2]|0;if(!d){n=m;i=e;return n|0}c[l>>2]=-1;a[k]=0;n=m;i=e;return n|0}m=c[b+44>>2]|0;k=(m|0)>1?m:1;a:do{if((k|0)>0){m=b+32|0;l=0;while(1){o=Nb(c[m>>2]|0)|0;if((o|0)==-1){n=-1;break}a[f+l|0]=o;l=l+1|0;if((l|0)>=(k|0)){break a}}i=e;return n|0}}while(0);b:do{if((a[b+53|0]|0)==0){l=b+40|0;m=b+36|0;o=g+4|0;p=b+32|0;q=k;while(1){r=c[l>>2]|0;s=r;t=c[s>>2]|0;u=c[s+4>>2]|0;s=c[m>>2]|0;v=f+q|0;w=rc[c[(c[s>>2]|0)+16>>2]&15](s,r,f,v,h,g,o,j)|0;if((w|0)==3){x=14;break}else if((w|0)==2){n=-1;x=22;break}else if((w|0)!=1){y=q;break b}w=c[l>>2]|0;c[w>>2]=t;c[w+4>>2]=u;if((q|0)==8){n=-1;x=22;break}u=Nb(c[p>>2]|0)|0;if((u|0)==-1){n=-1;x=22;break}a[v]=u;q=q+1|0}if((x|0)==14){c[g>>2]=a[f]|0;y=q;break}else if((x|0)==22){i=e;return n|0}}else{c[g>>2]=a[f]|0;y=k}}while(0);if(d){d=c[g>>2]|0;c[b+48>>2]=d;n=d;i=e;return n|0}d=b+32|0;b=y;while(1){if((b|0)<=0){break}y=b+ -1|0;if((Va(a[f+y|0]|0,c[d>>2]|0)|0)==-1){n=-1;x=22;break}else{b=y}}if((x|0)==22){i=e;return n|0}n=c[g>>2]|0;i=e;return n|0}function Ld(a){a=a|0;var b=0;b=i;c[a>>2]=5528;_j(a+4|0);i=b;return}function Md(a){a=a|0;var b=0;b=i;c[a>>2]=5528;_j(a+4|0);Um(a);i=b;return}function Nd(b,d){b=b|0;d=d|0;var e=0,f=0;e=i;mc[c[(c[b>>2]|0)+24>>2]&63](b)|0;f=$j(d,8320)|0;c[b+36>>2]=f;a[b+44|0]=(mc[c[(c[f>>2]|0)+28>>2]&63](f)|0)&1;i=e;return}function Od(a){a=a|0;var b=0,d=0,e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0;b=i;i=i+16|0;d=b+8|0;e=b;f=a+36|0;g=a+40|0;h=d+8|0;j=d;k=a+32|0;while(1){a=c[f>>2]|0;l=wc[c[(c[a>>2]|0)+20>>2]&15](a,c[g>>2]|0,d,h,e)|0;a=(c[e>>2]|0)-j|0;if((Wb(d|0,1,a|0,c[k>>2]|0)|0)!=(a|0)){m=-1;n=5;break}if((l|0)==2){m=-1;n=5;break}else if((l|0)!=1){n=4;break}}if((n|0)==4){m=((_a(c[k>>2]|0)|0)!=0)<<31>>31;i=b;return m|0}else if((n|0)==5){i=b;return m|0}return 0}function Pd(b,e,f){b=b|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0;g=i;if((a[b+44|0]|0)!=0){h=Wb(e|0,1,f|0,c[b+32>>2]|0)|0;i=g;return h|0}if((f|0)>0){j=e;k=0}else{h=0;i=g;return h|0}while(1){if((vc[c[(c[b>>2]|0)+52>>2]&15](b,d[j]|0)|0)==-1){h=k;l=6;break}e=k+1|0;if((e|0)<(f|0)){j=j+1|0;k=e}else{h=e;l=6;break}}if((l|0)==6){i=g;return h|0}return 0}function Qd(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0;e=i;i=i+32|0;f=e+16|0;g=e+8|0;h=e;j=e+4|0;k=(d|0)==-1;a:do{if(!k){a[g]=d;if((a[b+44|0]|0)!=0){if((Wb(g|0,1,1,c[b+32>>2]|0)|0)==1){break}else{l=-1}i=e;return l|0}c[h>>2]=f;m=g+1|0;n=b+36|0;o=b+40|0;p=f+8|0;q=f;r=b+32|0;s=g;while(1){t=c[n>>2]|0;u=rc[c[(c[t>>2]|0)+12>>2]&15](t,c[o>>2]|0,s,m,j,f,p,h)|0;if((c[j>>2]|0)==(s|0)){l=-1;v=12;break}if((u|0)==3){v=7;break}t=(u|0)==1;if(!(u>>>0<2)){l=-1;v=12;break}u=(c[h>>2]|0)-q|0;if((Wb(f|0,1,u|0,c[r>>2]|0)|0)!=(u|0)){l=-1;v=12;break}if(t){s=t?c[j>>2]|0:s}else{break a}}if((v|0)==7){if((Wb(s|0,1,1,c[r>>2]|0)|0)==1){break}else{l=-1}i=e;return l|0}else if((v|0)==12){i=e;return l|0}}}while(0);l=k?0:d;i=e;return l|0}function Rd(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0;f=i;i=i+16|0;g=f;c[b>>2]=5528;h=b+4|0;Yj(h);j=b+8|0;c[j+0>>2]=0;c[j+4>>2]=0;c[j+8>>2]=0;c[j+12>>2]=0;c[j+16>>2]=0;c[j+20>>2]=0;c[b>>2]=4960;c[b+32>>2]=d;c[b+40>>2]=e;c[b+48>>2]=-1;a[b+52|0]=0;Zj(g,h);h=$j(g,8320)|0;e=b+36|0;c[e>>2]=h;d=b+44|0;c[d>>2]=mc[c[(c[h>>2]|0)+24>>2]&63](h)|0;h=c[e>>2]|0;a[b+53|0]=(mc[c[(c[h>>2]|0)+28>>2]&63](h)|0)&1;if((c[d>>2]|0)>8){jj(4800)}else{_j(g);i=f;return}}function Sd(a){a=a|0;var b=0;b=i;c[a>>2]=5528;_j(a+4|0);i=b;return}function Td(a){a=a|0;var b=0;b=i;c[a>>2]=5528;_j(a+4|0);Um(a);i=b;return}function Ud(b,d){b=b|0;d=d|0;var e=0,f=0,g=0;e=i;f=$j(d,8320)|0;d=b+36|0;c[d>>2]=f;g=b+44|0;c[g>>2]=mc[c[(c[f>>2]|0)+24>>2]&63](f)|0;f=c[d>>2]|0;a[b+53|0]=(mc[c[(c[f>>2]|0)+28>>2]&63](f)|0)&1;if((c[g>>2]|0)>8){jj(4800)}else{i=e;return}}function Vd(a){a=a|0;var b=0,c=0;b=i;c=Yd(a,0)|0;i=b;return c|0}function Wd(a){a=a|0;var b=0,c=0;b=i;c=Yd(a,1)|0;i=b;return c|0}function Xd(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;e=i;i=i+32|0;f=e+16|0;g=e;h=e+8|0;j=e+4|0;k=b+52|0;l=(a[k]|0)!=0;if((d|0)==-1){if(l){m=-1;i=e;return m|0}n=c[b+48>>2]|0;a[k]=(n|0)!=-1|0;m=n;i=e;return m|0}n=b+48|0;a:do{if(l){a[h]=c[n>>2];o=c[b+36>>2]|0;p=rc[c[(c[o>>2]|0)+12>>2]&15](o,c[b+40>>2]|0,h,h+1|0,j,f,f+8|0,g)|0;if((p|0)==1|(p|0)==2){m=-1;i=e;return m|0}else if((p|0)==3){a[f]=c[n>>2];c[g>>2]=f+1}p=b+32|0;while(1){o=c[g>>2]|0;if(!(o>>>0>f>>>0)){break a}q=o+ -1|0;c[g>>2]=q;if((Va(a[q]|0,c[p>>2]|0)|0)==-1){m=-1;break}}i=e;return m|0}}while(0);c[n>>2]=d;a[k]=1;m=d;i=e;return m|0}function Yd(b,e){b=b|0;e=e|0;var f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0;f=i;i=i+32|0;g=f+16|0;h=f+8|0;j=f;k=f+4|0;l=b+52|0;if((a[l]|0)!=0){m=b+48|0;n=c[m>>2]|0;if(!e){o=n;i=f;return o|0}c[m>>2]=-1;a[l]=0;o=n;i=f;return o|0}n=c[b+44>>2]|0;l=(n|0)>1?n:1;a:do{if((l|0)>0){n=b+32|0;m=0;while(1){p=Nb(c[n>>2]|0)|0;if((p|0)==-1){o=-1;break}a[g+m|0]=p;m=m+1|0;if((m|0)>=(l|0)){break a}}i=f;return o|0}}while(0);b:do{if((a[b+53|0]|0)==0){m=b+40|0;n=b+36|0;p=h+1|0;q=b+32|0;r=l;while(1){s=c[m>>2]|0;t=s;u=c[t>>2]|0;v=c[t+4>>2]|0;t=c[n>>2]|0;w=g+r|0;x=rc[c[(c[t>>2]|0)+16>>2]&15](t,s,g,w,j,h,p,k)|0;if((x|0)==2){o=-1;y=23;break}else if((x|0)==3){y=14;break}else if((x|0)!=1){z=r;break b}x=c[m>>2]|0;c[x>>2]=u;c[x+4>>2]=v;if((r|0)==8){o=-1;y=23;break}v=Nb(c[q>>2]|0)|0;if((v|0)==-1){o=-1;y=23;break}a[w]=v;r=r+1|0}if((y|0)==14){a[h]=a[g]|0;z=r;break}else if((y|0)==23){i=f;return o|0}}else{a[h]=a[g]|0;z=l}}while(0);do{if(!e){l=b+32|0;k=z;while(1){if((k|0)<=0){y=21;break}j=k+ -1|0;if((Va(d[g+j|0]|0,c[l>>2]|0)|0)==-1){o=-1;y=23;break}else{k=j}}if((y|0)==21){A=a[h]|0;break}else if((y|0)==23){i=f;return o|0}}else{k=a[h]|0;c[b+48>>2]=k&255;A=k}}while(0);o=A&255;i=f;return o|0}function Zd(){var a=0;a=i;vd(0);Ra(111,4224,p|0)|0;i=a;return}function _d(a){a=a|0;return}function $d(a){a=a|0;var b=0;b=a+4|0;c[b>>2]=(c[b>>2]|0)+1;return}function ae(a){a=a|0;var b=0,d=0,e=0,f=0;b=i;d=a+4|0;e=c[d>>2]|0;c[d>>2]=e+ -1;if((e|0)!=0){f=0;i=b;return f|0}jc[c[(c[a>>2]|0)+8>>2]&127](a);f=1;i=b;return f|0}function be(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,g=0;d=i;c[a>>2]=5104;e=mn(b|0)|0;f=Tm(e+13|0)|0;c[f+4>>2]=e;c[f>>2]=e;g=f+12|0;c[a+4>>2]=g;c[f+8>>2]=0;nn(g|0,b|0,e+1|0)|0;i=d;return}function ce(a){a=a|0;var b=0,d=0,e=0,f=0;b=i;c[a>>2]=5104;d=a+4|0;e=(c[d>>2]|0)+ -4|0;f=c[e>>2]|0;c[e>>2]=f+ -1;if((f+ -1|0)<0){Vm((c[d>>2]|0)+ -12|0)}Mb(a|0);Um(a);i=b;return}function de(a){a=a|0;var b=0,d=0,e=0,f=0;b=i;c[a>>2]=5104;d=a+4|0;e=(c[d>>2]|0)+ -4|0;f=c[e>>2]|0;c[e>>2]=f+ -1;if((f+ -1|0)>=0){Mb(a|0);i=b;return}Vm((c[d>>2]|0)+ -12|0);Mb(a|0);i=b;return}function ee(a){a=a|0;return c[a+4>>2]|0}function fe(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0;e=i;c[b>>2]=5128;if((a[d]&1)==0){f=d+1|0}else{f=c[d+8>>2]|0}d=mn(f|0)|0;g=Tm(d+13|0)|0;c[g+4>>2]=d;c[g>>2]=d;h=g+12|0;c[b+4>>2]=h;c[g+8>>2]=0;nn(h|0,f|0,d+1|0)|0;i=e;return}function ge(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,g=0;d=i;c[a>>2]=5128;e=mn(b|0)|0;f=Tm(e+13|0)|0;c[f+4>>2]=e;c[f>>2]=e;g=f+12|0;c[a+4>>2]=g;c[f+8>>2]=0;nn(g|0,b|0,e+1|0)|0;i=d;return}function he(a){a=a|0;var b=0,d=0,e=0,f=0;b=i;c[a>>2]=5128;d=a+4|0;e=(c[d>>2]|0)+ -4|0;f=c[e>>2]|0;c[e>>2]=f+ -1;if((f+ -1|0)<0){Vm((c[d>>2]|0)+ -12|0)}Mb(a|0);Um(a);i=b;return}function ie(a){a=a|0;var b=0,d=0,e=0,f=0;b=i;c[a>>2]=5128;d=a+4|0;e=(c[d>>2]|0)+ -4|0;f=c[e>>2]|0;c[e>>2]=f+ -1;if((f+ -1|0)>=0){Mb(a|0);i=b;return}Vm((c[d>>2]|0)+ -12|0);Mb(a|0);i=b;return}function je(a){a=a|0;return c[a+4>>2]|0}function ke(a){a=a|0;var b=0,d=0,e=0,f=0;b=i;c[a>>2]=5104;d=a+4|0;e=(c[d>>2]|0)+ -4|0;f=c[e>>2]|0;c[e>>2]=f+ -1;if((f+ -1|0)<0){Vm((c[d>>2]|0)+ -12|0)}Mb(a|0);Um(a);i=b;return}function le(a){a=a|0;return}function me(a,b,d){a=a|0;b=b|0;d=d|0;c[a>>2]=d;c[a+4>>2]=b;return}function ne(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,g=0;e=i;i=i+16|0;f=e;oc[c[(c[a>>2]|0)+12>>2]&3](f,a,b);if((c[f+4>>2]|0)!=(c[d+4>>2]|0)){g=0;i=e;return g|0}g=(c[f>>2]|0)==(c[d>>2]|0);i=e;return g|0}function oe(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0;e=i;if((c[b+4>>2]|0)!=(a|0)){f=0;i=e;return f|0}f=(c[b>>2]|0)==(d|0);i=e;return f|0}function pe(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0,k=0;d=i;f=Zb(e|0)|0;e=mn(f|0)|0;if(e>>>0>4294967279){ve(0)}if(e>>>0<11){a[b]=e<<1;g=b+1|0;nn(g|0,f|0,e|0)|0;h=g+e|0;a[h]=0;i=d;return}else{j=e+16&-16;k=Sm(j)|0;c[b+8>>2]=k;c[b>>2]=j|1;c[b+4>>2]=e;g=k;nn(g|0,f|0,e|0)|0;h=g+e|0;a[h]=0;i=d;return}}function qe(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0;f=i;i=i+16|0;g=f;h=c[d>>2]|0;if((h|0)!=0){j=a[e]|0;if((j&1)==0){k=(j&255)>>>1}else{k=c[e+4>>2]|0}if((k|0)==0){l=h}else{Fe(e,5280,2)|0;l=c[d>>2]|0}h=c[d+4>>2]|0;oc[c[(c[h>>2]|0)+24>>2]&3](g,h,l);l=a[g]|0;if((l&1)==0){m=g+1|0;n=(l&255)>>>1}else{m=c[g+8>>2]|0;n=c[g+4>>2]|0}Fe(e,m,n)|0;if(!((a[g]&1)==0)){Um(c[g+8>>2]|0)}}c[b+0>>2]=c[e+0>>2];c[b+4>>2]=c[e+4>>2];c[b+8>>2]=c[e+8>>2];c[e+0>>2]=0;c[e+4>>2]=0;c[e+8>>2]=0;i=f;return}function re(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0,k=0,l=0,m=0;f=i;i=i+32|0;g=f+12|0;h=f;j=mn(e|0)|0;if(j>>>0>4294967279){ve(0)}if(j>>>0<11){a[h]=j<<1;k=h+1|0}else{l=j+16&-16;m=Sm(l)|0;c[h+8>>2]=m;c[h>>2]=l|1;c[h+4>>2]=j;k=m}nn(k|0,e|0,j|0)|0;a[k+j|0]=0;qe(g,d,h);fe(b,g);if(!((a[g]&1)==0)){Um(c[g+8>>2]|0)}if(!((a[h]&1)==0)){Um(c[h+8>>2]|0)}c[b>>2]=5296;h=d;d=c[h+4>>2]|0;g=b+8|0;c[g>>2]=c[h>>2];c[g+4>>2]=d;i=f;return}function se(a){a=a|0;var b=0;b=i;ie(a);Um(a);i=b;return}function te(a){a=a|0;var b=0;b=i;ie(a);i=b;return}function ue(a,b,d){a=a|0;b=b|0;d=d|0;var e=0;e=i;Yb(5432)|0;if((c[a>>2]|0)==1){do{Da(5456,5432)|0}while((c[a>>2]|0)==1)}if((c[a>>2]|0)==0){c[a>>2]=1;Ab(5432)|0;jc[d&127](b);Yb(5432)|0;c[a>>2]=-1;Ab(5432)|0;Cb(5456)|0;i=e;return}else{Ab(5432)|0;i=e;return}}function ve(a){a=a|0;a=Wa(8)|0;be(a,5504);c[a>>2]=5184;xb(a|0,5224,9)}function we(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0;e=i;if((a[d]&1)==0){c[b+0>>2]=c[d+0>>2];c[b+4>>2]=c[d+4>>2];c[b+8>>2]=c[d+8>>2];i=e;return}f=c[d+8>>2]|0;g=c[d+4>>2]|0;if(g>>>0>4294967279){ve(0)}if(g>>>0<11){a[b]=g<<1;h=b+1|0}else{d=g+16&-16;j=Sm(d)|0;c[b+8>>2]=j;c[b>>2]=d|1;c[b+4>>2]=g;h=j}nn(h|0,f|0,g|0)|0;a[h+g|0]=0;i=e;return}function xe(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0;f=i;if(e>>>0>4294967279){ve(0)}if(e>>>0<11){a[b]=e<<1;g=b+1|0}else{h=e+16&-16;j=Sm(h)|0;c[b+8>>2]=j;c[b>>2]=h|1;c[b+4>>2]=e;g=j}nn(g|0,d|0,e|0)|0;a[g+e|0]=0;i=f;return}function ye(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0;f=i;if(d>>>0>4294967279){ve(0)}if(d>>>0<11){a[b]=d<<1;g=b+1|0}else{h=d+16&-16;j=Sm(h)|0;c[b+8>>2]=j;c[b>>2]=h|1;c[b+4>>2]=d;g=j}pn(g|0,e|0,d|0)|0;a[g+d|0]=0;i=f;return}function ze(b){b=b|0;var d=0;d=i;if((a[b]&1)==0){i=d;return}Um(c[b+8>>2]|0);i=d;return}function Ae(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0;e=i;f=mn(d|0)|0;g=a[b]|0;if((g&1)==0){h=g;j=10}else{g=c[b>>2]|0;h=g&255;j=(g&-2)+ -1|0}g=(h&1)==0;if(j>>>0<f>>>0){if(g){k=(h&255)>>>1}else{k=c[b+4>>2]|0}Ge(b,j,f-j|0,k,0,k,f,d);i=e;return b|0}if(g){l=b+1|0}else{l=c[b+8>>2]|0}on(l|0,d|0,f|0)|0;a[l+f|0]=0;if((a[b]&1)==0){a[b]=f<<1;i=e;return b|0}else{c[b+4>>2]=f;i=e;return b|0}return 0}function Be(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0;f=i;g=a[b]|0;h=(g&1)==0;if(h){j=(g&255)>>>1}else{j=c[b+4>>2]|0}if(j>>>0<d>>>0){Ce(b,d-j|0,e)|0;i=f;return}if(h){a[b+d+1|0]=0;a[b]=d<<1;i=f;return}else{a[(c[b+8>>2]|0)+d|0]=0;c[b+4>>2]=d;i=f;return}}function Ce(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0,k=0,l=0,m=0;f=i;if((d|0)==0){i=f;return b|0}g=a[b]|0;if((g&1)==0){h=10;j=g}else{g=c[b>>2]|0;h=(g&-2)+ -1|0;j=g&255}if((j&1)==0){k=(j&255)>>>1}else{k=c[b+4>>2]|0}if((h-k|0)>>>0<d>>>0){He(b,h,d-h+k|0,k,k,0,0);l=a[b]|0}else{l=j}if((l&1)==0){m=b+1|0}else{m=c[b+8>>2]|0}pn(m+k|0,e|0,d|0)|0;e=k+d|0;if((a[b]&1)==0){a[b]=e<<1}else{c[b+4>>2]=e}a[m+e|0]=0;i=f;return b|0}function De(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;e=i;if(d>>>0>4294967279){ve(0)}f=a[b]|0;if((f&1)==0){g=10;h=f}else{f=c[b>>2]|0;g=(f&-2)+ -1|0;h=f&255}if((h&1)==0){j=(h&255)>>>1}else{j=c[b+4>>2]|0}f=j>>>0>d>>>0?j:d;if(f>>>0<11){k=10}else{k=(f+16&-16)+ -1|0}if((k|0)==(g|0)){i=e;return}do{if((k|0)!=10){f=k+1|0;if(k>>>0>g>>>0){l=Sm(f)|0}else{l=Sm(f)|0}if((h&1)==0){m=l;n=1;o=b+1|0;p=0;break}else{m=l;n=1;o=c[b+8>>2]|0;p=1;break}}else{m=b+1|0;n=0;o=c[b+8>>2]|0;p=1}}while(0);if((h&1)==0){q=(h&255)>>>1}else{q=c[b+4>>2]|0}nn(m|0,o|0,q+1|0)|0;if(p){Um(o)}if(n){c[b>>2]=k+1|1;c[b+4>>2]=j;c[b+8>>2]=m;i=e;return}else{a[b]=j<<1;i=e;return}}function Ee(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0;e=i;f=a[b]|0;g=(f&1)!=0;if(g){h=(c[b>>2]&-2)+ -1|0;j=c[b+4>>2]|0}else{h=10;j=(f&255)>>>1}if((j|0)==(h|0)){He(b,h,1,h,h,0,0);if((a[b]&1)==0){k=7}else{k=8}}else{if(g){k=8}else{k=7}}if((k|0)==7){a[b]=(j<<1)+2;l=b+1|0;m=j+1|0;n=l+j|0;a[n]=d;o=l+m|0;a[o]=0;i=e;return}else if((k|0)==8){k=c[b+8>>2]|0;g=j+1|0;c[b+4>>2]=g;l=k;m=g;n=l+j|0;a[n]=d;o=l+m|0;a[o]=0;i=e;return}}function Fe(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0,k=0,l=0;f=i;g=a[b]|0;if((g&1)==0){h=10;j=g}else{g=c[b>>2]|0;h=(g&-2)+ -1|0;j=g&255}if((j&1)==0){k=(j&255)>>>1}else{k=c[b+4>>2]|0}if((h-k|0)>>>0<e>>>0){Ge(b,h,e-h+k|0,k,k,0,e,d);i=f;return b|0}if((e|0)==0){i=f;return b|0}if((j&1)==0){l=b+1|0}else{l=c[b+8>>2]|0}nn(l+k|0,d|0,e|0)|0;d=k+e|0;if((a[b]&1)==0){a[b]=d<<1}else{c[b+4>>2]=d}a[l+d|0]=0;i=f;return b|0}function Ge(b,d,e,f,g,h,j,k){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0;l=i;if((-18-d|0)>>>0<e>>>0){ve(0)}if((a[b]&1)==0){m=b+1|0}else{m=c[b+8>>2]|0}if(d>>>0<2147483623){n=e+d|0;e=d<<1;o=n>>>0<e>>>0?e:n;if(o>>>0<11){p=11}else{p=o+16&-16}}else{p=-17}o=Sm(p)|0;if((g|0)!=0){nn(o|0,m|0,g|0)|0}if((j|0)!=0){nn(o+g|0,k|0,j|0)|0}k=f-h|0;if((k|0)!=(g|0)){nn(o+(j+g)|0,m+(h+g)|0,k-g|0)|0}if((d|0)==10){q=b+8|0;c[q>>2]=o;r=p|1;c[b>>2]=r;s=k+j|0;t=b+4|0;c[t>>2]=s;u=o+s|0;a[u]=0;i=l;return}Um(m);q=b+8|0;c[q>>2]=o;r=p|1;c[b>>2]=r;s=k+j|0;t=b+4|0;c[t>>2]=s;u=o+s|0;a[u]=0;i=l;return}function He(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0;k=i;if((-17-d|0)>>>0<e>>>0){ve(0)}if((a[b]&1)==0){l=b+1|0}else{l=c[b+8>>2]|0}if(d>>>0<2147483623){m=e+d|0;e=d<<1;n=m>>>0<e>>>0?e:m;if(n>>>0<11){o=11}else{o=n+16&-16}}else{o=-17}n=Sm(o)|0;if((g|0)!=0){nn(n|0,l|0,g|0)|0}m=f-h|0;if((m|0)!=(g|0)){nn(n+(j+g)|0,l+(h+g)|0,m-g|0)|0}if((d|0)==10){p=b+8|0;c[p>>2]=n;q=o|1;c[b>>2]=q;i=k;return}Um(l);p=b+8|0;c[p>>2]=n;q=o|1;c[b>>2]=q;i=k;return}function Ie(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0;f=i;if(e>>>0>1073741807){ve(0)}if(e>>>0<2){a[b]=e<<1;g=b+4|0}else{h=e+4&-4;j=Sm(h<<2)|0;c[b+8>>2]=j;c[b>>2]=h|1;c[b+4>>2]=e;g=j}mm(g,d,e)|0;c[g+(e<<2)>>2]=0;i=f;return}function Je(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0;f=i;if(d>>>0>1073741807){ve(0)}if(d>>>0<2){a[b]=d<<1;g=b+4|0}else{h=d+4&-4;j=Sm(h<<2)|0;c[b+8>>2]=j;c[b>>2]=h|1;c[b+4>>2]=d;g=j}om(g,e,d)|0;c[g+(d<<2)>>2]=0;i=f;return}function Ke(b){b=b|0;var d=0;d=i;if((a[b]&1)==0){i=d;return}Um(c[b+8>>2]|0);i=d;return}function Le(a,b){a=a|0;b=b|0;var c=0,d=0;c=i;d=Me(a,b,lm(b)|0)|0;i=c;return d|0}function Me(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0,k=0,l=0;f=i;g=a[b]|0;if((g&1)==0){h=1;j=g}else{g=c[b>>2]|0;h=(g&-2)+ -1|0;j=g&255}g=(j&1)==0;if(h>>>0<e>>>0){if(g){k=(j&255)>>>1}else{k=c[b+4>>2]|0}Pe(b,h,e-h|0,k,0,k,e,d);i=f;return b|0}if(g){l=b+4|0}else{l=c[b+8>>2]|0}nm(l,d,e)|0;c[l+(e<<2)>>2]=0;if((a[b]&1)==0){a[b]=e<<1;i=f;return b|0}else{c[b+4>>2]=e;i=f;return b|0}return 0}function Ne(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;e=i;if(d>>>0>1073741807){ve(0)}f=a[b]|0;if((f&1)==0){g=1;h=f}else{f=c[b>>2]|0;g=(f&-2)+ -1|0;h=f&255}if((h&1)==0){j=(h&255)>>>1}else{j=c[b+4>>2]|0}f=j>>>0>d>>>0?j:d;if(f>>>0<2){k=1}else{k=(f+4&-4)+ -1|0}if((k|0)==(g|0)){i=e;return}do{if((k|0)!=1){f=(k<<2)+4|0;if(k>>>0>g>>>0){l=Sm(f)|0}else{l=Sm(f)|0}if((h&1)==0){m=l;n=1;o=b+4|0;p=0;break}else{m=l;n=1;o=c[b+8>>2]|0;p=1;break}}else{m=b+4|0;n=0;o=c[b+8>>2]|0;p=1}}while(0);if((h&1)==0){q=(h&255)>>>1}else{q=c[b+4>>2]|0}mm(m,o,q+1|0)|0;if(p){Um(o)}if(n){c[b>>2]=k+1|1;c[b+4>>2]=j;c[b+8>>2]=m;i=e;return}else{a[b]=j<<1;i=e;return}}function Oe(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0;e=i;f=a[b]|0;g=(f&1)!=0;if(g){h=(c[b>>2]&-2)+ -1|0;j=c[b+4>>2]|0}else{h=1;j=(f&255)>>>1}if((j|0)==(h|0)){Qe(b,h,1,h,h,0,0);if((a[b]&1)==0){k=7}else{k=8}}else{if(g){k=8}else{k=7}}if((k|0)==7){a[b]=(j<<1)+2;l=b+4|0;m=j+1|0;n=l+(j<<2)|0;c[n>>2]=d;o=l+(m<<2)|0;c[o>>2]=0;i=e;return}else if((k|0)==8){k=c[b+8>>2]|0;g=j+1|0;c[b+4>>2]=g;l=k;m=g;n=l+(j<<2)|0;c[n>>2]=d;o=l+(m<<2)|0;c[o>>2]=0;i=e;return}}function Pe(b,d,e,f,g,h,j,k){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0;l=i;if((1073741806-d|0)>>>0<e>>>0){ve(0)}if((a[b]&1)==0){m=b+4|0}else{m=c[b+8>>2]|0}if(d>>>0<536870887){n=e+d|0;e=d<<1;o=n>>>0<e>>>0?e:n;if(o>>>0<2){p=2}else{p=o+4&-4}}else{p=1073741807}o=Sm(p<<2)|0;if((g|0)!=0){mm(o,m,g)|0}if((j|0)!=0){mm(o+(g<<2)|0,k,j)|0}k=f-h|0;if((k|0)!=(g|0)){mm(o+(j+g<<2)|0,m+(h+g<<2)|0,k-g|0)|0}if((d|0)==1){q=b+8|0;c[q>>2]=o;r=p|1;c[b>>2]=r;s=k+j|0;t=b+4|0;c[t>>2]=s;u=o+(s<<2)|0;c[u>>2]=0;i=l;return}Um(m);q=b+8|0;c[q>>2]=o;r=p|1;c[b>>2]=r;s=k+j|0;t=b+4|0;c[t>>2]=s;u=o+(s<<2)|0;c[u>>2]=0;i=l;return}function Qe(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0;k=i;if((1073741807-d|0)>>>0<e>>>0){ve(0)}if((a[b]&1)==0){l=b+4|0}else{l=c[b+8>>2]|0}if(d>>>0<536870887){m=e+d|0;e=d<<1;n=m>>>0<e>>>0?e:m;if(n>>>0<2){o=2}else{o=n+4&-4}}else{o=1073741807}n=Sm(o<<2)|0;if((g|0)!=0){mm(n,l,g)|0}m=f-h|0;if((m|0)!=(g|0)){mm(n+(j+g<<2)|0,l+(h+g<<2)|0,m-g|0)|0}if((d|0)==1){p=b+8|0;c[p>>2]=n;q=o|1;c[b>>2]=q;i=k;return}Um(l);p=b+8|0;c[p>>2]=n;q=o|1;c[b>>2]=q;i=k;return}function Re(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0;e=i;i=i+16|0;f=e+8|0;g=e;h=(c[b+24>>2]|0)==0;if(h){c[b+16>>2]=d|1}else{c[b+16>>2]=d}if(((h&1|d)&c[b+20>>2]|0)==0){i=e;return}e=Wa(16)|0;if((a[5872]|0)==0?(Oa(5872)|0)!=0:0){c[1466]=6568;Ra(40,5864,p|0)|0;$a(5872)}b=g;c[b>>2]=1;c[b+4>>2]=5864;c[f+0>>2]=c[g+0>>2];c[f+4>>2]=c[g+4>>2];re(e,f,5920);c[e>>2]=5888;xb(e|0,5968,36)}function Se(a){a=a|0;var b=0,d=0,e=0,f=0,g=0;b=i;c[a>>2]=5912;d=c[a+40>>2]|0;e=a+32|0;f=a+36|0;if((d|0)!=0){g=d;do{g=g+ -1|0;oc[c[(c[e>>2]|0)+(g<<2)>>2]&3](0,a,c[(c[f>>2]|0)+(g<<2)>>2]|0)}while((g|0)!=0)}_j(a+28|0);Mm(c[e>>2]|0);Mm(c[f>>2]|0);Mm(c[a+48>>2]|0);Mm(c[a+60>>2]|0);i=b;return}function Te(a,b){a=a|0;b=b|0;var c=0;c=i;Zj(a,b+28|0);i=c;return}function Ue(a,b){a=a|0;b=b|0;var d=0,e=0;d=i;c[a+24>>2]=b;c[a+16>>2]=(b|0)==0;c[a+20>>2]=0;c[a+4>>2]=4098;c[a+12>>2]=0;c[a+8>>2]=6;b=a+28|0;e=a+32|0;a=e+40|0;do{c[e>>2]=0;e=e+4|0}while((e|0)<(a|0));Yj(b);i=d;return}function Ve(a){a=a|0;var b=0;b=i;c[a>>2]=5528;_j(a+4|0);Um(a);i=b;return}function We(a){a=a|0;var b=0;b=i;c[a>>2]=5528;_j(a+4|0);i=b;return}function Xe(a,b){a=a|0;b=b|0;return}function Ye(a,b,c){a=a|0;b=b|0;c=c|0;return a|0}function Ze(a,b,d,e,f,g){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;g=a;c[g>>2]=0;c[g+4>>2]=0;g=a+8|0;c[g>>2]=-1;c[g+4>>2]=-1;return}function _e(a,b,d,e){a=a|0;b=b|0;d=d|0;e=e|0;e=a;c[e>>2]=0;c[e+4>>2]=0;e=a+8|0;c[e>>2]=-1;c[e+4>>2]=-1;return}function $e(a){a=a|0;return 0}function af(a){a=a|0;return 0}function bf(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0;f=i;if((e|0)<=0){g=0;i=f;return g|0}h=b+12|0;j=b+16|0;k=d;d=0;while(1){l=c[h>>2]|0;if(l>>>0<(c[j>>2]|0)>>>0){c[h>>2]=l+1;m=a[l]|0}else{l=mc[c[(c[b>>2]|0)+40>>2]&63](b)|0;if((l|0)==-1){g=d;n=8;break}m=l&255}a[k]=m;l=d+1|0;if((l|0)<(e|0)){k=k+1|0;d=l}else{g=l;n=8;break}}if((n|0)==8){i=f;return g|0}return 0}function cf(a){a=a|0;return-1}function df(a){a=a|0;var b=0,e=0,f=0;b=i;if((mc[c[(c[a>>2]|0)+36>>2]&63](a)|0)==-1){e=-1;i=b;return e|0}f=a+12|0;a=c[f>>2]|0;c[f>>2]=a+1;e=d[a]|0;i=b;return e|0}function ef(a,b){a=a|0;b=b|0;return-1}function ff(b,e,f){b=b|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0;g=i;if((f|0)<=0){h=0;i=g;return h|0}j=b+24|0;k=b+28|0;l=e;e=0;while(1){m=c[j>>2]|0;if(!(m>>>0<(c[k>>2]|0)>>>0)){if((vc[c[(c[b>>2]|0)+52>>2]&15](b,d[l]|0)|0)==-1){h=e;n=7;break}}else{o=a[l]|0;c[j>>2]=m+1;a[m]=o}o=e+1|0;if((o|0)<(f|0)){l=l+1|0;e=o}else{h=o;n=7;break}}if((n|0)==7){i=g;return h|0}return 0}function gf(a,b){a=a|0;b=b|0;return-1}function hf(a){a=a|0;var b=0;b=i;c[a>>2]=5592;_j(a+4|0);Um(a);i=b;return}function jf(a){a=a|0;var b=0;b=i;c[a>>2]=5592;_j(a+4|0);i=b;return}function kf(a,b){a=a|0;b=b|0;return}function lf(a,b,c){a=a|0;b=b|0;c=c|0;return a|0}function mf(a,b,d,e,f,g){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;g=a;c[g>>2]=0;c[g+4>>2]=0;g=a+8|0;c[g>>2]=-1;c[g+4>>2]=-1;return}function nf(a,b,d,e){a=a|0;b=b|0;d=d|0;e=e|0;e=a;c[e>>2]=0;c[e+4>>2]=0;e=a+8|0;c[e>>2]=-1;c[e+4>>2]=-1;return}function of(a){a=a|0;return 0}function pf(a){a=a|0;return 0}function qf(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0;e=i;if((d|0)<=0){f=0;i=e;return f|0}g=a+12|0;h=a+16|0;j=b;b=0;while(1){k=c[g>>2]|0;if(!(k>>>0<(c[h>>2]|0)>>>0)){l=mc[c[(c[a>>2]|0)+40>>2]&63](a)|0;if((l|0)==-1){f=b;m=8;break}else{n=l}}else{c[g>>2]=k+4;n=c[k>>2]|0}c[j>>2]=n;k=b+1|0;if((k|0)>=(d|0)){f=k;m=8;break}j=j+4|0;b=k}if((m|0)==8){i=e;return f|0}return 0}function rf(a){a=a|0;return-1}function sf(a){a=a|0;var b=0,d=0,e=0;b=i;if((mc[c[(c[a>>2]|0)+36>>2]&63](a)|0)==-1){d=-1;i=b;return d|0}e=a+12|0;a=c[e>>2]|0;c[e>>2]=a+4;d=c[a>>2]|0;i=b;return d|0}function tf(a,b){a=a|0;b=b|0;return-1}function uf(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0;e=i;if((d|0)<=0){f=0;i=e;return f|0}g=a+24|0;h=a+28|0;j=b;b=0;while(1){k=c[g>>2]|0;if(!(k>>>0<(c[h>>2]|0)>>>0)){if((vc[c[(c[a>>2]|0)+52>>2]&15](a,c[j>>2]|0)|0)==-1){f=b;l=8;break}}else{m=c[j>>2]|0;c[g>>2]=k+4;c[k>>2]=m}m=b+1|0;if((m|0)>=(d|0)){f=m;l=8;break}j=j+4|0;b=m}if((l|0)==8){i=e;return f|0}return 0}function vf(a,b){a=a|0;b=b|0;return-1}function wf(a){a=a|0;var b=0;b=i;Se(a+8|0);Um(a);i=b;return}function xf(a){a=a|0;var b=0;b=i;Se(a+8|0);i=b;return}function yf(a){a=a|0;var b=0,d=0;b=i;d=c[(c[a>>2]|0)+ -12>>2]|0;Se(a+(d+8)|0);Um(a+d|0);i=b;return}function zf(a){a=a|0;var b=0;b=i;Se(a+((c[(c[a>>2]|0)+ -12>>2]|0)+8)|0);i=b;return}function Af(b){b=b|0;var d=0,e=0,f=0,g=0,h=0;d=i;i=i+16|0;e=d;f=c[(c[b>>2]|0)+ -12>>2]|0;if((c[b+(f+24)>>2]|0)==0){i=d;return b|0}a[e]=0;c[e+4>>2]=b;if((c[b+(f+16)>>2]|0)==0){g=c[b+(f+72)>>2]|0;if((g|0)==0){h=f}else{Af(g)|0;h=c[(c[b>>2]|0)+ -12>>2]|0}a[e]=1;g=c[b+(h+24)>>2]|0;if((mc[c[(c[g>>2]|0)+24>>2]&63](g)|0)==-1){g=c[(c[b>>2]|0)+ -12>>2]|0;Re(b+g|0,c[b+(g+16)>>2]|1)}}Kf(e);i=d;return b|0}function Bf(a){a=a|0;var b=0;b=i;Se(a+8|0);Um(a);i=b;return}function Cf(a){a=a|0;var b=0;b=i;Se(a+8|0);i=b;return}function Df(a){a=a|0;var b=0,d=0;b=i;d=c[(c[a>>2]|0)+ -12>>2]|0;Se(a+(d+8)|0);Um(a+d|0);i=b;return}function Ef(a){a=a|0;var b=0;b=i;Se(a+((c[(c[a>>2]|0)+ -12>>2]|0)+8)|0);i=b;return}function Ff(b){b=b|0;var d=0,e=0,f=0,g=0,h=0;d=i;i=i+16|0;e=d;f=c[(c[b>>2]|0)+ -12>>2]|0;if((c[b+(f+24)>>2]|0)==0){i=d;return b|0}a[e]=0;c[e+4>>2]=b;if((c[b+(f+16)>>2]|0)==0){g=c[b+(f+72)>>2]|0;if((g|0)==0){h=f}else{Ff(g)|0;h=c[(c[b>>2]|0)+ -12>>2]|0}a[e]=1;g=c[b+(h+24)>>2]|0;if((mc[c[(c[g>>2]|0)+24>>2]&63](g)|0)==-1){g=c[(c[b>>2]|0)+ -12>>2]|0;Re(b+g|0,c[b+(g+16)>>2]|1)}}Pf(e);i=d;return b|0}function Gf(a){a=a|0;var b=0;b=i;Se(a+4|0);Um(a);i=b;return}function Hf(a){a=a|0;var b=0;b=i;Se(a+4|0);i=b;return}function If(a){a=a|0;var b=0,d=0;b=i;d=c[(c[a>>2]|0)+ -12>>2]|0;Se(a+(d+4)|0);Um(a+d|0);i=b;return}function Jf(a){a=a|0;var b=0;b=i;Se(a+((c[(c[a>>2]|0)+ -12>>2]|0)+4)|0);i=b;return}function Kf(a){a=a|0;var b=0,d=0,e=0;b=i;d=a+4|0;a=c[d>>2]|0;e=c[(c[a>>2]|0)+ -12>>2]|0;if((c[a+(e+24)>>2]|0)==0){i=b;return}if((c[a+(e+16)>>2]|0)!=0){i=b;return}if((c[a+(e+4)>>2]&8192|0)==0){i=b;return}if(Xa()|0){i=b;return}e=c[d>>2]|0;a=c[e+((c[(c[e>>2]|0)+ -12>>2]|0)+24)>>2]|0;if(!((mc[c[(c[a>>2]|0)+24>>2]&63](a)|0)==-1)){i=b;return}a=c[d>>2]|0;d=c[(c[a>>2]|0)+ -12>>2]|0;Re(a+d|0,c[a+(d+16)>>2]|1);i=b;return}function Lf(a){a=a|0;var b=0;b=i;Se(a+4|0);Um(a);i=b;return}function Mf(a){a=a|0;var b=0;b=i;Se(a+4|0);i=b;return}function Nf(a){a=a|0;var b=0,d=0;b=i;d=c[(c[a>>2]|0)+ -12>>2]|0;Se(a+(d+4)|0);Um(a+d|0);i=b;return}function Of(a){a=a|0;var b=0;b=i;Se(a+((c[(c[a>>2]|0)+ -12>>2]|0)+4)|0);i=b;return}function Pf(a){a=a|0;var b=0,d=0,e=0;b=i;d=a+4|0;a=c[d>>2]|0;e=c[(c[a>>2]|0)+ -12>>2]|0;if((c[a+(e+24)>>2]|0)==0){i=b;return}if((c[a+(e+16)>>2]|0)!=0){i=b;return}if((c[a+(e+4)>>2]&8192|0)==0){i=b;return}if(Xa()|0){i=b;return}e=c[d>>2]|0;a=c[e+((c[(c[e>>2]|0)+ -12>>2]|0)+24)>>2]|0;if(!((mc[c[(c[a>>2]|0)+24>>2]&63](a)|0)==-1)){i=b;return}a=c[d>>2]|0;d=c[(c[a>>2]|0)+ -12>>2]|0;Re(a+d|0,c[a+(d+16)>>2]|1);i=b;return}function Qf(a){a=a|0;return 5808}function Rf(a,b,c){a=a|0;b=b|0;c=c|0;var d=0;d=i;if((c|0)==1){xe(a,5824,35);i=d;return}else{pe(a,b,c);i=d;return}}function Sf(a){a=a|0;return}function Tf(a){a=a|0;var b=0;b=i;te(a);Um(a);i=b;return}function Uf(a){a=a|0;var b=0;b=i;te(a);i=b;return}function Vf(a){a=a|0;var b=0;b=i;Se(a);Um(a);i=b;return}function Wf(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function Xf(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function Yf(a){a=a|0;return}function Zf(a){a=a|0;return}function _f(b,c,d,e,f){b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0;b=i;a:do{if((e|0)==(f|0)){g=c;h=6}else{j=e;k=c;while(1){if((k|0)==(d|0)){l=-1;break a}m=a[k]|0;n=a[j]|0;if(m<<24>>24<n<<24>>24){l=-1;break a}if(n<<24>>24<m<<24>>24){l=1;break a}m=k+1|0;n=j+1|0;if((n|0)==(f|0)){g=m;h=6;break}else{j=n;k=m}}}}while(0);if((h|0)==6){l=(g|0)!=(d|0)|0}i=b;return l|0}function $f(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0;d=i;g=e;h=f-g|0;if(h>>>0>4294967279){ve(b)}if(h>>>0<11){a[b]=h<<1;j=b+1|0}else{k=h+16&-16;l=Sm(k)|0;c[b+8>>2]=l;c[b>>2]=k|1;c[b+4>>2]=h;j=l}if((e|0)==(f|0)){m=j;a[m]=0;i=d;return}else{n=e;o=j}while(1){a[o]=a[n]|0;n=n+1|0;if((n|0)==(f|0)){break}else{o=o+1|0}}m=j+(f+(0-g))|0;a[m]=0;i=d;return}function ag(b,c,d){b=b|0;c=c|0;d=d|0;var e=0,f=0,g=0,h=0,j=0;b=i;if((c|0)==(d|0)){e=0;i=b;return e|0}else{f=0;g=c}while(1){c=(a[g]|0)+(f<<4)|0;h=c&-268435456;j=(h>>>24|h)^c;c=g+1|0;if((c|0)==(d|0)){e=j;break}else{f=j;g=c}}i=b;return e|0}function bg(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function cg(a){a=a|0;return}function dg(a,b,d,e,f){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0;a=i;a:do{if((e|0)==(f|0)){g=b;h=6}else{j=e;k=b;while(1){if((k|0)==(d|0)){l=-1;break a}m=c[k>>2]|0;n=c[j>>2]|0;if((m|0)<(n|0)){l=-1;break a}if((n|0)<(m|0)){l=1;break a}m=k+4|0;n=j+4|0;if((n|0)==(f|0)){g=m;h=6;break}else{j=n;k=m}}}}while(0);if((h|0)==6){l=(g|0)!=(d|0)|0}i=a;return l|0}function eg(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0;d=i;g=e;h=f-g|0;j=h>>2;if(j>>>0>1073741807){ve(b)}if(j>>>0<2){a[b]=h>>>1;k=b+4|0}else{h=j+4&-4;l=Sm(h<<2)|0;c[b+8>>2]=l;c[b>>2]=h|1;c[b+4>>2]=j;k=l}if((e|0)==(f|0)){m=k;c[m>>2]=0;i=d;return}l=f+ -4+(0-g)|0;g=e;e=k;while(1){c[e>>2]=c[g>>2];g=g+4|0;if((g|0)==(f|0)){break}else{e=e+4|0}}m=k+((l>>>2)+1<<2)|0;c[m>>2]=0;i=d;return}function fg(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0;a=i;if((b|0)==(d|0)){e=0;i=a;return e|0}else{f=0;g=b}while(1){b=(c[g>>2]|0)+(f<<4)|0;h=b&-268435456;j=(h>>>24|h)^b;b=g+4|0;if((b|0)==(d|0)){e=j;break}else{f=j;g=b}}i=a;return e|0}function gg(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function hg(a){a=a|0;return}function ig(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0;k=i;i=i+80|0;l=k;m=k+64|0;n=k+60|0;o=k+56|0;p=k+52|0;q=k+68|0;r=k+16|0;s=k+12|0;t=k+24|0;u=k+48|0;if((c[g+4>>2]&1|0)==0){c[n>>2]=-1;v=c[(c[d>>2]|0)+16>>2]|0;c[p>>2]=c[e>>2];c[q>>2]=c[f>>2];c[m+0>>2]=c[p+0>>2];c[l+0>>2]=c[q+0>>2];hc[v&63](o,d,m,l,g,h,n);m=c[o>>2]|0;c[e>>2]=m;o=c[n>>2]|0;if((o|0)==1){a[j]=1}else if((o|0)==0){a[j]=0}else{a[j]=1;c[h>>2]=4}c[b>>2]=m;i=k;return}Te(r,g);m=c[r>>2]|0;if(!((c[2064]|0)==-1)){c[l>>2]=8256;c[l+4>>2]=112;c[l+8>>2]=0;ue(8256,l,113)}o=(c[8260>>2]|0)+ -1|0;n=c[m+8>>2]|0;if(!((c[m+12>>2]|0)-n>>2>>>0>o>>>0)){w=Wa(4)|0;qm(w);xb(w|0,16216,101)}m=c[n+(o<<2)>>2]|0;if((m|0)==0){w=Wa(4)|0;qm(w);xb(w|0,16216,101)}ae(c[r>>2]|0)|0;Te(s,g);g=c[s>>2]|0;if(!((c[2100]|0)==-1)){c[l>>2]=8400;c[l+4>>2]=112;c[l+8>>2]=0;ue(8400,l,113)}r=(c[8404>>2]|0)+ -1|0;w=c[g+8>>2]|0;if(!((c[g+12>>2]|0)-w>>2>>>0>r>>>0)){x=Wa(4)|0;qm(x);xb(x|0,16216,101)}g=c[w+(r<<2)>>2]|0;if((g|0)==0){x=Wa(4)|0;qm(x);xb(x|0,16216,101)}ae(c[s>>2]|0)|0;kc[c[(c[g>>2]|0)+24>>2]&63](t,g);kc[c[(c[g>>2]|0)+28>>2]&63](t+12|0,g);c[u>>2]=c[f>>2];f=t+24|0;c[l+0>>2]=c[u+0>>2];a[j]=(jg(e,l,t,f,m,h,1)|0)==(t|0)|0;c[b>>2]=c[e>>2];ze(t+12|0);ze(t);i=k;return}function jg(b,e,f,g,h,j,k){b=b|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0,X=0,Y=0,Z=0,_=0;l=i;i=i+112|0;m=l;n=(g-f|0)/12|0;if(n>>>0>100){o=Lm(n)|0;if((o|0)==0){_m()}else{p=o;q=o}}else{p=0;q=m}m=(f|0)==(g|0);if(m){r=0;s=n}else{o=f;t=0;u=n;n=q;while(1){v=a[o]|0;if((v&1)==0){w=(v&255)>>>1}else{w=c[o+4>>2]|0}if((w|0)==0){a[n]=2;x=t+1|0;y=u+ -1|0}else{a[n]=1;x=t;y=u}v=o+12|0;if((v|0)==(g|0)){r=x;s=y;break}else{o=v;t=x;u=y;n=n+1|0}}}n=0;y=r;r=s;a:while(1){s=c[b>>2]|0;do{if((s|0)!=0){if((c[s+12>>2]|0)==(c[s+16>>2]|0)){if((mc[c[(c[s>>2]|0)+36>>2]&63](s)|0)==-1){c[b>>2]=0;z=0;break}else{z=c[b>>2]|0;break}}else{z=s}}else{z=0}}while(0);s=(z|0)==0;u=c[e>>2]|0;if((u|0)!=0){if((c[u+12>>2]|0)==(c[u+16>>2]|0)?(mc[c[(c[u>>2]|0)+36>>2]&63](u)|0)==-1:0){c[e>>2]=0;A=0}else{A=u}}else{A=0}B=(A|0)==0;C=c[b>>2]|0;if(!((s^B)&(r|0)!=0)){break}s=c[C+12>>2]|0;if((s|0)==(c[C+16>>2]|0)){D=mc[c[(c[C>>2]|0)+36>>2]&63](C)|0}else{D=d[s]|0}s=D&255;if(k){E=s}else{E=vc[c[(c[h>>2]|0)+12>>2]&15](h,s)|0}s=n+1|0;if(m){n=s;continue}b:do{if(k){u=0;x=f;t=y;o=r;w=q;while(1){do{if((a[w]|0)==1){v=a[x]|0;F=(v&1)==0;if(F){G=x+1|0}else{G=c[x+8>>2]|0}if(!(E<<24>>24==(a[G+n|0]|0))){a[w]=0;H=u;I=t;J=o+ -1|0;break}if(F){K=(v&255)>>>1}else{K=c[x+4>>2]|0}if((K|0)==(s|0)){a[w]=2;H=1;I=t+1|0;J=o+ -1|0}else{H=1;I=t;J=o}}else{H=u;I=t;J=o}}while(0);v=x+12|0;if((v|0)==(g|0)){L=H;M=I;N=J;break b}u=H;x=v;t=I;o=J;w=w+1|0}}else{w=0;o=f;t=y;x=r;u=q;while(1){do{if((a[u]|0)==1){if((a[o]&1)==0){O=o+1|0}else{O=c[o+8>>2]|0}if(!(E<<24>>24==(vc[c[(c[h>>2]|0)+12>>2]&15](h,a[O+n|0]|0)|0)<<24>>24)){a[u]=0;P=w;Q=t;R=x+ -1|0;break}v=a[o]|0;if((v&1)==0){S=(v&255)>>>1}else{S=c[o+4>>2]|0}if((S|0)==(s|0)){a[u]=2;P=1;Q=t+1|0;R=x+ -1|0}else{P=1;Q=t;R=x}}else{P=w;Q=t;R=x}}while(0);v=o+12|0;if((v|0)==(g|0)){L=P;M=Q;N=R;break b}w=P;o=v;t=Q;x=R;u=u+1|0}}}while(0);if(!L){n=s;y=M;r=N;continue}u=c[b>>2]|0;x=u+12|0;t=c[x>>2]|0;if((t|0)==(c[u+16>>2]|0)){mc[c[(c[u>>2]|0)+40>>2]&63](u)|0}else{c[x>>2]=t+1}if((N+M|0)>>>0<2){n=s;y=M;r=N;continue}else{T=f;U=M;V=q}while(1){if((a[V]|0)==2){t=a[T]|0;if((t&1)==0){W=(t&255)>>>1}else{W=c[T+4>>2]|0}if((W|0)!=(s|0)){a[V]=0;X=U+ -1|0}else{X=U}}else{X=U}t=T+12|0;if((t|0)==(g|0)){n=s;y=X;r=N;continue a}else{T=t;U=X;V=V+1|0}}}do{if((C|0)!=0){if((c[C+12>>2]|0)==(c[C+16>>2]|0)){if((mc[c[(c[C>>2]|0)+36>>2]&63](C)|0)==-1){c[b>>2]=0;Y=0;break}else{Y=c[b>>2]|0;break}}else{Y=C}}else{Y=0}}while(0);C=(Y|0)==0;do{if(!B){if((c[A+12>>2]|0)!=(c[A+16>>2]|0)){if(C){break}else{Z=80;break}}if(!((mc[c[(c[A>>2]|0)+36>>2]&63](A)|0)==-1)){if(C){break}else{Z=80;break}}else{c[e>>2]=0;Z=78;break}}else{Z=78}}while(0);if((Z|0)==78?C:0){Z=80}if((Z|0)==80){c[j>>2]=c[j>>2]|2}c:do{if(!m){if((a[q]|0)==2){_=f}else{C=f;e=q;while(1){A=C+12|0;B=e+1|0;if((A|0)==(g|0)){Z=85;break c}if((a[B]|0)==2){_=A;break}else{C=A;e=B}}}}else{Z=85}}while(0);if((Z|0)==85){c[j>>2]=c[j>>2]|4;_=g}if((p|0)==0){i=l;return _|0}Mm(p);i=l;return _|0}function kg(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];lg(a,0,k,j,f,g,h);i=b;return}function lg(b,e,f,g,h,j,k){b=b|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0;e=i;i=i+224|0;l=e+198|0;m=e+196|0;n=e+4|0;o=e+16|0;p=e+28|0;q=e+32|0;r=e;s=e+192|0;t=c[h+4>>2]&74;if((t|0)==0){u=0}else if((t|0)==8){u=16}else if((t|0)==64){u=8}else{u=10}bh(n,h,l,m);c[o+0>>2]=0;c[o+4>>2]=0;c[o+8>>2]=0;Be(o,10,0);if((a[o]&1)==0){h=o+1|0;v=h;w=o+8|0;x=h}else{h=o+8|0;v=o+1|0;w=h;x=c[h>>2]|0}c[p>>2]=x;c[r>>2]=q;c[s>>2]=0;h=o+4|0;t=a[m]|0;m=c[f>>2]|0;y=x;a:while(1){if((m|0)!=0){if((c[m+12>>2]|0)==(c[m+16>>2]|0)?(mc[c[(c[m>>2]|0)+36>>2]&63](m)|0)==-1:0){c[f>>2]=0;z=0}else{z=m}}else{z=0}x=(z|0)==0;A=c[g>>2]|0;do{if((A|0)!=0){if((c[A+12>>2]|0)!=(c[A+16>>2]|0)){if(x){B=A;break}else{C=A;D=y;break a}}if(!((mc[c[(c[A>>2]|0)+36>>2]&63](A)|0)==-1)){if(x){B=A;break}else{C=A;D=y;break a}}else{c[g>>2]=0;E=18;break}}else{E=18}}while(0);if((E|0)==18){E=0;if(x){C=0;D=y;break}else{B=0}}A=a[o]|0;F=(A&1)==0;if(F){G=(A&255)>>>1}else{G=c[h>>2]|0}if(((c[p>>2]|0)-y|0)==(G|0)){if(F){H=(A&255)>>>1;I=(A&255)>>>1}else{A=c[h>>2]|0;H=A;I=A}Be(o,H<<1,0);if((a[o]&1)==0){J=10}else{J=(c[o>>2]&-2)+ -1|0}Be(o,J,0);if((a[o]&1)==0){K=v}else{K=c[w>>2]|0}c[p>>2]=K+I;L=K}else{L=y}A=z+12|0;F=c[A>>2]|0;M=z+16|0;if((F|0)==(c[M>>2]|0)){N=mc[c[(c[z>>2]|0)+36>>2]&63](z)|0}else{N=d[F]|0}if((Dg(N&255,u,L,p,s,t,n,q,r,l)|0)!=0){C=B;D=L;break}F=c[A>>2]|0;if((F|0)==(c[M>>2]|0)){mc[c[(c[z>>2]|0)+40>>2]&63](z)|0;m=z;y=L;continue}else{c[A>>2]=F+1;m=z;y=L;continue}}L=a[n]|0;if((L&1)==0){O=(L&255)>>>1}else{O=c[n+4>>2]|0}if((O|0)!=0?(O=c[r>>2]|0,(O-q|0)<160):0){L=c[s>>2]|0;c[r>>2]=O+4;c[O>>2]=L}c[k>>2]=Wl(D,c[p>>2]|0,j,u)|0;nj(n,q,c[r>>2]|0,j);if((z|0)!=0){if((c[z+12>>2]|0)==(c[z+16>>2]|0)?(mc[c[(c[z>>2]|0)+36>>2]&63](z)|0)==-1:0){c[f>>2]=0;P=0}else{P=z}}else{P=0}z=(P|0)==0;do{if((C|0)!=0){if((c[C+12>>2]|0)!=(c[C+16>>2]|0)){if(!z){break}c[b>>2]=P;ze(o);ze(n);i=e;return}if((mc[c[(c[C>>2]|0)+36>>2]&63](C)|0)==-1){c[g>>2]=0;E=54;break}if(z^(C|0)==0){c[b>>2]=P;ze(o);ze(n);i=e;return}}else{E=54}}while(0);if((E|0)==54?!z:0){c[b>>2]=P;ze(o);ze(n);i=e;return}c[j>>2]=c[j>>2]|2;c[b>>2]=P;ze(o);ze(n);i=e;return}function mg(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];ng(a,0,k,j,f,g,h);i=b;return}function ng(b,e,f,g,h,j,k){b=b|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0;e=i;i=i+224|0;l=e+198|0;m=e+196|0;n=e+4|0;o=e+16|0;p=e+28|0;q=e+32|0;r=e;s=e+192|0;t=c[h+4>>2]&74;if((t|0)==8){u=16}else if((t|0)==0){u=0}else if((t|0)==64){u=8}else{u=10}bh(n,h,l,m);c[o+0>>2]=0;c[o+4>>2]=0;c[o+8>>2]=0;Be(o,10,0);if((a[o]&1)==0){h=o+1|0;v=h;w=o+8|0;x=h}else{h=o+8|0;v=o+1|0;w=h;x=c[h>>2]|0}c[p>>2]=x;c[r>>2]=q;c[s>>2]=0;h=o+4|0;t=a[m]|0;m=c[f>>2]|0;y=x;a:while(1){if((m|0)!=0){if((c[m+12>>2]|0)==(c[m+16>>2]|0)?(mc[c[(c[m>>2]|0)+36>>2]&63](m)|0)==-1:0){c[f>>2]=0;z=0}else{z=m}}else{z=0}x=(z|0)==0;A=c[g>>2]|0;do{if((A|0)!=0){if((c[A+12>>2]|0)!=(c[A+16>>2]|0)){if(x){B=A;break}else{C=A;D=y;break a}}if(!((mc[c[(c[A>>2]|0)+36>>2]&63](A)|0)==-1)){if(x){B=A;break}else{C=A;D=y;break a}}else{c[g>>2]=0;E=18;break}}else{E=18}}while(0);if((E|0)==18){E=0;if(x){C=0;D=y;break}else{B=0}}A=a[o]|0;F=(A&1)==0;if(F){G=(A&255)>>>1}else{G=c[h>>2]|0}if(((c[p>>2]|0)-y|0)==(G|0)){if(F){H=(A&255)>>>1;J=(A&255)>>>1}else{A=c[h>>2]|0;H=A;J=A}Be(o,H<<1,0);if((a[o]&1)==0){K=10}else{K=(c[o>>2]&-2)+ -1|0}Be(o,K,0);if((a[o]&1)==0){L=v}else{L=c[w>>2]|0}c[p>>2]=L+J;M=L}else{M=y}A=z+12|0;F=c[A>>2]|0;N=z+16|0;if((F|0)==(c[N>>2]|0)){O=mc[c[(c[z>>2]|0)+36>>2]&63](z)|0}else{O=d[F]|0}if((Dg(O&255,u,M,p,s,t,n,q,r,l)|0)!=0){C=B;D=M;break}F=c[A>>2]|0;if((F|0)==(c[N>>2]|0)){mc[c[(c[z>>2]|0)+40>>2]&63](z)|0;m=z;y=M;continue}else{c[A>>2]=F+1;m=z;y=M;continue}}M=a[n]|0;if((M&1)==0){P=(M&255)>>>1}else{P=c[n+4>>2]|0}if((P|0)!=0?(P=c[r>>2]|0,(P-q|0)<160):0){M=c[s>>2]|0;c[r>>2]=P+4;c[P>>2]=M}M=Vl(D,c[p>>2]|0,j,u)|0;u=k;c[u>>2]=M;c[u+4>>2]=I;nj(n,q,c[r>>2]|0,j);if((z|0)!=0){if((c[z+12>>2]|0)==(c[z+16>>2]|0)?(mc[c[(c[z>>2]|0)+36>>2]&63](z)|0)==-1:0){c[f>>2]=0;Q=0}else{Q=z}}else{Q=0}z=(Q|0)==0;do{if((C|0)!=0){if((c[C+12>>2]|0)!=(c[C+16>>2]|0)){if(!z){break}c[b>>2]=Q;ze(o);ze(n);i=e;return}if((mc[c[(c[C>>2]|0)+36>>2]&63](C)|0)==-1){c[g>>2]=0;E=54;break}if(z^(C|0)==0){c[b>>2]=Q;ze(o);ze(n);i=e;return}}else{E=54}}while(0);if((E|0)==54?!z:0){c[b>>2]=Q;ze(o);ze(n);i=e;return}c[j>>2]=c[j>>2]|2;c[b>>2]=Q;ze(o);ze(n);i=e;return}function og(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];pg(a,0,k,j,f,g,h);i=b;return}function pg(e,f,g,h,j,k,l){e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;var m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0;f=i;i=i+224|0;m=f+198|0;n=f+196|0;o=f+4|0;p=f+16|0;q=f+28|0;r=f+32|0;s=f;t=f+192|0;u=c[j+4>>2]&74;if((u|0)==0){v=0}else if((u|0)==8){v=16}else if((u|0)==64){v=8}else{v=10}bh(o,j,m,n);c[p+0>>2]=0;c[p+4>>2]=0;c[p+8>>2]=0;Be(p,10,0);if((a[p]&1)==0){j=p+1|0;w=j;x=p+8|0;y=j}else{j=p+8|0;w=p+1|0;x=j;y=c[j>>2]|0}c[q>>2]=y;c[s>>2]=r;c[t>>2]=0;j=p+4|0;u=a[n]|0;n=c[g>>2]|0;z=y;a:while(1){if((n|0)!=0){if((c[n+12>>2]|0)==(c[n+16>>2]|0)?(mc[c[(c[n>>2]|0)+36>>2]&63](n)|0)==-1:0){c[g>>2]=0;A=0}else{A=n}}else{A=0}y=(A|0)==0;B=c[h>>2]|0;do{if((B|0)!=0){if((c[B+12>>2]|0)!=(c[B+16>>2]|0)){if(y){C=B;break}else{D=B;E=z;break a}}if(!((mc[c[(c[B>>2]|0)+36>>2]&63](B)|0)==-1)){if(y){C=B;break}else{D=B;E=z;break a}}else{c[h>>2]=0;F=18;break}}else{F=18}}while(0);if((F|0)==18){F=0;if(y){D=0;E=z;break}else{C=0}}B=a[p]|0;G=(B&1)==0;if(G){H=(B&255)>>>1}else{H=c[j>>2]|0}if(((c[q>>2]|0)-z|0)==(H|0)){if(G){I=(B&255)>>>1;J=(B&255)>>>1}else{B=c[j>>2]|0;I=B;J=B}Be(p,I<<1,0);if((a[p]&1)==0){K=10}else{K=(c[p>>2]&-2)+ -1|0}Be(p,K,0);if((a[p]&1)==0){L=w}else{L=c[x>>2]|0}c[q>>2]=L+J;M=L}else{M=z}B=A+12|0;G=c[B>>2]|0;N=A+16|0;if((G|0)==(c[N>>2]|0)){O=mc[c[(c[A>>2]|0)+36>>2]&63](A)|0}else{O=d[G]|0}if((Dg(O&255,v,M,q,t,u,o,r,s,m)|0)!=0){D=C;E=M;break}G=c[B>>2]|0;if((G|0)==(c[N>>2]|0)){mc[c[(c[A>>2]|0)+40>>2]&63](A)|0;n=A;z=M;continue}else{c[B>>2]=G+1;n=A;z=M;continue}}M=a[o]|0;if((M&1)==0){P=(M&255)>>>1}else{P=c[o+4>>2]|0}if((P|0)!=0?(P=c[s>>2]|0,(P-r|0)<160):0){M=c[t>>2]|0;c[s>>2]=P+4;c[P>>2]=M}b[l>>1]=Ul(E,c[q>>2]|0,k,v)|0;nj(o,r,c[s>>2]|0,k);if((A|0)!=0){if((c[A+12>>2]|0)==(c[A+16>>2]|0)?(mc[c[(c[A>>2]|0)+36>>2]&63](A)|0)==-1:0){c[g>>2]=0;Q=0}else{Q=A}}else{Q=0}A=(Q|0)==0;do{if((D|0)!=0){if((c[D+12>>2]|0)!=(c[D+16>>2]|0)){if(!A){break}c[e>>2]=Q;ze(p);ze(o);i=f;return}if((mc[c[(c[D>>2]|0)+36>>2]&63](D)|0)==-1){c[h>>2]=0;F=54;break}if(A^(D|0)==0){c[e>>2]=Q;ze(p);ze(o);i=f;return}}else{F=54}}while(0);if((F|0)==54?!A:0){c[e>>2]=Q;ze(p);ze(o);i=f;return}c[k>>2]=c[k>>2]|2;c[e>>2]=Q;ze(p);ze(o);i=f;return}function qg(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];rg(a,0,k,j,f,g,h);i=b;return}function rg(b,e,f,g,h,j,k){b=b|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0;e=i;i=i+224|0;l=e+198|0;m=e+196|0;n=e+4|0;o=e+16|0;p=e+28|0;q=e+32|0;r=e;s=e+192|0;t=c[h+4>>2]&74;if((t|0)==64){u=8}else if((t|0)==8){u=16}else if((t|0)==0){u=0}else{u=10}bh(n,h,l,m);c[o+0>>2]=0;c[o+4>>2]=0;c[o+8>>2]=0;Be(o,10,0);if((a[o]&1)==0){h=o+1|0;v=h;w=o+8|0;x=h}else{h=o+8|0;v=o+1|0;w=h;x=c[h>>2]|0}c[p>>2]=x;c[r>>2]=q;c[s>>2]=0;h=o+4|0;t=a[m]|0;m=c[f>>2]|0;y=x;a:while(1){if((m|0)!=0){if((c[m+12>>2]|0)==(c[m+16>>2]|0)?(mc[c[(c[m>>2]|0)+36>>2]&63](m)|0)==-1:0){c[f>>2]=0;z=0}else{z=m}}else{z=0}x=(z|0)==0;A=c[g>>2]|0;do{if((A|0)!=0){if((c[A+12>>2]|0)!=(c[A+16>>2]|0)){if(x){B=A;break}else{C=A;D=y;break a}}if(!((mc[c[(c[A>>2]|0)+36>>2]&63](A)|0)==-1)){if(x){B=A;break}else{C=A;D=y;break a}}else{c[g>>2]=0;E=18;break}}else{E=18}}while(0);if((E|0)==18){E=0;if(x){C=0;D=y;break}else{B=0}}A=a[o]|0;F=(A&1)==0;if(F){G=(A&255)>>>1}else{G=c[h>>2]|0}if(((c[p>>2]|0)-y|0)==(G|0)){if(F){H=(A&255)>>>1;I=(A&255)>>>1}else{A=c[h>>2]|0;H=A;I=A}Be(o,H<<1,0);if((a[o]&1)==0){J=10}else{J=(c[o>>2]&-2)+ -1|0}Be(o,J,0);if((a[o]&1)==0){K=v}else{K=c[w>>2]|0}c[p>>2]=K+I;L=K}else{L=y}A=z+12|0;F=c[A>>2]|0;M=z+16|0;if((F|0)==(c[M>>2]|0)){N=mc[c[(c[z>>2]|0)+36>>2]&63](z)|0}else{N=d[F]|0}if((Dg(N&255,u,L,p,s,t,n,q,r,l)|0)!=0){C=B;D=L;break}F=c[A>>2]|0;if((F|0)==(c[M>>2]|0)){mc[c[(c[z>>2]|0)+40>>2]&63](z)|0;m=z;y=L;continue}else{c[A>>2]=F+1;m=z;y=L;continue}}L=a[n]|0;if((L&1)==0){O=(L&255)>>>1}else{O=c[n+4>>2]|0}if((O|0)!=0?(O=c[r>>2]|0,(O-q|0)<160):0){L=c[s>>2]|0;c[r>>2]=O+4;c[O>>2]=L}c[k>>2]=Tl(D,c[p>>2]|0,j,u)|0;nj(n,q,c[r>>2]|0,j);if((z|0)!=0){if((c[z+12>>2]|0)==(c[z+16>>2]|0)?(mc[c[(c[z>>2]|0)+36>>2]&63](z)|0)==-1:0){c[f>>2]=0;P=0}else{P=z}}else{P=0}z=(P|0)==0;do{if((C|0)!=0){if((c[C+12>>2]|0)!=(c[C+16>>2]|0)){if(!z){break}c[b>>2]=P;ze(o);ze(n);i=e;return}if((mc[c[(c[C>>2]|0)+36>>2]&63](C)|0)==-1){c[g>>2]=0;E=54;break}if(z^(C|0)==0){c[b>>2]=P;ze(o);ze(n);i=e;return}}else{E=54}}while(0);if((E|0)==54?!z:0){c[b>>2]=P;ze(o);ze(n);i=e;return}c[j>>2]=c[j>>2]|2;c[b>>2]=P;ze(o);ze(n);i=e;return}function sg(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];tg(a,0,k,j,f,g,h);i=b;return}function tg(b,e,f,g,h,j,k){b=b|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0;e=i;i=i+224|0;l=e+198|0;m=e+196|0;n=e+4|0;o=e+16|0;p=e+28|0;q=e+32|0;r=e;s=e+192|0;t=c[h+4>>2]&74;if((t|0)==0){u=0}else if((t|0)==8){u=16}else if((t|0)==64){u=8}else{u=10}bh(n,h,l,m);c[o+0>>2]=0;c[o+4>>2]=0;c[o+8>>2]=0;Be(o,10,0);if((a[o]&1)==0){h=o+1|0;v=h;w=o+8|0;x=h}else{h=o+8|0;v=o+1|0;w=h;x=c[h>>2]|0}c[p>>2]=x;c[r>>2]=q;c[s>>2]=0;h=o+4|0;t=a[m]|0;m=c[f>>2]|0;y=x;a:while(1){if((m|0)!=0){if((c[m+12>>2]|0)==(c[m+16>>2]|0)?(mc[c[(c[m>>2]|0)+36>>2]&63](m)|0)==-1:0){c[f>>2]=0;z=0}else{z=m}}else{z=0}x=(z|0)==0;A=c[g>>2]|0;do{if((A|0)!=0){if((c[A+12>>2]|0)!=(c[A+16>>2]|0)){if(x){B=A;break}else{C=A;D=y;break a}}if(!((mc[c[(c[A>>2]|0)+36>>2]&63](A)|0)==-1)){if(x){B=A;break}else{C=A;D=y;break a}}else{c[g>>2]=0;E=18;break}}else{E=18}}while(0);if((E|0)==18){E=0;if(x){C=0;D=y;break}else{B=0}}A=a[o]|0;F=(A&1)==0;if(F){G=(A&255)>>>1}else{G=c[h>>2]|0}if(((c[p>>2]|0)-y|0)==(G|0)){if(F){H=(A&255)>>>1;I=(A&255)>>>1}else{A=c[h>>2]|0;H=A;I=A}Be(o,H<<1,0);if((a[o]&1)==0){J=10}else{J=(c[o>>2]&-2)+ -1|0}Be(o,J,0);if((a[o]&1)==0){K=v}else{K=c[w>>2]|0}c[p>>2]=K+I;L=K}else{L=y}A=z+12|0;F=c[A>>2]|0;M=z+16|0;if((F|0)==(c[M>>2]|0)){N=mc[c[(c[z>>2]|0)+36>>2]&63](z)|0}else{N=d[F]|0}if((Dg(N&255,u,L,p,s,t,n,q,r,l)|0)!=0){C=B;D=L;break}F=c[A>>2]|0;if((F|0)==(c[M>>2]|0)){mc[c[(c[z>>2]|0)+40>>2]&63](z)|0;m=z;y=L;continue}else{c[A>>2]=F+1;m=z;y=L;continue}}L=a[n]|0;if((L&1)==0){O=(L&255)>>>1}else{O=c[n+4>>2]|0}if((O|0)!=0?(O=c[r>>2]|0,(O-q|0)<160):0){L=c[s>>2]|0;c[r>>2]=O+4;c[O>>2]=L}c[k>>2]=Sl(D,c[p>>2]|0,j,u)|0;nj(n,q,c[r>>2]|0,j);if((z|0)!=0){if((c[z+12>>2]|0)==(c[z+16>>2]|0)?(mc[c[(c[z>>2]|0)+36>>2]&63](z)|0)==-1:0){c[f>>2]=0;P=0}else{P=z}}else{P=0}z=(P|0)==0;do{if((C|0)!=0){if((c[C+12>>2]|0)!=(c[C+16>>2]|0)){if(!z){break}c[b>>2]=P;ze(o);ze(n);i=e;return}if((mc[c[(c[C>>2]|0)+36>>2]&63](C)|0)==-1){c[g>>2]=0;E=54;break}if(z^(C|0)==0){c[b>>2]=P;ze(o);ze(n);i=e;return}}else{E=54}}while(0);if((E|0)==54?!z:0){c[b>>2]=P;ze(o);ze(n);i=e;return}c[j>>2]=c[j>>2]|2;c[b>>2]=P;ze(o);ze(n);i=e;return}function ug(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];vg(a,0,k,j,f,g,h);i=b;return}function vg(b,e,f,g,h,j,k){b=b|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0;e=i;i=i+224|0;l=e+198|0;m=e+196|0;n=e+4|0;o=e+16|0;p=e+28|0;q=e+32|0;r=e;s=e+192|0;t=c[h+4>>2]&74;if((t|0)==64){u=8}else if((t|0)==0){u=0}else if((t|0)==8){u=16}else{u=10}bh(n,h,l,m);c[o+0>>2]=0;c[o+4>>2]=0;c[o+8>>2]=0;Be(o,10,0);if((a[o]&1)==0){h=o+1|0;v=h;w=o+8|0;x=h}else{h=o+8|0;v=o+1|0;w=h;x=c[h>>2]|0}c[p>>2]=x;c[r>>2]=q;c[s>>2]=0;h=o+4|0;t=a[m]|0;m=c[f>>2]|0;y=x;a:while(1){if((m|0)!=0){if((c[m+12>>2]|0)==(c[m+16>>2]|0)?(mc[c[(c[m>>2]|0)+36>>2]&63](m)|0)==-1:0){c[f>>2]=0;z=0}else{z=m}}else{z=0}x=(z|0)==0;A=c[g>>2]|0;do{if((A|0)!=0){if((c[A+12>>2]|0)!=(c[A+16>>2]|0)){if(x){B=A;break}else{C=A;D=y;break a}}if(!((mc[c[(c[A>>2]|0)+36>>2]&63](A)|0)==-1)){if(x){B=A;break}else{C=A;D=y;break a}}else{c[g>>2]=0;E=18;break}}else{E=18}}while(0);if((E|0)==18){E=0;if(x){C=0;D=y;break}else{B=0}}A=a[o]|0;F=(A&1)==0;if(F){G=(A&255)>>>1}else{G=c[h>>2]|0}if(((c[p>>2]|0)-y|0)==(G|0)){if(F){H=(A&255)>>>1;J=(A&255)>>>1}else{A=c[h>>2]|0;H=A;J=A}Be(o,H<<1,0);if((a[o]&1)==0){K=10}else{K=(c[o>>2]&-2)+ -1|0}Be(o,K,0);if((a[o]&1)==0){L=v}else{L=c[w>>2]|0}c[p>>2]=L+J;M=L}else{M=y}A=z+12|0;F=c[A>>2]|0;N=z+16|0;if((F|0)==(c[N>>2]|0)){O=mc[c[(c[z>>2]|0)+36>>2]&63](z)|0}else{O=d[F]|0}if((Dg(O&255,u,M,p,s,t,n,q,r,l)|0)!=0){C=B;D=M;break}F=c[A>>2]|0;if((F|0)==(c[N>>2]|0)){mc[c[(c[z>>2]|0)+40>>2]&63](z)|0;m=z;y=M;continue}else{c[A>>2]=F+1;m=z;y=M;continue}}M=a[n]|0;if((M&1)==0){P=(M&255)>>>1}else{P=c[n+4>>2]|0}if((P|0)!=0?(P=c[r>>2]|0,(P-q|0)<160):0){M=c[s>>2]|0;c[r>>2]=P+4;c[P>>2]=M}M=Rl(D,c[p>>2]|0,j,u)|0;u=k;c[u>>2]=M;c[u+4>>2]=I;nj(n,q,c[r>>2]|0,j);if((z|0)!=0){if((c[z+12>>2]|0)==(c[z+16>>2]|0)?(mc[c[(c[z>>2]|0)+36>>2]&63](z)|0)==-1:0){c[f>>2]=0;Q=0}else{Q=z}}else{Q=0}z=(Q|0)==0;do{if((C|0)!=0){if((c[C+12>>2]|0)!=(c[C+16>>2]|0)){if(!z){break}c[b>>2]=Q;ze(o);ze(n);i=e;return}if((mc[c[(c[C>>2]|0)+36>>2]&63](C)|0)==-1){c[g>>2]=0;E=54;break}if(z^(C|0)==0){c[b>>2]=Q;ze(o);ze(n);i=e;return}}else{E=54}}while(0);if((E|0)==54?!z:0){c[b>>2]=Q;ze(o);ze(n);i=e;return}c[j>>2]=c[j>>2]|2;c[b>>2]=Q;ze(o);ze(n);i=e;return}function wg(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];xg(a,0,k,j,f,g,h);i=b;return}function xg(b,e,f,h,j,k,l){b=b|0;e=e|0;f=f|0;h=h|0;j=j|0;k=k|0;l=l|0;var m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0;e=i;i=i+240|0;m=e+200|0;n=e+199|0;o=e+198|0;p=e+8|0;q=e+20|0;r=e+192|0;s=e+32|0;t=e;u=e+4|0;v=e+197|0;w=e+196|0;ch(p,j,m,n,o);c[q+0>>2]=0;c[q+4>>2]=0;c[q+8>>2]=0;Be(q,10,0);if((a[q]&1)==0){j=q+1|0;x=j;y=q+8|0;z=j}else{j=q+8|0;x=q+1|0;y=j;z=c[j>>2]|0}c[r>>2]=z;c[t>>2]=s;c[u>>2]=0;a[v]=1;a[w]=69;j=q+4|0;A=a[n]|0;n=a[o]|0;o=c[f>>2]|0;B=z;a:while(1){if((o|0)!=0){if((c[o+12>>2]|0)==(c[o+16>>2]|0)?(mc[c[(c[o>>2]|0)+36>>2]&63](o)|0)==-1:0){c[f>>2]=0;C=0}else{C=o}}else{C=0}z=(C|0)==0;D=c[h>>2]|0;do{if((D|0)!=0){if((c[D+12>>2]|0)!=(c[D+16>>2]|0)){if(z){E=D;break}else{F=D;G=B;break a}}if(!((mc[c[(c[D>>2]|0)+36>>2]&63](D)|0)==-1)){if(z){E=D;break}else{F=D;G=B;break a}}else{c[h>>2]=0;H=14;break}}else{H=14}}while(0);if((H|0)==14){H=0;if(z){F=0;G=B;break}else{E=0}}D=a[q]|0;I=(D&1)==0;if(I){J=(D&255)>>>1}else{J=c[j>>2]|0}if(((c[r>>2]|0)-B|0)==(J|0)){if(I){K=(D&255)>>>1;L=(D&255)>>>1}else{D=c[j>>2]|0;K=D;L=D}Be(q,K<<1,0);if((a[q]&1)==0){M=10}else{M=(c[q>>2]&-2)+ -1|0}Be(q,M,0);if((a[q]&1)==0){N=x}else{N=c[y>>2]|0}c[r>>2]=N+L;O=N}else{O=B}D=C+12|0;I=c[D>>2]|0;P=C+16|0;if((I|0)==(c[P>>2]|0)){Q=mc[c[(c[C>>2]|0)+36>>2]&63](C)|0}else{Q=d[I]|0}if((dh(Q&255,v,w,O,r,A,n,p,s,t,u,m)|0)!=0){F=E;G=O;break}I=c[D>>2]|0;if((I|0)==(c[P>>2]|0)){mc[c[(c[C>>2]|0)+40>>2]&63](C)|0;o=C;B=O;continue}else{c[D>>2]=I+1;o=C;B=O;continue}}O=a[p]|0;if((O&1)==0){R=(O&255)>>>1}else{R=c[p+4>>2]|0}if(((R|0)!=0?(a[v]|0)!=0:0)?(v=c[t>>2]|0,(v-s|0)<160):0){R=c[u>>2]|0;c[t>>2]=v+4;c[v>>2]=R}g[l>>2]=+Ql(G,c[r>>2]|0,k);nj(p,s,c[t>>2]|0,k);if((C|0)!=0){if((c[C+12>>2]|0)==(c[C+16>>2]|0)?(mc[c[(c[C>>2]|0)+36>>2]&63](C)|0)==-1:0){c[f>>2]=0;S=0}else{S=C}}else{S=0}C=(S|0)==0;do{if((F|0)!=0){if((c[F+12>>2]|0)!=(c[F+16>>2]|0)){if(!C){break}c[b>>2]=S;ze(q);ze(p);i=e;return}if((mc[c[(c[F>>2]|0)+36>>2]&63](F)|0)==-1){c[h>>2]=0;H=51;break}if(C^(F|0)==0){c[b>>2]=S;ze(q);ze(p);i=e;return}}else{H=51}}while(0);if((H|0)==51?!C:0){c[b>>2]=S;ze(q);ze(p);i=e;return}c[k>>2]=c[k>>2]|2;c[b>>2]=S;ze(q);ze(p);i=e;return}function yg(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];zg(a,0,k,j,f,g,h);i=b;return}function zg(b,e,f,g,j,k,l){b=b|0;e=e|0;f=f|0;g=g|0;j=j|0;k=k|0;l=l|0;var m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0;e=i;i=i+240|0;m=e+200|0;n=e+199|0;o=e+198|0;p=e+8|0;q=e+20|0;r=e+192|0;s=e+32|0;t=e;u=e+4|0;v=e+197|0;w=e+196|0;ch(p,j,m,n,o);c[q+0>>2]=0;c[q+4>>2]=0;c[q+8>>2]=0;Be(q,10,0);if((a[q]&1)==0){j=q+1|0;x=j;y=q+8|0;z=j}else{j=q+8|0;x=q+1|0;y=j;z=c[j>>2]|0}c[r>>2]=z;c[t>>2]=s;c[u>>2]=0;a[v]=1;a[w]=69;j=q+4|0;A=a[n]|0;n=a[o]|0;o=c[f>>2]|0;B=z;a:while(1){if((o|0)!=0){if((c[o+12>>2]|0)==(c[o+16>>2]|0)?(mc[c[(c[o>>2]|0)+36>>2]&63](o)|0)==-1:0){c[f>>2]=0;C=0}else{C=o}}else{C=0}z=(C|0)==0;D=c[g>>2]|0;do{if((D|0)!=0){if((c[D+12>>2]|0)!=(c[D+16>>2]|0)){if(z){E=D;break}else{F=D;G=B;break a}}if(!((mc[c[(c[D>>2]|0)+36>>2]&63](D)|0)==-1)){if(z){E=D;break}else{F=D;G=B;break a}}else{c[g>>2]=0;H=14;break}}else{H=14}}while(0);if((H|0)==14){H=0;if(z){F=0;G=B;break}else{E=0}}D=a[q]|0;I=(D&1)==0;if(I){J=(D&255)>>>1}else{J=c[j>>2]|0}if(((c[r>>2]|0)-B|0)==(J|0)){if(I){K=(D&255)>>>1;L=(D&255)>>>1}else{D=c[j>>2]|0;K=D;L=D}Be(q,K<<1,0);if((a[q]&1)==0){M=10}else{M=(c[q>>2]&-2)+ -1|0}Be(q,M,0);if((a[q]&1)==0){N=x}else{N=c[y>>2]|0}c[r>>2]=N+L;O=N}else{O=B}D=C+12|0;I=c[D>>2]|0;P=C+16|0;if((I|0)==(c[P>>2]|0)){Q=mc[c[(c[C>>2]|0)+36>>2]&63](C)|0}else{Q=d[I]|0}if((dh(Q&255,v,w,O,r,A,n,p,s,t,u,m)|0)!=0){F=E;G=O;break}I=c[D>>2]|0;if((I|0)==(c[P>>2]|0)){mc[c[(c[C>>2]|0)+40>>2]&63](C)|0;o=C;B=O;continue}else{c[D>>2]=I+1;o=C;B=O;continue}}O=a[p]|0;if((O&1)==0){R=(O&255)>>>1}else{R=c[p+4>>2]|0}if(((R|0)!=0?(a[v]|0)!=0:0)?(v=c[t>>2]|0,(v-s|0)<160):0){R=c[u>>2]|0;c[t>>2]=v+4;c[v>>2]=R}h[l>>3]=+Pl(G,c[r>>2]|0,k);nj(p,s,c[t>>2]|0,k);if((C|0)!=0){if((c[C+12>>2]|0)==(c[C+16>>2]|0)?(mc[c[(c[C>>2]|0)+36>>2]&63](C)|0)==-1:0){c[f>>2]=0;S=0}else{S=C}}else{S=0}C=(S|0)==0;do{if((F|0)!=0){if((c[F+12>>2]|0)!=(c[F+16>>2]|0)){if(!C){break}c[b>>2]=S;ze(q);ze(p);i=e;return}if((mc[c[(c[F>>2]|0)+36>>2]&63](F)|0)==-1){c[g>>2]=0;H=51;break}if(C^(F|0)==0){c[b>>2]=S;ze(q);ze(p);i=e;return}}else{H=51}}while(0);if((H|0)==51?!C:0){c[b>>2]=S;ze(q);ze(p);i=e;return}c[k>>2]=c[k>>2]|2;c[b>>2]=S;ze(q);ze(p);i=e;return}function Ag(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];Bg(a,0,k,j,f,g,h);i=b;return}function Bg(b,e,f,g,j,k,l){b=b|0;e=e|0;f=f|0;g=g|0;j=j|0;k=k|0;l=l|0;var m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0;e=i;i=i+240|0;m=e+200|0;n=e+199|0;o=e+198|0;p=e+8|0;q=e+20|0;r=e+192|0;s=e+32|0;t=e;u=e+4|0;v=e+197|0;w=e+196|0;ch(p,j,m,n,o);c[q+0>>2]=0;c[q+4>>2]=0;c[q+8>>2]=0;Be(q,10,0);if((a[q]&1)==0){j=q+1|0;x=j;y=q+8|0;z=j}else{j=q+8|0;x=q+1|0;y=j;z=c[j>>2]|0}c[r>>2]=z;c[t>>2]=s;c[u>>2]=0;a[v]=1;a[w]=69;j=q+4|0;A=a[n]|0;n=a[o]|0;o=c[f>>2]|0;B=z;a:while(1){if((o|0)!=0){if((c[o+12>>2]|0)==(c[o+16>>2]|0)?(mc[c[(c[o>>2]|0)+36>>2]&63](o)|0)==-1:0){c[f>>2]=0;C=0}else{C=o}}else{C=0}z=(C|0)==0;D=c[g>>2]|0;do{if((D|0)!=0){if((c[D+12>>2]|0)!=(c[D+16>>2]|0)){if(z){E=D;break}else{F=D;G=B;break a}}if(!((mc[c[(c[D>>2]|0)+36>>2]&63](D)|0)==-1)){if(z){E=D;break}else{F=D;G=B;break a}}else{c[g>>2]=0;H=14;break}}else{H=14}}while(0);if((H|0)==14){H=0;if(z){F=0;G=B;break}else{E=0}}D=a[q]|0;I=(D&1)==0;if(I){J=(D&255)>>>1}else{J=c[j>>2]|0}if(((c[r>>2]|0)-B|0)==(J|0)){if(I){K=(D&255)>>>1;L=(D&255)>>>1}else{D=c[j>>2]|0;K=D;L=D}Be(q,K<<1,0);if((a[q]&1)==0){M=10}else{M=(c[q>>2]&-2)+ -1|0}Be(q,M,0);if((a[q]&1)==0){N=x}else{N=c[y>>2]|0}c[r>>2]=N+L;O=N}else{O=B}D=C+12|0;I=c[D>>2]|0;P=C+16|0;if((I|0)==(c[P>>2]|0)){Q=mc[c[(c[C>>2]|0)+36>>2]&63](C)|0}else{Q=d[I]|0}if((dh(Q&255,v,w,O,r,A,n,p,s,t,u,m)|0)!=0){F=E;G=O;break}I=c[D>>2]|0;if((I|0)==(c[P>>2]|0)){mc[c[(c[C>>2]|0)+40>>2]&63](C)|0;o=C;B=O;continue}else{c[D>>2]=I+1;o=C;B=O;continue}}O=a[p]|0;if((O&1)==0){R=(O&255)>>>1}else{R=c[p+4>>2]|0}if(((R|0)!=0?(a[v]|0)!=0:0)?(v=c[t>>2]|0,(v-s|0)<160):0){R=c[u>>2]|0;c[t>>2]=v+4;c[v>>2]=R}h[l>>3]=+Ol(G,c[r>>2]|0,k);nj(p,s,c[t>>2]|0,k);if((C|0)!=0){if((c[C+12>>2]|0)==(c[C+16>>2]|0)?(mc[c[(c[C>>2]|0)+36>>2]&63](C)|0)==-1:0){c[f>>2]=0;S=0}else{S=C}}else{S=0}C=(S|0)==0;do{if((F|0)!=0){if((c[F+12>>2]|0)!=(c[F+16>>2]|0)){if(!C){break}c[b>>2]=S;ze(q);ze(p);i=e;return}if((mc[c[(c[F>>2]|0)+36>>2]&63](F)|0)==-1){c[g>>2]=0;H=51;break}if(C^(F|0)==0){c[b>>2]=S;ze(q);ze(p);i=e;return}}else{H=51}}while(0);if((H|0)==51?!C:0){c[b>>2]=S;ze(q);ze(p);i=e;return}c[k>>2]=c[k>>2]|2;c[b>>2]=S;ze(q);ze(p);i=e;return}function Cg(b,e,f,g,h,j,k){b=b|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0,X=0;e=i;i=i+240|0;l=e;m=e+200|0;n=e+12|0;o=e+24|0;p=e+28|0;q=e+40|0;c[n+0>>2]=0;c[n+4>>2]=0;c[n+8>>2]=0;Te(o,h);h=c[o>>2]|0;if(!((c[2064]|0)==-1)){c[l>>2]=8256;c[l+4>>2]=112;c[l+8>>2]=0;ue(8256,l,113)}r=(c[8260>>2]|0)+ -1|0;s=c[h+8>>2]|0;if(!((c[h+12>>2]|0)-s>>2>>>0>r>>>0)){t=Wa(4)|0;qm(t);xb(t|0,16216,101)}h=c[s+(r<<2)>>2]|0;if((h|0)==0){t=Wa(4)|0;qm(t);xb(t|0,16216,101)}sc[c[(c[h>>2]|0)+32>>2]&7](h,6800,6826|0,m)|0;ae(c[o>>2]|0)|0;c[p+0>>2]=0;c[p+4>>2]=0;c[p+8>>2]=0;Be(p,10,0);if((a[p]&1)==0){o=p+1|0;u=o;v=p+8|0;w=o}else{o=p+8|0;u=p+1|0;v=o;w=c[o>>2]|0}o=p+4|0;h=m+24|0;t=m+25|0;r=q;s=m+26|0;x=m;y=n+4|0;z=c[f>>2]|0;A=q;q=0;B=w;C=w;a:while(1){if((z|0)!=0){if((c[z+12>>2]|0)==(c[z+16>>2]|0)?(mc[c[(c[z>>2]|0)+36>>2]&63](z)|0)==-1:0){c[f>>2]=0;D=0}else{D=z}}else{D=0}w=(D|0)==0;E=c[g>>2]|0;do{if((E|0)!=0){if((c[E+12>>2]|0)!=(c[E+16>>2]|0)){if(w){break}else{F=C;break a}}if(!((mc[c[(c[E>>2]|0)+36>>2]&63](E)|0)==-1)){if(w){break}else{F=C;break a}}else{c[g>>2]=0;G=19;break}}else{G=19}}while(0);if((G|0)==19?(G=0,w):0){F=C;break}E=a[p]|0;H=(E&1)==0;if(H){I=(E&255)>>>1}else{I=c[o>>2]|0}if((B-C|0)==(I|0)){if(H){J=(E&255)>>>1;K=(E&255)>>>1}else{E=c[o>>2]|0;J=E;K=E}Be(p,J<<1,0);if((a[p]&1)==0){L=10}else{L=(c[p>>2]&-2)+ -1|0}Be(p,L,0);if((a[p]&1)==0){M=u}else{M=c[v>>2]|0}N=M+K|0;O=M}else{N=B;O=C}E=c[D+12>>2]|0;if((E|0)==(c[D+16>>2]|0)){P=mc[c[(c[D>>2]|0)+36>>2]&63](D)|0}else{P=d[E]|0}E=P&255;H=(N|0)==(O|0);do{if(H){Q=(a[h]|0)==E<<24>>24;if(!Q?!((a[t]|0)==E<<24>>24):0){G=40;break}a[N]=Q?43:45;R=N+1|0;S=A;T=0}else{G=40}}while(0);do{if((G|0)==40){G=0;w=a[n]|0;if((w&1)==0){U=(w&255)>>>1}else{U=c[y>>2]|0}if((U|0)!=0&E<<24>>24==0){if((A-r|0)>=160){R=N;S=A;T=q;break}c[A>>2]=q;R=N;S=A+4|0;T=0;break}else{V=m}while(1){w=V+1|0;if((a[V]|0)==E<<24>>24){W=V;break}if((w|0)==(s|0)){W=s;break}else{V=w}}w=W-x|0;if((w|0)>23){F=O;break a}if((w|0)<22){a[N]=a[6800+w|0]|0;R=N+1|0;S=A;T=q+1|0;break}if(H){F=N;break a}if((N-O|0)>=3){F=O;break a}if((a[N+ -1|0]|0)!=48){F=O;break a}a[N]=a[6800+w|0]|0;R=N+1|0;S=A;T=0}}while(0);H=c[f>>2]|0;E=H+12|0;w=c[E>>2]|0;if((w|0)==(c[H+16>>2]|0)){mc[c[(c[H>>2]|0)+40>>2]&63](H)|0;z=H;A=S;q=T;B=R;C=O;continue}else{c[E>>2]=w+1;z=H;A=S;q=T;B=R;C=O;continue}}a[F+3|0]=0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}O=c[2036]|0;c[l>>2]=k;if((Eg(F,O,6840,l)|0)!=1){c[j>>2]=4}l=c[f>>2]|0;if((l|0)!=0){if((c[l+12>>2]|0)==(c[l+16>>2]|0)?(mc[c[(c[l>>2]|0)+36>>2]&63](l)|0)==-1:0){c[f>>2]=0;X=0}else{X=l}}else{X=0}l=(X|0)==0;f=c[g>>2]|0;do{if((f|0)!=0){if((c[f+12>>2]|0)!=(c[f+16>>2]|0)){if(!l){break}c[b>>2]=X;ze(p);ze(n);i=e;return}if((mc[c[(c[f>>2]|0)+36>>2]&63](f)|0)==-1){c[g>>2]=0;G=72;break}if(l^(f|0)==0){c[b>>2]=X;ze(p);ze(n);i=e;return}}else{G=72}}while(0);if((G|0)==72?!l:0){c[b>>2]=X;ze(p);ze(n);i=e;return}c[j>>2]=c[j>>2]|2;c[b>>2]=X;ze(p);ze(n);i=e;return}function Dg(b,d,e,f,g,h,j,k,l,m){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;m=m|0;var n=0,o=0,p=0,q=0,r=0,s=0,t=0;n=i;o=c[f>>2]|0;p=(o|0)==(e|0);do{if(p){q=(a[m+24|0]|0)==b<<24>>24;if(!q?!((a[m+25|0]|0)==b<<24>>24):0){break}c[f>>2]=e+1;a[e]=q?43:45;c[g>>2]=0;r=0;i=n;return r|0}}while(0);q=a[j]|0;if((q&1)==0){s=(q&255)>>>1}else{s=c[j+4>>2]|0}if((s|0)!=0&b<<24>>24==h<<24>>24){h=c[l>>2]|0;if((h-k|0)>=160){r=0;i=n;return r|0}k=c[g>>2]|0;c[l>>2]=h+4;c[h>>2]=k;c[g>>2]=0;r=0;i=n;return r|0}k=m+26|0;h=m;while(1){l=h+1|0;if((a[h]|0)==b<<24>>24){t=h;break}if((l|0)==(k|0)){t=k;break}else{h=l}}h=t-m|0;if((h|0)>23){r=-1;i=n;return r|0}if((d|0)==10|(d|0)==8){if((h|0)>=(d|0)){r=-1;i=n;return r|0}}else if((d|0)==16?(h|0)>=22:0){if(p){r=-1;i=n;return r|0}if((o-e|0)>=3){r=-1;i=n;return r|0}if((a[o+ -1|0]|0)!=48){r=-1;i=n;return r|0}c[g>>2]=0;e=a[6800+h|0]|0;c[f>>2]=o+1;a[o]=e;r=0;i=n;return r|0}e=a[6800+h|0]|0;c[f>>2]=o+1;a[o]=e;c[g>>2]=(c[g>>2]|0)+1;r=0;i=n;return r|0}function Eg(a,b,d,e){a=a|0;b=b|0;d=d|0;e=e|0;var f=0,g=0;f=i;i=i+16|0;g=f;c[g>>2]=e;e=jb(b|0)|0;b=Ka(a|0,d|0,g|0)|0;if((e|0)==0){i=f;return b|0}jb(e|0)|0;i=f;return b|0}function Fg(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function Gg(a){a=a|0;return}



function Hg(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0;k=i;i=i+80|0;l=k;m=k+64|0;n=k+60|0;o=k+56|0;p=k+52|0;q=k+68|0;r=k+16|0;s=k+12|0;t=k+24|0;u=k+48|0;if((c[g+4>>2]&1|0)==0){c[n>>2]=-1;v=c[(c[d>>2]|0)+16>>2]|0;c[p>>2]=c[e>>2];c[q>>2]=c[f>>2];c[m+0>>2]=c[p+0>>2];c[l+0>>2]=c[q+0>>2];hc[v&63](o,d,m,l,g,h,n);m=c[o>>2]|0;c[e>>2]=m;o=c[n>>2]|0;if((o|0)==1){a[j]=1}else if((o|0)==0){a[j]=0}else{a[j]=1;c[h>>2]=4}c[b>>2]=m;i=k;return}Te(r,g);m=c[r>>2]|0;if(!((c[2062]|0)==-1)){c[l>>2]=8248;c[l+4>>2]=112;c[l+8>>2]=0;ue(8248,l,113)}o=(c[8252>>2]|0)+ -1|0;n=c[m+8>>2]|0;if(!((c[m+12>>2]|0)-n>>2>>>0>o>>>0)){w=Wa(4)|0;qm(w);xb(w|0,16216,101)}m=c[n+(o<<2)>>2]|0;if((m|0)==0){w=Wa(4)|0;qm(w);xb(w|0,16216,101)}ae(c[r>>2]|0)|0;Te(s,g);g=c[s>>2]|0;if(!((c[2102]|0)==-1)){c[l>>2]=8408;c[l+4>>2]=112;c[l+8>>2]=0;ue(8408,l,113)}r=(c[8412>>2]|0)+ -1|0;w=c[g+8>>2]|0;if(!((c[g+12>>2]|0)-w>>2>>>0>r>>>0)){x=Wa(4)|0;qm(x);xb(x|0,16216,101)}g=c[w+(r<<2)>>2]|0;if((g|0)==0){x=Wa(4)|0;qm(x);xb(x|0,16216,101)}ae(c[s>>2]|0)|0;kc[c[(c[g>>2]|0)+24>>2]&63](t,g);kc[c[(c[g>>2]|0)+28>>2]&63](t+12|0,g);c[u>>2]=c[f>>2];f=t+24|0;c[l+0>>2]=c[u+0>>2];a[j]=(Ig(e,l,t,f,m,h,1)|0)==(t|0)|0;c[b>>2]=c[e>>2];Ke(t+12|0);Ke(t);i=k;return}function Ig(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0,X=0,Y=0,Z=0,_=0,$=0,aa=0,ba=0;k=i;i=i+112|0;l=k;m=(f-e|0)/12|0;if(m>>>0>100){n=Lm(m)|0;if((n|0)==0){_m()}else{o=n;p=n}}else{o=0;p=l}l=(e|0)==(f|0);if(l){q=0;r=m}else{n=e;s=0;t=m;m=p;while(1){u=a[n]|0;if((u&1)==0){v=(u&255)>>>1}else{v=c[n+4>>2]|0}if((v|0)==0){a[m]=2;w=s+1|0;x=t+ -1|0}else{a[m]=1;w=s;x=t}u=n+12|0;if((u|0)==(f|0)){q=w;r=x;break}else{n=u;s=w;t=x;m=m+1|0}}}m=0;x=q;q=r;a:while(1){r=c[b>>2]|0;do{if((r|0)!=0){t=c[r+12>>2]|0;if((t|0)==(c[r+16>>2]|0)){y=mc[c[(c[r>>2]|0)+36>>2]&63](r)|0}else{y=c[t>>2]|0}if((y|0)==-1){c[b>>2]=0;z=1;break}else{z=(c[b>>2]|0)==0;break}}else{z=1}}while(0);r=c[d>>2]|0;if((r|0)!=0){t=c[r+12>>2]|0;if((t|0)==(c[r+16>>2]|0)){A=mc[c[(c[r>>2]|0)+36>>2]&63](r)|0}else{A=c[t>>2]|0}if((A|0)==-1){c[d>>2]=0;B=0;C=1}else{B=r;C=0}}else{B=0;C=1}D=c[b>>2]|0;if(!((z^C)&(q|0)!=0)){break}r=c[D+12>>2]|0;if((r|0)==(c[D+16>>2]|0)){E=mc[c[(c[D>>2]|0)+36>>2]&63](D)|0}else{E=c[r>>2]|0}if(j){F=E}else{F=vc[c[(c[g>>2]|0)+28>>2]&15](g,E)|0}r=m+1|0;if(l){m=r;continue}b:do{if(j){t=0;w=e;s=x;n=q;v=p;while(1){do{if((a[v]|0)==1){u=a[w]|0;G=(u&1)==0;if(G){H=w+4|0}else{H=c[w+8>>2]|0}if((F|0)!=(c[H+(m<<2)>>2]|0)){a[v]=0;I=t;J=s;K=n+ -1|0;break}if(G){L=(u&255)>>>1}else{L=c[w+4>>2]|0}if((L|0)==(r|0)){a[v]=2;I=1;J=s+1|0;K=n+ -1|0}else{I=1;J=s;K=n}}else{I=t;J=s;K=n}}while(0);u=w+12|0;if((u|0)==(f|0)){M=I;N=J;O=K;break b}t=I;w=u;s=J;n=K;v=v+1|0}}else{v=0;n=e;s=x;w=q;t=p;while(1){do{if((a[t]|0)==1){if((a[n]&1)==0){P=n+4|0}else{P=c[n+8>>2]|0}if((F|0)!=(vc[c[(c[g>>2]|0)+28>>2]&15](g,c[P+(m<<2)>>2]|0)|0)){a[t]=0;Q=v;R=s;S=w+ -1|0;break}u=a[n]|0;if((u&1)==0){T=(u&255)>>>1}else{T=c[n+4>>2]|0}if((T|0)==(r|0)){a[t]=2;Q=1;R=s+1|0;S=w+ -1|0}else{Q=1;R=s;S=w}}else{Q=v;R=s;S=w}}while(0);u=n+12|0;if((u|0)==(f|0)){M=Q;N=R;O=S;break b}v=Q;n=u;s=R;w=S;t=t+1|0}}}while(0);if(!M){m=r;x=N;q=O;continue}t=c[b>>2]|0;w=t+12|0;s=c[w>>2]|0;if((s|0)==(c[t+16>>2]|0)){mc[c[(c[t>>2]|0)+40>>2]&63](t)|0}else{c[w>>2]=s+4}if((O+N|0)>>>0<2){m=r;x=N;q=O;continue}else{U=e;V=N;W=p}while(1){if((a[W]|0)==2){s=a[U]|0;if((s&1)==0){X=(s&255)>>>1}else{X=c[U+4>>2]|0}if((X|0)!=(r|0)){a[W]=0;Y=V+ -1|0}else{Y=V}}else{Y=V}s=U+12|0;if((s|0)==(f|0)){m=r;x=Y;q=O;continue a}else{U=s;V=Y;W=W+1|0}}}do{if((D|0)!=0){W=c[D+12>>2]|0;if((W|0)==(c[D+16>>2]|0)){Z=mc[c[(c[D>>2]|0)+36>>2]&63](D)|0}else{Z=c[W>>2]|0}if((Z|0)==-1){c[b>>2]=0;_=1;break}else{_=(c[b>>2]|0)==0;break}}else{_=1}}while(0);do{if((B|0)!=0){b=c[B+12>>2]|0;if((b|0)==(c[B+16>>2]|0)){$=mc[c[(c[B>>2]|0)+36>>2]&63](B)|0}else{$=c[b>>2]|0}if(!(($|0)==-1)){if(_){break}else{aa=87;break}}else{c[d>>2]=0;aa=85;break}}else{aa=85}}while(0);if((aa|0)==85?_:0){aa=87}if((aa|0)==87){c[h>>2]=c[h>>2]|2}c:do{if(!l){if((a[p]|0)==2){ba=e}else{_=e;d=p;while(1){$=_+12|0;B=d+1|0;if(($|0)==(f|0)){aa=92;break c}if((a[B]|0)==2){ba=$;break}else{_=$;d=B}}}}else{aa=92}}while(0);if((aa|0)==92){c[h>>2]=c[h>>2]|4;ba=f}if((o|0)==0){i=k;return ba|0}Mm(o);i=k;return ba|0}function Jg(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];Kg(a,0,k,j,f,g,h);i=b;return}function Kg(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0;d=i;i=i+304|0;k=d+200|0;l=d;m=d+4|0;n=d+16|0;o=d+28|0;p=d+32|0;q=d+192|0;r=d+196|0;s=c[g+4>>2]&74;if((s|0)==0){t=0}else if((s|0)==8){t=16}else if((s|0)==64){t=8}else{t=10}eh(m,g,k,l);c[n+0>>2]=0;c[n+4>>2]=0;c[n+8>>2]=0;Be(n,10,0);if((a[n]&1)==0){g=n+1|0;u=g;v=n+8|0;w=g}else{g=n+8|0;u=n+1|0;v=g;w=c[g>>2]|0}c[o>>2]=w;c[q>>2]=p;c[r>>2]=0;g=n+4|0;s=c[l>>2]|0;l=c[e>>2]|0;x=w;a:while(1){if((l|0)!=0){w=c[l+12>>2]|0;if((w|0)==(c[l+16>>2]|0)){y=mc[c[(c[l>>2]|0)+36>>2]&63](l)|0}else{y=c[w>>2]|0}if((y|0)==-1){c[e>>2]=0;z=1;A=0}else{z=0;A=l}}else{z=1;A=0}w=c[f>>2]|0;do{if((w|0)!=0){B=c[w+12>>2]|0;if((B|0)==(c[w+16>>2]|0)){C=mc[c[(c[w>>2]|0)+36>>2]&63](w)|0}else{C=c[B>>2]|0}if(!((C|0)==-1)){if(z){D=w;break}else{E=w;F=x;break a}}else{c[f>>2]=0;G=21;break}}else{G=21}}while(0);if((G|0)==21){G=0;if(z){E=0;F=x;break}else{D=0}}w=a[n]|0;B=(w&1)==0;if(B){H=(w&255)>>>1}else{H=c[g>>2]|0}if(((c[o>>2]|0)-x|0)==(H|0)){if(B){I=(w&255)>>>1;J=(w&255)>>>1}else{w=c[g>>2]|0;I=w;J=w}Be(n,I<<1,0);if((a[n]&1)==0){K=10}else{K=(c[n>>2]&-2)+ -1|0}Be(n,K,0);if((a[n]&1)==0){L=u}else{L=c[v>>2]|0}c[o>>2]=L+J;M=L}else{M=x}w=A+12|0;B=c[w>>2]|0;N=A+16|0;if((B|0)==(c[N>>2]|0)){O=mc[c[(c[A>>2]|0)+36>>2]&63](A)|0}else{O=c[B>>2]|0}if((ah(O,t,M,o,r,s,m,p,q,k)|0)!=0){E=D;F=M;break}B=c[w>>2]|0;if((B|0)==(c[N>>2]|0)){mc[c[(c[A>>2]|0)+40>>2]&63](A)|0;l=A;x=M;continue}else{c[w>>2]=B+4;l=A;x=M;continue}}M=a[m]|0;if((M&1)==0){P=(M&255)>>>1}else{P=c[m+4>>2]|0}if((P|0)!=0?(P=c[q>>2]|0,(P-p|0)<160):0){M=c[r>>2]|0;c[q>>2]=P+4;c[P>>2]=M}c[j>>2]=Wl(F,c[o>>2]|0,h,t)|0;nj(m,p,c[q>>2]|0,h);if((A|0)!=0){q=c[A+12>>2]|0;if((q|0)==(c[A+16>>2]|0)){Q=mc[c[(c[A>>2]|0)+36>>2]&63](A)|0}else{Q=c[q>>2]|0}if((Q|0)==-1){c[e>>2]=0;R=0;S=1}else{R=A;S=0}}else{R=0;S=1}do{if((E|0)!=0){A=c[E+12>>2]|0;if((A|0)==(c[E+16>>2]|0)){T=mc[c[(c[E>>2]|0)+36>>2]&63](E)|0}else{T=c[A>>2]|0}if((T|0)==-1){c[f>>2]=0;G=60;break}if(S){c[b>>2]=R;ze(n);ze(m);i=d;return}}else{G=60}}while(0);if((G|0)==60?!S:0){c[b>>2]=R;ze(n);ze(m);i=d;return}c[h>>2]=c[h>>2]|2;c[b>>2]=R;ze(n);ze(m);i=d;return}function Lg(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];Mg(a,0,k,j,f,g,h);i=b;return}function Mg(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0;d=i;i=i+304|0;k=d+200|0;l=d;m=d+4|0;n=d+16|0;o=d+28|0;p=d+32|0;q=d+192|0;r=d+196|0;s=c[g+4>>2]&74;if((s|0)==64){t=8}else if((s|0)==0){t=0}else if((s|0)==8){t=16}else{t=10}eh(m,g,k,l);c[n+0>>2]=0;c[n+4>>2]=0;c[n+8>>2]=0;Be(n,10,0);if((a[n]&1)==0){g=n+1|0;u=g;v=n+8|0;w=g}else{g=n+8|0;u=n+1|0;v=g;w=c[g>>2]|0}c[o>>2]=w;c[q>>2]=p;c[r>>2]=0;g=n+4|0;s=c[l>>2]|0;l=c[e>>2]|0;x=w;a:while(1){if((l|0)!=0){w=c[l+12>>2]|0;if((w|0)==(c[l+16>>2]|0)){y=mc[c[(c[l>>2]|0)+36>>2]&63](l)|0}else{y=c[w>>2]|0}if((y|0)==-1){c[e>>2]=0;z=1;A=0}else{z=0;A=l}}else{z=1;A=0}w=c[f>>2]|0;do{if((w|0)!=0){B=c[w+12>>2]|0;if((B|0)==(c[w+16>>2]|0)){C=mc[c[(c[w>>2]|0)+36>>2]&63](w)|0}else{C=c[B>>2]|0}if(!((C|0)==-1)){if(z){D=w;break}else{E=w;F=x;break a}}else{c[f>>2]=0;G=21;break}}else{G=21}}while(0);if((G|0)==21){G=0;if(z){E=0;F=x;break}else{D=0}}w=a[n]|0;B=(w&1)==0;if(B){H=(w&255)>>>1}else{H=c[g>>2]|0}if(((c[o>>2]|0)-x|0)==(H|0)){if(B){J=(w&255)>>>1;K=(w&255)>>>1}else{w=c[g>>2]|0;J=w;K=w}Be(n,J<<1,0);if((a[n]&1)==0){L=10}else{L=(c[n>>2]&-2)+ -1|0}Be(n,L,0);if((a[n]&1)==0){M=u}else{M=c[v>>2]|0}c[o>>2]=M+K;N=M}else{N=x}w=A+12|0;B=c[w>>2]|0;O=A+16|0;if((B|0)==(c[O>>2]|0)){P=mc[c[(c[A>>2]|0)+36>>2]&63](A)|0}else{P=c[B>>2]|0}if((ah(P,t,N,o,r,s,m,p,q,k)|0)!=0){E=D;F=N;break}B=c[w>>2]|0;if((B|0)==(c[O>>2]|0)){mc[c[(c[A>>2]|0)+40>>2]&63](A)|0;l=A;x=N;continue}else{c[w>>2]=B+4;l=A;x=N;continue}}N=a[m]|0;if((N&1)==0){Q=(N&255)>>>1}else{Q=c[m+4>>2]|0}if((Q|0)!=0?(Q=c[q>>2]|0,(Q-p|0)<160):0){N=c[r>>2]|0;c[q>>2]=Q+4;c[Q>>2]=N}N=Vl(F,c[o>>2]|0,h,t)|0;t=j;c[t>>2]=N;c[t+4>>2]=I;nj(m,p,c[q>>2]|0,h);if((A|0)!=0){q=c[A+12>>2]|0;if((q|0)==(c[A+16>>2]|0)){R=mc[c[(c[A>>2]|0)+36>>2]&63](A)|0}else{R=c[q>>2]|0}if((R|0)==-1){c[e>>2]=0;S=0;T=1}else{S=A;T=0}}else{S=0;T=1}do{if((E|0)!=0){A=c[E+12>>2]|0;if((A|0)==(c[E+16>>2]|0)){U=mc[c[(c[E>>2]|0)+36>>2]&63](E)|0}else{U=c[A>>2]|0}if((U|0)==-1){c[f>>2]=0;G=60;break}if(T){c[b>>2]=S;ze(n);ze(m);i=d;return}}else{G=60}}while(0);if((G|0)==60?!T:0){c[b>>2]=S;ze(n);ze(m);i=d;return}c[h>>2]=c[h>>2]|2;c[b>>2]=S;ze(n);ze(m);i=d;return}function Ng(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];Og(a,0,k,j,f,g,h);i=b;return}function Og(d,e,f,g,h,j,k){d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0;e=i;i=i+304|0;l=e+200|0;m=e;n=e+4|0;o=e+16|0;p=e+28|0;q=e+32|0;r=e+192|0;s=e+196|0;t=c[h+4>>2]&74;if((t|0)==0){u=0}else if((t|0)==64){u=8}else if((t|0)==8){u=16}else{u=10}eh(n,h,l,m);c[o+0>>2]=0;c[o+4>>2]=0;c[o+8>>2]=0;Be(o,10,0);if((a[o]&1)==0){h=o+1|0;v=h;w=o+8|0;x=h}else{h=o+8|0;v=o+1|0;w=h;x=c[h>>2]|0}c[p>>2]=x;c[r>>2]=q;c[s>>2]=0;h=o+4|0;t=c[m>>2]|0;m=c[f>>2]|0;y=x;a:while(1){if((m|0)!=0){x=c[m+12>>2]|0;if((x|0)==(c[m+16>>2]|0)){z=mc[c[(c[m>>2]|0)+36>>2]&63](m)|0}else{z=c[x>>2]|0}if((z|0)==-1){c[f>>2]=0;A=1;B=0}else{A=0;B=m}}else{A=1;B=0}x=c[g>>2]|0;do{if((x|0)!=0){C=c[x+12>>2]|0;if((C|0)==(c[x+16>>2]|0)){D=mc[c[(c[x>>2]|0)+36>>2]&63](x)|0}else{D=c[C>>2]|0}if(!((D|0)==-1)){if(A){E=x;break}else{F=x;G=y;break a}}else{c[g>>2]=0;H=21;break}}else{H=21}}while(0);if((H|0)==21){H=0;if(A){F=0;G=y;break}else{E=0}}x=a[o]|0;C=(x&1)==0;if(C){I=(x&255)>>>1}else{I=c[h>>2]|0}if(((c[p>>2]|0)-y|0)==(I|0)){if(C){J=(x&255)>>>1;K=(x&255)>>>1}else{x=c[h>>2]|0;J=x;K=x}Be(o,J<<1,0);if((a[o]&1)==0){L=10}else{L=(c[o>>2]&-2)+ -1|0}Be(o,L,0);if((a[o]&1)==0){M=v}else{M=c[w>>2]|0}c[p>>2]=M+K;N=M}else{N=y}x=B+12|0;C=c[x>>2]|0;O=B+16|0;if((C|0)==(c[O>>2]|0)){P=mc[c[(c[B>>2]|0)+36>>2]&63](B)|0}else{P=c[C>>2]|0}if((ah(P,u,N,p,s,t,n,q,r,l)|0)!=0){F=E;G=N;break}C=c[x>>2]|0;if((C|0)==(c[O>>2]|0)){mc[c[(c[B>>2]|0)+40>>2]&63](B)|0;m=B;y=N;continue}else{c[x>>2]=C+4;m=B;y=N;continue}}N=a[n]|0;if((N&1)==0){Q=(N&255)>>>1}else{Q=c[n+4>>2]|0}if((Q|0)!=0?(Q=c[r>>2]|0,(Q-q|0)<160):0){N=c[s>>2]|0;c[r>>2]=Q+4;c[Q>>2]=N}b[k>>1]=Ul(G,c[p>>2]|0,j,u)|0;nj(n,q,c[r>>2]|0,j);if((B|0)!=0){r=c[B+12>>2]|0;if((r|0)==(c[B+16>>2]|0)){R=mc[c[(c[B>>2]|0)+36>>2]&63](B)|0}else{R=c[r>>2]|0}if((R|0)==-1){c[f>>2]=0;S=0;T=1}else{S=B;T=0}}else{S=0;T=1}do{if((F|0)!=0){B=c[F+12>>2]|0;if((B|0)==(c[F+16>>2]|0)){U=mc[c[(c[F>>2]|0)+36>>2]&63](F)|0}else{U=c[B>>2]|0}if((U|0)==-1){c[g>>2]=0;H=60;break}if(T){c[d>>2]=S;ze(o);ze(n);i=e;return}}else{H=60}}while(0);if((H|0)==60?!T:0){c[d>>2]=S;ze(o);ze(n);i=e;return}c[j>>2]=c[j>>2]|2;c[d>>2]=S;ze(o);ze(n);i=e;return}function Pg(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];Qg(a,0,k,j,f,g,h);i=b;return}function Qg(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0;d=i;i=i+304|0;k=d+200|0;l=d;m=d+4|0;n=d+16|0;o=d+28|0;p=d+32|0;q=d+192|0;r=d+196|0;s=c[g+4>>2]&74;if((s|0)==64){t=8}else if((s|0)==0){t=0}else if((s|0)==8){t=16}else{t=10}eh(m,g,k,l);c[n+0>>2]=0;c[n+4>>2]=0;c[n+8>>2]=0;Be(n,10,0);if((a[n]&1)==0){g=n+1|0;u=g;v=n+8|0;w=g}else{g=n+8|0;u=n+1|0;v=g;w=c[g>>2]|0}c[o>>2]=w;c[q>>2]=p;c[r>>2]=0;g=n+4|0;s=c[l>>2]|0;l=c[e>>2]|0;x=w;a:while(1){if((l|0)!=0){w=c[l+12>>2]|0;if((w|0)==(c[l+16>>2]|0)){y=mc[c[(c[l>>2]|0)+36>>2]&63](l)|0}else{y=c[w>>2]|0}if((y|0)==-1){c[e>>2]=0;z=1;A=0}else{z=0;A=l}}else{z=1;A=0}w=c[f>>2]|0;do{if((w|0)!=0){B=c[w+12>>2]|0;if((B|0)==(c[w+16>>2]|0)){C=mc[c[(c[w>>2]|0)+36>>2]&63](w)|0}else{C=c[B>>2]|0}if(!((C|0)==-1)){if(z){D=w;break}else{E=w;F=x;break a}}else{c[f>>2]=0;G=21;break}}else{G=21}}while(0);if((G|0)==21){G=0;if(z){E=0;F=x;break}else{D=0}}w=a[n]|0;B=(w&1)==0;if(B){H=(w&255)>>>1}else{H=c[g>>2]|0}if(((c[o>>2]|0)-x|0)==(H|0)){if(B){I=(w&255)>>>1;J=(w&255)>>>1}else{w=c[g>>2]|0;I=w;J=w}Be(n,I<<1,0);if((a[n]&1)==0){K=10}else{K=(c[n>>2]&-2)+ -1|0}Be(n,K,0);if((a[n]&1)==0){L=u}else{L=c[v>>2]|0}c[o>>2]=L+J;M=L}else{M=x}w=A+12|0;B=c[w>>2]|0;N=A+16|0;if((B|0)==(c[N>>2]|0)){O=mc[c[(c[A>>2]|0)+36>>2]&63](A)|0}else{O=c[B>>2]|0}if((ah(O,t,M,o,r,s,m,p,q,k)|0)!=0){E=D;F=M;break}B=c[w>>2]|0;if((B|0)==(c[N>>2]|0)){mc[c[(c[A>>2]|0)+40>>2]&63](A)|0;l=A;x=M;continue}else{c[w>>2]=B+4;l=A;x=M;continue}}M=a[m]|0;if((M&1)==0){P=(M&255)>>>1}else{P=c[m+4>>2]|0}if((P|0)!=0?(P=c[q>>2]|0,(P-p|0)<160):0){M=c[r>>2]|0;c[q>>2]=P+4;c[P>>2]=M}c[j>>2]=Tl(F,c[o>>2]|0,h,t)|0;nj(m,p,c[q>>2]|0,h);if((A|0)!=0){q=c[A+12>>2]|0;if((q|0)==(c[A+16>>2]|0)){Q=mc[c[(c[A>>2]|0)+36>>2]&63](A)|0}else{Q=c[q>>2]|0}if((Q|0)==-1){c[e>>2]=0;R=0;S=1}else{R=A;S=0}}else{R=0;S=1}do{if((E|0)!=0){A=c[E+12>>2]|0;if((A|0)==(c[E+16>>2]|0)){T=mc[c[(c[E>>2]|0)+36>>2]&63](E)|0}else{T=c[A>>2]|0}if((T|0)==-1){c[f>>2]=0;G=60;break}if(S){c[b>>2]=R;ze(n);ze(m);i=d;return}}else{G=60}}while(0);if((G|0)==60?!S:0){c[b>>2]=R;ze(n);ze(m);i=d;return}c[h>>2]=c[h>>2]|2;c[b>>2]=R;ze(n);ze(m);i=d;return}function Rg(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];Sg(a,0,k,j,f,g,h);i=b;return}function Sg(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0;d=i;i=i+304|0;k=d+200|0;l=d;m=d+4|0;n=d+16|0;o=d+28|0;p=d+32|0;q=d+192|0;r=d+196|0;s=c[g+4>>2]&74;if((s|0)==0){t=0}else if((s|0)==8){t=16}else if((s|0)==64){t=8}else{t=10}eh(m,g,k,l);c[n+0>>2]=0;c[n+4>>2]=0;c[n+8>>2]=0;Be(n,10,0);if((a[n]&1)==0){g=n+1|0;u=g;v=n+8|0;w=g}else{g=n+8|0;u=n+1|0;v=g;w=c[g>>2]|0}c[o>>2]=w;c[q>>2]=p;c[r>>2]=0;g=n+4|0;s=c[l>>2]|0;l=c[e>>2]|0;x=w;a:while(1){if((l|0)!=0){w=c[l+12>>2]|0;if((w|0)==(c[l+16>>2]|0)){y=mc[c[(c[l>>2]|0)+36>>2]&63](l)|0}else{y=c[w>>2]|0}if((y|0)==-1){c[e>>2]=0;z=1;A=0}else{z=0;A=l}}else{z=1;A=0}w=c[f>>2]|0;do{if((w|0)!=0){B=c[w+12>>2]|0;if((B|0)==(c[w+16>>2]|0)){C=mc[c[(c[w>>2]|0)+36>>2]&63](w)|0}else{C=c[B>>2]|0}if(!((C|0)==-1)){if(z){D=w;break}else{E=w;F=x;break a}}else{c[f>>2]=0;G=21;break}}else{G=21}}while(0);if((G|0)==21){G=0;if(z){E=0;F=x;break}else{D=0}}w=a[n]|0;B=(w&1)==0;if(B){H=(w&255)>>>1}else{H=c[g>>2]|0}if(((c[o>>2]|0)-x|0)==(H|0)){if(B){I=(w&255)>>>1;J=(w&255)>>>1}else{w=c[g>>2]|0;I=w;J=w}Be(n,I<<1,0);if((a[n]&1)==0){K=10}else{K=(c[n>>2]&-2)+ -1|0}Be(n,K,0);if((a[n]&1)==0){L=u}else{L=c[v>>2]|0}c[o>>2]=L+J;M=L}else{M=x}w=A+12|0;B=c[w>>2]|0;N=A+16|0;if((B|0)==(c[N>>2]|0)){O=mc[c[(c[A>>2]|0)+36>>2]&63](A)|0}else{O=c[B>>2]|0}if((ah(O,t,M,o,r,s,m,p,q,k)|0)!=0){E=D;F=M;break}B=c[w>>2]|0;if((B|0)==(c[N>>2]|0)){mc[c[(c[A>>2]|0)+40>>2]&63](A)|0;l=A;x=M;continue}else{c[w>>2]=B+4;l=A;x=M;continue}}M=a[m]|0;if((M&1)==0){P=(M&255)>>>1}else{P=c[m+4>>2]|0}if((P|0)!=0?(P=c[q>>2]|0,(P-p|0)<160):0){M=c[r>>2]|0;c[q>>2]=P+4;c[P>>2]=M}c[j>>2]=Sl(F,c[o>>2]|0,h,t)|0;nj(m,p,c[q>>2]|0,h);if((A|0)!=0){q=c[A+12>>2]|0;if((q|0)==(c[A+16>>2]|0)){Q=mc[c[(c[A>>2]|0)+36>>2]&63](A)|0}else{Q=c[q>>2]|0}if((Q|0)==-1){c[e>>2]=0;R=0;S=1}else{R=A;S=0}}else{R=0;S=1}do{if((E|0)!=0){A=c[E+12>>2]|0;if((A|0)==(c[E+16>>2]|0)){T=mc[c[(c[E>>2]|0)+36>>2]&63](E)|0}else{T=c[A>>2]|0}if((T|0)==-1){c[f>>2]=0;G=60;break}if(S){c[b>>2]=R;ze(n);ze(m);i=d;return}}else{G=60}}while(0);if((G|0)==60?!S:0){c[b>>2]=R;ze(n);ze(m);i=d;return}c[h>>2]=c[h>>2]|2;c[b>>2]=R;ze(n);ze(m);i=d;return}function Tg(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];Ug(a,0,k,j,f,g,h);i=b;return}function Ug(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0;d=i;i=i+304|0;k=d+200|0;l=d;m=d+4|0;n=d+16|0;o=d+28|0;p=d+32|0;q=d+192|0;r=d+196|0;s=c[g+4>>2]&74;if((s|0)==8){t=16}else if((s|0)==0){t=0}else if((s|0)==64){t=8}else{t=10}eh(m,g,k,l);c[n+0>>2]=0;c[n+4>>2]=0;c[n+8>>2]=0;Be(n,10,0);if((a[n]&1)==0){g=n+1|0;u=g;v=n+8|0;w=g}else{g=n+8|0;u=n+1|0;v=g;w=c[g>>2]|0}c[o>>2]=w;c[q>>2]=p;c[r>>2]=0;g=n+4|0;s=c[l>>2]|0;l=c[e>>2]|0;x=w;a:while(1){if((l|0)!=0){w=c[l+12>>2]|0;if((w|0)==(c[l+16>>2]|0)){y=mc[c[(c[l>>2]|0)+36>>2]&63](l)|0}else{y=c[w>>2]|0}if((y|0)==-1){c[e>>2]=0;z=1;A=0}else{z=0;A=l}}else{z=1;A=0}w=c[f>>2]|0;do{if((w|0)!=0){B=c[w+12>>2]|0;if((B|0)==(c[w+16>>2]|0)){C=mc[c[(c[w>>2]|0)+36>>2]&63](w)|0}else{C=c[B>>2]|0}if(!((C|0)==-1)){if(z){D=w;break}else{E=w;F=x;break a}}else{c[f>>2]=0;G=21;break}}else{G=21}}while(0);if((G|0)==21){G=0;if(z){E=0;F=x;break}else{D=0}}w=a[n]|0;B=(w&1)==0;if(B){H=(w&255)>>>1}else{H=c[g>>2]|0}if(((c[o>>2]|0)-x|0)==(H|0)){if(B){J=(w&255)>>>1;K=(w&255)>>>1}else{w=c[g>>2]|0;J=w;K=w}Be(n,J<<1,0);if((a[n]&1)==0){L=10}else{L=(c[n>>2]&-2)+ -1|0}Be(n,L,0);if((a[n]&1)==0){M=u}else{M=c[v>>2]|0}c[o>>2]=M+K;N=M}else{N=x}w=A+12|0;B=c[w>>2]|0;O=A+16|0;if((B|0)==(c[O>>2]|0)){P=mc[c[(c[A>>2]|0)+36>>2]&63](A)|0}else{P=c[B>>2]|0}if((ah(P,t,N,o,r,s,m,p,q,k)|0)!=0){E=D;F=N;break}B=c[w>>2]|0;if((B|0)==(c[O>>2]|0)){mc[c[(c[A>>2]|0)+40>>2]&63](A)|0;l=A;x=N;continue}else{c[w>>2]=B+4;l=A;x=N;continue}}N=a[m]|0;if((N&1)==0){Q=(N&255)>>>1}else{Q=c[m+4>>2]|0}if((Q|0)!=0?(Q=c[q>>2]|0,(Q-p|0)<160):0){N=c[r>>2]|0;c[q>>2]=Q+4;c[Q>>2]=N}N=Rl(F,c[o>>2]|0,h,t)|0;t=j;c[t>>2]=N;c[t+4>>2]=I;nj(m,p,c[q>>2]|0,h);if((A|0)!=0){q=c[A+12>>2]|0;if((q|0)==(c[A+16>>2]|0)){R=mc[c[(c[A>>2]|0)+36>>2]&63](A)|0}else{R=c[q>>2]|0}if((R|0)==-1){c[e>>2]=0;S=0;T=1}else{S=A;T=0}}else{S=0;T=1}do{if((E|0)!=0){A=c[E+12>>2]|0;if((A|0)==(c[E+16>>2]|0)){U=mc[c[(c[E>>2]|0)+36>>2]&63](E)|0}else{U=c[A>>2]|0}if((U|0)==-1){c[f>>2]=0;G=60;break}if(T){c[b>>2]=S;ze(n);ze(m);i=d;return}}else{G=60}}while(0);if((G|0)==60?!T:0){c[b>>2]=S;ze(n);ze(m);i=d;return}c[h>>2]=c[h>>2]|2;c[b>>2]=S;ze(n);ze(m);i=d;return}function Vg(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];Wg(a,0,k,j,f,g,h);i=b;return}function Wg(b,d,e,f,h,j,k){b=b|0;d=d|0;e=e|0;f=f|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0;d=i;i=i+352|0;l=d+176|0;m=d+332|0;n=d+328|0;o=d+316|0;p=d+304|0;q=d+168|0;r=d+8|0;s=d+4|0;t=d;u=d+337|0;v=d+336|0;fh(o,h,l,m,n);c[p+0>>2]=0;c[p+4>>2]=0;c[p+8>>2]=0;Be(p,10,0);if((a[p]&1)==0){h=p+1|0;w=h;x=p+8|0;y=h}else{h=p+8|0;w=p+1|0;x=h;y=c[h>>2]|0}c[q>>2]=y;c[s>>2]=r;c[t>>2]=0;a[u]=1;a[v]=69;h=p+4|0;z=c[m>>2]|0;m=c[n>>2]|0;n=c[e>>2]|0;A=y;a:while(1){if((n|0)!=0){y=c[n+12>>2]|0;if((y|0)==(c[n+16>>2]|0)){B=mc[c[(c[n>>2]|0)+36>>2]&63](n)|0}else{B=c[y>>2]|0}if((B|0)==-1){c[e>>2]=0;C=1;D=0}else{C=0;D=n}}else{C=1;D=0}y=c[f>>2]|0;do{if((y|0)!=0){E=c[y+12>>2]|0;if((E|0)==(c[y+16>>2]|0)){F=mc[c[(c[y>>2]|0)+36>>2]&63](y)|0}else{F=c[E>>2]|0}if(!((F|0)==-1)){if(C){G=y;break}else{H=y;I=A;break a}}else{c[f>>2]=0;J=17;break}}else{J=17}}while(0);if((J|0)==17){J=0;if(C){H=0;I=A;break}else{G=0}}y=a[p]|0;E=(y&1)==0;if(E){K=(y&255)>>>1}else{K=c[h>>2]|0}if(((c[q>>2]|0)-A|0)==(K|0)){if(E){L=(y&255)>>>1;M=(y&255)>>>1}else{y=c[h>>2]|0;L=y;M=y}Be(p,L<<1,0);if((a[p]&1)==0){N=10}else{N=(c[p>>2]&-2)+ -1|0}Be(p,N,0);if((a[p]&1)==0){O=w}else{O=c[x>>2]|0}c[q>>2]=O+M;P=O}else{P=A}y=D+12|0;E=c[y>>2]|0;Q=D+16|0;if((E|0)==(c[Q>>2]|0)){R=mc[c[(c[D>>2]|0)+36>>2]&63](D)|0}else{R=c[E>>2]|0}if((gh(R,u,v,P,q,z,m,o,r,s,t,l)|0)!=0){H=G;I=P;break}E=c[y>>2]|0;if((E|0)==(c[Q>>2]|0)){mc[c[(c[D>>2]|0)+40>>2]&63](D)|0;n=D;A=P;continue}else{c[y>>2]=E+4;n=D;A=P;continue}}P=a[o]|0;if((P&1)==0){S=(P&255)>>>1}else{S=c[o+4>>2]|0}if(((S|0)!=0?(a[u]|0)!=0:0)?(u=c[s>>2]|0,(u-r|0)<160):0){S=c[t>>2]|0;c[s>>2]=u+4;c[u>>2]=S}g[k>>2]=+Ql(I,c[q>>2]|0,j);nj(o,r,c[s>>2]|0,j);if((D|0)!=0){s=c[D+12>>2]|0;if((s|0)==(c[D+16>>2]|0)){T=mc[c[(c[D>>2]|0)+36>>2]&63](D)|0}else{T=c[s>>2]|0}if((T|0)==-1){c[e>>2]=0;U=0;V=1}else{U=D;V=0}}else{U=0;V=1}do{if((H|0)!=0){D=c[H+12>>2]|0;if((D|0)==(c[H+16>>2]|0)){W=mc[c[(c[H>>2]|0)+36>>2]&63](H)|0}else{W=c[D>>2]|0}if((W|0)==-1){c[f>>2]=0;J=57;break}if(V){c[b>>2]=U;ze(p);ze(o);i=d;return}}else{J=57}}while(0);if((J|0)==57?!V:0){c[b>>2]=U;ze(p);ze(o);i=d;return}c[j>>2]=c[j>>2]|2;c[b>>2]=U;ze(p);ze(o);i=d;return}function Xg(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];Yg(a,0,k,j,f,g,h);i=b;return}function Yg(b,d,e,f,g,j,k){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0;d=i;i=i+352|0;l=d+176|0;m=d+332|0;n=d+328|0;o=d+316|0;p=d+304|0;q=d+168|0;r=d+8|0;s=d+4|0;t=d;u=d+337|0;v=d+336|0;fh(o,g,l,m,n);c[p+0>>2]=0;c[p+4>>2]=0;c[p+8>>2]=0;Be(p,10,0);if((a[p]&1)==0){g=p+1|0;w=g;x=p+8|0;y=g}else{g=p+8|0;w=p+1|0;x=g;y=c[g>>2]|0}c[q>>2]=y;c[s>>2]=r;c[t>>2]=0;a[u]=1;a[v]=69;g=p+4|0;z=c[m>>2]|0;m=c[n>>2]|0;n=c[e>>2]|0;A=y;a:while(1){if((n|0)!=0){y=c[n+12>>2]|0;if((y|0)==(c[n+16>>2]|0)){B=mc[c[(c[n>>2]|0)+36>>2]&63](n)|0}else{B=c[y>>2]|0}if((B|0)==-1){c[e>>2]=0;C=1;D=0}else{C=0;D=n}}else{C=1;D=0}y=c[f>>2]|0;do{if((y|0)!=0){E=c[y+12>>2]|0;if((E|0)==(c[y+16>>2]|0)){F=mc[c[(c[y>>2]|0)+36>>2]&63](y)|0}else{F=c[E>>2]|0}if(!((F|0)==-1)){if(C){G=y;break}else{H=y;I=A;break a}}else{c[f>>2]=0;J=17;break}}else{J=17}}while(0);if((J|0)==17){J=0;if(C){H=0;I=A;break}else{G=0}}y=a[p]|0;E=(y&1)==0;if(E){K=(y&255)>>>1}else{K=c[g>>2]|0}if(((c[q>>2]|0)-A|0)==(K|0)){if(E){L=(y&255)>>>1;M=(y&255)>>>1}else{y=c[g>>2]|0;L=y;M=y}Be(p,L<<1,0);if((a[p]&1)==0){N=10}else{N=(c[p>>2]&-2)+ -1|0}Be(p,N,0);if((a[p]&1)==0){O=w}else{O=c[x>>2]|0}c[q>>2]=O+M;P=O}else{P=A}y=D+12|0;E=c[y>>2]|0;Q=D+16|0;if((E|0)==(c[Q>>2]|0)){R=mc[c[(c[D>>2]|0)+36>>2]&63](D)|0}else{R=c[E>>2]|0}if((gh(R,u,v,P,q,z,m,o,r,s,t,l)|0)!=0){H=G;I=P;break}E=c[y>>2]|0;if((E|0)==(c[Q>>2]|0)){mc[c[(c[D>>2]|0)+40>>2]&63](D)|0;n=D;A=P;continue}else{c[y>>2]=E+4;n=D;A=P;continue}}P=a[o]|0;if((P&1)==0){S=(P&255)>>>1}else{S=c[o+4>>2]|0}if(((S|0)!=0?(a[u]|0)!=0:0)?(u=c[s>>2]|0,(u-r|0)<160):0){S=c[t>>2]|0;c[s>>2]=u+4;c[u>>2]=S}h[k>>3]=+Pl(I,c[q>>2]|0,j);nj(o,r,c[s>>2]|0,j);if((D|0)!=0){s=c[D+12>>2]|0;if((s|0)==(c[D+16>>2]|0)){T=mc[c[(c[D>>2]|0)+36>>2]&63](D)|0}else{T=c[s>>2]|0}if((T|0)==-1){c[e>>2]=0;U=0;V=1}else{U=D;V=0}}else{U=0;V=1}do{if((H|0)!=0){D=c[H+12>>2]|0;if((D|0)==(c[H+16>>2]|0)){W=mc[c[(c[H>>2]|0)+36>>2]&63](H)|0}else{W=c[D>>2]|0}if((W|0)==-1){c[f>>2]=0;J=57;break}if(V){c[b>>2]=U;ze(p);ze(o);i=d;return}}else{J=57}}while(0);if((J|0)==57?!V:0){c[b>>2]=U;ze(p);ze(o);i=d;return}c[j>>2]=c[j>>2]|2;c[b>>2]=U;ze(p);ze(o);i=d;return}function Zg(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0;b=i;i=i+16|0;j=b+12|0;k=b;l=b+4|0;m=b+8|0;c[l>>2]=c[d>>2];c[m>>2]=c[e>>2];c[k+0>>2]=c[l+0>>2];c[j+0>>2]=c[m+0>>2];_g(a,0,k,j,f,g,h);i=b;return}function _g(b,d,e,f,g,j,k){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0;d=i;i=i+352|0;l=d+176|0;m=d+332|0;n=d+328|0;o=d+316|0;p=d+304|0;q=d+168|0;r=d+8|0;s=d+4|0;t=d;u=d+337|0;v=d+336|0;fh(o,g,l,m,n);c[p+0>>2]=0;c[p+4>>2]=0;c[p+8>>2]=0;Be(p,10,0);if((a[p]&1)==0){g=p+1|0;w=g;x=p+8|0;y=g}else{g=p+8|0;w=p+1|0;x=g;y=c[g>>2]|0}c[q>>2]=y;c[s>>2]=r;c[t>>2]=0;a[u]=1;a[v]=69;g=p+4|0;z=c[m>>2]|0;m=c[n>>2]|0;n=c[e>>2]|0;A=y;a:while(1){if((n|0)!=0){y=c[n+12>>2]|0;if((y|0)==(c[n+16>>2]|0)){B=mc[c[(c[n>>2]|0)+36>>2]&63](n)|0}else{B=c[y>>2]|0}if((B|0)==-1){c[e>>2]=0;C=1;D=0}else{C=0;D=n}}else{C=1;D=0}y=c[f>>2]|0;do{if((y|0)!=0){E=c[y+12>>2]|0;if((E|0)==(c[y+16>>2]|0)){F=mc[c[(c[y>>2]|0)+36>>2]&63](y)|0}else{F=c[E>>2]|0}if(!((F|0)==-1)){if(C){G=y;break}else{H=y;I=A;break a}}else{c[f>>2]=0;J=17;break}}else{J=17}}while(0);if((J|0)==17){J=0;if(C){H=0;I=A;break}else{G=0}}y=a[p]|0;E=(y&1)==0;if(E){K=(y&255)>>>1}else{K=c[g>>2]|0}if(((c[q>>2]|0)-A|0)==(K|0)){if(E){L=(y&255)>>>1;M=(y&255)>>>1}else{y=c[g>>2]|0;L=y;M=y}Be(p,L<<1,0);if((a[p]&1)==0){N=10}else{N=(c[p>>2]&-2)+ -1|0}Be(p,N,0);if((a[p]&1)==0){O=w}else{O=c[x>>2]|0}c[q>>2]=O+M;P=O}else{P=A}y=D+12|0;E=c[y>>2]|0;Q=D+16|0;if((E|0)==(c[Q>>2]|0)){R=mc[c[(c[D>>2]|0)+36>>2]&63](D)|0}else{R=c[E>>2]|0}if((gh(R,u,v,P,q,z,m,o,r,s,t,l)|0)!=0){H=G;I=P;break}E=c[y>>2]|0;if((E|0)==(c[Q>>2]|0)){mc[c[(c[D>>2]|0)+40>>2]&63](D)|0;n=D;A=P;continue}else{c[y>>2]=E+4;n=D;A=P;continue}}P=a[o]|0;if((P&1)==0){S=(P&255)>>>1}else{S=c[o+4>>2]|0}if(((S|0)!=0?(a[u]|0)!=0:0)?(u=c[s>>2]|0,(u-r|0)<160):0){S=c[t>>2]|0;c[s>>2]=u+4;c[u>>2]=S}h[k>>3]=+Ol(I,c[q>>2]|0,j);nj(o,r,c[s>>2]|0,j);if((D|0)!=0){s=c[D+12>>2]|0;if((s|0)==(c[D+16>>2]|0)){T=mc[c[(c[D>>2]|0)+36>>2]&63](D)|0}else{T=c[s>>2]|0}if((T|0)==-1){c[e>>2]=0;U=0;V=1}else{U=D;V=0}}else{U=0;V=1}do{if((H|0)!=0){D=c[H+12>>2]|0;if((D|0)==(c[H+16>>2]|0)){W=mc[c[(c[H>>2]|0)+36>>2]&63](H)|0}else{W=c[D>>2]|0}if((W|0)==-1){c[f>>2]=0;J=57;break}if(V){c[b>>2]=U;ze(p);ze(o);i=d;return}}else{J=57}}while(0);if((J|0)==57?!V:0){c[b>>2]=U;ze(p);ze(o);i=d;return}c[j>>2]=c[j>>2]|2;c[b>>2]=U;ze(p);ze(o);i=d;return}function $g(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0,X=0,Y=0,Z=0,_=0,$=0;d=i;i=i+304|0;k=d;l=d+200|0;m=d+12|0;n=d+24|0;o=d+28|0;p=d+40|0;c[m+0>>2]=0;c[m+4>>2]=0;c[m+8>>2]=0;Te(n,g);g=c[n>>2]|0;if(!((c[2062]|0)==-1)){c[k>>2]=8248;c[k+4>>2]=112;c[k+8>>2]=0;ue(8248,k,113)}q=(c[8252>>2]|0)+ -1|0;r=c[g+8>>2]|0;if(!((c[g+12>>2]|0)-r>>2>>>0>q>>>0)){s=Wa(4)|0;qm(s);xb(s|0,16216,101)}g=c[r+(q<<2)>>2]|0;if((g|0)==0){s=Wa(4)|0;qm(s);xb(s|0,16216,101)}sc[c[(c[g>>2]|0)+48>>2]&7](g,6800,6826|0,l)|0;ae(c[n>>2]|0)|0;c[o+0>>2]=0;c[o+4>>2]=0;c[o+8>>2]=0;Be(o,10,0);if((a[o]&1)==0){n=o+1|0;t=n;u=o+8|0;v=n}else{n=o+8|0;t=o+1|0;u=n;v=c[n>>2]|0}n=o+4|0;g=l+96|0;s=l+100|0;q=p;r=l+104|0;w=l;x=m+4|0;y=c[e>>2]|0;z=p;p=0;A=v;B=v;a:while(1){if((y|0)!=0){v=c[y+12>>2]|0;if((v|0)==(c[y+16>>2]|0)){C=mc[c[(c[y>>2]|0)+36>>2]&63](y)|0}else{C=c[v>>2]|0}if((C|0)==-1){c[e>>2]=0;D=1;E=0}else{D=0;E=y}}else{D=1;E=0}v=c[f>>2]|0;do{if((v|0)!=0){F=c[v+12>>2]|0;if((F|0)==(c[v+16>>2]|0)){G=mc[c[(c[v>>2]|0)+36>>2]&63](v)|0}else{G=c[F>>2]|0}if(!((G|0)==-1)){if(D){break}else{H=B;break a}}else{c[f>>2]=0;I=22;break}}else{I=22}}while(0);if((I|0)==22?(I=0,D):0){H=B;break}v=a[o]|0;F=(v&1)==0;if(F){J=(v&255)>>>1}else{J=c[n>>2]|0}if((A-B|0)==(J|0)){if(F){K=(v&255)>>>1;L=(v&255)>>>1}else{v=c[n>>2]|0;K=v;L=v}Be(o,K<<1,0);if((a[o]&1)==0){M=10}else{M=(c[o>>2]&-2)+ -1|0}Be(o,M,0);if((a[o]&1)==0){N=t}else{N=c[u>>2]|0}O=N+L|0;P=N}else{O=A;P=B}v=c[E+12>>2]|0;if((v|0)==(c[E+16>>2]|0)){Q=mc[c[(c[E>>2]|0)+36>>2]&63](E)|0}else{Q=c[v>>2]|0}v=(O|0)==(P|0);do{if(v){F=(c[g>>2]|0)==(Q|0);if(!F?(c[s>>2]|0)!=(Q|0):0){I=43;break}a[O]=F?43:45;R=O+1|0;S=z;T=0}else{I=43}}while(0);do{if((I|0)==43){I=0;F=a[m]|0;if((F&1)==0){U=(F&255)>>>1}else{U=c[x>>2]|0}if((U|0)!=0&(Q|0)==0){if((z-q|0)>=160){R=O;S=z;T=p;break}c[z>>2]=p;R=O;S=z+4|0;T=0;break}else{V=l}while(1){F=V+4|0;if((c[V>>2]|0)==(Q|0)){W=V;break}if((F|0)==(r|0)){W=r;break}else{V=F}}F=W-w|0;X=F>>2;if((F|0)>92){H=P;break a}if((F|0)<88){a[O]=a[6800+X|0]|0;R=O+1|0;S=z;T=p+1|0;break}if(v){H=O;break a}if((O-P|0)>=3){H=P;break a}if((a[O+ -1|0]|0)!=48){H=P;break a}a[O]=a[6800+X|0]|0;R=O+1|0;S=z;T=0}}while(0);v=c[e>>2]|0;X=v+12|0;F=c[X>>2]|0;if((F|0)==(c[v+16>>2]|0)){mc[c[(c[v>>2]|0)+40>>2]&63](v)|0;y=v;z=S;p=T;A=R;B=P;continue}else{c[X>>2]=F+4;y=v;z=S;p=T;A=R;B=P;continue}}a[H+3|0]=0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}P=c[2036]|0;c[k>>2]=j;if((Eg(H,P,6840,k)|0)!=1){c[h>>2]=4}k=c[e>>2]|0;if((k|0)!=0){P=c[k+12>>2]|0;if((P|0)==(c[k+16>>2]|0)){Y=mc[c[(c[k>>2]|0)+36>>2]&63](k)|0}else{Y=c[P>>2]|0}if((Y|0)==-1){c[e>>2]=0;Z=0;_=1}else{Z=k;_=0}}else{Z=0;_=1}k=c[f>>2]|0;do{if((k|0)!=0){e=c[k+12>>2]|0;if((e|0)==(c[k+16>>2]|0)){$=mc[c[(c[k>>2]|0)+36>>2]&63](k)|0}else{$=c[e>>2]|0}if(($|0)==-1){c[f>>2]=0;I=78;break}if(_){c[b>>2]=Z;ze(o);ze(m);i=d;return}}else{I=78}}while(0);if((I|0)==78?!_:0){c[b>>2]=Z;ze(o);ze(m);i=d;return}c[h>>2]=c[h>>2]|2;c[b>>2]=Z;ze(o);ze(m);i=d;return}function ah(b,d,e,f,g,h,j,k,l,m){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;m=m|0;var n=0,o=0,p=0,q=0,r=0,s=0,t=0;n=i;o=c[f>>2]|0;p=(o|0)==(e|0);do{if(p){q=(c[m+96>>2]|0)==(b|0);if(!q?(c[m+100>>2]|0)!=(b|0):0){break}c[f>>2]=e+1;a[e]=q?43:45;c[g>>2]=0;r=0;i=n;return r|0}}while(0);q=a[j]|0;if((q&1)==0){s=(q&255)>>>1}else{s=c[j+4>>2]|0}if((s|0)!=0&(b|0)==(h|0)){h=c[l>>2]|0;if((h-k|0)>=160){r=0;i=n;return r|0}k=c[g>>2]|0;c[l>>2]=h+4;c[h>>2]=k;c[g>>2]=0;r=0;i=n;return r|0}k=m+104|0;h=m;while(1){l=h+4|0;if((c[h>>2]|0)==(b|0)){t=h;break}if((l|0)==(k|0)){t=k;break}else{h=l}}h=t-m|0;m=h>>2;if((h|0)>92){r=-1;i=n;return r|0}if((d|0)==10|(d|0)==8){if((m|0)>=(d|0)){r=-1;i=n;return r|0}}else if((d|0)==16?(h|0)>=88:0){if(p){r=-1;i=n;return r|0}if((o-e|0)>=3){r=-1;i=n;return r|0}if((a[o+ -1|0]|0)!=48){r=-1;i=n;return r|0}c[g>>2]=0;e=a[6800+m|0]|0;c[f>>2]=o+1;a[o]=e;r=0;i=n;return r|0}e=a[6800+m|0]|0;c[f>>2]=o+1;a[o]=e;c[g>>2]=(c[g>>2]|0)+1;r=0;i=n;return r|0}function bh(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0;g=i;i=i+16|0;h=g;j=g+12|0;Te(j,d);d=c[j>>2]|0;if(!((c[2064]|0)==-1)){c[h>>2]=8256;c[h+4>>2]=112;c[h+8>>2]=0;ue(8256,h,113)}k=(c[8260>>2]|0)+ -1|0;l=c[d+8>>2]|0;if(!((c[d+12>>2]|0)-l>>2>>>0>k>>>0)){m=Wa(4)|0;qm(m);xb(m|0,16216,101)}d=c[l+(k<<2)>>2]|0;if((d|0)==0){m=Wa(4)|0;qm(m);xb(m|0,16216,101)}sc[c[(c[d>>2]|0)+32>>2]&7](d,6800,6826|0,e)|0;e=c[j>>2]|0;if(!((c[2100]|0)==-1)){c[h>>2]=8400;c[h+4>>2]=112;c[h+8>>2]=0;ue(8400,h,113)}h=(c[8404>>2]|0)+ -1|0;d=c[e+8>>2]|0;if(!((c[e+12>>2]|0)-d>>2>>>0>h>>>0)){n=Wa(4)|0;qm(n);xb(n|0,16216,101)}e=c[d+(h<<2)>>2]|0;if((e|0)==0){n=Wa(4)|0;qm(n);xb(n|0,16216,101)}else{a[f]=mc[c[(c[e>>2]|0)+16>>2]&63](e)|0;kc[c[(c[e>>2]|0)+20>>2]&63](b,e);ae(c[j>>2]|0)|0;i=g;return}}function ch(b,d,e,f,g){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;var h=0,j=0,k=0,l=0,m=0,n=0,o=0;h=i;i=i+16|0;j=h;k=h+12|0;Te(k,d);d=c[k>>2]|0;if(!((c[2064]|0)==-1)){c[j>>2]=8256;c[j+4>>2]=112;c[j+8>>2]=0;ue(8256,j,113)}l=(c[8260>>2]|0)+ -1|0;m=c[d+8>>2]|0;if(!((c[d+12>>2]|0)-m>>2>>>0>l>>>0)){n=Wa(4)|0;qm(n);xb(n|0,16216,101)}d=c[m+(l<<2)>>2]|0;if((d|0)==0){n=Wa(4)|0;qm(n);xb(n|0,16216,101)}sc[c[(c[d>>2]|0)+32>>2]&7](d,6800,6832|0,e)|0;e=c[k>>2]|0;if(!((c[2100]|0)==-1)){c[j>>2]=8400;c[j+4>>2]=112;c[j+8>>2]=0;ue(8400,j,113)}j=(c[8404>>2]|0)+ -1|0;d=c[e+8>>2]|0;if(!((c[e+12>>2]|0)-d>>2>>>0>j>>>0)){o=Wa(4)|0;qm(o);xb(o|0,16216,101)}e=c[d+(j<<2)>>2]|0;if((e|0)==0){o=Wa(4)|0;qm(o);xb(o|0,16216,101)}else{a[f]=mc[c[(c[e>>2]|0)+12>>2]&63](e)|0;a[g]=mc[c[(c[e>>2]|0)+16>>2]&63](e)|0;kc[c[(c[e>>2]|0)+20>>2]&63](b,e);ae(c[k>>2]|0)|0;i=h;return}}function dh(b,d,e,f,g,h,j,k,l,m,n,o){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;m=m|0;n=n|0;o=o|0;var p=0,q=0,r=0,s=0,t=0,u=0;p=i;if(b<<24>>24==h<<24>>24){if((a[d]|0)==0){q=-1;i=p;return q|0}a[d]=0;h=c[g>>2]|0;c[g>>2]=h+1;a[h]=46;h=a[k]|0;if((h&1)==0){r=(h&255)>>>1}else{r=c[k+4>>2]|0}if((r|0)==0){q=0;i=p;return q|0}r=c[m>>2]|0;if((r-l|0)>=160){q=0;i=p;return q|0}h=c[n>>2]|0;c[m>>2]=r+4;c[r>>2]=h;q=0;i=p;return q|0}if(b<<24>>24==j<<24>>24){j=a[k]|0;if((j&1)==0){s=(j&255)>>>1}else{s=c[k+4>>2]|0}if((s|0)!=0){if((a[d]|0)==0){q=-1;i=p;return q|0}s=c[m>>2]|0;if((s-l|0)>=160){q=0;i=p;return q|0}j=c[n>>2]|0;c[m>>2]=s+4;c[s>>2]=j;c[n>>2]=0;q=0;i=p;return q|0}}j=o+32|0;s=o;while(1){h=s+1|0;if((a[s]|0)==b<<24>>24){t=s;break}if((h|0)==(j|0)){t=j;break}else{s=h}}s=t-o|0;if((s|0)>31){q=-1;i=p;return q|0}o=a[6800+s|0]|0;if((s|0)==23|(s|0)==22){a[e]=80;t=c[g>>2]|0;c[g>>2]=t+1;a[t]=o;q=0;i=p;return q|0}else if((s|0)==24|(s|0)==25){t=c[g>>2]|0;if((t|0)!=(f|0)?(a[t+ -1|0]&95|0)!=(a[e]&127|0):0){q=-1;i=p;return q|0}c[g>>2]=t+1;a[t]=o;q=0;i=p;return q|0}else{t=o&95;if((t|0)==(a[e]|0)?(a[e]=t|128,(a[d]|0)!=0):0){a[d]=0;d=a[k]|0;if((d&1)==0){u=(d&255)>>>1}else{u=c[k+4>>2]|0}if((u|0)!=0?(u=c[m>>2]|0,(u-l|0)<160):0){l=c[n>>2]|0;c[m>>2]=u+4;c[u>>2]=l}}l=c[g>>2]|0;c[g>>2]=l+1;a[l]=o;if((s|0)>21){q=0;i=p;return q|0}c[n>>2]=(c[n>>2]|0)+1;q=0;i=p;return q|0}return 0}function eh(a,b,d,e){a=a|0;b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0,k=0,l=0,m=0;f=i;i=i+16|0;g=f;h=f+12|0;Te(h,b);b=c[h>>2]|0;if(!((c[2062]|0)==-1)){c[g>>2]=8248;c[g+4>>2]=112;c[g+8>>2]=0;ue(8248,g,113)}j=(c[8252>>2]|0)+ -1|0;k=c[b+8>>2]|0;if(!((c[b+12>>2]|0)-k>>2>>>0>j>>>0)){l=Wa(4)|0;qm(l);xb(l|0,16216,101)}b=c[k+(j<<2)>>2]|0;if((b|0)==0){l=Wa(4)|0;qm(l);xb(l|0,16216,101)}sc[c[(c[b>>2]|0)+48>>2]&7](b,6800,6826|0,d)|0;d=c[h>>2]|0;if(!((c[2102]|0)==-1)){c[g>>2]=8408;c[g+4>>2]=112;c[g+8>>2]=0;ue(8408,g,113)}g=(c[8412>>2]|0)+ -1|0;b=c[d+8>>2]|0;if(!((c[d+12>>2]|0)-b>>2>>>0>g>>>0)){m=Wa(4)|0;qm(m);xb(m|0,16216,101)}d=c[b+(g<<2)>>2]|0;if((d|0)==0){m=Wa(4)|0;qm(m);xb(m|0,16216,101)}else{c[e>>2]=mc[c[(c[d>>2]|0)+16>>2]&63](d)|0;kc[c[(c[d>>2]|0)+20>>2]&63](a,d);ae(c[h>>2]|0)|0;i=f;return}}function fh(a,b,d,e,f){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0;g=i;i=i+16|0;h=g;j=g+12|0;Te(j,b);b=c[j>>2]|0;if(!((c[2062]|0)==-1)){c[h>>2]=8248;c[h+4>>2]=112;c[h+8>>2]=0;ue(8248,h,113)}k=(c[8252>>2]|0)+ -1|0;l=c[b+8>>2]|0;if(!((c[b+12>>2]|0)-l>>2>>>0>k>>>0)){m=Wa(4)|0;qm(m);xb(m|0,16216,101)}b=c[l+(k<<2)>>2]|0;if((b|0)==0){m=Wa(4)|0;qm(m);xb(m|0,16216,101)}sc[c[(c[b>>2]|0)+48>>2]&7](b,6800,6832|0,d)|0;d=c[j>>2]|0;if(!((c[2102]|0)==-1)){c[h>>2]=8408;c[h+4>>2]=112;c[h+8>>2]=0;ue(8408,h,113)}h=(c[8412>>2]|0)+ -1|0;b=c[d+8>>2]|0;if(!((c[d+12>>2]|0)-b>>2>>>0>h>>>0)){n=Wa(4)|0;qm(n);xb(n|0,16216,101)}d=c[b+(h<<2)>>2]|0;if((d|0)==0){n=Wa(4)|0;qm(n);xb(n|0,16216,101)}else{c[e>>2]=mc[c[(c[d>>2]|0)+12>>2]&63](d)|0;c[f>>2]=mc[c[(c[d>>2]|0)+16>>2]&63](d)|0;kc[c[(c[d>>2]|0)+20>>2]&63](a,d);ae(c[j>>2]|0)|0;i=g;return}}function gh(b,d,e,f,g,h,j,k,l,m,n,o){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;m=m|0;n=n|0;o=o|0;var p=0,q=0,r=0,s=0,t=0,u=0;p=i;if((b|0)==(h|0)){if((a[d]|0)==0){q=-1;i=p;return q|0}a[d]=0;h=c[g>>2]|0;c[g>>2]=h+1;a[h]=46;h=a[k]|0;if((h&1)==0){r=(h&255)>>>1}else{r=c[k+4>>2]|0}if((r|0)==0){q=0;i=p;return q|0}r=c[m>>2]|0;if((r-l|0)>=160){q=0;i=p;return q|0}h=c[n>>2]|0;c[m>>2]=r+4;c[r>>2]=h;q=0;i=p;return q|0}if((b|0)==(j|0)){j=a[k]|0;if((j&1)==0){s=(j&255)>>>1}else{s=c[k+4>>2]|0}if((s|0)!=0){if((a[d]|0)==0){q=-1;i=p;return q|0}s=c[m>>2]|0;if((s-l|0)>=160){q=0;i=p;return q|0}j=c[n>>2]|0;c[m>>2]=s+4;c[s>>2]=j;c[n>>2]=0;q=0;i=p;return q|0}}j=o+128|0;s=o;while(1){h=s+4|0;if((c[s>>2]|0)==(b|0)){t=s;break}if((h|0)==(j|0)){t=j;break}else{s=h}}s=t-o|0;o=s>>2;if((s|0)>124){q=-1;i=p;return q|0}t=a[6800+o|0]|0;if((o|0)==24|(o|0)==25){j=c[g>>2]|0;if((j|0)!=(f|0)?(a[j+ -1|0]&95|0)!=(a[e]&127|0):0){q=-1;i=p;return q|0}c[g>>2]=j+1;a[j]=t;q=0;i=p;return q|0}else if(!((o|0)==23|(o|0)==22)){o=t&95;if((o|0)==(a[e]|0)?(a[e]=o|128,(a[d]|0)!=0):0){a[d]=0;d=a[k]|0;if((d&1)==0){u=(d&255)>>>1}else{u=c[k+4>>2]|0}if((u|0)!=0?(u=c[m>>2]|0,(u-l|0)<160):0){l=c[n>>2]|0;c[m>>2]=u+4;c[u>>2]=l}}}else{a[e]=80}e=c[g>>2]|0;c[g>>2]=e+1;a[e]=t;if((s|0)>84){q=0;i=p;return q|0}c[n>>2]=(c[n>>2]|0)+1;q=0;i=p;return q|0}function hh(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function ih(a){a=a|0;return}function jh(b,d,e,f,g,h){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0;j=i;i=i+32|0;k=j;l=j+28|0;m=j+12|0;n=j+16|0;if((c[f+4>>2]&1|0)==0){o=c[(c[d>>2]|0)+24>>2]|0;c[l>>2]=c[e>>2];p=h&1;c[k+0>>2]=c[l+0>>2];uc[o&15](b,d,k,f,g,p);i=j;return}Te(m,f);f=c[m>>2]|0;if(!((c[2100]|0)==-1)){c[k>>2]=8400;c[k+4>>2]=112;c[k+8>>2]=0;ue(8400,k,113)}k=(c[8404>>2]|0)+ -1|0;p=c[f+8>>2]|0;if(!((c[f+12>>2]|0)-p>>2>>>0>k>>>0)){q=Wa(4)|0;qm(q);xb(q|0,16216,101)}f=c[p+(k<<2)>>2]|0;if((f|0)==0){q=Wa(4)|0;qm(q);xb(q|0,16216,101)}ae(c[m>>2]|0)|0;m=c[f>>2]|0;if(h){kc[c[m+24>>2]&63](n,f)}else{kc[c[m+28>>2]&63](n,f)}f=a[n]|0;if((f&1)==0){m=n+1|0;r=m;s=m;t=n+8|0}else{m=n+8|0;r=c[m>>2]|0;s=n+1|0;t=m}m=n+4|0;h=f;f=r;while(1){if((h&1)==0){u=s;v=(h&255)>>>1}else{u=c[t>>2]|0;v=c[m>>2]|0}if((f|0)==(u+v|0)){break}r=a[f]|0;q=c[e>>2]|0;do{if((q|0)!=0){k=q+24|0;p=c[k>>2]|0;if((p|0)!=(c[q+28>>2]|0)){c[k>>2]=p+1;a[p]=r;break}if((vc[c[(c[q>>2]|0)+52>>2]&15](q,r&255)|0)==-1){c[e>>2]=0}}}while(0);h=a[n]|0;f=f+1|0}c[b>>2]=c[e>>2];ze(n);i=j;return}function kh(b,d,e,f,g,h){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0;d=i;i=i+64|0;j=d;k=d+56|0;l=d+44|0;m=d+20|0;n=d+12|0;o=d+8|0;p=d+4|0;q=d+16|0;a[k+0|0]=a[7048|0]|0;a[k+1|0]=a[7049|0]|0;a[k+2|0]=a[7050|0]|0;a[k+3|0]=a[7051|0]|0;a[k+4|0]=a[7052|0]|0;a[k+5|0]=a[7053|0]|0;r=k+1|0;s=f+4|0;t=c[s>>2]|0;if((t&2048|0)==0){u=r}else{a[r]=43;u=k+2|0}if((t&512|0)==0){v=u}else{a[u]=35;v=u+1|0}a[v]=108;u=v+1|0;v=t&74;do{if((v|0)==8){if((t&16384|0)==0){a[u]=120;break}else{a[u]=88;break}}else if((v|0)==64){a[u]=111}else{a[u]=100}}while(0);if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}u=c[2036]|0;c[j>>2]=h;h=lh(l,12,u,k,j)|0;k=l+h|0;u=c[s>>2]&176;do{if((u|0)==16){s=a[l]|0;if(s<<24>>24==43|s<<24>>24==45){w=l+1|0;break}if((h|0)>1&s<<24>>24==48?(s=a[l+1|0]|0,s<<24>>24==88|s<<24>>24==120):0){w=l+2|0}else{x=20}}else if((u|0)==32){w=k}else{x=20}}while(0);if((x|0)==20){w=l}Te(p,f);mh(l,w,k,m,n,o,p);ae(c[p>>2]|0)|0;c[q>>2]=c[e>>2];e=c[n>>2]|0;n=c[o>>2]|0;c[j+0>>2]=c[q+0>>2];nh(b,j,m,e,n,f,g);i=d;return}function lh(a,b,d,e,f){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0;g=i;i=i+16|0;h=g;c[h>>2]=f;f=jb(d|0)|0;d=Lb(a|0,b|0,e|0,h|0)|0;if((f|0)==0){i=g;return d|0}jb(f|0)|0;i=g;return d|0}function mh(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0;k=i;i=i+32|0;l=k;m=k+12|0;n=c[j>>2]|0;if(!((c[2064]|0)==-1)){c[l>>2]=8256;c[l+4>>2]=112;c[l+8>>2]=0;ue(8256,l,113)}o=(c[8260>>2]|0)+ -1|0;p=c[n+8>>2]|0;if(!((c[n+12>>2]|0)-p>>2>>>0>o>>>0)){q=Wa(4)|0;qm(q);xb(q|0,16216,101)}n=c[p+(o<<2)>>2]|0;if((n|0)==0){q=Wa(4)|0;qm(q);xb(q|0,16216,101)}q=c[j>>2]|0;if(!((c[2100]|0)==-1)){c[l>>2]=8400;c[l+4>>2]=112;c[l+8>>2]=0;ue(8400,l,113)}l=(c[8404>>2]|0)+ -1|0;j=c[q+8>>2]|0;if(!((c[q+12>>2]|0)-j>>2>>>0>l>>>0)){r=Wa(4)|0;qm(r);xb(r|0,16216,101)}q=c[j+(l<<2)>>2]|0;if((q|0)==0){r=Wa(4)|0;qm(r);xb(r|0,16216,101)}kc[c[(c[q>>2]|0)+20>>2]&63](m,q);r=a[m]|0;if((r&1)==0){s=(r&255)>>>1}else{s=c[m+4>>2]|0}if((s|0)!=0){c[h>>2]=f;s=a[b]|0;if(s<<24>>24==43|s<<24>>24==45){r=vc[c[(c[n>>2]|0)+28>>2]&15](n,s)|0;s=c[h>>2]|0;c[h>>2]=s+1;a[s]=r;t=b+1|0}else{t=b}if(((e-t|0)>1?(a[t]|0)==48:0)?(r=t+1|0,s=a[r]|0,s<<24>>24==88|s<<24>>24==120):0){s=vc[c[(c[n>>2]|0)+28>>2]&15](n,48)|0;l=c[h>>2]|0;c[h>>2]=l+1;a[l]=s;s=vc[c[(c[n>>2]|0)+28>>2]&15](n,a[r]|0)|0;r=c[h>>2]|0;c[h>>2]=r+1;a[r]=s;u=t+2|0}else{u=t}if((u|0)!=(e|0)?(t=e+ -1|0,t>>>0>u>>>0):0){s=u;r=t;do{t=a[s]|0;a[s]=a[r]|0;a[r]=t;s=s+1|0;r=r+ -1|0}while(s>>>0<r>>>0)}r=mc[c[(c[q>>2]|0)+16>>2]&63](q)|0;if(u>>>0<e>>>0){q=m+1|0;s=m+4|0;t=m+8|0;l=0;j=0;o=u;while(1){p=(a[m]&1)==0;if((a[(p?q:c[t>>2]|0)+j|0]|0)!=0?(l|0)==(a[(p?q:c[t>>2]|0)+j|0]|0):0){p=c[h>>2]|0;c[h>>2]=p+1;a[p]=r;p=a[m]|0;if((p&1)==0){v=(p&255)>>>1}else{v=c[s>>2]|0}w=0;x=(j>>>0<(v+ -1|0)>>>0)+j|0}else{w=l;x=j}p=vc[c[(c[n>>2]|0)+28>>2]&15](n,a[o]|0)|0;y=c[h>>2]|0;c[h>>2]=y+1;a[y]=p;o=o+1|0;if(!(o>>>0<e>>>0)){break}else{l=w+1|0;j=x}}}x=f+(u-b)|0;u=c[h>>2]|0;if((x|0)!=(u|0)?(j=u+ -1|0,j>>>0>x>>>0):0){u=x;x=j;do{j=a[u]|0;a[u]=a[x]|0;a[x]=j;u=u+1|0;x=x+ -1|0}while(u>>>0<x>>>0)}}else{sc[c[(c[n>>2]|0)+32>>2]&7](n,b,e,f)|0;c[h>>2]=f+(e-b)}if((d|0)==(e|0)){z=c[h>>2]|0;c[g>>2]=z;ze(m);i=k;return}else{z=f+(d-b)|0;c[g>>2]=z;ze(m);i=k;return}}function nh(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0;k=i;i=i+16|0;l=k;m=c[d>>2]|0;if((m|0)==0){c[b>>2]=0;i=k;return}n=g;g=e;o=n-g|0;p=h+12|0;h=c[p>>2]|0;q=(h|0)>(o|0)?h-o|0:0;o=f;h=o-g|0;if((h|0)>0?(gc[c[(c[m>>2]|0)+48>>2]&31](m,e,h)|0)!=(h|0):0){c[d>>2]=0;c[b>>2]=0;i=k;return}do{if((q|0)>0){ye(l,q,j);if((a[l]&1)==0){r=l+1|0}else{r=c[l+8>>2]|0}if((gc[c[(c[m>>2]|0)+48>>2]&31](m,r,q)|0)==(q|0)){ze(l);break}c[d>>2]=0;c[b>>2]=0;ze(l);i=k;return}}while(0);l=n-o|0;if((l|0)>0?(gc[c[(c[m>>2]|0)+48>>2]&31](m,f,l)|0)!=(l|0):0){c[d>>2]=0;c[b>>2]=0;i=k;return}c[p>>2]=0;c[b>>2]=m;i=k;return}function oh(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0;d=i;i=i+96|0;k=d+8|0;l=d;m=d+74|0;n=d+32|0;o=d+20|0;p=d+24|0;q=d+16|0;r=d+28|0;s=l;c[s>>2]=37;c[s+4>>2]=0;s=l+1|0;t=f+4|0;u=c[t>>2]|0;if((u&2048|0)==0){v=s}else{a[s]=43;v=l+2|0}if((u&512|0)==0){w=v}else{a[v]=35;w=v+1|0}v=w+2|0;a[w]=108;a[w+1|0]=108;w=u&74;do{if((w|0)==64){a[v]=111}else if((w|0)==8){if((u&16384|0)==0){a[v]=120;break}else{a[v]=88;break}}else{a[v]=100}}while(0);if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}v=c[2036]|0;u=k;c[u>>2]=h;c[u+4>>2]=j;j=lh(m,22,v,l,k)|0;l=m+j|0;v=c[t>>2]&176;do{if((v|0)==16){t=a[m]|0;if(t<<24>>24==43|t<<24>>24==45){x=m+1|0;break}if((j|0)>1&t<<24>>24==48?(t=a[m+1|0]|0,t<<24>>24==88|t<<24>>24==120):0){x=m+2|0}else{y=20}}else if((v|0)==32){x=l}else{y=20}}while(0);if((y|0)==20){x=m}Te(q,f);mh(m,x,l,n,o,p,q);ae(c[q>>2]|0)|0;c[r>>2]=c[e>>2];e=c[o>>2]|0;o=c[p>>2]|0;c[k+0>>2]=c[r+0>>2];nh(b,k,n,e,o,f,g);i=d;return}function ph(b,d,e,f,g,h){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0;d=i;i=i+64|0;j=d;k=d+56|0;l=d+44|0;m=d+20|0;n=d+12|0;o=d+8|0;p=d+4|0;q=d+16|0;a[k+0|0]=a[7048|0]|0;a[k+1|0]=a[7049|0]|0;a[k+2|0]=a[7050|0]|0;a[k+3|0]=a[7051|0]|0;a[k+4|0]=a[7052|0]|0;a[k+5|0]=a[7053|0]|0;r=k+1|0;s=f+4|0;t=c[s>>2]|0;if((t&2048|0)==0){u=r}else{a[r]=43;u=k+2|0}if((t&512|0)==0){v=u}else{a[u]=35;v=u+1|0}a[v]=108;u=v+1|0;v=t&74;do{if((v|0)==8){if((t&16384|0)==0){a[u]=120;break}else{a[u]=88;break}}else if((v|0)==64){a[u]=111}else{a[u]=117}}while(0);if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}u=c[2036]|0;c[j>>2]=h;h=lh(l,12,u,k,j)|0;k=l+h|0;u=c[s>>2]&176;do{if((u|0)==32){w=k}else if((u|0)==16){s=a[l]|0;if(s<<24>>24==43|s<<24>>24==45){w=l+1|0;break}if((h|0)>1&s<<24>>24==48?(s=a[l+1|0]|0,s<<24>>24==88|s<<24>>24==120):0){w=l+2|0}else{x=20}}else{x=20}}while(0);if((x|0)==20){w=l}Te(p,f);mh(l,w,k,m,n,o,p);ae(c[p>>2]|0)|0;c[q>>2]=c[e>>2];e=c[n>>2]|0;n=c[o>>2]|0;c[j+0>>2]=c[q+0>>2];nh(b,j,m,e,n,f,g);i=d;return}function qh(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0;d=i;i=i+112|0;k=d+8|0;l=d;m=d+75|0;n=d+32|0;o=d+20|0;p=d+24|0;q=d+16|0;r=d+28|0;s=l;c[s>>2]=37;c[s+4>>2]=0;s=l+1|0;t=f+4|0;u=c[t>>2]|0;if((u&2048|0)==0){v=s}else{a[s]=43;v=l+2|0}if((u&512|0)==0){w=v}else{a[v]=35;w=v+1|0}v=w+2|0;a[w]=108;a[w+1|0]=108;w=u&74;do{if((w|0)==8){if((u&16384|0)==0){a[v]=120;break}else{a[v]=88;break}}else if((w|0)==64){a[v]=111}else{a[v]=117}}while(0);if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}v=c[2036]|0;w=k;c[w>>2]=h;c[w+4>>2]=j;j=lh(m,23,v,l,k)|0;l=m+j|0;v=c[t>>2]&176;do{if((v|0)==32){x=l}else if((v|0)==16){t=a[m]|0;if(t<<24>>24==43|t<<24>>24==45){x=m+1|0;break}if((j|0)>1&t<<24>>24==48?(t=a[m+1|0]|0,t<<24>>24==88|t<<24>>24==120):0){x=m+2|0}else{y=20}}else{y=20}}while(0);if((y|0)==20){x=m}Te(q,f);mh(m,x,l,n,o,p,q);ae(c[q>>2]|0)|0;c[r>>2]=c[e>>2];e=c[o>>2]|0;o=c[p>>2]|0;c[k+0>>2]=c[r+0>>2];nh(b,k,n,e,o,f,g);i=d;return}function rh(b,d,e,f,g,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;j=+j;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0;d=i;i=i+144|0;l=d+8|0;m=d;n=d+44|0;o=d+36|0;p=d+74|0;q=d+40|0;r=d+24|0;s=d+20|0;t=d+28|0;u=d+32|0;v=m;c[v>>2]=37;c[v+4>>2]=0;v=m+1|0;w=f+4|0;x=c[w>>2]|0;if((x&2048|0)==0){y=v}else{a[v]=43;y=m+2|0}if((x&1024|0)==0){z=y}else{a[y]=35;z=y+1|0}y=x&260;v=x>>>14;do{if((y|0)==260){if((v&1|0)==0){a[z]=97;A=0;break}else{a[z]=65;A=0;break}}else{a[z]=46;x=z+2|0;a[z+1|0]=42;if((y|0)==256){if((v&1|0)==0){a[x]=101;A=1;break}else{a[x]=69;A=1;break}}else if((y|0)==4){if((v&1|0)==0){a[x]=102;A=1;break}else{a[x]=70;A=1;break}}else{if((v&1|0)==0){a[x]=103;A=1;break}else{a[x]=71;A=1;break}}}}while(0);c[o>>2]=n;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}v=c[2036]|0;if(A){c[l>>2]=c[f+8>>2];y=l+4|0;h[k>>3]=j;c[y>>2]=c[k>>2];c[y+4>>2]=c[k+4>>2];B=lh(n,30,v,m,l)|0}else{h[k>>3]=j;c[l>>2]=c[k>>2];c[l+4>>2]=c[k+4>>2];B=lh(n,30,v,m,l)|0}if((B|0)>29){v=(a[8152]|0)==0;if(A){if(v?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}A=c[2036]|0;c[l>>2]=c[f+8>>2];y=l+4|0;h[k>>3]=j;c[y>>2]=c[k>>2];c[y+4>>2]=c[k+4>>2];C=sh(o,A,m,l)|0}else{if(v?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}v=c[2036]|0;c[l>>2]=c[f+8>>2];A=l+4|0;h[k>>3]=j;c[A>>2]=c[k>>2];c[A+4>>2]=c[k+4>>2];C=sh(o,v,m,l)|0}m=c[o>>2]|0;if((m|0)==0){_m()}else{D=m;E=m;F=C}}else{D=c[o>>2]|0;E=0;F=B}B=D+F|0;o=c[w>>2]&176;do{if((o|0)==16){w=a[D]|0;if(w<<24>>24==43|w<<24>>24==45){G=D+1|0;break}if((F|0)>1&w<<24>>24==48?(w=a[D+1|0]|0,w<<24>>24==88|w<<24>>24==120):0){G=D+2|0}else{H=44}}else if((o|0)==32){G=B}else{H=44}}while(0);if((H|0)==44){G=D}if((D|0)!=(n|0)){H=Lm(F<<1)|0;if((H|0)==0){_m()}else{I=D;J=H;K=H}}else{I=n;J=0;K=p}Te(s,f);th(I,G,B,K,q,r,s);ae(c[s>>2]|0)|0;c[u>>2]=c[e>>2];s=c[q>>2]|0;q=c[r>>2]|0;c[l+0>>2]=c[u+0>>2];nh(t,l,K,s,q,f,g);g=c[t>>2]|0;c[e>>2]=g;c[b>>2]=g;if((J|0)!=0){Mm(J)}if((E|0)==0){i=d;return}Mm(E);i=d;return}function sh(a,b,d,e){a=a|0;b=b|0;d=d|0;e=e|0;var f=0,g=0;f=i;i=i+16|0;g=f;c[g>>2]=e;e=jb(b|0)|0;b=Gb(a|0,d|0,g|0)|0;if((e|0)==0){i=f;return b|0}jb(e|0)|0;i=f;return b|0}function th(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0;k=i;i=i+32|0;l=k;m=k+12|0;n=c[j>>2]|0;if(!((c[2064]|0)==-1)){c[l>>2]=8256;c[l+4>>2]=112;c[l+8>>2]=0;ue(8256,l,113)}o=(c[8260>>2]|0)+ -1|0;p=c[n+8>>2]|0;if(!((c[n+12>>2]|0)-p>>2>>>0>o>>>0)){q=Wa(4)|0;qm(q);xb(q|0,16216,101)}n=c[p+(o<<2)>>2]|0;if((n|0)==0){q=Wa(4)|0;qm(q);xb(q|0,16216,101)}q=c[j>>2]|0;if(!((c[2100]|0)==-1)){c[l>>2]=8400;c[l+4>>2]=112;c[l+8>>2]=0;ue(8400,l,113)}l=(c[8404>>2]|0)+ -1|0;j=c[q+8>>2]|0;if(!((c[q+12>>2]|0)-j>>2>>>0>l>>>0)){r=Wa(4)|0;qm(r);xb(r|0,16216,101)}q=c[j+(l<<2)>>2]|0;if((q|0)==0){r=Wa(4)|0;qm(r);xb(r|0,16216,101)}kc[c[(c[q>>2]|0)+20>>2]&63](m,q);c[h>>2]=f;r=a[b]|0;if(r<<24>>24==43|r<<24>>24==45){l=vc[c[(c[n>>2]|0)+28>>2]&15](n,r)|0;r=c[h>>2]|0;c[h>>2]=r+1;a[r]=l;s=b+1|0}else{s=b}l=e;a:do{if(((l-s|0)>1?(a[s]|0)==48:0)?(r=s+1|0,j=a[r]|0,j<<24>>24==88|j<<24>>24==120):0){j=vc[c[(c[n>>2]|0)+28>>2]&15](n,48)|0;o=c[h>>2]|0;c[h>>2]=o+1;a[o]=j;j=s+2|0;o=vc[c[(c[n>>2]|0)+28>>2]&15](n,a[r]|0)|0;r=c[h>>2]|0;c[h>>2]=r+1;a[r]=o;if(j>>>0<e>>>0){o=j;while(1){r=a[o]|0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}p=o+1|0;if((Ya(r<<24>>24|0,c[2036]|0)|0)==0){t=j;u=o;break a}if(p>>>0<e>>>0){o=p}else{t=j;u=p;break}}}else{t=j;u=j}}else{v=14}}while(0);b:do{if((v|0)==14){if(s>>>0<e>>>0){o=s;while(1){p=a[o]|0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}r=o+1|0;if((ub(p<<24>>24|0,c[2036]|0)|0)==0){t=s;u=o;break b}if(r>>>0<e>>>0){o=r}else{t=s;u=r;break}}}else{t=s;u=s}}}while(0);s=a[m]|0;if((s&1)==0){w=(s&255)>>>1}else{w=c[m+4>>2]|0}if((w|0)!=0){if((t|0)!=(u|0)?(w=u+ -1|0,w>>>0>t>>>0):0){s=t;v=w;do{w=a[s]|0;a[s]=a[v]|0;a[v]=w;s=s+1|0;v=v+ -1|0}while(s>>>0<v>>>0)}v=mc[c[(c[q>>2]|0)+16>>2]&63](q)|0;if(t>>>0<u>>>0){s=m+1|0;w=m+4|0;o=m+8|0;j=0;r=0;p=t;while(1){x=(a[m]&1)==0;if((a[(x?s:c[o>>2]|0)+r|0]|0)>0?(j|0)==(a[(x?s:c[o>>2]|0)+r|0]|0):0){x=c[h>>2]|0;c[h>>2]=x+1;a[x]=v;x=a[m]|0;if((x&1)==0){y=(x&255)>>>1}else{y=c[w>>2]|0}z=0;A=(r>>>0<(y+ -1|0)>>>0)+r|0}else{z=j;A=r}x=vc[c[(c[n>>2]|0)+28>>2]&15](n,a[p]|0)|0;B=c[h>>2]|0;c[h>>2]=B+1;a[B]=x;p=p+1|0;if(!(p>>>0<u>>>0)){break}else{j=z+1|0;r=A}}}A=f+(t-b)|0;r=c[h>>2]|0;if((A|0)!=(r|0)?(z=r+ -1|0,z>>>0>A>>>0):0){r=A;A=z;do{z=a[r]|0;a[r]=a[A]|0;a[A]=z;r=r+1|0;A=A+ -1|0}while(r>>>0<A>>>0)}}else{sc[c[(c[n>>2]|0)+32>>2]&7](n,t,u,c[h>>2]|0)|0;c[h>>2]=(c[h>>2]|0)+(u-t)}c:do{if(u>>>0<e>>>0){t=u;while(1){A=a[t]|0;if(A<<24>>24==46){break}r=vc[c[(c[n>>2]|0)+28>>2]&15](n,A)|0;A=c[h>>2]|0;c[h>>2]=A+1;a[A]=r;r=t+1|0;if(r>>>0<e>>>0){t=r}else{C=r;break c}}r=mc[c[(c[q>>2]|0)+12>>2]&63](q)|0;A=c[h>>2]|0;c[h>>2]=A+1;a[A]=r;C=t+1|0}else{C=u}}while(0);sc[c[(c[n>>2]|0)+32>>2]&7](n,C,e,c[h>>2]|0)|0;n=(c[h>>2]|0)+(l-C)|0;c[h>>2]=n;if((d|0)==(e|0)){D=n;c[g>>2]=D;ze(m);i=k;return}D=f+(d-b)|0;c[g>>2]=D;ze(m);i=k;return}function uh(b,d,e,f,g,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;j=+j;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0;d=i;i=i+144|0;l=d+8|0;m=d;n=d+44|0;o=d+36|0;p=d+74|0;q=d+40|0;r=d+24|0;s=d+20|0;t=d+28|0;u=d+32|0;v=m;c[v>>2]=37;c[v+4>>2]=0;v=m+1|0;w=f+4|0;x=c[w>>2]|0;if((x&2048|0)==0){y=v}else{a[v]=43;y=m+2|0}if((x&1024|0)==0){z=y}else{a[y]=35;z=y+1|0}y=x&260;v=x>>>14;do{if((y|0)==260){a[z]=76;x=z+1|0;if((v&1|0)==0){a[x]=97;A=0;break}else{a[x]=65;A=0;break}}else{a[z]=46;a[z+1|0]=42;a[z+2|0]=76;x=z+3|0;if((y|0)==4){if((v&1|0)==0){a[x]=102;A=1;break}else{a[x]=70;A=1;break}}else if((y|0)==256){if((v&1|0)==0){a[x]=101;A=1;break}else{a[x]=69;A=1;break}}else{if((v&1|0)==0){a[x]=103;A=1;break}else{a[x]=71;A=1;break}}}}while(0);c[o>>2]=n;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}v=c[2036]|0;if(A){c[l>>2]=c[f+8>>2];y=l+4|0;h[k>>3]=j;c[y>>2]=c[k>>2];c[y+4>>2]=c[k+4>>2];B=lh(n,30,v,m,l)|0}else{h[k>>3]=j;c[l>>2]=c[k>>2];c[l+4>>2]=c[k+4>>2];B=lh(n,30,v,m,l)|0}if((B|0)>29){v=(a[8152]|0)==0;if(A){if(v?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}A=c[2036]|0;c[l>>2]=c[f+8>>2];y=l+4|0;h[k>>3]=j;c[y>>2]=c[k>>2];c[y+4>>2]=c[k+4>>2];C=sh(o,A,m,l)|0}else{if(v?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}v=c[2036]|0;h[k>>3]=j;c[l>>2]=c[k>>2];c[l+4>>2]=c[k+4>>2];C=sh(o,v,m,l)|0}m=c[o>>2]|0;if((m|0)==0){_m()}else{D=m;E=m;F=C}}else{D=c[o>>2]|0;E=0;F=B}B=D+F|0;o=c[w>>2]&176;do{if((o|0)==16){w=a[D]|0;if(w<<24>>24==43|w<<24>>24==45){G=D+1|0;break}if((F|0)>1&w<<24>>24==48?(w=a[D+1|0]|0,w<<24>>24==88|w<<24>>24==120):0){G=D+2|0}else{H=44}}else if((o|0)==32){G=B}else{H=44}}while(0);if((H|0)==44){G=D}if((D|0)!=(n|0)){H=Lm(F<<1)|0;if((H|0)==0){_m()}else{I=D;J=H;K=H}}else{I=n;J=0;K=p}Te(s,f);th(I,G,B,K,q,r,s);ae(c[s>>2]|0)|0;c[u>>2]=c[e>>2];s=c[q>>2]|0;q=c[r>>2]|0;c[l+0>>2]=c[u+0>>2];nh(t,l,K,s,q,f,g);g=c[t>>2]|0;c[e>>2]=g;c[b>>2]=g;if((J|0)!=0){Mm(J)}if((E|0)==0){i=d;return}Mm(E);i=d;return}function vh(b,d,e,f,g,h){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0;d=i;i=i+96|0;j=d;k=d+80|0;l=d+60|0;m=d+20|0;n=d+12|0;o=d+16|0;a[k+0|0]=a[7056|0]|0;a[k+1|0]=a[7057|0]|0;a[k+2|0]=a[7058|0]|0;a[k+3|0]=a[7059|0]|0;a[k+4|0]=a[7060|0]|0;a[k+5|0]=a[7061|0]|0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}p=c[2036]|0;c[j>>2]=h;h=lh(l,20,p,k,j)|0;k=l+h|0;p=c[f+4>>2]&176;do{if((p|0)==16){q=a[l]|0;if(q<<24>>24==43|q<<24>>24==45){r=l+1|0;break}if((h|0)>1&q<<24>>24==48?(q=a[l+1|0]|0,q<<24>>24==88|q<<24>>24==120):0){r=l+2|0}else{s=10}}else if((p|0)==32){r=k}else{s=10}}while(0);if((s|0)==10){r=l}Te(n,f);s=c[n>>2]|0;if(!((c[2064]|0)==-1)){c[j>>2]=8256;c[j+4>>2]=112;c[j+8>>2]=0;ue(8256,j,113)}p=(c[8260>>2]|0)+ -1|0;q=c[s+8>>2]|0;if(!((c[s+12>>2]|0)-q>>2>>>0>p>>>0)){t=Wa(4)|0;qm(t);xb(t|0,16216,101)}s=c[q+(p<<2)>>2]|0;if((s|0)==0){t=Wa(4)|0;qm(t);xb(t|0,16216,101)}ae(c[n>>2]|0)|0;sc[c[(c[s>>2]|0)+32>>2]&7](s,l,k,m)|0;s=m+h|0;if((r|0)==(k|0)){u=s;v=c[e>>2]|0;c[o>>2]=v;c[j+0>>2]=c[o+0>>2];nh(b,j,m,u,s,f,g);i=d;return}u=m+(r-l)|0;v=c[e>>2]|0;c[o>>2]=v;c[j+0>>2]=c[o+0>>2];nh(b,j,m,u,s,f,g);i=d;return}function wh(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function xh(a){a=a|0;return}function yh(b,d,e,f,g,h){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0;j=i;i=i+32|0;k=j;l=j+28|0;m=j+12|0;n=j+16|0;if((c[f+4>>2]&1|0)==0){o=c[(c[d>>2]|0)+24>>2]|0;c[l>>2]=c[e>>2];p=h&1;c[k+0>>2]=c[l+0>>2];uc[o&15](b,d,k,f,g,p);i=j;return}Te(m,f);f=c[m>>2]|0;if(!((c[2102]|0)==-1)){c[k>>2]=8408;c[k+4>>2]=112;c[k+8>>2]=0;ue(8408,k,113)}k=(c[8412>>2]|0)+ -1|0;p=c[f+8>>2]|0;if(!((c[f+12>>2]|0)-p>>2>>>0>k>>>0)){q=Wa(4)|0;qm(q);xb(q|0,16216,101)}f=c[p+(k<<2)>>2]|0;if((f|0)==0){q=Wa(4)|0;qm(q);xb(q|0,16216,101)}ae(c[m>>2]|0)|0;m=c[f>>2]|0;if(h){kc[c[m+24>>2]&63](n,f)}else{kc[c[m+28>>2]&63](n,f)}f=a[n]|0;if((f&1)==0){m=n+4|0;r=m;s=n+8|0;t=m}else{m=n+8|0;r=c[m>>2]|0;s=m;t=n+4|0}m=f;f=r;while(1){if((m&1)==0){u=t;v=(m&255)>>>1}else{u=c[s>>2]|0;v=c[t>>2]|0}if((f|0)==(u+(v<<2)|0)){break}r=c[f>>2]|0;h=c[e>>2]|0;if((h|0)!=0){q=h+24|0;k=c[q>>2]|0;if((k|0)==(c[h+28>>2]|0)){w=vc[c[(c[h>>2]|0)+52>>2]&15](h,r)|0}else{c[q>>2]=k+4;c[k>>2]=r;w=r}if((w|0)==-1){c[e>>2]=0}}m=a[n]|0;f=f+4|0}c[b>>2]=c[e>>2];Ke(n);i=j;return}function zh(b,d,e,f,g,h){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0;d=i;i=i+128|0;j=d;k=d+116|0;l=d+104|0;m=d+8|0;n=d+92|0;o=d+96|0;p=d+4|0;q=d+100|0;a[k+0|0]=a[7048|0]|0;a[k+1|0]=a[7049|0]|0;a[k+2|0]=a[7050|0]|0;a[k+3|0]=a[7051|0]|0;a[k+4|0]=a[7052|0]|0;a[k+5|0]=a[7053|0]|0;r=k+1|0;s=f+4|0;t=c[s>>2]|0;if((t&2048|0)==0){u=r}else{a[r]=43;u=k+2|0}if((t&512|0)==0){v=u}else{a[u]=35;v=u+1|0}a[v]=108;u=v+1|0;v=t&74;do{if((v|0)==64){a[u]=111}else if((v|0)==8){if((t&16384|0)==0){a[u]=120;break}else{a[u]=88;break}}else{a[u]=100}}while(0);if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}u=c[2036]|0;c[j>>2]=h;h=lh(l,12,u,k,j)|0;k=l+h|0;u=c[s>>2]&176;do{if((u|0)==16){s=a[l]|0;if(s<<24>>24==43|s<<24>>24==45){w=l+1|0;break}if((h|0)>1&s<<24>>24==48?(s=a[l+1|0]|0,s<<24>>24==88|s<<24>>24==120):0){w=l+2|0}else{x=20}}else if((u|0)==32){w=k}else{x=20}}while(0);if((x|0)==20){w=l}Te(p,f);Ah(l,w,k,m,n,o,p);ae(c[p>>2]|0)|0;c[q>>2]=c[e>>2];e=c[n>>2]|0;n=c[o>>2]|0;c[j+0>>2]=c[q+0>>2];Bh(b,j,m,e,n,f,g);i=d;return}function Ah(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0;k=i;i=i+32|0;l=k;m=k+12|0;n=c[j>>2]|0;if(!((c[2062]|0)==-1)){c[l>>2]=8248;c[l+4>>2]=112;c[l+8>>2]=0;ue(8248,l,113)}o=(c[8252>>2]|0)+ -1|0;p=c[n+8>>2]|0;if(!((c[n+12>>2]|0)-p>>2>>>0>o>>>0)){q=Wa(4)|0;qm(q);xb(q|0,16216,101)}n=c[p+(o<<2)>>2]|0;if((n|0)==0){q=Wa(4)|0;qm(q);xb(q|0,16216,101)}q=c[j>>2]|0;if(!((c[2102]|0)==-1)){c[l>>2]=8408;c[l+4>>2]=112;c[l+8>>2]=0;ue(8408,l,113)}l=(c[8412>>2]|0)+ -1|0;j=c[q+8>>2]|0;if(!((c[q+12>>2]|0)-j>>2>>>0>l>>>0)){r=Wa(4)|0;qm(r);xb(r|0,16216,101)}q=c[j+(l<<2)>>2]|0;if((q|0)==0){r=Wa(4)|0;qm(r);xb(r|0,16216,101)}kc[c[(c[q>>2]|0)+20>>2]&63](m,q);r=a[m]|0;if((r&1)==0){s=(r&255)>>>1}else{s=c[m+4>>2]|0}if((s|0)!=0){c[h>>2]=f;s=a[b]|0;if(s<<24>>24==43|s<<24>>24==45){r=vc[c[(c[n>>2]|0)+44>>2]&15](n,s)|0;s=c[h>>2]|0;c[h>>2]=s+4;c[s>>2]=r;t=b+1|0}else{t=b}if(((e-t|0)>1?(a[t]|0)==48:0)?(r=t+1|0,s=a[r]|0,s<<24>>24==88|s<<24>>24==120):0){s=vc[c[(c[n>>2]|0)+44>>2]&15](n,48)|0;l=c[h>>2]|0;c[h>>2]=l+4;c[l>>2]=s;s=vc[c[(c[n>>2]|0)+44>>2]&15](n,a[r]|0)|0;r=c[h>>2]|0;c[h>>2]=r+4;c[r>>2]=s;u=t+2|0}else{u=t}if((u|0)!=(e|0)?(t=e+ -1|0,t>>>0>u>>>0):0){s=u;r=t;do{t=a[s]|0;a[s]=a[r]|0;a[r]=t;s=s+1|0;r=r+ -1|0}while(s>>>0<r>>>0)}r=mc[c[(c[q>>2]|0)+16>>2]&63](q)|0;if(u>>>0<e>>>0){q=m+1|0;s=m+4|0;t=m+8|0;l=0;j=0;o=u;while(1){p=(a[m]&1)==0;if((a[(p?q:c[t>>2]|0)+j|0]|0)!=0?(l|0)==(a[(p?q:c[t>>2]|0)+j|0]|0):0){p=c[h>>2]|0;c[h>>2]=p+4;c[p>>2]=r;p=a[m]|0;if((p&1)==0){v=(p&255)>>>1}else{v=c[s>>2]|0}w=0;x=(j>>>0<(v+ -1|0)>>>0)+j|0}else{w=l;x=j}p=vc[c[(c[n>>2]|0)+44>>2]&15](n,a[o]|0)|0;y=c[h>>2]|0;z=y+4|0;c[h>>2]=z;c[y>>2]=p;p=o+1|0;if(p>>>0<e>>>0){l=w+1|0;j=x;o=p}else{A=z;break}}}else{A=c[h>>2]|0}o=f+(u-b<<2)|0;if((o|0)!=(A|0)?(u=A+ -4|0,u>>>0>o>>>0):0){x=o;o=u;while(1){u=c[x>>2]|0;c[x>>2]=c[o>>2];c[o>>2]=u;u=x+4|0;j=o+ -4|0;if(u>>>0<j>>>0){x=u;o=j}else{B=A;break}}}else{B=A}}else{sc[c[(c[n>>2]|0)+48>>2]&7](n,b,e,f)|0;n=f+(e-b<<2)|0;c[h>>2]=n;B=n}if((d|0)==(e|0)){C=B;c[g>>2]=C;ze(m);i=k;return}C=f+(d-b<<2)|0;c[g>>2]=C;ze(m);i=k;return}function Bh(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0;k=i;i=i+16|0;l=k;m=c[d>>2]|0;if((m|0)==0){c[b>>2]=0;i=k;return}n=g;g=e;o=n-g>>2;p=h+12|0;h=c[p>>2]|0;q=(h|0)>(o|0)?h-o|0:0;o=f;h=o-g|0;g=h>>2;if((h|0)>0?(gc[c[(c[m>>2]|0)+48>>2]&31](m,e,g)|0)!=(g|0):0){c[d>>2]=0;c[b>>2]=0;i=k;return}do{if((q|0)>0){Je(l,q,j);if((a[l]&1)==0){r=l+4|0}else{r=c[l+8>>2]|0}if((gc[c[(c[m>>2]|0)+48>>2]&31](m,r,q)|0)==(q|0)){Ke(l);break}c[d>>2]=0;c[b>>2]=0;Ke(l);i=k;return}}while(0);l=n-o|0;o=l>>2;if((l|0)>0?(gc[c[(c[m>>2]|0)+48>>2]&31](m,f,o)|0)!=(o|0):0){c[d>>2]=0;c[b>>2]=0;i=k;return}c[p>>2]=0;c[b>>2]=m;i=k;return}function Ch(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0;d=i;i=i+224|0;k=d+8|0;l=d;m=d+196|0;n=d+16|0;o=d+180|0;p=d+184|0;q=d+188|0;r=d+192|0;s=l;c[s>>2]=37;c[s+4>>2]=0;s=l+1|0;t=f+4|0;u=c[t>>2]|0;if((u&2048|0)==0){v=s}else{a[s]=43;v=l+2|0}if((u&512|0)==0){w=v}else{a[v]=35;w=v+1|0}v=w+2|0;a[w]=108;a[w+1|0]=108;w=u&74;do{if((w|0)==64){a[v]=111}else if((w|0)==8){if((u&16384|0)==0){a[v]=120;break}else{a[v]=88;break}}else{a[v]=100}}while(0);if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}v=c[2036]|0;u=k;c[u>>2]=h;c[u+4>>2]=j;j=lh(m,22,v,l,k)|0;l=m+j|0;v=c[t>>2]&176;do{if((v|0)==32){x=l}else if((v|0)==16){t=a[m]|0;if(t<<24>>24==43|t<<24>>24==45){x=m+1|0;break}if((j|0)>1&t<<24>>24==48?(t=a[m+1|0]|0,t<<24>>24==88|t<<24>>24==120):0){x=m+2|0}else{y=20}}else{y=20}}while(0);if((y|0)==20){x=m}Te(q,f);Ah(m,x,l,n,o,p,q);ae(c[q>>2]|0)|0;c[r>>2]=c[e>>2];e=c[o>>2]|0;o=c[p>>2]|0;c[k+0>>2]=c[r+0>>2];Bh(b,k,n,e,o,f,g);i=d;return}function Dh(b,d,e,f,g,h){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0;d=i;i=i+128|0;j=d;k=d+116|0;l=d+104|0;m=d+8|0;n=d+92|0;o=d+96|0;p=d+4|0;q=d+100|0;a[k+0|0]=a[7048|0]|0;a[k+1|0]=a[7049|0]|0;a[k+2|0]=a[7050|0]|0;a[k+3|0]=a[7051|0]|0;a[k+4|0]=a[7052|0]|0;a[k+5|0]=a[7053|0]|0;r=k+1|0;s=f+4|0;t=c[s>>2]|0;if((t&2048|0)==0){u=r}else{a[r]=43;u=k+2|0}if((t&512|0)==0){v=u}else{a[u]=35;v=u+1|0}a[v]=108;u=v+1|0;v=t&74;do{if((v|0)==64){a[u]=111}else if((v|0)==8){if((t&16384|0)==0){a[u]=120;break}else{a[u]=88;break}}else{a[u]=117}}while(0);if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}u=c[2036]|0;c[j>>2]=h;h=lh(l,12,u,k,j)|0;k=l+h|0;u=c[s>>2]&176;do{if((u|0)==16){s=a[l]|0;if(s<<24>>24==43|s<<24>>24==45){w=l+1|0;break}if((h|0)>1&s<<24>>24==48?(s=a[l+1|0]|0,s<<24>>24==88|s<<24>>24==120):0){w=l+2|0}else{x=20}}else if((u|0)==32){w=k}else{x=20}}while(0);if((x|0)==20){w=l}Te(p,f);Ah(l,w,k,m,n,o,p);ae(c[p>>2]|0)|0;c[q>>2]=c[e>>2];e=c[n>>2]|0;n=c[o>>2]|0;c[j+0>>2]=c[q+0>>2];Bh(b,j,m,e,n,f,g);i=d;return}function Eh(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0;d=i;i=i+240|0;k=d+8|0;l=d;m=d+204|0;n=d+16|0;o=d+188|0;p=d+192|0;q=d+196|0;r=d+200|0;s=l;c[s>>2]=37;c[s+4>>2]=0;s=l+1|0;t=f+4|0;u=c[t>>2]|0;if((u&2048|0)==0){v=s}else{a[s]=43;v=l+2|0}if((u&512|0)==0){w=v}else{a[v]=35;w=v+1|0}v=w+2|0;a[w]=108;a[w+1|0]=108;w=u&74;do{if((w|0)==64){a[v]=111}else if((w|0)==8){if((u&16384|0)==0){a[v]=120;break}else{a[v]=88;break}}else{a[v]=117}}while(0);if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}v=c[2036]|0;u=k;c[u>>2]=h;c[u+4>>2]=j;j=lh(m,23,v,l,k)|0;l=m+j|0;v=c[t>>2]&176;do{if((v|0)==32){x=l}else if((v|0)==16){t=a[m]|0;if(t<<24>>24==43|t<<24>>24==45){x=m+1|0;break}if((j|0)>1&t<<24>>24==48?(t=a[m+1|0]|0,t<<24>>24==88|t<<24>>24==120):0){x=m+2|0}else{y=20}}else{y=20}}while(0);if((y|0)==20){x=m}Te(q,f);Ah(m,x,l,n,o,p,q);ae(c[q>>2]|0)|0;c[r>>2]=c[e>>2];e=c[o>>2]|0;o=c[p>>2]|0;c[k+0>>2]=c[r+0>>2];Bh(b,k,n,e,o,f,g);i=d;return}function Fh(b,d,e,f,g,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;j=+j;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0;d=i;i=i+304|0;l=d+8|0;m=d;n=d+272|0;o=d+264|0;p=d+36|0;q=d+268|0;r=d+24|0;s=d+20|0;t=d+28|0;u=d+32|0;v=m;c[v>>2]=37;c[v+4>>2]=0;v=m+1|0;w=f+4|0;x=c[w>>2]|0;if((x&2048|0)==0){y=v}else{a[v]=43;y=m+2|0}if((x&1024|0)==0){z=y}else{a[y]=35;z=y+1|0}y=x&260;v=x>>>14;do{if((y|0)==260){if((v&1|0)==0){a[z]=97;A=0;break}else{a[z]=65;A=0;break}}else{a[z]=46;x=z+2|0;a[z+1|0]=42;if((y|0)==256){if((v&1|0)==0){a[x]=101;A=1;break}else{a[x]=69;A=1;break}}else if((y|0)==4){if((v&1|0)==0){a[x]=102;A=1;break}else{a[x]=70;A=1;break}}else{if((v&1|0)==0){a[x]=103;A=1;break}else{a[x]=71;A=1;break}}}}while(0);c[o>>2]=n;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}v=c[2036]|0;if(A){c[l>>2]=c[f+8>>2];y=l+4|0;h[k>>3]=j;c[y>>2]=c[k>>2];c[y+4>>2]=c[k+4>>2];B=lh(n,30,v,m,l)|0}else{h[k>>3]=j;c[l>>2]=c[k>>2];c[l+4>>2]=c[k+4>>2];B=lh(n,30,v,m,l)|0}if((B|0)>29){v=(a[8152]|0)==0;if(A){if(v?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}A=c[2036]|0;c[l>>2]=c[f+8>>2];y=l+4|0;h[k>>3]=j;c[y>>2]=c[k>>2];c[y+4>>2]=c[k+4>>2];C=sh(o,A,m,l)|0}else{if(v?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}v=c[2036]|0;c[l>>2]=c[f+8>>2];A=l+4|0;h[k>>3]=j;c[A>>2]=c[k>>2];c[A+4>>2]=c[k+4>>2];C=sh(o,v,m,l)|0}m=c[o>>2]|0;if((m|0)==0){_m()}else{D=m;E=m;F=C}}else{D=c[o>>2]|0;E=0;F=B}B=D+F|0;o=c[w>>2]&176;do{if((o|0)==16){w=a[D]|0;if(w<<24>>24==43|w<<24>>24==45){G=D+1|0;break}if((F|0)>1&w<<24>>24==48?(w=a[D+1|0]|0,w<<24>>24==88|w<<24>>24==120):0){G=D+2|0}else{H=44}}else if((o|0)==32){G=B}else{H=44}}while(0);if((H|0)==44){G=D}if((D|0)!=(n|0)){H=Lm(F<<3)|0;if((H|0)==0){_m()}else{I=D;J=H;K=H}}else{I=n;J=0;K=p}Te(s,f);Gh(I,G,B,K,q,r,s);ae(c[s>>2]|0)|0;c[u>>2]=c[e>>2];s=c[q>>2]|0;q=c[r>>2]|0;c[l+0>>2]=c[u+0>>2];Bh(t,l,K,s,q,f,g);g=c[t>>2]|0;c[e>>2]=g;c[b>>2]=g;if((J|0)!=0){Mm(J)}if((E|0)==0){i=d;return}Mm(E);i=d;return}function Gh(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0;k=i;i=i+32|0;l=k;m=k+12|0;n=c[j>>2]|0;if(!((c[2062]|0)==-1)){c[l>>2]=8248;c[l+4>>2]=112;c[l+8>>2]=0;ue(8248,l,113)}o=(c[8252>>2]|0)+ -1|0;p=c[n+8>>2]|0;if(!((c[n+12>>2]|0)-p>>2>>>0>o>>>0)){q=Wa(4)|0;qm(q);xb(q|0,16216,101)}n=c[p+(o<<2)>>2]|0;if((n|0)==0){q=Wa(4)|0;qm(q);xb(q|0,16216,101)}q=c[j>>2]|0;if(!((c[2102]|0)==-1)){c[l>>2]=8408;c[l+4>>2]=112;c[l+8>>2]=0;ue(8408,l,113)}l=(c[8412>>2]|0)+ -1|0;j=c[q+8>>2]|0;if(!((c[q+12>>2]|0)-j>>2>>>0>l>>>0)){r=Wa(4)|0;qm(r);xb(r|0,16216,101)}q=c[j+(l<<2)>>2]|0;if((q|0)==0){r=Wa(4)|0;qm(r);xb(r|0,16216,101)}kc[c[(c[q>>2]|0)+20>>2]&63](m,q);c[h>>2]=f;r=a[b]|0;if(r<<24>>24==43|r<<24>>24==45){l=vc[c[(c[n>>2]|0)+44>>2]&15](n,r)|0;r=c[h>>2]|0;c[h>>2]=r+4;c[r>>2]=l;s=b+1|0}else{s=b}l=e;a:do{if(((l-s|0)>1?(a[s]|0)==48:0)?(r=s+1|0,j=a[r]|0,j<<24>>24==88|j<<24>>24==120):0){j=vc[c[(c[n>>2]|0)+44>>2]&15](n,48)|0;o=c[h>>2]|0;c[h>>2]=o+4;c[o>>2]=j;j=s+2|0;o=vc[c[(c[n>>2]|0)+44>>2]&15](n,a[r]|0)|0;r=c[h>>2]|0;c[h>>2]=r+4;c[r>>2]=o;if(j>>>0<e>>>0){o=j;while(1){r=a[o]|0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}p=o+1|0;if((Ya(r<<24>>24|0,c[2036]|0)|0)==0){t=j;u=o;break a}if(p>>>0<e>>>0){o=p}else{t=j;u=p;break}}}else{t=j;u=j}}else{v=14}}while(0);b:do{if((v|0)==14){if(s>>>0<e>>>0){o=s;while(1){p=a[o]|0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}r=o+1|0;if((ub(p<<24>>24|0,c[2036]|0)|0)==0){t=s;u=o;break b}if(r>>>0<e>>>0){o=r}else{t=s;u=r;break}}}else{t=s;u=s}}}while(0);s=a[m]|0;if((s&1)==0){w=(s&255)>>>1}else{w=c[m+4>>2]|0}if((w|0)!=0){if((t|0)!=(u|0)?(w=u+ -1|0,w>>>0>t>>>0):0){s=t;v=w;do{w=a[s]|0;a[s]=a[v]|0;a[v]=w;s=s+1|0;v=v+ -1|0}while(s>>>0<v>>>0)}v=mc[c[(c[q>>2]|0)+16>>2]&63](q)|0;if(t>>>0<u>>>0){s=m+1|0;w=m+4|0;o=m+8|0;j=0;r=0;p=t;while(1){x=(a[m]&1)==0;if((a[(x?s:c[o>>2]|0)+r|0]|0)>0?(j|0)==(a[(x?s:c[o>>2]|0)+r|0]|0):0){x=c[h>>2]|0;c[h>>2]=x+4;c[x>>2]=v;x=a[m]|0;if((x&1)==0){y=(x&255)>>>1}else{y=c[w>>2]|0}z=0;A=(r>>>0<(y+ -1|0)>>>0)+r|0}else{z=j;A=r}x=vc[c[(c[n>>2]|0)+44>>2]&15](n,a[p]|0)|0;B=c[h>>2]|0;C=B+4|0;c[h>>2]=C;c[B>>2]=x;x=p+1|0;if(x>>>0<u>>>0){j=z+1|0;r=A;p=x}else{D=C;break}}}else{D=c[h>>2]|0}p=f+(t-b<<2)|0;if((p|0)!=(D|0)?(A=D+ -4|0,A>>>0>p>>>0):0){r=p;p=A;while(1){A=c[r>>2]|0;c[r>>2]=c[p>>2];c[p>>2]=A;A=r+4|0;z=p+ -4|0;if(A>>>0<z>>>0){r=A;p=z}else{E=D;break}}}else{E=D}}else{sc[c[(c[n>>2]|0)+48>>2]&7](n,t,u,c[h>>2]|0)|0;D=(c[h>>2]|0)+(u-t<<2)|0;c[h>>2]=D;E=D}c:do{if(u>>>0<e>>>0){D=u;while(1){t=a[D]|0;if(t<<24>>24==46){break}p=vc[c[(c[n>>2]|0)+44>>2]&15](n,t)|0;t=c[h>>2]|0;r=t+4|0;c[h>>2]=r;c[t>>2]=p;p=D+1|0;if(p>>>0<e>>>0){D=p}else{F=r;G=p;break c}}p=mc[c[(c[q>>2]|0)+12>>2]&63](q)|0;r=c[h>>2]|0;t=r+4|0;c[h>>2]=t;c[r>>2]=p;F=t;G=D+1|0}else{F=E;G=u}}while(0);sc[c[(c[n>>2]|0)+48>>2]&7](n,G,e,F)|0;F=(c[h>>2]|0)+(l-G<<2)|0;c[h>>2]=F;if((d|0)==(e|0)){H=F;c[g>>2]=H;ze(m);i=k;return}H=f+(d-b<<2)|0;c[g>>2]=H;ze(m);i=k;return}function Hh(b,d,e,f,g,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;j=+j;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0;d=i;i=i+304|0;l=d+8|0;m=d;n=d+272|0;o=d+264|0;p=d+36|0;q=d+268|0;r=d+24|0;s=d+20|0;t=d+28|0;u=d+32|0;v=m;c[v>>2]=37;c[v+4>>2]=0;v=m+1|0;w=f+4|0;x=c[w>>2]|0;if((x&2048|0)==0){y=v}else{a[v]=43;y=m+2|0}if((x&1024|0)==0){z=y}else{a[y]=35;z=y+1|0}y=x&260;v=x>>>14;do{if((y|0)==260){a[z]=76;x=z+1|0;if((v&1|0)==0){a[x]=97;A=0;break}else{a[x]=65;A=0;break}}else{a[z]=46;a[z+1|0]=42;a[z+2|0]=76;x=z+3|0;if((y|0)==256){if((v&1|0)==0){a[x]=101;A=1;break}else{a[x]=69;A=1;break}}else if((y|0)==4){if((v&1|0)==0){a[x]=102;A=1;break}else{a[x]=70;A=1;break}}else{if((v&1|0)==0){a[x]=103;A=1;break}else{a[x]=71;A=1;break}}}}while(0);c[o>>2]=n;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}v=c[2036]|0;if(A){c[l>>2]=c[f+8>>2];y=l+4|0;h[k>>3]=j;c[y>>2]=c[k>>2];c[y+4>>2]=c[k+4>>2];B=lh(n,30,v,m,l)|0}else{h[k>>3]=j;c[l>>2]=c[k>>2];c[l+4>>2]=c[k+4>>2];B=lh(n,30,v,m,l)|0}if((B|0)>29){v=(a[8152]|0)==0;if(A){if(v?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}A=c[2036]|0;c[l>>2]=c[f+8>>2];y=l+4|0;h[k>>3]=j;c[y>>2]=c[k>>2];c[y+4>>2]=c[k+4>>2];C=sh(o,A,m,l)|0}else{if(v?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}v=c[2036]|0;h[k>>3]=j;c[l>>2]=c[k>>2];c[l+4>>2]=c[k+4>>2];C=sh(o,v,m,l)|0}m=c[o>>2]|0;if((m|0)==0){_m()}else{D=m;E=m;F=C}}else{D=c[o>>2]|0;E=0;F=B}B=D+F|0;o=c[w>>2]&176;do{if((o|0)==16){w=a[D]|0;if(w<<24>>24==43|w<<24>>24==45){G=D+1|0;break}if((F|0)>1&w<<24>>24==48?(w=a[D+1|0]|0,w<<24>>24==88|w<<24>>24==120):0){G=D+2|0}else{H=44}}else if((o|0)==32){G=B}else{H=44}}while(0);if((H|0)==44){G=D}if((D|0)!=(n|0)){H=Lm(F<<3)|0;if((H|0)==0){_m()}else{I=D;J=H;K=H}}else{I=n;J=0;K=p}Te(s,f);Gh(I,G,B,K,q,r,s);ae(c[s>>2]|0)|0;c[u>>2]=c[e>>2];s=c[q>>2]|0;q=c[r>>2]|0;c[l+0>>2]=c[u+0>>2];Bh(t,l,K,s,q,f,g);g=c[t>>2]|0;c[e>>2]=g;c[b>>2]=g;if((J|0)!=0){Mm(J)}if((E|0)==0){i=d;return}Mm(E);i=d;return}function Ih(b,d,e,f,g,h){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0;d=i;i=i+208|0;j=d;k=d+188|0;l=d+168|0;m=d+16|0;n=d+12|0;o=d+164|0;a[k+0|0]=a[7056|0]|0;a[k+1|0]=a[7057|0]|0;a[k+2|0]=a[7058|0]|0;a[k+3|0]=a[7059|0]|0;a[k+4|0]=a[7060|0]|0;a[k+5|0]=a[7061|0]|0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}p=c[2036]|0;c[j>>2]=h;h=lh(l,20,p,k,j)|0;k=l+h|0;p=c[f+4>>2]&176;do{if((p|0)==16){q=a[l]|0;if(q<<24>>24==43|q<<24>>24==45){r=l+1|0;break}if((h|0)>1&q<<24>>24==48?(q=a[l+1|0]|0,q<<24>>24==88|q<<24>>24==120):0){r=l+2|0}else{s=10}}else if((p|0)==32){r=k}else{s=10}}while(0);if((s|0)==10){r=l}Te(n,f);s=c[n>>2]|0;if(!((c[2062]|0)==-1)){c[j>>2]=8248;c[j+4>>2]=112;c[j+8>>2]=0;ue(8248,j,113)}p=(c[8252>>2]|0)+ -1|0;q=c[s+8>>2]|0;if(!((c[s+12>>2]|0)-q>>2>>>0>p>>>0)){t=Wa(4)|0;qm(t);xb(t|0,16216,101)}s=c[q+(p<<2)>>2]|0;if((s|0)==0){t=Wa(4)|0;qm(t);xb(t|0,16216,101)}ae(c[n>>2]|0)|0;sc[c[(c[s>>2]|0)+48>>2]&7](s,l,k,m)|0;s=m+(h<<2)|0;if((r|0)==(k|0)){u=s;v=c[e>>2]|0;c[o>>2]=v;c[j+0>>2]=c[o+0>>2];Bh(b,j,m,u,s,f,g);i=d;return}u=m+(r-l<<2)|0;v=c[e>>2]|0;c[o>>2]=v;c[j+0>>2]=c[o+0>>2];Bh(b,j,m,u,s,f,g);i=d;return}function Jh(e,f,g,h,j,k,l,m,n){e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;m=m|0;n=n|0;var o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0,X=0,Y=0;o=i;i=i+32|0;p=o;q=o+28|0;r=o+12|0;s=o+16|0;t=o+20|0;u=o+24|0;Te(r,j);v=c[r>>2]|0;if(!((c[2064]|0)==-1)){c[p>>2]=8256;c[p+4>>2]=112;c[p+8>>2]=0;ue(8256,p,113)}w=(c[8260>>2]|0)+ -1|0;x=c[v+8>>2]|0;if(!((c[v+12>>2]|0)-x>>2>>>0>w>>>0)){y=Wa(4)|0;qm(y);xb(y|0,16216,101)}v=c[x+(w<<2)>>2]|0;if((v|0)==0){y=Wa(4)|0;qm(y);xb(y|0,16216,101)}ae(c[r>>2]|0)|0;c[k>>2]=0;a:do{if((m|0)!=(n|0)){r=v+8|0;y=m;w=0;b:while(1){x=w;while(1){if((x|0)!=0){z=65;break a}A=c[g>>2]|0;if((A|0)!=0){if((c[A+12>>2]|0)==(c[A+16>>2]|0)?(mc[c[(c[A>>2]|0)+36>>2]&63](A)|0)==-1:0){c[g>>2]=0;B=0}else{B=A}}else{B=0}A=(B|0)==0;C=c[h>>2]|0;do{if((C|0)!=0){if((c[C+12>>2]|0)==(c[C+16>>2]|0)?(mc[c[(c[C>>2]|0)+36>>2]&63](C)|0)==-1:0){c[h>>2]=0;z=19;break}if(A){D=C}else{z=20;break b}}else{z=19}}while(0);if((z|0)==19){z=0;if(A){z=20;break b}else{D=0}}if((gc[c[(c[v>>2]|0)+36>>2]&31](v,a[y]|0,0)|0)<<24>>24==37){z=22;break}C=a[y]|0;if(C<<24>>24>-1?(E=c[r>>2]|0,!((b[E+(C<<24>>24<<1)>>1]&8192)==0)):0){F=y;z=33;break}G=B+12|0;C=c[G>>2]|0;H=B+16|0;if((C|0)==(c[H>>2]|0)){I=mc[c[(c[B>>2]|0)+36>>2]&63](B)|0}else{I=d[C]|0}C=vc[c[(c[v>>2]|0)+12>>2]&15](v,I&255)|0;if(C<<24>>24==(vc[c[(c[v>>2]|0)+12>>2]&15](v,a[y]|0)|0)<<24>>24){z=60;break}c[k>>2]=4;x=4}c:do{if((z|0)==22){z=0;x=y+1|0;if((x|0)==(n|0)){z=23;break b}C=gc[c[(c[v>>2]|0)+36>>2]&31](v,a[x]|0,0)|0;if(C<<24>>24==48|C<<24>>24==69){J=y+2|0;if((J|0)==(n|0)){z=26;break b}K=J;L=gc[c[(c[v>>2]|0)+36>>2]&31](v,a[J]|0,0)|0;M=C}else{K=x;L=C;M=0}C=c[(c[f>>2]|0)+36>>2]|0;c[t>>2]=B;c[u>>2]=D;c[q+0>>2]=c[t+0>>2];c[p+0>>2]=c[u+0>>2];lc[C&3](s,f,q,p,j,k,l,L,M);c[g>>2]=c[s>>2];N=K+1|0}else if((z|0)==33){while(1){z=0;C=F+1|0;if((C|0)==(n|0)){O=n;break}x=a[C]|0;if(!(x<<24>>24>-1)){O=C;break}if((b[E+(x<<24>>24<<1)>>1]&8192)==0){O=C;break}else{F=C;z=33}}A=B;C=D;x=D;while(1){if((A|0)!=0){if((c[A+12>>2]|0)==(c[A+16>>2]|0)?(mc[c[(c[A>>2]|0)+36>>2]&63](A)|0)==-1:0){c[g>>2]=0;P=0}else{P=A}}else{P=0}J=(P|0)==0;do{if((x|0)!=0){if((c[x+12>>2]|0)!=(c[x+16>>2]|0)){if(J){Q=C;R=x;break}else{N=O;break c}}if(!((mc[c[(c[x>>2]|0)+36>>2]&63](x)|0)==-1)){if(J^(C|0)==0){Q=C;R=C;break}else{N=O;break c}}else{c[h>>2]=0;S=0;z=46;break}}else{S=C;z=46}}while(0);if((z|0)==46){z=0;if(J){N=O;break c}else{Q=S;R=0}}T=P+12|0;U=c[T>>2]|0;V=P+16|0;if((U|0)==(c[V>>2]|0)){W=mc[c[(c[P>>2]|0)+36>>2]&63](P)|0}else{W=d[U]|0}if(!((W&255)<<24>>24>-1)){N=O;break c}if((b[(c[r>>2]|0)+(W<<24>>24<<1)>>1]&8192)==0){N=O;break c}U=c[T>>2]|0;if((U|0)==(c[V>>2]|0)){mc[c[(c[P>>2]|0)+40>>2]&63](P)|0;A=P;C=Q;x=R;continue}else{c[T>>2]=U+1;A=P;C=Q;x=R;continue}}}else if((z|0)==60){z=0;x=c[G>>2]|0;if((x|0)==(c[H>>2]|0)){mc[c[(c[B>>2]|0)+40>>2]&63](B)|0}else{c[G>>2]=x+1}N=y+1|0}}while(0);if((N|0)==(n|0)){z=65;break a}y=N;w=c[k>>2]|0}if((z|0)==20){c[k>>2]=4;X=B;break}else if((z|0)==23){c[k>>2]=4;X=B;break}else if((z|0)==26){c[k>>2]=4;X=B;break}}else{z=65}}while(0);if((z|0)==65){X=c[g>>2]|0}if((X|0)!=0){if((c[X+12>>2]|0)==(c[X+16>>2]|0)?(mc[c[(c[X>>2]|0)+36>>2]&63](X)|0)==-1:0){c[g>>2]=0;Y=0}else{Y=X}}else{Y=0}X=(Y|0)==0;g=c[h>>2]|0;do{if((g|0)!=0){if((c[g+12>>2]|0)==(c[g+16>>2]|0)?(mc[c[(c[g>>2]|0)+36>>2]&63](g)|0)==-1:0){c[h>>2]=0;z=75;break}if(X){c[e>>2]=Y;i=o;return}}else{z=75}}while(0);if((z|0)==75?!X:0){c[e>>2]=Y;i=o;return}c[k>>2]=c[k>>2]|2;c[e>>2]=Y;i=o;return}function Kh(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function Lh(a){a=a|0;return}function Mh(a){a=a|0;return 2}function Nh(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0;j=i;i=i+16|0;k=j+12|0;l=j;m=j+4|0;n=j+8|0;c[m>>2]=c[d>>2];c[n>>2]=c[e>>2];c[l+0>>2]=c[m+0>>2];c[k+0>>2]=c[n+0>>2];Jh(a,b,l,k,f,g,h,7160,7168|0);i=j;return}function Oh(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0;k=i;i=i+16|0;l=k+12|0;m=k;n=k+4|0;o=k+8|0;p=d+8|0;q=mc[c[(c[p>>2]|0)+20>>2]&63](p)|0;c[n>>2]=c[e>>2];c[o>>2]=c[f>>2];f=a[q]|0;if((f&1)==0){r=q+1|0;s=(f&255)>>>1;t=q+1|0}else{f=c[q+8>>2]|0;r=f;s=c[q+4>>2]|0;t=f}f=r+s|0;c[m+0>>2]=c[n+0>>2];c[l+0>>2]=c[o+0>>2];Jh(b,d,m,l,g,h,j,t,f);i=k;return}function Ph(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;j=i;i=i+32|0;k=j;l=j+16|0;m=j+12|0;Te(m,f);f=c[m>>2]|0;if(!((c[2064]|0)==-1)){c[k>>2]=8256;c[k+4>>2]=112;c[k+8>>2]=0;ue(8256,k,113)}n=(c[8260>>2]|0)+ -1|0;o=c[f+8>>2]|0;if(!((c[f+12>>2]|0)-o>>2>>>0>n>>>0)){p=Wa(4)|0;qm(p);xb(p|0,16216,101)}f=c[o+(n<<2)>>2]|0;if((f|0)==0){p=Wa(4)|0;qm(p);xb(p|0,16216,101)}ae(c[m>>2]|0)|0;m=c[e>>2]|0;e=b+8|0;b=mc[c[c[e>>2]>>2]&63](e)|0;c[l>>2]=m;m=b+168|0;c[k+0>>2]=c[l+0>>2];l=(jg(d,k,b,m,f,g,0)|0)-b|0;if((l|0)>=168){q=c[d>>2]|0;c[a>>2]=q;i=j;return}c[h+24>>2]=((l|0)/12|0|0)%7|0;q=c[d>>2]|0;c[a>>2]=q;i=j;return}function Qh(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;j=i;i=i+32|0;k=j;l=j+16|0;m=j+12|0;Te(m,f);f=c[m>>2]|0;if(!((c[2064]|0)==-1)){c[k>>2]=8256;c[k+4>>2]=112;c[k+8>>2]=0;ue(8256,k,113)}n=(c[8260>>2]|0)+ -1|0;o=c[f+8>>2]|0;if(!((c[f+12>>2]|0)-o>>2>>>0>n>>>0)){p=Wa(4)|0;qm(p);xb(p|0,16216,101)}f=c[o+(n<<2)>>2]|0;if((f|0)==0){p=Wa(4)|0;qm(p);xb(p|0,16216,101)}ae(c[m>>2]|0)|0;m=c[e>>2]|0;e=b+8|0;b=mc[c[(c[e>>2]|0)+4>>2]&63](e)|0;c[l>>2]=m;m=b+288|0;c[k+0>>2]=c[l+0>>2];l=(jg(d,k,b,m,f,g,0)|0)-b|0;if((l|0)>=288){q=c[d>>2]|0;c[a>>2]=q;i=j;return}c[h+16>>2]=((l|0)/12|0|0)%12|0;q=c[d>>2]|0;c[a>>2]=q;i=j;return}function Rh(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;b=i;i=i+32|0;j=b;k=b+16|0;l=b+12|0;Te(l,f);f=c[l>>2]|0;if(!((c[2064]|0)==-1)){c[j>>2]=8256;c[j+4>>2]=112;c[j+8>>2]=0;ue(8256,j,113)}m=(c[8260>>2]|0)+ -1|0;n=c[f+8>>2]|0;if(!((c[f+12>>2]|0)-n>>2>>>0>m>>>0)){o=Wa(4)|0;qm(o);xb(o|0,16216,101)}f=c[n+(m<<2)>>2]|0;if((f|0)==0){o=Wa(4)|0;qm(o);xb(o|0,16216,101)}ae(c[l>>2]|0)|0;l=h+20|0;c[k>>2]=c[e>>2];c[j+0>>2]=c[k+0>>2];k=Vh(d,j,g,f,4)|0;if((c[g>>2]&4|0)!=0){p=c[d>>2]|0;c[a>>2]=p;i=b;return}if((k|0)<69){q=k+2e3|0}else{q=(k+ -69|0)>>>0<31?k+1900|0:k}c[l>>2]=q+ -1900;p=c[d>>2]|0;c[a>>2]=p;i=b;return}function Sh(b,d,e,f,g,h,j,k,l){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;var m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0,X=0,Y=0,Z=0,_=0,$=0,aa=0,ba=0,ca=0,da=0,ea=0,fa=0,ga=0,ha=0,ia=0,ja=0,ka=0;l=i;i=i+176|0;m=l;n=l+160|0;o=l+156|0;p=l+152|0;q=l+148|0;r=l+144|0;s=l+140|0;t=l+136|0;u=l+132|0;v=l+128|0;w=l+124|0;x=l+120|0;y=l+116|0;z=l+112|0;A=l+108|0;B=l+104|0;C=l+100|0;D=l+96|0;E=l+92|0;F=l+88|0;G=l+164|0;H=l+44|0;I=l+36|0;J=l+32|0;K=l+28|0;L=l+40|0;M=l+16|0;N=l+12|0;O=l+20|0;P=l+24|0;Q=l+80|0;R=l+48|0;S=l+52|0;T=l+56|0;U=l+60|0;V=l+64|0;W=l+68|0;X=l+72|0;Y=l+76|0;Z=l+84|0;c[h>>2]=0;Te(A,g);_=c[A>>2]|0;if(!((c[2064]|0)==-1)){c[m>>2]=8256;c[m+4>>2]=112;c[m+8>>2]=0;ue(8256,m,113)}$=(c[8260>>2]|0)+ -1|0;aa=c[_+8>>2]|0;if(!((c[_+12>>2]|0)-aa>>2>>>0>$>>>0)){ba=Wa(4)|0;qm(ba);xb(ba|0,16216,101)}_=c[aa+($<<2)>>2]|0;if((_|0)==0){ba=Wa(4)|0;qm(ba);xb(ba|0,16216,101)}ae(c[A>>2]|0)|0;a:do{switch(k<<24>>24|0){case 89:{c[n>>2]=c[f>>2];c[m+0>>2]=c[n+0>>2];A=Vh(e,m,h,_,4)|0;if((c[h>>2]&4|0)==0){c[j+20>>2]=A+ -1900}break};case 37:{c[Z>>2]=c[f>>2];c[m+0>>2]=c[Z+0>>2];Uh(0,e,m,h,_);break};case 73:{A=j+8|0;c[v>>2]=c[f>>2];c[m+0>>2]=c[v+0>>2];ba=Vh(e,m,h,_,2)|0;$=c[h>>2]|0;if(($&4|0)==0?(ba+ -1|0)>>>0<12:0){c[A>>2]=ba;break a}c[h>>2]=$|4;break};case 109:{c[t>>2]=c[f>>2];c[m+0>>2]=c[t+0>>2];$=Vh(e,m,h,_,2)|0;ba=c[h>>2]|0;if((ba&4|0)==0&($|0)<13){c[j+16>>2]=$+ -1;break a}else{c[h>>2]=ba|4;break a}break};case 77:{c[s>>2]=c[f>>2];c[m+0>>2]=c[s+0>>2];ba=Vh(e,m,h,_,2)|0;$=c[h>>2]|0;if(($&4|0)==0&(ba|0)<60){c[j+4>>2]=ba;break a}else{c[h>>2]=$|4;break a}break};case 114:{c[M>>2]=c[e>>2];c[N>>2]=c[f>>2];c[n+0>>2]=c[M+0>>2];c[m+0>>2]=c[N+0>>2];Jh(L,d,n,m,g,h,j,7184,7195|0);c[e>>2]=c[L>>2];break};case 82:{c[P>>2]=c[e>>2];c[Q>>2]=c[f>>2];c[n+0>>2]=c[P+0>>2];c[m+0>>2]=c[Q+0>>2];Jh(O,d,n,m,g,h,j,7200,7205|0);c[e>>2]=c[O>>2];break};case 119:{c[p>>2]=c[f>>2];c[m+0>>2]=c[p+0>>2];$=Vh(e,m,h,_,1)|0;ba=c[h>>2]|0;if((ba&4|0)==0&($|0)<7){c[j+24>>2]=$;break a}else{c[h>>2]=ba|4;break a}break};case 120:{ba=c[(c[d>>2]|0)+20>>2]|0;c[U>>2]=c[e>>2];c[V>>2]=c[f>>2];c[n+0>>2]=c[U+0>>2];c[m+0>>2]=c[V+0>>2];hc[ba&63](b,d,n,m,g,h,j);i=l;return};case 88:{ba=d+8|0;$=mc[c[(c[ba>>2]|0)+24>>2]&63](ba)|0;c[X>>2]=c[e>>2];c[Y>>2]=c[f>>2];ba=a[$]|0;if((ba&1)==0){ca=$+1|0;da=(ba&255)>>>1;ea=$+1|0}else{ba=c[$+8>>2]|0;ca=ba;da=c[$+4>>2]|0;ea=ba}c[n+0>>2]=c[X+0>>2];c[m+0>>2]=c[Y+0>>2];Jh(W,d,n,m,g,h,j,ea,ca+da|0);c[e>>2]=c[W>>2];break};case 83:{c[q>>2]=c[f>>2];c[m+0>>2]=c[q+0>>2];ba=Vh(e,m,h,_,2)|0;$=c[h>>2]|0;if(($&4|0)==0&(ba|0)<61){c[j>>2]=ba;break a}else{c[h>>2]=$|4;break a}break};case 104:case 66:case 98:{$=c[f>>2]|0;ba=d+8|0;A=mc[c[(c[ba>>2]|0)+4>>2]&63](ba)|0;c[y>>2]=$;c[m+0>>2]=c[y+0>>2];$=(jg(e,m,A,A+288|0,_,h,0)|0)-A|0;if(($|0)<288){c[j+16>>2]=(($|0)/12|0|0)%12|0}break};case 106:{c[u>>2]=c[f>>2];c[m+0>>2]=c[u+0>>2];$=Vh(e,m,h,_,3)|0;A=c[h>>2]|0;if((A&4|0)==0&($|0)<366){c[j+28>>2]=$;break a}else{c[h>>2]=A|4;break a}break};case 84:{c[S>>2]=c[e>>2];c[T>>2]=c[f>>2];c[n+0>>2]=c[S+0>>2];c[m+0>>2]=c[T+0>>2];Jh(R,d,n,m,g,h,j,7208,7216|0);c[e>>2]=c[R>>2];break};case 65:case 97:{A=c[f>>2]|0;$=d+8|0;ba=mc[c[c[$>>2]>>2]&63]($)|0;c[z>>2]=A;c[m+0>>2]=c[z+0>>2];A=(jg(e,m,ba,ba+168|0,_,h,0)|0)-ba|0;if((A|0)<168){c[j+24>>2]=((A|0)/12|0|0)%7|0}break};case 99:{A=d+8|0;ba=mc[c[(c[A>>2]|0)+12>>2]&63](A)|0;c[C>>2]=c[e>>2];c[D>>2]=c[f>>2];A=a[ba]|0;if((A&1)==0){fa=ba+1|0;ga=(A&255)>>>1;ha=ba+1|0}else{A=c[ba+8>>2]|0;fa=A;ga=c[ba+4>>2]|0;ha=A}c[n+0>>2]=c[C+0>>2];c[m+0>>2]=c[D+0>>2];Jh(B,d,n,m,g,h,j,ha,fa+ga|0);c[e>>2]=c[B>>2];break};case 72:{c[w>>2]=c[f>>2];c[m+0>>2]=c[w+0>>2];A=Vh(e,m,h,_,2)|0;ba=c[h>>2]|0;if((ba&4|0)==0&(A|0)<24){c[j+8>>2]=A;break a}else{c[h>>2]=ba|4;break a}break};case 70:{c[I>>2]=c[e>>2];c[J>>2]=c[f>>2];c[n+0>>2]=c[I+0>>2];c[m+0>>2]=c[J+0>>2];Jh(H,d,n,m,g,h,j,7176,7184|0);c[e>>2]=c[H>>2];break};case 101:case 100:{ba=j+12|0;c[x>>2]=c[f>>2];c[m+0>>2]=c[x+0>>2];A=Vh(e,m,h,_,2)|0;$=c[h>>2]|0;if(($&4|0)==0?(A+ -1|0)>>>0<31:0){c[ba>>2]=A;break a}c[h>>2]=$|4;break};case 116:case 110:{c[K>>2]=c[f>>2];c[m+0>>2]=c[K+0>>2];Th(0,e,m,h,_);break};case 68:{c[F>>2]=c[e>>2];c[G>>2]=c[f>>2];c[n+0>>2]=c[F+0>>2];c[m+0>>2]=c[G+0>>2];Jh(E,d,n,m,g,h,j,7168,7176|0);c[e>>2]=c[E>>2];break};case 112:{$=j+8|0;A=c[f>>2]|0;ba=d+8|0;aa=mc[c[(c[ba>>2]|0)+8>>2]&63](ba)|0;ba=a[aa]|0;if((ba&1)==0){ia=(ba&255)>>>1}else{ia=c[aa+4>>2]|0}ba=a[aa+12|0]|0;if((ba&1)==0){ja=(ba&255)>>>1}else{ja=c[aa+16>>2]|0}if((ia|0)==(0-ja|0)){c[h>>2]=c[h>>2]|4;break a}c[r>>2]=A;c[m+0>>2]=c[r+0>>2];A=jg(e,m,aa,aa+24|0,_,h,0)|0;ba=A-aa|0;if((A|0)==(aa|0)?(c[$>>2]|0)==12:0){c[$>>2]=0;break a}if((ba|0)==12?(ba=c[$>>2]|0,(ba|0)<12):0){c[$>>2]=ba+12}break};case 121:{ba=j+20|0;c[o>>2]=c[f>>2];c[m+0>>2]=c[o+0>>2];$=Vh(e,m,h,_,4)|0;if((c[h>>2]&4|0)==0){if(($|0)<69){ka=$+2e3|0}else{ka=($+ -69|0)>>>0<31?$+1900|0:$}c[ba>>2]=ka+ -1900}break};default:{c[h>>2]=c[h>>2]|4}}}while(0);c[b>>2]=c[e>>2];i=l;return}function Th(a,e,f,g,h){a=a|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0;a=i;j=h+8|0;a:while(1){h=c[e>>2]|0;do{if((h|0)!=0){if((c[h+12>>2]|0)==(c[h+16>>2]|0)){if((mc[c[(c[h>>2]|0)+36>>2]&63](h)|0)==-1){c[e>>2]=0;k=0;break}else{k=c[e>>2]|0;break}}else{k=h}}else{k=0}}while(0);h=(k|0)==0;l=c[f>>2]|0;do{if((l|0)!=0){if((c[l+12>>2]|0)!=(c[l+16>>2]|0)){if(h){m=l;break}else{n=l;break a}}if(!((mc[c[(c[l>>2]|0)+36>>2]&63](l)|0)==-1)){if(h){m=l;break}else{n=l;break a}}else{c[f>>2]=0;o=12;break}}else{o=12}}while(0);if((o|0)==12){o=0;if(h){n=0;break}else{m=0}}l=c[e>>2]|0;p=c[l+12>>2]|0;if((p|0)==(c[l+16>>2]|0)){q=mc[c[(c[l>>2]|0)+36>>2]&63](l)|0}else{q=d[p]|0}if(!((q&255)<<24>>24>-1)){n=m;break}if((b[(c[j>>2]|0)+(q<<24>>24<<1)>>1]&8192)==0){n=m;break}p=c[e>>2]|0;l=p+12|0;r=c[l>>2]|0;if((r|0)==(c[p+16>>2]|0)){mc[c[(c[p>>2]|0)+40>>2]&63](p)|0;continue}else{c[l>>2]=r+1;continue}}m=c[e>>2]|0;do{if((m|0)!=0){if((c[m+12>>2]|0)==(c[m+16>>2]|0)){if((mc[c[(c[m>>2]|0)+36>>2]&63](m)|0)==-1){c[e>>2]=0;s=0;break}else{s=c[e>>2]|0;break}}else{s=m}}else{s=0}}while(0);m=(s|0)==0;do{if((n|0)!=0){if((c[n+12>>2]|0)==(c[n+16>>2]|0)?(mc[c[(c[n>>2]|0)+36>>2]&63](n)|0)==-1:0){c[f>>2]=0;o=32;break}if(m){i=a;return}}else{o=32}}while(0);if((o|0)==32?!m:0){i=a;return}c[g>>2]=c[g>>2]|2;i=a;return}function Uh(a,b,e,f,g){a=a|0;b=b|0;e=e|0;f=f|0;g=g|0;var h=0,j=0,k=0,l=0,m=0,n=0;a=i;h=c[b>>2]|0;do{if((h|0)!=0){if((c[h+12>>2]|0)==(c[h+16>>2]|0)){if((mc[c[(c[h>>2]|0)+36>>2]&63](h)|0)==-1){c[b>>2]=0;j=0;break}else{j=c[b>>2]|0;break}}else{j=h}}else{j=0}}while(0);h=(j|0)==0;j=c[e>>2]|0;do{if((j|0)!=0){if((c[j+12>>2]|0)==(c[j+16>>2]|0)?(mc[c[(c[j>>2]|0)+36>>2]&63](j)|0)==-1:0){c[e>>2]=0;k=11;break}if(h){l=j}else{k=12}}else{k=11}}while(0);if((k|0)==11){if(h){k=12}else{l=0}}if((k|0)==12){c[f>>2]=c[f>>2]|6;i=a;return}h=c[b>>2]|0;j=c[h+12>>2]|0;if((j|0)==(c[h+16>>2]|0)){m=mc[c[(c[h>>2]|0)+36>>2]&63](h)|0}else{m=d[j]|0}if(!((gc[c[(c[g>>2]|0)+36>>2]&31](g,m&255,0)|0)<<24>>24==37)){c[f>>2]=c[f>>2]|4;i=a;return}m=c[b>>2]|0;g=m+12|0;j=c[g>>2]|0;if((j|0)==(c[m+16>>2]|0)){mc[c[(c[m>>2]|0)+40>>2]&63](m)|0}else{c[g>>2]=j+1}j=c[b>>2]|0;do{if((j|0)!=0){if((c[j+12>>2]|0)==(c[j+16>>2]|0)){if((mc[c[(c[j>>2]|0)+36>>2]&63](j)|0)==-1){c[b>>2]=0;n=0;break}else{n=c[b>>2]|0;break}}else{n=j}}else{n=0}}while(0);j=(n|0)==0;do{if((l|0)!=0){if((c[l+12>>2]|0)==(c[l+16>>2]|0)?(mc[c[(c[l>>2]|0)+36>>2]&63](l)|0)==-1:0){c[e>>2]=0;k=31;break}if(j){i=a;return}}else{k=31}}while(0);if((k|0)==31?!j:0){i=a;return}c[f>>2]=c[f>>2]|2;i=a;return}function Vh(a,e,f,g,h){a=a|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0;j=i;k=c[a>>2]|0;do{if((k|0)!=0){if((c[k+12>>2]|0)==(c[k+16>>2]|0)){if((mc[c[(c[k>>2]|0)+36>>2]&63](k)|0)==-1){c[a>>2]=0;l=0;break}else{l=c[a>>2]|0;break}}else{l=k}}else{l=0}}while(0);k=(l|0)==0;l=c[e>>2]|0;do{if((l|0)!=0){if((c[l+12>>2]|0)==(c[l+16>>2]|0)?(mc[c[(c[l>>2]|0)+36>>2]&63](l)|0)==-1:0){c[e>>2]=0;m=11;break}if(k){n=l}else{m=12}}else{m=11}}while(0);if((m|0)==11){if(k){m=12}else{n=0}}if((m|0)==12){c[f>>2]=c[f>>2]|6;o=0;i=j;return o|0}k=c[a>>2]|0;l=c[k+12>>2]|0;if((l|0)==(c[k+16>>2]|0)){p=mc[c[(c[k>>2]|0)+36>>2]&63](k)|0}else{p=d[l]|0}l=p&255;if(l<<24>>24>-1?(k=g+8|0,!((b[(c[k>>2]|0)+(p<<24>>24<<1)>>1]&2048)==0)):0){p=(gc[c[(c[g>>2]|0)+36>>2]&31](g,l,0)|0)<<24>>24;l=c[a>>2]|0;q=l+12|0;r=c[q>>2]|0;if((r|0)==(c[l+16>>2]|0)){mc[c[(c[l>>2]|0)+40>>2]&63](l)|0;s=h;t=n;u=n;v=p}else{c[q>>2]=r+1;s=h;t=n;u=n;v=p}while(1){w=v+ -48|0;p=s+ -1|0;n=c[a>>2]|0;do{if((n|0)!=0){if((c[n+12>>2]|0)==(c[n+16>>2]|0)){if((mc[c[(c[n>>2]|0)+36>>2]&63](n)|0)==-1){c[a>>2]=0;x=0;break}else{x=c[a>>2]|0;break}}else{x=n}}else{x=0}}while(0);n=(x|0)==0;if((u|0)!=0){if((c[u+12>>2]|0)==(c[u+16>>2]|0)){if((mc[c[(c[u>>2]|0)+36>>2]&63](u)|0)==-1){c[e>>2]=0;y=0;z=0}else{y=t;z=t}}else{y=t;z=u}}else{y=t;z=0}A=c[a>>2]|0;if(!((n^(z|0)==0)&(p|0)>0)){m=40;break}n=c[A+12>>2]|0;if((n|0)==(c[A+16>>2]|0)){B=mc[c[(c[A>>2]|0)+36>>2]&63](A)|0}else{B=d[n]|0}n=B&255;if(!(n<<24>>24>-1)){o=w;m=52;break}if((b[(c[k>>2]|0)+(B<<24>>24<<1)>>1]&2048)==0){o=w;m=52;break}h=((gc[c[(c[g>>2]|0)+36>>2]&31](g,n,0)|0)<<24>>24)+(w*10|0)|0;n=c[a>>2]|0;r=n+12|0;q=c[r>>2]|0;if((q|0)==(c[n+16>>2]|0)){mc[c[(c[n>>2]|0)+40>>2]&63](n)|0;s=p;t=y;u=z;v=h;continue}else{c[r>>2]=q+1;s=p;t=y;u=z;v=h;continue}}if((m|0)==40){do{if((A|0)!=0){if((c[A+12>>2]|0)==(c[A+16>>2]|0)){if((mc[c[(c[A>>2]|0)+36>>2]&63](A)|0)==-1){c[a>>2]=0;C=0;break}else{C=c[a>>2]|0;break}}else{C=A}}else{C=0}}while(0);A=(C|0)==0;do{if((y|0)!=0){if((c[y+12>>2]|0)==(c[y+16>>2]|0)?(mc[c[(c[y>>2]|0)+36>>2]&63](y)|0)==-1:0){c[e>>2]=0;m=50;break}if(A){o=w;i=j;return o|0}}else{m=50}}while(0);if((m|0)==50?!A:0){o=w;i=j;return o|0}c[f>>2]=c[f>>2]|2;o=w;i=j;return o|0}else if((m|0)==52){i=j;return o|0}}c[f>>2]=c[f>>2]|4;o=0;i=j;return o|0}function Wh(a,b,d,e,f,g,h,j,k){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0,X=0,Y=0,Z=0,_=0,$=0;l=i;i=i+32|0;m=l;n=l+28|0;o=l+12|0;p=l+16|0;q=l+20|0;r=l+24|0;Te(o,f);s=c[o>>2]|0;if(!((c[2062]|0)==-1)){c[m>>2]=8248;c[m+4>>2]=112;c[m+8>>2]=0;ue(8248,m,113)}t=(c[8252>>2]|0)+ -1|0;u=c[s+8>>2]|0;if(!((c[s+12>>2]|0)-u>>2>>>0>t>>>0)){v=Wa(4)|0;qm(v);xb(v|0,16216,101)}s=c[u+(t<<2)>>2]|0;if((s|0)==0){v=Wa(4)|0;qm(v);xb(v|0,16216,101)}ae(c[o>>2]|0)|0;c[g>>2]=0;a:do{if((j|0)!=(k|0)){o=j;v=0;b:while(1){t=v;while(1){if((t|0)!=0){w=69;break a}u=c[d>>2]|0;if((u|0)!=0){x=c[u+12>>2]|0;if((x|0)==(c[u+16>>2]|0)){y=mc[c[(c[u>>2]|0)+36>>2]&63](u)|0}else{y=c[x>>2]|0}if((y|0)==-1){c[d>>2]=0;z=1;A=0}else{z=0;A=u}}else{z=1;A=0}u=c[e>>2]|0;do{if((u|0)!=0){x=c[u+12>>2]|0;if((x|0)==(c[u+16>>2]|0)){B=mc[c[(c[u>>2]|0)+36>>2]&63](u)|0}else{B=c[x>>2]|0}if(!((B|0)==-1)){if(z){C=u;break}else{w=24;break b}}else{c[e>>2]=0;w=22;break}}else{w=22}}while(0);if((w|0)==22){w=0;if(z){w=24;break b}else{C=0}}if((gc[c[(c[s>>2]|0)+52>>2]&31](s,c[o>>2]|0,0)|0)<<24>>24==37){w=26;break}if(gc[c[(c[s>>2]|0)+12>>2]&31](s,8192,c[o>>2]|0)|0){D=o;w=36;break}E=A+12|0;u=c[E>>2]|0;F=A+16|0;if((u|0)==(c[F>>2]|0)){G=mc[c[(c[A>>2]|0)+36>>2]&63](A)|0}else{G=c[u>>2]|0}u=vc[c[(c[s>>2]|0)+28>>2]&15](s,G)|0;if((u|0)==(vc[c[(c[s>>2]|0)+28>>2]&15](s,c[o>>2]|0)|0)){w=64;break}c[g>>2]=4;t=4}c:do{if((w|0)==26){w=0;t=o+4|0;if((t|0)==(k|0)){w=27;break b}u=gc[c[(c[s>>2]|0)+52>>2]&31](s,c[t>>2]|0,0)|0;if(u<<24>>24==48|u<<24>>24==69){x=o+8|0;if((x|0)==(k|0)){w=30;break b}H=x;I=gc[c[(c[s>>2]|0)+52>>2]&31](s,c[x>>2]|0,0)|0;J=u}else{H=t;I=u;J=0}u=c[(c[b>>2]|0)+36>>2]|0;c[q>>2]=A;c[r>>2]=C;c[n+0>>2]=c[q+0>>2];c[m+0>>2]=c[r+0>>2];lc[u&3](p,b,n,m,f,g,h,I,J);c[d>>2]=c[p>>2];K=H+4|0}else if((w|0)==36){while(1){w=0;u=D+4|0;if((u|0)==(k|0)){L=k;break}if(gc[c[(c[s>>2]|0)+12>>2]&31](s,8192,c[u>>2]|0)|0){D=u;w=36}else{L=u;break}}u=A;t=C;x=C;while(1){if((u|0)!=0){M=c[u+12>>2]|0;if((M|0)==(c[u+16>>2]|0)){N=mc[c[(c[u>>2]|0)+36>>2]&63](u)|0}else{N=c[M>>2]|0}if((N|0)==-1){c[d>>2]=0;O=1;P=0}else{O=0;P=u}}else{O=1;P=0}do{if((x|0)!=0){M=c[x+12>>2]|0;if((M|0)==(c[x+16>>2]|0)){Q=mc[c[(c[x>>2]|0)+36>>2]&63](x)|0}else{Q=c[M>>2]|0}if(!((Q|0)==-1)){if(O^(t|0)==0){R=t;S=t;break}else{K=L;break c}}else{c[e>>2]=0;T=0;w=51;break}}else{T=t;w=51}}while(0);if((w|0)==51){w=0;if(O){K=L;break c}else{R=T;S=0}}M=P+12|0;U=c[M>>2]|0;V=P+16|0;if((U|0)==(c[V>>2]|0)){W=mc[c[(c[P>>2]|0)+36>>2]&63](P)|0}else{W=c[U>>2]|0}if(!(gc[c[(c[s>>2]|0)+12>>2]&31](s,8192,W)|0)){K=L;break c}U=c[M>>2]|0;if((U|0)==(c[V>>2]|0)){mc[c[(c[P>>2]|0)+40>>2]&63](P)|0;u=P;t=R;x=S;continue}else{c[M>>2]=U+4;u=P;t=R;x=S;continue}}}else if((w|0)==64){w=0;x=c[E>>2]|0;if((x|0)==(c[F>>2]|0)){mc[c[(c[A>>2]|0)+40>>2]&63](A)|0}else{c[E>>2]=x+4}K=o+4|0}}while(0);if((K|0)==(k|0)){w=69;break a}o=K;v=c[g>>2]|0}if((w|0)==24){c[g>>2]=4;X=A;break}else if((w|0)==27){c[g>>2]=4;X=A;break}else if((w|0)==30){c[g>>2]=4;X=A;break}}else{w=69}}while(0);if((w|0)==69){X=c[d>>2]|0}if((X|0)!=0){A=c[X+12>>2]|0;if((A|0)==(c[X+16>>2]|0)){Y=mc[c[(c[X>>2]|0)+36>>2]&63](X)|0}else{Y=c[A>>2]|0}if((Y|0)==-1){c[d>>2]=0;Z=0;_=1}else{Z=X;_=0}}else{Z=0;_=1}X=c[e>>2]|0;do{if((X|0)!=0){d=c[X+12>>2]|0;if((d|0)==(c[X+16>>2]|0)){$=mc[c[(c[X>>2]|0)+36>>2]&63](X)|0}else{$=c[d>>2]|0}if(($|0)==-1){c[e>>2]=0;w=82;break}if(_){c[a>>2]=Z;i=l;return}}else{w=82}}while(0);if((w|0)==82?!_:0){c[a>>2]=Z;i=l;return}c[g>>2]=c[g>>2]|2;c[a>>2]=Z;i=l;return}function Xh(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function Yh(a){a=a|0;return}function Zh(a){a=a|0;return 2}function _h(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0;j=i;i=i+16|0;k=j+12|0;l=j;m=j+4|0;n=j+8|0;c[m>>2]=c[d>>2];c[n>>2]=c[e>>2];c[l+0>>2]=c[m+0>>2];c[k+0>>2]=c[n+0>>2];Wh(a,b,l,k,f,g,h,7312,7344|0);i=j;return}function $h(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0;k=i;i=i+16|0;l=k+12|0;m=k;n=k+4|0;o=k+8|0;p=d+8|0;q=mc[c[(c[p>>2]|0)+20>>2]&63](p)|0;c[n>>2]=c[e>>2];c[o>>2]=c[f>>2];f=a[q]|0;if((f&1)==0){r=q+4|0;s=(f&255)>>>1;t=q+4|0}else{f=c[q+8>>2]|0;r=f;s=c[q+4>>2]|0;t=f}f=r+(s<<2)|0;c[m+0>>2]=c[n+0>>2];c[l+0>>2]=c[o+0>>2];Wh(b,d,m,l,g,h,j,t,f);i=k;return}function ai(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;j=i;i=i+32|0;k=j;l=j+16|0;m=j+12|0;Te(m,f);f=c[m>>2]|0;if(!((c[2062]|0)==-1)){c[k>>2]=8248;c[k+4>>2]=112;c[k+8>>2]=0;ue(8248,k,113)}n=(c[8252>>2]|0)+ -1|0;o=c[f+8>>2]|0;if(!((c[f+12>>2]|0)-o>>2>>>0>n>>>0)){p=Wa(4)|0;qm(p);xb(p|0,16216,101)}f=c[o+(n<<2)>>2]|0;if((f|0)==0){p=Wa(4)|0;qm(p);xb(p|0,16216,101)}ae(c[m>>2]|0)|0;m=c[e>>2]|0;e=b+8|0;b=mc[c[c[e>>2]>>2]&63](e)|0;c[l>>2]=m;m=b+168|0;c[k+0>>2]=c[l+0>>2];l=(Ig(d,k,b,m,f,g,0)|0)-b|0;if((l|0)>=168){q=c[d>>2]|0;c[a>>2]=q;i=j;return}c[h+24>>2]=((l|0)/12|0|0)%7|0;q=c[d>>2]|0;c[a>>2]=q;i=j;return}function bi(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;j=i;i=i+32|0;k=j;l=j+16|0;m=j+12|0;Te(m,f);f=c[m>>2]|0;if(!((c[2062]|0)==-1)){c[k>>2]=8248;c[k+4>>2]=112;c[k+8>>2]=0;ue(8248,k,113)}n=(c[8252>>2]|0)+ -1|0;o=c[f+8>>2]|0;if(!((c[f+12>>2]|0)-o>>2>>>0>n>>>0)){p=Wa(4)|0;qm(p);xb(p|0,16216,101)}f=c[o+(n<<2)>>2]|0;if((f|0)==0){p=Wa(4)|0;qm(p);xb(p|0,16216,101)}ae(c[m>>2]|0)|0;m=c[e>>2]|0;e=b+8|0;b=mc[c[(c[e>>2]|0)+4>>2]&63](e)|0;c[l>>2]=m;m=b+288|0;c[k+0>>2]=c[l+0>>2];l=(Ig(d,k,b,m,f,g,0)|0)-b|0;if((l|0)>=288){q=c[d>>2]|0;c[a>>2]=q;i=j;return}c[h+16>>2]=((l|0)/12|0|0)%12|0;q=c[d>>2]|0;c[a>>2]=q;i=j;return}function ci(a,b,d,e,f,g,h){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;b=i;i=i+32|0;j=b;k=b+16|0;l=b+12|0;Te(l,f);f=c[l>>2]|0;if(!((c[2062]|0)==-1)){c[j>>2]=8248;c[j+4>>2]=112;c[j+8>>2]=0;ue(8248,j,113)}m=(c[8252>>2]|0)+ -1|0;n=c[f+8>>2]|0;if(!((c[f+12>>2]|0)-n>>2>>>0>m>>>0)){o=Wa(4)|0;qm(o);xb(o|0,16216,101)}f=c[n+(m<<2)>>2]|0;if((f|0)==0){o=Wa(4)|0;qm(o);xb(o|0,16216,101)}ae(c[l>>2]|0)|0;l=h+20|0;c[k>>2]=c[e>>2];c[j+0>>2]=c[k+0>>2];k=gi(d,j,g,f,4)|0;if((c[g>>2]&4|0)!=0){p=c[d>>2]|0;c[a>>2]=p;i=b;return}if((k|0)<69){q=k+2e3|0}else{q=(k+ -69|0)>>>0<31?k+1900|0:k}c[l>>2]=q+ -1900;p=c[d>>2]|0;c[a>>2]=p;i=b;return}function di(b,d,e,f,g,h,j,k,l){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;var m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0,X=0,Y=0,Z=0,_=0,$=0,aa=0,ba=0,ca=0,da=0,ea=0,fa=0,ga=0,ha=0,ia=0,ja=0,ka=0;l=i;i=i+176|0;m=l;n=l+160|0;o=l+156|0;p=l+152|0;q=l+148|0;r=l+144|0;s=l+140|0;t=l+136|0;u=l+132|0;v=l+128|0;w=l+124|0;x=l+120|0;y=l+116|0;z=l+112|0;A=l+108|0;B=l+104|0;C=l+100|0;D=l+96|0;E=l+92|0;F=l+88|0;G=l+164|0;H=l+44|0;I=l+36|0;J=l+32|0;K=l+28|0;L=l+40|0;M=l+16|0;N=l+12|0;O=l+20|0;P=l+24|0;Q=l+80|0;R=l+48|0;S=l+52|0;T=l+56|0;U=l+60|0;V=l+64|0;W=l+68|0;X=l+72|0;Y=l+76|0;Z=l+84|0;c[h>>2]=0;Te(A,g);_=c[A>>2]|0;if(!((c[2062]|0)==-1)){c[m>>2]=8248;c[m+4>>2]=112;c[m+8>>2]=0;ue(8248,m,113)}$=(c[8252>>2]|0)+ -1|0;aa=c[_+8>>2]|0;if(!((c[_+12>>2]|0)-aa>>2>>>0>$>>>0)){ba=Wa(4)|0;qm(ba);xb(ba|0,16216,101)}_=c[aa+($<<2)>>2]|0;if((_|0)==0){ba=Wa(4)|0;qm(ba);xb(ba|0,16216,101)}ae(c[A>>2]|0)|0;a:do{switch(k<<24>>24|0){case 116:case 110:{c[K>>2]=c[f>>2];c[m+0>>2]=c[K+0>>2];ei(0,e,m,h,_);break};case 72:{c[w>>2]=c[f>>2];c[m+0>>2]=c[w+0>>2];A=gi(e,m,h,_,2)|0;ba=c[h>>2]|0;if((ba&4|0)==0&(A|0)<24){c[j+8>>2]=A;break a}else{c[h>>2]=ba|4;break a}break};case 73:{ba=j+8|0;c[v>>2]=c[f>>2];c[m+0>>2]=c[v+0>>2];A=gi(e,m,h,_,2)|0;$=c[h>>2]|0;if(($&4|0)==0?(A+ -1|0)>>>0<12:0){c[ba>>2]=A;break a}c[h>>2]=$|4;break};case 104:case 66:case 98:{$=c[f>>2]|0;A=d+8|0;ba=mc[c[(c[A>>2]|0)+4>>2]&63](A)|0;c[y>>2]=$;c[m+0>>2]=c[y+0>>2];$=(Ig(e,m,ba,ba+288|0,_,h,0)|0)-ba|0;if(($|0)<288){c[j+16>>2]=(($|0)/12|0|0)%12|0}break};case 106:{c[u>>2]=c[f>>2];c[m+0>>2]=c[u+0>>2];$=gi(e,m,h,_,3)|0;ba=c[h>>2]|0;if((ba&4|0)==0&($|0)<366){c[j+28>>2]=$;break a}else{c[h>>2]=ba|4;break a}break};case 77:{c[s>>2]=c[f>>2];c[m+0>>2]=c[s+0>>2];ba=gi(e,m,h,_,2)|0;$=c[h>>2]|0;if(($&4|0)==0&(ba|0)<60){c[j+4>>2]=ba;break a}else{c[h>>2]=$|4;break a}break};case 99:{$=d+8|0;ba=mc[c[(c[$>>2]|0)+12>>2]&63]($)|0;c[C>>2]=c[e>>2];c[D>>2]=c[f>>2];$=a[ba]|0;if(($&1)==0){ca=ba+4|0;da=($&255)>>>1;ea=ba+4|0}else{$=c[ba+8>>2]|0;ca=$;da=c[ba+4>>2]|0;ea=$}c[n+0>>2]=c[C+0>>2];c[m+0>>2]=c[D+0>>2];Wh(B,d,n,m,g,h,j,ea,ca+(da<<2)|0);c[e>>2]=c[B>>2];break};case 82:{c[P>>2]=c[e>>2];c[Q>>2]=c[f>>2];c[n+0>>2]=c[P+0>>2];c[m+0>>2]=c[Q+0>>2];Wh(O,d,n,m,g,h,j,7456,7476|0);c[e>>2]=c[O>>2];break};case 83:{c[q>>2]=c[f>>2];c[m+0>>2]=c[q+0>>2];$=gi(e,m,h,_,2)|0;ba=c[h>>2]|0;if((ba&4|0)==0&($|0)<61){c[j>>2]=$;break a}else{c[h>>2]=ba|4;break a}break};case 109:{c[t>>2]=c[f>>2];c[m+0>>2]=c[t+0>>2];ba=gi(e,m,h,_,2)|0;$=c[h>>2]|0;if(($&4|0)==0&(ba|0)<13){c[j+16>>2]=ba+ -1;break a}else{c[h>>2]=$|4;break a}break};case 84:{c[S>>2]=c[e>>2];c[T>>2]=c[f>>2];c[n+0>>2]=c[S+0>>2];c[m+0>>2]=c[T+0>>2];Wh(R,d,n,m,g,h,j,7480,7512|0);c[e>>2]=c[R>>2];break};case 119:{c[p>>2]=c[f>>2];c[m+0>>2]=c[p+0>>2];$=gi(e,m,h,_,1)|0;ba=c[h>>2]|0;if((ba&4|0)==0&($|0)<7){c[j+24>>2]=$;break a}else{c[h>>2]=ba|4;break a}break};case 68:{c[F>>2]=c[e>>2];c[G>>2]=c[f>>2];c[n+0>>2]=c[F+0>>2];c[m+0>>2]=c[G+0>>2];Wh(E,d,n,m,g,h,j,7344,7376|0);c[e>>2]=c[E>>2];break};case 70:{c[I>>2]=c[e>>2];c[J>>2]=c[f>>2];c[n+0>>2]=c[I+0>>2];c[m+0>>2]=c[J+0>>2];Wh(H,d,n,m,g,h,j,7376,7408|0);c[e>>2]=c[H>>2];break};case 101:case 100:{ba=j+12|0;c[x>>2]=c[f>>2];c[m+0>>2]=c[x+0>>2];$=gi(e,m,h,_,2)|0;A=c[h>>2]|0;if((A&4|0)==0?($+ -1|0)>>>0<31:0){c[ba>>2]=$;break a}c[h>>2]=A|4;break};case 114:{c[M>>2]=c[e>>2];c[N>>2]=c[f>>2];c[n+0>>2]=c[M+0>>2];c[m+0>>2]=c[N+0>>2];Wh(L,d,n,m,g,h,j,7408,7452|0);c[e>>2]=c[L>>2];break};case 112:{A=j+8|0;$=c[f>>2]|0;ba=d+8|0;aa=mc[c[(c[ba>>2]|0)+8>>2]&63](ba)|0;ba=a[aa]|0;if((ba&1)==0){fa=(ba&255)>>>1}else{fa=c[aa+4>>2]|0}ba=a[aa+12|0]|0;if((ba&1)==0){ga=(ba&255)>>>1}else{ga=c[aa+16>>2]|0}if((fa|0)==(0-ga|0)){c[h>>2]=c[h>>2]|4;break a}c[r>>2]=$;c[m+0>>2]=c[r+0>>2];$=Ig(e,m,aa,aa+24|0,_,h,0)|0;ba=$-aa|0;if(($|0)==(aa|0)?(c[A>>2]|0)==12:0){c[A>>2]=0;break a}if((ba|0)==12?(ba=c[A>>2]|0,(ba|0)<12):0){c[A>>2]=ba+12}break};case 65:case 97:{ba=c[f>>2]|0;A=d+8|0;aa=mc[c[c[A>>2]>>2]&63](A)|0;c[z>>2]=ba;c[m+0>>2]=c[z+0>>2];ba=(Ig(e,m,aa,aa+168|0,_,h,0)|0)-aa|0;if((ba|0)<168){c[j+24>>2]=((ba|0)/12|0|0)%7|0}break};case 121:{ba=j+20|0;c[o>>2]=c[f>>2];c[m+0>>2]=c[o+0>>2];aa=gi(e,m,h,_,4)|0;if((c[h>>2]&4|0)==0){if((aa|0)<69){ha=aa+2e3|0}else{ha=(aa+ -69|0)>>>0<31?aa+1900|0:aa}c[ba>>2]=ha+ -1900}break};case 120:{ba=c[(c[d>>2]|0)+20>>2]|0;c[U>>2]=c[e>>2];c[V>>2]=c[f>>2];c[n+0>>2]=c[U+0>>2];c[m+0>>2]=c[V+0>>2];hc[ba&63](b,d,n,m,g,h,j);i=l;return};case 89:{c[n>>2]=c[f>>2];c[m+0>>2]=c[n+0>>2];ba=gi(e,m,h,_,4)|0;if((c[h>>2]&4|0)==0){c[j+20>>2]=ba+ -1900}break};case 37:{c[Z>>2]=c[f>>2];c[m+0>>2]=c[Z+0>>2];fi(0,e,m,h,_);break};case 88:{ba=d+8|0;aa=mc[c[(c[ba>>2]|0)+24>>2]&63](ba)|0;c[X>>2]=c[e>>2];c[Y>>2]=c[f>>2];ba=a[aa]|0;if((ba&1)==0){ia=aa+4|0;ja=(ba&255)>>>1;ka=aa+4|0}else{ba=c[aa+8>>2]|0;ia=ba;ja=c[aa+4>>2]|0;ka=ba}c[n+0>>2]=c[X+0>>2];c[m+0>>2]=c[Y+0>>2];Wh(W,d,n,m,g,h,j,ka,ia+(ja<<2)|0);c[e>>2]=c[W>>2];break};default:{c[h>>2]=c[h>>2]|4}}}while(0);c[b>>2]=c[e>>2];i=l;return}function ei(a,b,d,e,f){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0;a=i;a:while(1){g=c[b>>2]|0;do{if((g|0)!=0){h=c[g+12>>2]|0;if((h|0)==(c[g+16>>2]|0)){j=mc[c[(c[g>>2]|0)+36>>2]&63](g)|0}else{j=c[h>>2]|0}if((j|0)==-1){c[b>>2]=0;k=1;break}else{k=(c[b>>2]|0)==0;break}}else{k=1}}while(0);g=c[d>>2]|0;do{if((g|0)!=0){h=c[g+12>>2]|0;if((h|0)==(c[g+16>>2]|0)){l=mc[c[(c[g>>2]|0)+36>>2]&63](g)|0}else{l=c[h>>2]|0}if(!((l|0)==-1)){if(k){m=g;break}else{n=g;break a}}else{c[d>>2]=0;o=15;break}}else{o=15}}while(0);if((o|0)==15){o=0;if(k){n=0;break}else{m=0}}g=c[b>>2]|0;h=c[g+12>>2]|0;if((h|0)==(c[g+16>>2]|0)){p=mc[c[(c[g>>2]|0)+36>>2]&63](g)|0}else{p=c[h>>2]|0}if(!(gc[c[(c[f>>2]|0)+12>>2]&31](f,8192,p)|0)){n=m;break}h=c[b>>2]|0;g=h+12|0;q=c[g>>2]|0;if((q|0)==(c[h+16>>2]|0)){mc[c[(c[h>>2]|0)+40>>2]&63](h)|0;continue}else{c[g>>2]=q+4;continue}}m=c[b>>2]|0;do{if((m|0)!=0){p=c[m+12>>2]|0;if((p|0)==(c[m+16>>2]|0)){r=mc[c[(c[m>>2]|0)+36>>2]&63](m)|0}else{r=c[p>>2]|0}if((r|0)==-1){c[b>>2]=0;s=1;break}else{s=(c[b>>2]|0)==0;break}}else{s=1}}while(0);do{if((n|0)!=0){b=c[n+12>>2]|0;if((b|0)==(c[n+16>>2]|0)){t=mc[c[(c[n>>2]|0)+36>>2]&63](n)|0}else{t=c[b>>2]|0}if((t|0)==-1){c[d>>2]=0;o=37;break}if(s){i=a;return}}else{o=37}}while(0);if((o|0)==37?!s:0){i=a;return}c[e>>2]=c[e>>2]|2;i=a;return}function fi(a,b,d,e,f){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0;a=i;g=c[b>>2]|0;do{if((g|0)!=0){h=c[g+12>>2]|0;if((h|0)==(c[g+16>>2]|0)){j=mc[c[(c[g>>2]|0)+36>>2]&63](g)|0}else{j=c[h>>2]|0}if((j|0)==-1){c[b>>2]=0;k=1;break}else{k=(c[b>>2]|0)==0;break}}else{k=1}}while(0);j=c[d>>2]|0;do{if((j|0)!=0){g=c[j+12>>2]|0;if((g|0)==(c[j+16>>2]|0)){l=mc[c[(c[j>>2]|0)+36>>2]&63](j)|0}else{l=c[g>>2]|0}if(!((l|0)==-1)){if(k){m=j;break}else{n=16;break}}else{c[d>>2]=0;n=14;break}}else{n=14}}while(0);if((n|0)==14){if(k){n=16}else{m=0}}if((n|0)==16){c[e>>2]=c[e>>2]|6;i=a;return}k=c[b>>2]|0;j=c[k+12>>2]|0;if((j|0)==(c[k+16>>2]|0)){o=mc[c[(c[k>>2]|0)+36>>2]&63](k)|0}else{o=c[j>>2]|0}if(!((gc[c[(c[f>>2]|0)+52>>2]&31](f,o,0)|0)<<24>>24==37)){c[e>>2]=c[e>>2]|4;i=a;return}o=c[b>>2]|0;f=o+12|0;j=c[f>>2]|0;if((j|0)==(c[o+16>>2]|0)){mc[c[(c[o>>2]|0)+40>>2]&63](o)|0}else{c[f>>2]=j+4}j=c[b>>2]|0;do{if((j|0)!=0){f=c[j+12>>2]|0;if((f|0)==(c[j+16>>2]|0)){p=mc[c[(c[j>>2]|0)+36>>2]&63](j)|0}else{p=c[f>>2]|0}if((p|0)==-1){c[b>>2]=0;q=1;break}else{q=(c[b>>2]|0)==0;break}}else{q=1}}while(0);do{if((m|0)!=0){b=c[m+12>>2]|0;if((b|0)==(c[m+16>>2]|0)){r=mc[c[(c[m>>2]|0)+36>>2]&63](m)|0}else{r=c[b>>2]|0}if((r|0)==-1){c[d>>2]=0;n=38;break}if(q){i=a;return}}else{n=38}}while(0);if((n|0)==38?!q:0){i=a;return}c[e>>2]=c[e>>2]|2;i=a;return}function gi(a,b,d,e,f){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0;g=i;h=c[a>>2]|0;do{if((h|0)!=0){j=c[h+12>>2]|0;if((j|0)==(c[h+16>>2]|0)){k=mc[c[(c[h>>2]|0)+36>>2]&63](h)|0}else{k=c[j>>2]|0}if((k|0)==-1){c[a>>2]=0;l=1;break}else{l=(c[a>>2]|0)==0;break}}else{l=1}}while(0);k=c[b>>2]|0;do{if((k|0)!=0){h=c[k+12>>2]|0;if((h|0)==(c[k+16>>2]|0)){m=mc[c[(c[k>>2]|0)+36>>2]&63](k)|0}else{m=c[h>>2]|0}if(!((m|0)==-1)){if(l){n=k;break}else{o=16;break}}else{c[b>>2]=0;o=14;break}}else{o=14}}while(0);if((o|0)==14){if(l){o=16}else{n=0}}if((o|0)==16){c[d>>2]=c[d>>2]|6;p=0;i=g;return p|0}l=c[a>>2]|0;k=c[l+12>>2]|0;if((k|0)==(c[l+16>>2]|0)){q=mc[c[(c[l>>2]|0)+36>>2]&63](l)|0}else{q=c[k>>2]|0}if(!(gc[c[(c[e>>2]|0)+12>>2]&31](e,2048,q)|0)){c[d>>2]=c[d>>2]|4;p=0;i=g;return p|0}k=(gc[c[(c[e>>2]|0)+52>>2]&31](e,q,0)|0)<<24>>24;q=c[a>>2]|0;l=q+12|0;m=c[l>>2]|0;if((m|0)==(c[q+16>>2]|0)){mc[c[(c[q>>2]|0)+40>>2]&63](q)|0;r=f;s=n;t=n;u=k}else{c[l>>2]=m+4;r=f;s=n;t=n;u=k}while(1){v=u+ -48|0;k=r+ -1|0;n=c[a>>2]|0;do{if((n|0)!=0){f=c[n+12>>2]|0;if((f|0)==(c[n+16>>2]|0)){w=mc[c[(c[n>>2]|0)+36>>2]&63](n)|0}else{w=c[f>>2]|0}if((w|0)==-1){c[a>>2]=0;x=1;break}else{x=(c[a>>2]|0)==0;break}}else{x=1}}while(0);do{if((t|0)!=0){n=c[t+12>>2]|0;if((n|0)==(c[t+16>>2]|0)){y=mc[c[(c[t>>2]|0)+36>>2]&63](t)|0}else{y=c[n>>2]|0}if((y|0)==-1){c[b>>2]=0;z=0;A=0;B=1;break}else{z=s;A=s;B=(s|0)==0;break}}else{z=s;A=0;B=1}}while(0);C=c[a>>2]|0;if(!((x^B)&(k|0)>0)){break}n=c[C+12>>2]|0;if((n|0)==(c[C+16>>2]|0)){D=mc[c[(c[C>>2]|0)+36>>2]&63](C)|0}else{D=c[n>>2]|0}if(!(gc[c[(c[e>>2]|0)+12>>2]&31](e,2048,D)|0)){p=v;o=63;break}n=((gc[c[(c[e>>2]|0)+52>>2]&31](e,D,0)|0)<<24>>24)+(v*10|0)|0;f=c[a>>2]|0;m=f+12|0;l=c[m>>2]|0;if((l|0)==(c[f+16>>2]|0)){mc[c[(c[f>>2]|0)+40>>2]&63](f)|0;r=k;s=z;t=A;u=n;continue}else{c[m>>2]=l+4;r=k;s=z;t=A;u=n;continue}}if((o|0)==63){i=g;return p|0}do{if((C|0)!=0){u=c[C+12>>2]|0;if((u|0)==(c[C+16>>2]|0)){E=mc[c[(c[C>>2]|0)+36>>2]&63](C)|0}else{E=c[u>>2]|0}if((E|0)==-1){c[a>>2]=0;F=1;break}else{F=(c[a>>2]|0)==0;break}}else{F=1}}while(0);do{if((z|0)!=0){a=c[z+12>>2]|0;if((a|0)==(c[z+16>>2]|0)){G=mc[c[(c[z>>2]|0)+36>>2]&63](z)|0}else{G=c[a>>2]|0}if((G|0)==-1){c[b>>2]=0;o=60;break}if(F){p=v;i=g;return p|0}}else{o=60}}while(0);if((o|0)==60?!F:0){p=v;i=g;return p|0}c[d>>2]=c[d>>2]|2;p=v;i=g;return p|0}function hi(b){b=b|0;var d=0,e=0,f=0;d=i;e=b+8|0;f=c[e>>2]|0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}if((f|0)==(c[2036]|0)){Um(b);i=d;return}Ea(c[e>>2]|0);Um(b);i=d;return}function ii(b){b=b|0;var d=0,e=0;d=i;e=b+8|0;b=c[e>>2]|0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}if((b|0)==(c[2036]|0)){i=d;return}Ea(c[e>>2]|0);i=d;return}function ji(b,d,e,f,g,h,j,k){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0;g=i;i=i+112|0;f=g+100|0;l=g;a[f]=37;m=f+1|0;a[m]=j;n=f+2|0;a[n]=k;a[f+3|0]=0;if(!(k<<24>>24==0)){a[m]=k;a[n]=j}j=fb(l|0,100,f|0,h|0,c[d+8>>2]|0)|0;d=l+j|0;h=c[e>>2]|0;if((j|0)==0){o=h;c[b>>2]=o;i=g;return}else{p=l;q=h;r=h}while(1){h=a[p]|0;do{if((q|0)!=0){l=q+24|0;j=c[l>>2]|0;if((j|0)==(c[q+28>>2]|0)){e=(vc[c[(c[q>>2]|0)+52>>2]&15](q,h&255)|0)==-1;s=e?0:r;t=e?0:q;break}else{c[l>>2]=j+1;a[j]=h;s=r;t=q;break}}else{s=r;t=0}}while(0);h=p+1|0;if((h|0)==(d|0)){o=s;break}else{p=h;q=t;r=s}}c[b>>2]=o;i=g;return}function ki(b){b=b|0;var d=0,e=0,f=0;d=i;e=b+8|0;f=c[e>>2]|0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}if((f|0)==(c[2036]|0)){Um(b);i=d;return}Ea(c[e>>2]|0);Um(b);i=d;return}function li(b){b=b|0;var d=0,e=0;d=i;e=b+8|0;b=c[e>>2]|0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}if((b|0)==(c[2036]|0)){i=d;return}Ea(c[e>>2]|0);i=d;return}function mi(a,b,d,e,f,g,h,j){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0;f=i;i=i+416|0;e=f+8|0;k=f;c[k>>2]=e+400;ni(b+8|0,e,k,g,h,j);j=c[k>>2]|0;k=c[d>>2]|0;if((e|0)==(j|0)){l=k;c[a>>2]=l;i=f;return}else{m=e;n=k;o=k}while(1){k=c[m>>2]|0;if((o|0)==0){p=n;q=0}else{e=o+24|0;d=c[e>>2]|0;if((d|0)==(c[o+28>>2]|0)){r=vc[c[(c[o>>2]|0)+52>>2]&15](o,k)|0}else{c[e>>2]=d+4;c[d>>2]=k;r=k}k=(r|0)==-1;p=k?0:n;q=k?0:o}k=m+4|0;if((k|0)==(j|0)){l=p;break}else{m=k;n=p;o=q}}c[a>>2]=l;i=f;return}function ni(b,d,e,f,g,h){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0;j=i;i=i+128|0;k=j+112|0;l=j+12|0;m=j;n=j+8|0;a[k]=37;o=k+1|0;a[o]=g;p=k+2|0;a[p]=h;a[k+3|0]=0;if(!(h<<24>>24==0)){a[o]=h;a[p]=g}fb(l|0,100,k|0,f|0,c[b>>2]|0)|0;f=m;c[f>>2]=0;c[f+4>>2]=0;c[n>>2]=l;l=(c[e>>2]|0)-d>>2;f=jb(c[b>>2]|0)|0;b=gm(d,n,l,m)|0;if((f|0)!=0){jb(f|0)|0}if((b|0)==-1){jj(9136)}else{c[e>>2]=d+(b<<2);i=j;return}}function oi(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function pi(a){a=a|0;return}function qi(a){a=a|0;return 127}function ri(a){a=a|0;return 127}function si(a,b){a=a|0;b=b|0;b=i;c[a+0>>2]=0;c[a+4>>2]=0;c[a+8>>2]=0;i=b;return}function ti(a,b){a=a|0;b=b|0;b=i;c[a+0>>2]=0;c[a+4>>2]=0;c[a+8>>2]=0;i=b;return}function ui(a,b){a=a|0;b=b|0;b=i;c[a+0>>2]=0;c[a+4>>2]=0;c[a+8>>2]=0;i=b;return}function vi(a,b){a=a|0;b=b|0;b=i;ye(a,1,45);i=b;return}function wi(a){a=a|0;return 0}function xi(b,c){b=b|0;c=c|0;a[b]=67109634;a[b+1|0]=262147;a[b+2|0]=1024;a[b+3|0]=4;return}function yi(b,c){b=b|0;c=c|0;a[b]=67109634;a[b+1|0]=262147;a[b+2|0]=1024;a[b+3|0]=4;return}function zi(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function Ai(a){a=a|0;return}function Bi(a){a=a|0;return 127}function Ci(a){a=a|0;return 127}function Di(a,b){a=a|0;b=b|0;b=i;c[a+0>>2]=0;c[a+4>>2]=0;c[a+8>>2]=0;i=b;return}function Ei(a,b){a=a|0;b=b|0;b=i;c[a+0>>2]=0;c[a+4>>2]=0;c[a+8>>2]=0;i=b;return}function Fi(a,b){a=a|0;b=b|0;b=i;c[a+0>>2]=0;c[a+4>>2]=0;c[a+8>>2]=0;i=b;return}function Gi(a,b){a=a|0;b=b|0;b=i;ye(a,1,45);i=b;return}function Hi(a){a=a|0;return 0}function Ii(b,c){b=b|0;c=c|0;a[b]=67109634;a[b+1|0]=262147;a[b+2|0]=1024;a[b+3|0]=4;return}function Ji(b,c){b=b|0;c=c|0;a[b]=67109634;a[b+1|0]=262147;a[b+2|0]=1024;a[b+3|0]=4;return}function Ki(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function Li(a){a=a|0;return}function Mi(a){a=a|0;return 2147483647}function Ni(a){a=a|0;return 2147483647}function Oi(a,b){a=a|0;b=b|0;b=i;c[a+0>>2]=0;c[a+4>>2]=0;c[a+8>>2]=0;i=b;return}function Pi(a,b){a=a|0;b=b|0;b=i;c[a+0>>2]=0;c[a+4>>2]=0;c[a+8>>2]=0;i=b;return}function Qi(a,b){a=a|0;b=b|0;b=i;c[a+0>>2]=0;c[a+4>>2]=0;c[a+8>>2]=0;i=b;return}function Ri(a,b){a=a|0;b=b|0;b=i;Je(a,1,45);i=b;return}function Si(a){a=a|0;return 0}function Ti(b,c){b=b|0;c=c|0;a[b]=67109634;a[b+1|0]=262147;a[b+2|0]=1024;a[b+3|0]=4;return}function Ui(b,c){b=b|0;c=c|0;a[b]=67109634;a[b+1|0]=262147;a[b+2|0]=1024;a[b+3|0]=4;return}function Vi(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function Wi(a){a=a|0;return}function Xi(a){a=a|0;return 2147483647}function Yi(a){a=a|0;return 2147483647}function Zi(a,b){a=a|0;b=b|0;b=i;c[a+0>>2]=0;c[a+4>>2]=0;c[a+8>>2]=0;i=b;return}function _i(a,b){a=a|0;b=b|0;b=i;c[a+0>>2]=0;c[a+4>>2]=0;c[a+8>>2]=0;i=b;return}function $i(a,b){a=a|0;b=b|0;b=i;c[a+0>>2]=0;c[a+4>>2]=0;c[a+8>>2]=0;i=b;return}function aj(a,b){a=a|0;b=b|0;b=i;Je(a,1,45);i=b;return}function bj(a){a=a|0;return 0}function cj(b,c){b=b|0;c=c|0;a[b]=67109634;a[b+1|0]=262147;a[b+2|0]=1024;a[b+3|0]=4;return}function dj(b,c){b=b|0;c=c|0;a[b]=67109634;a[b+1|0]=262147;a[b+2|0]=1024;a[b+3|0]=4;return}function ej(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function fj(a){a=a|0;return}function gj(b,d,e,f,g,h,j,k){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0;d=i;i=i+256|0;l=d;m=d+32|0;n=d+24|0;o=d+16|0;p=d+12|0;q=d+244|0;r=d+20|0;s=d+132|0;t=d+144|0;c[n>>2]=m;u=n+4|0;c[u>>2]=114;v=m+100|0;Te(p,h);m=c[p>>2]|0;if(!((c[2064]|0)==-1)){c[l>>2]=8256;c[l+4>>2]=112;c[l+8>>2]=0;ue(8256,l,113)}w=(c[8260>>2]|0)+ -1|0;x=c[m+8>>2]|0;if(!((c[m+12>>2]|0)-x>>2>>>0>w>>>0)){y=Wa(4)|0;qm(y);xb(y|0,16216,101)}m=c[x+(w<<2)>>2]|0;if((m|0)==0){y=Wa(4)|0;qm(y);xb(y|0,16216,101)}a[q]=0;c[r>>2]=c[f>>2];y=c[h+4>>2]|0;c[l+0>>2]=c[r+0>>2];if(ij(e,l,g,p,y,j,q,m,n,o,v)|0){sc[c[(c[m>>2]|0)+32>>2]&7](m,7872,7882|0,s)|0;m=c[o>>2]|0;v=c[n>>2]|0;y=m-v|0;if((y|0)>98){g=Lm(y+2|0)|0;if((g|0)==0){_m()}else{z=g;A=g}}else{z=0;A=t}if((a[q]|0)==0){B=A}else{a[A]=45;B=A+1|0}if(v>>>0<m>>>0){m=s+10|0;A=s;q=B;g=v;while(1){v=a[g]|0;y=s;while(1){r=y+1|0;if((a[y]|0)==v<<24>>24){C=y;break}if((r|0)==(m|0)){C=m;break}else{y=r}}a[q]=a[7872+(C-A)|0]|0;y=g+1|0;v=q+1|0;if(y>>>0<(c[o>>2]|0)>>>0){q=v;g=y}else{D=v;break}}}else{D=B}a[D]=0;c[l>>2]=k;if((gb(t|0,7888,l|0)|0)!=1){l=Wa(8)|0;ge(l,7896);xb(l|0,5264,11)}if((z|0)!=0){Mm(z)}}z=c[e>>2]|0;if((z|0)!=0){if((c[z+12>>2]|0)==(c[z+16>>2]|0)?(mc[c[(c[z>>2]|0)+36>>2]&63](z)|0)==-1:0){c[e>>2]=0;E=0}else{E=z}}else{E=0}z=(E|0)==0;e=c[f>>2]|0;do{if((e|0)!=0){if((c[e+12>>2]|0)!=(c[e+16>>2]|0)){if(z){break}else{F=33;break}}if(!((mc[c[(c[e>>2]|0)+36>>2]&63](e)|0)==-1)){if(z){break}else{F=33;break}}else{c[f>>2]=0;F=31;break}}else{F=31}}while(0);if((F|0)==31?z:0){F=33}if((F|0)==33){c[j>>2]=c[j>>2]|2}c[b>>2]=E;ae(c[p>>2]|0)|0;p=c[n>>2]|0;c[n>>2]=0;if((p|0)==0){i=d;return}jc[c[u>>2]&127](p);i=d;return}function hj(a){a=a|0;return}function ij(e,f,g,h,j,k,l,m,n,o,p){e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;m=m|0;n=n|0;o=o|0;p=p|0;var q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0,X=0,Y=0,Z=0,_=0,$=0,aa=0,ba=0,ca=0,da=0,ea=0,fa=0,ga=0,ha=0,ia=0,ja=0,ka=0,la=0,ma=0,na=0,oa=0,pa=0,qa=0,ra=0,sa=0,ta=0,ua=0,va=0,wa=0,xa=0,ya=0,za=0,Aa=0,Ba=0,Ca=0,Da=0,Ea=0,Fa=0,Ga=0,Ha=0,Ia=0,Ja=0,Ka=0,La=0,Ma=0,Na=0,Oa=0,Pa=0,Qa=0,Ra=0,Sa=0,Ta=0,Ua=0,Va=0,Wa=0,Xa=0,Ya=0,Za=0,_a=0,$a=0,ab=0,bb=0,cb=0,db=0,eb=0,fb=0,gb=0,hb=0,ib=0,jb=0,kb=0,lb=0,mb=0,nb=0,ob=0,pb=0,qb=0,rb=0,sb=0,tb=0,ub=0,vb=0,wb=0,xb=0,yb=0,zb=0,Ab=0,Bb=0,Cb=0,Db=0,Eb=0,Fb=0,Gb=0,Hb=0,Ib=0,Jb=0,Kb=0,Lb=0,Mb=0,Nb=0,Ob=0,Pb=0,Qb=0,Rb=0,Sb=0,Tb=0,Ub=0,Vb=0,Wb=0,Xb=0,Yb=0,Zb=0,_b=0,$b=0,ac=0,bc=0,cc=0,dc=0,ec=0,fc=0,gc=0,hc=0,ic=0,kc=0,lc=0;q=i;i=i+480|0;r=q+72|0;s=q+68|0;t=q+472|0;u=q+473|0;v=q+56|0;w=q+44|0;x=q+28|0;y=q+16|0;z=q+4|0;A=q;B=q+40|0;c[s>>2]=0;c[v+0>>2]=0;c[v+4>>2]=0;c[v+8>>2]=0;c[w+0>>2]=0;c[w+4>>2]=0;c[w+8>>2]=0;c[x+0>>2]=0;c[x+4>>2]=0;c[x+8>>2]=0;c[y+0>>2]=0;c[y+4>>2]=0;c[y+8>>2]=0;c[z+0>>2]=0;c[z+4>>2]=0;c[z+8>>2]=0;mj(g,h,s,t,u,v,w,x,y,A);c[o>>2]=c[n>>2];h=m+8|0;m=y+1|0;g=y+4|0;C=y+8|0;D=x+1|0;E=x+4|0;F=x+8|0;G=(j&512|0)!=0;j=w+1|0;H=w+8|0;I=w+4|0;J=z+1|0;K=z+8|0;L=z+4|0;M=s+3|0;N=n+4|0;O=v+4|0;P=r+400|0;Q=r;R=r;r=p;p=0;S=0;T=114;a:while(1){U=c[e>>2]|0;do{if((U|0)!=0){if((c[U+12>>2]|0)==(c[U+16>>2]|0)){if((mc[c[(c[U>>2]|0)+36>>2]&63](U)|0)==-1){c[e>>2]=0;V=0;break}else{V=c[e>>2]|0;break}}else{V=U}}else{V=0}}while(0);U=(V|0)==0;W=c[f>>2]|0;do{if((W|0)!=0){if((c[W+12>>2]|0)!=(c[W+16>>2]|0)){if(U){X=W;break}else{Y=R;Z=Q;_=S;$=T;aa=269;break a}}if(!((mc[c[(c[W>>2]|0)+36>>2]&63](W)|0)==-1)){if(U){X=W;break}else{Y=R;Z=Q;_=S;$=T;aa=269;break a}}else{c[f>>2]=0;aa=12;break}}else{aa=12}}while(0);if((aa|0)==12){aa=0;if(U){Y=R;Z=Q;_=S;$=T;aa=269;break}else{X=0}}b:do{switch(a[s+p|0]|0){case 0:{aa=26;break};case 1:{if((p|0)==3){Y=R;Z=Q;_=S;$=T;aa=269;break a}W=c[e>>2]|0;ba=c[W+12>>2]|0;if((ba|0)==(c[W+16>>2]|0)){ca=mc[c[(c[W>>2]|0)+36>>2]&63](W)|0}else{ca=d[ba]|0}if(!((ca&255)<<24>>24>-1)){aa=25;break a}if((b[(c[h>>2]|0)+(ca<<24>>24<<1)>>1]&8192)==0){aa=25;break a}ba=c[e>>2]|0;W=ba+12|0;da=c[W>>2]|0;if((da|0)==(c[ba+16>>2]|0)){ea=mc[c[(c[ba>>2]|0)+40>>2]&63](ba)|0}else{c[W>>2]=da+1;ea=d[da]|0}Ee(z,ea&255);aa=26;break};case 3:{da=a[x]|0;W=(da&1)==0;if(W){fa=(da&255)>>>1}else{fa=c[E>>2]|0}ba=a[y]|0;ga=(ba&1)==0;if(ga){ha=(ba&255)>>>1}else{ha=c[g>>2]|0}if((fa|0)==(0-ha|0)){ia=r;ja=R;ka=Q;la=P;ma=S;na=T}else{if(W){oa=(da&255)>>>1}else{oa=c[E>>2]|0}if((oa|0)!=0){if(ga){pa=(ba&255)>>>1}else{pa=c[g>>2]|0}if((pa|0)!=0){ga=c[e>>2]|0;qa=c[ga+12>>2]|0;ra=c[ga+16>>2]|0;if((qa|0)==(ra|0)){sa=mc[c[(c[ga>>2]|0)+36>>2]&63](ga)|0;ta=c[e>>2]|0;ua=sa;va=a[x]|0;wa=ta;xa=c[ta+12>>2]|0;ya=c[ta+16>>2]|0}else{ua=d[qa]|0;va=da;wa=ga;xa=qa;ya=ra}ra=wa+12|0;qa=(xa|0)==(ya|0);if((ua&255)<<24>>24==(a[(va&1)==0?D:c[F>>2]|0]|0)){if(qa){mc[c[(c[wa>>2]|0)+40>>2]&63](wa)|0}else{c[ra>>2]=xa+1}ra=a[x]|0;if((ra&1)==0){za=(ra&255)>>>1}else{za=c[E>>2]|0}ia=r;ja=R;ka=Q;la=P;ma=za>>>0>1?x:S;na=T;break b}if(qa){Aa=mc[c[(c[wa>>2]|0)+36>>2]&63](wa)|0}else{Aa=d[xa]|0}if(!((Aa&255)<<24>>24==(a[(a[y]&1)==0?m:c[C>>2]|0]|0))){aa=112;break a}qa=c[e>>2]|0;ra=qa+12|0;ga=c[ra>>2]|0;if((ga|0)==(c[qa+16>>2]|0)){mc[c[(c[qa>>2]|0)+40>>2]&63](qa)|0}else{c[ra>>2]=ga+1}a[l]=1;ga=a[y]|0;if((ga&1)==0){Ba=(ga&255)>>>1}else{Ba=c[g>>2]|0}ia=r;ja=R;ka=Q;la=P;ma=Ba>>>0>1?y:S;na=T;break b}}if(W){Ca=(da&255)>>>1}else{Ca=c[E>>2]|0}W=c[e>>2]|0;ga=c[W+12>>2]|0;ra=(ga|0)==(c[W+16>>2]|0);if((Ca|0)==0){if(ra){qa=mc[c[(c[W>>2]|0)+36>>2]&63](W)|0;Da=qa;Ea=a[y]|0}else{Da=d[ga]|0;Ea=ba}if(!((Da&255)<<24>>24==(a[(Ea&1)==0?m:c[C>>2]|0]|0))){ia=r;ja=R;ka=Q;la=P;ma=S;na=T;break b}ba=c[e>>2]|0;qa=ba+12|0;ta=c[qa>>2]|0;if((ta|0)==(c[ba+16>>2]|0)){mc[c[(c[ba>>2]|0)+40>>2]&63](ba)|0}else{c[qa>>2]=ta+1}a[l]=1;ta=a[y]|0;if((ta&1)==0){Fa=(ta&255)>>>1}else{Fa=c[g>>2]|0}ia=r;ja=R;ka=Q;la=P;ma=Fa>>>0>1?y:S;na=T;break b}if(ra){ra=mc[c[(c[W>>2]|0)+36>>2]&63](W)|0;Ga=ra;Ha=a[x]|0}else{Ga=d[ga]|0;Ha=da}if(!((Ga&255)<<24>>24==(a[(Ha&1)==0?D:c[F>>2]|0]|0))){a[l]=1;ia=r;ja=R;ka=Q;la=P;ma=S;na=T;break b}da=c[e>>2]|0;ga=da+12|0;ra=c[ga>>2]|0;if((ra|0)==(c[da+16>>2]|0)){mc[c[(c[da>>2]|0)+40>>2]&63](da)|0}else{c[ga>>2]=ra+1}ra=a[x]|0;if((ra&1)==0){Ia=(ra&255)>>>1}else{Ia=c[E>>2]|0}ia=r;ja=R;ka=Q;la=P;ma=Ia>>>0>1?x:S;na=T}break};case 2:{if(!((S|0)!=0|p>>>0<2)){if((p|0)==2){Ja=(a[M]|0)!=0}else{Ja=0}if(!(G|Ja)){ia=r;ja=R;ka=Q;la=P;ma=0;na=T;break b}}ra=a[w]|0;ga=(ra&1)==0;da=ga?j:c[H>>2]|0;c:do{if((p|0)!=0?(d[s+(p+ -1)|0]|0)<2:0){W=da+(ga?(ra&255)>>>1:c[I>>2]|0)|0;ta=da;while(1){if((ta|0)==(W|0)){Ka=W;break}qa=a[ta]|0;if(!(qa<<24>>24>-1)){Ka=ta;break}if((b[(c[h>>2]|0)+(qa<<24>>24<<1)>>1]&8192)==0){Ka=ta;break}else{ta=ta+1|0}}ta=Ka-da|0;W=a[z]|0;qa=(W&1)==0;if(qa){La=(W&255)>>>1}else{La=c[L>>2]|0}if(!(ta>>>0>La>>>0)){if(qa){qa=(W&255)>>>1;Ma=J;Na=qa;Oa=z+(qa-ta)+1|0}else{qa=c[K>>2]|0;W=c[L>>2]|0;Ma=qa;Na=W;Oa=qa+(W-ta)|0}ta=Ma+Na|0;if((Oa|0)==(ta|0)){Pa=X;Qa=ra;Ra=Ka;Sa=X}else{W=Oa;qa=da;while(1){if((a[W]|0)!=(a[qa]|0)){Pa=X;Qa=ra;Ra=da;Sa=X;break c}ba=W+1|0;if((ba|0)==(ta|0)){Pa=X;Qa=ra;Ra=Ka;Sa=X;break}else{W=ba;qa=qa+1|0}}}}else{Pa=X;Qa=ra;Ra=da;Sa=X}}else{Pa=X;Qa=ra;Ra=da;Sa=X}}while(0);d:while(1){if((Qa&1)==0){Ta=j;Ua=(Qa&255)>>>1}else{Ta=c[H>>2]|0;Ua=c[I>>2]|0}if((Ra|0)==(Ta+Ua|0)){break}da=c[e>>2]|0;do{if((da|0)!=0){if((c[da+12>>2]|0)==(c[da+16>>2]|0)){if((mc[c[(c[da>>2]|0)+36>>2]&63](da)|0)==-1){c[e>>2]=0;Va=0;break}else{Va=c[e>>2]|0;break}}else{Va=da}}else{Va=0}}while(0);da=(Va|0)==0;do{if((Sa|0)!=0){if((c[Sa+12>>2]|0)!=(c[Sa+16>>2]|0)){if(da){Wa=Pa;Xa=Sa;break}else{break d}}if(!((mc[c[(c[Sa>>2]|0)+36>>2]&63](Sa)|0)==-1)){if(da^(Pa|0)==0){Wa=Pa;Xa=Pa;break}else{break d}}else{c[f>>2]=0;Ya=0;aa=147;break}}else{Ya=Pa;aa=147}}while(0);if((aa|0)==147){aa=0;if(da){break}else{Wa=Ya;Xa=0}}ra=c[e>>2]|0;ga=c[ra+12>>2]|0;if((ga|0)==(c[ra+16>>2]|0)){Za=mc[c[(c[ra>>2]|0)+36>>2]&63](ra)|0}else{Za=d[ga]|0}if(!((Za&255)<<24>>24==(a[Ra]|0))){break}ga=c[e>>2]|0;ra=ga+12|0;qa=c[ra>>2]|0;if((qa|0)==(c[ga+16>>2]|0)){mc[c[(c[ga>>2]|0)+40>>2]&63](ga)|0}else{c[ra>>2]=qa+1}Pa=Wa;Qa=a[w]|0;Ra=Ra+1|0;Sa=Xa}if(G){qa=a[w]|0;if((qa&1)==0){_a=j;$a=(qa&255)>>>1}else{_a=c[H>>2]|0;$a=c[I>>2]|0}if((Ra|0)!=(_a+$a|0)){aa=162;break a}else{ia=r;ja=R;ka=Q;la=P;ma=S;na=T}}else{ia=r;ja=R;ka=Q;la=P;ma=S;na=T}break};case 4:{qa=r;ra=Q;ga=P;W=R;ta=0;ba=T;e:while(1){sa=c[e>>2]|0;do{if((sa|0)!=0){if((c[sa+12>>2]|0)==(c[sa+16>>2]|0)){if((mc[c[(c[sa>>2]|0)+36>>2]&63](sa)|0)==-1){c[e>>2]=0;ab=0;break}else{ab=c[e>>2]|0;break}}else{ab=sa}}else{ab=0}}while(0);sa=(ab|0)==0;da=c[f>>2]|0;do{if((da|0)!=0){if((c[da+12>>2]|0)!=(c[da+16>>2]|0)){if(sa){break}else{break e}}if(!((mc[c[(c[da>>2]|0)+36>>2]&63](da)|0)==-1)){if(sa){break}else{break e}}else{c[f>>2]=0;aa=173;break}}else{aa=173}}while(0);if((aa|0)==173?(aa=0,sa):0){break}da=c[e>>2]|0;bb=c[da+12>>2]|0;if((bb|0)==(c[da+16>>2]|0)){cb=mc[c[(c[da>>2]|0)+36>>2]&63](da)|0}else{cb=d[bb]|0}bb=cb&255;if(bb<<24>>24>-1?!((b[(c[h>>2]|0)+(cb<<24>>24<<1)>>1]&2048)==0):0){da=c[o>>2]|0;if((da|0)==(qa|0)){db=(c[N>>2]|0)!=114;eb=c[n>>2]|0;fb=qa-eb|0;gb=fb>>>0<2147483647?fb<<1:-1;hb=Nm(db?eb:0,gb)|0;if((hb|0)==0){aa=182;break a}if(!db){db=c[n>>2]|0;c[n>>2]=hb;if((db|0)==0){ib=hb}else{jc[c[N>>2]&127](db);ib=c[n>>2]|0}}else{c[n>>2]=hb;ib=hb}c[N>>2]=115;hb=ib+fb|0;c[o>>2]=hb;jb=hb;kb=(c[n>>2]|0)+gb|0}else{jb=da;kb=qa}c[o>>2]=jb+1;a[jb]=bb;lb=kb;mb=W;nb=ra;ob=ga;pb=ta+1|0;qb=ba}else{da=a[v]|0;if((da&1)==0){rb=(da&255)>>>1}else{rb=c[O>>2]|0}if((rb|0)==0|(ta|0)==0){break}if(!(bb<<24>>24==(a[u]|0))){break}if((ra|0)==(ga|0)){bb=ra-W|0;da=bb>>>0<2147483647?bb<<1:-1;if((ba|0)==114){sb=0}else{sb=W}gb=Nm(sb,da)|0;if((gb|0)==0){aa=198;break a}tb=gb+(bb>>2<<2)|0;ub=gb;vb=gb+(da>>>2<<2)|0;wb=115}else{tb=ra;ub=W;vb=ga;wb=ba}c[tb>>2]=ta;lb=qa;mb=ub;nb=tb+4|0;ob=vb;pb=0;qb=wb}da=c[e>>2]|0;gb=da+12|0;bb=c[gb>>2]|0;if((bb|0)==(c[da+16>>2]|0)){mc[c[(c[da>>2]|0)+40>>2]&63](da)|0;qa=lb;ra=nb;ga=ob;W=mb;ta=pb;ba=qb;continue}else{c[gb>>2]=bb+1;qa=lb;ra=nb;ga=ob;W=mb;ta=pb;ba=qb;continue}}if((W|0)==(ra|0)|(ta|0)==0){xb=W;yb=ra;zb=ga;Ab=ba}else{if((ra|0)==(ga|0)){bb=ra-W|0;gb=bb>>>0<2147483647?bb<<1:-1;if((ba|0)==114){Bb=0}else{Bb=W}da=Nm(Bb,gb)|0;if((da|0)==0){aa=209;break a}Cb=da+(bb>>2<<2)|0;Db=da;Eb=da+(gb>>>2<<2)|0;Fb=115}else{Cb=ra;Db=W;Eb=ga;Fb=ba}c[Cb>>2]=ta;xb=Db;yb=Cb+4|0;zb=Eb;Ab=Fb}gb=c[A>>2]|0;if((gb|0)>0){da=c[e>>2]|0;do{if((da|0)!=0){if((c[da+12>>2]|0)==(c[da+16>>2]|0)){if((mc[c[(c[da>>2]|0)+36>>2]&63](da)|0)==-1){c[e>>2]=0;Gb=0;break}else{Gb=c[e>>2]|0;break}}else{Gb=da}}else{Gb=0}}while(0);da=(Gb|0)==0;ta=c[f>>2]|0;do{if((ta|0)!=0){if((c[ta+12>>2]|0)!=(c[ta+16>>2]|0)){if(da){Hb=ta;break}else{aa=229;break a}}if(!((mc[c[(c[ta>>2]|0)+36>>2]&63](ta)|0)==-1)){if(da){Hb=ta;break}else{aa=229;break a}}else{c[f>>2]=0;aa=223;break}}else{aa=223}}while(0);if((aa|0)==223){aa=0;if(da){aa=229;break a}else{Hb=0}}ta=c[e>>2]|0;ba=c[ta+12>>2]|0;if((ba|0)==(c[ta+16>>2]|0)){Ib=mc[c[(c[ta>>2]|0)+36>>2]&63](ta)|0}else{Ib=d[ba]|0}if(!((Ib&255)<<24>>24==(a[t]|0))){aa=229;break a}ba=c[e>>2]|0;ta=ba+12|0;ga=c[ta>>2]|0;if((ga|0)==(c[ba+16>>2]|0)){mc[c[(c[ba>>2]|0)+40>>2]&63](ba)|0;Jb=Hb;Kb=Hb;Lb=qa;Mb=gb}else{c[ta>>2]=ga+1;Jb=Hb;Kb=Hb;Lb=qa;Mb=gb}while(1){ga=c[e>>2]|0;do{if((ga|0)!=0){if((c[ga+12>>2]|0)==(c[ga+16>>2]|0)){if((mc[c[(c[ga>>2]|0)+36>>2]&63](ga)|0)==-1){c[e>>2]=0;Nb=0;break}else{Nb=c[e>>2]|0;break}}else{Nb=ga}}else{Nb=0}}while(0);ga=(Nb|0)==0;do{if((Kb|0)!=0){if((c[Kb+12>>2]|0)!=(c[Kb+16>>2]|0)){if(ga){Ob=Jb;Pb=Kb;break}else{aa=250;break a}}if(!((mc[c[(c[Kb>>2]|0)+36>>2]&63](Kb)|0)==-1)){if(ga^(Jb|0)==0){Ob=Jb;Pb=Jb;break}else{aa=250;break a}}else{c[f>>2]=0;Qb=0;aa=243;break}}else{Qb=Jb;aa=243}}while(0);if((aa|0)==243){aa=0;if(ga){aa=250;break a}else{Ob=Qb;Pb=0}}sa=c[e>>2]|0;ta=c[sa+12>>2]|0;if((ta|0)==(c[sa+16>>2]|0)){Rb=mc[c[(c[sa>>2]|0)+36>>2]&63](sa)|0}else{Rb=d[ta]|0}if(!((Rb&255)<<24>>24>-1)){aa=250;break a}if((b[(c[h>>2]|0)+(Rb<<24>>24<<1)>>1]&2048)==0){aa=250;break a}ta=c[o>>2]|0;if((ta|0)==(Lb|0)){sa=(c[N>>2]|0)!=114;ba=c[n>>2]|0;W=Lb-ba|0;ra=W>>>0<2147483647?W<<1:-1;bb=Nm(sa?ba:0,ra)|0;if((bb|0)==0){aa=253;break a}if(!sa){sa=c[n>>2]|0;c[n>>2]=bb;if((sa|0)==0){Sb=bb}else{jc[c[N>>2]&127](sa);Sb=c[n>>2]|0}}else{c[n>>2]=bb;Sb=bb}c[N>>2]=115;bb=Sb+W|0;c[o>>2]=bb;Tb=bb;Ub=(c[n>>2]|0)+ra|0}else{Tb=ta;Ub=Lb}ta=c[e>>2]|0;ra=c[ta+12>>2]|0;if((ra|0)==(c[ta+16>>2]|0)){bb=mc[c[(c[ta>>2]|0)+36>>2]&63](ta)|0;Vb=bb;Wb=c[o>>2]|0}else{Vb=d[ra]|0;Wb=Tb}c[o>>2]=Wb+1;a[Wb]=Vb;ra=Mb+ -1|0;c[A>>2]=ra;bb=c[e>>2]|0;ta=bb+12|0;W=c[ta>>2]|0;if((W|0)==(c[bb+16>>2]|0)){mc[c[(c[bb>>2]|0)+40>>2]&63](bb)|0}else{c[ta>>2]=W+1}if((ra|0)>0){Jb=Ob;Kb=Pb;Lb=Ub;Mb=ra}else{Xb=Ub;break}}}else{Xb=qa}if((c[o>>2]|0)==(c[n>>2]|0)){aa=267;break a}else{ia=Xb;ja=xb;ka=yb;la=zb;ma=S;na=Ab}break};default:{ia=r;ja=R;ka=Q;la=P;ma=S;na=T}}}while(0);f:do{if((aa|0)==26){aa=0;if((p|0)==3){Y=R;Z=Q;_=S;$=T;aa=269;break a}else{Yb=X;Zb=X}while(1){U=c[e>>2]|0;do{if((U|0)!=0){if((c[U+12>>2]|0)==(c[U+16>>2]|0)){if((mc[c[(c[U>>2]|0)+36>>2]&63](U)|0)==-1){c[e>>2]=0;_b=0;break}else{_b=c[e>>2]|0;break}}else{_b=U}}else{_b=0}}while(0);U=(_b|0)==0;do{if((Zb|0)!=0){if((c[Zb+12>>2]|0)!=(c[Zb+16>>2]|0)){if(U){$b=Yb;ac=Zb;break}else{ia=r;ja=R;ka=Q;la=P;ma=S;na=T;break f}}if(!((mc[c[(c[Zb>>2]|0)+36>>2]&63](Zb)|0)==-1)){if(U^(Yb|0)==0){$b=Yb;ac=Yb;break}else{ia=r;ja=R;ka=Q;la=P;ma=S;na=T;break f}}else{c[f>>2]=0;bc=0;aa=37;break}}else{bc=Yb;aa=37}}while(0);if((aa|0)==37){aa=0;if(U){ia=r;ja=R;ka=Q;la=P;ma=S;na=T;break f}else{$b=bc;ac=0}}ga=c[e>>2]|0;gb=c[ga+12>>2]|0;if((gb|0)==(c[ga+16>>2]|0)){cc=mc[c[(c[ga>>2]|0)+36>>2]&63](ga)|0}else{cc=d[gb]|0}if(!((cc&255)<<24>>24>-1)){ia=r;ja=R;ka=Q;la=P;ma=S;na=T;break f}if((b[(c[h>>2]|0)+(cc<<24>>24<<1)>>1]&8192)==0){ia=r;ja=R;ka=Q;la=P;ma=S;na=T;break f}gb=c[e>>2]|0;ga=gb+12|0;da=c[ga>>2]|0;if((da|0)==(c[gb+16>>2]|0)){dc=mc[c[(c[gb>>2]|0)+40>>2]&63](gb)|0}else{c[ga>>2]=da+1;dc=d[da]|0}Ee(z,dc&255);Yb=$b;Zb=ac}}}while(0);qa=p+1|0;if(qa>>>0<4){P=la;Q=ka;R=ja;r=ia;p=qa;S=ma;T=na}else{Y=ja;Z=ka;_=ma;$=na;aa=269;break}}g:do{if((aa|0)==25){c[k>>2]=c[k>>2]|4;ec=0;fc=R;gc=T}else if((aa|0)==112){c[k>>2]=c[k>>2]|4;ec=0;fc=R;gc=T}else if((aa|0)==162){c[k>>2]=c[k>>2]|4;ec=0;fc=R;gc=T}else if((aa|0)==182){_m()}else if((aa|0)==198){_m()}else if((aa|0)==209){_m()}else if((aa|0)==229){c[k>>2]=c[k>>2]|4;ec=0;fc=xb;gc=Ab}else if((aa|0)==250){c[k>>2]=c[k>>2]|4;ec=0;fc=xb;gc=Ab}else if((aa|0)==253){_m()}else if((aa|0)==267){c[k>>2]=c[k>>2]|4;ec=0;fc=xb;gc=Ab}else if((aa|0)==269){h:do{if((_|0)!=0){na=_+1|0;ma=_+8|0;ka=_+4|0;ja=1;i:while(1){S=a[_]|0;if((S&1)==0){hc=(S&255)>>>1}else{hc=c[ka>>2]|0}if(!(ja>>>0<hc>>>0)){break h}S=c[e>>2]|0;do{if((S|0)!=0){if((c[S+12>>2]|0)==(c[S+16>>2]|0)){if((mc[c[(c[S>>2]|0)+36>>2]&63](S)|0)==-1){c[e>>2]=0;ic=0;break}else{ic=c[e>>2]|0;break}}else{ic=S}}else{ic=0}}while(0);S=(ic|0)==0;U=c[f>>2]|0;do{if((U|0)!=0){if((c[U+12>>2]|0)!=(c[U+16>>2]|0)){if(S){break}else{break i}}if(!((mc[c[(c[U>>2]|0)+36>>2]&63](U)|0)==-1)){if(S){break}else{break i}}else{c[f>>2]=0;aa=285;break}}else{aa=285}}while(0);if((aa|0)==285?(aa=0,S):0){break}U=c[e>>2]|0;p=c[U+12>>2]|0;if((p|0)==(c[U+16>>2]|0)){kc=mc[c[(c[U>>2]|0)+36>>2]&63](U)|0}else{kc=d[p]|0}if((a[_]&1)==0){lc=na}else{lc=c[ma>>2]|0}if(!((kc&255)<<24>>24==(a[lc+ja|0]|0))){break}p=ja+1|0;U=c[e>>2]|0;ia=U+12|0;r=c[ia>>2]|0;if((r|0)==(c[U+16>>2]|0)){mc[c[(c[U>>2]|0)+40>>2]&63](U)|0;ja=p;continue}else{c[ia>>2]=r+1;ja=p;continue}}c[k>>2]=c[k>>2]|4;ec=0;fc=Y;gc=$;break g}}while(0);if((Y|0)!=(Z|0)){c[B>>2]=0;nj(v,Y,Z,B);if((c[B>>2]|0)==0){ec=1;fc=Y;gc=$}else{c[k>>2]=c[k>>2]|4;ec=0;fc=Y;gc=$}}else{ec=1;fc=Z;gc=$}}}while(0);ze(z);ze(y);ze(x);ze(w);ze(v);if((fc|0)==0){i=q;return ec|0}jc[gc&127](fc);i=q;return ec|0}function jj(a){a=a|0;var b=0;b=Wa(8)|0;ge(b,a);xb(b|0,5264,11)}function kj(b,d,e,f,g,h,j,k){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0;d=i;i=i+144|0;l=d;m=d+40|0;n=d+16|0;o=d+24|0;p=d+28|0;q=d+36|0;r=d+32|0;c[n>>2]=m;s=n+4|0;c[s>>2]=114;t=m+100|0;Te(p,h);m=c[p>>2]|0;if(!((c[2064]|0)==-1)){c[l>>2]=8256;c[l+4>>2]=112;c[l+8>>2]=0;ue(8256,l,113)}u=(c[8260>>2]|0)+ -1|0;v=c[m+8>>2]|0;if(!((c[m+12>>2]|0)-v>>2>>>0>u>>>0)){w=Wa(4)|0;qm(w);xb(w|0,16216,101)}m=c[v+(u<<2)>>2]|0;if((m|0)==0){w=Wa(4)|0;qm(w);xb(w|0,16216,101)}a[q]=0;w=c[f>>2]|0;c[r>>2]=w;u=c[h+4>>2]|0;c[l+0>>2]=c[r+0>>2];if(ij(e,l,g,p,u,j,q,m,n,o,t)|0){if((a[k]&1)==0){a[k+1|0]=0;a[k]=0}else{a[c[k+8>>2]|0]=0;c[k+4>>2]=0}if((a[q]|0)!=0){Ee(k,vc[c[(c[m>>2]|0)+28>>2]&15](m,45)|0)}q=vc[c[(c[m>>2]|0)+28>>2]&15](m,48)|0;m=c[n>>2]|0;t=c[o>>2]|0;o=t+ -1|0;a:do{if(m>>>0<o>>>0){u=m;while(1){g=u+1|0;if(!((a[u]|0)==q<<24>>24)){x=u;break a}if(g>>>0<o>>>0){u=g}else{x=g;break}}}else{x=m}}while(0);lj(k,x,t)|0}t=c[e>>2]|0;if((t|0)!=0){if((c[t+12>>2]|0)==(c[t+16>>2]|0)?(mc[c[(c[t>>2]|0)+36>>2]&63](t)|0)==-1:0){c[e>>2]=0;y=0}else{y=t}}else{y=0}t=(y|0)==0;do{if((w|0)!=0){if((c[w+12>>2]|0)!=(c[w+16>>2]|0)){if(t){break}else{z=27;break}}if(!((mc[c[(c[w>>2]|0)+36>>2]&63](w)|0)==-1)){if(t^(w|0)==0){break}else{z=27;break}}else{c[f>>2]=0;z=25;break}}else{z=25}}while(0);if((z|0)==25?t:0){z=27}if((z|0)==27){c[j>>2]=c[j>>2]|2}c[b>>2]=y;ae(c[p>>2]|0)|0;p=c[n>>2]|0;c[n>>2]=0;if((p|0)==0){i=d;return}jc[c[s>>2]&127](p);i=d;return}function lj(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0;f=i;g=d;h=a[b]|0;if((h&1)==0){j=10;k=(h&255)>>>1;l=h}else{h=c[b>>2]|0;j=(h&-2)+ -1|0;k=c[b+4>>2]|0;l=h&255}h=e-g|0;if((e|0)==(d|0)){i=f;return b|0}if((j-k|0)>>>0<h>>>0){He(b,j,k+h-j|0,k,k,0,0);m=a[b]|0}else{m=l}if((m&1)==0){n=b+1|0}else{n=c[b+8>>2]|0}m=e+(k-g)|0;g=d;d=n+k|0;while(1){a[d]=a[g]|0;g=g+1|0;if((g|0)==(e|0)){break}else{d=d+1|0}}a[n+m|0]=0;m=k+h|0;if((a[b]&1)==0){a[b]=m<<1;i=f;return b|0}else{c[b+4>>2]=m;i=f;return b|0}return 0}



function mj(b,d,e,f,g,h,j,k,l,m){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;m=m|0;var n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0;n=i;i=i+128|0;o=n;p=n+100|0;q=n+88|0;r=n+76|0;s=n+64|0;t=n+104|0;u=n+36|0;v=n+24|0;w=n+12|0;x=n+40|0;y=n+52|0;if(b){b=c[d>>2]|0;if(!((c[1924]|0)==-1)){c[o>>2]=7696;c[o+4>>2]=112;c[o+8>>2]=0;ue(7696,o,113)}z=(c[7700>>2]|0)+ -1|0;A=c[b+8>>2]|0;if(!((c[b+12>>2]|0)-A>>2>>>0>z>>>0)){B=Wa(4)|0;qm(B);xb(B|0,16216,101)}b=c[A+(z<<2)>>2]|0;if((b|0)==0){B=Wa(4)|0;qm(B);xb(B|0,16216,101)}kc[c[(c[b>>2]|0)+44>>2]&63](p,b);B=c[p>>2]|0;a[e]=B;a[e+1|0]=B>>8;a[e+2|0]=B>>16;a[e+3|0]=B>>24;kc[c[(c[b>>2]|0)+32>>2]&63](q,b);if((a[l]&1)==0){a[l+1|0]=0;a[l]=0}else{a[c[l+8>>2]|0]=0;c[l+4>>2]=0}De(l,0);c[l+0>>2]=c[q+0>>2];c[l+4>>2]=c[q+4>>2];c[l+8>>2]=c[q+8>>2];c[q+0>>2]=0;c[q+4>>2]=0;c[q+8>>2]=0;ze(q);kc[c[(c[b>>2]|0)+28>>2]&63](r,b);if((a[k]&1)==0){a[k+1|0]=0;a[k]=0}else{a[c[k+8>>2]|0]=0;c[k+4>>2]=0}De(k,0);c[k+0>>2]=c[r+0>>2];c[k+4>>2]=c[r+4>>2];c[k+8>>2]=c[r+8>>2];c[r+0>>2]=0;c[r+4>>2]=0;c[r+8>>2]=0;ze(r);a[f]=mc[c[(c[b>>2]|0)+12>>2]&63](b)|0;a[g]=mc[c[(c[b>>2]|0)+16>>2]&63](b)|0;kc[c[(c[b>>2]|0)+20>>2]&63](s,b);if((a[h]&1)==0){a[h+1|0]=0;a[h]=0}else{a[c[h+8>>2]|0]=0;c[h+4>>2]=0}De(h,0);c[h+0>>2]=c[s+0>>2];c[h+4>>2]=c[s+4>>2];c[h+8>>2]=c[s+8>>2];c[s+0>>2]=0;c[s+4>>2]=0;c[s+8>>2]=0;ze(s);kc[c[(c[b>>2]|0)+24>>2]&63](t,b);if((a[j]&1)==0){a[j+1|0]=0;a[j]=0}else{a[c[j+8>>2]|0]=0;c[j+4>>2]=0}De(j,0);c[j+0>>2]=c[t+0>>2];c[j+4>>2]=c[t+4>>2];c[j+8>>2]=c[t+8>>2];c[t+0>>2]=0;c[t+4>>2]=0;c[t+8>>2]=0;ze(t);C=mc[c[(c[b>>2]|0)+36>>2]&63](b)|0;c[m>>2]=C;i=n;return}else{b=c[d>>2]|0;if(!((c[1908]|0)==-1)){c[o>>2]=7632;c[o+4>>2]=112;c[o+8>>2]=0;ue(7632,o,113)}o=(c[7636>>2]|0)+ -1|0;d=c[b+8>>2]|0;if(!((c[b+12>>2]|0)-d>>2>>>0>o>>>0)){D=Wa(4)|0;qm(D);xb(D|0,16216,101)}b=c[d+(o<<2)>>2]|0;if((b|0)==0){D=Wa(4)|0;qm(D);xb(D|0,16216,101)}kc[c[(c[b>>2]|0)+44>>2]&63](u,b);D=c[u>>2]|0;a[e]=D;a[e+1|0]=D>>8;a[e+2|0]=D>>16;a[e+3|0]=D>>24;kc[c[(c[b>>2]|0)+32>>2]&63](v,b);if((a[l]&1)==0){a[l+1|0]=0;a[l]=0}else{a[c[l+8>>2]|0]=0;c[l+4>>2]=0}De(l,0);c[l+0>>2]=c[v+0>>2];c[l+4>>2]=c[v+4>>2];c[l+8>>2]=c[v+8>>2];c[v+0>>2]=0;c[v+4>>2]=0;c[v+8>>2]=0;ze(v);kc[c[(c[b>>2]|0)+28>>2]&63](w,b);if((a[k]&1)==0){a[k+1|0]=0;a[k]=0}else{a[c[k+8>>2]|0]=0;c[k+4>>2]=0}De(k,0);c[k+0>>2]=c[w+0>>2];c[k+4>>2]=c[w+4>>2];c[k+8>>2]=c[w+8>>2];c[w+0>>2]=0;c[w+4>>2]=0;c[w+8>>2]=0;ze(w);a[f]=mc[c[(c[b>>2]|0)+12>>2]&63](b)|0;a[g]=mc[c[(c[b>>2]|0)+16>>2]&63](b)|0;kc[c[(c[b>>2]|0)+20>>2]&63](x,b);if((a[h]&1)==0){a[h+1|0]=0;a[h]=0}else{a[c[h+8>>2]|0]=0;c[h+4>>2]=0}De(h,0);c[h+0>>2]=c[x+0>>2];c[h+4>>2]=c[x+4>>2];c[h+8>>2]=c[x+8>>2];c[x+0>>2]=0;c[x+4>>2]=0;c[x+8>>2]=0;ze(x);kc[c[(c[b>>2]|0)+24>>2]&63](y,b);if((a[j]&1)==0){a[j+1|0]=0;a[j]=0}else{a[c[j+8>>2]|0]=0;c[j+4>>2]=0}De(j,0);c[j+0>>2]=c[y+0>>2];c[j+4>>2]=c[y+4>>2];c[j+8>>2]=c[y+8>>2];c[y+0>>2]=0;c[y+4>>2]=0;c[y+8>>2]=0;ze(y);C=mc[c[(c[b>>2]|0)+36>>2]&63](b)|0;c[m>>2]=C;i=n;return}}function nj(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0;g=i;h=a[b]|0;if((h&1)==0){j=(h&255)>>>1}else{j=c[b+4>>2]|0}if((j|0)==0){i=g;return}if((d|0)!=(e|0)?(j=e+ -4|0,j>>>0>d>>>0):0){k=d;l=j;do{j=c[k>>2]|0;c[k>>2]=c[l>>2];c[l>>2]=j;k=k+4|0;l=l+ -4|0}while(k>>>0<l>>>0);m=a[b]|0}else{m=h}if((m&1)==0){n=b+1|0;o=(m&255)>>>1}else{n=c[b+8>>2]|0;o=c[b+4>>2]|0}b=e+ -4|0;e=a[n]|0;m=e<<24>>24<1|e<<24>>24==127;a:do{if(b>>>0>d>>>0){h=n+o|0;l=e;k=n;j=d;p=m;while(1){if(!p?(l<<24>>24|0)!=(c[j>>2]|0):0){break}q=(h-k|0)>1?k+1|0:k;r=j+4|0;s=a[q]|0;t=s<<24>>24<1|s<<24>>24==127;if(r>>>0<b>>>0){l=s;k=q;j=r;p=t}else{u=s;v=t;break a}}c[f>>2]=4;i=g;return}else{u=e;v=m}}while(0);if(v){i=g;return}v=c[b>>2]|0;if(!(u<<24>>24>>>0<v>>>0|(v|0)==0)){i=g;return}c[f>>2]=4;i=g;return}function oj(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function pj(a){a=a|0;return}function qj(b,d,e,f,g,h,j,k){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0;d=i;i=i+576|0;l=d;m=d+72|0;n=d+64|0;o=d+60|0;p=d+56|0;q=d+572|0;r=d+12|0;s=d+16|0;t=d+472|0;c[n>>2]=m;u=n+4|0;c[u>>2]=114;v=m+400|0;Te(p,h);m=c[p>>2]|0;if(!((c[2062]|0)==-1)){c[l>>2]=8248;c[l+4>>2]=112;c[l+8>>2]=0;ue(8248,l,113)}w=(c[8252>>2]|0)+ -1|0;x=c[m+8>>2]|0;if(!((c[m+12>>2]|0)-x>>2>>>0>w>>>0)){y=Wa(4)|0;qm(y);xb(y|0,16216,101)}m=c[x+(w<<2)>>2]|0;if((m|0)==0){y=Wa(4)|0;qm(y);xb(y|0,16216,101)}a[q]=0;c[r>>2]=c[f>>2];y=c[h+4>>2]|0;c[l+0>>2]=c[r+0>>2];if(rj(e,l,g,p,y,j,q,m,n,o,v)|0){sc[c[(c[m>>2]|0)+48>>2]&7](m,7952,7962|0,s)|0;m=c[o>>2]|0;v=c[n>>2]|0;y=m-v|0;if((y|0)>392){g=Lm((y>>2)+2|0)|0;if((g|0)==0){_m()}else{z=g;A=g}}else{z=0;A=t}if((a[q]|0)==0){B=A}else{a[A]=45;B=A+1|0}if(v>>>0<m>>>0){m=s+40|0;A=s;q=B;g=v;while(1){v=c[g>>2]|0;y=s;while(1){r=y+4|0;if((c[y>>2]|0)==(v|0)){C=y;break}if((r|0)==(m|0)){C=m;break}else{y=r}}a[q]=a[7952+(C-A>>2)|0]|0;y=g+4|0;v=q+1|0;if(y>>>0<(c[o>>2]|0)>>>0){q=v;g=y}else{D=v;break}}}else{D=B}a[D]=0;c[l>>2]=k;if((gb(t|0,7888,l|0)|0)!=1){l=Wa(8)|0;ge(l,7896);xb(l|0,5264,11)}if((z|0)!=0){Mm(z)}}z=c[e>>2]|0;do{if((z|0)!=0){l=c[z+12>>2]|0;if((l|0)==(c[z+16>>2]|0)){E=mc[c[(c[z>>2]|0)+36>>2]&63](z)|0}else{E=c[l>>2]|0}if((E|0)==-1){c[e>>2]=0;F=1;break}else{F=(c[e>>2]|0)==0;break}}else{F=1}}while(0);E=c[f>>2]|0;do{if((E|0)!=0){z=c[E+12>>2]|0;if((z|0)==(c[E+16>>2]|0)){G=mc[c[(c[E>>2]|0)+36>>2]&63](E)|0}else{G=c[z>>2]|0}if(!((G|0)==-1)){if(F){break}else{H=37;break}}else{c[f>>2]=0;H=35;break}}else{H=35}}while(0);if((H|0)==35?F:0){H=37}if((H|0)==37){c[j>>2]=c[j>>2]|2}c[b>>2]=c[e>>2];ae(c[p>>2]|0)|0;p=c[n>>2]|0;c[n>>2]=0;if((p|0)==0){i=d;return}jc[c[u>>2]&127](p);i=d;return}function rj(b,e,f,g,h,j,k,l,m,n,o){b=b|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;m=m|0;n=n|0;o=o|0;var p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0,X=0,Y=0,Z=0,_=0,$=0,aa=0,ba=0,ca=0,da=0,ea=0,fa=0,ga=0,ha=0,ia=0,ja=0,ka=0,la=0,ma=0,na=0,oa=0,pa=0,qa=0,ra=0,sa=0,ta=0,ua=0,va=0,wa=0,xa=0,ya=0,za=0,Aa=0,Ba=0,Ca=0,Da=0,Ea=0,Fa=0,Ga=0,Ha=0,Ia=0,Ja=0,Ka=0,La=0,Ma=0,Na=0,Oa=0,Pa=0,Qa=0,Ra=0,Sa=0,Ta=0,Ua=0,Va=0,Wa=0,Xa=0,Ya=0,Za=0,_a=0,$a=0,ab=0,bb=0,cb=0,db=0,eb=0,fb=0,gb=0,hb=0,ib=0,jb=0,kb=0,lb=0,mb=0,nb=0,ob=0,pb=0,qb=0,rb=0,sb=0,tb=0,ub=0,vb=0,wb=0,xb=0,yb=0,zb=0,Ab=0,Bb=0,Cb=0,Db=0,Eb=0,Fb=0,Gb=0,Hb=0,Ib=0,Jb=0,Kb=0,Lb=0,Mb=0,Nb=0,Ob=0,Pb=0,Qb=0,Rb=0,Sb=0,Tb=0,Ub=0,Vb=0,Wb=0,Xb=0,Yb=0,Zb=0,_b=0,$b=0,ac=0,bc=0,cc=0,dc=0,ec=0,fc=0,hc=0,ic=0,kc=0,lc=0,nc=0,oc=0,pc=0,qc=0,rc=0,sc=0,tc=0,uc=0,vc=0,wc=0;p=i;i=i+480|0;q=p+80|0;r=p+76|0;s=p+72|0;t=p+68|0;u=p+56|0;v=p+44|0;w=p+28|0;x=p+16|0;y=p+4|0;z=p;A=p+40|0;c[r>>2]=0;c[u+0>>2]=0;c[u+4>>2]=0;c[u+8>>2]=0;c[v+0>>2]=0;c[v+4>>2]=0;c[v+8>>2]=0;c[w+0>>2]=0;c[w+4>>2]=0;c[w+8>>2]=0;c[x+0>>2]=0;c[x+4>>2]=0;c[x+8>>2]=0;c[y+0>>2]=0;c[y+4>>2]=0;c[y+8>>2]=0;uj(f,g,r,s,t,u,v,w,x,z);c[n>>2]=c[m>>2];g=x+4|0;f=x+8|0;B=w+4|0;C=w+8|0;D=(h&512|0)!=0;h=v+4|0;E=v+8|0;F=y+4|0;G=y+8|0;H=r+3|0;I=m+4|0;J=u+4|0;K=q+400|0;L=q;M=q;q=o;o=0;N=0;O=114;a:while(1){P=c[b>>2]|0;do{if((P|0)!=0){Q=c[P+12>>2]|0;if((Q|0)==(c[P+16>>2]|0)){R=mc[c[(c[P>>2]|0)+36>>2]&63](P)|0}else{R=c[Q>>2]|0}if((R|0)==-1){c[b>>2]=0;S=1;break}else{S=(c[b>>2]|0)==0;break}}else{S=1}}while(0);P=c[e>>2]|0;do{if((P|0)!=0){Q=c[P+12>>2]|0;if((Q|0)==(c[P+16>>2]|0)){T=mc[c[(c[P>>2]|0)+36>>2]&63](P)|0}else{T=c[Q>>2]|0}if(!((T|0)==-1)){if(S){U=P;break}else{V=M;W=L;X=N;Y=O;Z=292;break a}}else{c[e>>2]=0;Z=15;break}}else{Z=15}}while(0);if((Z|0)==15){Z=0;if(S){V=M;W=L;X=N;Y=O;Z=292;break}else{U=0}}b:do{switch(a[r+o|0]|0){case 1:{if((o|0)==3){V=M;W=L;X=N;Y=O;Z=292;break a}P=c[b>>2]|0;Q=c[P+12>>2]|0;if((Q|0)==(c[P+16>>2]|0)){_=mc[c[(c[P>>2]|0)+36>>2]&63](P)|0}else{_=c[Q>>2]|0}if(!(gc[c[(c[l>>2]|0)+12>>2]&31](l,8192,_)|0)){Z=27;break a}Q=c[b>>2]|0;P=Q+12|0;$=c[P>>2]|0;if(($|0)==(c[Q+16>>2]|0)){aa=mc[c[(c[Q>>2]|0)+40>>2]&63](Q)|0}else{c[P>>2]=$+4;aa=c[$>>2]|0}Oe(y,aa);Z=28;break};case 3:{$=a[w]|0;P=($&1)==0;if(P){ba=($&255)>>>1}else{ba=c[B>>2]|0}Q=a[x]|0;ca=(Q&1)==0;if(ca){da=(Q&255)>>>1}else{da=c[g>>2]|0}if((ba|0)==(0-da|0)){ea=q;fa=M;ga=L;ha=K;ia=N;ja=O}else{if(P){ka=($&255)>>>1}else{ka=c[B>>2]|0}if((ka|0)!=0){if(ca){la=(Q&255)>>>1}else{la=c[g>>2]|0}if((la|0)!=0){ca=c[b>>2]|0;ma=c[ca+12>>2]|0;if((ma|0)==(c[ca+16>>2]|0)){na=mc[c[(c[ca>>2]|0)+36>>2]&63](ca)|0;oa=na;pa=a[w]|0}else{oa=c[ma>>2]|0;pa=$}ma=c[b>>2]|0;na=ma+12|0;ca=c[na>>2]|0;qa=(ca|0)==(c[ma+16>>2]|0);if((oa|0)==(c[((pa&1)==0?B:c[C>>2]|0)>>2]|0)){if(qa){mc[c[(c[ma>>2]|0)+40>>2]&63](ma)|0}else{c[na>>2]=ca+4}na=a[w]|0;if((na&1)==0){ra=(na&255)>>>1}else{ra=c[B>>2]|0}ea=q;fa=M;ga=L;ha=K;ia=ra>>>0>1?w:N;ja=O;break b}if(qa){sa=mc[c[(c[ma>>2]|0)+36>>2]&63](ma)|0}else{sa=c[ca>>2]|0}if((sa|0)!=(c[((a[x]&1)==0?g:c[f>>2]|0)>>2]|0)){Z=116;break a}ca=c[b>>2]|0;ma=ca+12|0;qa=c[ma>>2]|0;if((qa|0)==(c[ca+16>>2]|0)){mc[c[(c[ca>>2]|0)+40>>2]&63](ca)|0}else{c[ma>>2]=qa+4}a[k]=1;qa=a[x]|0;if((qa&1)==0){ta=(qa&255)>>>1}else{ta=c[g>>2]|0}ea=q;fa=M;ga=L;ha=K;ia=ta>>>0>1?x:N;ja=O;break b}}if(P){ua=($&255)>>>1}else{ua=c[B>>2]|0}P=c[b>>2]|0;qa=c[P+12>>2]|0;ma=(qa|0)==(c[P+16>>2]|0);if((ua|0)==0){if(ma){ca=mc[c[(c[P>>2]|0)+36>>2]&63](P)|0;va=ca;wa=a[x]|0}else{va=c[qa>>2]|0;wa=Q}if((va|0)!=(c[((wa&1)==0?g:c[f>>2]|0)>>2]|0)){ea=q;fa=M;ga=L;ha=K;ia=N;ja=O;break b}Q=c[b>>2]|0;ca=Q+12|0;na=c[ca>>2]|0;if((na|0)==(c[Q+16>>2]|0)){mc[c[(c[Q>>2]|0)+40>>2]&63](Q)|0}else{c[ca>>2]=na+4}a[k]=1;na=a[x]|0;if((na&1)==0){xa=(na&255)>>>1}else{xa=c[g>>2]|0}ea=q;fa=M;ga=L;ha=K;ia=xa>>>0>1?x:N;ja=O;break b}if(ma){ma=mc[c[(c[P>>2]|0)+36>>2]&63](P)|0;ya=ma;za=a[w]|0}else{ya=c[qa>>2]|0;za=$}if((ya|0)!=(c[((za&1)==0?B:c[C>>2]|0)>>2]|0)){a[k]=1;ea=q;fa=M;ga=L;ha=K;ia=N;ja=O;break b}$=c[b>>2]|0;qa=$+12|0;ma=c[qa>>2]|0;if((ma|0)==(c[$+16>>2]|0)){mc[c[(c[$>>2]|0)+40>>2]&63]($)|0}else{c[qa>>2]=ma+4}ma=a[w]|0;if((ma&1)==0){Aa=(ma&255)>>>1}else{Aa=c[B>>2]|0}ea=q;fa=M;ga=L;ha=K;ia=Aa>>>0>1?w:N;ja=O}break};case 0:{Z=28;break};case 2:{if(!((N|0)!=0|o>>>0<2)){if((o|0)==2){Ba=(a[H]|0)!=0}else{Ba=0}if(!(D|Ba)){ea=q;fa=M;ga=L;ha=K;ia=0;ja=O;break b}}ma=a[v]|0;qa=(ma&1)==0?h:c[E>>2]|0;c:do{if((o|0)!=0?(d[r+(o+ -1)|0]|0)<2:0){$=ma;P=qa;while(1){if(($&1)==0){Ca=h;Da=($&255)>>>1}else{Ca=c[E>>2]|0;Da=c[h>>2]|0}if((P|0)==(Ca+(Da<<2)|0)){Ea=$;break}if(!(gc[c[(c[l>>2]|0)+12>>2]&31](l,8192,c[P>>2]|0)|0)){Z=129;break}$=a[v]|0;P=P+4|0}if((Z|0)==129){Z=0;Ea=a[v]|0}$=(Ea&1)==0;na=P-($?h:c[E>>2]|0)>>2;ca=a[y]|0;Q=(ca&1)==0;if(Q){Fa=(ca&255)>>>1}else{Fa=c[F>>2]|0}d:do{if(!(na>>>0>Fa>>>0)){if(Q){Ga=F;Ha=(ca&255)>>>1;Ia=F+(((ca&255)>>>1)-na<<2)|0}else{Ja=c[G>>2]|0;Ka=c[F>>2]|0;Ga=Ja;Ha=Ka;Ia=Ja+(Ka-na<<2)|0}Ka=Ga+(Ha<<2)|0;if((Ia|0)==(Ka|0)){La=U;Ma=Ea;Na=P;Oa=U;break c}else{Pa=Ia;Qa=$?h:c[E>>2]|0}while(1){if((c[Pa>>2]|0)!=(c[Qa>>2]|0)){break d}Ja=Pa+4|0;if((Ja|0)==(Ka|0)){La=U;Ma=Ea;Na=P;Oa=U;break c}Pa=Ja;Qa=Qa+4|0}}}while(0);La=U;Ma=Ea;Na=$?h:c[E>>2]|0;Oa=U}else{La=U;Ma=ma;Na=qa;Oa=U}}while(0);e:while(1){if((Ma&1)==0){Ra=h;Sa=(Ma&255)>>>1}else{Ra=c[E>>2]|0;Sa=c[h>>2]|0}if((Na|0)==(Ra+(Sa<<2)|0)){break}qa=c[b>>2]|0;do{if((qa|0)!=0){ma=c[qa+12>>2]|0;if((ma|0)==(c[qa+16>>2]|0)){Ta=mc[c[(c[qa>>2]|0)+36>>2]&63](qa)|0}else{Ta=c[ma>>2]|0}if((Ta|0)==-1){c[b>>2]=0;Ua=1;break}else{Ua=(c[b>>2]|0)==0;break}}else{Ua=1}}while(0);do{if((Oa|0)!=0){qa=c[Oa+12>>2]|0;if((qa|0)==(c[Oa+16>>2]|0)){Va=mc[c[(c[Oa>>2]|0)+36>>2]&63](Oa)|0}else{Va=c[qa>>2]|0}if(!((Va|0)==-1)){if(Ua^(La|0)==0){Wa=La;Xa=La;break}else{break e}}else{c[e>>2]=0;Ya=0;Z=159;break}}else{Ya=La;Z=159}}while(0);if((Z|0)==159){Z=0;if(Ua){break}else{Wa=Ya;Xa=0}}qa=c[b>>2]|0;$=c[qa+12>>2]|0;if(($|0)==(c[qa+16>>2]|0)){Za=mc[c[(c[qa>>2]|0)+36>>2]&63](qa)|0}else{Za=c[$>>2]|0}if((Za|0)!=(c[Na>>2]|0)){break}$=c[b>>2]|0;qa=$+12|0;ma=c[qa>>2]|0;if((ma|0)==(c[$+16>>2]|0)){mc[c[(c[$>>2]|0)+40>>2]&63]($)|0}else{c[qa>>2]=ma+4}La=Wa;Ma=a[v]|0;Na=Na+4|0;Oa=Xa}if(D){ma=a[v]|0;if((ma&1)==0){_a=h;$a=(ma&255)>>>1}else{_a=c[E>>2]|0;$a=c[h>>2]|0}if((Na|0)!=(_a+($a<<2)|0)){Z=174;break a}else{ea=q;fa=M;ga=L;ha=K;ia=N;ja=O}}else{ea=q;fa=M;ga=L;ha=K;ia=N;ja=O}break};case 4:{ma=q;qa=L;$=K;P=M;na=0;ca=O;f:while(1){Q=c[b>>2]|0;do{if((Q|0)!=0){Ka=c[Q+12>>2]|0;if((Ka|0)==(c[Q+16>>2]|0)){ab=mc[c[(c[Q>>2]|0)+36>>2]&63](Q)|0}else{ab=c[Ka>>2]|0}if((ab|0)==-1){c[b>>2]=0;bb=1;break}else{bb=(c[b>>2]|0)==0;break}}else{bb=1}}while(0);Q=c[e>>2]|0;do{if((Q|0)!=0){Ka=c[Q+12>>2]|0;if((Ka|0)==(c[Q+16>>2]|0)){cb=mc[c[(c[Q>>2]|0)+36>>2]&63](Q)|0}else{cb=c[Ka>>2]|0}if(!((cb|0)==-1)){if(bb){break}else{break f}}else{c[e>>2]=0;Z=188;break}}else{Z=188}}while(0);if((Z|0)==188?(Z=0,bb):0){break}Q=c[b>>2]|0;Ka=c[Q+12>>2]|0;if((Ka|0)==(c[Q+16>>2]|0)){db=mc[c[(c[Q>>2]|0)+36>>2]&63](Q)|0}else{db=c[Ka>>2]|0}if(gc[c[(c[l>>2]|0)+12>>2]&31](l,2048,db)|0){Ka=c[n>>2]|0;if((Ka|0)==(ma|0)){Q=(c[I>>2]|0)!=114;Ja=c[m>>2]|0;eb=ma-Ja|0;fb=eb>>>0<2147483647?eb<<1:-1;gb=eb>>2;if(Q){hb=Ja}else{hb=0}Ja=Nm(hb,fb)|0;if((Ja|0)==0){Z=198;break a}if(!Q){Q=c[m>>2]|0;c[m>>2]=Ja;if((Q|0)==0){ib=Ja}else{jc[c[I>>2]&127](Q);ib=c[m>>2]|0}}else{c[m>>2]=Ja;ib=Ja}c[I>>2]=115;Ja=ib+(gb<<2)|0;c[n>>2]=Ja;jb=Ja;kb=(c[m>>2]|0)+(fb>>>2<<2)|0}else{jb=Ka;kb=ma}c[n>>2]=jb+4;c[jb>>2]=db;lb=kb;mb=P;nb=qa;ob=$;pb=na+1|0;qb=ca}else{Ka=a[u]|0;if((Ka&1)==0){rb=(Ka&255)>>>1}else{rb=c[J>>2]|0}if((rb|0)==0|(na|0)==0){break}if((db|0)!=(c[t>>2]|0)){break}if((qa|0)==($|0)){Ka=qa-P|0;fb=Ka>>>0<2147483647?Ka<<1:-1;if((ca|0)!=114){sb=P}else{sb=0}Ja=Nm(sb,fb)|0;if((Ja|0)==0){Z=214;break a}tb=Ja+(Ka>>2<<2)|0;ub=Ja;vb=Ja+(fb>>>2<<2)|0;wb=115}else{tb=qa;ub=P;vb=$;wb=ca}c[tb>>2]=na;lb=ma;mb=ub;nb=tb+4|0;ob=vb;pb=0;qb=wb}fb=c[b>>2]|0;Ja=fb+12|0;Ka=c[Ja>>2]|0;if((Ka|0)==(c[fb+16>>2]|0)){mc[c[(c[fb>>2]|0)+40>>2]&63](fb)|0;ma=lb;qa=nb;$=ob;P=mb;na=pb;ca=qb;continue}else{c[Ja>>2]=Ka+4;ma=lb;qa=nb;$=ob;P=mb;na=pb;ca=qb;continue}}if((P|0)==(qa|0)|(na|0)==0){xb=P;yb=qa;zb=$;Ab=ca}else{if((qa|0)==($|0)){Ka=qa-P|0;Ja=Ka>>>0<2147483647?Ka<<1:-1;if((ca|0)!=114){Bb=P}else{Bb=0}fb=Nm(Bb,Ja)|0;if((fb|0)==0){Z=225;break a}Cb=fb+(Ka>>2<<2)|0;Db=fb;Eb=fb+(Ja>>>2<<2)|0;Fb=115}else{Cb=qa;Db=P;Eb=$;Fb=ca}c[Cb>>2]=na;xb=Db;yb=Cb+4|0;zb=Eb;Ab=Fb}Ja=c[z>>2]|0;if((Ja|0)>0){fb=c[b>>2]|0;do{if((fb|0)!=0){Ka=c[fb+12>>2]|0;if((Ka|0)==(c[fb+16>>2]|0)){Gb=mc[c[(c[fb>>2]|0)+36>>2]&63](fb)|0}else{Gb=c[Ka>>2]|0}if((Gb|0)==-1){c[b>>2]=0;Hb=1;break}else{Hb=(c[b>>2]|0)==0;break}}else{Hb=1}}while(0);fb=c[e>>2]|0;do{if((fb|0)!=0){na=c[fb+12>>2]|0;if((na|0)==(c[fb+16>>2]|0)){Ib=mc[c[(c[fb>>2]|0)+36>>2]&63](fb)|0}else{Ib=c[na>>2]|0}if(!((Ib|0)==-1)){if(Hb){Jb=fb;break}else{Z=248;break a}}else{c[e>>2]=0;Z=242;break}}else{Z=242}}while(0);if((Z|0)==242){Z=0;if(Hb){Z=248;break a}else{Jb=0}}fb=c[b>>2]|0;na=c[fb+12>>2]|0;if((na|0)==(c[fb+16>>2]|0)){Kb=mc[c[(c[fb>>2]|0)+36>>2]&63](fb)|0}else{Kb=c[na>>2]|0}if((Kb|0)!=(c[s>>2]|0)){Z=248;break a}na=c[b>>2]|0;fb=na+12|0;ca=c[fb>>2]|0;if((ca|0)==(c[na+16>>2]|0)){mc[c[(c[na>>2]|0)+40>>2]&63](na)|0;Lb=Jb;Mb=Jb;Nb=ma;Ob=Ja}else{c[fb>>2]=ca+4;Lb=Jb;Mb=Jb;Nb=ma;Ob=Ja}while(1){ca=c[b>>2]|0;do{if((ca|0)!=0){fb=c[ca+12>>2]|0;if((fb|0)==(c[ca+16>>2]|0)){Pb=mc[c[(c[ca>>2]|0)+36>>2]&63](ca)|0}else{Pb=c[fb>>2]|0}if((Pb|0)==-1){c[b>>2]=0;Qb=1;break}else{Qb=(c[b>>2]|0)==0;break}}else{Qb=1}}while(0);do{if((Mb|0)!=0){ca=c[Mb+12>>2]|0;if((ca|0)==(c[Mb+16>>2]|0)){Rb=mc[c[(c[Mb>>2]|0)+36>>2]&63](Mb)|0}else{Rb=c[ca>>2]|0}if(!((Rb|0)==-1)){if(Qb^(Lb|0)==0){Sb=Lb;Tb=Lb;break}else{Z=271;break a}}else{c[e>>2]=0;Ub=0;Z=265;break}}else{Ub=Lb;Z=265}}while(0);if((Z|0)==265){Z=0;if(Qb){Z=271;break a}else{Sb=Ub;Tb=0}}ca=c[b>>2]|0;fb=c[ca+12>>2]|0;if((fb|0)==(c[ca+16>>2]|0)){Vb=mc[c[(c[ca>>2]|0)+36>>2]&63](ca)|0}else{Vb=c[fb>>2]|0}if(!(gc[c[(c[l>>2]|0)+12>>2]&31](l,2048,Vb)|0)){Z=271;break a}fb=c[n>>2]|0;if((fb|0)==(Nb|0)){ca=(c[I>>2]|0)!=114;na=c[m>>2]|0;$=Nb-na|0;P=$>>>0<2147483647?$<<1:-1;qa=$>>2;if(ca){Wb=na}else{Wb=0}na=Nm(Wb,P)|0;if((na|0)==0){Z=276;break a}if(!ca){ca=c[m>>2]|0;c[m>>2]=na;if((ca|0)==0){Xb=na}else{jc[c[I>>2]&127](ca);Xb=c[m>>2]|0}}else{c[m>>2]=na;Xb=na}c[I>>2]=115;na=Xb+(qa<<2)|0;c[n>>2]=na;Yb=na;Zb=(c[m>>2]|0)+(P>>>2<<2)|0}else{Yb=fb;Zb=Nb}fb=c[b>>2]|0;P=c[fb+12>>2]|0;if((P|0)==(c[fb+16>>2]|0)){na=mc[c[(c[fb>>2]|0)+36>>2]&63](fb)|0;_b=na;$b=c[n>>2]|0}else{_b=c[P>>2]|0;$b=Yb}c[n>>2]=$b+4;c[$b>>2]=_b;P=Ob+ -1|0;c[z>>2]=P;na=c[b>>2]|0;fb=na+12|0;qa=c[fb>>2]|0;if((qa|0)==(c[na+16>>2]|0)){mc[c[(c[na>>2]|0)+40>>2]&63](na)|0}else{c[fb>>2]=qa+4}if((P|0)>0){Lb=Sb;Mb=Tb;Nb=Zb;Ob=P}else{ac=Zb;break}}}else{ac=ma}if((c[n>>2]|0)==(c[m>>2]|0)){Z=290;break a}else{ea=ac;fa=xb;ga=yb;ha=zb;ia=N;ja=Ab}break};default:{ea=q;fa=M;ga=L;ha=K;ia=N;ja=O}}}while(0);g:do{if((Z|0)==28){Z=0;if((o|0)==3){V=M;W=L;X=N;Y=O;Z=292;break a}else{bc=U;cc=U}while(1){Ja=c[b>>2]|0;do{if((Ja|0)!=0){P=c[Ja+12>>2]|0;if((P|0)==(c[Ja+16>>2]|0)){dc=mc[c[(c[Ja>>2]|0)+36>>2]&63](Ja)|0}else{dc=c[P>>2]|0}if((dc|0)==-1){c[b>>2]=0;ec=1;break}else{ec=(c[b>>2]|0)==0;break}}else{ec=1}}while(0);do{if((cc|0)!=0){Ja=c[cc+12>>2]|0;if((Ja|0)==(c[cc+16>>2]|0)){fc=mc[c[(c[cc>>2]|0)+36>>2]&63](cc)|0}else{fc=c[Ja>>2]|0}if(!((fc|0)==-1)){if(ec^(bc|0)==0){hc=bc;ic=bc;break}else{ea=q;fa=M;ga=L;ha=K;ia=N;ja=O;break g}}else{c[e>>2]=0;kc=0;Z=42;break}}else{kc=bc;Z=42}}while(0);if((Z|0)==42){Z=0;if(ec){ea=q;fa=M;ga=L;ha=K;ia=N;ja=O;break g}else{hc=kc;ic=0}}Ja=c[b>>2]|0;P=c[Ja+12>>2]|0;if((P|0)==(c[Ja+16>>2]|0)){lc=mc[c[(c[Ja>>2]|0)+36>>2]&63](Ja)|0}else{lc=c[P>>2]|0}if(!(gc[c[(c[l>>2]|0)+12>>2]&31](l,8192,lc)|0)){ea=q;fa=M;ga=L;ha=K;ia=N;ja=O;break g}P=c[b>>2]|0;Ja=P+12|0;qa=c[Ja>>2]|0;if((qa|0)==(c[P+16>>2]|0)){nc=mc[c[(c[P>>2]|0)+40>>2]&63](P)|0}else{c[Ja>>2]=qa+4;nc=c[qa>>2]|0}Oe(y,nc);bc=hc;cc=ic}}}while(0);ma=o+1|0;if(ma>>>0<4){K=ha;L=ga;M=fa;q=ea;o=ma;N=ia;O=ja}else{V=fa;W=ga;X=ia;Y=ja;Z=292;break}}h:do{if((Z|0)==27){c[j>>2]=c[j>>2]|4;oc=0;pc=M;qc=O}else if((Z|0)==116){c[j>>2]=c[j>>2]|4;oc=0;pc=M;qc=O}else if((Z|0)==174){c[j>>2]=c[j>>2]|4;oc=0;pc=M;qc=O}else if((Z|0)==198){_m()}else if((Z|0)==214){_m()}else if((Z|0)==225){_m()}else if((Z|0)==248){c[j>>2]=c[j>>2]|4;oc=0;pc=xb;qc=Ab}else if((Z|0)==271){c[j>>2]=c[j>>2]|4;oc=0;pc=xb;qc=Ab}else if((Z|0)==276){_m()}else if((Z|0)==290){c[j>>2]=c[j>>2]|4;oc=0;pc=xb;qc=Ab}else if((Z|0)==292){i:do{if((X|0)!=0){ja=X+4|0;ia=X+8|0;ga=1;j:while(1){fa=a[X]|0;if((fa&1)==0){rc=(fa&255)>>>1}else{rc=c[ja>>2]|0}if(!(ga>>>0<rc>>>0)){break i}fa=c[b>>2]|0;do{if((fa|0)!=0){N=c[fa+12>>2]|0;if((N|0)==(c[fa+16>>2]|0)){sc=mc[c[(c[fa>>2]|0)+36>>2]&63](fa)|0}else{sc=c[N>>2]|0}if((sc|0)==-1){c[b>>2]=0;tc=1;break}else{tc=(c[b>>2]|0)==0;break}}else{tc=1}}while(0);fa=c[e>>2]|0;do{if((fa|0)!=0){N=c[fa+12>>2]|0;if((N|0)==(c[fa+16>>2]|0)){uc=mc[c[(c[fa>>2]|0)+36>>2]&63](fa)|0}else{uc=c[N>>2]|0}if(!((uc|0)==-1)){if(tc){break}else{break j}}else{c[e>>2]=0;Z=311;break}}else{Z=311}}while(0);if((Z|0)==311?(Z=0,tc):0){break}fa=c[b>>2]|0;N=c[fa+12>>2]|0;if((N|0)==(c[fa+16>>2]|0)){vc=mc[c[(c[fa>>2]|0)+36>>2]&63](fa)|0}else{vc=c[N>>2]|0}if((a[X]&1)==0){wc=ja}else{wc=c[ia>>2]|0}if((vc|0)!=(c[wc+(ga<<2)>>2]|0)){break}N=ga+1|0;fa=c[b>>2]|0;o=fa+12|0;ea=c[o>>2]|0;if((ea|0)==(c[fa+16>>2]|0)){mc[c[(c[fa>>2]|0)+40>>2]&63](fa)|0;ga=N;continue}else{c[o>>2]=ea+4;ga=N;continue}}c[j>>2]=c[j>>2]|4;oc=0;pc=V;qc=Y;break h}}while(0);if((V|0)!=(W|0)){c[A>>2]=0;nj(u,V,W,A);if((c[A>>2]|0)==0){oc=1;pc=V;qc=Y}else{c[j>>2]=c[j>>2]|4;oc=0;pc=V;qc=Y}}else{oc=1;pc=W;qc=Y}}}while(0);Ke(y);Ke(x);Ke(w);Ke(v);ze(u);if((pc|0)==0){i=p;return oc|0}jc[qc&127](pc);i=p;return oc|0}function sj(b,d,e,f,g,h,j,k){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0;d=i;i=i+448|0;l=d;m=d+40|0;n=d+16|0;o=d+24|0;p=d+28|0;q=d+440|0;r=d+32|0;c[n>>2]=m;s=n+4|0;c[s>>2]=114;t=m+400|0;Te(p,h);m=c[p>>2]|0;if(!((c[2062]|0)==-1)){c[l>>2]=8248;c[l+4>>2]=112;c[l+8>>2]=0;ue(8248,l,113)}u=(c[8252>>2]|0)+ -1|0;v=c[m+8>>2]|0;if(!((c[m+12>>2]|0)-v>>2>>>0>u>>>0)){w=Wa(4)|0;qm(w);xb(w|0,16216,101)}m=c[v+(u<<2)>>2]|0;if((m|0)==0){w=Wa(4)|0;qm(w);xb(w|0,16216,101)}a[q]=0;w=c[f>>2]|0;c[r>>2]=w;u=c[h+4>>2]|0;c[l+0>>2]=c[r+0>>2];if(rj(e,l,g,p,u,j,q,m,n,o,t)|0){if((a[k]&1)==0){c[k+4>>2]=0;a[k]=0}else{c[c[k+8>>2]>>2]=0;c[k+4>>2]=0}if((a[q]|0)!=0){Oe(k,vc[c[(c[m>>2]|0)+44>>2]&15](m,45)|0)}q=vc[c[(c[m>>2]|0)+44>>2]&15](m,48)|0;m=c[n>>2]|0;t=c[o>>2]|0;o=t+ -4|0;a:do{if(m>>>0<o>>>0){u=m;while(1){g=u+4|0;if((c[u>>2]|0)!=(q|0)){x=u;break a}if(g>>>0<o>>>0){u=g}else{x=g;break}}}else{x=m}}while(0);tj(k,x,t)|0}t=c[e>>2]|0;do{if((t|0)!=0){x=c[t+12>>2]|0;if((x|0)==(c[t+16>>2]|0)){y=mc[c[(c[t>>2]|0)+36>>2]&63](t)|0}else{y=c[x>>2]|0}if((y|0)==-1){c[e>>2]=0;z=1;break}else{z=(c[e>>2]|0)==0;break}}else{z=1}}while(0);do{if((w|0)!=0){y=c[w+12>>2]|0;if((y|0)==(c[w+16>>2]|0)){A=mc[c[(c[w>>2]|0)+36>>2]&63](w)|0}else{A=c[y>>2]|0}if(!((A|0)==-1)){if(z){break}else{B=31;break}}else{c[f>>2]=0;B=29;break}}else{B=29}}while(0);if((B|0)==29?z:0){B=31}if((B|0)==31){c[j>>2]=c[j>>2]|2}c[b>>2]=c[e>>2];ae(c[p>>2]|0)|0;p=c[n>>2]|0;c[n>>2]=0;if((p|0)==0){i=d;return}jc[c[s>>2]&127](p);i=d;return}function tj(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0;f=i;g=d;h=a[b]|0;if((h&1)==0){j=1;k=(h&255)>>>1;l=h}else{h=c[b>>2]|0;j=(h&-2)+ -1|0;k=c[b+4>>2]|0;l=h&255}h=e-g>>2;if((h|0)==0){i=f;return b|0}if((j-k|0)>>>0<h>>>0){Qe(b,j,k+h-j|0,k,k,0,0);m=a[b]|0}else{m=l}if((m&1)==0){n=b+4|0}else{n=c[b+8>>2]|0}m=n+(k<<2)|0;if((d|0)==(e|0)){o=m}else{l=k+((e+ -4+(0-g)|0)>>>2)+1|0;g=d;d=m;while(1){c[d>>2]=c[g>>2];g=g+4|0;if((g|0)==(e|0)){break}else{d=d+4|0}}o=n+(l<<2)|0}c[o>>2]=0;o=k+h|0;if((a[b]&1)==0){a[b]=o<<1;i=f;return b|0}else{c[b+4>>2]=o;i=f;return b|0}return 0}function uj(b,d,e,f,g,h,j,k,l,m){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;m=m|0;var n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0;n=i;i=i+128|0;o=n;p=n+100|0;q=n+88|0;r=n+76|0;s=n+64|0;t=n+104|0;u=n+36|0;v=n+24|0;w=n+12|0;x=n+40|0;y=n+52|0;if(b){b=c[d>>2]|0;if(!((c[1956]|0)==-1)){c[o>>2]=7824;c[o+4>>2]=112;c[o+8>>2]=0;ue(7824,o,113)}z=(c[7828>>2]|0)+ -1|0;A=c[b+8>>2]|0;if(!((c[b+12>>2]|0)-A>>2>>>0>z>>>0)){B=Wa(4)|0;qm(B);xb(B|0,16216,101)}b=c[A+(z<<2)>>2]|0;if((b|0)==0){B=Wa(4)|0;qm(B);xb(B|0,16216,101)}kc[c[(c[b>>2]|0)+44>>2]&63](p,b);B=c[p>>2]|0;a[e]=B;a[e+1|0]=B>>8;a[e+2|0]=B>>16;a[e+3|0]=B>>24;kc[c[(c[b>>2]|0)+32>>2]&63](q,b);if((a[l]&1)==0){c[l+4>>2]=0;a[l]=0}else{c[c[l+8>>2]>>2]=0;c[l+4>>2]=0}Ne(l,0);c[l+0>>2]=c[q+0>>2];c[l+4>>2]=c[q+4>>2];c[l+8>>2]=c[q+8>>2];c[q+0>>2]=0;c[q+4>>2]=0;c[q+8>>2]=0;Ke(q);kc[c[(c[b>>2]|0)+28>>2]&63](r,b);if((a[k]&1)==0){c[k+4>>2]=0;a[k]=0}else{c[c[k+8>>2]>>2]=0;c[k+4>>2]=0}Ne(k,0);c[k+0>>2]=c[r+0>>2];c[k+4>>2]=c[r+4>>2];c[k+8>>2]=c[r+8>>2];c[r+0>>2]=0;c[r+4>>2]=0;c[r+8>>2]=0;Ke(r);c[f>>2]=mc[c[(c[b>>2]|0)+12>>2]&63](b)|0;c[g>>2]=mc[c[(c[b>>2]|0)+16>>2]&63](b)|0;kc[c[(c[b>>2]|0)+20>>2]&63](s,b);if((a[h]&1)==0){a[h+1|0]=0;a[h]=0}else{a[c[h+8>>2]|0]=0;c[h+4>>2]=0}De(h,0);c[h+0>>2]=c[s+0>>2];c[h+4>>2]=c[s+4>>2];c[h+8>>2]=c[s+8>>2];c[s+0>>2]=0;c[s+4>>2]=0;c[s+8>>2]=0;ze(s);kc[c[(c[b>>2]|0)+24>>2]&63](t,b);if((a[j]&1)==0){c[j+4>>2]=0;a[j]=0}else{c[c[j+8>>2]>>2]=0;c[j+4>>2]=0}Ne(j,0);c[j+0>>2]=c[t+0>>2];c[j+4>>2]=c[t+4>>2];c[j+8>>2]=c[t+8>>2];c[t+0>>2]=0;c[t+4>>2]=0;c[t+8>>2]=0;Ke(t);C=mc[c[(c[b>>2]|0)+36>>2]&63](b)|0;c[m>>2]=C;i=n;return}else{b=c[d>>2]|0;if(!((c[1940]|0)==-1)){c[o>>2]=7760;c[o+4>>2]=112;c[o+8>>2]=0;ue(7760,o,113)}o=(c[7764>>2]|0)+ -1|0;d=c[b+8>>2]|0;if(!((c[b+12>>2]|0)-d>>2>>>0>o>>>0)){D=Wa(4)|0;qm(D);xb(D|0,16216,101)}b=c[d+(o<<2)>>2]|0;if((b|0)==0){D=Wa(4)|0;qm(D);xb(D|0,16216,101)}kc[c[(c[b>>2]|0)+44>>2]&63](u,b);D=c[u>>2]|0;a[e]=D;a[e+1|0]=D>>8;a[e+2|0]=D>>16;a[e+3|0]=D>>24;kc[c[(c[b>>2]|0)+32>>2]&63](v,b);if((a[l]&1)==0){c[l+4>>2]=0;a[l]=0}else{c[c[l+8>>2]>>2]=0;c[l+4>>2]=0}Ne(l,0);c[l+0>>2]=c[v+0>>2];c[l+4>>2]=c[v+4>>2];c[l+8>>2]=c[v+8>>2];c[v+0>>2]=0;c[v+4>>2]=0;c[v+8>>2]=0;Ke(v);kc[c[(c[b>>2]|0)+28>>2]&63](w,b);if((a[k]&1)==0){c[k+4>>2]=0;a[k]=0}else{c[c[k+8>>2]>>2]=0;c[k+4>>2]=0}Ne(k,0);c[k+0>>2]=c[w+0>>2];c[k+4>>2]=c[w+4>>2];c[k+8>>2]=c[w+8>>2];c[w+0>>2]=0;c[w+4>>2]=0;c[w+8>>2]=0;Ke(w);c[f>>2]=mc[c[(c[b>>2]|0)+12>>2]&63](b)|0;c[g>>2]=mc[c[(c[b>>2]|0)+16>>2]&63](b)|0;kc[c[(c[b>>2]|0)+20>>2]&63](x,b);if((a[h]&1)==0){a[h+1|0]=0;a[h]=0}else{a[c[h+8>>2]|0]=0;c[h+4>>2]=0}De(h,0);c[h+0>>2]=c[x+0>>2];c[h+4>>2]=c[x+4>>2];c[h+8>>2]=c[x+8>>2];c[x+0>>2]=0;c[x+4>>2]=0;c[x+8>>2]=0;ze(x);kc[c[(c[b>>2]|0)+24>>2]&63](y,b);if((a[j]&1)==0){c[j+4>>2]=0;a[j]=0}else{c[c[j+8>>2]>>2]=0;c[j+4>>2]=0}Ne(j,0);c[j+0>>2]=c[y+0>>2];c[j+4>>2]=c[y+4>>2];c[j+8>>2]=c[y+8>>2];c[y+0>>2]=0;c[y+4>>2]=0;c[y+8>>2]=0;Ke(y);C=mc[c[(c[b>>2]|0)+36>>2]&63](b)|0;c[m>>2]=C;i=n;return}}function vj(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function wj(a){a=a|0;return}function xj(b,d,e,f,g,j,l){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;j=j|0;l=+l;var m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0;d=i;i=i+384|0;m=d;n=d+276|0;o=d+52|0;p=d+176|0;q=d+60|0;r=d+56|0;s=d+376|0;t=d+377|0;u=d+64|0;v=d+32|0;w=d+20|0;x=d+16|0;y=d+76|0;z=d+12|0;A=d+44|0;B=d+48|0;c[o>>2]=n;h[k>>3]=l;c[m>>2]=c[k>>2];c[m+4>>2]=c[k+4>>2];C=mb(n|0,100,8008,m|0)|0;if(C>>>0>99){if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}n=c[2036]|0;h[k>>3]=l;c[m>>2]=c[k>>2];c[m+4>>2]=c[k+4>>2];D=sh(o,n,8008,m)|0;n=c[o>>2]|0;if((n|0)==0){_m()}E=Lm(D)|0;if((E|0)==0){_m()}else{F=E;G=n;H=E;I=D}}else{F=0;G=0;H=p;I=C}Te(q,g);C=c[q>>2]|0;if(!((c[2064]|0)==-1)){c[m>>2]=8256;c[m+4>>2]=112;c[m+8>>2]=0;ue(8256,m,113)}p=(c[8260>>2]|0)+ -1|0;D=c[C+8>>2]|0;if(!((c[C+12>>2]|0)-D>>2>>>0>p>>>0)){J=Wa(4)|0;qm(J);xb(J|0,16216,101)}C=c[D+(p<<2)>>2]|0;if((C|0)==0){J=Wa(4)|0;qm(J);xb(J|0,16216,101)}J=c[o>>2]|0;sc[c[(c[C>>2]|0)+32>>2]&7](C,J,J+I|0,H)|0;if((I|0)==0){K=0}else{K=(a[c[o>>2]|0]|0)==45}c[r>>2]=0;c[u+0>>2]=0;c[u+4>>2]=0;c[u+8>>2]=0;c[v+0>>2]=0;c[v+4>>2]=0;c[v+8>>2]=0;c[w+0>>2]=0;c[w+4>>2]=0;c[w+8>>2]=0;yj(f,K,q,r,s,t,u,v,w,x);f=c[x>>2]|0;if((I|0)>(f|0)){x=a[w]|0;if((x&1)==0){L=(x&255)>>>1}else{L=c[w+4>>2]|0}x=a[v]|0;if((x&1)==0){M=(x&255)>>>1}else{M=c[v+4>>2]|0}N=L+(I-f<<1|1)+M|0}else{M=a[w]|0;if((M&1)==0){O=(M&255)>>>1}else{O=c[w+4>>2]|0}M=a[v]|0;if((M&1)==0){P=(M&255)>>>1}else{P=c[v+4>>2]|0}N=O+2+P|0}P=N+f|0;if(P>>>0>100){N=Lm(P)|0;if((N|0)==0){_m()}else{Q=N;R=N}}else{Q=0;R=y}zj(R,z,A,c[g+4>>2]|0,H,H+I|0,C,K,r,a[s]|0,a[t]|0,u,v,w,f);c[B>>2]=c[e>>2];e=c[z>>2]|0;z=c[A>>2]|0;c[m+0>>2]=c[B+0>>2];nh(b,m,R,e,z,g,j);if((Q|0)!=0){Mm(Q)}ze(w);ze(v);ze(u);ae(c[q>>2]|0)|0;if((F|0)!=0){Mm(F)}if((G|0)==0){i=d;return}Mm(G);i=d;return}function yj(b,d,e,f,g,h,j,k,l,m){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;m=m|0;var n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0;n=i;i=i+128|0;o=n;p=n+108|0;q=n+96|0;r=n+92|0;s=n+80|0;t=n+68|0;u=n+112|0;v=n+40|0;w=n+28|0;x=n+24|0;y=n+12|0;z=n+44|0;A=n+56|0;B=c[e>>2]|0;if(b){if(!((c[1924]|0)==-1)){c[o>>2]=7696;c[o+4>>2]=112;c[o+8>>2]=0;ue(7696,o,113)}b=(c[7700>>2]|0)+ -1|0;e=c[B+8>>2]|0;if(!((c[B+12>>2]|0)-e>>2>>>0>b>>>0)){C=Wa(4)|0;qm(C);xb(C|0,16216,101)}D=c[e+(b<<2)>>2]|0;if((D|0)==0){C=Wa(4)|0;qm(C);xb(C|0,16216,101)}C=c[D>>2]|0;if(d){kc[c[C+44>>2]&63](p,D);b=c[p>>2]|0;a[f]=b;a[f+1|0]=b>>8;a[f+2|0]=b>>16;a[f+3|0]=b>>24;kc[c[(c[D>>2]|0)+32>>2]&63](q,D);if((a[l]&1)==0){a[l+1|0]=0;a[l]=0}else{a[c[l+8>>2]|0]=0;c[l+4>>2]=0}De(l,0);c[l+0>>2]=c[q+0>>2];c[l+4>>2]=c[q+4>>2];c[l+8>>2]=c[q+8>>2];c[q+0>>2]=0;c[q+4>>2]=0;c[q+8>>2]=0;ze(q)}else{kc[c[C+40>>2]&63](r,D);C=c[r>>2]|0;a[f]=C;a[f+1|0]=C>>8;a[f+2|0]=C>>16;a[f+3|0]=C>>24;kc[c[(c[D>>2]|0)+28>>2]&63](s,D);if((a[l]&1)==0){a[l+1|0]=0;a[l]=0}else{a[c[l+8>>2]|0]=0;c[l+4>>2]=0}De(l,0);c[l+0>>2]=c[s+0>>2];c[l+4>>2]=c[s+4>>2];c[l+8>>2]=c[s+8>>2];c[s+0>>2]=0;c[s+4>>2]=0;c[s+8>>2]=0;ze(s)}a[g]=mc[c[(c[D>>2]|0)+12>>2]&63](D)|0;a[h]=mc[c[(c[D>>2]|0)+16>>2]&63](D)|0;kc[c[(c[D>>2]|0)+20>>2]&63](t,D);if((a[j]&1)==0){a[j+1|0]=0;a[j]=0}else{a[c[j+8>>2]|0]=0;c[j+4>>2]=0}De(j,0);c[j+0>>2]=c[t+0>>2];c[j+4>>2]=c[t+4>>2];c[j+8>>2]=c[t+8>>2];c[t+0>>2]=0;c[t+4>>2]=0;c[t+8>>2]=0;ze(t);kc[c[(c[D>>2]|0)+24>>2]&63](u,D);if((a[k]&1)==0){a[k+1|0]=0;a[k]=0}else{a[c[k+8>>2]|0]=0;c[k+4>>2]=0}De(k,0);c[k+0>>2]=c[u+0>>2];c[k+4>>2]=c[u+4>>2];c[k+8>>2]=c[u+8>>2];c[u+0>>2]=0;c[u+4>>2]=0;c[u+8>>2]=0;ze(u);E=mc[c[(c[D>>2]|0)+36>>2]&63](D)|0;c[m>>2]=E;i=n;return}else{if(!((c[1908]|0)==-1)){c[o>>2]=7632;c[o+4>>2]=112;c[o+8>>2]=0;ue(7632,o,113)}o=(c[7636>>2]|0)+ -1|0;D=c[B+8>>2]|0;if(!((c[B+12>>2]|0)-D>>2>>>0>o>>>0)){F=Wa(4)|0;qm(F);xb(F|0,16216,101)}B=c[D+(o<<2)>>2]|0;if((B|0)==0){F=Wa(4)|0;qm(F);xb(F|0,16216,101)}F=c[B>>2]|0;if(d){kc[c[F+44>>2]&63](v,B);d=c[v>>2]|0;a[f]=d;a[f+1|0]=d>>8;a[f+2|0]=d>>16;a[f+3|0]=d>>24;kc[c[(c[B>>2]|0)+32>>2]&63](w,B);if((a[l]&1)==0){a[l+1|0]=0;a[l]=0}else{a[c[l+8>>2]|0]=0;c[l+4>>2]=0}De(l,0);c[l+0>>2]=c[w+0>>2];c[l+4>>2]=c[w+4>>2];c[l+8>>2]=c[w+8>>2];c[w+0>>2]=0;c[w+4>>2]=0;c[w+8>>2]=0;ze(w)}else{kc[c[F+40>>2]&63](x,B);F=c[x>>2]|0;a[f]=F;a[f+1|0]=F>>8;a[f+2|0]=F>>16;a[f+3|0]=F>>24;kc[c[(c[B>>2]|0)+28>>2]&63](y,B);if((a[l]&1)==0){a[l+1|0]=0;a[l]=0}else{a[c[l+8>>2]|0]=0;c[l+4>>2]=0}De(l,0);c[l+0>>2]=c[y+0>>2];c[l+4>>2]=c[y+4>>2];c[l+8>>2]=c[y+8>>2];c[y+0>>2]=0;c[y+4>>2]=0;c[y+8>>2]=0;ze(y)}a[g]=mc[c[(c[B>>2]|0)+12>>2]&63](B)|0;a[h]=mc[c[(c[B>>2]|0)+16>>2]&63](B)|0;kc[c[(c[B>>2]|0)+20>>2]&63](z,B);if((a[j]&1)==0){a[j+1|0]=0;a[j]=0}else{a[c[j+8>>2]|0]=0;c[j+4>>2]=0}De(j,0);c[j+0>>2]=c[z+0>>2];c[j+4>>2]=c[z+4>>2];c[j+8>>2]=c[z+8>>2];c[z+0>>2]=0;c[z+4>>2]=0;c[z+8>>2]=0;ze(z);kc[c[(c[B>>2]|0)+24>>2]&63](A,B);if((a[k]&1)==0){a[k+1|0]=0;a[k]=0}else{a[c[k+8>>2]|0]=0;c[k+4>>2]=0}De(k,0);c[k+0>>2]=c[A+0>>2];c[k+4>>2]=c[A+4>>2];c[k+8>>2]=c[A+8>>2];c[A+0>>2]=0;c[A+4>>2]=0;c[A+8>>2]=0;ze(A);E=mc[c[(c[B>>2]|0)+36>>2]&63](B)|0;c[m>>2]=E;i=n;return}}function zj(d,e,f,g,h,j,k,l,m,n,o,p,q,r,s){d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;m=m|0;n=n|0;o=o|0;p=p|0;q=q|0;r=r|0;s=s|0;var t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0,X=0,Y=0,Z=0,_=0,$=0,aa=0,ba=0,ca=0,da=0,ea=0,fa=0,ga=0,ha=0,ia=0,ja=0,ka=0,la=0,ma=0,na=0,oa=0,pa=0,qa=0,ra=0,sa=0,ta=0,ua=0,va=0;t=i;c[f>>2]=d;u=r+1|0;v=r+8|0;w=r+4|0;x=(g&512|0)==0;y=q+1|0;z=q+8|0;A=q+4|0;B=(s|0)>0;C=p+1|0;D=p+8|0;E=p+4|0;F=k+8|0;G=0-s|0;H=h;h=0;while(1){switch(a[m+h|0]|0){case 0:{c[e>>2]=c[f>>2];I=H;break};case 4:{J=c[f>>2]|0;K=l?H+1|0:H;a:do{if(K>>>0<j>>>0){L=K;while(1){M=a[L]|0;if(!(M<<24>>24>-1)){N=L;break a}O=L+1|0;if((b[(c[F>>2]|0)+(M<<24>>24<<1)>>1]&2048)==0){N=L;break a}if(O>>>0<j>>>0){L=O}else{N=O;break}}}else{N=K}}while(0);L=N;if(B){if(N>>>0>K>>>0){O=K+(0-L)|0;L=O>>>0<G>>>0?G:O;O=L+s|0;M=J;P=N;Q=s;while(1){R=P+ -1|0;S=a[R]|0;c[f>>2]=M+1;a[M]=S;S=Q+ -1|0;T=(S|0)>0;if(!(R>>>0>K>>>0&T)){break}M=c[f>>2]|0;P=R;Q=S}Q=N+L|0;if(T){U=Q;V=O;W=32}else{X=0;Y=Q;Z=O}}else{U=N;V=s;W=32}if((W|0)==32){W=0;X=vc[c[(c[k>>2]|0)+28>>2]&15](k,48)|0;Y=U;Z=V}Q=c[f>>2]|0;c[f>>2]=Q+1;if((Z|0)>0){P=Q;M=Z;while(1){a[P]=X;S=M+ -1|0;R=c[f>>2]|0;c[f>>2]=R+1;if((S|0)>0){P=R;M=S}else{_=R;break}}}else{_=Q}a[_]=n;$=Y}else{$=N}if(($|0)==(K|0)){M=vc[c[(c[k>>2]|0)+28>>2]&15](k,48)|0;P=c[f>>2]|0;c[f>>2]=P+1;a[P]=M}else{M=a[p]|0;P=(M&1)==0;if(P){aa=(M&255)>>>1}else{aa=c[E>>2]|0}if((aa|0)==0){ba=$;ca=-1;da=0;ea=0}else{if(P){fa=C}else{fa=c[D>>2]|0}ba=$;ca=a[fa]|0;da=0;ea=0}while(1){if((ea|0)==(ca|0)){P=c[f>>2]|0;c[f>>2]=P+1;a[P]=o;P=da+1|0;M=a[p]|0;O=(M&1)==0;if(O){ga=(M&255)>>>1}else{ga=c[E>>2]|0}if(P>>>0<ga>>>0){if(O){ha=C}else{ha=c[D>>2]|0}if((a[ha+P|0]|0)==127){ia=-1;ja=P;ka=0}else{if(O){la=C}else{la=c[D>>2]|0}ia=a[la+P|0]|0;ja=P;ka=0}}else{ia=ca;ja=P;ka=0}}else{ia=ca;ja=da;ka=ea}ba=ba+ -1|0;P=a[ba]|0;O=c[f>>2]|0;c[f>>2]=O+1;a[O]=P;if((ba|0)==(K|0)){break}else{ca=ia;da=ja;ea=ka+1|0}}}Q=c[f>>2]|0;if((J|0)!=(Q|0)?(P=Q+ -1|0,P>>>0>J>>>0):0){Q=J;O=P;while(1){P=a[Q]|0;a[Q]=a[O]|0;a[O]=P;P=Q+1|0;M=O+ -1|0;if(P>>>0<M>>>0){Q=P;O=M}else{I=K;break}}}else{I=K}break};case 2:{O=a[q]|0;Q=(O&1)==0;if(Q){ma=(O&255)>>>1}else{ma=c[A>>2]|0}if((ma|0)==0|x){I=H}else{if(Q){na=y;oa=(O&255)>>>1}else{na=c[z>>2]|0;oa=c[A>>2]|0}O=na+oa|0;Q=c[f>>2]|0;if((na|0)==(O|0)){pa=Q}else{J=Q;Q=na;while(1){a[J]=a[Q]|0;M=Q+1|0;P=J+1|0;if((M|0)==(O|0)){pa=P;break}else{J=P;Q=M}}}c[f>>2]=pa;I=H}break};case 1:{c[e>>2]=c[f>>2];Q=vc[c[(c[k>>2]|0)+28>>2]&15](k,32)|0;J=c[f>>2]|0;c[f>>2]=J+1;a[J]=Q;I=H;break};case 3:{Q=a[r]|0;J=(Q&1)==0;if(J){qa=(Q&255)>>>1}else{qa=c[w>>2]|0}if((qa|0)==0){I=H}else{if(J){ra=u}else{ra=c[v>>2]|0}J=a[ra]|0;Q=c[f>>2]|0;c[f>>2]=Q+1;a[Q]=J;I=H}break};default:{I=H}}h=h+1|0;if((h|0)==4){break}else{H=I}}I=a[r]|0;r=(I&1)==0;if(r){sa=(I&255)>>>1}else{sa=c[w>>2]|0}if(sa>>>0>1){if(r){ta=u;ua=(I&255)>>>1}else{ta=c[v>>2]|0;ua=c[w>>2]|0}w=ta+1|0;v=ta+ua|0;ua=c[f>>2]|0;if((w|0)==(v|0)){va=ua}else{ta=ua;ua=w;while(1){a[ta]=a[ua]|0;w=ua+1|0;I=ta+1|0;if((w|0)==(v|0)){va=I;break}else{ta=I;ua=w}}}c[f>>2]=va}va=g&176;if((va|0)==32){c[e>>2]=c[f>>2];i=t;return}else if((va|0)==16){i=t;return}else{c[e>>2]=d;i=t;return}}function Aj(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0;d=i;i=i+176|0;k=d;l=d+52|0;m=d+56|0;n=d+172|0;o=d+173|0;p=d+40|0;q=d+60|0;r=d+20|0;s=d+16|0;t=d+72|0;u=d+12|0;v=d+32|0;w=d+36|0;Te(l,g);x=c[l>>2]|0;if(!((c[2064]|0)==-1)){c[k>>2]=8256;c[k+4>>2]=112;c[k+8>>2]=0;ue(8256,k,113)}y=(c[8260>>2]|0)+ -1|0;z=c[x+8>>2]|0;if(!((c[x+12>>2]|0)-z>>2>>>0>y>>>0)){A=Wa(4)|0;qm(A);xb(A|0,16216,101)}x=c[z+(y<<2)>>2]|0;if((x|0)==0){A=Wa(4)|0;qm(A);xb(A|0,16216,101)}A=a[j]|0;y=(A&1)==0;if(y){B=(A&255)>>>1}else{B=c[j+4>>2]|0}if((B|0)==0){C=0}else{if(y){D=j+1|0}else{D=c[j+8>>2]|0}y=a[D]|0;C=y<<24>>24==(vc[c[(c[x>>2]|0)+28>>2]&15](x,45)|0)<<24>>24}c[m>>2]=0;c[p+0>>2]=0;c[p+4>>2]=0;c[p+8>>2]=0;c[q+0>>2]=0;c[q+4>>2]=0;c[q+8>>2]=0;c[r+0>>2]=0;c[r+4>>2]=0;c[r+8>>2]=0;yj(f,C,l,m,n,o,p,q,r,s);f=a[j]|0;y=(f&1)==0;if(y){E=(f&255)>>>1}else{E=c[j+4>>2]|0}D=c[s>>2]|0;if((E|0)>(D|0)){if(y){F=(f&255)>>>1}else{F=c[j+4>>2]|0}y=a[r]|0;if((y&1)==0){G=(y&255)>>>1}else{G=c[r+4>>2]|0}y=a[q]|0;if((y&1)==0){H=(y&255)>>>1}else{H=c[q+4>>2]|0}I=G+(F-D<<1|1)+H|0}else{H=a[r]|0;if((H&1)==0){J=(H&255)>>>1}else{J=c[r+4>>2]|0}H=a[q]|0;if((H&1)==0){K=(H&255)>>>1}else{K=c[q+4>>2]|0}I=J+2+K|0}K=I+D|0;if(K>>>0>100){I=Lm(K)|0;if((I|0)==0){_m()}else{L=I;M=I}}else{L=0;M=t}if((f&1)==0){N=j+1|0;O=(f&255)>>>1}else{N=c[j+8>>2]|0;O=c[j+4>>2]|0}zj(M,u,v,c[g+4>>2]|0,N,N+O|0,x,C,m,a[n]|0,a[o]|0,p,q,r,D);c[w>>2]=c[e>>2];e=c[u>>2]|0;u=c[v>>2]|0;c[k+0>>2]=c[w+0>>2];nh(b,k,M,e,u,g,h);if((L|0)==0){ze(r);ze(q);ze(p);P=c[l>>2]|0;ae(P)|0;i=d;return}Mm(L);ze(r);ze(q);ze(p);P=c[l>>2]|0;ae(P)|0;i=d;return}function Bj(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function Cj(a){a=a|0;return}function Dj(b,d,e,f,g,j,l){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;j=j|0;l=+l;var m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0;d=i;i=i+992|0;m=d;n=d+888|0;o=d+872|0;p=d+472|0;q=d+464|0;r=d+460|0;s=d+456|0;t=d+452|0;u=d+876|0;v=d+432|0;w=d+420|0;x=d+416|0;y=d+16|0;z=d+12|0;A=d+444|0;B=d+448|0;c[o>>2]=n;h[k>>3]=l;c[m>>2]=c[k>>2];c[m+4>>2]=c[k+4>>2];C=mb(n|0,100,8008,m|0)|0;if(C>>>0>99){if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}n=c[2036]|0;h[k>>3]=l;c[m>>2]=c[k>>2];c[m+4>>2]=c[k+4>>2];D=sh(o,n,8008,m)|0;n=c[o>>2]|0;if((n|0)==0){_m()}E=Lm(D<<2)|0;if((E|0)==0){_m()}else{F=E;G=n;H=E;I=D}}else{F=0;G=0;H=p;I=C}Te(q,g);C=c[q>>2]|0;if(!((c[2062]|0)==-1)){c[m>>2]=8248;c[m+4>>2]=112;c[m+8>>2]=0;ue(8248,m,113)}p=(c[8252>>2]|0)+ -1|0;D=c[C+8>>2]|0;if(!((c[C+12>>2]|0)-D>>2>>>0>p>>>0)){J=Wa(4)|0;qm(J);xb(J|0,16216,101)}C=c[D+(p<<2)>>2]|0;if((C|0)==0){J=Wa(4)|0;qm(J);xb(J|0,16216,101)}J=c[o>>2]|0;sc[c[(c[C>>2]|0)+48>>2]&7](C,J,J+I|0,H)|0;if((I|0)==0){K=0}else{K=(a[c[o>>2]|0]|0)==45}c[r>>2]=0;c[u+0>>2]=0;c[u+4>>2]=0;c[u+8>>2]=0;c[v+0>>2]=0;c[v+4>>2]=0;c[v+8>>2]=0;c[w+0>>2]=0;c[w+4>>2]=0;c[w+8>>2]=0;Ej(f,K,q,r,s,t,u,v,w,x);f=c[x>>2]|0;if((I|0)>(f|0)){x=a[w]|0;if((x&1)==0){L=(x&255)>>>1}else{L=c[w+4>>2]|0}x=a[v]|0;if((x&1)==0){M=(x&255)>>>1}else{M=c[v+4>>2]|0}N=L+(I-f<<1|1)+M|0}else{M=a[w]|0;if((M&1)==0){O=(M&255)>>>1}else{O=c[w+4>>2]|0}M=a[v]|0;if((M&1)==0){P=(M&255)>>>1}else{P=c[v+4>>2]|0}N=O+2+P|0}P=N+f|0;if(P>>>0>100){N=Lm(P<<2)|0;if((N|0)==0){_m()}else{Q=N;R=N}}else{Q=0;R=y}Fj(R,z,A,c[g+4>>2]|0,H,H+(I<<2)|0,C,K,r,c[s>>2]|0,c[t>>2]|0,u,v,w,f);c[B>>2]=c[e>>2];e=c[z>>2]|0;z=c[A>>2]|0;c[m+0>>2]=c[B+0>>2];Bh(b,m,R,e,z,g,j);if((Q|0)!=0){Mm(Q)}Ke(w);Ke(v);ze(u);ae(c[q>>2]|0)|0;if((F|0)!=0){Mm(F)}if((G|0)==0){i=d;return}Mm(G);i=d;return}function Ej(b,d,e,f,g,h,j,k,l,m){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;m=m|0;var n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0;n=i;i=i+128|0;o=n;p=n+108|0;q=n+96|0;r=n+92|0;s=n+80|0;t=n+68|0;u=n+112|0;v=n+40|0;w=n+28|0;x=n+24|0;y=n+12|0;z=n+44|0;A=n+56|0;B=c[e>>2]|0;if(b){if(!((c[1956]|0)==-1)){c[o>>2]=7824;c[o+4>>2]=112;c[o+8>>2]=0;ue(7824,o,113)}b=(c[7828>>2]|0)+ -1|0;e=c[B+8>>2]|0;if(!((c[B+12>>2]|0)-e>>2>>>0>b>>>0)){C=Wa(4)|0;qm(C);xb(C|0,16216,101)}D=c[e+(b<<2)>>2]|0;if((D|0)==0){C=Wa(4)|0;qm(C);xb(C|0,16216,101)}C=c[D>>2]|0;if(d){kc[c[C+44>>2]&63](p,D);b=c[p>>2]|0;a[f]=b;a[f+1|0]=b>>8;a[f+2|0]=b>>16;a[f+3|0]=b>>24;kc[c[(c[D>>2]|0)+32>>2]&63](q,D);if((a[l]&1)==0){c[l+4>>2]=0;a[l]=0}else{c[c[l+8>>2]>>2]=0;c[l+4>>2]=0}Ne(l,0);c[l+0>>2]=c[q+0>>2];c[l+4>>2]=c[q+4>>2];c[l+8>>2]=c[q+8>>2];c[q+0>>2]=0;c[q+4>>2]=0;c[q+8>>2]=0;Ke(q)}else{kc[c[C+40>>2]&63](r,D);C=c[r>>2]|0;a[f]=C;a[f+1|0]=C>>8;a[f+2|0]=C>>16;a[f+3|0]=C>>24;kc[c[(c[D>>2]|0)+28>>2]&63](s,D);if((a[l]&1)==0){c[l+4>>2]=0;a[l]=0}else{c[c[l+8>>2]>>2]=0;c[l+4>>2]=0}Ne(l,0);c[l+0>>2]=c[s+0>>2];c[l+4>>2]=c[s+4>>2];c[l+8>>2]=c[s+8>>2];c[s+0>>2]=0;c[s+4>>2]=0;c[s+8>>2]=0;Ke(s)}c[g>>2]=mc[c[(c[D>>2]|0)+12>>2]&63](D)|0;c[h>>2]=mc[c[(c[D>>2]|0)+16>>2]&63](D)|0;kc[c[(c[D>>2]|0)+20>>2]&63](t,D);if((a[j]&1)==0){a[j+1|0]=0;a[j]=0}else{a[c[j+8>>2]|0]=0;c[j+4>>2]=0}De(j,0);c[j+0>>2]=c[t+0>>2];c[j+4>>2]=c[t+4>>2];c[j+8>>2]=c[t+8>>2];c[t+0>>2]=0;c[t+4>>2]=0;c[t+8>>2]=0;ze(t);kc[c[(c[D>>2]|0)+24>>2]&63](u,D);if((a[k]&1)==0){c[k+4>>2]=0;a[k]=0}else{c[c[k+8>>2]>>2]=0;c[k+4>>2]=0}Ne(k,0);c[k+0>>2]=c[u+0>>2];c[k+4>>2]=c[u+4>>2];c[k+8>>2]=c[u+8>>2];c[u+0>>2]=0;c[u+4>>2]=0;c[u+8>>2]=0;Ke(u);E=mc[c[(c[D>>2]|0)+36>>2]&63](D)|0;c[m>>2]=E;i=n;return}else{if(!((c[1940]|0)==-1)){c[o>>2]=7760;c[o+4>>2]=112;c[o+8>>2]=0;ue(7760,o,113)}o=(c[7764>>2]|0)+ -1|0;D=c[B+8>>2]|0;if(!((c[B+12>>2]|0)-D>>2>>>0>o>>>0)){F=Wa(4)|0;qm(F);xb(F|0,16216,101)}B=c[D+(o<<2)>>2]|0;if((B|0)==0){F=Wa(4)|0;qm(F);xb(F|0,16216,101)}F=c[B>>2]|0;if(d){kc[c[F+44>>2]&63](v,B);d=c[v>>2]|0;a[f]=d;a[f+1|0]=d>>8;a[f+2|0]=d>>16;a[f+3|0]=d>>24;kc[c[(c[B>>2]|0)+32>>2]&63](w,B);if((a[l]&1)==0){c[l+4>>2]=0;a[l]=0}else{c[c[l+8>>2]>>2]=0;c[l+4>>2]=0}Ne(l,0);c[l+0>>2]=c[w+0>>2];c[l+4>>2]=c[w+4>>2];c[l+8>>2]=c[w+8>>2];c[w+0>>2]=0;c[w+4>>2]=0;c[w+8>>2]=0;Ke(w)}else{kc[c[F+40>>2]&63](x,B);F=c[x>>2]|0;a[f]=F;a[f+1|0]=F>>8;a[f+2|0]=F>>16;a[f+3|0]=F>>24;kc[c[(c[B>>2]|0)+28>>2]&63](y,B);if((a[l]&1)==0){c[l+4>>2]=0;a[l]=0}else{c[c[l+8>>2]>>2]=0;c[l+4>>2]=0}Ne(l,0);c[l+0>>2]=c[y+0>>2];c[l+4>>2]=c[y+4>>2];c[l+8>>2]=c[y+8>>2];c[y+0>>2]=0;c[y+4>>2]=0;c[y+8>>2]=0;Ke(y)}c[g>>2]=mc[c[(c[B>>2]|0)+12>>2]&63](B)|0;c[h>>2]=mc[c[(c[B>>2]|0)+16>>2]&63](B)|0;kc[c[(c[B>>2]|0)+20>>2]&63](z,B);if((a[j]&1)==0){a[j+1|0]=0;a[j]=0}else{a[c[j+8>>2]|0]=0;c[j+4>>2]=0}De(j,0);c[j+0>>2]=c[z+0>>2];c[j+4>>2]=c[z+4>>2];c[j+8>>2]=c[z+8>>2];c[z+0>>2]=0;c[z+4>>2]=0;c[z+8>>2]=0;ze(z);kc[c[(c[B>>2]|0)+24>>2]&63](A,B);if((a[k]&1)==0){c[k+4>>2]=0;a[k]=0}else{c[c[k+8>>2]>>2]=0;c[k+4>>2]=0}Ne(k,0);c[k+0>>2]=c[A+0>>2];c[k+4>>2]=c[A+4>>2];c[k+8>>2]=c[A+8>>2];c[A+0>>2]=0;c[A+4>>2]=0;c[A+8>>2]=0;Ke(A);E=mc[c[(c[B>>2]|0)+36>>2]&63](B)|0;c[m>>2]=E;i=n;return}}function Fj(b,d,e,f,g,h,j,k,l,m,n,o,p,q,r){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;m=m|0;n=n|0;o=o|0;p=p|0;q=q|0;r=r|0;var s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0,X=0,Y=0,Z=0,_=0,$=0,aa=0,ba=0,ca=0,da=0,ea=0,fa=0,ga=0,ha=0,ia=0,ja=0,ka=0,la=0,ma=0,na=0,oa=0,pa=0,qa=0,ra=0,sa=0;s=i;c[e>>2]=b;t=q+4|0;u=q+8|0;v=(f&512|0)==0;w=p+4|0;x=p+8|0;y=(r|0)>0;z=o+1|0;A=o+8|0;B=o+4|0;C=g;g=0;while(1){switch(a[l+g|0]|0){case 1:{c[d>>2]=c[e>>2];D=vc[c[(c[j>>2]|0)+44>>2]&15](j,32)|0;E=c[e>>2]|0;c[e>>2]=E+4;c[E>>2]=D;F=C;break};case 0:{c[d>>2]=c[e>>2];F=C;break};case 4:{D=c[e>>2]|0;E=k?C+4|0:C;a:do{if(E>>>0<h>>>0){G=E;while(1){H=G+4|0;if(!(gc[c[(c[j>>2]|0)+12>>2]&31](j,2048,c[G>>2]|0)|0)){I=G;break a}if(H>>>0<h>>>0){G=H}else{I=H;break}}}else{I=E}}while(0);if(y){if(I>>>0>E>>>0){G=c[e>>2]|0;H=I;J=r;while(1){K=H+ -4|0;L=G+4|0;c[G>>2]=c[K>>2];M=J+ -1|0;N=(M|0)>0;if(K>>>0>E>>>0&N){G=L;H=K;J=M}else{break}}c[e>>2]=L;if(N){O=K;P=M;Q=34}else{J=c[e>>2]|0;c[e>>2]=J+4;R=J;S=K}}else{O=I;P=r;Q=34}if((Q|0)==34){Q=0;J=vc[c[(c[j>>2]|0)+44>>2]&15](j,48)|0;H=c[e>>2]|0;G=H+4|0;c[e>>2]=G;if((P|0)>0){T=H;U=G;G=P;while(1){c[T>>2]=J;G=G+ -1|0;if((G|0)<=0){break}else{V=U;U=U+4|0;T=V}}c[e>>2]=H+(P+1<<2);R=H+(P<<2)|0;S=O}else{R=H;S=O}}c[R>>2]=m;W=S}else{W=I}if((W|0)==(E|0)){T=vc[c[(c[j>>2]|0)+44>>2]&15](j,48)|0;U=c[e>>2]|0;G=U+4|0;c[e>>2]=G;c[U>>2]=T;X=G}else{G=a[o]|0;T=(G&1)==0;if(T){Y=(G&255)>>>1}else{Y=c[B>>2]|0}if((Y|0)==0){Z=W;_=-1;$=0;aa=0}else{if(T){ba=z}else{ba=c[A>>2]|0}Z=W;_=a[ba]|0;$=0;aa=0}while(1){T=c[e>>2]|0;if((aa|0)==(_|0)){G=T+4|0;c[e>>2]=G;c[T>>2]=n;U=$+1|0;J=a[o]|0;V=(J&1)==0;if(V){ca=(J&255)>>>1}else{ca=c[B>>2]|0}if(U>>>0<ca>>>0){if(V){da=z}else{da=c[A>>2]|0}if((a[da+U|0]|0)==127){ea=G;fa=-1;ga=U;ha=0}else{if(V){ia=z}else{ia=c[A>>2]|0}ea=G;fa=a[ia+U|0]|0;ga=U;ha=0}}else{ea=G;fa=_;ga=U;ha=0}}else{ea=T;fa=_;ga=$;ha=aa}T=Z+ -4|0;U=c[T>>2]|0;G=ea+4|0;c[e>>2]=G;c[ea>>2]=U;if((T|0)==(E|0)){X=G;break}else{Z=T;_=fa;$=ga;aa=ha+1|0}}}if((D|0)!=(X|0)?(H=X+ -4|0,H>>>0>D>>>0):0){T=D;G=H;while(1){H=c[T>>2]|0;c[T>>2]=c[G>>2];c[G>>2]=H;H=T+4|0;U=G+ -4|0;if(H>>>0<U>>>0){T=H;G=U}else{F=E;break}}}else{F=E}break};case 2:{G=a[p]|0;T=(G&1)==0;if(T){ja=(G&255)>>>1}else{ja=c[w>>2]|0}if((ja|0)==0|v){F=C}else{if(T){ka=w;la=(G&255)>>>1}else{ka=c[x>>2]|0;la=c[w>>2]|0}G=ka+(la<<2)|0;T=c[e>>2]|0;if((ka|0)==(G|0)){ma=T}else{D=(ka+(la+ -1<<2)+(0-ka)|0)>>>2;U=T;H=ka;while(1){c[U>>2]=c[H>>2];V=H+4|0;if((V|0)==(G|0)){break}U=U+4|0;H=V}ma=T+(D+1<<2)|0}c[e>>2]=ma;F=C}break};case 3:{H=a[q]|0;U=(H&1)==0;if(U){na=(H&255)>>>1}else{na=c[t>>2]|0}if((na|0)==0){F=C}else{if(U){oa=t}else{oa=c[u>>2]|0}U=c[oa>>2]|0;H=c[e>>2]|0;c[e>>2]=H+4;c[H>>2]=U;F=C}break};default:{F=C}}g=g+1|0;if((g|0)==4){break}else{C=F}}F=a[q]|0;q=(F&1)==0;if(q){pa=(F&255)>>>1}else{pa=c[t>>2]|0}if(pa>>>0>1){if(q){qa=t;ra=(F&255)>>>1}else{qa=c[u>>2]|0;ra=c[t>>2]|0}t=qa+4|0;u=qa+(ra<<2)|0;F=c[e>>2]|0;if((t|0)==(u|0)){sa=F}else{q=(qa+(ra+ -1<<2)+(0-t)|0)>>>2;ra=F;qa=t;while(1){c[ra>>2]=c[qa>>2];qa=qa+4|0;if((qa|0)==(u|0)){break}else{ra=ra+4|0}}sa=F+(q+1<<2)|0}c[e>>2]=sa}sa=f&176;if((sa|0)==16){i=s;return}else if((sa|0)==32){c[d>>2]=c[e>>2];i=s;return}else{c[d>>2]=b;i=s;return}}function Gj(b,d,e,f,g,h,j){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0;d=i;i=i+480|0;k=d;l=d+464|0;m=d+460|0;n=d+456|0;o=d+452|0;p=d+440|0;q=d+468|0;r=d+420|0;s=d+416|0;t=d+16|0;u=d+12|0;v=d+432|0;w=d+436|0;Te(l,g);x=c[l>>2]|0;if(!((c[2062]|0)==-1)){c[k>>2]=8248;c[k+4>>2]=112;c[k+8>>2]=0;ue(8248,k,113)}y=(c[8252>>2]|0)+ -1|0;z=c[x+8>>2]|0;if(!((c[x+12>>2]|0)-z>>2>>>0>y>>>0)){A=Wa(4)|0;qm(A);xb(A|0,16216,101)}x=c[z+(y<<2)>>2]|0;if((x|0)==0){A=Wa(4)|0;qm(A);xb(A|0,16216,101)}A=a[j]|0;y=(A&1)==0;if(y){B=(A&255)>>>1}else{B=c[j+4>>2]|0}if((B|0)==0){C=0}else{if(y){D=j+4|0}else{D=c[j+8>>2]|0}y=c[D>>2]|0;C=(y|0)==(vc[c[(c[x>>2]|0)+44>>2]&15](x,45)|0)}c[m>>2]=0;c[p+0>>2]=0;c[p+4>>2]=0;c[p+8>>2]=0;c[q+0>>2]=0;c[q+4>>2]=0;c[q+8>>2]=0;c[r+0>>2]=0;c[r+4>>2]=0;c[r+8>>2]=0;Ej(f,C,l,m,n,o,p,q,r,s);f=a[j]|0;y=(f&1)==0;if(y){E=(f&255)>>>1}else{E=c[j+4>>2]|0}D=c[s>>2]|0;if((E|0)>(D|0)){if(y){F=(f&255)>>>1}else{F=c[j+4>>2]|0}y=a[r]|0;if((y&1)==0){G=(y&255)>>>1}else{G=c[r+4>>2]|0}y=a[q]|0;if((y&1)==0){H=(y&255)>>>1}else{H=c[q+4>>2]|0}I=G+(F-D<<1|1)+H|0}else{H=a[r]|0;if((H&1)==0){J=(H&255)>>>1}else{J=c[r+4>>2]|0}H=a[q]|0;if((H&1)==0){K=(H&255)>>>1}else{K=c[q+4>>2]|0}I=J+2+K|0}K=I+D|0;if(K>>>0>100){I=Lm(K<<2)|0;if((I|0)==0){_m()}else{L=I;M=I}}else{L=0;M=t}if((f&1)==0){N=j+4|0;O=(f&255)>>>1}else{N=c[j+8>>2]|0;O=c[j+4>>2]|0}Fj(M,u,v,c[g+4>>2]|0,N,N+(O<<2)|0,x,C,m,c[n>>2]|0,c[o>>2]|0,p,q,r,D);c[w>>2]=c[e>>2];e=c[u>>2]|0;u=c[v>>2]|0;c[k+0>>2]=c[w+0>>2];Bh(b,k,M,e,u,g,h);if((L|0)==0){Ke(r);Ke(q);ze(p);P=c[l>>2]|0;ae(P)|0;i=d;return}Mm(L);Ke(r);Ke(q);ze(p);P=c[l>>2]|0;ae(P)|0;i=d;return}function Hj(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function Ij(a){a=a|0;return}function Jj(b,d,e){b=b|0;d=d|0;e=e|0;var f=0;e=i;if((a[d]&1)==0){f=d+1|0}else{f=c[d+8>>2]|0}d=bc(f|0,1)|0;i=e;return d>>>((d|0)!=(-1|0)|0)|0}function Kj(b,d,e,f,g,h){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0;d=i;i=i+16|0;j=d;c[j+0>>2]=0;c[j+4>>2]=0;c[j+8>>2]=0;k=a[h]|0;if((k&1)==0){l=h+1|0;m=(k&255)>>>1;n=h+1|0}else{k=c[h+8>>2]|0;l=k;m=c[h+4>>2]|0;n=k}k=l+m|0;if(n>>>0<k>>>0){m=n;do{Ee(j,a[m]|0);m=m+1|0}while((m|0)!=(k|0));k=(e|0)==-1?-1:e<<1;if((a[j]&1)==0){o=k;p=9}else{q=k;r=c[j+8>>2]|0}}else{o=(e|0)==-1?-1:e<<1;p=9}if((p|0)==9){q=o;r=j+1|0}o=kb(q|0,f|0,g|0,r|0)|0;c[b+0>>2]=0;c[b+4>>2]=0;c[b+8>>2]=0;r=mn(o|0)|0;g=o+r|0;if((r|0)>0){s=o}else{ze(j);i=d;return}do{Ee(b,a[s]|0);s=s+1|0}while((s|0)!=(g|0));ze(j);i=d;return}function Lj(a,b){a=a|0;b=b|0;a=i;Vb(((b|0)==-1?-1:b<<1)|0)|0;i=a;return}function Mj(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function Nj(a){a=a|0;return}function Oj(b,d,e){b=b|0;d=d|0;e=e|0;var f=0;e=i;if((a[d]&1)==0){f=d+1|0}else{f=c[d+8>>2]|0}d=bc(f|0,1)|0;i=e;return d>>>((d|0)!=(-1|0)|0)|0}function Pj(b,d,e,f,g,h){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0;d=i;i=i+176|0;j=d;k=d+48|0;l=d+8|0;m=d+12|0;n=d+16|0;o=d+32|0;p=d+40|0;c[n+0>>2]=0;c[n+4>>2]=0;c[n+8>>2]=0;c[o+4>>2]=0;c[o>>2]=9912;q=a[h]|0;if((q&1)==0){r=h+4|0;s=(q&255)>>>1;t=h+4|0}else{q=c[h+8>>2]|0;r=q;s=c[h+4>>2]|0;t=q}q=r+(s<<2)|0;s=j;c[s>>2]=0;c[s+4>>2]=0;a:do{if(t>>>0<q>>>0){s=k+32|0;r=t;h=9912|0;while(1){c[m>>2]=r;u=(rc[c[h+12>>2]&15](o,j,r,q,m,k,s,l)|0)==2;v=c[m>>2]|0;if(u|(v|0)==(r|0)){break}if(k>>>0<(c[l>>2]|0)>>>0){u=k;do{Ee(n,a[u]|0);u=u+1|0}while(u>>>0<(c[l>>2]|0)>>>0);w=c[m>>2]|0}else{w=v}if(!(w>>>0<q>>>0)){break a}r=w;h=c[o>>2]|0}jj(9136)}}while(0);if((a[n]&1)==0){x=n+1|0}else{x=c[n+8>>2]|0}o=kb(((e|0)==-1?-1:e<<1)|0,f|0,g|0,x|0)|0;c[b+0>>2]=0;c[b+4>>2]=0;c[b+8>>2]=0;c[p+4>>2]=0;c[p>>2]=10016;x=mn(o|0)|0;g=o+x|0;f=j;c[f>>2]=0;c[f+4>>2]=0;if((x|0)<=0){ze(n);i=d;return}x=g;f=k+128|0;e=o;o=10016|0;while(1){c[m>>2]=e;w=(rc[c[o+16>>2]&15](p,j,e,(x-e|0)>32?e+32|0:g,m,k,f,l)|0)==2;q=c[m>>2]|0;if(w|(q|0)==(e|0)){y=20;break}if(k>>>0<(c[l>>2]|0)>>>0){w=k;do{Oe(b,c[w>>2]|0);w=w+4|0}while(w>>>0<(c[l>>2]|0)>>>0);z=c[m>>2]|0}else{z=q}if(!(z>>>0<g>>>0)){y=25;break}e=z;o=c[p>>2]|0}if((y|0)==20){jj(9136)}else if((y|0)==25){ze(n);i=d;return}}function Qj(a,b){a=a|0;b=b|0;a=i;Vb(((b|0)==-1?-1:b<<1)|0)|0;i=a;return}function Rj(b){b=b|0;var d=0,e=0;d=i;c[b>>2]=8344;e=b+8|0;b=c[e>>2]|0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}if((b|0)==(c[2036]|0)){i=d;return}Ea(c[e>>2]|0);i=d;return}function Sj(a){a=a|0;a=Wa(8)|0;be(a,8136);c[a>>2]=5184;xb(a|0,5224,9)}function Tj(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0;e=i;i=i+16|0;f=e;c[b+4>>2]=d+ -1;c[b>>2]=8176;d=b+8|0;g=b+12|0;h=b+136|0;j=b+24|0;a[h]=1;c[g>>2]=j;c[d>>2]=j;c[b+16>>2]=h;h=28;k=j;do{if((k|0)==0){l=0}else{c[k>>2]=0;l=c[g>>2]|0}k=l+4|0;c[g>>2]=k;h=h+ -1|0}while((h|0)!=0);xe(b+144|0,8160,1);h=c[d>>2]|0;d=c[g>>2]|0;if((d|0)!=(h|0)){c[g>>2]=d+(~((d+ -4+(0-h)|0)>>>2)<<2)}c[13068>>2]=0;c[3266]=6656;if(!((c[1670]|0)==-1)){c[f>>2]=6680;c[f+4>>2]=112;c[f+8>>2]=0;ue(6680,f,113)}Uj(b,13064,(c[6684>>2]|0)+ -1|0);c[13060>>2]=0;c[3264]=6696;if(!((c[1680]|0)==-1)){c[f>>2]=6720;c[f+4>>2]=112;c[f+8>>2]=0;ue(6720,f,113)}Uj(b,13056,(c[6724>>2]|0)+ -1|0);c[13044>>2]=0;c[3260]=8272;c[13048>>2]=0;a[13052|0]=0;c[13048>>2]=c[(La()|0)>>2];if(!((c[2064]|0)==-1)){c[f>>2]=8256;c[f+4>>2]=112;c[f+8>>2]=0;ue(8256,f,113)}Uj(b,13040,(c[8260>>2]|0)+ -1|0);c[13036>>2]=0;c[3258]=9232;if(!((c[2062]|0)==-1)){c[f>>2]=8248;c[f+4>>2]=112;c[f+8>>2]=0;ue(8248,f,113)}Uj(b,13032,(c[8252>>2]|0)+ -1|0);c[13028>>2]=0;c[3256]=9448;if(!((c[2080]|0)==-1)){c[f>>2]=8320;c[f+4>>2]=112;c[f+8>>2]=0;ue(8320,f,113)}Uj(b,13024,(c[8324>>2]|0)+ -1|0);c[13012>>2]=0;c[3252]=8344;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}c[13016>>2]=c[2036];if(!((c[2082]|0)==-1)){c[f>>2]=8328;c[f+4>>2]=112;c[f+8>>2]=0;ue(8328,f,113)}Uj(b,13008,(c[8332>>2]|0)+ -1|0);c[13004>>2]=0;c[3250]=9672;if(!((c[2096]|0)==-1)){c[f>>2]=8384;c[f+4>>2]=112;c[f+8>>2]=0;ue(8384,f,113)}Uj(b,13e3,(c[8388>>2]|0)+ -1|0);c[12996>>2]=0;c[3248]=9792;if(!((c[2098]|0)==-1)){c[f>>2]=8392;c[f+4>>2]=112;c[f+8>>2]=0;ue(8392,f,113)}Uj(b,12992,(c[8396>>2]|0)+ -1|0);c[12972>>2]=0;c[3242]=8424;a[12976|0]=46;a[12977|0]=44;c[12980>>2]=0;c[12984>>2]=0;c[12988>>2]=0;if(!((c[2100]|0)==-1)){c[f>>2]=8400;c[f+4>>2]=112;c[f+8>>2]=0;ue(8400,f,113)}Uj(b,12968,(c[8404>>2]|0)+ -1|0);c[12940>>2]=0;c[3234]=8464;c[12944>>2]=46;c[12948>>2]=44;c[12952>>2]=0;c[12956>>2]=0;c[12960>>2]=0;if(!((c[2102]|0)==-1)){c[f>>2]=8408;c[f+4>>2]=112;c[f+8>>2]=0;ue(8408,f,113)}Uj(b,12936,(c[8412>>2]|0)+ -1|0);c[12932>>2]=0;c[3232]=6736;if(!((c[1698]|0)==-1)){c[f>>2]=6792;c[f+4>>2]=112;c[f+8>>2]=0;ue(6792,f,113)}Uj(b,12928,(c[6796>>2]|0)+ -1|0);c[12924>>2]=0;c[3230]=6856;if(!((c[1728]|0)==-1)){c[f>>2]=6912;c[f+4>>2]=112;c[f+8>>2]=0;ue(6912,f,113)}Uj(b,12920,(c[6916>>2]|0)+ -1|0);c[12916>>2]=0;c[3228]=6928;if(!((c[1744]|0)==-1)){c[f>>2]=6976;c[f+4>>2]=112;c[f+8>>2]=0;ue(6976,f,113)}Uj(b,12912,(c[6980>>2]|0)+ -1|0);c[12908>>2]=0;c[3226]=6992;if(!((c[1760]|0)==-1)){c[f>>2]=7040;c[f+4>>2]=112;c[f+8>>2]=0;ue(7040,f,113)}Uj(b,12904,(c[7044>>2]|0)+ -1|0);c[12900>>2]=0;c[3224]=7584;if(!((c[1908]|0)==-1)){c[f>>2]=7632;c[f+4>>2]=112;c[f+8>>2]=0;ue(7632,f,113)}Uj(b,12896,(c[7636>>2]|0)+ -1|0);c[12892>>2]=0;c[3222]=7648;if(!((c[1924]|0)==-1)){c[f>>2]=7696;c[f+4>>2]=112;c[f+8>>2]=0;ue(7696,f,113)}Uj(b,12888,(c[7700>>2]|0)+ -1|0);c[12884>>2]=0;c[3220]=7712;if(!((c[1940]|0)==-1)){c[f>>2]=7760;c[f+4>>2]=112;c[f+8>>2]=0;ue(7760,f,113)}Uj(b,12880,(c[7764>>2]|0)+ -1|0);c[12876>>2]=0;c[3218]=7776;if(!((c[1956]|0)==-1)){c[f>>2]=7824;c[f+4>>2]=112;c[f+8>>2]=0;ue(7824,f,113)}Uj(b,12872,(c[7828>>2]|0)+ -1|0);c[12868>>2]=0;c[3216]=7840;if(!((c[1966]|0)==-1)){c[f>>2]=7864;c[f+4>>2]=112;c[f+8>>2]=0;ue(7864,f,113)}Uj(b,12864,(c[7868>>2]|0)+ -1|0);c[12860>>2]=0;c[3214]=7920;if(!((c[1986]|0)==-1)){c[f>>2]=7944;c[f+4>>2]=112;c[f+8>>2]=0;ue(7944,f,113)}Uj(b,12856,(c[7948>>2]|0)+ -1|0);c[12852>>2]=0;c[3212]=7976;if(!((c[2e3]|0)==-1)){c[f>>2]=8e3;c[f+4>>2]=112;c[f+8>>2]=0;ue(8e3,f,113)}Uj(b,12848,(c[8004>>2]|0)+ -1|0);c[12844>>2]=0;c[3210]=8024;if(!((c[2012]|0)==-1)){c[f>>2]=8048;c[f+4>>2]=112;c[f+8>>2]=0;ue(8048,f,113)}Uj(b,12840,(c[8052>>2]|0)+ -1|0);c[12828>>2]=0;c[3206]=7072;c[12832>>2]=7120;if(!((c[1788]|0)==-1)){c[f>>2]=7152;c[f+4>>2]=112;c[f+8>>2]=0;ue(7152,f,113)}Uj(b,12824,(c[7156>>2]|0)+ -1|0);c[12812>>2]=0;c[3202]=7224;c[12816>>2]=7272;if(!((c[1826]|0)==-1)){c[f>>2]=7304;c[f+4>>2]=112;c[f+8>>2]=0;ue(7304,f,113)}Uj(b,12808,(c[7308>>2]|0)+ -1|0);c[12796>>2]=0;c[3198]=9168;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}c[12800>>2]=c[2036];c[3198]=7520;if(!((c[1884]|0)==-1)){c[f>>2]=7536;c[f+4>>2]=112;c[f+8>>2]=0;ue(7536,f,113)}Uj(b,12792,(c[7540>>2]|0)+ -1|0);c[12780>>2]=0;c[3194]=9168;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}c[12784>>2]=c[2036];c[3194]=7552;if(!((c[1892]|0)==-1)){c[f>>2]=7568;c[f+4>>2]=112;c[f+8>>2]=0;ue(7568,f,113)}Uj(b,12776,(c[7572>>2]|0)+ -1|0);c[12772>>2]=0;c[3192]=8064;if(!((c[2022]|0)==-1)){c[f>>2]=8088;c[f+4>>2]=112;c[f+8>>2]=0;ue(8088,f,113)}Uj(b,12768,(c[8092>>2]|0)+ -1|0);c[12764>>2]=0;c[3190]=8104;if((c[2032]|0)==-1){m=c[8132>>2]|0;n=m+ -1|0;Uj(b,12760,n);i=e;return}c[f>>2]=8128;c[f+4>>2]=112;c[f+8>>2]=0;ue(8128,f,113);m=c[8132>>2]|0;n=m+ -1|0;Uj(b,12760,n);i=e;return}function Uj(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0;e=i;$d(b);f=a+8|0;g=a+12|0;a=c[g>>2]|0;h=c[f>>2]|0;j=a-h>>2;do{if(!(j>>>0>d>>>0)){k=d+1|0;if(j>>>0<k>>>0){Yl(f,k-j|0);l=c[f>>2]|0;break}if(j>>>0>k>>>0?(m=h+(k<<2)|0,(a|0)!=(m|0)):0){c[g>>2]=a+(~((a+ -4+(0-m)|0)>>>2)<<2);l=h}else{l=h}}else{l=h}}while(0);h=c[l+(d<<2)>>2]|0;if((h|0)==0){n=l;o=n+(d<<2)|0;c[o>>2]=b;i=e;return}ae(h)|0;n=c[f>>2]|0;o=n+(d<<2)|0;c[o>>2]=b;i=e;return}function Vj(a){a=a|0;var b=0;b=i;Wj(a);Um(a);i=b;return}function Wj(b){b=b|0;var d=0,e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0;d=i;c[b>>2]=8176;e=b+12|0;f=c[e>>2]|0;g=b+8|0;h=c[g>>2]|0;if((f|0)!=(h|0)){j=f;f=h;h=0;while(1){k=c[f+(h<<2)>>2]|0;if((k|0)==0){l=j;m=f}else{ae(k)|0;l=c[e>>2]|0;m=c[g>>2]|0}h=h+1|0;if(!(h>>>0<l-m>>2>>>0)){break}else{j=l;f=m}}}ze(b+144|0);m=c[g>>2]|0;if((m|0)==0){i=d;return}g=c[e>>2]|0;if((g|0)!=(m|0)){c[e>>2]=g+(~((g+ -4+(0-m)|0)>>>2)<<2)}if((b+24|0)==(m|0)){a[b+136|0]=0;i=d;return}else{Um(m);i=d;return}}function Xj(){var b=0,d=0,e=0;b=i;if((a[8232]|0)!=0){d=c[2056]|0;i=b;return d|0}if((Oa(8232)|0)==0){d=c[2056]|0;i=b;return d|0}if((a[8208]|0)==0?(Oa(8208)|0)!=0:0){Tj(12600,1);c[2048]=12600;c[2050]=8192;$a(8208)}e=c[c[2050]>>2]|0;c[2054]=e;$d(e);c[2056]=8216;$a(8232);d=c[2056]|0;i=b;return d|0}function Yj(a){a=a|0;var b=0,d=0;b=i;d=c[(Xj()|0)>>2]|0;c[a>>2]=d;$d(d);i=b;return}function Zj(a,b){a=a|0;b=b|0;var d=0,e=0;d=i;e=c[b>>2]|0;c[a>>2]=e;$d(e);i=d;return}function _j(a){a=a|0;var b=0;b=i;ae(c[a>>2]|0)|0;i=b;return}function $j(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,g=0;d=i;i=i+16|0;e=d;f=c[a>>2]|0;if(!((c[b>>2]|0)==-1)){c[e>>2]=b;c[e+4>>2]=112;c[e+8>>2]=0;ue(b,e,113)}e=(c[b+4>>2]|0)+ -1|0;b=c[f+8>>2]|0;if(!((c[f+12>>2]|0)-b>>2>>>0>e>>>0)){g=Wa(4)|0;qm(g);xb(g|0,16216,101)}f=c[b+(e<<2)>>2]|0;if((f|0)==0){g=Wa(4)|0;qm(g);xb(g|0,16216,101)}else{i=d;return f|0}return 0}function ak(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function bk(a){a=a|0;var b=0;b=i;if((a|0)==0){i=b;return}jc[c[(c[a>>2]|0)+4>>2]&127](a);i=b;return}function ck(a){a=a|0;var b=0;b=c[2060]|0;c[2060]=b+1;c[a+4>>2]=b+1;return}function dk(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function ek(a,d,e){a=a|0;d=d|0;e=e|0;var f=0;a=i;if(!(e>>>0<128)){f=0;i=a;return f|0}f=(b[(c[(La()|0)>>2]|0)+(e<<1)>>1]&d)<<16>>16!=0;i=a;return f|0}function fk(a,d,e,f){a=a|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0;a=i;if((d|0)==(e|0)){g=d;i=a;return g|0}else{h=d;j=f}while(1){f=c[h>>2]|0;if(f>>>0<128){k=b[(c[(La()|0)>>2]|0)+(f<<1)>>1]|0}else{k=0}b[j>>1]=k;f=h+4|0;if((f|0)==(e|0)){g=e;break}else{h=f;j=j+2|0}}i=a;return g|0}function gk(a,d,e,f){a=a|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0;a=i;a:do{if((e|0)==(f|0)){g=e}else{h=e;while(1){j=c[h>>2]|0;if(j>>>0<128?!((b[(c[(La()|0)>>2]|0)+(j<<1)>>1]&d)<<16>>16==0):0){g=h;break a}j=h+4|0;if((j|0)==(f|0)){g=f;break}else{h=j}}}}while(0);i=a;return g|0}function hk(a,d,e,f){a=a|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0;a=i;a:do{if((e|0)==(f|0)){g=e}else{h=e;while(1){j=c[h>>2]|0;if(!(j>>>0<128)){g=h;break a}k=h+4|0;if((b[(c[(La()|0)>>2]|0)+(j<<1)>>1]&d)<<16>>16==0){g=h;break a}if((k|0)==(f|0)){g=f;break}else{h=k}}}}while(0);i=a;return g|0}function ik(a,b){a=a|0;b=b|0;var d=0;a=i;if(!(b>>>0<128)){d=b;i=a;return d|0}d=c[(c[(Za()|0)>>2]|0)+(b<<2)>>2]|0;i=a;return d|0}function jk(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,g=0;a=i;if((b|0)==(d|0)){e=b;i=a;return e|0}else{f=b}while(1){b=c[f>>2]|0;if(b>>>0<128){g=c[(c[(Za()|0)>>2]|0)+(b<<2)>>2]|0}else{g=b}c[f>>2]=g;b=f+4|0;if((b|0)==(d|0)){e=d;break}else{f=b}}i=a;return e|0}function kk(a,b){a=a|0;b=b|0;var d=0;a=i;if(!(b>>>0<128)){d=b;i=a;return d|0}d=c[(c[(yb()|0)>>2]|0)+(b<<2)>>2]|0;i=a;return d|0}function lk(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,g=0;a=i;if((b|0)==(d|0)){e=b;i=a;return e|0}else{f=b}while(1){b=c[f>>2]|0;if(b>>>0<128){g=c[(c[(yb()|0)>>2]|0)+(b<<2)>>2]|0}else{g=b}c[f>>2]=g;b=f+4|0;if((b|0)==(d|0)){e=d;break}else{f=b}}i=a;return e|0}function mk(a,b){a=a|0;b=b|0;return b<<24>>24|0}function nk(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0;b=i;if((d|0)==(e|0)){g=d;i=b;return g|0}else{h=d;j=f}while(1){c[j>>2]=a[h]|0;f=h+1|0;if((f|0)==(e|0)){g=e;break}else{h=f;j=j+4|0}}i=b;return g|0}function ok(a,b,c){a=a|0;b=b|0;c=c|0;return(b>>>0<128?b&255:c)|0}function pk(b,d,e,f,g){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;var h=0,j=0,k=0,l=0;b=i;if((d|0)==(e|0)){h=d;i=b;return h|0}j=((e+ -4+(0-d)|0)>>>2)+1|0;k=d;l=g;while(1){g=c[k>>2]|0;a[l]=g>>>0<128?g&255:f;k=k+4|0;if((k|0)==(e|0)){break}else{l=l+1|0}}h=d+(j<<2)|0;i=b;return h|0}function qk(b){b=b|0;var d=0,e=0;d=i;c[b>>2]=8272;e=c[b+8>>2]|0;if((e|0)!=0?(a[b+12|0]|0)!=0:0){Vm(e)}Um(b);i=d;return}function rk(b){b=b|0;var d=0,e=0;d=i;c[b>>2]=8272;e=c[b+8>>2]|0;if((e|0)!=0?(a[b+12|0]|0)!=0:0){Vm(e)}i=d;return}function sk(a,b){a=a|0;b=b|0;var d=0;a=i;if(!(b<<24>>24>-1)){d=b;i=a;return d|0}d=c[(c[(Za()|0)>>2]|0)+((b&255)<<2)>>2]&255;i=a;return d|0}function tk(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0;b=i;if((d|0)==(e|0)){f=d;i=b;return f|0}else{g=d}while(1){d=a[g]|0;if(d<<24>>24>-1){h=c[(c[(Za()|0)>>2]|0)+(d<<24>>24<<2)>>2]&255}else{h=d}a[g]=h;d=g+1|0;if((d|0)==(e|0)){f=e;break}else{g=d}}i=b;return f|0}function uk(a,b){a=a|0;b=b|0;var d=0;a=i;if(!(b<<24>>24>-1)){d=b;i=a;return d|0}d=c[(c[(yb()|0)>>2]|0)+(b<<24>>24<<2)>>2]&255;i=a;return d|0}function vk(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0;b=i;if((d|0)==(e|0)){f=d;i=b;return f|0}else{g=d}while(1){d=a[g]|0;if(d<<24>>24>-1){h=c[(c[(yb()|0)>>2]|0)+(d<<24>>24<<2)>>2]&255}else{h=d}a[g]=h;d=g+1|0;if((d|0)==(e|0)){f=e;break}else{g=d}}i=b;return f|0}function wk(a,b){a=a|0;b=b|0;return b|0}function xk(b,c,d,e){b=b|0;c=c|0;d=d|0;e=e|0;var f=0,g=0;b=i;if((c|0)==(d|0)){f=c}else{g=c;c=e;while(1){a[c]=a[g]|0;e=g+1|0;if((e|0)==(d|0)){f=d;break}else{g=e;c=c+1|0}}}i=b;return f|0}function yk(a,b,c){a=a|0;b=b|0;c=c|0;return(b<<24>>24>-1?b:c)|0}function zk(b,c,d,e,f){b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0;b=i;if((c|0)==(d|0)){g=c;i=b;return g|0}else{h=c;j=f}while(1){f=a[h]|0;a[j]=f<<24>>24>-1?f:e;f=h+1|0;if((f|0)==(d|0)){g=d;break}else{h=f;j=j+1|0}}i=b;return g|0}function Ak(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function Bk(a,b,d,e,f,g,h,i){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;i=i|0;c[f>>2]=d;c[i>>2]=g;return 3}function Ck(a,b,d,e,f,g,h,i){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;i=i|0;c[f>>2]=d;c[i>>2]=g;return 3}function Dk(a,b,d,e,f){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;c[f>>2]=d;return 3}function Ek(a){a=a|0;return 1}function Fk(a){a=a|0;return 1}function Gk(a,b,c,d,e){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;b=d-c|0;return(b>>>0<e>>>0?b:e)|0}function Hk(a){a=a|0;return 1}function Ik(a){a=a|0;var b=0;b=i;Rj(a);Um(a);i=b;return}function Jk(b,d,e,f,g,h,j,k){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0;l=i;i=i+16|0;m=l;n=l+8|0;o=(e|0)==(f|0);a:do{if(!o){p=e;while(1){q=p+4|0;if((c[p>>2]|0)==0){r=p;break}if((q|0)==(f|0)){r=f;break}else{p=q}}c[k>>2]=h;c[g>>2]=e;if(!(o|(h|0)==(j|0))){p=j;q=b+8|0;s=e;t=h;u=r;while(1){v=d;w=c[v+4>>2]|0;x=m;c[x>>2]=c[v>>2];c[x+4>>2]=w;w=jb(c[q>>2]|0)|0;x=jm(t,g,u-s>>2,p-t|0,d)|0;if((w|0)!=0){jb(w|0)|0}if((x|0)==0){y=1;z=33;break}else if((x|0)==-1){z=10;break}w=(c[k>>2]|0)+x|0;c[k>>2]=w;if((w|0)==(j|0)){z=31;break}if((u|0)==(f|0)){A=c[g>>2]|0;B=w;C=f}else{w=jb(c[q>>2]|0)|0;x=im(n,0,d)|0;if((w|0)!=0){jb(w|0)|0}if((x|0)==-1){y=2;z=33;break}w=c[k>>2]|0;if(x>>>0>(p-w|0)>>>0){y=1;z=33;break}b:do{if((x|0)!=0){v=w;D=x;E=n;while(1){F=a[E]|0;c[k>>2]=v+1;a[v]=F;F=D+ -1|0;if((F|0)==0){break b}v=c[k>>2]|0;D=F;E=E+1|0}}}while(0);x=(c[g>>2]|0)+4|0;c[g>>2]=x;c:do{if((x|0)==(f|0)){G=f}else{w=x;while(1){E=w+4|0;if((c[w>>2]|0)==0){G=w;break c}if((E|0)==(f|0)){G=f;break}else{w=E}}}}while(0);A=x;B=c[k>>2]|0;C=G}if((A|0)==(f|0)|(B|0)==(j|0)){H=A;break a}else{s=A;t=B;u=C}}if((z|0)==10){c[k>>2]=t;d:do{if((s|0)==(c[g>>2]|0)){I=s}else{u=s;p=t;while(1){w=c[u>>2]|0;E=jb(c[q>>2]|0)|0;D=im(p,w,m)|0;if((E|0)!=0){jb(E|0)|0}if((D|0)==-1){I=u;break d}E=(c[k>>2]|0)+D|0;c[k>>2]=E;D=u+4|0;if((D|0)==(c[g>>2]|0)){I=D;break}else{u=D;p=E}}}}while(0);c[g>>2]=I;y=2;i=l;return y|0}else if((z|0)==31){H=c[g>>2]|0;break}else if((z|0)==33){i=l;return y|0}}else{H=e}}else{c[k>>2]=h;c[g>>2]=e;H=e}}while(0);y=(H|0)!=(f|0)|0;i=l;return y|0}function Kk(b,d,e,f,g,h,j,k){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0;l=i;i=i+16|0;m=l;n=(e|0)==(f|0);a:do{if(!n){o=e;while(1){p=o+1|0;if((a[o]|0)==0){q=o;break}if((p|0)==(f|0)){q=f;break}else{o=p}}c[k>>2]=h;c[g>>2]=e;if(!(n|(h|0)==(j|0))){o=j;p=b+8|0;r=e;s=h;t=q;while(1){u=d;v=c[u+4>>2]|0;w=m;c[w>>2]=c[u>>2];c[w+4>>2]=v;x=t;v=jb(c[p>>2]|0)|0;w=fm(s,g,x-r|0,o-s>>2,d)|0;if((v|0)!=0){jb(v|0)|0}if((w|0)==0){y=2;z=32;break}else if((w|0)==-1){z=10;break}v=(c[k>>2]|0)+(w<<2)|0;c[k>>2]=v;if((v|0)==(j|0)){z=30;break}w=c[g>>2]|0;if((t|0)==(f|0)){A=w;B=v;C=f}else{u=jb(c[p>>2]|0)|0;D=em(v,w,1,d)|0;if((u|0)!=0){jb(u|0)|0}if((D|0)!=0){y=2;z=32;break}c[k>>2]=(c[k>>2]|0)+4;D=(c[g>>2]|0)+1|0;c[g>>2]=D;b:do{if((D|0)==(f|0)){E=f}else{u=D;while(1){w=u+1|0;if((a[u]|0)==0){E=u;break b}if((w|0)==(f|0)){E=f;break}else{u=w}}}}while(0);A=D;B=c[k>>2]|0;C=E}if((A|0)==(f|0)|(B|0)==(j|0)){F=A;break a}else{r=A;s=B;t=C}}if((z|0)==10){c[k>>2]=s;c:do{if((r|0)!=(c[g>>2]|0)){t=r;o=s;while(1){u=jb(c[p>>2]|0)|0;w=em(o,t,x-t|0,m)|0;if((u|0)!=0){jb(u|0)|0}if((w|0)==-2){z=16;break}else if((w|0)==0){G=t+1|0}else if((w|0)==-1){z=15;break}else{G=t+w|0}w=(c[k>>2]|0)+4|0;c[k>>2]=w;if((G|0)==(c[g>>2]|0)){H=G;break c}else{t=G;o=w}}if((z|0)==15){c[g>>2]=t;y=2;i=l;return y|0}else if((z|0)==16){c[g>>2]=t;y=1;i=l;return y|0}}else{H=r}}while(0);c[g>>2]=H;y=(H|0)!=(f|0)|0;i=l;return y|0}else if((z|0)==30){F=c[g>>2]|0;break}else if((z|0)==32){i=l;return y|0}}else{F=e}}else{c[k>>2]=h;c[g>>2]=e;F=e}}while(0);y=(F|0)!=(f|0)|0;i=l;return y|0}function Lk(b,d,e,f,g){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;var h=0,j=0,k=0,l=0,m=0,n=0;h=i;i=i+16|0;j=h;c[g>>2]=e;e=jb(c[b+8>>2]|0)|0;b=im(j,0,d)|0;if((e|0)!=0){jb(e|0)|0}if((b|0)==0|(b|0)==-1){k=2;i=h;return k|0}e=b+ -1|0;b=c[g>>2]|0;if(e>>>0>(f-b|0)>>>0){k=1;i=h;return k|0}if((e|0)==0){k=0;i=h;return k|0}else{l=b;m=e;n=j}while(1){j=a[n]|0;c[g>>2]=l+1;a[l]=j;j=m+ -1|0;if((j|0)==0){k=0;break}l=c[g>>2]|0;m=j;n=n+1|0}i=h;return k|0}function Mk(a){a=a|0;var b=0,d=0,e=0,f=0;b=i;d=a+8|0;a=jb(c[d>>2]|0)|0;e=hm(0,0,4)|0;if((a|0)!=0){jb(a|0)|0}if((e|0)==0){e=c[d>>2]|0;if((e|0)!=0){d=jb(e|0)|0;if((d|0)==0){f=0}else{jb(d|0)|0;f=0}}else{f=1}}else{f=-1}i=b;return f|0}function Nk(a){a=a|0;return 0}function Ok(a,b,d,e,f){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0;g=i;if((f|0)==0|(d|0)==(e|0)){h=0;i=g;return h|0}j=e;k=a+8|0;a=d;d=0;l=0;while(1){m=jb(c[k>>2]|0)|0;n=dm(a,j-a|0,b)|0;if((m|0)!=0){jb(m|0)|0}if((n|0)==-2|(n|0)==-1){h=d;o=9;break}else if((n|0)==0){p=a+1|0;q=1}else{p=a+n|0;q=n}n=q+d|0;m=l+1|0;if(m>>>0>=f>>>0|(p|0)==(e|0)){h=n;o=9;break}else{a=p;d=n;l=m}}if((o|0)==9){i=g;return h|0}return 0}function Pk(a){a=a|0;var b=0,d=0,e=0;b=i;d=c[a+8>>2]|0;if((d|0)!=0){a=jb(d|0)|0;if((a|0)==0){e=4}else{jb(a|0)|0;e=4}}else{e=1}i=b;return e|0}function Qk(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function Rk(a,b,d,e,f,g,h,j){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0;b=i;i=i+16|0;a=b+4|0;k=b;c[a>>2]=d;c[k>>2]=g;l=Sk(d,e,a,g,h,k,1114111,0)|0;c[f>>2]=d+((c[a>>2]|0)-d>>1<<1);c[j>>2]=g+((c[k>>2]|0)-g);i=b;return l|0}function Sk(d,f,g,h,j,k,l,m){d=d|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;m=m|0;var n=0,o=0,p=0,q=0,r=0,s=0,t=0;n=i;c[g>>2]=d;c[k>>2]=h;do{if((m&2|0)!=0){if((j-h|0)<3){o=1;i=n;return o|0}else{c[k>>2]=h+1;a[h]=-17;d=c[k>>2]|0;c[k>>2]=d+1;a[d]=-69;d=c[k>>2]|0;c[k>>2]=d+1;a[d]=-65;break}}}while(0);h=f;m=c[g>>2]|0;if(!(m>>>0<f>>>0)){o=0;i=n;return o|0}d=j;j=m;a:while(1){m=b[j>>1]|0;p=m&65535;if(p>>>0>l>>>0){o=2;q=26;break}do{if((m&65535)<128){r=c[k>>2]|0;if((d-r|0)<1){o=1;q=26;break a}c[k>>2]=r+1;a[r]=m}else{if((m&65535)<2048){r=c[k>>2]|0;if((d-r|0)<2){o=1;q=26;break a}c[k>>2]=r+1;a[r]=p>>>6|192;r=c[k>>2]|0;c[k>>2]=r+1;a[r]=p&63|128;break}if((m&65535)<55296){r=c[k>>2]|0;if((d-r|0)<3){o=1;q=26;break a}c[k>>2]=r+1;a[r]=p>>>12|224;r=c[k>>2]|0;c[k>>2]=r+1;a[r]=p>>>6&63|128;r=c[k>>2]|0;c[k>>2]=r+1;a[r]=p&63|128;break}if(!((m&65535)<56320)){if((m&65535)<57344){o=2;q=26;break a}r=c[k>>2]|0;if((d-r|0)<3){o=1;q=26;break a}c[k>>2]=r+1;a[r]=p>>>12|224;r=c[k>>2]|0;c[k>>2]=r+1;a[r]=p>>>6&63|128;r=c[k>>2]|0;c[k>>2]=r+1;a[r]=p&63|128;break}if((h-j|0)<4){o=1;q=26;break a}r=j+2|0;s=e[r>>1]|0;if((s&64512|0)!=56320){o=2;q=26;break a}if((d-(c[k>>2]|0)|0)<4){o=1;q=26;break a}t=p&960;if(((t<<10)+65536|p<<10&64512|s&1023)>>>0>l>>>0){o=2;q=26;break a}c[g>>2]=r;r=(t>>>6)+1|0;t=c[k>>2]|0;c[k>>2]=t+1;a[t]=r>>>2|240;t=c[k>>2]|0;c[k>>2]=t+1;a[t]=p>>>2&15|r<<4&48|128;r=c[k>>2]|0;c[k>>2]=r+1;a[r]=p<<4&48|s>>>6&15|128;r=c[k>>2]|0;c[k>>2]=r+1;a[r]=s&63|128}}while(0);p=(c[g>>2]|0)+2|0;c[g>>2]=p;if(p>>>0<f>>>0){j=p}else{o=0;q=26;break}}if((q|0)==26){i=n;return o|0}return 0}function Tk(a,b,d,e,f,g,h,j){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0;b=i;i=i+16|0;a=b+4|0;k=b;c[a>>2]=d;c[k>>2]=g;l=Uk(d,e,a,g,h,k,1114111,0)|0;c[f>>2]=d+((c[a>>2]|0)-d);c[j>>2]=g+((c[k>>2]|0)-g>>1<<1);i=b;return l|0}function Uk(e,f,g,h,j,k,l,m){e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;m=m|0;var n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0;n=i;c[g>>2]=e;c[k>>2]=h;h=c[g>>2]|0;if(((((m&4|0)!=0?(f-h|0)>2:0)?(a[h]|0)==-17:0)?(a[h+1|0]|0)==-69:0)?(a[h+2|0]|0)==-65:0){m=h+3|0;c[g>>2]=m;o=m}else{o=h}a:do{if(o>>>0<f>>>0){h=f;m=j;e=c[k>>2]|0;p=o;b:while(1){if(!(e>>>0<j>>>0)){q=p;break a}r=a[p]|0;s=r&255;if(s>>>0>l>>>0){t=2;u=41;break}do{if(r<<24>>24>-1){b[e>>1]=r&255;c[g>>2]=p+1}else{if((r&255)<194){t=2;u=41;break b}if((r&255)<224){if((h-p|0)<2){t=1;u=41;break b}v=d[p+1|0]|0;if((v&192|0)!=128){t=2;u=41;break b}w=v&63|s<<6&1984;if(w>>>0>l>>>0){t=2;u=41;break b}b[e>>1]=w;c[g>>2]=p+2;break}if((r&255)<240){if((h-p|0)<3){t=1;u=41;break b}w=a[p+1|0]|0;v=a[p+2|0]|0;if((s|0)==224){if(!((w&-32)<<24>>24==-96)){t=2;u=41;break b}}else if((s|0)==237){if(!((w&-32)<<24>>24==-128)){t=2;u=41;break b}}else{if(!((w&-64)<<24>>24==-128)){t=2;u=41;break b}}x=v&255;if((x&192|0)!=128){t=2;u=41;break b}v=(w&255)<<6&4032|s<<12|x&63;if((v&65535)>>>0>l>>>0){t=2;u=41;break b}b[e>>1]=v;c[g>>2]=p+3;break}if(!((r&255)<245)){t=2;u=41;break b}if((h-p|0)<4){t=1;u=41;break b}v=a[p+1|0]|0;x=a[p+2|0]|0;w=a[p+3|0]|0;if((s|0)==244){if(!((v&-16)<<24>>24==-128)){t=2;u=41;break b}}else if((s|0)==240){if(!((v+112<<24>>24&255)<48)){t=2;u=41;break b}}else{if(!((v&-64)<<24>>24==-128)){t=2;u=41;break b}}y=x&255;if((y&192|0)!=128){t=2;u=41;break b}x=w&255;if((x&192|0)!=128){t=2;u=41;break b}if((m-e|0)<4){t=1;u=41;break b}w=s&7;z=v&255;v=y<<6;A=x&63;if((z<<12&258048|w<<18|v&4032|A)>>>0>l>>>0){t=2;u=41;break b}b[e>>1]=z<<2&60|y>>>4&3|((z>>>4&3|w<<2)<<6)+16320|55296;w=e+2|0;c[k>>2]=w;b[w>>1]=A|v&960|56320;c[g>>2]=(c[g>>2]|0)+4}}while(0);s=(c[k>>2]|0)+2|0;c[k>>2]=s;r=c[g>>2]|0;if(r>>>0<f>>>0){e=s;p=r}else{q=r;break a}}if((u|0)==41){i=n;return t|0}}else{q=o}}while(0);t=q>>>0<f>>>0|0;i=n;return t|0}function Vk(a,b,d,e,f){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;c[f>>2]=d;return 3}function Wk(a){a=a|0;return 0}function Xk(a){a=a|0;return 0}function Yk(a,b,c,d,e){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;b=i;a=Zk(c,d,e,1114111,0)|0;i=b;return a|0}function Zk(b,c,e,f,g){b=b|0;c=c|0;e=e|0;f=f|0;g=g|0;var h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0;h=i;if((((g&4|0)!=0?(c-b|0)>2:0)?(a[b]|0)==-17:0)?(a[b+1|0]|0)==-69:0){j=(a[b+2|0]|0)==-65?b+3|0:b}else{j=b}a:do{if(j>>>0<c>>>0&(e|0)!=0){g=c;k=j;l=0;b:while(1){m=a[k]|0;n=m&255;if(n>>>0>f>>>0){o=k;break a}do{if(m<<24>>24>-1){p=k+1|0;q=l}else{if((m&255)<194){o=k;break a}if((m&255)<224){if((g-k|0)<2){o=k;break a}r=d[k+1|0]|0;if((r&192|0)!=128){o=k;break a}if((r&63|n<<6&1984)>>>0>f>>>0){o=k;break a}p=k+2|0;q=l;break}if((m&255)<240){s=k;if((g-s|0)<3){o=k;break a}r=a[k+1|0]|0;t=a[k+2|0]|0;if((n|0)==224){if(!((r&-32)<<24>>24==-96)){u=21;break b}}else if((n|0)==237){if(!((r&-32)<<24>>24==-128)){u=23;break b}}else{if(!((r&-64)<<24>>24==-128)){u=25;break b}}v=t&255;if((v&192|0)!=128){o=k;break a}if(((r&255)<<6&4032|n<<12&61440|v&63)>>>0>f>>>0){o=k;break a}p=k+3|0;q=l;break}if(!((m&255)<245)){o=k;break a}w=k;if((g-w|0)<4){o=k;break a}if((e-l|0)>>>0<2){o=k;break a}v=a[k+1|0]|0;r=a[k+2|0]|0;t=a[k+3|0]|0;if((n|0)==244){if(!((v&-16)<<24>>24==-128)){u=36;break b}}else if((n|0)==240){if(!((v+112<<24>>24&255)<48)){u=34;break b}}else{if(!((v&-64)<<24>>24==-128)){u=38;break b}}x=r&255;if((x&192|0)!=128){o=k;break a}r=t&255;if((r&192|0)!=128){o=k;break a}if(((v&255)<<12&258048|n<<18&1835008|x<<6&4032|r&63)>>>0>f>>>0){o=k;break a}p=k+4|0;q=l+1|0}}while(0);n=q+1|0;if(p>>>0<c>>>0&n>>>0<e>>>0){k=p;l=n}else{o=p;break a}}if((u|0)==21){y=s-b|0;i=h;return y|0}else if((u|0)==23){y=s-b|0;i=h;return y|0}else if((u|0)==25){y=s-b|0;i=h;return y|0}else if((u|0)==34){y=w-b|0;i=h;return y|0}else if((u|0)==36){y=w-b|0;i=h;return y|0}else if((u|0)==38){y=w-b|0;i=h;return y|0}}else{o=j}}while(0);y=o-b|0;i=h;return y|0}function _k(a){a=a|0;return 4}function $k(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function al(a,b,d,e,f,g,h,j){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0;b=i;i=i+16|0;a=b+4|0;k=b;c[a>>2]=d;c[k>>2]=g;l=bl(d,e,a,g,h,k,1114111,0)|0;c[f>>2]=d+((c[a>>2]|0)-d>>2<<2);c[j>>2]=g+((c[k>>2]|0)-g);i=b;return l|0}function bl(b,d,e,f,g,h,j,k){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;var l=0,m=0,n=0,o=0,p=0;l=i;c[e>>2]=b;c[h>>2]=f;do{if((k&2|0)!=0){if((g-f|0)<3){m=1;i=l;return m|0}else{c[h>>2]=f+1;a[f]=-17;b=c[h>>2]|0;c[h>>2]=b+1;a[b]=-69;b=c[h>>2]|0;c[h>>2]=b+1;a[b]=-65;break}}}while(0);f=c[e>>2]|0;if(!(f>>>0<d>>>0)){m=0;i=l;return m|0}k=g;g=f;a:while(1){f=c[g>>2]|0;if((f&-2048|0)==55296|f>>>0>j>>>0){m=2;n=19;break}do{if(!(f>>>0<128)){if(f>>>0<2048){b=c[h>>2]|0;if((k-b|0)<2){m=1;n=19;break a}c[h>>2]=b+1;a[b]=f>>>6|192;b=c[h>>2]|0;c[h>>2]=b+1;a[b]=f&63|128;break}b=c[h>>2]|0;o=k-b|0;if(f>>>0<65536){if((o|0)<3){m=1;n=19;break a}c[h>>2]=b+1;a[b]=f>>>12|224;p=c[h>>2]|0;c[h>>2]=p+1;a[p]=f>>>6&63|128;p=c[h>>2]|0;c[h>>2]=p+1;a[p]=f&63|128;break}else{if((o|0)<4){m=1;n=19;break a}c[h>>2]=b+1;a[b]=f>>>18|240;b=c[h>>2]|0;c[h>>2]=b+1;a[b]=f>>>12&63|128;b=c[h>>2]|0;c[h>>2]=b+1;a[b]=f>>>6&63|128;b=c[h>>2]|0;c[h>>2]=b+1;a[b]=f&63|128;break}}else{b=c[h>>2]|0;if((k-b|0)<1){m=1;n=19;break a}c[h>>2]=b+1;a[b]=f}}while(0);f=(c[e>>2]|0)+4|0;c[e>>2]=f;if(f>>>0<d>>>0){g=f}else{m=0;n=19;break}}if((n|0)==19){i=l;return m|0}return 0}function cl(a,b,d,e,f,g,h,j){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;var k=0,l=0;b=i;i=i+16|0;a=b+4|0;k=b;c[a>>2]=d;c[k>>2]=g;l=dl(d,e,a,g,h,k,1114111,0)|0;c[f>>2]=d+((c[a>>2]|0)-d);c[j>>2]=g+((c[k>>2]|0)-g>>2<<2);i=b;return l|0}function dl(b,e,f,g,h,j,k,l){b=b|0;e=e|0;f=f|0;g=g|0;h=h|0;j=j|0;k=k|0;l=l|0;var m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0;m=i;c[f>>2]=b;c[j>>2]=g;g=c[f>>2]|0;if(((((l&4|0)!=0?(e-g|0)>2:0)?(a[g]|0)==-17:0)?(a[g+1|0]|0)==-69:0)?(a[g+2|0]|0)==-65:0){l=g+3|0;c[f>>2]=l;n=l}else{n=g}a:do{if(n>>>0<e>>>0){g=e;l=c[j>>2]|0;b=n;while(1){if(!(l>>>0<h>>>0)){o=b;p=39;break a}q=a[b]|0;r=q&255;do{if(q<<24>>24>-1){if(r>>>0>k>>>0){s=2;break a}c[l>>2]=r;c[f>>2]=b+1}else{if((q&255)<194){s=2;break a}if((q&255)<224){if((g-b|0)<2){s=1;break a}t=d[b+1|0]|0;if((t&192|0)!=128){s=2;break a}u=t&63|r<<6&1984;if(u>>>0>k>>>0){s=2;break a}c[l>>2]=u;c[f>>2]=b+2;break}if((q&255)<240){if((g-b|0)<3){s=1;break a}u=a[b+1|0]|0;t=a[b+2|0]|0;if((r|0)==237){if(!((u&-32)<<24>>24==-128)){s=2;break a}}else if((r|0)==224){if(!((u&-32)<<24>>24==-96)){s=2;break a}}else{if(!((u&-64)<<24>>24==-128)){s=2;break a}}v=t&255;if((v&192|0)!=128){s=2;break a}t=(u&255)<<6&4032|r<<12&61440|v&63;if(t>>>0>k>>>0){s=2;break a}c[l>>2]=t;c[f>>2]=b+3;break}if(!((q&255)<245)){s=2;break a}if((g-b|0)<4){s=1;break a}t=a[b+1|0]|0;v=a[b+2|0]|0;u=a[b+3|0]|0;if((r|0)==240){if(!((t+112<<24>>24&255)<48)){s=2;break a}}else if((r|0)==244){if(!((t&-16)<<24>>24==-128)){s=2;break a}}else{if(!((t&-64)<<24>>24==-128)){s=2;break a}}w=v&255;if((w&192|0)!=128){s=2;break a}v=u&255;if((v&192|0)!=128){s=2;break a}u=(t&255)<<12&258048|r<<18&1835008|w<<6&4032|v&63;if(u>>>0>k>>>0){s=2;break a}c[l>>2]=u;c[f>>2]=b+4}}while(0);r=(c[j>>2]|0)+4|0;c[j>>2]=r;q=c[f>>2]|0;if(q>>>0<e>>>0){l=r;b=q}else{o=q;p=39;break}}}else{o=n;p=39}}while(0);if((p|0)==39){s=o>>>0<e>>>0|0}i=m;return s|0}function el(a,b,d,e,f){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;c[f>>2]=d;return 3}function fl(a){a=a|0;return 0}function gl(a){a=a|0;return 0}function hl(a,b,c,d,e){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;b=i;a=il(c,d,e,1114111,0)|0;i=b;return a|0}function il(b,c,e,f,g){b=b|0;c=c|0;e=e|0;f=f|0;g=g|0;var h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0;h=i;if((((g&4|0)!=0?(c-b|0)>2:0)?(a[b]|0)==-17:0)?(a[b+1|0]|0)==-69:0){j=(a[b+2|0]|0)==-65?b+3|0:b}else{j=b}a:do{if(j>>>0<c>>>0&(e|0)!=0){g=c;k=j;l=1;b:while(1){m=a[k]|0;n=m&255;do{if(m<<24>>24>-1){if(n>>>0>f>>>0){o=k;break a}p=k+1|0}else{if((m&255)<194){o=k;break a}if((m&255)<224){if((g-k|0)<2){o=k;break a}q=d[k+1|0]|0;if((q&192|0)!=128){o=k;break a}if((q&63|n<<6&1984)>>>0>f>>>0){o=k;break a}p=k+2|0;break}if((m&255)<240){r=k;if((g-r|0)<3){o=k;break a}q=a[k+1|0]|0;s=a[k+2|0]|0;if((n|0)==237){if(!((q&-32)<<24>>24==-128)){t=23;break b}}else if((n|0)==224){if(!((q&-32)<<24>>24==-96)){t=21;break b}}else{if(!((q&-64)<<24>>24==-128)){t=25;break b}}u=s&255;if((u&192|0)!=128){o=k;break a}if(((q&255)<<6&4032|n<<12&61440|u&63)>>>0>f>>>0){o=k;break a}p=k+3|0;break}if(!((m&255)<245)){o=k;break a}v=k;if((g-v|0)<4){o=k;break a}u=a[k+1|0]|0;q=a[k+2|0]|0;s=a[k+3|0]|0;if((n|0)==240){if(!((u+112<<24>>24&255)<48)){t=33;break b}}else if((n|0)==244){if(!((u&-16)<<24>>24==-128)){t=35;break b}}else{if(!((u&-64)<<24>>24==-128)){t=37;break b}}w=q&255;if((w&192|0)!=128){o=k;break a}q=s&255;if((q&192|0)!=128){o=k;break a}if(((u&255)<<12&258048|n<<18&1835008|w<<6&4032|q&63)>>>0>f>>>0){o=k;break a}p=k+4|0}}while(0);if(!(p>>>0<c>>>0&l>>>0<e>>>0)){o=p;break a}k=p;l=l+1|0}if((t|0)==21){x=r-b|0;i=h;return x|0}else if((t|0)==23){x=r-b|0;i=h;return x|0}else if((t|0)==25){x=r-b|0;i=h;return x|0}else if((t|0)==33){x=v-b|0;i=h;return x|0}else if((t|0)==35){x=v-b|0;i=h;return x|0}else if((t|0)==37){x=v-b|0;i=h;return x|0}}else{o=j}}while(0);x=o-b|0;i=h;return x|0}function jl(a){a=a|0;return 4}function kl(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function ll(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function ml(a){a=a|0;var b=0;b=i;c[a>>2]=8424;ze(a+12|0);Um(a);i=b;return}function nl(a){a=a|0;var b=0;b=i;c[a>>2]=8424;ze(a+12|0);i=b;return}function ol(a){a=a|0;var b=0;b=i;c[a>>2]=8464;ze(a+16|0);Um(a);i=b;return}function pl(a){a=a|0;var b=0;b=i;c[a>>2]=8464;ze(a+16|0);i=b;return}function ql(b){b=b|0;return a[b+8|0]|0}function rl(a){a=a|0;return c[a+8>>2]|0}function sl(b){b=b|0;return a[b+9|0]|0}function tl(a){a=a|0;return c[a+12>>2]|0}function ul(a,b){a=a|0;b=b|0;var c=0;c=i;we(a,b+12|0);i=c;return}function vl(a,b){a=a|0;b=b|0;var c=0;c=i;we(a,b+16|0);i=c;return}function wl(a,b){a=a|0;b=b|0;b=i;xe(a,8496,4);i=b;return}function xl(a,b){a=a|0;b=b|0;b=i;Ie(a,8504,lm(8504)|0);i=b;return}function yl(a,b){a=a|0;b=b|0;b=i;xe(a,8528,5);i=b;return}function zl(a,b){a=a|0;b=b|0;b=i;Ie(a,8536,lm(8536)|0);i=b;return}function Al(b){b=b|0;var d=0;b=i;if((a[8568]|0)!=0){d=c[2140]|0;i=b;return d|0}if((Oa(8568)|0)==0){d=c[2140]|0;i=b;return d|0}if((a[15768]|0)==0?(Oa(15768)|0)!=0:0){pn(15600,0,168)|0;Ra(116,0,p|0)|0;$a(15768)}Ae(15600,15776)|0;Ae(15612|0,15784)|0;Ae(15624|0,15792)|0;Ae(15636|0,15800)|0;Ae(15648|0,15816)|0;Ae(15660|0,15832)|0;Ae(15672|0,15840)|0;Ae(15684|0,15856)|0;Ae(15696|0,15864)|0;Ae(15708|0,15872)|0;Ae(15720|0,15880)|0;Ae(15732|0,15888)|0;Ae(15744|0,15896)|0;Ae(15756|0,15904)|0;c[2140]=15600;$a(8568);d=c[2140]|0;i=b;return d|0}function Bl(b){b=b|0;var d=0;b=i;if((a[8584]|0)!=0){d=c[2144]|0;i=b;return d|0}if((Oa(8584)|0)==0){d=c[2144]|0;i=b;return d|0}if((a[15232]|0)==0?(Oa(15232)|0)!=0:0){pn(15064,0,168)|0;Ra(117,0,p|0)|0;$a(15232)}Le(15064,15240)|0;Le(15076|0,15272)|0;Le(15088|0,15304)|0;Le(15100|0,15336)|0;Le(15112|0,15376)|0;Le(15124|0,15416)|0;Le(15136|0,15448)|0;Le(15148|0,15488)|0;Le(15160|0,15504)|0;Le(15172|0,15520)|0;Le(15184|0,15536)|0;Le(15196|0,15552)|0;Le(15208|0,15568)|0;Le(15220|0,15584)|0;c[2144]=15064;$a(8584);d=c[2144]|0;i=b;return d|0}function Cl(b){b=b|0;var d=0;b=i;if((a[8600]|0)!=0){d=c[2148]|0;i=b;return d|0}if((Oa(8600)|0)==0){d=c[2148]|0;i=b;return d|0}if((a[14840]|0)==0?(Oa(14840)|0)!=0:0){pn(14552,0,288)|0;Ra(118,0,p|0)|0;$a(14840)}Ae(14552,14848)|0;Ae(14564|0,14856)|0;Ae(14576|0,14872)|0;Ae(14588|0,14880)|0;Ae(14600|0,14888)|0;Ae(14612|0,14896)|0;Ae(14624|0,14904)|0;Ae(14636|0,14912)|0;Ae(14648|0,14920)|0;Ae(14660|0,14936)|0;Ae(14672|0,14944)|0;Ae(14684|0,14960)|0;Ae(14696|0,14976)|0;Ae(14708|0,14984)|0;Ae(14720|0,14992)|0;Ae(14732|0,15e3)|0;Ae(14744|0,14888)|0;Ae(14756|0,15008)|0;Ae(14768|0,15016)|0;Ae(14780|0,15024)|0;Ae(14792|0,15032)|0;Ae(14804|0,15040)|0;Ae(14816|0,15048)|0;Ae(14828|0,15056)|0;c[2148]=14552;$a(8600);d=c[2148]|0;i=b;return d|0}function Dl(b){b=b|0;var d=0;b=i;if((a[8616]|0)!=0){d=c[2152]|0;i=b;return d|0}if((Oa(8616)|0)==0){d=c[2152]|0;i=b;return d|0}if((a[14e3]|0)==0?(Oa(14e3)|0)!=0:0){pn(13712,0,288)|0;Ra(119,0,p|0)|0;$a(14e3)}Le(13712,14008)|0;Le(13724|0,14040)|0;Le(13736|0,14080)|0;Le(13748|0,14104)|0;Le(13760|0,14424)|0;Le(13772|0,14128)|0;Le(13784|0,14152)|0;Le(13796|0,14176)|0;Le(13808|0,14208)|0;Le(13820|0,14248)|0;Le(13832|0,14280)|0;Le(13844|0,14320)|0;Le(13856|0,14360)|0;Le(13868|0,14376)|0;Le(13880|0,14392)|0;Le(13892|0,14408)|0;Le(13904|0,14424)|0;Le(13916|0,14440)|0;Le(13928|0,14456)|0;Le(13940|0,14472)|0;Le(13952|0,14488)|0;Le(13964|0,14504)|0;Le(13976|0,14520)|0;Le(13988|0,14536)|0;c[2152]=13712;$a(8616);d=c[2152]|0;i=b;return d|0}function El(b){b=b|0;var d=0;b=i;if((a[8632]|0)!=0){d=c[2156]|0;i=b;return d|0}if((Oa(8632)|0)==0){d=c[2156]|0;i=b;return d|0}if((a[13688]|0)==0?(Oa(13688)|0)!=0:0){pn(13400,0,288)|0;Ra(120,0,p|0)|0;$a(13688)}Ae(13400,13696)|0;Ae(13412|0,13704)|0;c[2156]=13400;$a(8632);d=c[2156]|0;i=b;return d|0}function Fl(b){b=b|0;var d=0;b=i;if((a[8648]|0)!=0){d=c[2160]|0;i=b;return d|0}if((Oa(8648)|0)==0){d=c[2160]|0;i=b;return d|0}if((a[13360]|0)==0?(Oa(13360)|0)!=0:0){pn(13072,0,288)|0;Ra(121,0,p|0)|0;$a(13360)}Le(13072,13368)|0;Le(13084|0,13384)|0;c[2160]=13072;$a(8648);d=c[2160]|0;i=b;return d|0}function Gl(b){b=b|0;b=i;if((a[8672]|0)==0?(Oa(8672)|0)!=0:0){xe(8656,8680,8);Ra(122,8656,p|0)|0;$a(8672)}i=b;return 8656}function Hl(b){b=b|0;b=i;if((a[8712]|0)!=0){i=b;return 8696}if((Oa(8712)|0)==0){i=b;return 8696}Ie(8696,8720,lm(8720)|0);Ra(123,8696,p|0)|0;$a(8712);i=b;return 8696}function Il(b){b=b|0;b=i;if((a[8776]|0)==0?(Oa(8776)|0)!=0:0){xe(8760,8784,8);Ra(122,8760,p|0)|0;$a(8776)}i=b;return 8760}function Jl(b){b=b|0;b=i;if((a[8816]|0)!=0){i=b;return 8800}if((Oa(8816)|0)==0){i=b;return 8800}Ie(8800,8824,lm(8824)|0);Ra(123,8800,p|0)|0;$a(8816);i=b;return 8800}function Kl(b){b=b|0;b=i;if((a[8880]|0)==0?(Oa(8880)|0)!=0:0){xe(8864,8888,20);Ra(122,8864,p|0)|0;$a(8880)}i=b;return 8864}function Ll(b){b=b|0;b=i;if((a[8928]|0)!=0){i=b;return 8912}if((Oa(8928)|0)==0){i=b;return 8912}Ie(8912,8936,lm(8936)|0);Ra(123,8912,p|0)|0;$a(8928);i=b;return 8912}function Ml(b){b=b|0;b=i;if((a[9040]|0)==0?(Oa(9040)|0)!=0:0){xe(9024,9048,11);Ra(122,9024,p|0)|0;$a(9040)}i=b;return 9024}function Nl(b){b=b|0;b=i;if((a[9080]|0)!=0){i=b;return 9064}if((Oa(9080)|0)==0){i=b;return 9064}Ie(9064,9088,lm(9088)|0);Ra(123,9064,p|0)|0;$a(9080);i=b;return 9064}function Ol(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0.0,j=0,k=0,l=0.0;f=i;i=i+16|0;g=f;if((b|0)==(d|0)){c[e>>2]=4;h=0.0;i=f;return+h}j=bb()|0;k=c[j>>2]|0;c[j>>2]=0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}l=+hn(b,g,c[2036]|0);b=c[j>>2]|0;if((b|0)==0){c[j>>2]=k}if((c[g>>2]|0)!=(d|0)){c[e>>2]=4;h=0.0;i=f;return+h}if((b|0)!=34){h=l;i=f;return+h}c[e>>2]=4;h=l;i=f;return+h}function Pl(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0.0,j=0,k=0,l=0.0;f=i;i=i+16|0;g=f;if((b|0)==(d|0)){c[e>>2]=4;h=0.0;i=f;return+h}j=bb()|0;k=c[j>>2]|0;c[j>>2]=0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}l=+hn(b,g,c[2036]|0);b=c[j>>2]|0;if((b|0)==0){c[j>>2]=k}if((c[g>>2]|0)!=(d|0)){c[e>>2]=4;h=0.0;i=f;return+h}if((b|0)!=34){h=l;i=f;return+h}c[e>>2]=4;h=l;i=f;return+h}function Ql(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0.0,j=0,k=0,l=0.0;f=i;i=i+16|0;g=f;if((b|0)==(d|0)){c[e>>2]=4;h=0.0;i=f;return+h}j=bb()|0;k=c[j>>2]|0;c[j>>2]=0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}l=+hn(b,g,c[2036]|0);b=c[j>>2]|0;if((b|0)==0){c[j>>2]=k}if((c[g>>2]|0)!=(d|0)){c[e>>2]=4;h=0.0;i=f;return+h}if((b|0)==34){c[e>>2]=4}h=l;i=f;return+h}function Rl(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0;g=i;i=i+16|0;h=g;do{if((b|0)!=(d|0)){if((a[b]|0)==45){c[e>>2]=4;j=0;k=0;break}l=bb()|0;m=c[l>>2]|0;c[l>>2]=0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}n=Ob(b|0,h|0,f|0,c[2036]|0)|0;o=c[l>>2]|0;if((o|0)==0){c[l>>2]=m}if((c[h>>2]|0)!=(d|0)){c[e>>2]=4;j=0;k=0;break}if((o|0)==34){c[e>>2]=4;j=-1;k=-1}else{j=I;k=n}}else{c[e>>2]=4;j=0;k=0}}while(0);I=j;i=g;return k|0}function Sl(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0;g=i;i=i+16|0;h=g;if((b|0)==(d|0)){c[e>>2]=4;j=0;i=g;return j|0}if((a[b]|0)==45){c[e>>2]=4;j=0;i=g;return j|0}k=bb()|0;l=c[k>>2]|0;c[k>>2]=0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}m=Ob(b|0,h|0,f|0,c[2036]|0)|0;f=I;b=c[k>>2]|0;if((b|0)==0){c[k>>2]=l}if((c[h>>2]|0)!=(d|0)){c[e>>2]=4;j=0;i=g;return j|0}if((b|0)==34|(f>>>0>0|(f|0)==0&m>>>0>4294967295)){c[e>>2]=4;j=-1;i=g;return j|0}else{j=m;i=g;return j|0}return 0}function Tl(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0;g=i;i=i+16|0;h=g;if((b|0)==(d|0)){c[e>>2]=4;j=0;i=g;return j|0}if((a[b]|0)==45){c[e>>2]=4;j=0;i=g;return j|0}k=bb()|0;l=c[k>>2]|0;c[k>>2]=0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}m=Ob(b|0,h|0,f|0,c[2036]|0)|0;f=I;b=c[k>>2]|0;if((b|0)==0){c[k>>2]=l}if((c[h>>2]|0)!=(d|0)){c[e>>2]=4;j=0;i=g;return j|0}if((b|0)==34|(f>>>0>0|(f|0)==0&m>>>0>4294967295)){c[e>>2]=4;j=-1;i=g;return j|0}else{j=m;i=g;return j|0}return 0}function Ul(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0;g=i;i=i+16|0;h=g;if((b|0)==(d|0)){c[e>>2]=4;j=0;i=g;return j|0}if((a[b]|0)==45){c[e>>2]=4;j=0;i=g;return j|0}k=bb()|0;l=c[k>>2]|0;c[k>>2]=0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}m=Ob(b|0,h|0,f|0,c[2036]|0)|0;f=I;b=c[k>>2]|0;if((b|0)==0){c[k>>2]=l}if((c[h>>2]|0)!=(d|0)){c[e>>2]=4;j=0;i=g;return j|0}if((b|0)==34|(f>>>0>0|(f|0)==0&m>>>0>65535)){c[e>>2]=4;j=-1;i=g;return j|0}else{j=m&65535;i=g;return j|0}return 0}function Vl(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0;g=i;i=i+16|0;h=g;if((b|0)==(d|0)){c[e>>2]=4;j=0;k=0;I=j;i=g;return k|0}l=bb()|0;m=c[l>>2]|0;c[l>>2]=0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}n=Ja(b|0,h|0,f|0,c[2036]|0)|0;f=I;b=c[l>>2]|0;if((b|0)==0){c[l>>2]=m}if((c[h>>2]|0)!=(d|0)){c[e>>2]=4;j=0;k=0;I=j;i=g;return k|0}if((b|0)==34){c[e>>2]=4;e=(f|0)>0|(f|0)==0&n>>>0>0;I=e?2147483647:-2147483648;i=g;return(e?-1:0)|0}else{j=f;k=n;I=j;i=g;return k|0}return 0}function Wl(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0;g=i;i=i+16|0;h=g;if((b|0)==(d|0)){c[e>>2]=4;j=0;i=g;return j|0}k=bb()|0;l=c[k>>2]|0;c[k>>2]=0;if((a[8152]|0)==0?(Oa(8152)|0)!=0:0){c[2036]=lb(2147483647,8160,0)|0;$a(8152)}m=Ja(b|0,h|0,f|0,c[2036]|0)|0;f=I;b=c[k>>2]|0;if((b|0)==0){c[k>>2]=l}if((c[h>>2]|0)!=(d|0)){c[e>>2]=4;j=0;i=g;return j|0}do{if((b|0)==34){c[e>>2]=4;if((f|0)>0|(f|0)==0&m>>>0>0){j=2147483647;i=g;return j|0}}else{if((f|0)<-1|(f|0)==-1&m>>>0<2147483648){c[e>>2]=4;break}if((f|0)>0|(f|0)==0&m>>>0>2147483647){c[e>>2]=4;j=2147483647;i=g;return j|0}else{j=m;i=g;return j|0}}}while(0);j=-2147483648;i=g;return j|0}function Xl(a){a=a|0;var b=0,e=0,f=0,g=0,h=0;b=i;e=a+4|0;f=d[e]|d[e+1|0]<<8|d[e+2|0]<<16|d[e+3|0]<<24;g=e+4|0;e=d[g]|d[g+1|0]<<8|d[g+2|0]<<16|d[g+3|0]<<24;g=(c[a>>2]|0)+(e>>1)|0;if((e&1|0)==0){h=f;jc[h&127](g);i=b;return}else{h=c[(c[g>>2]|0)+f>>2]|0;jc[h&127](g);i=b;return}}function Yl(b,d){b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0;e=i;f=b+8|0;g=b+4|0;h=c[g>>2]|0;j=c[f>>2]|0;k=h;if(!(j-k>>2>>>0<d>>>0)){l=d;m=h;do{if((m|0)==0){n=0}else{c[m>>2]=0;n=c[g>>2]|0}m=n+4|0;c[g>>2]=m;l=l+ -1|0}while((l|0)!=0);i=e;return}l=b+16|0;m=c[b>>2]|0;n=k-m>>2;k=n+d|0;if(k>>>0>1073741823){Sj(0)}h=j-m|0;if(h>>2>>>0<536870911){m=h>>1;h=m>>>0<k>>>0?k:m;if((h|0)!=0){m=b+128|0;if((a[m]|0)==0&h>>>0<29){a[m]=1;o=h;p=l}else{q=h;r=11}}else{o=0;p=0}}else{q=1073741823;r=11}if((r|0)==11){o=q;p=Sm(q<<2)|0}q=d;d=p+(n<<2)|0;do{if((d|0)==0){s=0}else{c[d>>2]=0;s=d}d=s+4|0;q=q+ -1|0}while((q|0)!=0);q=c[b>>2]|0;s=(c[g>>2]|0)-q|0;r=p+(n-(s>>2)<<2)|0;nn(r|0,q|0,s|0)|0;c[b>>2]=r;c[g>>2]=d;c[f>>2]=p+(o<<2);if((q|0)==0){i=e;return}if((l|0)==(q|0)){a[b+128|0]=0;i=e;return}else{Um(q);i=e;return}}function Zl(a){a=a|0;a=i;Ke(13348|0);Ke(13336|0);Ke(13324|0);Ke(13312|0);Ke(13300|0);Ke(13288|0);Ke(13276|0);Ke(13264|0);Ke(13252|0);Ke(13240|0);Ke(13228|0);Ke(13216|0);Ke(13204|0);Ke(13192|0);Ke(13180|0);Ke(13168|0);Ke(13156|0);Ke(13144|0);Ke(13132|0);Ke(13120|0);Ke(13108|0);Ke(13096|0);Ke(13084|0);Ke(13072);i=a;return}function _l(a){a=a|0;a=i;ze(13676|0);ze(13664|0);ze(13652|0);ze(13640|0);ze(13628|0);ze(13616|0);ze(13604|0);ze(13592|0);ze(13580|0);ze(13568|0);ze(13556|0);ze(13544|0);ze(13532|0);ze(13520|0);ze(13508|0);ze(13496|0);ze(13484|0);ze(13472|0);ze(13460|0);ze(13448|0);ze(13436|0);ze(13424|0);ze(13412|0);ze(13400);i=a;return}function $l(a){a=a|0;a=i;Ke(13988|0);Ke(13976|0);Ke(13964|0);Ke(13952|0);Ke(13940|0);Ke(13928|0);Ke(13916|0);Ke(13904|0);Ke(13892|0);Ke(13880|0);Ke(13868|0);Ke(13856|0);Ke(13844|0);Ke(13832|0);Ke(13820|0);Ke(13808|0);Ke(13796|0);Ke(13784|0);Ke(13772|0);Ke(13760|0);Ke(13748|0);Ke(13736|0);Ke(13724|0);Ke(13712);i=a;return}function am(a){a=a|0;a=i;ze(14828|0);ze(14816|0);ze(14804|0);ze(14792|0);ze(14780|0);ze(14768|0);ze(14756|0);ze(14744|0);ze(14732|0);ze(14720|0);ze(14708|0);ze(14696|0);ze(14684|0);ze(14672|0);ze(14660|0);ze(14648|0);ze(14636|0);ze(14624|0);ze(14612|0);ze(14600|0);ze(14588|0);ze(14576|0);ze(14564|0);ze(14552);i=a;return}function bm(a){a=a|0;a=i;Ke(15220|0);Ke(15208|0);Ke(15196|0);Ke(15184|0);Ke(15172|0);Ke(15160|0);Ke(15148|0);Ke(15136|0);Ke(15124|0);Ke(15112|0);Ke(15100|0);Ke(15088|0);Ke(15076|0);Ke(15064);i=a;return}function cm(a){a=a|0;a=i;ze(15756|0);ze(15744|0);ze(15732|0);ze(15720|0);ze(15708|0);ze(15696|0);ze(15684|0);ze(15672|0);ze(15660|0);ze(15648|0);ze(15636|0);ze(15624|0);ze(15612|0);ze(15600);i=a;return}function dm(a,b,c){a=a|0;b=b|0;c=c|0;var d=0,e=0;d=i;e=em(0,a,b,(c|0)!=0?c:16120)|0;i=d;return e|0}function em(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0;g=i;i=i+16|0;h=g;c[h>>2]=b;j=(f|0)==0?16128:f;f=c[j>>2]|0;a:do{if((d|0)==0){if((f|0)==0){k=0;i=g;return k|0}}else{if((b|0)==0){c[h>>2]=h;l=h}else{l=b}if((e|0)==0){k=-2;i=g;return k|0}do{if((f|0)==0){m=a[d]|0;n=m&255;if(m<<24>>24>-1){c[l>>2]=n;k=m<<24>>24!=0|0;i=g;return k|0}else{m=n+ -194|0;if(m>>>0>50){break a}o=e+ -1|0;p=c[15912+(m<<2)>>2]|0;q=d+1|0;break}}else{o=e;p=f;q=d}}while(0);b:do{if((o|0)==0){r=p}else{m=a[q]|0;n=(m&255)>>>3;if((n+ -16|n+(p>>26))>>>0>7){break a}else{s=o;t=m;u=p;v=q}while(1){v=v+1|0;u=(t&255)+ -128|u<<6;s=s+ -1|0;if((u|0)>=0){break}if((s|0)==0){r=u;break b}t=a[v]|0;if(((t&255)+ -128|0)>>>0>63){break a}}c[j>>2]=0;c[l>>2]=u;k=e-s|0;i=g;return k|0}}while(0);c[j>>2]=r;k=-2;i=g;return k|0}}while(0);c[j>>2]=0;c[(bb()|0)>>2]=84;k=-1;i=g;return k|0}function fm(a,b,d,e,f){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0;g=i;i=i+1040|0;h=g+8|0;j=g;k=c[b>>2]|0;c[j>>2]=k;l=(a|0)!=0;m=l?e:256;e=l?a:h;a:do{if((k|0)==0|(m|0)==0){n=d;o=m;p=k;q=0;r=e}else{a=d;s=m;t=k;u=0;v=e;while(1){w=a>>>2;x=w>>>0>=s>>>0;if(!(x|a>>>0>131)){n=a;o=s;p=t;q=u;r=v;break a}y=x?s:w;z=a-y|0;w=gm(v,j,y,f)|0;if((w|0)==-1){break}if((v|0)==(h|0)){A=s;B=h}else{A=s-w|0;B=v+(w<<2)|0}y=w+u|0;w=c[j>>2]|0;if((w|0)==0|(A|0)==0){n=z;o=A;p=w;q=y;r=B;break a}else{a=z;s=A;t=w;u=y;v=B}}n=z;o=0;p=c[j>>2]|0;q=-1;r=v}}while(0);b:do{if((p|0)!=0?!((o|0)==0|(n|0)==0):0){z=n;B=o;A=p;h=q;e=r;while(1){C=em(e,A,z,f)|0;if((C+2|0)>>>0<3){break}k=(c[j>>2]|0)+C|0;c[j>>2]=k;m=B+ -1|0;d=h+1|0;if((m|0)==0|(z|0)==(C|0)){D=d;break b}else{z=z-C|0;B=m;A=k;h=d;e=e+4|0}}if((C|0)==-1){D=-1;break}else if((C|0)==0){c[j>>2]=0;D=h;break}else{c[f>>2]=0;D=h;break}}else{D=q}}while(0);if(!l){i=g;return D|0}c[b>>2]=c[j>>2];i=g;return D|0}function gm(b,e,f,g){b=b|0;e=e|0;f=f|0;g=g|0;var h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0,X=0,Y=0,Z=0,_=0,$=0,aa=0,ba=0,ca=0,da=0,ea=0,fa=0,ga=0,ha=0,ia=0,ja=0;h=i;j=c[e>>2]|0;if((g|0)!=0?(k=c[g>>2]|0,(k|0)!=0):0){if((b|0)==0){l=f;m=k;n=j;o=16}else{c[g>>2]=0;p=b;q=f;r=k;s=j;o=36}}else{if((b|0)==0){t=f;u=j;o=7}else{v=b;w=f;x=j;o=6}}a:while(1){if((o|0)==6){o=0;if((w|0)==0){y=f;o=53;break}else{z=v;A=w;B=x}while(1){j=a[B]|0;do{if(((j&255)+ -1|0)>>>0<127?(B&3|0)==0&A>>>0>3:0){k=z;g=A;C=B;while(1){D=c[C>>2]|0;if(((D+ -16843009|D)&-2139062144|0)!=0){o=30;break}c[k>>2]=D&255;c[k+4>>2]=d[C+1|0]|0;c[k+8>>2]=d[C+2|0]|0;E=C+4|0;F=k+16|0;c[k+12>>2]=d[C+3|0]|0;G=g+ -4|0;if(G>>>0>3){k=F;g=G;C=E}else{o=31;break}}if((o|0)==30){o=0;H=k;I=g;J=D&255;K=C;break}else if((o|0)==31){o=0;H=F;I=G;J=a[E]|0;K=E;break}}else{H=z;I=A;J=j;K=B}}while(0);L=J&255;if(!((L+ -1|0)>>>0<127)){break}c[H>>2]=L;j=I+ -1|0;if((j|0)==0){y=f;o=53;break a}else{z=H+4|0;A=j;B=K+1|0}}j=L+ -194|0;if(j>>>0>50){M=H;N=I;O=K;o=47;break}p=H;q=I;r=c[15912+(j<<2)>>2]|0;s=K+1|0;o=36;continue}else if((o|0)==7){o=0;j=a[u]|0;if(((j&255)+ -1|0)>>>0<127?(u&3|0)==0:0){P=c[u>>2]|0;if(((P+ -16843009|P)&-2139062144|0)==0){Q=t;R=u;while(1){S=R+4|0;T=Q+ -4|0;U=c[S>>2]|0;if(((U+ -16843009|U)&-2139062144|0)==0){Q=T;R=S}else{V=T;W=U;X=S;break}}}else{V=t;W=P;X=u}Y=V;Z=W&255;_=X}else{Y=t;Z=j;_=u}R=Z&255;if((R+ -1|0)>>>0<127){t=Y+ -1|0;u=_+1|0;o=7;continue}Q=R+ -194|0;if(Q>>>0>50){M=b;N=Y;O=_;o=47;break}l=Y;m=c[15912+(Q<<2)>>2]|0;n=_+1|0;o=16;continue}else if((o|0)==16){o=0;Q=(d[n]|0)>>>3;if((Q+ -16|Q+(m>>26))>>>0>7){o=17;break}Q=n+1|0;if((m&33554432|0)!=0){if(((d[Q]|0)+ -128|0)>>>0>63){o=20;break}R=n+2|0;if((m&524288|0)==0){$=R}else{if(((d[R]|0)+ -128|0)>>>0>63){o=23;break}$=n+3|0}}else{$=Q}t=l+ -1|0;u=$;o=7;continue}else if((o|0)==36){o=0;Q=d[s]|0;R=Q>>>3;if((R+ -16|R+(r>>26))>>>0>7){o=37;break}R=s+1|0;aa=Q+ -128|r<<6;if((aa|0)<0){Q=(d[R]|0)+ -128|0;if(Q>>>0>63){o=40;break}S=s+2|0;ba=Q|aa<<6;if((ba|0)<0){Q=(d[S]|0)+ -128|0;if(Q>>>0>63){o=43;break}ca=Q|ba<<6;da=s+3|0}else{ca=ba;da=S}}else{ca=aa;da=R}c[p>>2]=ca;v=p+4|0;w=q+ -1|0;x=da;o=6;continue}}if((o|0)==17){ea=b;fa=l;ga=m;ha=n+ -1|0;o=46}else if((o|0)==20){ea=b;fa=l;ga=m;ha=n+ -1|0;o=46}else if((o|0)==23){ea=b;fa=l;ga=m;ha=n+ -1|0;o=46}else if((o|0)==37){ea=p;fa=q;ga=r;ha=s+ -1|0;o=46}else if((o|0)==40){ea=p;fa=q;ga=aa;ha=s+ -1|0;o=46}else if((o|0)==43){ea=p;fa=q;ga=ba;ha=s+ -1|0;o=46}else if((o|0)==53){i=h;return y|0}if((o|0)==46){if((ga|0)==0){M=ea;N=fa;O=ha;o=47}else{ia=ea;ja=ha}}if((o|0)==47){if((a[O]|0)==0){if((M|0)!=0){c[M>>2]=0;c[e>>2]=0}y=f-N|0;i=h;return y|0}else{ia=M;ja=O}}c[(bb()|0)>>2]=84;if((ia|0)==0){y=-1;i=h;return y|0}c[e>>2]=ja;y=-1;i=h;return y|0}function hm(b,e,f){b=b|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0;g=i;i=i+16|0;h=g;c[h>>2]=b;if((e|0)==0){j=0;i=g;return j|0}do{if((f|0)!=0){if((b|0)==0){c[h>>2]=h;k=h}else{k=b}l=a[e]|0;m=l&255;if(l<<24>>24>-1){c[k>>2]=m;j=l<<24>>24!=0|0;i=g;return j|0}l=m+ -194|0;if(!(l>>>0>50)){m=e+1|0;n=c[15912+(l<<2)>>2]|0;if(f>>>0<4?(n&-2147483648>>>((f*6|0)+ -6|0)|0)!=0:0){break}l=d[m]|0;m=l>>>3;if(!((m+ -16|m+(n>>26))>>>0>7)){m=l+ -128|n<<6;if((m|0)>=0){c[k>>2]=m;j=2;i=g;return j|0}n=(d[e+2|0]|0)+ -128|0;if(!(n>>>0>63)){l=n|m<<6;if((l|0)>=0){c[k>>2]=l;j=3;i=g;return j|0}m=(d[e+3|0]|0)+ -128|0;if(!(m>>>0>63)){c[k>>2]=m|l<<6;j=4;i=g;return j|0}}}}}}while(0);c[(bb()|0)>>2]=84;j=-1;i=g;return j|0}function im(b,d,e){b=b|0;d=d|0;e=e|0;var f=0;e=i;if((b|0)==0){f=1;i=e;return f|0}if(d>>>0<128){a[b]=d;f=1;i=e;return f|0}if(d>>>0<2048){a[b]=d>>>6|192;a[b+1|0]=d&63|128;f=2;i=e;return f|0}if(d>>>0<55296|(d+ -57344|0)>>>0<8192){a[b]=d>>>12|224;a[b+1|0]=d>>>6&63|128;a[b+2|0]=d&63|128;f=3;i=e;return f|0}if((d+ -65536|0)>>>0<1048576){a[b]=d>>>18|240;a[b+1|0]=d>>>12&63|128;a[b+2|0]=d>>>6&63|128;a[b+3|0]=d&63|128;f=4;i=e;return f|0}else{c[(bb()|0)>>2]=84;f=-1;i=e;return f|0}return 0}function jm(a,b,d,e,f){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0;f=i;i=i+272|0;g=f+8|0;h=f;j=c[b>>2]|0;c[h>>2]=j;k=(a|0)!=0;l=k?e:256;e=k?a:g;a:do{if((j|0)==0|(l|0)==0){m=d;n=l;o=j;p=0;q=e}else{a=d;r=l;s=j;t=0;u=e;while(1){v=a>>>0>=r>>>0;if(!(v|a>>>0>32)){m=a;n=r;o=s;p=t;q=u;break a}w=v?r:a;x=a-w|0;v=km(u,h,w,0)|0;if((v|0)==-1){break}if((u|0)==(g|0)){y=r;z=g}else{y=r-v|0;z=u+v|0}w=v+t|0;v=c[h>>2]|0;if((v|0)==0|(y|0)==0){m=x;n=y;o=v;p=w;q=z;break a}else{a=x;r=y;s=v;t=w;u=z}}m=x;n=0;o=c[h>>2]|0;p=-1;q=u}}while(0);b:do{if((o|0)!=0?!((n|0)==0|(m|0)==0):0){x=m;z=n;y=o;g=p;e=q;while(1){A=im(e,c[y>>2]|0,0)|0;if((A+1|0)>>>0<2){break}j=(c[h>>2]|0)+4|0;c[h>>2]=j;l=x+ -1|0;d=g+1|0;if((z|0)==(A|0)|(l|0)==0){B=d;break b}else{x=l;z=z-A|0;y=j;g=d;e=e+A|0}}if((A|0)==0){c[h>>2]=0;B=g}else{B=-1}}else{B=p}}while(0);if(!k){i=f;return B|0}c[b>>2]=c[h>>2];i=f;return B|0}function km(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0;f=i;i=i+16|0;g=f;if((b|0)==0){h=c[d>>2]|0;j=c[h>>2]|0;if((j|0)==0){k=0;i=f;return k|0}else{l=0;m=j;n=h}while(1){if(m>>>0>127){h=im(g,m,0)|0;if((h|0)==-1){k=-1;o=26;break}else{p=h}}else{p=1}h=p+l|0;j=n+4|0;q=c[j>>2]|0;if((q|0)==0){k=h;o=26;break}else{l=h;m=q;n=j}}if((o|0)==26){i=f;return k|0}}a:do{if(e>>>0>3){n=b;m=e;l=c[d>>2]|0;while(1){p=c[l>>2]|0;if((p|0)==0){r=n;s=m;break a}if(p>>>0>127){j=im(n,p,0)|0;if((j|0)==-1){k=-1;break}t=n+j|0;u=m-j|0;v=l}else{a[n]=p;t=n+1|0;u=m+ -1|0;v=c[d>>2]|0}p=v+4|0;c[d>>2]=p;if(u>>>0>3){n=t;m=u;l=p}else{r=t;s=u;break a}}i=f;return k|0}else{r=b;s=e}}while(0);b:do{if((s|0)!=0){b=r;u=s;t=c[d>>2]|0;while(1){v=c[t>>2]|0;if((v|0)==0){o=24;break}if(v>>>0>127){l=im(g,v,0)|0;if((l|0)==-1){k=-1;o=26;break}if(l>>>0>u>>>0){o=20;break}im(b,c[t>>2]|0,0)|0;w=b+l|0;x=u-l|0;y=t}else{a[b]=v;w=b+1|0;x=u+ -1|0;y=c[d>>2]|0}v=y+4|0;c[d>>2]=v;if((x|0)==0){z=0;break b}else{b=w;u=x;t=v}}if((o|0)==20){k=e-u|0;i=f;return k|0}else if((o|0)==24){a[b]=0;z=u;break}else if((o|0)==26){i=f;return k|0}}else{z=0}}while(0);c[d>>2]=0;k=e-z|0;i=f;return k|0}function lm(a){a=a|0;var b=0,d=0;b=i;d=a;while(1){if((c[d>>2]|0)==0){break}else{d=d+4|0}}i=b;return d-a>>2|0}function mm(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,g=0,h=0;e=i;if((d|0)==0){i=e;return a|0}else{f=d;g=b;h=a}while(1){f=f+ -1|0;c[h>>2]=c[g>>2];if((f|0)==0){break}else{g=g+4|0;h=h+4|0}}i=e;return a|0}function nm(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,g=0;e=i;f=(d|0)==0;if(a-b>>2>>>0<d>>>0){if(!f){g=d;do{g=g+ -1|0;c[a+(g<<2)>>2]=c[b+(g<<2)>>2]}while((g|0)!=0)}}else{if(!f){f=b;b=a;g=d;while(1){g=g+ -1|0;c[b>>2]=c[f>>2];if((g|0)==0){break}else{f=f+4|0;b=b+4|0}}}}i=e;return a|0}function om(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0;e=i;if((d|0)!=0){f=d;d=a;while(1){f=f+ -1|0;c[d>>2]=b;if((f|0)==0){break}else{d=d+4|0}}}i=e;return a|0}function pm(a){a=a|0;return}function qm(a){a=a|0;c[a>>2]=16144;return}function rm(a){a=a|0;var b=0;b=i;Mb(a|0);Um(a);i=b;return}function sm(a){a=a|0;var b=0;b=i;Mb(a|0);i=b;return}function tm(a){a=a|0;return 16160}function um(a){a=a|0;return}function vm(a){a=a|0;return}function wm(a){a=a|0;return}function xm(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function ym(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function zm(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function Am(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0;e=i;i=i+64|0;f=e;if((a|0)==(b|0)){g=1;i=e;return g|0}if((b|0)==0){g=0;i=e;return g|0}h=Em(b,16272,16328,0)|0;if((h|0)==0){g=0;i=e;return g|0}b=f+0|0;j=b+56|0;do{c[b>>2]=0;b=b+4|0}while((b|0)<(j|0));c[f>>2]=h;c[f+8>>2]=a;c[f+12>>2]=-1;c[f+48>>2]=1;xc[c[(c[h>>2]|0)+28>>2]&7](h,f,c[d>>2]|0,1);if((c[f+24>>2]|0)!=1){g=0;i=e;return g|0}c[d>>2]=c[f+16>>2];g=1;i=e;return g|0}function Bm(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0;g=i;if((c[d+8>>2]|0)!=(b|0)){i=g;return}b=d+16|0;h=c[b>>2]|0;if((h|0)==0){c[b>>2]=e;c[d+24>>2]=f;c[d+36>>2]=1;i=g;return}if((h|0)!=(e|0)){e=d+36|0;c[e>>2]=(c[e>>2]|0)+1;c[d+24>>2]=2;a[d+54|0]=1;i=g;return}e=d+24|0;if((c[e>>2]|0)!=2){i=g;return}c[e>>2]=f;i=g;return}function Cm(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0;g=i;if((b|0)!=(c[d+8>>2]|0)){h=c[b+8>>2]|0;xc[c[(c[h>>2]|0)+28>>2]&7](h,d,e,f);i=g;return}h=d+16|0;b=c[h>>2]|0;if((b|0)==0){c[h>>2]=e;c[d+24>>2]=f;c[d+36>>2]=1;i=g;return}if((b|0)!=(e|0)){e=d+36|0;c[e>>2]=(c[e>>2]|0)+1;c[d+24>>2]=2;a[d+54|0]=1;i=g;return}e=d+24|0;if((c[e>>2]|0)!=2){i=g;return}c[e>>2]=f;i=g;return}function Dm(b,d,e,f){b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0;g=i;if((b|0)==(c[d+8>>2]|0)){h=d+16|0;j=c[h>>2]|0;if((j|0)==0){c[h>>2]=e;c[d+24>>2]=f;c[d+36>>2]=1;i=g;return}if((j|0)!=(e|0)){j=d+36|0;c[j>>2]=(c[j>>2]|0)+1;c[d+24>>2]=2;a[d+54|0]=1;i=g;return}j=d+24|0;if((c[j>>2]|0)!=2){i=g;return}c[j>>2]=f;i=g;return}j=c[b+12>>2]|0;h=b+(j<<3)+16|0;k=c[b+20>>2]|0;l=k>>8;if((k&1|0)==0){m=l}else{m=c[(c[e>>2]|0)+l>>2]|0}l=c[b+16>>2]|0;xc[c[(c[l>>2]|0)+28>>2]&7](l,d,e+m|0,(k&2|0)!=0?f:2);if((j|0)<=1){i=g;return}j=d+54|0;k=b+24|0;while(1){b=c[k+4>>2]|0;m=b>>8;if((b&1|0)==0){n=m}else{n=c[(c[e>>2]|0)+m>>2]|0}m=c[k>>2]|0;xc[c[(c[m>>2]|0)+28>>2]&7](m,d,e+n|0,(b&2|0)!=0?f:2);if((a[j]|0)!=0){o=16;break}b=k+8|0;if(b>>>0<h>>>0){k=b}else{o=16;break}}if((o|0)==16){i=g;return}}function Em(d,e,f,g){d=d|0;e=e|0;f=f|0;g=g|0;var h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0;h=i;i=i+64|0;j=h;k=c[d>>2]|0;l=d+(c[k+ -8>>2]|0)|0;m=c[k+ -4>>2]|0;c[j>>2]=f;c[j+4>>2]=d;c[j+8>>2]=e;c[j+12>>2]=g;g=j+16|0;e=j+20|0;d=j+24|0;k=j+28|0;n=j+32|0;o=j+40|0;p=(m|0)==(f|0);f=g+0|0;q=f+36|0;do{c[f>>2]=0;f=f+4|0}while((f|0)<(q|0));b[g+36>>1]=0;a[g+38|0]=0;if(p){c[j+48>>2]=1;uc[c[(c[m>>2]|0)+20>>2]&15](m,j,l,l,1,0);r=(c[d>>2]|0)==1?l:0;i=h;return r|0}ic[c[(c[m>>2]|0)+24>>2]&3](m,j,l,1,0);l=c[j+36>>2]|0;if((l|0)==1){if((c[d>>2]|0)!=1){if((c[o>>2]|0)!=0){r=0;i=h;return r|0}if((c[k>>2]|0)!=1){r=0;i=h;return r|0}if((c[n>>2]|0)!=1){r=0;i=h;return r|0}}r=c[g>>2]|0;i=h;return r|0}else if((l|0)==0){if((c[o>>2]|0)!=1){r=0;i=h;return r|0}if((c[k>>2]|0)!=1){r=0;i=h;return r|0}r=(c[n>>2]|0)==1?c[e>>2]|0:0;i=h;return r|0}else{r=0;i=h;return r|0}return 0}function Fm(b,d,e,f,g){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;var h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0;h=i;if((b|0)==(c[d+8>>2]|0)){if((c[d+4>>2]|0)!=(e|0)){i=h;return}j=d+28|0;if((c[j>>2]|0)==1){i=h;return}c[j>>2]=f;i=h;return}if((b|0)==(c[d>>2]|0)){if((c[d+16>>2]|0)!=(e|0)?(j=d+20|0,(c[j>>2]|0)!=(e|0)):0){c[d+32>>2]=f;k=d+44|0;if((c[k>>2]|0)==4){i=h;return}l=c[b+12>>2]|0;m=b+(l<<3)+16|0;a:do{if((l|0)>0){n=d+52|0;o=d+53|0;p=d+54|0;q=b+8|0;r=d+24|0;s=0;t=0;u=b+16|0;b:while(1){a[n]=0;a[o]=0;v=c[u+4>>2]|0;w=v>>8;if((v&1|0)==0){x=w}else{x=c[(c[e>>2]|0)+w>>2]|0}w=c[u>>2]|0;uc[c[(c[w>>2]|0)+20>>2]&15](w,d,e,e+x|0,2-(v>>>1&1)|0,g);if((a[p]|0)!=0){y=s;z=t;break}do{if((a[o]|0)!=0){if((a[n]|0)==0){if((c[q>>2]&1|0)==0){y=s;z=1;break b}else{A=s;B=1;break}}if((c[r>>2]|0)==1){C=27;break a}if((c[q>>2]&2|0)==0){C=27;break a}else{A=1;B=1}}else{A=s;B=t}}while(0);v=u+8|0;if(v>>>0<m>>>0){s=A;t=B;u=v}else{y=A;z=B;break}}if(y){D=z;C=26}else{E=z;C=23}}else{E=0;C=23}}while(0);if((C|0)==23){c[j>>2]=e;j=d+40|0;c[j>>2]=(c[j>>2]|0)+1;if((c[d+36>>2]|0)==1?(c[d+24>>2]|0)==2:0){a[d+54|0]=1;if(E){C=27}else{C=28}}else{D=E;C=26}}if((C|0)==26){if(D){C=27}else{C=28}}if((C|0)==27){c[k>>2]=3;i=h;return}else if((C|0)==28){c[k>>2]=4;i=h;return}}if((f|0)!=1){i=h;return}c[d+32>>2]=1;i=h;return}k=c[b+12>>2]|0;D=b+(k<<3)+16|0;E=c[b+20>>2]|0;j=E>>8;if((E&1|0)==0){F=j}else{F=c[(c[e>>2]|0)+j>>2]|0}j=c[b+16>>2]|0;ic[c[(c[j>>2]|0)+24>>2]&3](j,d,e+F|0,(E&2|0)!=0?f:2,g);E=b+24|0;if((k|0)<=1){i=h;return}k=c[b+8>>2]|0;if((k&2|0)==0?(b=d+36|0,(c[b>>2]|0)!=1):0){if((k&1|0)==0){k=d+54|0;F=E;while(1){if((a[k]|0)!=0){C=53;break}if((c[b>>2]|0)==1){C=53;break}j=c[F+4>>2]|0;z=j>>8;if((j&1|0)==0){G=z}else{G=c[(c[e>>2]|0)+z>>2]|0}z=c[F>>2]|0;ic[c[(c[z>>2]|0)+24>>2]&3](z,d,e+G|0,(j&2|0)!=0?f:2,g);j=F+8|0;if(j>>>0<D>>>0){F=j}else{C=53;break}}if((C|0)==53){i=h;return}}F=d+24|0;G=d+54|0;k=E;while(1){if((a[G]|0)!=0){C=53;break}if((c[b>>2]|0)==1?(c[F>>2]|0)==1:0){C=53;break}j=c[k+4>>2]|0;z=j>>8;if((j&1|0)==0){H=z}else{H=c[(c[e>>2]|0)+z>>2]|0}z=c[k>>2]|0;ic[c[(c[z>>2]|0)+24>>2]&3](z,d,e+H|0,(j&2|0)!=0?f:2,g);j=k+8|0;if(j>>>0<D>>>0){k=j}else{C=53;break}}if((C|0)==53){i=h;return}}k=d+54|0;H=E;while(1){if((a[k]|0)!=0){C=53;break}E=c[H+4>>2]|0;F=E>>8;if((E&1|0)==0){I=F}else{I=c[(c[e>>2]|0)+F>>2]|0}F=c[H>>2]|0;ic[c[(c[F>>2]|0)+24>>2]&3](F,d,e+I|0,(E&2|0)!=0?f:2,g);E=H+8|0;if(E>>>0<D>>>0){H=E}else{C=53;break}}if((C|0)==53){i=h;return}}function Gm(b,d,e,f,g){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;var h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0;h=i;if((b|0)==(c[d+8>>2]|0)){if((c[d+4>>2]|0)!=(e|0)){i=h;return}j=d+28|0;if((c[j>>2]|0)==1){i=h;return}c[j>>2]=f;i=h;return}if((b|0)!=(c[d>>2]|0)){j=c[b+8>>2]|0;ic[c[(c[j>>2]|0)+24>>2]&3](j,d,e,f,g);i=h;return}if((c[d+16>>2]|0)!=(e|0)?(j=d+20|0,(c[j>>2]|0)!=(e|0)):0){c[d+32>>2]=f;k=d+44|0;if((c[k>>2]|0)==4){i=h;return}l=d+52|0;a[l]=0;m=d+53|0;a[m]=0;n=c[b+8>>2]|0;uc[c[(c[n>>2]|0)+20>>2]&15](n,d,e,e,1,g);if((a[m]|0)!=0){if((a[l]|0)==0){o=1;p=13}}else{o=0;p=13}do{if((p|0)==13){c[j>>2]=e;l=d+40|0;c[l>>2]=(c[l>>2]|0)+1;if((c[d+36>>2]|0)==1?(c[d+24>>2]|0)==2:0){a[d+54|0]=1;if(o){break}}else{p=16}if((p|0)==16?o:0){break}c[k>>2]=4;i=h;return}}while(0);c[k>>2]=3;i=h;return}if((f|0)!=1){i=h;return}c[d+32>>2]=1;i=h;return}function Hm(b,d,e,f,g){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;var h=0;g=i;if((c[d+8>>2]|0)==(b|0)){if((c[d+4>>2]|0)!=(e|0)){i=g;return}h=d+28|0;if((c[h>>2]|0)==1){i=g;return}c[h>>2]=f;i=g;return}if((c[d>>2]|0)!=(b|0)){i=g;return}if((c[d+16>>2]|0)!=(e|0)?(b=d+20|0,(c[b>>2]|0)!=(e|0)):0){c[d+32>>2]=f;c[b>>2]=e;e=d+40|0;c[e>>2]=(c[e>>2]|0)+1;if((c[d+36>>2]|0)==1?(c[d+24>>2]|0)==2:0){a[d+54|0]=1}c[d+44>>2]=4;i=g;return}if((f|0)!=1){i=g;return}c[d+32>>2]=1;i=g;return}function Im(b,d,e,f,g,h){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0;j=i;if((b|0)!=(c[d+8>>2]|0)){k=d+52|0;l=a[k]|0;m=d+53|0;n=a[m]|0;o=c[b+12>>2]|0;p=b+(o<<3)+16|0;a[k]=0;a[m]=0;q=c[b+20>>2]|0;r=q>>8;if((q&1|0)==0){s=r}else{s=c[(c[f>>2]|0)+r>>2]|0}r=c[b+16>>2]|0;uc[c[(c[r>>2]|0)+20>>2]&15](r,d,e,f+s|0,(q&2|0)!=0?g:2,h);a:do{if((o|0)>1){q=d+24|0;s=b+8|0;r=d+54|0;t=b+24|0;do{if((a[r]|0)!=0){break a}if((a[k]|0)==0){if((a[m]|0)!=0?(c[s>>2]&1|0)==0:0){break a}}else{if((c[q>>2]|0)==1){break a}if((c[s>>2]&2|0)==0){break a}}a[k]=0;a[m]=0;u=c[t+4>>2]|0;v=u>>8;if((u&1|0)==0){w=v}else{w=c[(c[f>>2]|0)+v>>2]|0}v=c[t>>2]|0;uc[c[(c[v>>2]|0)+20>>2]&15](v,d,e,f+w|0,(u&2|0)!=0?g:2,h);t=t+8|0}while(t>>>0<p>>>0)}}while(0);a[k]=l;a[m]=n;i=j;return}a[d+53|0]=1;if((c[d+4>>2]|0)!=(f|0)){i=j;return}a[d+52|0]=1;f=d+16|0;n=c[f>>2]|0;if((n|0)==0){c[f>>2]=e;c[d+24>>2]=g;c[d+36>>2]=1;if(!((c[d+48>>2]|0)==1&(g|0)==1)){i=j;return}a[d+54|0]=1;i=j;return}if((n|0)!=(e|0)){e=d+36|0;c[e>>2]=(c[e>>2]|0)+1;a[d+54|0]=1;i=j;return}e=d+24|0;n=c[e>>2]|0;if((n|0)==2){c[e>>2]=g;x=g}else{x=n}if(!((c[d+48>>2]|0)==1&(x|0)==1)){i=j;return}a[d+54|0]=1;i=j;return}function Jm(b,d,e,f,g,h){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0,k=0,l=0;j=i;if((b|0)!=(c[d+8>>2]|0)){k=c[b+8>>2]|0;uc[c[(c[k>>2]|0)+20>>2]&15](k,d,e,f,g,h);i=j;return}a[d+53|0]=1;if((c[d+4>>2]|0)!=(f|0)){i=j;return}a[d+52|0]=1;f=d+16|0;h=c[f>>2]|0;if((h|0)==0){c[f>>2]=e;c[d+24>>2]=g;c[d+36>>2]=1;if(!((c[d+48>>2]|0)==1&(g|0)==1)){i=j;return}a[d+54|0]=1;i=j;return}if((h|0)!=(e|0)){e=d+36|0;c[e>>2]=(c[e>>2]|0)+1;a[d+54|0]=1;i=j;return}e=d+24|0;h=c[e>>2]|0;if((h|0)==2){c[e>>2]=g;l=g}else{l=h}if(!((c[d+48>>2]|0)==1&(l|0)==1)){i=j;return}a[d+54|0]=1;i=j;return}function Km(b,d,e,f,g,h){b=b|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;var j=0;h=i;if((c[d+8>>2]|0)!=(b|0)){i=h;return}a[d+53|0]=1;if((c[d+4>>2]|0)!=(f|0)){i=h;return}a[d+52|0]=1;f=d+16|0;b=c[f>>2]|0;if((b|0)==0){c[f>>2]=e;c[d+24>>2]=g;c[d+36>>2]=1;if(!((c[d+48>>2]|0)==1&(g|0)==1)){i=h;return}a[d+54|0]=1;i=h;return}if((b|0)!=(e|0)){e=d+36|0;c[e>>2]=(c[e>>2]|0)+1;a[d+54|0]=1;i=h;return}e=d+24|0;b=c[e>>2]|0;if((b|0)==2){c[e>>2]=g;j=g}else{j=b}if(!((c[d+48>>2]|0)==1&(j|0)==1)){i=h;return}a[d+54|0]=1;i=h;return}function Lm(a){a=a|0;var b=0,d=0,e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0,T=0,U=0,V=0,W=0,X=0,Y=0,Z=0,_=0,$=0,aa=0,ba=0,ca=0,da=0,ea=0,fa=0,ga=0,ha=0,ia=0,ja=0,ka=0,la=0,ma=0,na=0,oa=0,pa=0,qa=0,ra=0,sa=0,ta=0,ua=0,va=0,wa=0,xa=0,ya=0,za=0,Aa=0,Ba=0,Ca=0,Da=0,Ea=0,Fa=0,Ga=0,Ha=0,Ia=0,Ja=0,Ka=0;b=i;do{if(a>>>0<245){if(a>>>0<11){d=16}else{d=a+11&-8}e=d>>>3;f=c[4144]|0;g=f>>>e;if((g&3|0)!=0){h=(g&1^1)+e|0;j=h<<1;k=16616+(j<<2)|0;l=16616+(j+2<<2)|0;j=c[l>>2]|0;m=j+8|0;n=c[m>>2]|0;do{if((k|0)!=(n|0)){if(n>>>0<(c[16592>>2]|0)>>>0){Ub()}o=n+12|0;if((c[o>>2]|0)==(j|0)){c[o>>2]=k;c[l>>2]=n;break}else{Ub()}}else{c[4144]=f&~(1<<h)}}while(0);n=h<<3;c[j+4>>2]=n|3;l=j+(n|4)|0;c[l>>2]=c[l>>2]|1;p=m;i=b;return p|0}if(d>>>0>(c[16584>>2]|0)>>>0){if((g|0)!=0){l=2<<e;n=g<<e&(l|0-l);l=(n&0-n)+ -1|0;n=l>>>12&16;k=l>>>n;l=k>>>5&8;o=k>>>l;k=o>>>2&4;q=o>>>k;o=q>>>1&2;r=q>>>o;q=r>>>1&1;s=(l|n|k|o|q)+(r>>>q)|0;q=s<<1;r=16616+(q<<2)|0;o=16616+(q+2<<2)|0;q=c[o>>2]|0;k=q+8|0;n=c[k>>2]|0;do{if((r|0)!=(n|0)){if(n>>>0<(c[16592>>2]|0)>>>0){Ub()}l=n+12|0;if((c[l>>2]|0)==(q|0)){c[l>>2]=r;c[o>>2]=n;break}else{Ub()}}else{c[4144]=f&~(1<<s)}}while(0);f=s<<3;n=f-d|0;c[q+4>>2]=d|3;o=q+d|0;c[q+(d|4)>>2]=n|1;c[q+f>>2]=n;f=c[16584>>2]|0;if((f|0)!=0){r=c[16596>>2]|0;e=f>>>3;f=e<<1;g=16616+(f<<2)|0;m=c[4144]|0;j=1<<e;if((m&j|0)!=0){e=16616+(f+2<<2)|0;h=c[e>>2]|0;if(h>>>0<(c[16592>>2]|0)>>>0){Ub()}else{t=e;u=h}}else{c[4144]=m|j;t=16616+(f+2<<2)|0;u=g}c[t>>2]=r;c[u+12>>2]=r;c[r+8>>2]=u;c[r+12>>2]=g}c[16584>>2]=n;c[16596>>2]=o;p=k;i=b;return p|0}o=c[16580>>2]|0;if((o|0)!=0){n=(o&0-o)+ -1|0;o=n>>>12&16;g=n>>>o;n=g>>>5&8;r=g>>>n;g=r>>>2&4;f=r>>>g;r=f>>>1&2;j=f>>>r;f=j>>>1&1;m=c[16880+((n|o|g|r|f)+(j>>>f)<<2)>>2]|0;f=(c[m+4>>2]&-8)-d|0;j=m;r=m;while(1){m=c[j+16>>2]|0;if((m|0)==0){g=c[j+20>>2]|0;if((g|0)==0){break}else{v=g}}else{v=m}m=(c[v+4>>2]&-8)-d|0;g=m>>>0<f>>>0;f=g?m:f;j=v;r=g?v:r}j=c[16592>>2]|0;if(r>>>0<j>>>0){Ub()}k=r+d|0;if(!(r>>>0<k>>>0)){Ub()}q=c[r+24>>2]|0;s=c[r+12>>2]|0;do{if((s|0)==(r|0)){g=r+20|0;m=c[g>>2]|0;if((m|0)==0){o=r+16|0;n=c[o>>2]|0;if((n|0)==0){w=0;break}else{x=n;y=o}}else{x=m;y=g}while(1){g=x+20|0;m=c[g>>2]|0;if((m|0)!=0){x=m;y=g;continue}g=x+16|0;m=c[g>>2]|0;if((m|0)==0){break}else{x=m;y=g}}if(y>>>0<j>>>0){Ub()}else{c[y>>2]=0;w=x;break}}else{g=c[r+8>>2]|0;if(g>>>0<j>>>0){Ub()}m=g+12|0;if((c[m>>2]|0)!=(r|0)){Ub()}o=s+8|0;if((c[o>>2]|0)==(r|0)){c[m>>2]=s;c[o>>2]=g;w=s;break}else{Ub()}}}while(0);do{if((q|0)!=0){s=c[r+28>>2]|0;j=16880+(s<<2)|0;if((r|0)==(c[j>>2]|0)){c[j>>2]=w;if((w|0)==0){c[16580>>2]=c[16580>>2]&~(1<<s);break}}else{if(q>>>0<(c[16592>>2]|0)>>>0){Ub()}s=q+16|0;if((c[s>>2]|0)==(r|0)){c[s>>2]=w}else{c[q+20>>2]=w}if((w|0)==0){break}}if(w>>>0<(c[16592>>2]|0)>>>0){Ub()}c[w+24>>2]=q;s=c[r+16>>2]|0;do{if((s|0)!=0){if(s>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[w+16>>2]=s;c[s+24>>2]=w;break}}}while(0);s=c[r+20>>2]|0;if((s|0)!=0){if(s>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[w+20>>2]=s;c[s+24>>2]=w;break}}}}while(0);if(f>>>0<16){q=f+d|0;c[r+4>>2]=q|3;s=r+(q+4)|0;c[s>>2]=c[s>>2]|1}else{c[r+4>>2]=d|3;c[r+(d|4)>>2]=f|1;c[r+(f+d)>>2]=f;s=c[16584>>2]|0;if((s|0)!=0){q=c[16596>>2]|0;j=s>>>3;s=j<<1;g=16616+(s<<2)|0;o=c[4144]|0;m=1<<j;if((o&m|0)!=0){j=16616+(s+2<<2)|0;n=c[j>>2]|0;if(n>>>0<(c[16592>>2]|0)>>>0){Ub()}else{z=j;A=n}}else{c[4144]=o|m;z=16616+(s+2<<2)|0;A=g}c[z>>2]=q;c[A+12>>2]=q;c[q+8>>2]=A;c[q+12>>2]=g}c[16584>>2]=f;c[16596>>2]=k}p=r+8|0;i=b;return p|0}else{B=d}}else{B=d}}else{if(!(a>>>0>4294967231)){g=a+11|0;q=g&-8;s=c[16580>>2]|0;if((s|0)!=0){m=0-q|0;o=g>>>8;if((o|0)!=0){if(q>>>0>16777215){C=31}else{g=(o+1048320|0)>>>16&8;n=o<<g;o=(n+520192|0)>>>16&4;j=n<<o;n=(j+245760|0)>>>16&2;h=14-(o|g|n)+(j<<n>>>15)|0;C=q>>>(h+7|0)&1|h<<1}}else{C=0}h=c[16880+(C<<2)>>2]|0;a:do{if((h|0)==0){D=m;E=0;F=0}else{if((C|0)==31){G=0}else{G=25-(C>>>1)|0}n=m;j=0;g=q<<G;o=h;e=0;while(1){l=c[o+4>>2]&-8;H=l-q|0;if(H>>>0<n>>>0){if((l|0)==(q|0)){D=H;E=o;F=o;break a}else{I=H;J=o}}else{I=n;J=e}H=c[o+20>>2]|0;l=c[o+(g>>>31<<2)+16>>2]|0;K=(H|0)==0|(H|0)==(l|0)?j:H;if((l|0)==0){D=I;E=K;F=J;break}else{n=I;j=K;g=g<<1;o=l;e=J}}}}while(0);if((E|0)==0&(F|0)==0){h=2<<C;m=s&(h|0-h);if((m|0)==0){B=q;break}h=(m&0-m)+ -1|0;m=h>>>12&16;r=h>>>m;h=r>>>5&8;k=r>>>h;r=k>>>2&4;f=k>>>r;k=f>>>1&2;e=f>>>k;f=e>>>1&1;L=c[16880+((h|m|r|k|f)+(e>>>f)<<2)>>2]|0}else{L=E}if((L|0)==0){M=D;N=F}else{f=D;e=L;k=F;while(1){r=(c[e+4>>2]&-8)-q|0;m=r>>>0<f>>>0;h=m?r:f;r=m?e:k;m=c[e+16>>2]|0;if((m|0)!=0){f=h;e=m;k=r;continue}m=c[e+20>>2]|0;if((m|0)==0){M=h;N=r;break}else{f=h;e=m;k=r}}}if((N|0)!=0?M>>>0<((c[16584>>2]|0)-q|0)>>>0:0){k=c[16592>>2]|0;if(N>>>0<k>>>0){Ub()}e=N+q|0;if(!(N>>>0<e>>>0)){Ub()}f=c[N+24>>2]|0;s=c[N+12>>2]|0;do{if((s|0)==(N|0)){r=N+20|0;m=c[r>>2]|0;if((m|0)==0){h=N+16|0;o=c[h>>2]|0;if((o|0)==0){O=0;break}else{P=o;Q=h}}else{P=m;Q=r}while(1){r=P+20|0;m=c[r>>2]|0;if((m|0)!=0){P=m;Q=r;continue}r=P+16|0;m=c[r>>2]|0;if((m|0)==0){break}else{P=m;Q=r}}if(Q>>>0<k>>>0){Ub()}else{c[Q>>2]=0;O=P;break}}else{r=c[N+8>>2]|0;if(r>>>0<k>>>0){Ub()}m=r+12|0;if((c[m>>2]|0)!=(N|0)){Ub()}h=s+8|0;if((c[h>>2]|0)==(N|0)){c[m>>2]=s;c[h>>2]=r;O=s;break}else{Ub()}}}while(0);do{if((f|0)!=0){s=c[N+28>>2]|0;k=16880+(s<<2)|0;if((N|0)==(c[k>>2]|0)){c[k>>2]=O;if((O|0)==0){c[16580>>2]=c[16580>>2]&~(1<<s);break}}else{if(f>>>0<(c[16592>>2]|0)>>>0){Ub()}s=f+16|0;if((c[s>>2]|0)==(N|0)){c[s>>2]=O}else{c[f+20>>2]=O}if((O|0)==0){break}}if(O>>>0<(c[16592>>2]|0)>>>0){Ub()}c[O+24>>2]=f;s=c[N+16>>2]|0;do{if((s|0)!=0){if(s>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[O+16>>2]=s;c[s+24>>2]=O;break}}}while(0);s=c[N+20>>2]|0;if((s|0)!=0){if(s>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[O+20>>2]=s;c[s+24>>2]=O;break}}}}while(0);b:do{if(!(M>>>0<16)){c[N+4>>2]=q|3;c[N+(q|4)>>2]=M|1;c[N+(M+q)>>2]=M;f=M>>>3;if(M>>>0<256){s=f<<1;k=16616+(s<<2)|0;r=c[4144]|0;h=1<<f;if((r&h|0)!=0){f=16616+(s+2<<2)|0;m=c[f>>2]|0;if(m>>>0<(c[16592>>2]|0)>>>0){Ub()}else{R=f;S=m}}else{c[4144]=r|h;R=16616+(s+2<<2)|0;S=k}c[R>>2]=e;c[S+12>>2]=e;c[N+(q+8)>>2]=S;c[N+(q+12)>>2]=k;break}k=M>>>8;if((k|0)!=0){if(M>>>0>16777215){T=31}else{s=(k+1048320|0)>>>16&8;h=k<<s;k=(h+520192|0)>>>16&4;r=h<<k;h=(r+245760|0)>>>16&2;m=14-(k|s|h)+(r<<h>>>15)|0;T=M>>>(m+7|0)&1|m<<1}}else{T=0}m=16880+(T<<2)|0;c[N+(q+28)>>2]=T;c[N+(q+20)>>2]=0;c[N+(q+16)>>2]=0;h=c[16580>>2]|0;r=1<<T;if((h&r|0)==0){c[16580>>2]=h|r;c[m>>2]=e;c[N+(q+24)>>2]=m;c[N+(q+12)>>2]=e;c[N+(q+8)>>2]=e;break}r=c[m>>2]|0;if((T|0)==31){U=0}else{U=25-(T>>>1)|0}c:do{if((c[r+4>>2]&-8|0)!=(M|0)){m=M<<U;h=r;while(1){V=h+(m>>>31<<2)+16|0;s=c[V>>2]|0;if((s|0)==0){break}if((c[s+4>>2]&-8|0)==(M|0)){W=s;break c}else{m=m<<1;h=s}}if(V>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[V>>2]=e;c[N+(q+24)>>2]=h;c[N+(q+12)>>2]=e;c[N+(q+8)>>2]=e;break b}}else{W=r}}while(0);r=W+8|0;m=c[r>>2]|0;s=c[16592>>2]|0;if(W>>>0<s>>>0){Ub()}if(m>>>0<s>>>0){Ub()}else{c[m+12>>2]=e;c[r>>2]=e;c[N+(q+8)>>2]=m;c[N+(q+12)>>2]=W;c[N+(q+24)>>2]=0;break}}else{m=M+q|0;c[N+4>>2]=m|3;r=N+(m+4)|0;c[r>>2]=c[r>>2]|1}}while(0);p=N+8|0;i=b;return p|0}else{B=q}}else{B=q}}else{B=-1}}}while(0);N=c[16584>>2]|0;if(!(B>>>0>N>>>0)){M=N-B|0;W=c[16596>>2]|0;if(M>>>0>15){c[16596>>2]=W+B;c[16584>>2]=M;c[W+(B+4)>>2]=M|1;c[W+N>>2]=M;c[W+4>>2]=B|3}else{c[16584>>2]=0;c[16596>>2]=0;c[W+4>>2]=N|3;M=W+(N+4)|0;c[M>>2]=c[M>>2]|1}p=W+8|0;i=b;return p|0}W=c[16588>>2]|0;if(B>>>0<W>>>0){M=W-B|0;c[16588>>2]=M;W=c[16600>>2]|0;c[16600>>2]=W+B;c[W+(B+4)>>2]=M|1;c[W+4>>2]=B|3;p=W+8|0;i=b;return p|0}do{if((c[4262]|0)==0){W=sb(30)|0;if((W+ -1&W|0)==0){c[17056>>2]=W;c[17052>>2]=W;c[17060>>2]=-1;c[17064>>2]=-1;c[17068>>2]=0;c[17020>>2]=0;c[4262]=(Xb(0)|0)&-16^1431655768;break}else{Ub()}}}while(0);W=B+48|0;M=c[17056>>2]|0;N=B+47|0;V=M+N|0;U=0-M|0;M=V&U;if(!(M>>>0>B>>>0)){p=0;i=b;return p|0}T=c[17016>>2]|0;if((T|0)!=0?(S=c[17008>>2]|0,R=S+M|0,R>>>0<=S>>>0|R>>>0>T>>>0):0){p=0;i=b;return p|0}d:do{if((c[17020>>2]&4|0)==0){T=c[16600>>2]|0;e:do{if((T|0)!=0){R=17024|0;while(1){S=c[R>>2]|0;if(!(S>>>0>T>>>0)?(X=R+4|0,(S+(c[X>>2]|0)|0)>>>0>T>>>0):0){break}S=c[R+8>>2]|0;if((S|0)==0){Y=182;break e}else{R=S}}if((R|0)!=0){S=V-(c[16588>>2]|0)&U;if(S>>>0<2147483647){O=ib(S|0)|0;P=(O|0)==((c[R>>2]|0)+(c[X>>2]|0)|0);Z=O;_=S;$=P?O:-1;aa=P?S:0;Y=191}else{ba=0}}else{Y=182}}else{Y=182}}while(0);do{if((Y|0)==182){T=ib(0)|0;if((T|0)!=(-1|0)){q=T;S=c[17052>>2]|0;P=S+ -1|0;if((P&q|0)==0){ca=M}else{ca=M-q+(P+q&0-S)|0}S=c[17008>>2]|0;q=S+ca|0;if(ca>>>0>B>>>0&ca>>>0<2147483647){P=c[17016>>2]|0;if((P|0)!=0?q>>>0<=S>>>0|q>>>0>P>>>0:0){ba=0;break}P=ib(ca|0)|0;q=(P|0)==(T|0);Z=P;_=ca;$=q?T:-1;aa=q?ca:0;Y=191}else{ba=0}}else{ba=0}}}while(0);f:do{if((Y|0)==191){q=0-_|0;if(($|0)!=(-1|0)){da=$;ea=aa;Y=202;break d}do{if((Z|0)!=(-1|0)&_>>>0<2147483647&_>>>0<W>>>0?(T=c[17056>>2]|0,P=N-_+T&0-T,P>>>0<2147483647):0){if((ib(P|0)|0)==(-1|0)){ib(q|0)|0;ba=aa;break f}else{fa=P+_|0;break}}else{fa=_}}while(0);if((Z|0)==(-1|0)){ba=aa}else{da=Z;ea=fa;Y=202;break d}}}while(0);c[17020>>2]=c[17020>>2]|4;ga=ba;Y=199}else{ga=0;Y=199}}while(0);if((((Y|0)==199?M>>>0<2147483647:0)?(ba=ib(M|0)|0,M=ib(0)|0,(M|0)!=(-1|0)&(ba|0)!=(-1|0)&ba>>>0<M>>>0):0)?(fa=M-ba|0,M=fa>>>0>(B+40|0)>>>0,M):0){da=ba;ea=M?fa:ga;Y=202}if((Y|0)==202){ga=(c[17008>>2]|0)+ea|0;c[17008>>2]=ga;if(ga>>>0>(c[17012>>2]|0)>>>0){c[17012>>2]=ga}ga=c[16600>>2]|0;g:do{if((ga|0)!=0){fa=17024|0;while(1){ha=c[fa>>2]|0;ia=fa+4|0;ja=c[ia>>2]|0;if((da|0)==(ha+ja|0)){Y=214;break}M=c[fa+8>>2]|0;if((M|0)==0){break}else{fa=M}}if(((Y|0)==214?(c[fa+12>>2]&8|0)==0:0)?ga>>>0>=ha>>>0&ga>>>0<da>>>0:0){c[ia>>2]=ja+ea;M=(c[16588>>2]|0)+ea|0;ba=ga+8|0;if((ba&7|0)==0){ka=0}else{ka=0-ba&7}ba=M-ka|0;c[16600>>2]=ga+ka;c[16588>>2]=ba;c[ga+(ka+4)>>2]=ba|1;c[ga+(M+4)>>2]=40;c[16604>>2]=c[17064>>2];break}if(da>>>0<(c[16592>>2]|0)>>>0){c[16592>>2]=da}M=da+ea|0;ba=17024|0;while(1){if((c[ba>>2]|0)==(M|0)){Y=224;break}Z=c[ba+8>>2]|0;if((Z|0)==0){break}else{ba=Z}}if((Y|0)==224?(c[ba+12>>2]&8|0)==0:0){c[ba>>2]=da;M=ba+4|0;c[M>>2]=(c[M>>2]|0)+ea;M=da+8|0;if((M&7|0)==0){la=0}else{la=0-M&7}M=da+(ea+8)|0;if((M&7|0)==0){ma=0}else{ma=0-M&7}M=da+(ma+ea)|0;fa=la+B|0;Z=da+fa|0;aa=M-(da+la)-B|0;c[da+(la+4)>>2]=B|3;h:do{if((M|0)!=(c[16600>>2]|0)){if((M|0)==(c[16596>>2]|0)){_=(c[16584>>2]|0)+aa|0;c[16584>>2]=_;c[16596>>2]=Z;c[da+(fa+4)>>2]=_|1;c[da+(_+fa)>>2]=_;break}_=ea+4|0;N=c[da+(_+ma)>>2]|0;if((N&3|0)==1){W=N&-8;$=N>>>3;do{if(!(N>>>0<256)){ca=c[da+((ma|24)+ea)>>2]|0;X=c[da+(ea+12+ma)>>2]|0;do{if((X|0)==(M|0)){U=ma|16;V=da+(_+U)|0;q=c[V>>2]|0;if((q|0)==0){R=da+(U+ea)|0;U=c[R>>2]|0;if((U|0)==0){na=0;break}else{oa=U;pa=R}}else{oa=q;pa=V}while(1){V=oa+20|0;q=c[V>>2]|0;if((q|0)!=0){oa=q;pa=V;continue}V=oa+16|0;q=c[V>>2]|0;if((q|0)==0){break}else{oa=q;pa=V}}if(pa>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[pa>>2]=0;na=oa;break}}else{V=c[da+((ma|8)+ea)>>2]|0;if(V>>>0<(c[16592>>2]|0)>>>0){Ub()}q=V+12|0;if((c[q>>2]|0)!=(M|0)){Ub()}R=X+8|0;if((c[R>>2]|0)==(M|0)){c[q>>2]=X;c[R>>2]=V;na=X;break}else{Ub()}}}while(0);if((ca|0)!=0){X=c[da+(ea+28+ma)>>2]|0;h=16880+(X<<2)|0;if((M|0)==(c[h>>2]|0)){c[h>>2]=na;if((na|0)==0){c[16580>>2]=c[16580>>2]&~(1<<X);break}}else{if(ca>>>0<(c[16592>>2]|0)>>>0){Ub()}X=ca+16|0;if((c[X>>2]|0)==(M|0)){c[X>>2]=na}else{c[ca+20>>2]=na}if((na|0)==0){break}}if(na>>>0<(c[16592>>2]|0)>>>0){Ub()}c[na+24>>2]=ca;X=ma|16;h=c[da+(X+ea)>>2]|0;do{if((h|0)!=0){if(h>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[na+16>>2]=h;c[h+24>>2]=na;break}}}while(0);h=c[da+(_+X)>>2]|0;if((h|0)!=0){if(h>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[na+20>>2]=h;c[h+24>>2]=na;break}}}}else{h=c[da+((ma|8)+ea)>>2]|0;ca=c[da+(ea+12+ma)>>2]|0;V=16616+($<<1<<2)|0;if((h|0)!=(V|0)){if(h>>>0<(c[16592>>2]|0)>>>0){Ub()}if((c[h+12>>2]|0)!=(M|0)){Ub()}}if((ca|0)==(h|0)){c[4144]=c[4144]&~(1<<$);break}if((ca|0)!=(V|0)){if(ca>>>0<(c[16592>>2]|0)>>>0){Ub()}V=ca+8|0;if((c[V>>2]|0)==(M|0)){qa=V}else{Ub()}}else{qa=ca+8|0}c[h+12>>2]=ca;c[qa>>2]=h}}while(0);ra=da+((W|ma)+ea)|0;sa=W+aa|0}else{ra=M;sa=aa}$=ra+4|0;c[$>>2]=c[$>>2]&-2;c[da+(fa+4)>>2]=sa|1;c[da+(sa+fa)>>2]=sa;$=sa>>>3;if(sa>>>0<256){_=$<<1;N=16616+(_<<2)|0;h=c[4144]|0;ca=1<<$;if((h&ca|0)!=0){$=16616+(_+2<<2)|0;V=c[$>>2]|0;if(V>>>0<(c[16592>>2]|0)>>>0){Ub()}else{ta=$;ua=V}}else{c[4144]=h|ca;ta=16616+(_+2<<2)|0;ua=N}c[ta>>2]=Z;c[ua+12>>2]=Z;c[da+(fa+8)>>2]=ua;c[da+(fa+12)>>2]=N;break}N=sa>>>8;if((N|0)!=0){if(sa>>>0>16777215){va=31}else{_=(N+1048320|0)>>>16&8;ca=N<<_;N=(ca+520192|0)>>>16&4;h=ca<<N;ca=(h+245760|0)>>>16&2;V=14-(N|_|ca)+(h<<ca>>>15)|0;va=sa>>>(V+7|0)&1|V<<1}}else{va=0}V=16880+(va<<2)|0;c[da+(fa+28)>>2]=va;c[da+(fa+20)>>2]=0;c[da+(fa+16)>>2]=0;ca=c[16580>>2]|0;h=1<<va;if((ca&h|0)==0){c[16580>>2]=ca|h;c[V>>2]=Z;c[da+(fa+24)>>2]=V;c[da+(fa+12)>>2]=Z;c[da+(fa+8)>>2]=Z;break}h=c[V>>2]|0;if((va|0)==31){wa=0}else{wa=25-(va>>>1)|0}i:do{if((c[h+4>>2]&-8|0)!=(sa|0)){V=sa<<wa;ca=h;while(1){xa=ca+(V>>>31<<2)+16|0;_=c[xa>>2]|0;if((_|0)==0){break}if((c[_+4>>2]&-8|0)==(sa|0)){ya=_;break i}else{V=V<<1;ca=_}}if(xa>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[xa>>2]=Z;c[da+(fa+24)>>2]=ca;c[da+(fa+12)>>2]=Z;c[da+(fa+8)>>2]=Z;break h}}else{ya=h}}while(0);h=ya+8|0;W=c[h>>2]|0;V=c[16592>>2]|0;if(ya>>>0<V>>>0){Ub()}if(W>>>0<V>>>0){Ub()}else{c[W+12>>2]=Z;c[h>>2]=Z;c[da+(fa+8)>>2]=W;c[da+(fa+12)>>2]=ya;c[da+(fa+24)>>2]=0;break}}else{W=(c[16588>>2]|0)+aa|0;c[16588>>2]=W;c[16600>>2]=Z;c[da+(fa+4)>>2]=W|1}}while(0);p=da+(la|8)|0;i=b;return p|0}fa=17024|0;while(1){za=c[fa>>2]|0;if(!(za>>>0>ga>>>0)?(Aa=c[fa+4>>2]|0,Ba=za+Aa|0,Ba>>>0>ga>>>0):0){break}fa=c[fa+8>>2]|0}fa=za+(Aa+ -39)|0;if((fa&7|0)==0){Ca=0}else{Ca=0-fa&7}fa=za+(Aa+ -47+Ca)|0;Z=fa>>>0<(ga+16|0)>>>0?ga:fa;fa=Z+8|0;aa=da+8|0;if((aa&7|0)==0){Da=0}else{Da=0-aa&7}aa=ea+ -40-Da|0;c[16600>>2]=da+Da;c[16588>>2]=aa;c[da+(Da+4)>>2]=aa|1;c[da+(ea+ -36)>>2]=40;c[16604>>2]=c[17064>>2];c[Z+4>>2]=27;c[fa+0>>2]=c[17024>>2];c[fa+4>>2]=c[17028>>2];c[fa+8>>2]=c[17032>>2];c[fa+12>>2]=c[17036>>2];c[17024>>2]=da;c[17028>>2]=ea;c[17036>>2]=0;c[17032>>2]=fa;fa=Z+28|0;c[fa>>2]=7;if((Z+32|0)>>>0<Ba>>>0){aa=fa;while(1){fa=aa+4|0;c[fa>>2]=7;if((aa+8|0)>>>0<Ba>>>0){aa=fa}else{break}}}if((Z|0)!=(ga|0)){aa=Z-ga|0;fa=ga+(aa+4)|0;c[fa>>2]=c[fa>>2]&-2;c[ga+4>>2]=aa|1;c[ga+aa>>2]=aa;fa=aa>>>3;if(aa>>>0<256){M=fa<<1;ba=16616+(M<<2)|0;W=c[4144]|0;h=1<<fa;if((W&h|0)!=0){fa=16616+(M+2<<2)|0;V=c[fa>>2]|0;if(V>>>0<(c[16592>>2]|0)>>>0){Ub()}else{Ea=fa;Fa=V}}else{c[4144]=W|h;Ea=16616+(M+2<<2)|0;Fa=ba}c[Ea>>2]=ga;c[Fa+12>>2]=ga;c[ga+8>>2]=Fa;c[ga+12>>2]=ba;break}ba=aa>>>8;if((ba|0)!=0){if(aa>>>0>16777215){Ga=31}else{M=(ba+1048320|0)>>>16&8;h=ba<<M;ba=(h+520192|0)>>>16&4;W=h<<ba;h=(W+245760|0)>>>16&2;V=14-(ba|M|h)+(W<<h>>>15)|0;Ga=aa>>>(V+7|0)&1|V<<1}}else{Ga=0}V=16880+(Ga<<2)|0;c[ga+28>>2]=Ga;c[ga+20>>2]=0;c[ga+16>>2]=0;h=c[16580>>2]|0;W=1<<Ga;if((h&W|0)==0){c[16580>>2]=h|W;c[V>>2]=ga;c[ga+24>>2]=V;c[ga+12>>2]=ga;c[ga+8>>2]=ga;break}W=c[V>>2]|0;if((Ga|0)==31){Ha=0}else{Ha=25-(Ga>>>1)|0}j:do{if((c[W+4>>2]&-8|0)!=(aa|0)){V=aa<<Ha;h=W;while(1){Ia=h+(V>>>31<<2)+16|0;M=c[Ia>>2]|0;if((M|0)==0){break}if((c[M+4>>2]&-8|0)==(aa|0)){Ja=M;break j}else{V=V<<1;h=M}}if(Ia>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[Ia>>2]=ga;c[ga+24>>2]=h;c[ga+12>>2]=ga;c[ga+8>>2]=ga;break g}}else{Ja=W}}while(0);W=Ja+8|0;aa=c[W>>2]|0;Z=c[16592>>2]|0;if(Ja>>>0<Z>>>0){Ub()}if(aa>>>0<Z>>>0){Ub()}else{c[aa+12>>2]=ga;c[W>>2]=ga;c[ga+8>>2]=aa;c[ga+12>>2]=Ja;c[ga+24>>2]=0;break}}}else{aa=c[16592>>2]|0;if((aa|0)==0|da>>>0<aa>>>0){c[16592>>2]=da}c[17024>>2]=da;c[17028>>2]=ea;c[17036>>2]=0;c[16612>>2]=c[4262];c[16608>>2]=-1;aa=0;do{W=aa<<1;Z=16616+(W<<2)|0;c[16616+(W+3<<2)>>2]=Z;c[16616+(W+2<<2)>>2]=Z;aa=aa+1|0}while((aa|0)!=32);aa=da+8|0;if((aa&7|0)==0){Ka=0}else{Ka=0-aa&7}aa=ea+ -40-Ka|0;c[16600>>2]=da+Ka;c[16588>>2]=aa;c[da+(Ka+4)>>2]=aa|1;c[da+(ea+ -36)>>2]=40;c[16604>>2]=c[17064>>2]}}while(0);ea=c[16588>>2]|0;if(ea>>>0>B>>>0){da=ea-B|0;c[16588>>2]=da;ea=c[16600>>2]|0;c[16600>>2]=ea+B;c[ea+(B+4)>>2]=da|1;c[ea+4>>2]=B|3;p=ea+8|0;i=b;return p|0}}c[(bb()|0)>>2]=12;p=0;i=b;return p|0}



function Mm(a){a=a|0;var b=0,d=0,e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,I=0,J=0,K=0;b=i;if((a|0)==0){i=b;return}d=a+ -8|0;e=c[16592>>2]|0;if(d>>>0<e>>>0){Ub()}f=c[a+ -4>>2]|0;g=f&3;if((g|0)==1){Ub()}h=f&-8;j=a+(h+ -8)|0;do{if((f&1|0)==0){k=c[d>>2]|0;if((g|0)==0){i=b;return}l=-8-k|0;m=a+l|0;n=k+h|0;if(m>>>0<e>>>0){Ub()}if((m|0)==(c[16596>>2]|0)){o=a+(h+ -4)|0;if((c[o>>2]&3|0)!=3){p=m;q=n;break}c[16584>>2]=n;c[o>>2]=c[o>>2]&-2;c[a+(l+4)>>2]=n|1;c[j>>2]=n;i=b;return}o=k>>>3;if(k>>>0<256){k=c[a+(l+8)>>2]|0;r=c[a+(l+12)>>2]|0;s=16616+(o<<1<<2)|0;if((k|0)!=(s|0)){if(k>>>0<e>>>0){Ub()}if((c[k+12>>2]|0)!=(m|0)){Ub()}}if((r|0)==(k|0)){c[4144]=c[4144]&~(1<<o);p=m;q=n;break}if((r|0)!=(s|0)){if(r>>>0<e>>>0){Ub()}s=r+8|0;if((c[s>>2]|0)==(m|0)){t=s}else{Ub()}}else{t=r+8|0}c[k+12>>2]=r;c[t>>2]=k;p=m;q=n;break}k=c[a+(l+24)>>2]|0;r=c[a+(l+12)>>2]|0;do{if((r|0)==(m|0)){s=a+(l+20)|0;o=c[s>>2]|0;if((o|0)==0){u=a+(l+16)|0;v=c[u>>2]|0;if((v|0)==0){w=0;break}else{x=v;y=u}}else{x=o;y=s}while(1){s=x+20|0;o=c[s>>2]|0;if((o|0)!=0){x=o;y=s;continue}s=x+16|0;o=c[s>>2]|0;if((o|0)==0){break}else{x=o;y=s}}if(y>>>0<e>>>0){Ub()}else{c[y>>2]=0;w=x;break}}else{s=c[a+(l+8)>>2]|0;if(s>>>0<e>>>0){Ub()}o=s+12|0;if((c[o>>2]|0)!=(m|0)){Ub()}u=r+8|0;if((c[u>>2]|0)==(m|0)){c[o>>2]=r;c[u>>2]=s;w=r;break}else{Ub()}}}while(0);if((k|0)!=0){r=c[a+(l+28)>>2]|0;s=16880+(r<<2)|0;if((m|0)==(c[s>>2]|0)){c[s>>2]=w;if((w|0)==0){c[16580>>2]=c[16580>>2]&~(1<<r);p=m;q=n;break}}else{if(k>>>0<(c[16592>>2]|0)>>>0){Ub()}r=k+16|0;if((c[r>>2]|0)==(m|0)){c[r>>2]=w}else{c[k+20>>2]=w}if((w|0)==0){p=m;q=n;break}}if(w>>>0<(c[16592>>2]|0)>>>0){Ub()}c[w+24>>2]=k;r=c[a+(l+16)>>2]|0;do{if((r|0)!=0){if(r>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[w+16>>2]=r;c[r+24>>2]=w;break}}}while(0);r=c[a+(l+20)>>2]|0;if((r|0)!=0){if(r>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[w+20>>2]=r;c[r+24>>2]=w;p=m;q=n;break}}else{p=m;q=n}}else{p=m;q=n}}else{p=d;q=h}}while(0);if(!(p>>>0<j>>>0)){Ub()}d=a+(h+ -4)|0;w=c[d>>2]|0;if((w&1|0)==0){Ub()}if((w&2|0)==0){if((j|0)==(c[16600>>2]|0)){e=(c[16588>>2]|0)+q|0;c[16588>>2]=e;c[16600>>2]=p;c[p+4>>2]=e|1;if((p|0)!=(c[16596>>2]|0)){i=b;return}c[16596>>2]=0;c[16584>>2]=0;i=b;return}if((j|0)==(c[16596>>2]|0)){e=(c[16584>>2]|0)+q|0;c[16584>>2]=e;c[16596>>2]=p;c[p+4>>2]=e|1;c[p+e>>2]=e;i=b;return}e=(w&-8)+q|0;x=w>>>3;do{if(!(w>>>0<256)){y=c[a+(h+16)>>2]|0;t=c[a+(h|4)>>2]|0;do{if((t|0)==(j|0)){g=a+(h+12)|0;f=c[g>>2]|0;if((f|0)==0){r=a+(h+8)|0;k=c[r>>2]|0;if((k|0)==0){z=0;break}else{A=k;B=r}}else{A=f;B=g}while(1){g=A+20|0;f=c[g>>2]|0;if((f|0)!=0){A=f;B=g;continue}g=A+16|0;f=c[g>>2]|0;if((f|0)==0){break}else{A=f;B=g}}if(B>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[B>>2]=0;z=A;break}}else{g=c[a+h>>2]|0;if(g>>>0<(c[16592>>2]|0)>>>0){Ub()}f=g+12|0;if((c[f>>2]|0)!=(j|0)){Ub()}r=t+8|0;if((c[r>>2]|0)==(j|0)){c[f>>2]=t;c[r>>2]=g;z=t;break}else{Ub()}}}while(0);if((y|0)!=0){t=c[a+(h+20)>>2]|0;n=16880+(t<<2)|0;if((j|0)==(c[n>>2]|0)){c[n>>2]=z;if((z|0)==0){c[16580>>2]=c[16580>>2]&~(1<<t);break}}else{if(y>>>0<(c[16592>>2]|0)>>>0){Ub()}t=y+16|0;if((c[t>>2]|0)==(j|0)){c[t>>2]=z}else{c[y+20>>2]=z}if((z|0)==0){break}}if(z>>>0<(c[16592>>2]|0)>>>0){Ub()}c[z+24>>2]=y;t=c[a+(h+8)>>2]|0;do{if((t|0)!=0){if(t>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[z+16>>2]=t;c[t+24>>2]=z;break}}}while(0);t=c[a+(h+12)>>2]|0;if((t|0)!=0){if(t>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[z+20>>2]=t;c[t+24>>2]=z;break}}}}else{t=c[a+h>>2]|0;y=c[a+(h|4)>>2]|0;n=16616+(x<<1<<2)|0;if((t|0)!=(n|0)){if(t>>>0<(c[16592>>2]|0)>>>0){Ub()}if((c[t+12>>2]|0)!=(j|0)){Ub()}}if((y|0)==(t|0)){c[4144]=c[4144]&~(1<<x);break}if((y|0)!=(n|0)){if(y>>>0<(c[16592>>2]|0)>>>0){Ub()}n=y+8|0;if((c[n>>2]|0)==(j|0)){C=n}else{Ub()}}else{C=y+8|0}c[t+12>>2]=y;c[C>>2]=t}}while(0);c[p+4>>2]=e|1;c[p+e>>2]=e;if((p|0)==(c[16596>>2]|0)){c[16584>>2]=e;i=b;return}else{D=e}}else{c[d>>2]=w&-2;c[p+4>>2]=q|1;c[p+q>>2]=q;D=q}q=D>>>3;if(D>>>0<256){w=q<<1;d=16616+(w<<2)|0;e=c[4144]|0;C=1<<q;if((e&C|0)!=0){q=16616+(w+2<<2)|0;j=c[q>>2]|0;if(j>>>0<(c[16592>>2]|0)>>>0){Ub()}else{E=q;F=j}}else{c[4144]=e|C;E=16616+(w+2<<2)|0;F=d}c[E>>2]=p;c[F+12>>2]=p;c[p+8>>2]=F;c[p+12>>2]=d;i=b;return}d=D>>>8;if((d|0)!=0){if(D>>>0>16777215){G=31}else{F=(d+1048320|0)>>>16&8;E=d<<F;d=(E+520192|0)>>>16&4;w=E<<d;E=(w+245760|0)>>>16&2;C=14-(d|F|E)+(w<<E>>>15)|0;G=D>>>(C+7|0)&1|C<<1}}else{G=0}C=16880+(G<<2)|0;c[p+28>>2]=G;c[p+20>>2]=0;c[p+16>>2]=0;E=c[16580>>2]|0;w=1<<G;a:do{if((E&w|0)!=0){F=c[C>>2]|0;if((G|0)==31){H=0}else{H=25-(G>>>1)|0}b:do{if((c[F+4>>2]&-8|0)!=(D|0)){d=D<<H;e=F;while(1){I=e+(d>>>31<<2)+16|0;j=c[I>>2]|0;if((j|0)==0){break}if((c[j+4>>2]&-8|0)==(D|0)){J=j;break b}else{d=d<<1;e=j}}if(I>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[I>>2]=p;c[p+24>>2]=e;c[p+12>>2]=p;c[p+8>>2]=p;break a}}else{J=F}}while(0);F=J+8|0;d=c[F>>2]|0;j=c[16592>>2]|0;if(J>>>0<j>>>0){Ub()}if(d>>>0<j>>>0){Ub()}else{c[d+12>>2]=p;c[F>>2]=p;c[p+8>>2]=d;c[p+12>>2]=J;c[p+24>>2]=0;break}}else{c[16580>>2]=E|w;c[C>>2]=p;c[p+24>>2]=C;c[p+12>>2]=p;c[p+8>>2]=p}}while(0);p=(c[16608>>2]|0)+ -1|0;c[16608>>2]=p;if((p|0)==0){K=17032|0}else{i=b;return}while(1){p=c[K>>2]|0;if((p|0)==0){break}else{K=p+8|0}}c[16608>>2]=-1;i=b;return}function Nm(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,g=0,h=0,j=0;d=i;do{if((a|0)!=0){if(b>>>0>4294967231){c[(bb()|0)>>2]=12;e=0;break}if(b>>>0<11){f=16}else{f=b+11&-8}g=Om(a+ -8|0,f)|0;if((g|0)!=0){e=g+8|0;break}g=Lm(b)|0;if((g|0)==0){e=0}else{h=c[a+ -4>>2]|0;j=(h&-8)-((h&3|0)==0?8:4)|0;nn(g|0,a|0,(j>>>0<b>>>0?j:b)|0)|0;Mm(a);e=g}}else{e=Lm(b)|0}}while(0);i=d;return e|0}function Om(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0;d=i;e=a+4|0;f=c[e>>2]|0;g=f&-8;h=a+g|0;j=c[16592>>2]|0;if(a>>>0<j>>>0){Ub()}k=f&3;if(!((k|0)!=1&a>>>0<h>>>0)){Ub()}l=a+(g|4)|0;m=c[l>>2]|0;if((m&1|0)==0){Ub()}if((k|0)==0){if(b>>>0<256){n=0;i=d;return n|0}if(!(g>>>0<(b+4|0)>>>0)?!((g-b|0)>>>0>c[17056>>2]<<1>>>0):0){n=a;i=d;return n|0}n=0;i=d;return n|0}if(!(g>>>0<b>>>0)){k=g-b|0;if(!(k>>>0>15)){n=a;i=d;return n|0}c[e>>2]=f&1|b|2;c[a+(b+4)>>2]=k|3;c[l>>2]=c[l>>2]|1;Rm(a+b|0,k);n=a;i=d;return n|0}if((h|0)==(c[16600>>2]|0)){k=(c[16588>>2]|0)+g|0;if(!(k>>>0>b>>>0)){n=0;i=d;return n|0}l=k-b|0;c[e>>2]=f&1|b|2;c[a+(b+4)>>2]=l|1;c[16600>>2]=a+b;c[16588>>2]=l;n=a;i=d;return n|0}if((h|0)==(c[16596>>2]|0)){l=(c[16584>>2]|0)+g|0;if(l>>>0<b>>>0){n=0;i=d;return n|0}k=l-b|0;if(k>>>0>15){c[e>>2]=f&1|b|2;c[a+(b+4)>>2]=k|1;c[a+l>>2]=k;o=a+(l+4)|0;c[o>>2]=c[o>>2]&-2;p=a+b|0;q=k}else{c[e>>2]=f&1|l|2;f=a+(l+4)|0;c[f>>2]=c[f>>2]|1;p=0;q=0}c[16584>>2]=q;c[16596>>2]=p;n=a;i=d;return n|0}if((m&2|0)!=0){n=0;i=d;return n|0}p=(m&-8)+g|0;if(p>>>0<b>>>0){n=0;i=d;return n|0}q=p-b|0;f=m>>>3;do{if(!(m>>>0<256)){l=c[a+(g+24)>>2]|0;k=c[a+(g+12)>>2]|0;do{if((k|0)==(h|0)){o=a+(g+20)|0;r=c[o>>2]|0;if((r|0)==0){s=a+(g+16)|0;t=c[s>>2]|0;if((t|0)==0){u=0;break}else{v=t;w=s}}else{v=r;w=o}while(1){o=v+20|0;r=c[o>>2]|0;if((r|0)!=0){v=r;w=o;continue}o=v+16|0;r=c[o>>2]|0;if((r|0)==0){break}else{v=r;w=o}}if(w>>>0<j>>>0){Ub()}else{c[w>>2]=0;u=v;break}}else{o=c[a+(g+8)>>2]|0;if(o>>>0<j>>>0){Ub()}r=o+12|0;if((c[r>>2]|0)!=(h|0)){Ub()}s=k+8|0;if((c[s>>2]|0)==(h|0)){c[r>>2]=k;c[s>>2]=o;u=k;break}else{Ub()}}}while(0);if((l|0)!=0){k=c[a+(g+28)>>2]|0;o=16880+(k<<2)|0;if((h|0)==(c[o>>2]|0)){c[o>>2]=u;if((u|0)==0){c[16580>>2]=c[16580>>2]&~(1<<k);break}}else{if(l>>>0<(c[16592>>2]|0)>>>0){Ub()}k=l+16|0;if((c[k>>2]|0)==(h|0)){c[k>>2]=u}else{c[l+20>>2]=u}if((u|0)==0){break}}if(u>>>0<(c[16592>>2]|0)>>>0){Ub()}c[u+24>>2]=l;k=c[a+(g+16)>>2]|0;do{if((k|0)!=0){if(k>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[u+16>>2]=k;c[k+24>>2]=u;break}}}while(0);k=c[a+(g+20)>>2]|0;if((k|0)!=0){if(k>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[u+20>>2]=k;c[k+24>>2]=u;break}}}}else{k=c[a+(g+8)>>2]|0;l=c[a+(g+12)>>2]|0;o=16616+(f<<1<<2)|0;if((k|0)!=(o|0)){if(k>>>0<j>>>0){Ub()}if((c[k+12>>2]|0)!=(h|0)){Ub()}}if((l|0)==(k|0)){c[4144]=c[4144]&~(1<<f);break}if((l|0)!=(o|0)){if(l>>>0<j>>>0){Ub()}o=l+8|0;if((c[o>>2]|0)==(h|0)){x=o}else{Ub()}}else{x=l+8|0}c[k+12>>2]=l;c[x>>2]=k}}while(0);if(q>>>0<16){c[e>>2]=p|c[e>>2]&1|2;x=a+(p|4)|0;c[x>>2]=c[x>>2]|1;n=a;i=d;return n|0}else{c[e>>2]=c[e>>2]&1|b|2;c[a+(b+4)>>2]=q|3;e=a+(p|4)|0;c[e>>2]=c[e>>2]|1;Rm(a+b|0,q);n=a;i=d;return n|0}return 0}function Pm(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0;d=i;e=a>>>0<16?16:a;if((e+ -1&e|0)==0){f=e}else{a=16;while(1){if(a>>>0<e>>>0){a=a<<1}else{f=a;break}}}if(!((-64-f|0)>>>0>b>>>0)){c[(bb()|0)>>2]=12;g=0;i=d;return g|0}if(b>>>0<11){h=16}else{h=b+11&-8}b=Lm(f+12+h|0)|0;if((b|0)==0){g=0;i=d;return g|0}a=b+ -8|0;e=f+ -1|0;do{if((b&e|0)!=0){j=b+e&0-f;k=j+ -8|0;l=a;if((k-l|0)>>>0>15){m=k}else{m=j+(f+ -8)|0}j=m-l|0;l=b+ -4|0;k=c[l>>2]|0;n=(k&-8)-j|0;if((k&3|0)==0){c[m>>2]=(c[a>>2]|0)+j;c[m+4>>2]=n;o=m;break}else{k=m+4|0;c[k>>2]=n|c[k>>2]&1|2;k=m+(n+4)|0;c[k>>2]=c[k>>2]|1;c[l>>2]=j|c[l>>2]&1|2;l=b+(j+ -4)|0;c[l>>2]=c[l>>2]|1;Rm(a,j);o=m;break}}else{o=a}}while(0);a=o+4|0;m=c[a>>2]|0;if((m&3|0)!=0?(b=m&-8,b>>>0>(h+16|0)>>>0):0){f=b-h|0;c[a>>2]=h|m&1|2;c[o+(h|4)>>2]=f|3;m=o+(b|4)|0;c[m>>2]=c[m>>2]|1;Rm(o+h|0,f)}g=o+8|0;i=d;return g|0}function Qm(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,g=0,h=0,j=0;e=i;if((b|0)!=8){f=b>>>2;if(!((b&3|0)!=0|(f|0)==0)?(f+1073741823&f|0)==0:0){if((-64-b|0)>>>0<d>>>0){g=12}else{h=Pm(b>>>0<16?16:b,d)|0;j=7}}else{g=22}}else{h=Lm(d)|0;j=7}if((j|0)==7){if((h|0)==0){g=12}else{c[a>>2]=h;g=0}}i=e;return g|0}function Rm(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0;d=i;e=a+b|0;f=c[a+4>>2]|0;do{if((f&1|0)==0){g=c[a>>2]|0;if((f&3|0)==0){i=d;return}h=a+(0-g)|0;j=g+b|0;k=c[16592>>2]|0;if(h>>>0<k>>>0){Ub()}if((h|0)==(c[16596>>2]|0)){l=a+(b+4)|0;if((c[l>>2]&3|0)!=3){m=h;n=j;break}c[16584>>2]=j;c[l>>2]=c[l>>2]&-2;c[a+(4-g)>>2]=j|1;c[e>>2]=j;i=d;return}l=g>>>3;if(g>>>0<256){o=c[a+(8-g)>>2]|0;p=c[a+(12-g)>>2]|0;q=16616+(l<<1<<2)|0;if((o|0)!=(q|0)){if(o>>>0<k>>>0){Ub()}if((c[o+12>>2]|0)!=(h|0)){Ub()}}if((p|0)==(o|0)){c[4144]=c[4144]&~(1<<l);m=h;n=j;break}if((p|0)!=(q|0)){if(p>>>0<k>>>0){Ub()}q=p+8|0;if((c[q>>2]|0)==(h|0)){r=q}else{Ub()}}else{r=p+8|0}c[o+12>>2]=p;c[r>>2]=o;m=h;n=j;break}o=c[a+(24-g)>>2]|0;p=c[a+(12-g)>>2]|0;do{if((p|0)==(h|0)){q=16-g|0;l=a+(q+4)|0;s=c[l>>2]|0;if((s|0)==0){t=a+q|0;q=c[t>>2]|0;if((q|0)==0){u=0;break}else{v=q;w=t}}else{v=s;w=l}while(1){l=v+20|0;s=c[l>>2]|0;if((s|0)!=0){v=s;w=l;continue}l=v+16|0;s=c[l>>2]|0;if((s|0)==0){break}else{v=s;w=l}}if(w>>>0<k>>>0){Ub()}else{c[w>>2]=0;u=v;break}}else{l=c[a+(8-g)>>2]|0;if(l>>>0<k>>>0){Ub()}s=l+12|0;if((c[s>>2]|0)!=(h|0)){Ub()}t=p+8|0;if((c[t>>2]|0)==(h|0)){c[s>>2]=p;c[t>>2]=l;u=p;break}else{Ub()}}}while(0);if((o|0)!=0){p=c[a+(28-g)>>2]|0;k=16880+(p<<2)|0;if((h|0)==(c[k>>2]|0)){c[k>>2]=u;if((u|0)==0){c[16580>>2]=c[16580>>2]&~(1<<p);m=h;n=j;break}}else{if(o>>>0<(c[16592>>2]|0)>>>0){Ub()}p=o+16|0;if((c[p>>2]|0)==(h|0)){c[p>>2]=u}else{c[o+20>>2]=u}if((u|0)==0){m=h;n=j;break}}if(u>>>0<(c[16592>>2]|0)>>>0){Ub()}c[u+24>>2]=o;p=16-g|0;k=c[a+p>>2]|0;do{if((k|0)!=0){if(k>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[u+16>>2]=k;c[k+24>>2]=u;break}}}while(0);k=c[a+(p+4)>>2]|0;if((k|0)!=0){if(k>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[u+20>>2]=k;c[k+24>>2]=u;m=h;n=j;break}}else{m=h;n=j}}else{m=h;n=j}}else{m=a;n=b}}while(0);u=c[16592>>2]|0;if(e>>>0<u>>>0){Ub()}v=a+(b+4)|0;w=c[v>>2]|0;if((w&2|0)==0){if((e|0)==(c[16600>>2]|0)){r=(c[16588>>2]|0)+n|0;c[16588>>2]=r;c[16600>>2]=m;c[m+4>>2]=r|1;if((m|0)!=(c[16596>>2]|0)){i=d;return}c[16596>>2]=0;c[16584>>2]=0;i=d;return}if((e|0)==(c[16596>>2]|0)){r=(c[16584>>2]|0)+n|0;c[16584>>2]=r;c[16596>>2]=m;c[m+4>>2]=r|1;c[m+r>>2]=r;i=d;return}r=(w&-8)+n|0;f=w>>>3;do{if(!(w>>>0<256)){k=c[a+(b+24)>>2]|0;g=c[a+(b+12)>>2]|0;do{if((g|0)==(e|0)){o=a+(b+20)|0;l=c[o>>2]|0;if((l|0)==0){t=a+(b+16)|0;s=c[t>>2]|0;if((s|0)==0){x=0;break}else{y=s;z=t}}else{y=l;z=o}while(1){o=y+20|0;l=c[o>>2]|0;if((l|0)!=0){y=l;z=o;continue}o=y+16|0;l=c[o>>2]|0;if((l|0)==0){break}else{y=l;z=o}}if(z>>>0<u>>>0){Ub()}else{c[z>>2]=0;x=y;break}}else{o=c[a+(b+8)>>2]|0;if(o>>>0<u>>>0){Ub()}l=o+12|0;if((c[l>>2]|0)!=(e|0)){Ub()}t=g+8|0;if((c[t>>2]|0)==(e|0)){c[l>>2]=g;c[t>>2]=o;x=g;break}else{Ub()}}}while(0);if((k|0)!=0){g=c[a+(b+28)>>2]|0;j=16880+(g<<2)|0;if((e|0)==(c[j>>2]|0)){c[j>>2]=x;if((x|0)==0){c[16580>>2]=c[16580>>2]&~(1<<g);break}}else{if(k>>>0<(c[16592>>2]|0)>>>0){Ub()}g=k+16|0;if((c[g>>2]|0)==(e|0)){c[g>>2]=x}else{c[k+20>>2]=x}if((x|0)==0){break}}if(x>>>0<(c[16592>>2]|0)>>>0){Ub()}c[x+24>>2]=k;g=c[a+(b+16)>>2]|0;do{if((g|0)!=0){if(g>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[x+16>>2]=g;c[g+24>>2]=x;break}}}while(0);g=c[a+(b+20)>>2]|0;if((g|0)!=0){if(g>>>0<(c[16592>>2]|0)>>>0){Ub()}else{c[x+20>>2]=g;c[g+24>>2]=x;break}}}}else{g=c[a+(b+8)>>2]|0;k=c[a+(b+12)>>2]|0;j=16616+(f<<1<<2)|0;if((g|0)!=(j|0)){if(g>>>0<u>>>0){Ub()}if((c[g+12>>2]|0)!=(e|0)){Ub()}}if((k|0)==(g|0)){c[4144]=c[4144]&~(1<<f);break}if((k|0)!=(j|0)){if(k>>>0<u>>>0){Ub()}j=k+8|0;if((c[j>>2]|0)==(e|0)){A=j}else{Ub()}}else{A=k+8|0}c[g+12>>2]=k;c[A>>2]=g}}while(0);c[m+4>>2]=r|1;c[m+r>>2]=r;if((m|0)==(c[16596>>2]|0)){c[16584>>2]=r;i=d;return}else{B=r}}else{c[v>>2]=w&-2;c[m+4>>2]=n|1;c[m+n>>2]=n;B=n}n=B>>>3;if(B>>>0<256){w=n<<1;v=16616+(w<<2)|0;r=c[4144]|0;A=1<<n;if((r&A|0)!=0){n=16616+(w+2<<2)|0;e=c[n>>2]|0;if(e>>>0<(c[16592>>2]|0)>>>0){Ub()}else{C=n;D=e}}else{c[4144]=r|A;C=16616+(w+2<<2)|0;D=v}c[C>>2]=m;c[D+12>>2]=m;c[m+8>>2]=D;c[m+12>>2]=v;i=d;return}v=B>>>8;if((v|0)!=0){if(B>>>0>16777215){E=31}else{D=(v+1048320|0)>>>16&8;C=v<<D;v=(C+520192|0)>>>16&4;w=C<<v;C=(w+245760|0)>>>16&2;A=14-(v|D|C)+(w<<C>>>15)|0;E=B>>>(A+7|0)&1|A<<1}}else{E=0}A=16880+(E<<2)|0;c[m+28>>2]=E;c[m+20>>2]=0;c[m+16>>2]=0;C=c[16580>>2]|0;w=1<<E;if((C&w|0)==0){c[16580>>2]=C|w;c[A>>2]=m;c[m+24>>2]=A;c[m+12>>2]=m;c[m+8>>2]=m;i=d;return}w=c[A>>2]|0;if((E|0)==31){F=0}else{F=25-(E>>>1)|0}a:do{if((c[w+4>>2]&-8|0)==(B|0)){G=w}else{E=B<<F;A=w;while(1){H=A+(E>>>31<<2)+16|0;C=c[H>>2]|0;if((C|0)==0){break}if((c[C+4>>2]&-8|0)==(B|0)){G=C;break a}else{E=E<<1;A=C}}if(H>>>0<(c[16592>>2]|0)>>>0){Ub()}c[H>>2]=m;c[m+24>>2]=A;c[m+12>>2]=m;c[m+8>>2]=m;i=d;return}}while(0);H=G+8|0;B=c[H>>2]|0;w=c[16592>>2]|0;if(G>>>0<w>>>0){Ub()}if(B>>>0<w>>>0){Ub()}c[B+12>>2]=m;c[H>>2]=m;c[m+8>>2]=B;c[m+12>>2]=G;c[m+24>>2]=0;i=d;return}function Sm(a){a=a|0;var b=0,d=0,e=0,f=0;b=i;d=(a|0)==0?1:a;while(1){e=Lm(d)|0;if((e|0)!=0){f=6;break}a=c[4268]|0;c[4268]=a+0;if((a|0)==0){f=5;break}qc[a&0]()}if((f|0)==5){d=Wa(4)|0;c[d>>2]=17088;xb(d|0,17136,109)}else if((f|0)==6){i=b;return e|0}return 0}function Tm(a){a=a|0;var b=0,c=0;b=i;c=Sm(a)|0;i=b;return c|0}function Um(a){a=a|0;var b=0;b=i;if((a|0)!=0){Mm(a)}i=b;return}function Vm(a){a=a|0;var b=0;b=i;Um(a);i=b;return}function Wm(a){a=a|0;c[a>>2]=17088;return}function Xm(a){a=a|0;var b=0;b=i;Mb(a|0);Um(a);i=b;return}function Ym(a){a=a|0;var b=0;b=i;Mb(a|0);i=b;return}function Zm(a){a=a|0;return 17104}function _m(){var a=0;a=Wa(4)|0;c[a>>2]=17088;xb(a|0,17136,109)}function $m(b,e,f){b=b|0;e=e|0;f=f|0;var g=0,h=0,j=0,k=0,l=0.0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,J=0,K=0,L=0,M=0,N=0,O=0,P=0,Q=0,R=0,S=0.0,U=0,V=0.0,W=0,X=0,Y=0,Z=0,_=0,$=0,aa=0,ba=0,ca=0,da=0,fa=0.0,ga=0,ha=0.0,ia=0,ja=0.0,ka=0,la=0.0,ma=0,na=0.0,oa=0,pa=0,qa=0,ra=0,sa=0,ta=0,ua=0,va=0.0,wa=0,xa=0.0,ya=0,za=0,Aa=0,Ba=0,Ca=0.0,Da=0,Ea=0.0,Fa=0.0,Ga=0,Ha=0.0,Ia=0,Ja=0,Ka=0,La=0,Ma=0,Oa=0,Pa=0,Qa=0,Ra=0,Ta=0,Ua=0,Va=0,Wa=0,Xa=0,Ya=0,Za=0,_a=0,$a=0,ab=0,cb=0,db=0,eb=0,fb=0,gb=0,hb=0,ib=0,jb=0,kb=0,lb=0,mb=0,nb=0,ob=0,pb=0,qb=0,rb=0,sb=0,tb=0,ub=0,vb=0,wb=0,xb=0,yb=0,Ab=0,Bb=0,Cb=0,Db=0,Eb=0,Fb=0,Gb=0,Hb=0,Ib=0,Jb=0,Kb=0,Lb=0,Mb=0,Nb=0,Ob=0,Pb=0,Qb=0,Rb=0,Sb=0,Tb=0,Ub=0,Vb=0,Wb=0,Xb=0,Yb=0,Zb=0,_b=0,$b=0,ac=0,bc=0,cc=0,dc=0,ec=0,fc=0,gc=0,hc=0,ic=0,jc=0,kc=0,lc=0,mc=0,nc=0,oc=0,pc=0,qc=0,rc=0,sc=0,tc=0,uc=0,vc=0.0,wc=0,xc=0,yc=0.0,zc=0.0,Ac=0.0,Bc=0.0,Cc=0.0,Dc=0.0,Ec=0,Fc=0,Gc=0.0,Hc=0,Ic=0.0;g=i;i=i+512|0;h=g;if((e|0)==1){j=53;k=-1074}else if((e|0)==2){j=53;k=-1074}else if((e|0)==0){j=24;k=-149}else{l=0.0;i=g;return+l}e=b+4|0;m=b+100|0;do{n=c[e>>2]|0;if(n>>>0<(c[m>>2]|0)>>>0){c[e>>2]=n+1;o=d[n]|0}else{o=cn(b)|0}}while((zb(o|0)|0)!=0);do{if((o|0)==43|(o|0)==45){n=1-(((o|0)==45)<<1)|0;p=c[e>>2]|0;if(p>>>0<(c[m>>2]|0)>>>0){c[e>>2]=p+1;q=d[p]|0;r=n;break}else{q=cn(b)|0;r=n;break}}else{q=o;r=1}}while(0);o=q;q=0;while(1){if((o|32|0)!=(a[17152+q|0]|0)){s=o;t=q;break}do{if(q>>>0<7){n=c[e>>2]|0;if(n>>>0<(c[m>>2]|0)>>>0){c[e>>2]=n+1;u=d[n]|0;break}else{u=cn(b)|0;break}}else{u=o}}while(0);n=q+1|0;if(n>>>0<8){o=u;q=n}else{s=u;t=n;break}}do{if((t|0)==3){v=23}else if((t|0)!=8){u=(f|0)==0;if(!(t>>>0<4|u)){if((t|0)==8){break}else{v=23;break}}a:do{if((t|0)==0){q=s;o=0;while(1){if((q|32|0)!=(a[17168+o|0]|0)){w=q;z=o;break a}do{if(o>>>0<2){n=c[e>>2]|0;if(n>>>0<(c[m>>2]|0)>>>0){c[e>>2]=n+1;A=d[n]|0;break}else{A=cn(b)|0;break}}else{A=q}}while(0);n=o+1|0;if(n>>>0<3){q=A;o=n}else{w=A;z=n;break}}}else{w=s;z=t}}while(0);if((z|0)==3){o=c[e>>2]|0;if(o>>>0<(c[m>>2]|0)>>>0){c[e>>2]=o+1;B=d[o]|0}else{B=cn(b)|0}if((B|0)==40){C=1}else{if((c[m>>2]|0)==0){l=x;i=g;return+l}c[e>>2]=(c[e>>2]|0)+ -1;l=x;i=g;return+l}while(1){o=c[e>>2]|0;if(o>>>0<(c[m>>2]|0)>>>0){c[e>>2]=o+1;D=d[o]|0}else{D=cn(b)|0}if(!((D+ -48|0)>>>0<10|(D+ -65|0)>>>0<26)?!((D+ -97|0)>>>0<26|(D|0)==95):0){break}C=C+1|0}if((D|0)==41){l=x;i=g;return+l}o=(c[m>>2]|0)==0;if(!o){c[e>>2]=(c[e>>2]|0)+ -1}if(u){c[(bb()|0)>>2]=22;bn(b,0);l=0.0;i=g;return+l}if((C|0)==0|o){l=x;i=g;return+l}else{E=C}while(1){o=E+ -1|0;c[e>>2]=(c[e>>2]|0)+ -1;if((o|0)==0){l=x;break}else{E=o}}i=g;return+l}else if((z|0)==0){do{if((w|0)==48){o=c[e>>2]|0;if(o>>>0<(c[m>>2]|0)>>>0){c[e>>2]=o+1;F=d[o]|0}else{F=cn(b)|0}if((F|32|0)!=120){if((c[m>>2]|0)==0){G=48;break}c[e>>2]=(c[e>>2]|0)+ -1;G=48;break}o=c[e>>2]|0;if(o>>>0<(c[m>>2]|0)>>>0){c[e>>2]=o+1;H=d[o]|0;J=0}else{H=cn(b)|0;J=0}while(1){if((H|0)==46){v=70;break}else if((H|0)!=48){K=0;L=0;M=0;N=0;O=H;P=J;Q=0;R=0;S=1.0;U=0;V=0.0;break}o=c[e>>2]|0;if(o>>>0<(c[m>>2]|0)>>>0){c[e>>2]=o+1;H=d[o]|0;J=1;continue}else{H=cn(b)|0;J=1;continue}}b:do{if((v|0)==70){o=c[e>>2]|0;if(o>>>0<(c[m>>2]|0)>>>0){c[e>>2]=o+1;W=d[o]|0}else{W=cn(b)|0}if((W|0)==48){o=-1;q=-1;while(1){n=c[e>>2]|0;if(n>>>0<(c[m>>2]|0)>>>0){c[e>>2]=n+1;X=d[n]|0}else{X=cn(b)|0}if((X|0)!=48){K=0;L=0;M=o;N=q;O=X;P=1;Q=1;R=0;S=1.0;U=0;V=0.0;break b}n=ln(o|0,q|0,-1,-1)|0;o=n;q=I}}else{K=0;L=0;M=0;N=0;O=W;P=J;Q=1;R=0;S=1.0;U=0;V=0.0}}}while(0);c:while(1){q=O+ -48|0;do{if(!(q>>>0<10)){o=O|32;n=(O|0)==46;if(!((o+ -97|0)>>>0<6|n)){Y=O;break c}if(n){if((Q|0)==0){Z=L;_=K;$=L;aa=K;ba=P;ca=1;da=R;fa=S;ga=U;ha=V;break}else{Y=46;break c}}else{ia=(O|0)>57?o+ -87|0:q;v=84;break}}else{ia=q;v=84}}while(0);if((v|0)==84){v=0;do{if(!((K|0)<0|(K|0)==0&L>>>0<8)){if((K|0)<0|(K|0)==0&L>>>0<14){ja=S*.0625;ka=R;la=ja;ma=U;na=V+ja*+(ia|0);break}if((ia|0)!=0&(R|0)==0){ka=1;la=S;ma=U;na=V+S*.5}else{ka=R;la=S;ma=U;na=V}}else{ka=R;la=S;ma=ia+(U<<4)|0;na=V}}while(0);q=ln(L|0,K|0,1,0)|0;Z=M;_=N;$=q;aa=I;ba=1;ca=Q;da=ka;fa=la;ga=ma;ha=na}q=c[e>>2]|0;if(q>>>0<(c[m>>2]|0)>>>0){c[e>>2]=q+1;K=aa;L=$;M=Z;N=_;O=d[q]|0;P=ba;Q=ca;R=da;S=fa;U=ga;V=ha;continue}else{K=aa;L=$;M=Z;N=_;O=cn(b)|0;P=ba;Q=ca;R=da;S=fa;U=ga;V=ha;continue}}if((P|0)==0){q=(c[m>>2]|0)==0;if(!q){c[e>>2]=(c[e>>2]|0)+ -1}if(!u){if(!q?(q=c[e>>2]|0,c[e>>2]=q+ -1,(Q|0)!=0):0){c[e>>2]=q+ -2}}else{bn(b,0)}l=+(r|0)*0.0;i=g;return+l}q=(Q|0)==0;o=q?L:M;n=q?K:N;if((K|0)<0|(K|0)==0&L>>>0<8){q=L;p=K;oa=U;while(1){pa=oa<<4;qa=ln(q|0,p|0,1,0)|0;ra=I;if((ra|0)<0|(ra|0)==0&qa>>>0<8){q=qa;p=ra;oa=pa}else{sa=pa;break}}}else{sa=U}do{if((Y|32|0)==112){oa=an(b,f)|0;p=I;if((oa|0)==0&(p|0)==-2147483648){if(u){bn(b,0);l=0.0;i=g;return+l}else{if((c[m>>2]|0)==0){ta=0;ua=0;break}c[e>>2]=(c[e>>2]|0)+ -1;ta=0;ua=0;break}}else{ta=oa;ua=p}}else{if((c[m>>2]|0)==0){ta=0;ua=0}else{c[e>>2]=(c[e>>2]|0)+ -1;ta=0;ua=0}}}while(0);p=qn(o|0,n|0,2)|0;oa=ln(p|0,I|0,-32,-1)|0;p=ln(oa|0,I|0,ta|0,ua|0)|0;oa=I;if((sa|0)==0){l=+(r|0)*0.0;i=g;return+l}if((oa|0)>0|(oa|0)==0&p>>>0>(0-k|0)>>>0){c[(bb()|0)>>2]=34;l=+(r|0)*1.7976931348623157e+308*1.7976931348623157e+308;i=g;return+l}q=k+ -106|0;pa=((q|0)<0)<<31>>31;if((oa|0)<(pa|0)|(oa|0)==(pa|0)&p>>>0<q>>>0){c[(bb()|0)>>2]=34;l=+(r|0)*2.2250738585072014e-308*2.2250738585072014e-308;i=g;return+l}if((sa|0)>-1){q=p;pa=oa;ra=sa;ja=V;while(1){qa=ra<<1;if(!(ja>=.5)){va=ja;wa=qa}else{va=ja+-1.0;wa=qa|1}xa=ja+va;qa=ln(q|0,pa|0,-1,-1)|0;ya=I;if((wa|0)>-1){q=qa;pa=ya;ra=wa;ja=xa}else{za=qa;Aa=ya;Ba=wa;Ca=xa;break}}}else{za=p;Aa=oa;Ba=sa;Ca=V}ra=kn(32,0,k|0,((k|0)<0)<<31>>31|0)|0;pa=ln(za|0,Aa|0,ra|0,I|0)|0;ra=I;if(0>(ra|0)|0==(ra|0)&j>>>0>pa>>>0){Da=(pa|0)<0?0:pa}else{Da=j}if((Da|0)<53){ja=+(r|0);xa=+Sa(+(+dn(1.0,84-Da|0)),+ja);if((Da|0)<32&Ca!=0.0){pa=Ba&1;Ea=ja;Fa=xa;Ga=(pa^1)+Ba|0;Ha=(pa|0)==0?0.0:Ca}else{Ea=ja;Fa=xa;Ga=Ba;Ha=Ca}}else{Ea=+(r|0);Fa=0.0;Ga=Ba;Ha=Ca}xa=Ea*Ha+(Fa+Ea*+(Ga>>>0))-Fa;if(!(xa!=0.0)){c[(bb()|0)>>2]=34}l=+en(xa,za);i=g;return+l}else{G=w}}while(0);pa=k+j|0;ra=0-pa|0;q=G;n=0;while(1){if((q|0)==46){v=139;break}else if((q|0)!=48){Ia=q;Ja=0;Ka=0;La=n;Ma=0;break}o=c[e>>2]|0;if(o>>>0<(c[m>>2]|0)>>>0){c[e>>2]=o+1;q=d[o]|0;n=1;continue}else{q=cn(b)|0;n=1;continue}}d:do{if((v|0)==139){q=c[e>>2]|0;if(q>>>0<(c[m>>2]|0)>>>0){c[e>>2]=q+1;Oa=d[q]|0}else{Oa=cn(b)|0}if((Oa|0)==48){q=-1;o=-1;while(1){ya=c[e>>2]|0;if(ya>>>0<(c[m>>2]|0)>>>0){c[e>>2]=ya+1;Pa=d[ya]|0}else{Pa=cn(b)|0}if((Pa|0)!=48){Ia=Pa;Ja=q;Ka=o;La=1;Ma=1;break d}ya=ln(q|0,o|0,-1,-1)|0;q=ya;o=I}}else{Ia=Oa;Ja=0;Ka=0;La=n;Ma=1}}}while(0);c[h>>2]=0;n=Ia+ -48|0;o=(Ia|0)==46;e:do{if(n>>>0<10|o){q=h+496|0;oa=Ia;p=0;ya=0;qa=o;Qa=n;Ra=Ja;Ta=Ka;Ua=La;Va=Ma;Wa=0;Xa=0;Ya=0;while(1){do{if(qa){if((Va|0)==0){Za=p;_a=ya;$a=p;ab=ya;cb=Ua;db=1;eb=Wa;fb=Xa;gb=Ya}else{hb=oa;ib=Ra;jb=Ta;kb=p;lb=ya;mb=Ua;nb=Wa;ob=Xa;pb=Ya;break e}}else{qb=ln(p|0,ya|0,1,0)|0;rb=I;sb=(oa|0)!=48;if((Xa|0)>=125){if(!sb){Za=Ra;_a=Ta;$a=qb;ab=rb;cb=Ua;db=Va;eb=Wa;fb=Xa;gb=Ya;break}c[q>>2]=c[q>>2]|1;Za=Ra;_a=Ta;$a=qb;ab=rb;cb=Ua;db=Va;eb=Wa;fb=Xa;gb=Ya;break}tb=h+(Xa<<2)|0;if((Wa|0)==0){ub=Qa}else{ub=oa+ -48+((c[tb>>2]|0)*10|0)|0}c[tb>>2]=ub;tb=Wa+1|0;vb=(tb|0)==9;Za=Ra;_a=Ta;$a=qb;ab=rb;cb=1;db=Va;eb=vb?0:tb;fb=(vb&1)+Xa|0;gb=sb?qb:Ya}}while(0);qb=c[e>>2]|0;if(qb>>>0<(c[m>>2]|0)>>>0){c[e>>2]=qb+1;wb=d[qb]|0}else{wb=cn(b)|0}qb=wb+ -48|0;sb=(wb|0)==46;if(qb>>>0<10|sb){oa=wb;p=$a;ya=ab;qa=sb;Qa=qb;Ra=Za;Ta=_a;Ua=cb;Va=db;Wa=eb;Xa=fb;Ya=gb}else{xb=wb;yb=$a;Ab=Za;Bb=ab;Cb=_a;Db=cb;Eb=db;Fb=eb;Gb=fb;Hb=gb;v=162;break}}}else{xb=Ia;yb=0;Ab=Ja;Bb=0;Cb=Ka;Db=La;Eb=Ma;Fb=0;Gb=0;Hb=0;v=162}}while(0);if((v|0)==162){n=(Eb|0)==0;hb=xb;ib=n?yb:Ab;jb=n?Bb:Cb;kb=yb;lb=Bb;mb=Db;nb=Fb;ob=Gb;pb=Hb}n=(mb|0)!=0;if(n?(hb|32|0)==101:0){o=an(b,f)|0;Ya=I;do{if((o|0)==0&(Ya|0)==-2147483648){if(u){bn(b,0);l=0.0;i=g;return+l}else{if((c[m>>2]|0)==0){Ib=0;Jb=0;break}c[e>>2]=(c[e>>2]|0)+ -1;Ib=0;Jb=0;break}}else{Ib=o;Jb=Ya}}while(0);Ya=ln(Ib|0,Jb|0,ib|0,jb|0)|0;Kb=Ya;Lb=I}else{if((hb|0)>-1?(c[m>>2]|0)!=0:0){c[e>>2]=(c[e>>2]|0)+ -1;Kb=ib;Lb=jb}else{Kb=ib;Lb=jb}}if(!n){c[(bb()|0)>>2]=22;bn(b,0);l=0.0;i=g;return+l}Ya=c[h>>2]|0;if((Ya|0)==0){l=+(r|0)*0.0;i=g;return+l}do{if((Kb|0)==(kb|0)&(Lb|0)==(lb|0)&((lb|0)<0|(lb|0)==0&kb>>>0<10)){if(!(j>>>0>30)?(Ya>>>j|0)!=0:0){break}l=+(r|0)*+(Ya>>>0);i=g;return+l}}while(0);Ya=(k|0)/-2|0;n=((Ya|0)<0)<<31>>31;if((Lb|0)>(n|0)|(Lb|0)==(n|0)&Kb>>>0>Ya>>>0){c[(bb()|0)>>2]=34;l=+(r|0)*1.7976931348623157e+308*1.7976931348623157e+308;i=g;return+l}Ya=k+ -106|0;n=((Ya|0)<0)<<31>>31;if((Lb|0)<(n|0)|(Lb|0)==(n|0)&Kb>>>0<Ya>>>0){c[(bb()|0)>>2]=34;l=+(r|0)*2.2250738585072014e-308*2.2250738585072014e-308;i=g;return+l}if((nb|0)==0){Mb=ob}else{if((nb|0)<9){Ya=h+(ob<<2)|0;n=c[Ya>>2]|0;o=nb;do{n=n*10|0;o=o+1|0}while((o|0)!=9);c[Ya>>2]=n}Mb=ob+1|0}do{if((pb|0)<9?(pb|0)<=(Kb|0)&(Kb|0)<18:0){if((Kb|0)==9){l=+(r|0)*+((c[h>>2]|0)>>>0);i=g;return+l}if((Kb|0)<9){l=+(r|0)*+((c[h>>2]|0)>>>0)/+(c[17184+(8-Kb<<2)>>2]|0);i=g;return+l}o=j+27+(ea(Kb,-3)|0)|0;u=c[h>>2]|0;if((o|0)<=30?(u>>>o|0)!=0:0){break}l=+(r|0)*+(u>>>0)*+(c[17184+(Kb+ -10<<2)>>2]|0);i=g;return+l}}while(0);n=(Kb|0)%9|0;if((n|0)==0){Nb=0;Ob=0;Pb=Kb;Qb=Mb}else{Ya=(Kb|0)>-1?n:n+9|0;n=c[17184+(8-Ya<<2)>>2]|0;if((Mb|0)!=0){u=1e9/(n|0)|0;o=0;Xa=0;Wa=0;Va=Kb;while(1){Ua=h+(Wa<<2)|0;Ta=c[Ua>>2]|0;Ra=((Ta>>>0)/(n>>>0)|0)+Xa|0;c[Ua>>2]=Ra;Rb=ea((Ta>>>0)%(n>>>0)|0,u)|0;Ta=Wa+1|0;if((Wa|0)==(o|0)&(Ra|0)==0){Sb=Ta&127;Tb=Va+ -9|0}else{Sb=o;Tb=Va}if((Ta|0)==(Mb|0)){break}else{o=Sb;Xa=Rb;Wa=Ta;Va=Tb}}if((Rb|0)==0){Ub=Sb;Vb=Tb;Wb=Mb}else{c[h+(Mb<<2)>>2]=Rb;Ub=Sb;Vb=Tb;Wb=Mb+1|0}}else{Ub=0;Vb=Kb;Wb=0}Nb=Ub;Ob=0;Pb=9-Ya+Vb|0;Qb=Wb}f:while(1){Va=h+(Nb<<2)|0;if((Pb|0)<18){Wa=Ob;Xa=Qb;while(1){o=0;u=Xa+127|0;n=Xa;while(1){Ta=u&127;Ra=h+(Ta<<2)|0;Ua=qn(c[Ra>>2]|0,0,29)|0;Qa=ln(Ua|0,I|0,o|0,0)|0;Ua=I;if(Ua>>>0>0|(Ua|0)==0&Qa>>>0>1e9){qa=zn(Qa|0,Ua|0,1e9,0)|0;ya=An(Qa|0,Ua|0,1e9,0)|0;Xb=ya;Yb=qa}else{Xb=Qa;Yb=0}c[Ra>>2]=Xb;Ra=(Ta|0)==(Nb|0);if((Ta|0)!=(n+127&127|0)|Ra){Zb=n}else{Zb=(Xb|0)==0?Ta:n}if(Ra){break}else{o=Yb;u=Ta+ -1|0;n=Zb}}n=Wa+ -29|0;if((Yb|0)==0){Wa=n;Xa=Zb}else{_b=n;$b=Yb;ac=Zb;break}}}else{if((Pb|0)==18){bc=Ob;cc=Qb}else{dc=Nb;ec=Ob;fc=Pb;gc=Qb;break}while(1){if(!((c[Va>>2]|0)>>>0<9007199)){dc=Nb;ec=bc;fc=18;gc=cc;break f}Xa=0;Wa=cc+127|0;n=cc;while(1){u=Wa&127;o=h+(u<<2)|0;Ta=qn(c[o>>2]|0,0,29)|0;Ra=ln(Ta|0,I|0,Xa|0,0)|0;Ta=I;if(Ta>>>0>0|(Ta|0)==0&Ra>>>0>1e9){Qa=zn(Ra|0,Ta|0,1e9,0)|0;qa=An(Ra|0,Ta|0,1e9,0)|0;hc=qa;ic=Qa}else{hc=Ra;ic=0}c[o>>2]=hc;o=(u|0)==(Nb|0);if((u|0)!=(n+127&127|0)|o){jc=n}else{jc=(hc|0)==0?u:n}if(o){break}else{Xa=ic;Wa=u+ -1|0;n=jc}}n=bc+ -29|0;if((ic|0)==0){bc=n;cc=jc}else{_b=n;$b=ic;ac=jc;break}}}Va=Nb+127&127;if((Va|0)==(ac|0)){n=ac+127&127;Wa=h+((ac+126&127)<<2)|0;c[Wa>>2]=c[Wa>>2]|c[h+(n<<2)>>2];kc=n}else{kc=ac}c[h+(Va<<2)>>2]=$b;Nb=Va;Ob=_b;Pb=Pb+9|0;Qb=kc}g:while(1){lc=gc+1&127;Ya=h+((gc+127&127)<<2)|0;Va=dc;n=ec;Wa=fc;while(1){Xa=(Wa|0)==18;u=(Wa|0)>27?9:1;mc=Va;nc=n;while(1){o=0;while(1){Ra=o+mc&127;if((Ra|0)==(gc|0)){oc=2;break}Qa=c[h+(Ra<<2)>>2]|0;Ra=c[17176+(o<<2)>>2]|0;if(Qa>>>0<Ra>>>0){oc=2;break}qa=o+1|0;if(Qa>>>0>Ra>>>0){oc=o;break}if((qa|0)<2){o=qa}else{oc=qa;break}}if((oc|0)==2&Xa){break g}pc=u+nc|0;if((mc|0)==(gc|0)){mc=gc;nc=pc}else{break}}Xa=(1<<u)+ -1|0;o=1e9>>>u;qc=mc;rc=0;qa=mc;sc=Wa;do{Ra=h+(qa<<2)|0;Qa=c[Ra>>2]|0;Ta=(Qa>>>u)+rc|0;c[Ra>>2]=Ta;rc=ea(Qa&Xa,o)|0;Qa=(qa|0)==(qc|0)&(Ta|0)==0;qa=qa+1&127;sc=Qa?sc+ -9|0:sc;qc=Qa?qa:qc}while((qa|0)!=(gc|0));if((rc|0)==0){Va=qc;n=pc;Wa=sc;continue}if((lc|0)!=(qc|0)){break}c[Ya>>2]=c[Ya>>2]|1;Va=qc;n=pc;Wa=sc}c[h+(gc<<2)>>2]=rc;dc=qc;ec=pc;fc=sc;gc=lc}Wa=mc&127;if((Wa|0)==(gc|0)){c[h+(lc+ -1<<2)>>2]=0;tc=lc}else{tc=gc}xa=+((c[h+(Wa<<2)>>2]|0)>>>0);Wa=mc+1&127;if((Wa|0)==(tc|0)){n=tc+1&127;c[h+(n+ -1<<2)>>2]=0;uc=n}else{uc=tc}ja=+(r|0);vc=ja*(xa*1.0e9+ +((c[h+(Wa<<2)>>2]|0)>>>0));Wa=nc+53|0;n=Wa-k|0;if((n|0)<(j|0)){wc=(n|0)<0?0:n;xc=1}else{wc=j;xc=0}if((wc|0)<53){xa=+Sa(+(+dn(1.0,105-wc|0)),+vc);yc=+Na(+vc,+(+dn(1.0,53-wc|0)));zc=xa;Ac=yc;Bc=xa+(vc-yc)}else{zc=0.0;Ac=0.0;Bc=vc}Va=mc+2&127;if((Va|0)!=(uc|0)){Ya=c[h+(Va<<2)>>2]|0;do{if(!(Ya>>>0<5e8)){if(Ya>>>0>5e8){Cc=ja*.75+Ac;break}if((mc+3&127|0)==(uc|0)){Cc=ja*.5+Ac;break}else{Cc=ja*.75+Ac;break}}else{if((Ya|0)==0?(mc+3&127|0)==(uc|0):0){Cc=Ac;break}Cc=ja*.25+Ac}}while(0);if((53-wc|0)>1?!(+Na(+Cc,1.0)!=0.0):0){Dc=Cc+1.0}else{Dc=Cc}}else{Dc=Ac}ja=Bc+Dc-zc;do{if((Wa&2147483647|0)>(-2-pa|0)){if(!(+T(+ja)>=9007199254740992.0)){Ec=xc;Fc=nc;Gc=ja}else{Ec=(xc|0)!=0&(wc|0)==(n|0)?0:xc;Fc=nc+1|0;Gc=ja*.5}if((Fc+50|0)<=(ra|0)?!((Ec|0)!=0&Dc!=0.0):0){Hc=Fc;Ic=Gc;break}c[(bb()|0)>>2]=34;Hc=Fc;Ic=Gc}else{Hc=nc;Ic=ja}}while(0);l=+en(Ic,Hc);i=g;return+l}else{if((c[m>>2]|0)!=0){c[e>>2]=(c[e>>2]|0)+ -1}c[(bb()|0)>>2]=22;bn(b,0);l=0.0;i=g;return+l}}}while(0);if((v|0)==23){v=(c[m>>2]|0)==0;if(!v){c[e>>2]=(c[e>>2]|0)+ -1}if(!(t>>>0<4|(f|0)==0|v)){v=t;do{c[e>>2]=(c[e>>2]|0)+ -1;v=v+ -1|0}while(v>>>0>3)}}l=+(r|0)*y;i=g;return+l}function an(a,b){a=a|0;b=b|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0;e=i;f=a+4|0;g=c[f>>2]|0;h=a+100|0;if(g>>>0<(c[h>>2]|0)>>>0){c[f>>2]=g+1;j=d[g]|0}else{j=cn(a)|0}if((j|0)==43|(j|0)==45){g=(j|0)==45|0;k=c[f>>2]|0;if(k>>>0<(c[h>>2]|0)>>>0){c[f>>2]=k+1;l=d[k]|0}else{l=cn(a)|0}if(!((l+ -48|0)>>>0<10|(b|0)==0)?(c[h>>2]|0)!=0:0){c[f>>2]=(c[f>>2]|0)+ -1;m=l;n=g}else{m=l;n=g}}else{m=j;n=0}if((m+ -48|0)>>>0>9){if((c[h>>2]|0)==0){o=-2147483648;p=0;I=o;i=e;return p|0}c[f>>2]=(c[f>>2]|0)+ -1;o=-2147483648;p=0;I=o;i=e;return p|0}else{q=m;r=0}while(1){s=q+ -48+r|0;m=c[f>>2]|0;if(m>>>0<(c[h>>2]|0)>>>0){c[f>>2]=m+1;t=d[m]|0}else{t=cn(a)|0}if(!((t+ -48|0)>>>0<10&(s|0)<214748364)){break}q=t;r=s*10|0}r=((s|0)<0)<<31>>31;if((t+ -48|0)>>>0<10){q=s;m=r;j=t;while(1){g=yn(q|0,m|0,10,0)|0;l=I;b=ln(j|0,((j|0)<0)<<31>>31|0,-48,-1)|0;k=ln(b|0,I|0,g|0,l|0)|0;l=I;g=c[f>>2]|0;if(g>>>0<(c[h>>2]|0)>>>0){c[f>>2]=g+1;u=d[g]|0}else{u=cn(a)|0}if((u+ -48|0)>>>0<10&((l|0)<21474836|(l|0)==21474836&k>>>0<2061584302)){q=k;m=l;j=u}else{v=k;w=l;x=u;break}}}else{v=s;w=r;x=t}if((x+ -48|0)>>>0<10){do{x=c[f>>2]|0;if(x>>>0<(c[h>>2]|0)>>>0){c[f>>2]=x+1;y=d[x]|0}else{y=cn(a)|0}}while((y+ -48|0)>>>0<10)}if((c[h>>2]|0)!=0){c[f>>2]=(c[f>>2]|0)+ -1}f=(n|0)!=0;n=kn(0,0,v|0,w|0)|0;o=f?I:w;p=f?n:v;I=o;i=e;return p|0}function bn(a,b){a=a|0;b=b|0;var d=0,e=0,f=0,g=0;d=i;c[a+104>>2]=b;e=c[a+8>>2]|0;f=c[a+4>>2]|0;g=e-f|0;c[a+108>>2]=g;if((b|0)!=0&(g|0)>(b|0)){c[a+100>>2]=f+b;i=d;return}else{c[a+100>>2]=e;i=d;return}}function cn(b){b=b|0;var e=0,f=0,g=0,h=0,j=0,k=0,l=0;e=i;f=b+104|0;g=c[f>>2]|0;if(!((g|0)!=0?(c[b+108>>2]|0)>=(g|0):0)){h=3}if((h|0)==3?(h=gn(b)|0,(h|0)>=0):0){g=c[f>>2]|0;f=c[b+8>>2]|0;if((g|0)!=0?(j=c[b+4>>2]|0,k=g-(c[b+108>>2]|0)+ -1|0,(f-j|0)>(k|0)):0){c[b+100>>2]=j+k}else{c[b+100>>2]=f}k=c[b+4>>2]|0;if((f|0)!=0){j=b+108|0;c[j>>2]=f+1-k+(c[j>>2]|0)}j=k+ -1|0;if((d[j]|0|0)==(h|0)){l=h;i=e;return l|0}a[j]=h;l=h;i=e;return l|0}c[b+100>>2]=0;l=-1;i=e;return l|0}function dn(a,b){a=+a;b=b|0;var d=0,e=0.0,f=0,g=0,j=0,l=0.0;d=i;if((b|0)>1023){e=a*8.98846567431158e+307;f=b+ -1023|0;if((f|0)>1023){g=b+ -2046|0;j=(g|0)>1023?1023:g;l=e*8.98846567431158e+307}else{j=f;l=e}}else{if((b|0)<-1022){e=a*2.2250738585072014e-308;f=b+1022|0;if((f|0)<-1022){g=b+2044|0;j=(g|0)<-1022?-1022:g;l=e*2.2250738585072014e-308}else{j=f;l=e}}else{j=b;l=a}}b=qn(j+1023|0,0,52)|0;j=I;c[k>>2]=b;c[k+4>>2]=j;a=l*+h[k>>3];i=d;return+a}function en(a,b){a=+a;b=b|0;var c=0,d=0.0;c=i;d=+dn(a,b);i=c;return+d}function fn(b){b=b|0;var d=0,e=0,f=0,g=0,h=0;d=i;e=b+74|0;f=a[e]|0;a[e]=f+255|f;f=b+20|0;e=b+44|0;if((c[f>>2]|0)>>>0>(c[e>>2]|0)>>>0){gc[c[b+36>>2]&31](b,0,0)|0}c[b+16>>2]=0;c[b+28>>2]=0;c[f>>2]=0;f=c[b>>2]|0;if((f&20|0)==0){g=c[e>>2]|0;c[b+8>>2]=g;c[b+4>>2]=g;h=0;i=d;return h|0}if((f&4|0)==0){h=-1;i=d;return h|0}c[b>>2]=f|32;h=-1;i=d;return h|0}function gn(a){a=a|0;var b=0,e=0,f=0;b=i;i=i+16|0;e=b;if((c[a+8>>2]|0)==0?(fn(a)|0)!=0:0){f=-1}else{if((gc[c[a+32>>2]&31](a,e,1)|0)==1){f=d[e]|0}else{f=-1}}i=b;return f|0}function hn(a,b,d){a=a|0;b=b|0;d=d|0;var e=0,f=0,g=0,h=0.0,j=0,k=0;d=i;i=i+112|0;e=d;f=e+0|0;g=f+112|0;do{c[f>>2]=0;f=f+4|0}while((f|0)<(g|0));f=e+4|0;c[f>>2]=a;g=e+8|0;c[g>>2]=-1;c[e+44>>2]=a;c[e+76>>2]=-1;bn(e,0);h=+$m(e,2,1);j=(c[f>>2]|0)-(c[g>>2]|0)+(c[e+108>>2]|0)|0;if((b|0)==0){i=d;return+h}if((j|0)==0){k=a}else{k=a+j|0}c[b>>2]=k;i=d;return+h}function jn(){c[1292]=o;c[1318]=o;c[4056]=o;c[4286]=o}function kn(a,b,c,d){a=a|0;b=b|0;c=c|0;d=d|0;var e=0;e=b-d>>>0;e=b-d-(c>>>0>a>>>0|0)>>>0;return(I=e,a-c>>>0|0)|0}function ln(a,b,c,d){a=a|0;b=b|0;c=c|0;d=d|0;var e=0;e=a+c>>>0;return(I=b+d+(e>>>0<a>>>0|0)>>>0,e|0)|0}function mn(b){b=b|0;var c=0;c=b;while(a[c]|0){c=c+1|0}return c-b|0}function nn(b,d,e){b=b|0;d=d|0;e=e|0;var f=0;if((e|0)>=4096)return ob(b|0,d|0,e|0)|0;f=b|0;if((b&3)==(d&3)){while(b&3){if((e|0)==0)return f|0;a[b]=a[d]|0;b=b+1|0;d=d+1|0;e=e-1|0}while((e|0)>=4){c[b>>2]=c[d>>2];b=b+4|0;d=d+4|0;e=e-4|0}}while((e|0)>0){a[b]=a[d]|0;b=b+1|0;d=d+1|0;e=e-1|0}return f|0}function on(b,c,d){b=b|0;c=c|0;d=d|0;var e=0;if((c|0)<(b|0)&(b|0)<(c+d|0)){e=b;c=c+d|0;b=b+d|0;while((d|0)>0){b=b-1|0;c=c-1|0;d=d-1|0;a[b]=a[c]|0}b=e}else{nn(b,c,d)|0}return b|0}function pn(b,d,e){b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,i=0;f=b+e|0;if((e|0)>=20){d=d&255;g=b&3;h=d|d<<8|d<<16|d<<24;i=f&~3;if(g){g=b+4-g|0;while((b|0)<(g|0)){a[b]=d;b=b+1|0}}while((b|0)<(i|0)){c[b>>2]=h;b=b+4|0}}while((b|0)<(f|0)){a[b]=d;b=b+1|0}return b-e|0}function qn(a,b,c){a=a|0;b=b|0;c=c|0;if((c|0)<32){I=b<<c|(a&(1<<c)-1<<32-c)>>>32-c;return a<<c}I=a<<c-32;return 0}function rn(a,b,c){a=a|0;b=b|0;c=c|0;if((c|0)<32){I=b>>>c;return a>>>c|(b&(1<<c)-1)<<32-c}I=0;return b>>>c-32|0}function sn(a,b,c){a=a|0;b=b|0;c=c|0;if((c|0)<32){I=b>>c;return a>>>c|(b&(1<<c)-1)<<32-c}I=(b|0)<0?-1:0;return b>>c-32|0}function tn(b){b=b|0;var c=0;c=a[n+(b>>>24)|0]|0;if((c|0)<8)return c|0;c=a[n+(b>>16&255)|0]|0;if((c|0)<8)return c+8|0;c=a[n+(b>>8&255)|0]|0;if((c|0)<8)return c+16|0;return(a[n+(b&255)|0]|0)+24|0}function un(b){b=b|0;var c=0;c=a[m+(b&255)|0]|0;if((c|0)<8)return c|0;c=a[m+(b>>8&255)|0]|0;if((c|0)<8)return c+8|0;c=a[m+(b>>16&255)|0]|0;if((c|0)<8)return c+16|0;return(a[m+(b>>>24)|0]|0)+24|0}function vn(a,b){a=a|0;b=b|0;var c=0,d=0,e=0,f=0;c=a&65535;d=b&65535;e=ea(d,c)|0;f=a>>>16;a=(e>>>16)+(ea(d,f)|0)|0;d=b>>>16;b=ea(d,c)|0;return(I=(a>>>16)+(ea(d,f)|0)+(((a&65535)+b|0)>>>16)|0,a+b<<16|e&65535|0)|0}function wn(a,b,c,d){a=a|0;b=b|0;c=c|0;d=d|0;var e=0,f=0,g=0,h=0,i=0;e=b>>31|((b|0)<0?-1:0)<<1;f=((b|0)<0?-1:0)>>31|((b|0)<0?-1:0)<<1;g=d>>31|((d|0)<0?-1:0)<<1;h=((d|0)<0?-1:0)>>31|((d|0)<0?-1:0)<<1;i=kn(e^a,f^b,e,f)|0;b=I;a=g^e;e=h^f;f=kn((Bn(i,b,kn(g^c,h^d,g,h)|0,I,0)|0)^a,I^e,a,e)|0;return f|0}function xn(a,b,d,e){a=a|0;b=b|0;d=d|0;e=e|0;var f=0,g=0,h=0,j=0,k=0,l=0,m=0;f=i;i=i+8|0;g=f|0;h=b>>31|((b|0)<0?-1:0)<<1;j=((b|0)<0?-1:0)>>31|((b|0)<0?-1:0)<<1;k=e>>31|((e|0)<0?-1:0)<<1;l=((e|0)<0?-1:0)>>31|((e|0)<0?-1:0)<<1;m=kn(h^a,j^b,h,j)|0;b=I;Bn(m,b,kn(k^d,l^e,k,l)|0,I,g)|0;l=kn(c[g>>2]^h,c[g+4>>2]^j,h,j)|0;j=I;i=f;return(I=j,l)|0}function yn(a,b,c,d){a=a|0;b=b|0;c=c|0;d=d|0;var e=0,f=0;e=a;a=c;c=vn(e,a)|0;f=I;return(I=(ea(b,a)|0)+(ea(d,e)|0)+f|f&0,c|0|0)|0}function zn(a,b,c,d){a=a|0;b=b|0;c=c|0;d=d|0;var e=0;e=Bn(a,b,c,d,0)|0;return e|0}function An(a,b,d,e){a=a|0;b=b|0;d=d|0;e=e|0;var f=0,g=0;f=i;i=i+8|0;g=f|0;Bn(a,b,d,e,g)|0;i=f;return(I=c[g+4>>2]|0,c[g>>2]|0)|0}function Bn(a,b,d,e,f){a=a|0;b=b|0;d=d|0;e=e|0;f=f|0;var g=0,h=0,i=0,j=0,k=0,l=0,m=0,n=0,o=0,p=0,q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,A=0,B=0,C=0,D=0,E=0,F=0,G=0,H=0,J=0,K=0,L=0,M=0;g=a;h=b;i=h;j=d;k=e;l=k;if((i|0)==0){m=(f|0)!=0;if((l|0)==0){if(m){c[f>>2]=(g>>>0)%(j>>>0);c[f+4>>2]=0}n=0;o=(g>>>0)/(j>>>0)>>>0;return(I=n,o)|0}else{if(!m){n=0;o=0;return(I=n,o)|0}c[f>>2]=a|0;c[f+4>>2]=b&0;n=0;o=0;return(I=n,o)|0}}m=(l|0)==0;do{if((j|0)!=0){if(!m){p=(tn(l|0)|0)-(tn(i|0)|0)|0;if(p>>>0<=31){q=p+1|0;r=31-p|0;s=p-31>>31;t=q;u=g>>>(q>>>0)&s|i<<r;v=i>>>(q>>>0)&s;w=0;x=g<<r;break}if((f|0)==0){n=0;o=0;return(I=n,o)|0}c[f>>2]=a|0;c[f+4>>2]=h|b&0;n=0;o=0;return(I=n,o)|0}r=j-1|0;if((r&j|0)!=0){s=(tn(j|0)|0)+33-(tn(i|0)|0)|0;q=64-s|0;p=32-s|0;y=p>>31;z=s-32|0;A=z>>31;t=s;u=p-1>>31&i>>>(z>>>0)|(i<<p|g>>>(s>>>0))&A;v=A&i>>>(s>>>0);w=g<<q&y;x=(i<<q|g>>>(z>>>0))&y|g<<p&s-33>>31;break}if((f|0)!=0){c[f>>2]=r&g;c[f+4>>2]=0}if((j|0)==1){n=h|b&0;o=a|0|0;return(I=n,o)|0}else{r=un(j|0)|0;n=i>>>(r>>>0)|0;o=i<<32-r|g>>>(r>>>0)|0;return(I=n,o)|0}}else{if(m){if((f|0)!=0){c[f>>2]=(i>>>0)%(j>>>0);c[f+4>>2]=0}n=0;o=(i>>>0)/(j>>>0)>>>0;return(I=n,o)|0}if((g|0)==0){if((f|0)!=0){c[f>>2]=0;c[f+4>>2]=(i>>>0)%(l>>>0)}n=0;o=(i>>>0)/(l>>>0)>>>0;return(I=n,o)|0}r=l-1|0;if((r&l|0)==0){if((f|0)!=0){c[f>>2]=a|0;c[f+4>>2]=r&i|b&0}n=0;o=i>>>((un(l|0)|0)>>>0);return(I=n,o)|0}r=(tn(l|0)|0)-(tn(i|0)|0)|0;if(r>>>0<=30){s=r+1|0;p=31-r|0;t=s;u=i<<p|g>>>(s>>>0);v=i>>>(s>>>0);w=0;x=g<<p;break}if((f|0)==0){n=0;o=0;return(I=n,o)|0}c[f>>2]=a|0;c[f+4>>2]=h|b&0;n=0;o=0;return(I=n,o)|0}}while(0);if((t|0)==0){B=x;C=w;D=v;E=u;F=0;G=0}else{b=d|0|0;d=k|e&0;e=ln(b,d,-1,-1)|0;k=I;h=x;x=w;w=v;v=u;u=t;t=0;while(1){H=x>>>31|h<<1;J=t|x<<1;a=v<<1|h>>>31|0;g=v>>>31|w<<1|0;kn(e,k,a,g)|0;i=I;l=i>>31|((i|0)<0?-1:0)<<1;K=l&1;L=kn(a,g,l&b,(((i|0)<0?-1:0)>>31|((i|0)<0?-1:0)<<1)&d)|0;M=I;i=u-1|0;if((i|0)==0){break}else{h=H;x=J;w=M;v=L;u=i;t=K}}B=H;C=J;D=M;E=L;F=0;G=K}K=C;C=0;if((f|0)!=0){c[f>>2]=E;c[f+4>>2]=D}n=(K|0)>>>31|(B|C)<<1|(C<<1|K>>>31)&0|F;o=(K<<1|0>>>31)&-2|G;return(I=n,o)|0}function Cn(a,b,c,d){a=a|0;b=b|0;c=c|0;d=d|0;return gc[a&31](b|0,c|0,d|0)|0}function Dn(a,b,c,d,e,f,g,h){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;hc[a&63](b|0,c|0,d|0,e|0,f|0,g|0,h|0)}function En(a,b,c,d,e,f){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;ic[a&3](b|0,c|0,d|0,e|0,f|0)}function Fn(a,b){a=a|0;b=b|0;jc[a&127](b|0)}function Gn(a,b,c){a=a|0;b=b|0;c=c|0;kc[a&63](b|0,c|0)}function Hn(a,b,c,d,e,f,g,h,i,j){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;i=i|0;j=j|0;lc[a&3](b|0,c|0,d|0,e|0,f|0,g|0,h|0,i|0,j|0)}function In(a,b){a=a|0;b=b|0;return mc[a&63](b|0)|0}function Jn(a,b,c,d,e,f,g,h){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;g=g|0;h=+h;nc[a&3](b|0,c|0,d|0,e|0,f|0,g|0,+h)}function Kn(a,b,c,d){a=a|0;b=b|0;c=c|0;d=d|0;oc[a&3](b|0,c|0,d|0)}function Ln(a,b,c,d,e,f,g){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;g=+g;pc[a&7](b|0,c|0,d|0,e|0,f|0,+g)}function Mn(a){a=a|0;qc[a&0]()}function Nn(a,b,c,d,e,f,g,h,i){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;i=i|0;return rc[a&15](b|0,c|0,d|0,e|0,f|0,g|0,h|0,i|0)|0}function On(a,b,c,d,e){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;return sc[a&7](b|0,c|0,d|0,e|0)|0}function Pn(a,b,c,d,e,f,g,h,i){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;i=i|0;tc[a&7](b|0,c|0,d|0,e|0,f|0,g|0,h|0,i|0)}function Qn(a,b,c,d,e,f,g){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;g=g|0;uc[a&15](b|0,c|0,d|0,e|0,f|0,g|0)}function Rn(a,b,c){a=a|0;b=b|0;c=c|0;return vc[a&15](b|0,c|0)|0}function Sn(a,b,c,d,e,f){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;return wc[a&15](b|0,c|0,d|0,e|0,f|0)|0}function Tn(a,b,c,d,e){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;xc[a&7](b|0,c|0,d|0,e|0)}function Un(a,b,c){a=a|0;b=b|0;c=c|0;fa(0);return 0}function Vn(a,b,c,d,e,f,g){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;g=g|0;fa(1)}function Wn(a,b,c,d,e){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;fa(2)}function Xn(a){a=a|0;fa(3)}function Yn(a,b){a=a|0;b=b|0;fa(4)}function Zn(a,b,c,d,e,f,g,h,i){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;i=i|0;fa(5)}function _n(a){a=a|0;fa(6);return 0}function $n(a,b,c,d,e,f,g){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;g=+g;fa(7)}function ao(a,b,c){a=a|0;b=b|0;c=c|0;fa(8)}function bo(a,b,c,d,e,f){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=+f;fa(9)}function co(){fa(10)}function eo(a,b,c,d,e,f,g,h){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;fa(11);return 0}function fo(a,b,c,d){a=a|0;b=b|0;c=c|0;d=d|0;fa(12);return 0}function go(a,b,c,d,e,f,g,h){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;g=g|0;h=h|0;fa(13)}function ho(a,b,c,d,e,f){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;f=f|0;fa(14)}function io(a,b){a=a|0;b=b|0;fa(15);return 0}function jo(a,b,c,d,e){a=a|0;b=b|0;c=c|0;d=d|0;e=e|0;fa(16);return 0}function ko(a,b,c,d){a=a|0;b=b|0;c=c|0;d=d|0;fa(17)}




// EMSCRIPTEN_END_FUNCS
var gc=[Un,lf,qf,Bd,uf,Ye,bf,Pd,ff,ne,oe,ag,fg,Jj,Oj,tk,vk,yk,ek,jk,lk,ok,Am,Un,Un,Un,Un,Un,Un,Un,Un,Un];var hc=[Vn,ig,kg,mg,og,qg,sg,ug,wg,yg,Ag,Cg,Hg,Jg,Lg,Ng,Pg,Rg,Tg,Vg,Xg,Zg,$g,oh,qh,Ch,Eh,Nh,Oh,Ph,Qh,Rh,_h,$h,ai,bi,ci,Aj,Gj,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn,Vn];var ic=[Wn,Hm,Gm,Fm];var jc=[Xn,xd,yd,Ed,Fd,Ld,Md,Sd,Td,de,ce,ie,he,ke,te,se,We,Ve,jf,hf,xf,wf,zf,yf,Cf,Bf,Ef,Df,Hf,Gf,Jf,If,Mf,Lf,Of,Nf,Uf,Tf,Se,Vf,Sf,Wf,Yf,Xf,bk,cg,bg,hg,gg,Gg,Fg,ih,hh,xh,wh,Lh,Kh,Yh,Xh,ii,hi,li,ki,pi,oi,Ai,zi,Li,Ki,Wi,Vi,fj,ej,pj,oj,wj,vj,Cj,Bj,Ij,Hj,Nj,Mj,Wj,Vj,rk,qk,Rj,Ik,nl,ml,pl,ol,Zf,ak,dk,Ak,Qk,$k,kl,ll,sm,rm,um,xm,vm,wm,ym,zm,Ym,Xm,wd,ck,Xl,hj,Mm,cm,bm,am,$l,_l,Zl,ze,Ke,Xn,Xn,Xn,Xn];var kc=[Yn,zd,Gd,Nd,Ud,Xe,kf,si,ti,ui,vi,xi,yi,Di,Ei,Fi,Gi,Ii,Ji,Oi,Pi,Qi,Ri,Ti,Ui,Zi,_i,$i,aj,cj,dj,Lj,Qj,ul,wl,yl,vl,xl,zl,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn,Yn];var lc=[Zn,Sh,di,Zn];var mc=[_n,Ad,pf,rf,sf,of,Hd,Id,Od,af,cf,df,$e,Vd,Wd,ee,je,Qf,Mh,Al,Cl,El,Kl,Ml,Gl,Il,Zh,Bl,Dl,Fl,Ll,Nl,Hl,Jl,qi,ri,wi,Bi,Ci,Hi,Mi,Ni,Si,Xi,Yi,bj,Mk,Nk,Pk,ql,sl,rl,tl,Ek,Fk,Hk,Wk,Xk,_k,fl,gl,jl,tm,Zm];var nc=[$n,xj,Dj,$n];var oc=[ao,me,Rf,ao];var pc=[bo,rh,uh,Fh,Hh,bo,bo,bo];var qc=[co];var rc=[eo,Jk,Kk,Bk,Ck,Rk,Tk,al,cl,eo,eo,eo,eo,eo,eo,eo];var sc=[fo,xk,fk,gk,hk,nk,fo,fo];var tc=[go,ji,mi,gj,kj,qj,sj,go];var uc=[ho,mf,Ze,jh,kh,ph,vh,yh,zh,Dh,Ih,Kj,Pj,Km,Jm,Im];var vc=[io,tf,Cd,Jd,vf,ef,Qd,Xd,gf,sk,uk,wk,ik,kk,mk,io];var wc=[jo,_f,dg,zk,Lk,Ok,pk,Dk,Gk,Vk,Yk,el,hl,jo,jo,jo];var xc=[ko,nf,_e,$f,eg,Bm,Cm,Dm];return{_i64Subtract:kn,_free:Mm,_realloc:Nm,_i64Add:ln,_memmove:on,_strlen:mn,_process:Oc,_memset:pn,_malloc:Lm,_memcpy:nn,_bitshift64Shl:qn,__GLOBAL__I_a:Zd,runPostSets:jn,stackAlloc:yc,stackSave:zc,stackRestore:Ac,setThrew:Bc,setTempRet0:Ec,setTempRet1:Fc,setTempRet2:Gc,setTempRet3:Hc,setTempRet4:Ic,setTempRet5:Jc,setTempRet6:Kc,setTempRet7:Lc,setTempRet8:Mc,setTempRet9:Nc,dynCall_iiii:Cn,dynCall_viiiiiii:Dn,dynCall_viiiii:En,dynCall_vi:Fn,dynCall_vii:Gn,dynCall_viiiiiiiii:Hn,dynCall_ii:In,dynCall_viiiiiid:Jn,dynCall_viii:Kn,dynCall_viiiiid:Ln,dynCall_v:Mn,dynCall_iiiiiiiii:Nn,dynCall_iiiii:On,dynCall_viiiiiiii:Pn,dynCall_viiiiii:Qn,dynCall_iii:Rn,dynCall_iiiiii:Sn,dynCall_viiii:Tn}})


// EMSCRIPTEN_END_ASM
({ "Math": Math, "Int8Array": Int8Array, "Int16Array": Int16Array, "Int32Array": Int32Array, "Uint8Array": Uint8Array, "Uint16Array": Uint16Array, "Uint32Array": Uint32Array, "Float32Array": Float32Array, "Float64Array": Float64Array }, { "abort": abort, "assert": assert, "asmPrintInt": asmPrintInt, "asmPrintFloat": asmPrintFloat, "min": Math_min, "invoke_iiii": invoke_iiii, "invoke_viiiiiii": invoke_viiiiiii, "invoke_viiiii": invoke_viiiii, "invoke_vi": invoke_vi, "invoke_vii": invoke_vii, "invoke_viiiiiiiii": invoke_viiiiiiiii, "invoke_ii": invoke_ii, "invoke_viiiiiid": invoke_viiiiiid, "invoke_viii": invoke_viii, "invoke_viiiiid": invoke_viiiiid, "invoke_v": invoke_v, "invoke_iiiiiiiii": invoke_iiiiiiiii, "invoke_iiiii": invoke_iiiii, "invoke_viiiiiiii": invoke_viiiiiiii, "invoke_viiiiii": invoke_viiiiii, "invoke_iii": invoke_iii, "invoke_iiiiii": invoke_iiiiii, "invoke_viiii": invoke_viiii, "_fabs": _fabs, "_pthread_cond_wait": _pthread_cond_wait, "_freelocale": _freelocale, "__formatString": __formatString, "_asprintf": _asprintf, "__reallyNegative": __reallyNegative, "_send": _send, "_strtoll_l": _strtoll_l, "_vsscanf": _vsscanf, "___ctype_b_loc": ___ctype_b_loc, "__ZSt9terminatev": __ZSt9terminatev, "_fmod": _fmod, "___cxa_guard_acquire": ___cxa_guard_acquire, "___setErrNo": ___setErrNo, "___cxa_is_number_type": ___cxa_is_number_type, "_atexit": _atexit, "_copysign": _copysign, "_emscripten_asm_const_int": _emscripten_asm_const_int, "___assert_fail": ___assert_fail, "_ungetc": _ungetc, "___cxa_allocate_exception": ___cxa_allocate_exception, "__ZSt18uncaught_exceptionv": __ZSt18uncaught_exceptionv, "_isxdigit_l": _isxdigit_l, "___ctype_toupper_loc": ___ctype_toupper_loc, "_fflush": _fflush, "___cxa_guard_release": ___cxa_guard_release, "__addDays": __addDays, "___errno_location": ___errno_location, "_strtoll": _strtoll, "_strerror_r": _strerror_r, "_emscripten_asm_const_double": _emscripten_asm_const_double, "_strftime_l": _strftime_l, "_sscanf": _sscanf, "_fabsf": _fabsf, "_sbrk": _sbrk, "_uselocale": _uselocale, "_catgets": _catgets, "_newlocale": _newlocale, "_snprintf": _snprintf, "___cxa_begin_catch": ___cxa_begin_catch, "_emscripten_memcpy_big": _emscripten_memcpy_big, "_fileno": _fileno, "___resumeException": ___resumeException, "___cxa_find_matching_catch": ___cxa_find_matching_catch, "_sysconf": _sysconf, "_strtoull": _strtoull, "_isdigit_l": _isdigit_l, "_strftime": _strftime, "__arraySum": __arraySum, "___cxa_throw": ___cxa_throw, "___ctype_tolower_loc": ___ctype_tolower_loc, "_isspace": _isspace, "_pthread_mutex_unlock": _pthread_mutex_unlock, "_fread": _fread, "_pthread_cond_broadcast": _pthread_cond_broadcast, "_isxdigit": _isxdigit, "_sprintf": _sprintf, "_sqrtf": _sqrtf, "_vasprintf": _vasprintf, "_write": _write, "__isLeapYear": __isLeapYear, "__scanString": __scanString, "_recv": _recv, "_vsnprintf": _vsnprintf, "__ZNSt9exceptionD2Ev": __ZNSt9exceptionD2Ev, "_fgetc": _fgetc, "_strtoull_l": _strtoull_l, "_mkport": _mkport, "___cxa_does_inherit": ___cxa_does_inherit, "__exit": __exit, "_read": _read, "__parseInt64": __parseInt64, "_abort": _abort, "_catclose": _catclose, "_fwrite": _fwrite, "_time": _time, "_pthread_mutex_lock": _pthread_mutex_lock, "_strerror": _strerror, "_pread": _pread, "_pwrite": _pwrite, "_sqrt": _sqrt, "_catopen": _catopen, "_exit": _exit, "_isdigit": _isdigit, "__getFloat": __getFloat, "STACKTOP": STACKTOP, "STACK_MAX": STACK_MAX, "tempDoublePtr": tempDoublePtr, "ABORT": ABORT, "cttz_i8": cttz_i8, "ctlz_i8": ctlz_i8, "NaN": NaN, "Infinity": Infinity, "__ZTISt9exception": __ZTISt9exception, "___dso_handle": ___dso_handle, "_stderr": _stderr, "_stdin": _stdin, "_stdout": _stdout }, buffer);
var _i64Subtract = Module["_i64Subtract"] = asm["_i64Subtract"];
var _free = Module["_free"] = asm["_free"];
var _realloc = Module["_realloc"] = asm["_realloc"];
var _i64Add = Module["_i64Add"] = asm["_i64Add"];
var _memmove = Module["_memmove"] = asm["_memmove"];
var _strlen = Module["_strlen"] = asm["_strlen"];
var _process = Module["_process"] = asm["_process"];
var _memset = Module["_memset"] = asm["_memset"];
var _malloc = Module["_malloc"] = asm["_malloc"];
var _memcpy = Module["_memcpy"] = asm["_memcpy"];
var _bitshift64Shl = Module["_bitshift64Shl"] = asm["_bitshift64Shl"];
var __GLOBAL__I_a = Module["__GLOBAL__I_a"] = asm["__GLOBAL__I_a"];
var runPostSets = Module["runPostSets"] = asm["runPostSets"];
var dynCall_iiii = Module["dynCall_iiii"] = asm["dynCall_iiii"];
var dynCall_viiiiiii = Module["dynCall_viiiiiii"] = asm["dynCall_viiiiiii"];
var dynCall_viiiii = Module["dynCall_viiiii"] = asm["dynCall_viiiii"];
var dynCall_vi = Module["dynCall_vi"] = asm["dynCall_vi"];
var dynCall_vii = Module["dynCall_vii"] = asm["dynCall_vii"];
var dynCall_viiiiiiiii = Module["dynCall_viiiiiiiii"] = asm["dynCall_viiiiiiiii"];
var dynCall_ii = Module["dynCall_ii"] = asm["dynCall_ii"];
var dynCall_viiiiiid = Module["dynCall_viiiiiid"] = asm["dynCall_viiiiiid"];
var dynCall_viii = Module["dynCall_viii"] = asm["dynCall_viii"];
var dynCall_viiiiid = Module["dynCall_viiiiid"] = asm["dynCall_viiiiid"];
var dynCall_v = Module["dynCall_v"] = asm["dynCall_v"];
var dynCall_iiiiiiiii = Module["dynCall_iiiiiiiii"] = asm["dynCall_iiiiiiiii"];
var dynCall_iiiii = Module["dynCall_iiiii"] = asm["dynCall_iiiii"];
var dynCall_viiiiiiii = Module["dynCall_viiiiiiii"] = asm["dynCall_viiiiiiii"];
var dynCall_viiiiii = Module["dynCall_viiiiii"] = asm["dynCall_viiiiii"];
var dynCall_iii = Module["dynCall_iii"] = asm["dynCall_iii"];
var dynCall_iiiiii = Module["dynCall_iiiiii"] = asm["dynCall_iiiiii"];
var dynCall_viiii = Module["dynCall_viiii"] = asm["dynCall_viiii"];

Runtime.stackAlloc = function(size) { return asm['stackAlloc'](size) };
Runtime.stackSave = function() { return asm['stackSave']() };
Runtime.stackRestore = function(top) { asm['stackRestore'](top) };


// TODO: strip out parts of this we do not need

//======= begin closure i64 code =======

// Copyright 2009 The Closure Library Authors. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS-IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

/**
 * @fileoverview Defines a Long class for representing a 64-bit two's-complement
 * integer value, which faithfully simulates the behavior of a Java "long". This
 * implementation is derived from LongLib in GWT.
 *
 */

var i64Math = (function() { // Emscripten wrapper
  var goog = { math: {} };


  /**
   * Constructs a 64-bit two's-complement integer, given its low and high 32-bit
   * values as *signed* integers.  See the from* functions below for more
   * convenient ways of constructing Longs.
   *
   * The internal representation of a long is the two given signed, 32-bit values.
   * We use 32-bit pieces because these are the size of integers on which
   * Javascript performs bit-operations.  For operations like addition and
   * multiplication, we split each number into 16-bit pieces, which can easily be
   * multiplied within Javascript's floating-point representation without overflow
   * or change in sign.
   *
   * In the algorithms below, we frequently reduce the negative case to the
   * positive case by negating the input(s) and then post-processing the result.
   * Note that we must ALWAYS check specially whether those values are MIN_VALUE
   * (-2^63) because -MIN_VALUE == MIN_VALUE (since 2^63 cannot be represented as
   * a positive number, it overflows back into a negative).  Not handling this
   * case would often result in infinite recursion.
   *
   * @param {number} low  The low (signed) 32 bits of the long.
   * @param {number} high  The high (signed) 32 bits of the long.
   * @constructor
   */
  goog.math.Long = function(low, high) {
    /**
     * @type {number}
     * @private
     */
    this.low_ = low | 0;  // force into 32 signed bits.

    /**
     * @type {number}
     * @private
     */
    this.high_ = high | 0;  // force into 32 signed bits.
  };


  // NOTE: Common constant values ZERO, ONE, NEG_ONE, etc. are defined below the
  // from* methods on which they depend.


  /**
   * A cache of the Long representations of small integer values.
   * @type {!Object}
   * @private
   */
  goog.math.Long.IntCache_ = {};


  /**
   * Returns a Long representing the given (32-bit) integer value.
   * @param {number} value The 32-bit integer in question.
   * @return {!goog.math.Long} The corresponding Long value.
   */
  goog.math.Long.fromInt = function(value) {
    if (-128 <= value && value < 128) {
      var cachedObj = goog.math.Long.IntCache_[value];
      if (cachedObj) {
        return cachedObj;
      }
    }

    var obj = new goog.math.Long(value | 0, value < 0 ? -1 : 0);
    if (-128 <= value && value < 128) {
      goog.math.Long.IntCache_[value] = obj;
    }
    return obj;
  };


  /**
   * Returns a Long representing the given value, provided that it is a finite
   * number.  Otherwise, zero is returned.
   * @param {number} value The number in question.
   * @return {!goog.math.Long} The corresponding Long value.
   */
  goog.math.Long.fromNumber = function(value) {
    if (isNaN(value) || !isFinite(value)) {
      return goog.math.Long.ZERO;
    } else if (value <= -goog.math.Long.TWO_PWR_63_DBL_) {
      return goog.math.Long.MIN_VALUE;
    } else if (value + 1 >= goog.math.Long.TWO_PWR_63_DBL_) {
      return goog.math.Long.MAX_VALUE;
    } else if (value < 0) {
      return goog.math.Long.fromNumber(-value).negate();
    } else {
      return new goog.math.Long(
          (value % goog.math.Long.TWO_PWR_32_DBL_) | 0,
          (value / goog.math.Long.TWO_PWR_32_DBL_) | 0);
    }
  };


  /**
   * Returns a Long representing the 64-bit integer that comes by concatenating
   * the given high and low bits.  Each is assumed to use 32 bits.
   * @param {number} lowBits The low 32-bits.
   * @param {number} highBits The high 32-bits.
   * @return {!goog.math.Long} The corresponding Long value.
   */
  goog.math.Long.fromBits = function(lowBits, highBits) {
    return new goog.math.Long(lowBits, highBits);
  };


  /**
   * Returns a Long representation of the given string, written using the given
   * radix.
   * @param {string} str The textual representation of the Long.
   * @param {number=} opt_radix The radix in which the text is written.
   * @return {!goog.math.Long} The corresponding Long value.
   */
  goog.math.Long.fromString = function(str, opt_radix) {
    if (str.length == 0) {
      throw Error('number format error: empty string');
    }

    var radix = opt_radix || 10;
    if (radix < 2 || 36 < radix) {
      throw Error('radix out of range: ' + radix);
    }

    if (str.charAt(0) == '-') {
      return goog.math.Long.fromString(str.substring(1), radix).negate();
    } else if (str.indexOf('-') >= 0) {
      throw Error('number format error: interior "-" character: ' + str);
    }

    // Do several (8) digits each time through the loop, so as to
    // minimize the calls to the very expensive emulated div.
    var radixToPower = goog.math.Long.fromNumber(Math.pow(radix, 8));

    var result = goog.math.Long.ZERO;
    for (var i = 0; i < str.length; i += 8) {
      var size = Math.min(8, str.length - i);
      var value = parseInt(str.substring(i, i + size), radix);
      if (size < 8) {
        var power = goog.math.Long.fromNumber(Math.pow(radix, size));
        result = result.multiply(power).add(goog.math.Long.fromNumber(value));
      } else {
        result = result.multiply(radixToPower);
        result = result.add(goog.math.Long.fromNumber(value));
      }
    }
    return result;
  };


  // NOTE: the compiler should inline these constant values below and then remove
  // these variables, so there should be no runtime penalty for these.


  /**
   * Number used repeated below in calculations.  This must appear before the
   * first call to any from* function below.
   * @type {number}
   * @private
   */
  goog.math.Long.TWO_PWR_16_DBL_ = 1 << 16;


  /**
   * @type {number}
   * @private
   */
  goog.math.Long.TWO_PWR_24_DBL_ = 1 << 24;


  /**
   * @type {number}
   * @private
   */
  goog.math.Long.TWO_PWR_32_DBL_ =
      goog.math.Long.TWO_PWR_16_DBL_ * goog.math.Long.TWO_PWR_16_DBL_;


  /**
   * @type {number}
   * @private
   */
  goog.math.Long.TWO_PWR_31_DBL_ =
      goog.math.Long.TWO_PWR_32_DBL_ / 2;


  /**
   * @type {number}
   * @private
   */
  goog.math.Long.TWO_PWR_48_DBL_ =
      goog.math.Long.TWO_PWR_32_DBL_ * goog.math.Long.TWO_PWR_16_DBL_;


  /**
   * @type {number}
   * @private
   */
  goog.math.Long.TWO_PWR_64_DBL_ =
      goog.math.Long.TWO_PWR_32_DBL_ * goog.math.Long.TWO_PWR_32_DBL_;


  /**
   * @type {number}
   * @private
   */
  goog.math.Long.TWO_PWR_63_DBL_ =
      goog.math.Long.TWO_PWR_64_DBL_ / 2;


  /** @type {!goog.math.Long} */
  goog.math.Long.ZERO = goog.math.Long.fromInt(0);


  /** @type {!goog.math.Long} */
  goog.math.Long.ONE = goog.math.Long.fromInt(1);


  /** @type {!goog.math.Long} */
  goog.math.Long.NEG_ONE = goog.math.Long.fromInt(-1);


  /** @type {!goog.math.Long} */
  goog.math.Long.MAX_VALUE =
      goog.math.Long.fromBits(0xFFFFFFFF | 0, 0x7FFFFFFF | 0);


  /** @type {!goog.math.Long} */
  goog.math.Long.MIN_VALUE = goog.math.Long.fromBits(0, 0x80000000 | 0);


  /**
   * @type {!goog.math.Long}
   * @private
   */
  goog.math.Long.TWO_PWR_24_ = goog.math.Long.fromInt(1 << 24);


  /** @return {number} The value, assuming it is a 32-bit integer. */
  goog.math.Long.prototype.toInt = function() {
    return this.low_;
  };


  /** @return {number} The closest floating-point representation to this value. */
  goog.math.Long.prototype.toNumber = function() {
    return this.high_ * goog.math.Long.TWO_PWR_32_DBL_ +
           this.getLowBitsUnsigned();
  };


  /**
   * @param {number=} opt_radix The radix in which the text should be written.
   * @return {string} The textual representation of this value.
   */
  goog.math.Long.prototype.toString = function(opt_radix) {
    var radix = opt_radix || 10;
    if (radix < 2 || 36 < radix) {
      throw Error('radix out of range: ' + radix);
    }

    if (this.isZero()) {
      return '0';
    }

    if (this.isNegative()) {
      if (this.equals(goog.math.Long.MIN_VALUE)) {
        // We need to change the Long value before it can be negated, so we remove
        // the bottom-most digit in this base and then recurse to do the rest.
        var radixLong = goog.math.Long.fromNumber(radix);
        var div = this.div(radixLong);
        var rem = div.multiply(radixLong).subtract(this);
        return div.toString(radix) + rem.toInt().toString(radix);
      } else {
        return '-' + this.negate().toString(radix);
      }
    }

    // Do several (6) digits each time through the loop, so as to
    // minimize the calls to the very expensive emulated div.
    var radixToPower = goog.math.Long.fromNumber(Math.pow(radix, 6));

    var rem = this;
    var result = '';
    while (true) {
      var remDiv = rem.div(radixToPower);
      var intval = rem.subtract(remDiv.multiply(radixToPower)).toInt();
      var digits = intval.toString(radix);

      rem = remDiv;
      if (rem.isZero()) {
        return digits + result;
      } else {
        while (digits.length < 6) {
          digits = '0' + digits;
        }
        result = '' + digits + result;
      }
    }
  };


  /** @return {number} The high 32-bits as a signed value. */
  goog.math.Long.prototype.getHighBits = function() {
    return this.high_;
  };


  /** @return {number} The low 32-bits as a signed value. */
  goog.math.Long.prototype.getLowBits = function() {
    return this.low_;
  };


  /** @return {number} The low 32-bits as an unsigned value. */
  goog.math.Long.prototype.getLowBitsUnsigned = function() {
    return (this.low_ >= 0) ?
        this.low_ : goog.math.Long.TWO_PWR_32_DBL_ + this.low_;
  };


  /**
   * @return {number} Returns the number of bits needed to represent the absolute
   *     value of this Long.
   */
  goog.math.Long.prototype.getNumBitsAbs = function() {
    if (this.isNegative()) {
      if (this.equals(goog.math.Long.MIN_VALUE)) {
        return 64;
      } else {
        return this.negate().getNumBitsAbs();
      }
    } else {
      var val = this.high_ != 0 ? this.high_ : this.low_;
      for (var bit = 31; bit > 0; bit--) {
        if ((val & (1 << bit)) != 0) {
          break;
        }
      }
      return this.high_ != 0 ? bit + 33 : bit + 1;
    }
  };


  /** @return {boolean} Whether this value is zero. */
  goog.math.Long.prototype.isZero = function() {
    return this.high_ == 0 && this.low_ == 0;
  };


  /** @return {boolean} Whether this value is negative. */
  goog.math.Long.prototype.isNegative = function() {
    return this.high_ < 0;
  };


  /** @return {boolean} Whether this value is odd. */
  goog.math.Long.prototype.isOdd = function() {
    return (this.low_ & 1) == 1;
  };


  /**
   * @param {goog.math.Long} other Long to compare against.
   * @return {boolean} Whether this Long equals the other.
   */
  goog.math.Long.prototype.equals = function(other) {
    return (this.high_ == other.high_) && (this.low_ == other.low_);
  };


  /**
   * @param {goog.math.Long} other Long to compare against.
   * @return {boolean} Whether this Long does not equal the other.
   */
  goog.math.Long.prototype.notEquals = function(other) {
    return (this.high_ != other.high_) || (this.low_ != other.low_);
  };


  /**
   * @param {goog.math.Long} other Long to compare against.
   * @return {boolean} Whether this Long is less than the other.
   */
  goog.math.Long.prototype.lessThan = function(other) {
    return this.compare(other) < 0;
  };


  /**
   * @param {goog.math.Long} other Long to compare against.
   * @return {boolean} Whether this Long is less than or equal to the other.
   */
  goog.math.Long.prototype.lessThanOrEqual = function(other) {
    return this.compare(other) <= 0;
  };


  /**
   * @param {goog.math.Long} other Long to compare against.
   * @return {boolean} Whether this Long is greater than the other.
   */
  goog.math.Long.prototype.greaterThan = function(other) {
    return this.compare(other) > 0;
  };


  /**
   * @param {goog.math.Long} other Long to compare against.
   * @return {boolean} Whether this Long is greater than or equal to the other.
   */
  goog.math.Long.prototype.greaterThanOrEqual = function(other) {
    return this.compare(other) >= 0;
  };


  /**
   * Compares this Long with the given one.
   * @param {goog.math.Long} other Long to compare against.
   * @return {number} 0 if they are the same, 1 if the this is greater, and -1
   *     if the given one is greater.
   */
  goog.math.Long.prototype.compare = function(other) {
    if (this.equals(other)) {
      return 0;
    }

    var thisNeg = this.isNegative();
    var otherNeg = other.isNegative();
    if (thisNeg && !otherNeg) {
      return -1;
    }
    if (!thisNeg && otherNeg) {
      return 1;
    }

    // at this point, the signs are the same, so subtraction will not overflow
    if (this.subtract(other).isNegative()) {
      return -1;
    } else {
      return 1;
    }
  };


  /** @return {!goog.math.Long} The negation of this value. */
  goog.math.Long.prototype.negate = function() {
    if (this.equals(goog.math.Long.MIN_VALUE)) {
      return goog.math.Long.MIN_VALUE;
    } else {
      return this.not().add(goog.math.Long.ONE);
    }
  };


  /**
   * Returns the sum of this and the given Long.
   * @param {goog.math.Long} other Long to add to this one.
   * @return {!goog.math.Long} The sum of this and the given Long.
   */
  goog.math.Long.prototype.add = function(other) {
    // Divide each number into 4 chunks of 16 bits, and then sum the chunks.

    var a48 = this.high_ >>> 16;
    var a32 = this.high_ & 0xFFFF;
    var a16 = this.low_ >>> 16;
    var a00 = this.low_ & 0xFFFF;

    var b48 = other.high_ >>> 16;
    var b32 = other.high_ & 0xFFFF;
    var b16 = other.low_ >>> 16;
    var b00 = other.low_ & 0xFFFF;

    var c48 = 0, c32 = 0, c16 = 0, c00 = 0;
    c00 += a00 + b00;
    c16 += c00 >>> 16;
    c00 &= 0xFFFF;
    c16 += a16 + b16;
    c32 += c16 >>> 16;
    c16 &= 0xFFFF;
    c32 += a32 + b32;
    c48 += c32 >>> 16;
    c32 &= 0xFFFF;
    c48 += a48 + b48;
    c48 &= 0xFFFF;
    return goog.math.Long.fromBits((c16 << 16) | c00, (c48 << 16) | c32);
  };


  /**
   * Returns the difference of this and the given Long.
   * @param {goog.math.Long} other Long to subtract from this.
   * @return {!goog.math.Long} The difference of this and the given Long.
   */
  goog.math.Long.prototype.subtract = function(other) {
    return this.add(other.negate());
  };


  /**
   * Returns the product of this and the given long.
   * @param {goog.math.Long} other Long to multiply with this.
   * @return {!goog.math.Long} The product of this and the other.
   */
  goog.math.Long.prototype.multiply = function(other) {
    if (this.isZero()) {
      return goog.math.Long.ZERO;
    } else if (other.isZero()) {
      return goog.math.Long.ZERO;
    }

    if (this.equals(goog.math.Long.MIN_VALUE)) {
      return other.isOdd() ? goog.math.Long.MIN_VALUE : goog.math.Long.ZERO;
    } else if (other.equals(goog.math.Long.MIN_VALUE)) {
      return this.isOdd() ? goog.math.Long.MIN_VALUE : goog.math.Long.ZERO;
    }

    if (this.isNegative()) {
      if (other.isNegative()) {
        return this.negate().multiply(other.negate());
      } else {
        return this.negate().multiply(other).negate();
      }
    } else if (other.isNegative()) {
      return this.multiply(other.negate()).negate();
    }

    // If both longs are small, use float multiplication
    if (this.lessThan(goog.math.Long.TWO_PWR_24_) &&
        other.lessThan(goog.math.Long.TWO_PWR_24_)) {
      return goog.math.Long.fromNumber(this.toNumber() * other.toNumber());
    }

    // Divide each long into 4 chunks of 16 bits, and then add up 4x4 products.
    // We can skip products that would overflow.

    var a48 = this.high_ >>> 16;
    var a32 = this.high_ & 0xFFFF;
    var a16 = this.low_ >>> 16;
    var a00 = this.low_ & 0xFFFF;

    var b48 = other.high_ >>> 16;
    var b32 = other.high_ & 0xFFFF;
    var b16 = other.low_ >>> 16;
    var b00 = other.low_ & 0xFFFF;

    var c48 = 0, c32 = 0, c16 = 0, c00 = 0;
    c00 += a00 * b00;
    c16 += c00 >>> 16;
    c00 &= 0xFFFF;
    c16 += a16 * b00;
    c32 += c16 >>> 16;
    c16 &= 0xFFFF;
    c16 += a00 * b16;
    c32 += c16 >>> 16;
    c16 &= 0xFFFF;
    c32 += a32 * b00;
    c48 += c32 >>> 16;
    c32 &= 0xFFFF;
    c32 += a16 * b16;
    c48 += c32 >>> 16;
    c32 &= 0xFFFF;
    c32 += a00 * b32;
    c48 += c32 >>> 16;
    c32 &= 0xFFFF;
    c48 += a48 * b00 + a32 * b16 + a16 * b32 + a00 * b48;
    c48 &= 0xFFFF;
    return goog.math.Long.fromBits((c16 << 16) | c00, (c48 << 16) | c32);
  };


  /**
   * Returns this Long divided by the given one.
   * @param {goog.math.Long} other Long by which to divide.
   * @return {!goog.math.Long} This Long divided by the given one.
   */
  goog.math.Long.prototype.div = function(other) {
    if (other.isZero()) {
      throw Error('division by zero');
    } else if (this.isZero()) {
      return goog.math.Long.ZERO;
    }

    if (this.equals(goog.math.Long.MIN_VALUE)) {
      if (other.equals(goog.math.Long.ONE) ||
          other.equals(goog.math.Long.NEG_ONE)) {
        return goog.math.Long.MIN_VALUE;  // recall that -MIN_VALUE == MIN_VALUE
      } else if (other.equals(goog.math.Long.MIN_VALUE)) {
        return goog.math.Long.ONE;
      } else {
        // At this point, we have |other| >= 2, so |this/other| < |MIN_VALUE|.
        var halfThis = this.shiftRight(1);
        var approx = halfThis.div(other).shiftLeft(1);
        if (approx.equals(goog.math.Long.ZERO)) {
          return other.isNegative() ? goog.math.Long.ONE : goog.math.Long.NEG_ONE;
        } else {
          var rem = this.subtract(other.multiply(approx));
          var result = approx.add(rem.div(other));
          return result;
        }
      }
    } else if (other.equals(goog.math.Long.MIN_VALUE)) {
      return goog.math.Long.ZERO;
    }

    if (this.isNegative()) {
      if (other.isNegative()) {
        return this.negate().div(other.negate());
      } else {
        return this.negate().div(other).negate();
      }
    } else if (other.isNegative()) {
      return this.div(other.negate()).negate();
    }

    // Repeat the following until the remainder is less than other:  find a
    // floating-point that approximates remainder / other *from below*, add this
    // into the result, and subtract it from the remainder.  It is critical that
    // the approximate value is less than or equal to the real value so that the
    // remainder never becomes negative.
    var res = goog.math.Long.ZERO;
    var rem = this;
    while (rem.greaterThanOrEqual(other)) {
      // Approximate the result of division. This may be a little greater or
      // smaller than the actual value.
      var approx = Math.max(1, Math.floor(rem.toNumber() / other.toNumber()));

      // We will tweak the approximate result by changing it in the 48-th digit or
      // the smallest non-fractional digit, whichever is larger.
      var log2 = Math.ceil(Math.log(approx) / Math.LN2);
      var delta = (log2 <= 48) ? 1 : Math.pow(2, log2 - 48);

      // Decrease the approximation until it is smaller than the remainder.  Note
      // that if it is too large, the product overflows and is negative.
      var approxRes = goog.math.Long.fromNumber(approx);
      var approxRem = approxRes.multiply(other);
      while (approxRem.isNegative() || approxRem.greaterThan(rem)) {
        approx -= delta;
        approxRes = goog.math.Long.fromNumber(approx);
        approxRem = approxRes.multiply(other);
      }

      // We know the answer can't be zero... and actually, zero would cause
      // infinite recursion since we would make no progress.
      if (approxRes.isZero()) {
        approxRes = goog.math.Long.ONE;
      }

      res = res.add(approxRes);
      rem = rem.subtract(approxRem);
    }
    return res;
  };


  /**
   * Returns this Long modulo the given one.
   * @param {goog.math.Long} other Long by which to mod.
   * @return {!goog.math.Long} This Long modulo the given one.
   */
  goog.math.Long.prototype.modulo = function(other) {
    return this.subtract(this.div(other).multiply(other));
  };


  /** @return {!goog.math.Long} The bitwise-NOT of this value. */
  goog.math.Long.prototype.not = function() {
    return goog.math.Long.fromBits(~this.low_, ~this.high_);
  };


  /**
   * Returns the bitwise-AND of this Long and the given one.
   * @param {goog.math.Long} other The Long with which to AND.
   * @return {!goog.math.Long} The bitwise-AND of this and the other.
   */
  goog.math.Long.prototype.and = function(other) {
    return goog.math.Long.fromBits(this.low_ & other.low_,
                                   this.high_ & other.high_);
  };


  /**
   * Returns the bitwise-OR of this Long and the given one.
   * @param {goog.math.Long} other The Long with which to OR.
   * @return {!goog.math.Long} The bitwise-OR of this and the other.
   */
  goog.math.Long.prototype.or = function(other) {
    return goog.math.Long.fromBits(this.low_ | other.low_,
                                   this.high_ | other.high_);
  };


  /**
   * Returns the bitwise-XOR of this Long and the given one.
   * @param {goog.math.Long} other The Long with which to XOR.
   * @return {!goog.math.Long} The bitwise-XOR of this and the other.
   */
  goog.math.Long.prototype.xor = function(other) {
    return goog.math.Long.fromBits(this.low_ ^ other.low_,
                                   this.high_ ^ other.high_);
  };


  /**
   * Returns this Long with bits shifted to the left by the given amount.
   * @param {number} numBits The number of bits by which to shift.
   * @return {!goog.math.Long} This shifted to the left by the given amount.
   */
  goog.math.Long.prototype.shiftLeft = function(numBits) {
    numBits &= 63;
    if (numBits == 0) {
      return this;
    } else {
      var low = this.low_;
      if (numBits < 32) {
        var high = this.high_;
        return goog.math.Long.fromBits(
            low << numBits,
            (high << numBits) | (low >>> (32 - numBits)));
      } else {
        return goog.math.Long.fromBits(0, low << (numBits - 32));
      }
    }
  };


  /**
   * Returns this Long with bits shifted to the right by the given amount.
   * @param {number} numBits The number of bits by which to shift.
   * @return {!goog.math.Long} This shifted to the right by the given amount.
   */
  goog.math.Long.prototype.shiftRight = function(numBits) {
    numBits &= 63;
    if (numBits == 0) {
      return this;
    } else {
      var high = this.high_;
      if (numBits < 32) {
        var low = this.low_;
        return goog.math.Long.fromBits(
            (low >>> numBits) | (high << (32 - numBits)),
            high >> numBits);
      } else {
        return goog.math.Long.fromBits(
            high >> (numBits - 32),
            high >= 0 ? 0 : -1);
      }
    }
  };


  /**
   * Returns this Long with bits shifted to the right by the given amount, with
   * the new top bits matching the current sign bit.
   * @param {number} numBits The number of bits by which to shift.
   * @return {!goog.math.Long} This shifted to the right by the given amount, with
   *     zeros placed into the new leading bits.
   */
  goog.math.Long.prototype.shiftRightUnsigned = function(numBits) {
    numBits &= 63;
    if (numBits == 0) {
      return this;
    } else {
      var high = this.high_;
      if (numBits < 32) {
        var low = this.low_;
        return goog.math.Long.fromBits(
            (low >>> numBits) | (high << (32 - numBits)),
            high >>> numBits);
      } else if (numBits == 32) {
        return goog.math.Long.fromBits(high, 0);
      } else {
        return goog.math.Long.fromBits(high >>> (numBits - 32), 0);
      }
    }
  };

  //======= begin jsbn =======

  var navigator = { appName: 'Modern Browser' }; // polyfill a little

  // Copyright (c) 2005  Tom Wu
  // All Rights Reserved.
  // http://www-cs-students.stanford.edu/~tjw/jsbn/

  /*
   * Copyright (c) 2003-2005  Tom Wu
   * All Rights Reserved.
   *
   * Permission is hereby granted, free of charge, to any person obtaining
   * a copy of this software and associated documentation files (the
   * "Software"), to deal in the Software without restriction, including
   * without limitation the rights to use, copy, modify, merge, publish,
   * distribute, sublicense, and/or sell copies of the Software, and to
   * permit persons to whom the Software is furnished to do so, subject to
   * the following conditions:
   *
   * The above copyright notice and this permission notice shall be
   * included in all copies or substantial portions of the Software.
   *
   * THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND, 
   * EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY 
   * WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  
   *
   * IN NO EVENT SHALL TOM WU BE LIABLE FOR ANY SPECIAL, INCIDENTAL,
   * INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND, OR ANY DAMAGES WHATSOEVER
   * RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER OR NOT ADVISED OF
   * THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF LIABILITY, ARISING OUT
   * OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
   *
   * In addition, the following condition applies:
   *
   * All redistributions must retain an intact copy of this copyright notice
   * and disclaimer.
   */

  // Basic JavaScript BN library - subset useful for RSA encryption.

  // Bits per digit
  var dbits;

  // JavaScript engine analysis
  var canary = 0xdeadbeefcafe;
  var j_lm = ((canary&0xffffff)==0xefcafe);

  // (public) Constructor
  function BigInteger(a,b,c) {
    if(a != null)
      if("number" == typeof a) this.fromNumber(a,b,c);
      else if(b == null && "string" != typeof a) this.fromString(a,256);
      else this.fromString(a,b);
  }

  // return new, unset BigInteger
  function nbi() { return new BigInteger(null); }

  // am: Compute w_j += (x*this_i), propagate carries,
  // c is initial carry, returns final carry.
  // c < 3*dvalue, x < 2*dvalue, this_i < dvalue
  // We need to select the fastest one that works in this environment.

  // am1: use a single mult and divide to get the high bits,
  // max digit bits should be 26 because
  // max internal value = 2*dvalue^2-2*dvalue (< 2^53)
  function am1(i,x,w,j,c,n) {
    while(--n >= 0) {
      var v = x*this[i++]+w[j]+c;
      c = Math.floor(v/0x4000000);
      w[j++] = v&0x3ffffff;
    }
    return c;
  }
  // am2 avoids a big mult-and-extract completely.
  // Max digit bits should be <= 30 because we do bitwise ops
  // on values up to 2*hdvalue^2-hdvalue-1 (< 2^31)
  function am2(i,x,w,j,c,n) {
    var xl = x&0x7fff, xh = x>>15;
    while(--n >= 0) {
      var l = this[i]&0x7fff;
      var h = this[i++]>>15;
      var m = xh*l+h*xl;
      l = xl*l+((m&0x7fff)<<15)+w[j]+(c&0x3fffffff);
      c = (l>>>30)+(m>>>15)+xh*h+(c>>>30);
      w[j++] = l&0x3fffffff;
    }
    return c;
  }
  // Alternately, set max digit bits to 28 since some
  // browsers slow down when dealing with 32-bit numbers.
  function am3(i,x,w,j,c,n) {
    var xl = x&0x3fff, xh = x>>14;
    while(--n >= 0) {
      var l = this[i]&0x3fff;
      var h = this[i++]>>14;
      var m = xh*l+h*xl;
      l = xl*l+((m&0x3fff)<<14)+w[j]+c;
      c = (l>>28)+(m>>14)+xh*h;
      w[j++] = l&0xfffffff;
    }
    return c;
  }
  if(j_lm && (navigator.appName == "Microsoft Internet Explorer")) {
    BigInteger.prototype.am = am2;
    dbits = 30;
  }
  else if(j_lm && (navigator.appName != "Netscape")) {
    BigInteger.prototype.am = am1;
    dbits = 26;
  }
  else { // Mozilla/Netscape seems to prefer am3
    BigInteger.prototype.am = am3;
    dbits = 28;
  }

  BigInteger.prototype.DB = dbits;
  BigInteger.prototype.DM = ((1<<dbits)-1);
  BigInteger.prototype.DV = (1<<dbits);

  var BI_FP = 52;
  BigInteger.prototype.FV = Math.pow(2,BI_FP);
  BigInteger.prototype.F1 = BI_FP-dbits;
  BigInteger.prototype.F2 = 2*dbits-BI_FP;

  // Digit conversions
  var BI_RM = "0123456789abcdefghijklmnopqrstuvwxyz";
  var BI_RC = new Array();
  var rr,vv;
  rr = "0".charCodeAt(0);
  for(vv = 0; vv <= 9; ++vv) BI_RC[rr++] = vv;
  rr = "a".charCodeAt(0);
  for(vv = 10; vv < 36; ++vv) BI_RC[rr++] = vv;
  rr = "A".charCodeAt(0);
  for(vv = 10; vv < 36; ++vv) BI_RC[rr++] = vv;

  function int2char(n) { return BI_RM.charAt(n); }
  function intAt(s,i) {
    var c = BI_RC[s.charCodeAt(i)];
    return (c==null)?-1:c;
  }

  // (protected) copy this to r
  function bnpCopyTo(r) {
    for(var i = this.t-1; i >= 0; --i) r[i] = this[i];
    r.t = this.t;
    r.s = this.s;
  }

  // (protected) set from integer value x, -DV <= x < DV
  function bnpFromInt(x) {
    this.t = 1;
    this.s = (x<0)?-1:0;
    if(x > 0) this[0] = x;
    else if(x < -1) this[0] = x+DV;
    else this.t = 0;
  }

  // return bigint initialized to value
  function nbv(i) { var r = nbi(); r.fromInt(i); return r; }

  // (protected) set from string and radix
  function bnpFromString(s,b) {
    var k;
    if(b == 16) k = 4;
    else if(b == 8) k = 3;
    else if(b == 256) k = 8; // byte array
    else if(b == 2) k = 1;
    else if(b == 32) k = 5;
    else if(b == 4) k = 2;
    else { this.fromRadix(s,b); return; }
    this.t = 0;
    this.s = 0;
    var i = s.length, mi = false, sh = 0;
    while(--i >= 0) {
      var x = (k==8)?s[i]&0xff:intAt(s,i);
      if(x < 0) {
        if(s.charAt(i) == "-") mi = true;
        continue;
      }
      mi = false;
      if(sh == 0)
        this[this.t++] = x;
      else if(sh+k > this.DB) {
        this[this.t-1] |= (x&((1<<(this.DB-sh))-1))<<sh;
        this[this.t++] = (x>>(this.DB-sh));
      }
      else
        this[this.t-1] |= x<<sh;
      sh += k;
      if(sh >= this.DB) sh -= this.DB;
    }
    if(k == 8 && (s[0]&0x80) != 0) {
      this.s = -1;
      if(sh > 0) this[this.t-1] |= ((1<<(this.DB-sh))-1)<<sh;
    }
    this.clamp();
    if(mi) BigInteger.ZERO.subTo(this,this);
  }

  // (protected) clamp off excess high words
  function bnpClamp() {
    var c = this.s&this.DM;
    while(this.t > 0 && this[this.t-1] == c) --this.t;
  }

  // (public) return string representation in given radix
  function bnToString(b) {
    if(this.s < 0) return "-"+this.negate().toString(b);
    var k;
    if(b == 16) k = 4;
    else if(b == 8) k = 3;
    else if(b == 2) k = 1;
    else if(b == 32) k = 5;
    else if(b == 4) k = 2;
    else return this.toRadix(b);
    var km = (1<<k)-1, d, m = false, r = "", i = this.t;
    var p = this.DB-(i*this.DB)%k;
    if(i-- > 0) {
      if(p < this.DB && (d = this[i]>>p) > 0) { m = true; r = int2char(d); }
      while(i >= 0) {
        if(p < k) {
          d = (this[i]&((1<<p)-1))<<(k-p);
          d |= this[--i]>>(p+=this.DB-k);
        }
        else {
          d = (this[i]>>(p-=k))&km;
          if(p <= 0) { p += this.DB; --i; }
        }
        if(d > 0) m = true;
        if(m) r += int2char(d);
      }
    }
    return m?r:"0";
  }

  // (public) -this
  function bnNegate() { var r = nbi(); BigInteger.ZERO.subTo(this,r); return r; }

  // (public) |this|
  function bnAbs() { return (this.s<0)?this.negate():this; }

  // (public) return + if this > a, - if this < a, 0 if equal
  function bnCompareTo(a) {
    var r = this.s-a.s;
    if(r != 0) return r;
    var i = this.t;
    r = i-a.t;
    if(r != 0) return (this.s<0)?-r:r;
    while(--i >= 0) if((r=this[i]-a[i]) != 0) return r;
    return 0;
  }

  // returns bit length of the integer x
  function nbits(x) {
    var r = 1, t;
    if((t=x>>>16) != 0) { x = t; r += 16; }
    if((t=x>>8) != 0) { x = t; r += 8; }
    if((t=x>>4) != 0) { x = t; r += 4; }
    if((t=x>>2) != 0) { x = t; r += 2; }
    if((t=x>>1) != 0) { x = t; r += 1; }
    return r;
  }

  // (public) return the number of bits in "this"
  function bnBitLength() {
    if(this.t <= 0) return 0;
    return this.DB*(this.t-1)+nbits(this[this.t-1]^(this.s&this.DM));
  }

  // (protected) r = this << n*DB
  function bnpDLShiftTo(n,r) {
    var i;
    for(i = this.t-1; i >= 0; --i) r[i+n] = this[i];
    for(i = n-1; i >= 0; --i) r[i] = 0;
    r.t = this.t+n;
    r.s = this.s;
  }

  // (protected) r = this >> n*DB
  function bnpDRShiftTo(n,r) {
    for(var i = n; i < this.t; ++i) r[i-n] = this[i];
    r.t = Math.max(this.t-n,0);
    r.s = this.s;
  }

  // (protected) r = this << n
  function bnpLShiftTo(n,r) {
    var bs = n%this.DB;
    var cbs = this.DB-bs;
    var bm = (1<<cbs)-1;
    var ds = Math.floor(n/this.DB), c = (this.s<<bs)&this.DM, i;
    for(i = this.t-1; i >= 0; --i) {
      r[i+ds+1] = (this[i]>>cbs)|c;
      c = (this[i]&bm)<<bs;
    }
    for(i = ds-1; i >= 0; --i) r[i] = 0;
    r[ds] = c;
    r.t = this.t+ds+1;
    r.s = this.s;
    r.clamp();
  }

  // (protected) r = this >> n
  function bnpRShiftTo(n,r) {
    r.s = this.s;
    var ds = Math.floor(n/this.DB);
    if(ds >= this.t) { r.t = 0; return; }
    var bs = n%this.DB;
    var cbs = this.DB-bs;
    var bm = (1<<bs)-1;
    r[0] = this[ds]>>bs;
    for(var i = ds+1; i < this.t; ++i) {
      r[i-ds-1] |= (this[i]&bm)<<cbs;
      r[i-ds] = this[i]>>bs;
    }
    if(bs > 0) r[this.t-ds-1] |= (this.s&bm)<<cbs;
    r.t = this.t-ds;
    r.clamp();
  }

  // (protected) r = this - a
  function bnpSubTo(a,r) {
    var i = 0, c = 0, m = Math.min(a.t,this.t);
    while(i < m) {
      c += this[i]-a[i];
      r[i++] = c&this.DM;
      c >>= this.DB;
    }
    if(a.t < this.t) {
      c -= a.s;
      while(i < this.t) {
        c += this[i];
        r[i++] = c&this.DM;
        c >>= this.DB;
      }
      c += this.s;
    }
    else {
      c += this.s;
      while(i < a.t) {
        c -= a[i];
        r[i++] = c&this.DM;
        c >>= this.DB;
      }
      c -= a.s;
    }
    r.s = (c<0)?-1:0;
    if(c < -1) r[i++] = this.DV+c;
    else if(c > 0) r[i++] = c;
    r.t = i;
    r.clamp();
  }

  // (protected) r = this * a, r != this,a (HAC 14.12)
  // "this" should be the larger one if appropriate.
  function bnpMultiplyTo(a,r) {
    var x = this.abs(), y = a.abs();
    var i = x.t;
    r.t = i+y.t;
    while(--i >= 0) r[i] = 0;
    for(i = 0; i < y.t; ++i) r[i+x.t] = x.am(0,y[i],r,i,0,x.t);
    r.s = 0;
    r.clamp();
    if(this.s != a.s) BigInteger.ZERO.subTo(r,r);
  }

  // (protected) r = this^2, r != this (HAC 14.16)
  function bnpSquareTo(r) {
    var x = this.abs();
    var i = r.t = 2*x.t;
    while(--i >= 0) r[i] = 0;
    for(i = 0; i < x.t-1; ++i) {
      var c = x.am(i,x[i],r,2*i,0,1);
      if((r[i+x.t]+=x.am(i+1,2*x[i],r,2*i+1,c,x.t-i-1)) >= x.DV) {
        r[i+x.t] -= x.DV;
        r[i+x.t+1] = 1;
      }
    }
    if(r.t > 0) r[r.t-1] += x.am(i,x[i],r,2*i,0,1);
    r.s = 0;
    r.clamp();
  }

  // (protected) divide this by m, quotient and remainder to q, r (HAC 14.20)
  // r != q, this != m.  q or r may be null.
  function bnpDivRemTo(m,q,r) {
    var pm = m.abs();
    if(pm.t <= 0) return;
    var pt = this.abs();
    if(pt.t < pm.t) {
      if(q != null) q.fromInt(0);
      if(r != null) this.copyTo(r);
      return;
    }
    if(r == null) r = nbi();
    var y = nbi(), ts = this.s, ms = m.s;
    var nsh = this.DB-nbits(pm[pm.t-1]);	// normalize modulus
    if(nsh > 0) { pm.lShiftTo(nsh,y); pt.lShiftTo(nsh,r); }
    else { pm.copyTo(y); pt.copyTo(r); }
    var ys = y.t;
    var y0 = y[ys-1];
    if(y0 == 0) return;
    var yt = y0*(1<<this.F1)+((ys>1)?y[ys-2]>>this.F2:0);
    var d1 = this.FV/yt, d2 = (1<<this.F1)/yt, e = 1<<this.F2;
    var i = r.t, j = i-ys, t = (q==null)?nbi():q;
    y.dlShiftTo(j,t);
    if(r.compareTo(t) >= 0) {
      r[r.t++] = 1;
      r.subTo(t,r);
    }
    BigInteger.ONE.dlShiftTo(ys,t);
    t.subTo(y,y);	// "negative" y so we can replace sub with am later
    while(y.t < ys) y[y.t++] = 0;
    while(--j >= 0) {
      // Estimate quotient digit
      var qd = (r[--i]==y0)?this.DM:Math.floor(r[i]*d1+(r[i-1]+e)*d2);
      if((r[i]+=y.am(0,qd,r,j,0,ys)) < qd) {	// Try it out
        y.dlShiftTo(j,t);
        r.subTo(t,r);
        while(r[i] < --qd) r.subTo(t,r);
      }
    }
    if(q != null) {
      r.drShiftTo(ys,q);
      if(ts != ms) BigInteger.ZERO.subTo(q,q);
    }
    r.t = ys;
    r.clamp();
    if(nsh > 0) r.rShiftTo(nsh,r);	// Denormalize remainder
    if(ts < 0) BigInteger.ZERO.subTo(r,r);
  }

  // (public) this mod a
  function bnMod(a) {
    var r = nbi();
    this.abs().divRemTo(a,null,r);
    if(this.s < 0 && r.compareTo(BigInteger.ZERO) > 0) a.subTo(r,r);
    return r;
  }

  // Modular reduction using "classic" algorithm
  function Classic(m) { this.m = m; }
  function cConvert(x) {
    if(x.s < 0 || x.compareTo(this.m) >= 0) return x.mod(this.m);
    else return x;
  }
  function cRevert(x) { return x; }
  function cReduce(x) { x.divRemTo(this.m,null,x); }
  function cMulTo(x,y,r) { x.multiplyTo(y,r); this.reduce(r); }
  function cSqrTo(x,r) { x.squareTo(r); this.reduce(r); }

  Classic.prototype.convert = cConvert;
  Classic.prototype.revert = cRevert;
  Classic.prototype.reduce = cReduce;
  Classic.prototype.mulTo = cMulTo;
  Classic.prototype.sqrTo = cSqrTo;

  // (protected) return "-1/this % 2^DB"; useful for Mont. reduction
  // justification:
  //         xy == 1 (mod m)
  //         xy =  1+km
  //   xy(2-xy) = (1+km)(1-km)
  // x[y(2-xy)] = 1-k^2m^2
  // x[y(2-xy)] == 1 (mod m^2)
  // if y is 1/x mod m, then y(2-xy) is 1/x mod m^2
  // should reduce x and y(2-xy) by m^2 at each step to keep size bounded.
  // JS multiply "overflows" differently from C/C++, so care is needed here.
  function bnpInvDigit() {
    if(this.t < 1) return 0;
    var x = this[0];
    if((x&1) == 0) return 0;
    var y = x&3;		// y == 1/x mod 2^2
    y = (y*(2-(x&0xf)*y))&0xf;	// y == 1/x mod 2^4
    y = (y*(2-(x&0xff)*y))&0xff;	// y == 1/x mod 2^8
    y = (y*(2-(((x&0xffff)*y)&0xffff)))&0xffff;	// y == 1/x mod 2^16
    // last step - calculate inverse mod DV directly;
    // assumes 16 < DB <= 32 and assumes ability to handle 48-bit ints
    y = (y*(2-x*y%this.DV))%this.DV;		// y == 1/x mod 2^dbits
    // we really want the negative inverse, and -DV < y < DV
    return (y>0)?this.DV-y:-y;
  }

  // Montgomery reduction
  function Montgomery(m) {
    this.m = m;
    this.mp = m.invDigit();
    this.mpl = this.mp&0x7fff;
    this.mph = this.mp>>15;
    this.um = (1<<(m.DB-15))-1;
    this.mt2 = 2*m.t;
  }

  // xR mod m
  function montConvert(x) {
    var r = nbi();
    x.abs().dlShiftTo(this.m.t,r);
    r.divRemTo(this.m,null,r);
    if(x.s < 0 && r.compareTo(BigInteger.ZERO) > 0) this.m.subTo(r,r);
    return r;
  }

  // x/R mod m
  function montRevert(x) {
    var r = nbi();
    x.copyTo(r);
    this.reduce(r);
    return r;
  }

  // x = x/R mod m (HAC 14.32)
  function montReduce(x) {
    while(x.t <= this.mt2)	// pad x so am has enough room later
      x[x.t++] = 0;
    for(var i = 0; i < this.m.t; ++i) {
      // faster way of calculating u0 = x[i]*mp mod DV
      var j = x[i]&0x7fff;
      var u0 = (j*this.mpl+(((j*this.mph+(x[i]>>15)*this.mpl)&this.um)<<15))&x.DM;
      // use am to combine the multiply-shift-add into one call
      j = i+this.m.t;
      x[j] += this.m.am(0,u0,x,i,0,this.m.t);
      // propagate carry
      while(x[j] >= x.DV) { x[j] -= x.DV; x[++j]++; }
    }
    x.clamp();
    x.drShiftTo(this.m.t,x);
    if(x.compareTo(this.m) >= 0) x.subTo(this.m,x);
  }

  // r = "x^2/R mod m"; x != r
  function montSqrTo(x,r) { x.squareTo(r); this.reduce(r); }

  // r = "xy/R mod m"; x,y != r
  function montMulTo(x,y,r) { x.multiplyTo(y,r); this.reduce(r); }

  Montgomery.prototype.convert = montConvert;
  Montgomery.prototype.revert = montRevert;
  Montgomery.prototype.reduce = montReduce;
  Montgomery.prototype.mulTo = montMulTo;
  Montgomery.prototype.sqrTo = montSqrTo;

  // (protected) true iff this is even
  function bnpIsEven() { return ((this.t>0)?(this[0]&1):this.s) == 0; }

  // (protected) this^e, e < 2^32, doing sqr and mul with "r" (HAC 14.79)
  function bnpExp(e,z) {
    if(e > 0xffffffff || e < 1) return BigInteger.ONE;
    var r = nbi(), r2 = nbi(), g = z.convert(this), i = nbits(e)-1;
    g.copyTo(r);
    while(--i >= 0) {
      z.sqrTo(r,r2);
      if((e&(1<<i)) > 0) z.mulTo(r2,g,r);
      else { var t = r; r = r2; r2 = t; }
    }
    return z.revert(r);
  }

  // (public) this^e % m, 0 <= e < 2^32
  function bnModPowInt(e,m) {
    var z;
    if(e < 256 || m.isEven()) z = new Classic(m); else z = new Montgomery(m);
    return this.exp(e,z);
  }

  // protected
  BigInteger.prototype.copyTo = bnpCopyTo;
  BigInteger.prototype.fromInt = bnpFromInt;
  BigInteger.prototype.fromString = bnpFromString;
  BigInteger.prototype.clamp = bnpClamp;
  BigInteger.prototype.dlShiftTo = bnpDLShiftTo;
  BigInteger.prototype.drShiftTo = bnpDRShiftTo;
  BigInteger.prototype.lShiftTo = bnpLShiftTo;
  BigInteger.prototype.rShiftTo = bnpRShiftTo;
  BigInteger.prototype.subTo = bnpSubTo;
  BigInteger.prototype.multiplyTo = bnpMultiplyTo;
  BigInteger.prototype.squareTo = bnpSquareTo;
  BigInteger.prototype.divRemTo = bnpDivRemTo;
  BigInteger.prototype.invDigit = bnpInvDigit;
  BigInteger.prototype.isEven = bnpIsEven;
  BigInteger.prototype.exp = bnpExp;

  // public
  BigInteger.prototype.toString = bnToString;
  BigInteger.prototype.negate = bnNegate;
  BigInteger.prototype.abs = bnAbs;
  BigInteger.prototype.compareTo = bnCompareTo;
  BigInteger.prototype.bitLength = bnBitLength;
  BigInteger.prototype.mod = bnMod;
  BigInteger.prototype.modPowInt = bnModPowInt;

  // "constants"
  BigInteger.ZERO = nbv(0);
  BigInteger.ONE = nbv(1);

  // jsbn2 stuff

  // (protected) convert from radix string
  function bnpFromRadix(s,b) {
    this.fromInt(0);
    if(b == null) b = 10;
    var cs = this.chunkSize(b);
    var d = Math.pow(b,cs), mi = false, j = 0, w = 0;
    for(var i = 0; i < s.length; ++i) {
      var x = intAt(s,i);
      if(x < 0) {
        if(s.charAt(i) == "-" && this.signum() == 0) mi = true;
        continue;
      }
      w = b*w+x;
      if(++j >= cs) {
        this.dMultiply(d);
        this.dAddOffset(w,0);
        j = 0;
        w = 0;
      }
    }
    if(j > 0) {
      this.dMultiply(Math.pow(b,j));
      this.dAddOffset(w,0);
    }
    if(mi) BigInteger.ZERO.subTo(this,this);
  }

  // (protected) return x s.t. r^x < DV
  function bnpChunkSize(r) { return Math.floor(Math.LN2*this.DB/Math.log(r)); }

  // (public) 0 if this == 0, 1 if this > 0
  function bnSigNum() {
    if(this.s < 0) return -1;
    else if(this.t <= 0 || (this.t == 1 && this[0] <= 0)) return 0;
    else return 1;
  }

  // (protected) this *= n, this >= 0, 1 < n < DV
  function bnpDMultiply(n) {
    this[this.t] = this.am(0,n-1,this,0,0,this.t);
    ++this.t;
    this.clamp();
  }

  // (protected) this += n << w words, this >= 0
  function bnpDAddOffset(n,w) {
    if(n == 0) return;
    while(this.t <= w) this[this.t++] = 0;
    this[w] += n;
    while(this[w] >= this.DV) {
      this[w] -= this.DV;
      if(++w >= this.t) this[this.t++] = 0;
      ++this[w];
    }
  }

  // (protected) convert to radix string
  function bnpToRadix(b) {
    if(b == null) b = 10;
    if(this.signum() == 0 || b < 2 || b > 36) return "0";
    var cs = this.chunkSize(b);
    var a = Math.pow(b,cs);
    var d = nbv(a), y = nbi(), z = nbi(), r = "";
    this.divRemTo(d,y,z);
    while(y.signum() > 0) {
      r = (a+z.intValue()).toString(b).substr(1) + r;
      y.divRemTo(d,y,z);
    }
    return z.intValue().toString(b) + r;
  }

  // (public) return value as integer
  function bnIntValue() {
    if(this.s < 0) {
      if(this.t == 1) return this[0]-this.DV;
      else if(this.t == 0) return -1;
    }
    else if(this.t == 1) return this[0];
    else if(this.t == 0) return 0;
    // assumes 16 < DB < 32
    return ((this[1]&((1<<(32-this.DB))-1))<<this.DB)|this[0];
  }

  // (protected) r = this + a
  function bnpAddTo(a,r) {
    var i = 0, c = 0, m = Math.min(a.t,this.t);
    while(i < m) {
      c += this[i]+a[i];
      r[i++] = c&this.DM;
      c >>= this.DB;
    }
    if(a.t < this.t) {
      c += a.s;
      while(i < this.t) {
        c += this[i];
        r[i++] = c&this.DM;
        c >>= this.DB;
      }
      c += this.s;
    }
    else {
      c += this.s;
      while(i < a.t) {
        c += a[i];
        r[i++] = c&this.DM;
        c >>= this.DB;
      }
      c += a.s;
    }
    r.s = (c<0)?-1:0;
    if(c > 0) r[i++] = c;
    else if(c < -1) r[i++] = this.DV+c;
    r.t = i;
    r.clamp();
  }

  BigInteger.prototype.fromRadix = bnpFromRadix;
  BigInteger.prototype.chunkSize = bnpChunkSize;
  BigInteger.prototype.signum = bnSigNum;
  BigInteger.prototype.dMultiply = bnpDMultiply;
  BigInteger.prototype.dAddOffset = bnpDAddOffset;
  BigInteger.prototype.toRadix = bnpToRadix;
  BigInteger.prototype.intValue = bnIntValue;
  BigInteger.prototype.addTo = bnpAddTo;

  //======= end jsbn =======

  // Emscripten wrapper
  var Wrapper = {
    abs: function(l, h) {
      var x = new goog.math.Long(l, h);
      var ret;
      if (x.isNegative()) {
        ret = x.negate();
      } else {
        ret = x;
      }
      HEAP32[tempDoublePtr>>2] = ret.low_;
      HEAP32[tempDoublePtr+4>>2] = ret.high_;
    },
    ensureTemps: function() {
      if (Wrapper.ensuredTemps) return;
      Wrapper.ensuredTemps = true;
      Wrapper.two32 = new BigInteger();
      Wrapper.two32.fromString('4294967296', 10);
      Wrapper.two64 = new BigInteger();
      Wrapper.two64.fromString('18446744073709551616', 10);
      Wrapper.temp1 = new BigInteger();
      Wrapper.temp2 = new BigInteger();
    },
    lh2bignum: function(l, h) {
      var a = new BigInteger();
      a.fromString(h.toString(), 10);
      var b = new BigInteger();
      a.multiplyTo(Wrapper.two32, b);
      var c = new BigInteger();
      c.fromString(l.toString(), 10);
      var d = new BigInteger();
      c.addTo(b, d);
      return d;
    },
    stringify: function(l, h, unsigned) {
      var ret = new goog.math.Long(l, h).toString();
      if (unsigned && ret[0] == '-') {
        // unsign slowly using jsbn bignums
        Wrapper.ensureTemps();
        var bignum = new BigInteger();
        bignum.fromString(ret, 10);
        ret = new BigInteger();
        Wrapper.two64.addTo(bignum, ret);
        ret = ret.toString(10);
      }
      return ret;
    },
    fromString: function(str, base, min, max, unsigned) {
      Wrapper.ensureTemps();
      var bignum = new BigInteger();
      bignum.fromString(str, base);
      var bigmin = new BigInteger();
      bigmin.fromString(min, 10);
      var bigmax = new BigInteger();
      bigmax.fromString(max, 10);
      if (unsigned && bignum.compareTo(BigInteger.ZERO) < 0) {
        var temp = new BigInteger();
        bignum.addTo(Wrapper.two64, temp);
        bignum = temp;
      }
      var error = false;
      if (bignum.compareTo(bigmin) < 0) {
        bignum = bigmin;
        error = true;
      } else if (bignum.compareTo(bigmax) > 0) {
        bignum = bigmax;
        error = true;
      }
      var ret = goog.math.Long.fromString(bignum.toString()); // min-max checks should have clamped this to a range goog.math.Long can handle well
      HEAP32[tempDoublePtr>>2] = ret.low_;
      HEAP32[tempDoublePtr+4>>2] = ret.high_;
      if (error) throw 'range error';
    }
  };
  return Wrapper;
})();

//======= end closure i64 code =======



// === Auto-generated postamble setup entry stuff ===

if (memoryInitializer) {
  if (ENVIRONMENT_IS_NODE || ENVIRONMENT_IS_SHELL) {
    var data = Module['readBinary'](memoryInitializer);
    HEAPU8.set(data, STATIC_BASE);
  } else {
    addRunDependency('memory initializer');
    Browser.asyncLoad(memoryInitializer, function(data) {
      HEAPU8.set(data, STATIC_BASE);
      removeRunDependency('memory initializer');
    }, function(data) {
      throw 'could not load memory initializer ' + memoryInitializer;
    });
  }
}

function ExitStatus(status) {
  this.name = "ExitStatus";
  this.message = "Program terminated with exit(" + status + ")";
  this.status = status;
};
ExitStatus.prototype = new Error();
ExitStatus.prototype.constructor = ExitStatus;

var initialStackTop;
var preloadStartTime = null;
var calledMain = false;

dependenciesFulfilled = function runCaller() {
  // If run has never been called, and we should call run (INVOKE_RUN is true, and Module.noInitialRun is not false)
  if (!Module['calledRun'] && shouldRunNow) run();
  if (!Module['calledRun']) dependenciesFulfilled = runCaller; // try this again later, after new deps are fulfilled
}

Module['callMain'] = Module.callMain = function callMain(args) {
  assert(runDependencies == 0, 'cannot call main when async dependencies remain! (listen on __ATMAIN__)');
  assert(__ATPRERUN__.length == 0, 'cannot call main when preRun functions remain to be called');

  args = args || [];

  ensureInitRuntime();

  var argc = args.length+1;
  function pad() {
    for (var i = 0; i < 4-1; i++) {
      argv.push(0);
    }
  }
  var argv = [allocate(intArrayFromString("/bin/this.program"), 'i8', ALLOC_NORMAL) ];
  pad();
  for (var i = 0; i < argc-1; i = i + 1) {
    argv.push(allocate(intArrayFromString(args[i]), 'i8', ALLOC_NORMAL));
    pad();
  }
  argv.push(0);
  argv = allocate(argv, 'i32', ALLOC_NORMAL);

  initialStackTop = STACKTOP;

  try {

    var ret = Module['_main'](argc, argv, 0);


    // if we're not running an evented main loop, it's time to exit
    if (!Module['noExitRuntime']) {
      exit(ret);
    }
  }
  catch(e) {
    if (e instanceof ExitStatus) {
      // exit() throws this once it's done to make sure execution
      // has been stopped completely
      return;
    } else if (e == 'SimulateInfiniteLoop') {
      // running an evented main loop, don't immediately exit
      Module['noExitRuntime'] = true;
      return;
    } else {
      if (e && typeof e === 'object' && e.stack) Module.printErr('exception thrown: ' + [e, e.stack]);
      throw e;
    }
  } finally {
    calledMain = true;
  }
}




function run(args) {
  args = args || Module['arguments'];

  if (preloadStartTime === null) preloadStartTime = Date.now();

  if (runDependencies > 0) {
    Module.printErr('run() called, but dependencies remain, so not running');
    return;
  }

  preRun();

  if (runDependencies > 0) return; // a preRun added a dependency, run will be called later
  if (Module['calledRun']) return; // run may have just been called through dependencies being fulfilled just in this very frame

  function doRun() {
    if (Module['calledRun']) return; // run may have just been called while the async setStatus time below was happening
    Module['calledRun'] = true;

    ensureInitRuntime();

    preMain();

    if (ENVIRONMENT_IS_WEB && preloadStartTime !== null) {
      Module.printErr('pre-main prep time: ' + (Date.now() - preloadStartTime) + ' ms');
    }

    if (Module['_main'] && shouldRunNow) {
      Module['callMain'](args);
    }

    postRun();
  }

  if (Module['setStatus']) {
    Module['setStatus']('Running...');
    setTimeout(function() {
      setTimeout(function() {
        Module['setStatus']('');
      }, 1);
      if (!ABORT) doRun();
    }, 1);
  } else {
    doRun();
  }
}
Module['run'] = Module.run = run;

function exit(status) {
  ABORT = true;
  EXITSTATUS = status;
  STACKTOP = initialStackTop;

  // exit the runtime
  exitRuntime();

  // TODO We should handle this differently based on environment.
  // In the browser, the best we can do is throw an exception
  // to halt execution, but in node we could process.exit and
  // I'd imagine SM shell would have something equivalent.
  // This would let us set a proper exit status (which
  // would be great for checking test exit statuses).
  // https://github.com/kripken/emscripten/issues/1371

  // throw an exception to halt the current execution
  throw new ExitStatus(status);
}
Module['exit'] = Module.exit = exit;

function abort(text) {
  if (text) {
    Module.print(text);
    Module.printErr(text);
  }

  ABORT = true;
  EXITSTATUS = 1;

  var extra = '\nIf this abort() is unexpected, build with -s ASSERTIONS=1 which can give more information.';

  throw 'abort() at ' + stackTrace() + extra;
}
Module['abort'] = Module.abort = abort;

// {{PRE_RUN_ADDITIONS}}

if (Module['preInit']) {
  if (typeof Module['preInit'] == 'function') Module['preInit'] = [Module['preInit']];
  while (Module['preInit'].length > 0) {
    Module['preInit'].pop()();
  }
}

// shouldRunNow refers to calling main(), not run().
var shouldRunNow = true;
if (Module['noInitialRun']) {
  shouldRunNow = false;
}


run();

// {{POST_RUN_ADDITIONS}}






// {{MODULE_ADDITIONS}}






