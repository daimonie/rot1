#from http://www.jesshamrick.com/2012/04/29/the-demise-of-for-loops/
from timeit import Timer
import numpy as np
import math

def timer(*funcs):
    # find the maximum function name length
    if len(funcs) > 1:
        maxlen = max(*[len(func) for func in funcs])
    elif len(funcs) == 1:
        maxlen = len(funcs[0])
    else:
        return

    # run each function 10000 times and print statistics
    times = []
    print "--"
    for func in funcs:
        timerfunc = Timer("%s()" % func, "from __main__ import %s" % func)
        runtime = timerfunc.repeat(repeat=10000, number=1)
        mtime = np.mean(runtime)
        stime = np.std(runtime)
        dfunc = func + (" " * (maxlen - len(func) + 1))
        print "%s: %.6f +/- %.6f seconds" % (dfunc, mtime, stime)
    
    
def numpy_linspace():
  total = np.linspace(0,100,1000);
def numpy_arange():
  total = np.arange(0,100,0.1);
  
  
def numpy_linspace2():
  total = np.linspace(0,100,100000);
def numpy_arange2():
  total = np.arange(0,100,0.001);
  
  
def numpy_linspace3():
  total = np.linspace(0,100,10000000);
def numpy_arange3():
  total = np.arange(0,100,0.00001);
  
timer("numpy_linspace","numpy_arange")
timer("numpy_linspace2","numpy_arange2")


def aranger():
	total = np.arange(0,100,0.25);
def ranger():
	total = range(0,100,400)
def aranger2():
	total = np.arange(0,100,0.025);
def ranger2():
	total = range(0,100,4000)
def aranger3():
	total = np.arange(0,100,0.0025);
timer("aranger","ranger");
timer("aranger2","ranger2");