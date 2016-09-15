#!/usr/bin/env python
from time import time
NN = 100000
t1 = time()
def fact(n):
    result = 1
    for i in xrange(1,n+1):
        result *= i
    return result

fact(NN)
print "python function takes %.7fs"%(time()-t1)

tm = time()
from ctypes import cdll
import os
print "importing takes       %.7fs"%(time()-tm)
## = = = = = use c recursion = = = = =
# - - - - [c_recursion.c] - - - -
# int factorial(int n){
#     if (n < 2)
#         return 1;
#     else
#         return factorial(n-1);
# }
# - - - -   end .c file  - - - -
# compile:
# $ gcc -c -fPIC c_recursion.c
# $ gcc -shared c_recursion.c -o c_recursion.so
#
# in windows, 
# > cl -LD c_recursion.c -c_recursion.dll
# and replace .so with .dll
t2 = time()
c_recursion = cdll.LoadLibrary(os.getcwd() + '/c_recursion.so')
a1 = c_recursion.factorial(NN)
print "c recursion takes     %.7fs"%(time() - t2)

## = = = = = use c for-loop = = = = =
# - - - - [c_forloop.c] - - - -
# int factor(int n){
#     int result = 1;
#     for (int i = 1; i <= n; ++i){
#     	result *= i;
#     }
#     return result;
# }
# - - - -   end .c file  - - - -
# compile:
# $ gcc -c -fPIC c_forloop.c
# $ gcc -shared c_forloop.c -o c_forloop.so
t3 = time()
c_forloop = cdll.LoadLibrary(os.getcwd() + '/c_forloop.so')
a2 = c_forloop.factor(NN)
print "c for-loop takes      %.7fs"%(time()-t3)

## result
# if NN = 100,000, output is:
# python function takes 2.8248279s
# importing takes       0.0025470s
# c recursion takes     0.0017569s
# c for-loop takes      0.0005419s
#
# if NN = 10,000, output is:
# python function takes 0.0223629s
# importing takes       0.0027299s
# c recursion takes     0.0003290s
# c for-loop takes      0.0002100s
#
# so c function is much faster than python.
# Larger the number is, faster it runs
# and for-loop is faster than recursion in large number case
