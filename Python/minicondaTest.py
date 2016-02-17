#!~/Program/miniconda/bin
from numpy import *
import matplotlib.pyplot as plt
x=linspace(0,2,101)
y=sin(x*2*pi)
plt.plot(x,y,'r--')
plt.show()
