from numpy import *
import matplotlib.pyplot as plt

def f(x,y):
	return (1-x/2+x**5+y**3)*exp(-x**2-y**2)

n=256
x=linspace(-3,3,n)
y=linspace(-3,3,n)
X,Y = meshgrid(x,y)

plt.figure()
plt.axes([0.025,0.025,0.95,0.95])

plt.contourf(X,Y,f(X,Y),8,alpha=0.75,cmap=plt.cm.hot)
C=plt.contour(X,Y,f(X,Y),8,colors='black',linewidth=0.5)
plt.clabel(C,inline=1,fontsize=10)

plt.xticks([]), plt.yticks([])
plt.savefig('./figures/contour_ex.png',dpi=60)
plt.show()