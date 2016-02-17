from numpy import *
from pylab import *

X=linspace(-pi,pi,256,endpoint=True)
C=cos(X)
S=sin(X)
figure(figsize=(8,6),dpi=80)
plot(X,C,'b-',linewidth=2.0,label='cos(x)')
plot(X,S,'r--',linewidth=2.0,label='sin(x)')
legend(loc='upper left')
xlim(-4.0,4.0)
# xticks(linspace(-4.0,4.0,9,endpoint=True))
xticks([-pi,-pi/2,0,pi/2,pi],
	[r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$',])
xlim(X.min()*1.1, X.max()*1.1)
ylim(-1.0,1.0)
# yticks(linspace(-1,1,5,endpoint=True))
yticks([-1,0,1],[r'$-1$',r'$0$',r'$+1$'])
ylim(C.min()*1.1,C.max()*1.1)

ax=gca()
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.spines['bottom'].set_position(('data',0))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))

t=2*pi/3
plot([t,t],[0,cos(t)],'b--',linewidth=1.5)
scatter([t,],[cos(t),],40,color='green')
annotate(r'$cos(\frac{2\pi}{3})=-\frac{1}{2}$',
	xy=(t,cos(t)),xycoords='data',
	xytext=(-90,-50),textcoords='offset points',fontsize=16,
	arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2"))

plot([t,t],[0,sin(t)],'r--',linewidth=1.5)
scatter([t,],[sin(t),],40,color='green')
annotate(r'$sin(\frac{2\pi}{3})=\frac{\sqrt{3}}{2}$',
	xy=(t,sin(t)),xycoords='data',
	xytext=(-10,+25),textcoords='offset points',fontsize=16,
	arrowprops=dict(arrowstyle="->",connectionstyle="arc3,rad=.2"))

for label in ax.get_xticklabels()+ax.get_yticklabels():
	label.set_fontsize(16)
	label.set_bbox(dict(facecolor='white',edgecolor='None',alpha=0.5))

show()
