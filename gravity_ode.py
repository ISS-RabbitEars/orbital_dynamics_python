import numpy as np
from scipy.integrate import odeint
from matplotlib import pyplot as plt
from matplotlib import animation


def a_alpha(rt, t, p):
	r,v,theta,omega=rt
	m,G=p
	return [v,r*(omega**2)-G*(m/(r**2)),omega,-2*v*omega/r]
mc=1
m=1
r=1
v=0.5
theta=0
omega=1
G=1

theta=(theta*np.pi/180)

p=[mc,G]
rt=[r,v,theta,omega]

tf = 240
nfps = 60
nframes = tf * nfps
t = np.linspace(0, tf, nframes)

rth = odeint(a_alpha, rt, t, args = (p,))

r1=rth[:,0]
th1=rth[:,2]

x=r1*np.cos(th1)
y=r1*np.sin(th1)

xmin=min(x)-0.2
xmax=max(x)+0.2
ymin=min(y)-0.2
ymax=max(y)+0.2


vr=rth[:,1]
w=rth[:,3]

ke=0.5*m*((vr**2)+((r1*w)**2))
pe=-G*m*mc/r1
E=ke+pe
Emax=abs(max(E))
ke/=Emax
pe/=Emax
E/=Emax
ke+=1
pe+=1
E+=1

fig, a=plt.subplots()
fig.tight_layout()

def run(frame):
	plt.clf()
	plt.subplot(211)
	circle=plt.Circle((x[frame],y[frame]),radius=0.05,fc='r')
	plt.gca().add_patch(circle)
	circle=plt.Circle((0,0),radius=0.05,fc='y')
	plt.gca().add_patch(circle)
	plt.title("Orbital Dynamics")
	ax=plt.gca()
	ax.set_aspect(1)
	plt.xlim([xmin,xmax])
	plt.ylim([ymin,ymax])
	ax.xaxis.set_ticklabels([])
	ax.yaxis.set_ticklabels([])
	ax.xaxis.set_ticks_position('none')
	ax.yaxis.set_ticks_position('none')
	ax.set_facecolor('xkcd:black')
	plt.subplot(212)
	plt.plot(t[0:frame],ke[0:frame],'r',lw=0.5)
	plt.plot(t[0:frame],pe[0:frame],'b',lw=0.5)
	plt.plot(t[0:frame],E[0:frame],'g',lw=0.5)
	plt.xlim([0,tf])
	plt.title("Energy (Rescaled and Shifted)")
	ax=plt.gca()
	ax.legend(['T','V','E'],labelcolor='w',frameon=False)
	ax.set_facecolor('xkcd:black')


ani=animation.FuncAnimation(fig,run,frames=nframes)
writervideo = animation.FFMpegWriter(fps=nfps)
ani.save('gravity_ode_wgphs.mp4', writer=writervideo)

plt.show()


