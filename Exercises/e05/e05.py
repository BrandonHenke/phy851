import matplotlib.pyplot as plt
import numpy as np
import math
import os
from pylab import *

NMAX=10000
WF=zeros(NMAX+1)
XMAX=40.0
DX=XMAX/NMAX
a=0.707
B=2.2
hbarc=197.327
mu=0.5*939.0
x=DX*arange(0,NMAX+1)

def GetQ2(x,V0):
	V=-V0/(1.0+exp((x-a)/a))
	q2= 2*mu*(V+B)/hbarc**2
	return q2,V

def GetIntercept(V0):
	q=sqrt(2.0*mu*B)/hbarc
	WF[NMAX]=exp(-q*NMAX*DX)
	WF[NMAX-1]=exp(-q*(NMAX-1)*DX)
	N=NMAX-1
	while N>0:
		qx2,V=GetQ2(N*DX,V0)
		WF[N-1]=2*WF[N]-WF[N+1]+DX**2*(qx2+V*2*mu/hbarc**2)*WF[N]
		N=N-1

	return WF[0]

def GetR2():
	R2bar=0.0
	norm=0.0
	N=0
	while N <= NMAX:
		r=N*DX
		WF2=WF[N]*WF[N]
		norm+=DX*WF2
		R2bar+=DX*WF2*r*r
		N=N+1
	R2bar=R2bar/norm
	return R2bar

V0 = float(input("Enter V0: "))
intercept=GetIntercept(V0)
print('Intercept=',intercept)
R2=GetR2()
print('R^2bar=',R2,', R=',np.sqrt(R2))

plt.figure(figsize=(6,5))
fig = plt.figure(1)
plt.plot(x,WF,linewidth=3,color='k',label='$\\beta a$=0.5')
plt.xlim(0.0,15)
plt.ylim(-1.0,1.0)
plt.xlabel('$x$ (fm)', fontsize=18, weight='normal')
plt.ylabel('$\psi$',fontsize=18)

# If you have a mac, this will open Preview
plt.savefig('exercise5_fig.pdf',format='pdf')
os.system('open -a Preview exercise5_fig.pdf')
# if you have something else, use this instead
#plt.show()
quit()
