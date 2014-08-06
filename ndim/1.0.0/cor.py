##!/usr/bin/python

import numpy as np
import pylab as plt
import matplotlib as mpl
import seaborn as sns

sns.set_context("poster")

mpl.rcParams['lines.linewidth'] = 2

data = np.genfromtxt(fname='cor.dat') 

ncols = data.shape[1]

#for x in range(1,ncols):
plt.plot(data[:,0],data[:,1],linewidth=2,label='$\Re(C_{xx})$')
plt.plot(data[:,0],data[:,2],linewidth=2,label='$\Im(C_{xx})$')
plt.plot(data[:,0],data[:,3],linewidth=2,label='$|C_{xx}|$')

    

#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
#pl.ylim(0,1)
plt.legend()
plt.xlabel('Time [a.u.]')
plt.ylabel('Positions')
plt.savefig('traj.pdf')
plt.show() 

