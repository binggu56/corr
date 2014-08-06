#!/usr/bin/python

import numpy as np
import pylab as pl
import seaborn as sns

#with open("traj.dat") as f:
#    data = f.read()
#
#    data = data.split('\n')
#
#    x = [row.split(' ')[0] for row in data]
#    y = [row.split(' ')[1] for row in data]
#
#    fig = plt.figure()
#
#    ax1 = fig.add_subplot(111)
#
#    ax1.set_title("Plot title...")    
#    ax1.set_xlabel('your x label..')
#    ax1.set_ylabel('your y label...')
#
#    ax1.plot(x,y, c='r', label='the data')
#
#    leg = ax1.legend()
#fig = plt.figure()
data = np.genfromtxt(fname='energy.dat')

pl.plot(data[:,0],data[:,1],linewidth=2,label='Kinetic')
pl.plot(data[:,0],data[:,2],linewidth=2,label='Potential')
pl.plot(data[:,0],data[:,3],linewidth=2,label='Quantum Potential')
pl.plot(data[:,0],data[:,4],linewidth=2,label='Energy')
#pl.legend(bbox_to_anchor=(0.5, 0.38, 0.42, .302), loc=3,ncol=1, mode="expand", borderaxespad=0.)
pl.legend()
pl.ylabel('Energy [hartree]')
pl.savefig('energy.pdf')
pl.show() 

