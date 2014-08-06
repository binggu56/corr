##!/usr/bin/python

import numpy as np
import pylab as pl

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
data0 = np.genfromtxt(fname='dat2/en.dat') 
data = np.genfromtxt(fname='dat3/en.dat')
dat = np.genfromtxt(fname='../1.0.1/dat1/en.dat')
dat1 = np.genfromtxt(fname='../1.0.2/dat4/en.dat')
#data1 = np.genfromtxt(fname='err.dat')
#data2 = np.genfromtxt(fname='../2d8/err.dat')
#dat1 = np.genfromtxt(fname='../2d7/en.dat')
#dat2 = np.genfromtxt(fname='../2d7/err.dat')

#data = np.loadtxt('traj.dat')
#for x in range(1,10):
pl.subplot(111)
pl.xlim(0,4)
#pl.title('$\gamma = 0.5$, $V(x,y)=x^2/2+x^4/4+y^2/2+y^4/4+\gamma xy$')
pl.ylabel('Energy [hartree]')
#pl.plot(data[:,0],data[:,2],'b-',linewidth=2,label='Potential')
#pl.plot(data[:,0],data[:,3],'g-',linewidth=2,label='Quantum Potential')
pl.plot(data[:,0],data[:,4],'k-',linewidth=2,label='One-step')
#pl.plot(data0[:,0],data0[:,2],'b--',linewidth=2,label='Potential')
#pl.plot(data0[:,0],data0[:,3],'g--',linewidth=2,label='Quantum Potential')
pl.plot(data0[:,0],data0[:,4],'g-.',linewidth=2,label='three steps')
pl.plot(dat[:,0],dat[:,4],'g-.',linewidth=2,label='9600, 12, 0.75')
pl.plot(dat1[:,0],dat1[:,4],'g-.',linewidth=2,label='9600, 12, 0.8')

#pl.plot(dat1[:,0],dat1[:,4],'r--',linewidth=2,label='two steps')
#pl.legend(bbox_to_anchor=(0.5, 0.38, 0.42, .302), loc=3,ncol=1, mode="expand", borderaxespad=0.)
pl.legend()

#pl.subplot(212)
#pl.xlim(0,4)
#pl.xlabel('time [a.u.]')
#pl.plot(data1[:,0],data1[:,1],'r-',linewidth=2,label='2s err($r_x$)')
##pl.plot(data1[:,0],data1[:,2],'g-',linewidth=2,label='2s err($r_y$)')
#
#pl.plot(data2[:,0],data2[:,1],'r--',linewidth=2,label='err($r_x$)')
##pl.plot(data2[:,0],data2[:,2],'g--',linewidth=2,label='err($r_y$)')
##pl.plot(dat2[:,0],dat2[:,1],'k--',linewidth=2,label='3s err($r_x$)')
pl.legend(loc=1)
pl.savefig('err.pdf')
pl.show() 

