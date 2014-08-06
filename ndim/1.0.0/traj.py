##!/usr/bin/python

import numpy as np
import pylab as plt
import matplotlib as mpl
import seaborn as sns

sns.set_context("poster")

mpl.rcParams['lines.linewidth'] = 2

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

data = np.genfromtxt(fname='traj.dat') 

ncol = data.shape[1]

for x in range(1,20):
    plt.plot(data[:,0],data[:,x],'k-',linewidth=1)
    

#plt.figure(1) 
#plt.plot(x,y1,'-')
#plt.plot(x,y2,'g-')
#pl.ylim(0,1)
plt.xlabel('Time [a.u.]')
plt.ylabel('Positions')
plt.savefig('traj.pdf')
plt.show() 

