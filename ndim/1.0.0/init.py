import sys

Ndim = 2
fi = open('ini','w')
fi.write("1000 " + " # number of trajectories" + '\n')
fi.write(str(Ndim) + " # number of degree of freedom"+'\n')
fi.write("50,0.01 " + " # time steps"+'\n')
fi.write("0d0,"*(Ndim-1) + '0d0'+ "# inital momentum"+'\n')
fi.write("1,"*(Ndim-1) + '1'

