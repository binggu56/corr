objects = mod.o qm.o fcl.o qpot.o  gasdev.o ran2.o mativs.o
MODULES = cdat.mod 
FC = gfortran
switch = -O3 -g 
EXEC = qm.x 
#LIB = /usr/lib64/atlas/liblapack.so.3.0
LIB = -llapack  
$(EXEC): $(objects) $(MODULES)
	$(FC) -o $(EXEC) $(LIB) $(switch) $(objects)

%.o: %.f90
	$(FC) -c $<
%.o: %.f
	$(FC) -c  $<
%.mod: %.f90 
	$(FC) -c $<
clean:
	rm $(objects) 
veryclean:
	rm *.o *.dat *.mod 
