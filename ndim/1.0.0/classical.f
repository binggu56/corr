        program main1
        implicit none
	real*8 :: t,ax,ay,k1,k2,m1,m2
     	real*8 :: kinetic_energy,potential_energy,total_energy,dv,dvdx,dvdy,dt,dt2,s
	integer :: k,kmax,i
	integer, parameter :: Ntraj=2
	
	common/pot_par/m1,m2
 	m1=10d0
 	m2=1d0
! the initial conditions
	real,allocatable :: x(:),y(:),px(:),py(:) ! N trajectories of x,y
	do i=1,Ntraj
	x(i) = i/10d0
	y(i) = 0d0
	px(i) = i/10d0
	py(i) = 0d0
	end do
	allocate(x(Ntraj),y(Ntraj),px(Ntraj),py(Ntraj))


        open(5,file='IN')
!	read(5,*) Ntraj
        read(5,*) kmax,dt
        read(5,*) k1,k2 ! constants for potential energy
!	read(5,*) h12
!	read(5,*) m1,m2 ! mass of different dimension
        close(5)
        write(*,*) 'Classical Trajectories'

! time step
        dt2=dt/2d0
        t=0d0
        
        
! x,px,y,py define the center of GWP        
!       write(*,*) 'INITIAL GWP WIDTH and CENTER'
!	write(*,*) x,px
!	write(*,*) y,py
!	write(*,*) m1,m2     



!	call derivs(x,y,potential_energy,dvdx,dvdy)
	open(10,file='output')
	
	potential_energy=0d0
	kinetic_energy=0d0
        do k=1,kmax

! increase t by dt
	t=t+dt	

! half-step increments of momenta
	do i=1,Ntraj
	call derivs(x(i),y(i),potential_energy,dvdx,dvdy)
        px(i)=px(i)-dvdx*dt2
        py(i)=py(i)-dvdy*dt2

! full step increment of positions
        x(i)=x(i)+px(i)*dt/m1
        y(i)=y(i)+py(i)*dt/m2
        
!        call derivs(x,y,v,dvdx,dvdy)

! half-step increments of momenta
	px(i)=px(i)-dvdx*dt2
	py(i)=py(i)-dvdy*dt2
	

! update potential, kinetic, and total energy	
        
	potential_energy=(k1*x(i)*x(i)+k2*y(i)*y(i))/2d0+potential_energy
	kinetic_energy=px(i)*px(i)/(2d0*m1)+py(i)*py(i)/(2d0*m2)+kinetic_energy
	end do


	total_energy=potential_energy+kinetic_energy
	

! print out the results
	write(10,*) t,total_energy,kinetic_energy,potential_energy
	
	
	end do
		
	end program main

        

! determine (vector) a=dv/dt, given(vector) x   
     
        subroutine derivs(x,y,potential_energy,dvdx,dvdy)
	implicit none
	real*8 :: x,y,potential_energy,dvdx,dvdy,m1,m2
	common/pot_par/m1,m2
        potential_energy=(k1*x*x+k2*y*y)/2d0
        dvdx=k1*x
        dvdy=k2*y
        return

        end
       
