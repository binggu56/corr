      program main
!**********************************************************************
!	N dimensional quantum trajectory code with linear quantum force
!	test with coupled oscilator for each DOF (3 dim)
!*********************************************************************
      use cdat 

      implicit real*8(a-h,o-z)
     
      integer*4 :: Ntraj,counter
      
      real*8,dimension(:,:),allocatable :: dvdx,x,p,qf
      
      integer*4,dimension(:),allocatable :: idum
      
      real*8,dimension(:),allocatable :: am,p0,w,alpha,x0,pe,ke,qp,cf
      
      real*8 :: gasdev

      complex*16, dimension(:), allocatable :: g,dg

      complex*16 :: cor 

! --- the initial conditions
      open(100,file='energy.dat')
      open(101,file='traj.dat')
      open(103, file='cor.dat',status='unknown',action='write')

      open(10,file='IN', status='old', action='read')
      
      read(10,*) Ntraj
      read(10,*) Ndim
      read(10,*) kmax,dt
      read(10,*) cf0
      read(10,*) am0 
      
      close(10)

! --- linear basis 
      nb = ndim+1

      allocate(dvdx(Ndim,Ntraj),pe(Ntraj),ke(Ntraj),   &
                x(Ndim,Ntraj),p(Ndim,Ntraj),           &
                qp(Ntraj),qf(Ndim,Ntraj),w(Ntraj))

      allocate(p0(Ndim),idum(Ndim),alpha(Ndim),am(Ndim),x0(Ndim),cf(Ndim)) 

      allocate(g(nb),dg(nb))

!	read(10,*) (p0(i),i=1,Ndim)
!	read(10,*) (idum(i),i=1,Ndim)
!	read(10,*) (alpha(i),i=1,Ndim)
!	read(10,*) (am(i),i=1,Ndim)
!	read(10,*) (x0(i),i=1,Ndim)


      dt2 = dt/2d0
      t = 0d0
      p0 = 0d0
      x0 = 0d0
      cf = cf0
      am = am0 
      
      sigma = 1d0

      alpha = 1d0/2d0/sigma**2

      call seed(idum,Ndim)

      write(6,1001) sigma,Ntraj,kmax,dt,am0,p0,q0
1001  format('Initial Conditions'/ ,                   &
            'intial gaussian width = ', f10.6/,        & 
            'Ntraj = ', i6/ ,          &
            'Kmax  = ', i6/ ,          &
            'dt    = ', f10.6/ ,       &
            'Mass  = ', f10.6/ ,       &
            'p0    = ', f10.6/ ,       &
            'q0    = ', f10.6/)

! initial grid points 
      do i=1,Ntraj
        do j=1,Ndim
          x(j,i)=gasdev(idum(j))
          x(j,i)=x(j,i)/sqrt(4d0*alpha(j))+x0(j)
        enddo
      enddo

! initial momentum 
	do i=1,Ntraj
	  do j=1,Ndim
          p(j,i)=p0(j)
        enddo
	enddo

      w=1d0/dble(Ntraj)
      do j=1,ndim 
        g(j) = (1d0,0d0)
      enddo 
      g(ndim+1) = (0d0,0d0)

!-----------------------------------------	

! ---  time propagation	
      time: do counter=1,kmax

! ----- increase t by dt
        
        t=t+dt

! ----- corrleation function         
        call corr(am,dt,nb,ndim,ntraj,w,x,p,g,cor)

! ----- propagate trajectories
        call leap_frog(dt,ndim,ntraj,cf,am,x,x0,p,w,env,enk,enq)

        write(103,10000) t,cor,abs(cor) 
      
        write(101,10000) t,(x(1,i),i=1,20)

!	  s = 0d0
        write(100,10000) t,enk,env,enq,(enk+env+enq)

      enddo time 

10000 format(100(e14.7,1x))
      
      stop 
      end program 

!---------------------------------------------
! average over trajectories y = sum(w(i)*x(i)) 
!---------------------------------------------
	subroutine aver(Ntraj,w,x,y)
	implicit real*8(a-h,o-z)
	real*8 :: x(Ntraj),w(Ntraj)

	y = 0d0

	do i=1,Ntraj
	  y = y+x(i)*w(i) 
	enddo
	
	return
	end subroutine
!----------------------------------------------
!     assign idum[Ndim]

	subroutine seed(idum,Ndim)
	implicit real*8(a-h,o-z)
	integer*4, intent(IN) :: Ndim
	integer*4, intent(OUT) :: idum(Ndim)
	do i=1,Ndim
	  idum(i) = 5 + i
	enddo

	return
	end subroutine
!----------------------------------------------

      subroutine prop_c(am,dt,nb,ndim,ntraj,w,x,p,c,dc)

      use cdat, only : im 

      implicit real*8 (a-h, o-z) 

      integer*4, intent(in) :: ndim,ntraj,nb
      
      real*8,    intent(in) :: dt,am(ndim),w(ntraj),x(ndim,ntraj), & 
                               p(ndim,ntraj)

      real*8 :: f(nb),df(ndim,nb),s(nb,nb),sinv(nb,nb)

      complex*16 :: cor,mat1(nb,nb),mat(nb,nb),mat2(nb,nb), c(nb),dc(nb)

      s = (0d0,0d0) 
      mat1 = (0d0,0d0) 
      mat = (0d0,0d0)
      am0 = am(1)

! --- get df(ndim,nb)
      df = 0d0 
      do k=1,ndim 
        df(k,k) = 1d0 
      enddo 

! --- matrix mat2 = <df|df> 
!      do i=1,nb 
!        do j=1,nb 
!          do k=1,ndim 
!            mat2(j,i) = df(k,j)*df(k,i)
!          enddo 
!        enddo 
!      enddo 
!      df_t = transpose(df)
      mat2 = matmul(transpose(df),df)

! --- matrix matrix mat = <f|f> 
      do i=1,ntraj 

        do j=1,ndim 
          f(j) = x(j,i)
        enddo 
        f(ndim+1) = 1d0 
        
        do k=1,nb
          do j=1,k  
            s(j,k) = s(j,k)+f(j)*f(k)*w(i)
          enddo 
        enddo 
      enddo 

      do k=1,nb 
        do j=1,k
          s(k,j) = s(j,k)
        enddo 
      enddo       

! --- matrix mat1 = <p*f|df>, simplified version for linear basis 
!     needs modification if using higher order basis 
      traj: do i=1,ntraj 
        do j=1,ndim 
          f(j) = x(j,i) 
        enddo 
        f(ndim+1) = 0d0 

        do j=1,nb 
          do k=1,nb 
            mat1(k,j) = mat1(k,j)+p(j,i)*f(k)*df(j,j)*w(i)
          enddo 
        enddo 
      enddo traj 

      mat = 2d0*mat1 + im*mat2
      mat = -mat/2d0/am0

      dc = (0d0,0d0) 
      dc = matmul(mat,c)

      call inverse(s,sinv,nb)

      dc = matmul(sinv,dc)
      
      return 
      end subroutine 

! --- subroutine to solve ODE problem S*dc/dt = h*c 

!      subroutine ode(n,s,h,c,dc)

!      implicit real*8(a-h, o-z) 

! --- subroutine 
      subroutine corr(am,dt,nb,ndim,ntraj,w,x,p,g,cor)
     
      implicit real*8 (a-h,o-z)

      real*8, intent(in) :: am(ndim),w(ntraj),x(ndim,ntraj),p(ndim,ntraj)

      integer*4 ,intent(in) :: nb,ntraj,ndim 

      real*8 :: f(nb)

      complex*16,intent(inout) :: g(nb)

      complex*16,intent(out) :: cor 
      
      complex*16 :: z0,dg(nb)

      cor = (0d0,0d0) 

      do i=1,ntraj 
        do j=1,ndim 
          f(j) = x(j,i)
        enddo 
        f(ndim+1) = 1d0 
        
        z0 = dot_product(f,g)

        do j=1,ndim
          do k=1,nb 
            cor = cor+ conjg(z0)*x(j,i)*w(i)
          enddo 
        enddo 

      enddo

! --- update c 
      call prop_c(am,dt,nb,ndim,ntraj,w,x,p,g,dg)

      do j=1,nb 
        g(j) = g(j)+dg(j)*dt 
      enddo 


      return 
      end subroutine 

! --- leap_frog subroutine to solve ode 
      subroutine leap_frog(dt,ndim,ntraj,cf,am,x,x0,p,w,env,enk,enq)

      use cdat

      implicit real*8(a-h, o-z)

      integer*4,intent(in)    :: ntraj,ndim
      real*8,   intent(in)    :: dt

      real*8, intent(in), dimension(ntraj)  :: w
      real*8, intent(in), dimension(ndim) :: am,cf,x0

      real*8, intent(inout), dimension(ndim,ntraj) :: x,p

      real*8, dimension(ndim,ntraj) :: qf,dv
!      real*8, dimension(ndim,)       :: dv

      real*8, dimension(ntraj)  :: pe,qp,ke

!      real*8 :: q(NATOM3)
     

      dt2 = dt/2d0
      pe = 0d0

      call derivs(x,x0,dv,Ndim,Ntraj,pe)
      call qpot(am,x,w,Ndim,Ntraj,qp,qf)

! half-step increments of momenta & full step increment of positions
      do i=1,Ntraj
        do j=1,Ndim
          p(j,i)=p(j,i)+(-dv(j,i)+qf(j,i)-cf(j)*p(j,i))*dt2
          x(j,i)=x(j,i)+p(j,i)*dt/am(j)
        enddo   
      enddo 
    
      call derivs(x,x0,dv,Ndim,Ntraj,pe)
      call qpot(am,x,w,Ndim,Ntraj,qp,qf)

! half-step increments of momenta
      do i=1,Ntraj
        do j=1,Ndim
          p(j,i)=p(j,i)+(-dv(j,i)+qf(j,i)-cf(j)*p(j,i))*dt2
        enddo
      enddo
!------------------------------------------------------------
! update potential, kinetic, and total energy
        ke = 0d0
        do i=1, Ntraj
          do j=1,Ndim
            ke(i)=p(j,i)**2/(2d0*am(j))+ke(i)
          enddo
        enddo

      call aver(Ntraj,w,pe,env)
      call aver(Ntraj,w,ke,enk)
      call aver(Ntraj,w,qp,enq)

      return
      end subroutine
