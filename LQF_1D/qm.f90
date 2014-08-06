      program main
      use cdat

      implicit real*8(a-h,o-z)

     	real*8,allocatable :: ke(:),v(:),x(:),y(:),p(:),w(:)
      real*8,allocatable :: qfx(:),qfy(:),qp(:),fx(:),ddv(:),s(:)
      complex*16,allocatable :: h(:,:),mat(:,:)
     	real*8 :: f(3),ki

      complex*16,allocatable :: c(:),c0(:),hc(:),g(:),dg(:)
      real*8 :: gasdev

      complex*16 :: psi,cor

      common/wave/al,q0,p0

      open(100,file='energy.dat')
      open(101,file = 'x.dat')
      open(102,file='aver.dat')  
      open(103,file='wf.dat')
      open(104,file='wf0.dat')
      open(105,file = 'cor.dat')

      open(5,file='IN')
      read(5,*) Ntraj
      read(5,*) kmax,dt
      read(5,*) am
      read(5,*) al
      read(5,*) q0
      read(5,*) p0
      read(5,*) idum1
      read(5,*) nb
      close(5)

      write(*,1002) 
1002  format(/'1D QTM code with LQF, version 1.0'//)

      allocate(ke(Ntraj),v(Ntraj),x(Ntraj),p(Ntraj),s(Ntraj),w(Ntraj),qfx(Ntraj), &
               qp(Ntraj),fx(Ntraj),ddv(Ntraj))

      allocate(mat(Ntraj,Ntraj),c(Ntraj),c0(Ntraj),h(Ntraj,Ntraj),hc(Ntraj))     

      allocate(g(nb),dg(nb))

      dt2=dt/2d0
      t=0d0

      g(1) = (0d0,0d0)
      g(2) = (1d0,0d0)

!      Nb = Ntraj

!     grid size
!      xmin = -4d0
!      xmax = 4d0
!      Np = 1000
!      dx = (xmax-xmin)/(np-1)
      xav = 0d0
	do i=1,Ntraj
          x(i) = gasdev(idum1)
          x(i)=x(i)/sqrt(4d0*al)+q0
          xav = xav+w(i)*x(i)
	enddo 
      p = p0

      write(*,1003) xav
1003  format('position average =', f13.6)

!     print out the initial conditions        
      write(*,1001) al,Ntraj,kmax,dt,am,p0,q0
1001  format('Initial Conditions'// ,  &
            'ax    = ', f10.6/,        & 
            'Ntraj = ', i6/ ,          &
            'Kmax  = ', i6/ ,          &
            'dt    = ', f10.6/ ,       &
            'Mass  = ', f10.6/ ,       &
            'p0    = ', f10.6/ ,       &
            'q0    = ', f10.6/)
       
      w=1d0/Ntraj

!     begin the time step loop
      do 10 k=1,kmax
        
        call  prop(am,dt,nb,ntraj,w,x,p,g,dg)

        do j=1,nb
          g(j) = g(j)+dg(j)*dt
        enddo

        cor = (0d0,0d0)

        do i=1,ntraj
          cor = cor+conjg(g(1)+g(2)*x(i))*x(i)*w(i)
        enddo

        write(105,1000) t,cor,abs(cor)

        t = t+dt2

!       propogate trajectories with LQF

        call derivs(x,Ntraj,fx,v,ddv)
        call qpot(am,x,w,Ntraj,qp,qfx)

!       half-step increments of moment, full step increment of positions
        do 11 i=1,Ntraj
          p(i)=p(i)+(fx(i)+qfx(i))*dt2
          x(i)=x(i)+p(i)*dt/am
11      enddo
        
        t = t+dt2

        call derivs(x,Ntraj,fx,v,ddv)
        call qpot(am,x,w,Ntraj,qp,qfx)

        do 12 i=1,Ntraj
!         half-step increments of momenta
          p(i)=p(i)+(fx(i)+qfx(i))*dt2
          ke(i)=p(i)**2/(2d0*am)
12      enddo


!     calculate action function
      
!	 do i = 1,Ntraj
!	  s(i)=s(i)+(ke(i)-pe(i)-qp(i))*dt
!	 end do

!       print out a random trajectory	
        write(101,10000) t,(x(i),i=1,20)


!       calculate the expectation values
      
       call aver(Ntraj,ke,w,enk)
       call aver(ntraj,v, w,env)
       call aver(ntraj,qp,w,enq)

       en = enk+env+enq
       write(100,1000) t,enk,env,enq,en
        x2 = 0d0
        do i=1,Ntraj
          x2 = x2+x(i)**2*w(i)
        enddo
          
        write(102,10000) t,x2


!     end of time loop
10    end do

10000 format(2000(e14.7,1x))
1000  format(20(e14.7,1x))
      end program main

        
!     ----------------------------------------------
!     average
!     ---------------------------------------------
      subroutine aver(Ntraj,x,w,xav)
      implicit real*8(a-h,o-z)
      integer*4,intent(in)::Ntraj
      real*8,intent(in)::x(Ntraj),w(Ntraj)
      real*8,intent(out)::xav

      xav = 0d0
      do i=1,Ntraj
        xav = xav+w(i)*x(i)
      enddo

      return
      end subroutine


      subroutine prop(am,dt,nb,ntraj,w,x,p,c,dc)

      use cdat, only: im

! --- propagate exp(-iHt)*x*psi(x,0) 
! --- approximated by [c0(t)+c1(t)*x]*psi(x,t) 
! --- get coefficients c 

      implicit real*8(a-h,o-z)

      integer*4, intent(IN) :: nb,ntraj

      real*8,    intent(in) :: dt,am,w(ntraj),x(ntraj),p(ntraj)
      
      real*8 :: f(nb),df(nb),s(nb,nb),sinv(nb,nb)
      
      complex*16 :: cor,m1(nb,nb),m(nb,nb),m2(nb,nb)

      complex*16 :: c(nb),dc(nb)

      df = (/0d0,1d0/)
      m1 = (0d0,0d0)
      s  = (0d0,0d0)

      do i=1,nb
        do j=1,nb
          m2(j,i) = df(j)*df(i)
        enddo
      enddo

      do i=1,ntraj
        f = (/1d0,x(i)/)
        do j=1,nb
          do k=1,nb
            m1(k,j) = m1(k,j)+p(i)*w(i)*f(k)*df(j)
            s(k,j)  = s(k,j)+w(i)*f(k)*f(j)
          enddo
        enddo
      enddo

!      do i=1,nb
!        do j=1,nb
!          m(j,i) = 2d0*m1(j,i) + im*m2(j,i)
!          m(j,i) = -1d0/2d0/am*m(j,i)
!        enddo
!      enddo
      m = 2d0*m1+im*m2 
      m = -m/2d0/am 

      dc = matmul(m,c)
!      do i=1,nb
!        dc(i) = (0d0,0d0)
!        do j=1,nb
!          dc(i) = dc(i)+m(i,j)*c(j)
!        enddo
!      enddo
      
      call inverse(s,sinv,nb)

      dc = matmul(sinv,dc)

!      call zposv('U', nb, 1, s, nb, dc, nb, INFO)

!      if (info .ne. 0) then 
!         write(*,*) 'ZPOSV fails. INFO =', info
!         stop
!       endif
      
      return 
      end subroutine 



      subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
      implicit none 
      integer n
      double precision a(n,n), c(n,n)
      double precision L(n,n), U(n,n), b(n), d(n), x(n)
      double precision coeff
      integer i, j, k
      
      ! step 0: initialization for matrices L and U and b
      ! Fortran 90/95 aloows such operations on matrices
      L=0.0
      U=0.0
      b=0.0
      
      ! step 1: forward elimination
      do k=1, n-1
         do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
               a(i,j) = a(i,j)-coeff*a(k,j)
            end do
         end do
      end do
      
      ! Step 2: prepare L and U matrices 
      ! L matrix is a matrix of the elimination coefficient
      ! + the diagonal elements are 1.0
      do i=1,n
        L(i,i) = 1.0
      end do
      ! U matrix is the upper triangular part of A
      do j=1,n
        do i=1,j
          U(i,j) = a(i,j)
        end do
      end do
      
      ! Step 3: compute columns of the inverse matrix C
      do k=1,n
        b(k)=1.0
        d(1) = b(1)
      ! Step 3a: Solve Ld=b using the forward substitution
        do i=2,n
          d(i)=b(i)
          do j=1,i-1
            d(i) = d(i) - L(i,j)*d(j)
          end do
        end do
      ! Step 3b: Solve Ux=d using the back substitution
        x(n)=d(n)/U(n,n)
        do i = n-1,1,-1
          x(i) = d(i)
          do j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
          end do
          x(i) = x(i)/u(i,i)
        end do
      ! Step 3c: fill the solutions x(n) into column k of C
        do i=1,n
          c(i,k) = x(i)
        end do
        b(k)=0.0
      end do
      end subroutine inverse


      
