        subroutine ph2(k,Ndim,Ntraj,x,v)
        implicit real*8(a-h,o-z)
        integer*4,intent(in)::Ndim,Ntraj
        real*8,intent(in)::x(Ntraj,Ndim)
        real*8,intent(out)::v
        real*8 :: r(Ndim/3,Ndim/3)
        
        Np = Ndim/3        
        c6 = 12.14d0
        c8 = 215.2d0
        c10 = 4813.9d0
        c9 = 143.1d0
        al = 1.713
        be = 1.5671
        rm = 3.44d0/0.52918d0
        

!        do k=1,Ntraj
          do i=1,Np
            do j=i+1,Np
              r = dsqrt((x(3*j-2,k)-x(3*i-2,k))**2+
     &                 (x(3*j-1,k)-x(3*i-1,k))**2+
     &                 (x(3*j,k)-x(3*i,k))**2)
              if(r < 1.28*rm) then
                fc = exp(-(1.28d0*rm/r-1d0)**2)
              else
                fc = 1.d0
              endif
              v = exp(al-be*r-ga*r**2)-
     &                (c6/r**6+c8/r**8+c10/r**10)*fc+
     &                c9/r**9*fc
            enddo
          enddo
!        enddo      
        return
        end subroutine
        
        subroutine cforce(m,n,Ndim,Ntraj,x,fx,v)
        implicit real*8(a-h,o-z)
        integer*4,intent(in)::m,n,Ndim,Ntraj
        real*8,intent(out)::fx,v
        real*8 :: y(Ndim,Ntraj)
        call ph2(m,Ndim,Ntraj,x,v)
        dx = 1d-3
        y = x
        y(n,m) = y(n,m)+dx
        call ph2(m,Ndim,Ntraj,y,v1)
        fx = -(v1-v)/dx

        
        return
        end subroutine        

        
!        do m=1,Ntraj
!          do i=1,Np
!            r = dsqrt((x(3*j-2,k)-x(3*i-2,k))**2+
!     &                 (x(3*j-1,k)-x(3*i-1,k))**2+
!     &                 (x(3*j,k)-x(3*i,k))**2)
!              if(r < 1.28*rm) then
!                fc = exp(-(1.28d0*rm/r-1d0)**2)
!                dfc = 2.56d0*(1.28d0*rm/r-1d0)*rm*fc
!              else
!                fc = 1.d0
!                dfc = 0d0
!              endif
!              f1 = exp(al-be*r-ga*r**2)
!              df1 = (-be-2d0*ga*r)*f1
!              f2 = c6/r**6+c8/r**8+c10/r**10
!
!
!            
!
!        subroutine dist(x,r,Np,Ntraj,Ndim)
!        implicit real*8(a-h,o-z)
!
!        r = 0d0
!        do i=1,Np
!          do j=i+1,Np
!            r(i,j) = dsqrt((x(3*j-2)-x(3*i-2))**2+
!     &               (x(3*j-1)-x(3*i-1))**2+
!     &               (x(3*j-1)-x(3*i))**2)
!          enddo
!        enddo
