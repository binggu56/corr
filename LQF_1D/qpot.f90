        subroutine qpot(m1,x,w,Ntraj,qp,qfx)
        integer*4,intent(IN) :: Ntraj
        real*8, intent(OUT) :: qfx(Ntraj),qp(Ntraj)
        integer*4 :: i,j
	integer INFO
        real*8,intent(IN) :: m1
        real*8,intent(IN) :: w(Ntraj),x(Ntraj)
        real*8 :: f(2),rx(Ntraj)
	real*8 ::  s(2,2),c(2,1)
! define c matrix 
        c(1,1)=-0.5d0
        c(2,1)=0d0

! quantum force
        s=0d0
        do i=1,Ntraj
          f=(/x(i),1d0/)
          do m=1,2
            do n=1,2
              s(m,n)=w(i)*f(m)*f(n)+s(m,n)
            end do
          end do
        end do

        !do m=1,3
        !print *, s(m,1),s(m,2),s(m,3)
        !enddo
        !STOP
        
! calculate matrix c(t)
        call DPOSV('U',2,1,s,2,c,2,INFO)
        if(INFO/=0) then
        print *, "info=",info
        print *, "matrix fails"
        stop
        end if

! the momentum operator r=cf
        do i=1,Ntraj
        rx(i)=c(1,1)*x(i)+c(2,1)

! calculate quantum potential
        qp(i)=-rx(i)**2/(2d0*m1)-c(1,1)/(2d0*m1)
        qfx(i)=rx(i)*c(1,1)/m1

        end do

        return
        end subroutine
