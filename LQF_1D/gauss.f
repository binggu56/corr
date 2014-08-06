C     
C     two Guassians multiplication <gi|gj> 

      program main
      use cdat
      implicit real*8(a-h,o-z)
      parameter (NP = 1000)
      parameter (Nb = 2)
      complex*16 :: gi,gj,c
      real*8,    dimension(Nb) :: q0(Nb),p0(Nb)
      complex*16,dimension(Nb,Nb) :: mat
      common/grid/xmin,xmax,dx

C      im = (0d0,1d0)
      xmin = -50d0
      xmax = 50d0
      dx = (xmax-xmin)/(Np-1)
      al = 0.5d0

c      do i=1,Nb
c        p0(i) = 2d0
c      enddo
      q0(1) = -1d0
      q0(2) = 1d0
      p0(1) = 0d0
      p0(2) = 0d0
      
C      do m = 1,Nb
C        do n = 1,Nb
C          c = (0d0,0d0)
C          do i=1,Np
C            x = xmin + dx*(i-1)
C            gi = exp(-al*(x-q0(m))**2+im*p0(m)*(x-q0(m)))
C            gj = exp(-al*(x-q0(n))**2+im*p0(n)*(x-q0(n)))
C            c = c + conjg(gi)*gj
C          enddo
C          mat(m,n) = c
C        enddo
C      enddo
C     
      call multi(np,nb,al,q0,p0,mat)

      print *,"matrix = "
      do i=1,Nb
        print *,(mat(i,j),j=1,Nb)
      enddo

      end program

C     ********************************
      subroutine multi(np,nb,al,q,p,mat)
      use cdat, only : pi, im
      implicit real*8(a-h,o-z)
      integer*4, intent(in) :: nb,np
      real*8   , intent(in) :: al,q(nb),p(nb)
      complex*16,intent(out) :: mat(nb,nb)
      complex*16 :: c,gi,gj
      common/grid/xmin,xmax,dx

C      im = (0d0,1d0)

      do m=1,nb
        do n=1,nb
          c = (0d0,0d0)
          do i=1,Np
            x = xmin + dx*(i-1)
            gi = dsqrt(dsqrt(2d0*al/pi))*exp(-al*(x-q(m))**2+
     &           im*p(m)*(x-q(m)))
            gj = dsqrt(dsqrt(2d0*al/pi))*exp(-al*(x-q(n))**2+
     &           im*p(n)*(x-q(n)))
            c = c + conjg(gi)*gj*dx
          enddo
          mat(m,n) = c
        enddo
      enddo
      end subroutine

