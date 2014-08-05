      subroutine derivs(x,x0,dv,Ndim,Ntraj,pe)

      implicit real*8(a-h,o-z)

      integer*4,intent(IN) :: Ntraj,Ndim
      real*8,intent(IN) :: x(Ndim,Ntraj),x0(Ndim)
      real*8,intent(OUT) :: dv(Ndim,Ntraj),pe(Ntraj)
      
      real*8 :: hess(Ndim,Ndim)

      call hessian(hess,Ndim)

      pi=4.d0*atan(1.d0)
      pe = 0d0
      dv = 0d0

      dv = matmul(hess,x)

      do i=1,Ntraj
        do j=1,Ndim
          do k=1,Ndim
!              dv(j,i)=hess(j,k)*x(k,i)+dvdx(j,i)
!        xi(i)=-exp(2d0*pi*x(2,i)/l)
!    &  B*xi(i))/((1-xi(i))**3*l)
            pe(i)=x(k,i)*hess(k,j)*x(j,i)/2d0+pe(i)
!     &  B*xi(i)/(1d0-xi(i))**2
          enddo
        enddo
      enddo

      return
      end subroutine
!---------------------------------------
! hessian matrix
!--------------------------------------
	subroutine hessian(hess,Ndim)
	
      implicit real*8(a-h,o-z)
	
      real*8 :: hess(Ndim,Ndim)
!	hess(1,1) = 1d0
!	hess(1,2) = 0.4d0
!	hess(2,2) = 1d0
!	hess(1,3) = 0.2d0
!	hess(2,3) = 0.4d0
!	hess(3,3) = 1d0
	hess = 0d0

	do j=1,Ndim
	  hess(j,j) = 1.d0
	enddo

	return
	end subroutine
