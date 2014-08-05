      module cdat

        implicit real*8(a-h,o-z)
      
!      integer*4, public, parameter  :: 
        real*8, public, parameter  :: pi = 4.d0*atan(1.d0)

        complex*16, public, parameter :: im=(0d0,1d0)

      save

      contains

        double precision function trace(n,A)

!        implicit real*8(a-h,o-z)

        real*8, dimension(n,n), intent(in) :: A

        trace = 0d0 

        do i=1,n
          trace = trace + A(i,i)
        enddo 

        end function trace  

      
      end module cdat
