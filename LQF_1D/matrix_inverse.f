
C  --------------------------------------------
C  Do matrix inversion for square matrix A(N,N)
C  Output will be written to A(N,N)
C  --------------------------------------------
      subroutine matrix_inverse(N,A)
      IMPLICIT real*8(a-h,o-z)
      integer*4,intent(in) :: N
      complex*16,intent(in) :: A(N,N)
      complex*16 :: WORK(N*N)
      INTEGER :: LWORK, LDA
      INTEGER :: IPIV(N)

C      INTEGER DeAllocateStatus
C      CHARACTER :: A_filename*20
C      external ZGETRF
C      external ZGETRI

      LDA = N
      LWORK = N*N
C      ALLOCATE (A(LDA,N))
C      ALLOCATE (Ainv(LDA,N))
C      ALLOCATE (WORK(LWORK))
C      ALLOCATE (IPIV(N))

C     DGETRF computes an LU factorization of a general M-by-N matrix A
C     using partial pivoting with row interchanges.

C     Store A in Ainv to prevent it from being overwritten by LAPACK
C     Ainv = A
C      print *,'Before inversion'
C      do i=1,N
C        print *,(A(i,j),j=1,N)
C      enddo

      CALL ZGETRF( N, N, A, LDA, IPIV, INFO )

C      IF(INFO.EQ.0)THEN
C         PRINT '(" LU decomposition successful ")'
C      ENDIF
      IF(INFO.LT.0)THEN
         PRINT '(" LU decomposition:  illegal value ")'
         STOP
      ENDIF
      IF(INFO.GT.0)THEN
         WRITE(*,35)INFO,INFO
35       FORMAT( 'LU decomposition: U(',I4,',',I4,') = 0 ')
      ENDIF

C  DGETRI computes the inverse of a matrix using the LU factorization
C  computed by DGETRF.

      CALL ZGETRI(N, A, N, IPIV, WORK, LWORK, INFO)
C      PRINT '(" ")'
      IF (info.NE.0) THEN
         stop 'Matrix inversion failed!'
      ENDIF


C      DEALLOCATE (A, STAT = DeAllocateStatus)
C      DEALLOCATE (Ainv, STAT = DeAllocateStatus)
C      DEALLOCATE (IPIV, STAT = DeAllocateStatus)
C      DEALLOCATE (WORK, STAT = DeAllocateStatus)

      return
      END
