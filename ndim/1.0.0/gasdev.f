      DOUBLE PRECISION FUNCTION GASDEV(IDUM)
C
C      USES ran1
C      Returns a normally distributed deviate with zero mean and unit
C      variance, using ran1(idum) as the source of uniform deviates.
C
C     .. Scalar Arguments ..
      INTEGER IDUM
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FAC,GSET,RSQ,V1,V2
      INTEGER ISET
C     ..
C     .. External Functions ..
      DOUBLE PRECISION RAN1
      EXTERNAL RAN1
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LOG,SQRT
C     ..
C     .. Save statement ..
      SAVE ISET,GSET
C     ..
C     .. Data statements ..
      DATA ISET/0/
C     ..
      IF (ISET.EQ.0) THEN
   10    V1 = 2.D0*RAN2(IDUM) - 1.D0
         V2 = 2.D0*RAN2(IDUM) - 1.D0
         RSQ = V1**2 + V2**2
         IF (RSQ.GE.1.D0 .OR. RSQ.EQ.0.D0) GO TO 10
         FAC = SQRT(-2.D0*LOG(RSQ)/RSQ)
         GSET = V1*FAC
         GASDEV = V2*FAC
         ISET = 1
      ELSE
         GASDEV = GSET
         ISET = 0
      END IF
      RETURN
      END

