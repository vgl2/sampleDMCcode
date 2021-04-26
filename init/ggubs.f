      SUBROUTINE GGUBS (DSEED,NR,R)
! get a vector of rand nums uniformly distributed 0-1
! dseed - input big number of order of 1.d10;
! defined as double precision
! but should be integer
! nr - input num of rand numbers desired
! r - vector of output random numbers
! PROGRAM ASSUMED TO BE SET TO DOUBLE PRECISION!!!
         implicit real*8 (a-h,o-z)
!                                  SPECIFICATIONS FOR ARGUMENTS
      INTEGER            NR
      REAL*8           R(NR)
      DOUBLE PRECISION   DSEED
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I
      REAL*8           S2P31,S2P31M,SEED
!                                  S2P31M = (2**31) - 1
!                                  S2P31 = (2**31)
      DATA               S2P31M/2147483647.E0/,S2P31/2147483648.E0/
!                                  FIRST EXECUTABLE STATEMENT
      SEED = DSEED
      DO 5 I=1,NR
         SEED = dMOD(16807.E0*SEED,S2P31M)
    5 R(I) = SEED / S2P31
      DSEED = SEED
      RETURN
      END 

