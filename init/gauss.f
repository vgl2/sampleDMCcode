        SUBROUTINE GAUSS(AM,SD,V,nr,dseed)
c       gaussian distr random nums
C       SD is the input variance of the gaussian distribution
c       half width (I belive) ABM 2/17/94
C       AM is the input mean of the gaussian distribution
c       dseed - input big number of order of 1.d10;
c       defined as double precision
c       but should be integer
c       nr - input num of rand numbers desired
c       V - vector of output random numbers
c       PROGRAM ASSUMED TO BE SET TO DOUBLE PRECISION!!!
        implicit real*8 (a-h,o-z)
        parameter (npts = 12)
        dimension v(nr),tmp(npts*nr)
        half = float(npts)*.5
        call ggubs(dseed,npts*nr,tmp)
c        call r8vec_uniform_01(npts*nr,dseed,tmp)
        ip = 0
        v = 0.
c       am = 0.
c       sd = 1.
        do n=1,npts
           do m = 1,nr
              ip = ip + 1
              v(m) = v(m) + tmp(ip)
           enddo
        enddo
        do n = 1,nr
C          v(n) = (v(n)/float(npts)-.5)*sd+am
           V(n)=(v(n)-half)*SD+AM
        enddo
        RETURN
        END
