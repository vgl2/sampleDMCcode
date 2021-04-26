        subroutine calc_r(n,coord,nmono,r)
        use omp_lib
        implicit real*8(a-h,o-z)
        parameter (nmax=70000)
        parameter (ndim=54)
        parameter (ndimw=9)
        parameter (nvibw=3)
        parameter (nvib=18)
        dimension coord(ndim,nmax),coord_mono(ndimw,nmono,nmax),
     1  r(nvibw,nmono,nmax)
        do i = 1,n
            ip = 0
            do j = 1,nmono
                do k = 1,ndimw
                    ip = ip + 1
                    coord_mono(k,j,i) = coord(ip,i)
                enddo
            enddo
        enddo
!$omp parallel
!$omp do
        do i = 1,n
            do j = 1,nmono
                r(1,j,i)=sqrt(((coord_mono(4,j,i)-coord_mono(1,j,i))**2)
     1          +((coord_mono(5,j,i)-coord_mono(2,j,i))**2)+
     1          ((coord_mono(6,j,i)-coord_mono(3,j,i))**2))
                r(2,j,i)=sqrt(((coord_mono(7,j,i)-coord_mono(1,j,i))**2)
     1          +((coord_mono(8,j,i)-coord_mono(2,j,i))**2)+
     1          ((coord_mono(9,j,i)-coord_mono(3,j,i))**2))
                r(3,j,i)=dacos((((coord_mono(4,j,i)-coord_mono(1,j,i))*
     1          (coord_mono(7,j,i)-coord_mono(1,j,i)))+
     1          ((coord_mono(5,j,i)-coord_mono(2,j,i))*
     1          (coord_mono(8,j,i)-coord_mono(2,j,i)))+
     1          ((coord_mono(6,j,i)-coord_mono(3,j,i))*
     1          (coord_mono(9,j,i)-coord_mono(3,j,i))))/(r(1,j,i)*
     1          r(2,j,i)))
            enddo
        enddo
!$omp end parallel
        return
        end subroutine

        subroutine calc_psi(n,r,sig_bend,nmono,psi,dr,first_mat,
     1  second_mat)
        use omp_lib
        implicit real*8(a-h,o-z)
        parameter (nmax=70000)
        parameter (nvibw=3)
        parameter (nvib=18)
        parameter (nparam=500)
        common/pcoef/alpha,c2,re,rea,c1
        common/wa1/roh(nparam),roh1(nparam),wave_oh(nparam),
     1  feed_oh(nparam)
        common/wa2/rod(nparam),rod1(nparam),wave_od(nparam),
     1  feed_od(nparam)
        dimension r(nvibw,nmono,nmax),wf(nvibw,nmono,nmax),psi(nmax),
     1  dr(nvibw,nmono,nmax),first(nvibw,nmono,nmax),
     1  second(nvibw,nmono,nmax),first_mat(nvib,nmax),
     1  second_mat(nvib,nmax),wf2(nvibw,nmono,nmax)
        pi = dacos(-1.d0)
        do i = 1,n
            psi(i) = 1.d0
        enddo
!$omp parallel
!$omp do
        do i = 1,n
            do j = 1,nmono
                do k = 1,nvibw
                    if (k.ne.3) then
                        call splint_derivatives(roh1(:),wave_oh(:),
     1                  feed_oh(:),nparam,r(k,j,i),wf(k,j,i),
     1                  first(k,j,i),second(k,j,i))
                        psi(i) = psi(i)*wf(k,j,i)
c                       psi(i) = psi(i)
                    else
                        wf(k,j,i)=(alpha/((sig_bend**2)*pi))**(0.25)*
     1                  exp(-0.5*alpha*(((r(k,j,i)-rea)**2)/
     1                  (sig_bend**2)))
                        psi(i) = psi(i)*wf(k,j,i)
                    endif
                enddo
            enddo
        enddo
!$omp end parallel
!$omp parallel
!$omp do
        do i = 1,n
            do j = 1,nmono
                do k = 1,nvibw
                    if (k.ne.3) then
                        first(k,j,i) = (first(k,j,i)*psi(i))/wf(k,j,i)
                        second(k,j,i) = (second(k,j,i)*psi(i))/wf(k,j,i)
                    endif
                enddo
            enddo
        enddo
!$omp end parallel
!$omp parallel
!$omp do
        do i = 1,n
            do j = 1,nmono
                do k = 1,nvibw
                    if (k.ne.3) then
                        dr(k,j,i) = (2*first(k,j,i))/psi(i)
                    else
                        dr(k,j,i)=(-2*alpha*(r(k,j,i)-rea)/
     1                  (sig_bend**2))
                        first(k,j,i) = dr(k,j,i)*psi(i)/2.d0
                        second(k,j,i) = (((alpha**2*((r(k,j,i)-rea)**2))
     1                  /(sig_bend**4))-(alpha/(sig_bend**2)))*psi(i)
                    endif
                enddo
            enddo
        enddo
!$omp end parallel
!$omp parallel
!$omp do
        do i = 1,n
            ip = 0
            do j = 1,nmono
                do k = 1,nvibw
                    ip = ip + 1
                    first_mat(ip,i) = first(k,j,i)
                    second_mat(ip,i) = second(k,j,i)
                enddo
            enddo
        enddo
!$omp end parallel
        return
        end subroutine

        subroutine calc_t(n,coord,psi,sig_bend,nmono,first_mat,
     1  second_mat,t)
        use omp_lib
        implicit real*8(a-h,o-z)
        parameter (nmax=70000)
        parameter (nvib=18)
        parameter (nparam=500)
        parameter (ndim=54)
        common/pcoef/alpha,c2,re,rea,c1
        common/kin/weight(ndim),d(ndim)
        common /wa1/roh(nparam),roh1(nparam),wave_oh(nparam),
     1  feed_oh(nparam)
        common /wa2/rod(nparam),rod1(nparam),wave_od(nparam),
     1  feed_od(nparam)
        common/der/drdx(nvib,ndim,nmax),dr2dx2(nvib,ndim,nmax)
        dimension coord(ndim,nmax),psi(nmax),t(nmax),second(nvib,nvib),
     1  tot_deriv(ndim,nmax),tot_deriv1(ndim,nmax),tot_deriv2(ndim,nmax)
     1  ,first_mat(nvib,nmax),second_mat(nvib,nmax)
        do i = 1,n
            t(i) = 0.d0
            do k = 1,ndim
                tot_deriv1(k,i) = 0.d0
                tot_deriv(k,i) = 0.d0
                tot_deriv2(k,i) = 0.d0
            enddo
        enddo
!$omp parallel
!$omp do
        do i = 1,n
            do k =1,ndim
                do l = 1,nvib
                    do j = 1,nvib
                        if (l.eq.j) then
                            tot_deriv1(k,i) = tot_deriv1(k,i)+
     1                      (second_mat(l,i)*drdx(l,k,i)*drdx(j,k,i))
                        else
                            tot_deriv1(k,i) = tot_deriv1(k,i)+
     1                      (((first_mat(l,i)*drdx(l,k,i)*first_mat(j,i)
     1                      *drdx(j,k,i)))/(psi(i)))
                        endif
                    enddo
                enddo
                do l = 1,nvib
                    tot_deriv2(k,i) = tot_deriv2(k,i)+(first_mat(l,i)*
     1              dr2dx2(l,k,i))
                enddo
                tot_deriv(k,i) = tot_deriv1(k,i) + tot_deriv2(k,i)
            enddo
        enddo
!$omp end parallel
        do i = 1,n
            do j =1,ndim
                t(i) = t(i) -(0.5*(tot_deriv(j,i)/weight(j)))
            enddo
        enddo
        do i = 1,n
            t(i) = t(i)/psi(i)
        enddo
        return
        end subroutine
            
C
c       calculate second derivatives        
        
        
        
        subroutine initial_roh_wave(sig)
        implicit real*8(a-h,o-z)
        parameter (nparam=500)
        common/wa1/roh(nparam),roh1(nparam),wave_oh(nparam),
     1  feed_oh(nparam)
c        common/pcoef/alpha,c2,re,rea,c1
        dimension harm(nparam),wave2(nparam)
        open(unit=200,file='roh_morse.dat',status='old')
        do i = 1,nparam
            read(200,*) roh1(i),harm(i),wave2(i),wave_oh(i)
        enddo
        close(200)
        maxwave = maxloc(wave_oh,1)
        roh = roh1
        re = roh1(maxwave)
        do i =1 ,nparam
           roh1(i) = roh1(i) - re
        enddo 
        do i = 1,nparam
            roh1(i) = roh1(i)*sig
        enddo
        do i = 1,nparam
            roh1(i) = roh1(i) + re
        enddo
        dx = roh1(2)-roh1(1)
        tot_wave = 0.d0
        do i = 1,nparam
            tot_wave = tot_wave + wave_oh(i)**2*dx
        enddo
        do i = 1,nparam
            wave_oh(i) = wave_oh(i)/sqrt(tot_wave)
        enddo
        pt1 = (wave_oh(2)-wave_oh(1))/(roh1(2)-roh1(1))
        ptn =(wave_oh(nparam)-wave_oh(nparam-1))/(roh1(nparam)-
     1  roh1(nparam-1))
        call spline(roh1,wave_oh,nparam,pt1,ptn,feed_oh)
        return
        end subroutine

        subroutine initial_rod_wave(sig)
        implicit real*8(a-h,o-z)
        parameter (nparam=500)
        common/wa2/rod(nparam),rod1(nparam),wave_od(nparam),
     1  feed_od(nparam)
c        common/pcoef/alpha,c2,re,rea,c1
        dimension harm(nparam),wave2(nparam)
        open(unit=300,file='rod_morse.dat',status='old')
        do i = 1,nparam
            read(300,*) rod1(i),harm(i),wave2(i),wave_od(i)
        enddo
        close(300)
        maxwave = maxloc(wave_od,1)
        re = rod1(maxwave)
        do i =1 ,nparam
           rod1(i) = rod1(i) - re
        enddo 
        do i = 1,nparam
            rod1(i) = rod1(i)*sig
        enddo
        do i = 1,nparam
            rod1(i) = rod1(i) + re
        enddo
        dx = rod1(2)-rod1(1)
        tot_wave = 0.d0
        do i = 1,nparam
            tot_wave = tot_wave + wave_od(i)**2*dx
        enddo
        do i = 1,nparam
            wave_od(i) = wave_od(i)/sqrt(tot_wave)
        enddo
        pt1 = (wave_od(2)-wave_od(1))/(rod1(2)-rod1(1))
        ptn =(wave_od(nparam)-wave_od(nparam-1))/(rod1(nparam)-
     1  rod1(nparam-1))
        call spline(rod1,wave_od,nparam,pt1,ptn,feed_od)
        return
        end subroutine

