        subroutine mc_sample(dseed,coord,ntot,ngood,freq,evec,it,psips)
        implicit real*8(a-h,o-z)
        parameter (natoms=9)
        parameter (nmax=10000)
        parameter (ndim=54)
        parameter (nvib=ndim-6)
        parameter (nvibw=3)
        parameter (nmono=6)
        parameter (nparam=500)
        common/wa1/roh(nparam),roh1(nparam),wave_oh(nparam),
     1  feed_oh(nparam)
        common/wa2/rod(nparam),rod1(nparam),wave_od(nparam),
     1  feed_od(nparam)
        common/pcoef/alpha,c2,re,rea,c1
        dimension coord(ndim),freq(nvib),evec(ndim,nvib),
     1  output(nmax*nvib),check1(nvib,nmax),psi(nmax),
     1  check2(ndim,nmax),evec_trans(nvib,ndim),psips(ndim,nmax),
     1  r(nvibw,nmono,nmax),psit(nmax),dr(nvibw,nmono,nmax),
     1  first(18,nmax),second(18,nmax)
        pi = dacos(-1.d0)
        sig = 2.d0
        bend = 1.d0
        o_mass = 15.99491461957d0*1822.88852962d0
        h_mass = 1.00782503223d0*1822.88852962d0
        d_mass = 2.01410177812d0*1822.88852962d0
        nslab = int(ngood/ntot)
        evec_trans = transpose(evec)
        psi_max = 1.d0
        do i = 1,nvib
            psi_max = psi_max*((freq(i)/(sig**2*pi))**(0.25))
        enddo
        it = 0
        do while (it.ne.nslab)
            call gauss(0.d0,1.d0,output,nmax*nvib,dseed)
            ip = 0
            do i =1,nmax
                do j = 1,nvib
                    ip = ip + 1
                    check1(j,i) = output(ip)/(sqrt(freq(j))/sig)
                enddo 
            enddo
            do i = 1,nmax
                check2(:,i) = matmul(check1(:,i),evec_trans(:,:))
            enddo
            do i = 1,nmax
                do j = 1,ndim
                    if (j.lt.4) then
                        check2(j,i) = check2(j,i)/sqrt(o_mass)
                    else if ((j.gt.9).and.(j.lt.13)) then
                        check2(j,i) = check2(j,i)/sqrt(o_mass)
                    else if ((j.gt.18).and.(j.lt.22)) then
                        check2(j,i) = check2(j,i)/sqrt(o_mass)
                    else if ((j.gt.27).and.(j.lt.31)) then
                        check2(j,i) = check2(j,i)/sqrt(o_mass)
                    else if ((j.gt.36).and.(j.lt.40)) then
                        check2(j,i) = check2(j,i)/sqrt(o_mass)
                    else if ((j.gt.45).and.(j.lt.49)) then
                        check2(j,i) = check2(j,i)/sqrt(o_mass)
                    else 
                        check2(j,i) = check2(j,i)/sqrt(h_mass)
                    endif
                enddo
            enddo
            do i =1,nmax
                do j = 1,ndim
                    check2(j,i) = check2(j,i) + coord(j)
                enddo
            enddo
            call calc_r(nmax,check2,nmono,r)
            call calc_psi(nmax,r,bend,nmono,psit,dr,first,second)
            do i = 1,nmax
                call calc_psi_ho(freq(:),sig,psi(i),check1(:,i))
                psi(i) = psi(i)/psi_max
            enddo
            do i = 1,nmax
                do j = 1,nmono
                    if ((r(1,j,i).ge.roh(nparam)).or.(r(2,j,i).ge.
     1              roh(nparam))) then
                        psi(i) = 0.d0
                        psit(i) = 0.d0
                    endif
                enddo
            enddo
            do i = 1,nmax
                if ((psi(i).gt.1e-10).and.(abs(psit(i)).gt.1e-12)) then
                    it = it + 1
                    psips(:,it) = check2(:,i)
                    if (it.eq.nslab) then
                        exit
                    endif
                endif
            enddo
        enddo
        return
        end subroutine

        subroutine calc_psi_ho(freq,sig,psi,q)
        implicit real*8(a-h,o-z)
        parameter (nvib=48)
        dimension freq(nvib),q(nvib)
        pi = dacos(-1.d0)
        psi = 1.d0
        do i = 1,nvib
            psi = psi*((freq(i)/(pi*sig**2))**(0.25)*exp(-0.5*freq(i)*
     1      (q(i)/sig)**2))
        enddo
        return
        end subroutine
