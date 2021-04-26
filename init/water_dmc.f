        program qmc
        use omp_lib
        implicit real*8 (a-h,o-z)
C
C

        parameter (ndim =54)
        parameter (nmax = 70000)
        parameter (nparam=500)
        parameter (nvib=18)
        parameter (nvibw=3)
        parameter (nmono=6)
        parameter (ntmax=50)
        parameter (nsymm=8)
        common/pts/psips(ndim,nmax,nsymm),t(nmax),v(nmax)
        common/kin/weight(ndim),d(ndim)
        common/pcoef/alpha,c2,re,rea,c1
        common/wa1/roh(nparam),roh1(nparam),wave_oh(nparam),
     1  feed_oh(nparam)
        common/wa2/rod(nparam),rod1(nparam),wave_od(nparam),
     1  feed_od(nparam)
        common/wav/psi_old(nmax),r(nvibw,nmono,nmax),vmin,sig_bend
        common/der/drdx(nvib,ndim,nmax),dr2dx2(nvib,ndim,nmax)
        dimension step(nmax*ndim),fate(nmax),
     1  q_force_new(nvibw,nmono,nmax),psi_new(nmax),drift(nmax*ndim),
     1  g_top(nmax),old_psips(ndim,nmax,nsymm),psi_sq_new(nmax),
     1  psi_sq_old(nmax),g_bottom(nmax),w(nmax),acceptance_ratio(nmax),
     1  acceptance(nmax),psi_old2(nmax),q_f_old(ndim,nmax),
     1  q_f_new(ndim,nmax),t_old(nmax),angle(nmax),v_old(nmax),
     1  r_old(nvibw,nmono,nmax),drdx_old(nvib,ndim,nmax),
     1  dr2dx2_old(nvib,ndim,nmax),q_force_old(nvibw,nmono,nmax),
     1  first(nvib,nmax),second(nvib,nmax),first_new(nvib,nmax),
     1  second_new(nvib,nmax),dseed(ntmax),step_chunk(ndim*nmax,ntmax),
     1  n0(nsymm),n(nsymm),ntot(nsymm)
        character (len=1024) energy,wavefunction(nsymm),wf
C
        read(*,*) dt
        read(*,*) endtime
        read(*,*) dseedi
        read(*,*) energy
        open (unit=12,file=energy,status='unknown')
C
C
C		This code is an importance sampled DMC code where we will read the 
C		timestep (dt), time propagated (endtime), name of file where the coordinates of 
C		the walkers will be stored (wavefunction) and the file where the energy of each 
C		timestep. (energy). Example below is shown for water hexamer (H2O)_6
C		Order of operation: simple details below but important details can be found in 
C		following publications: 
C		1. Reynolds, P. J.; Ceperley, D. M.; Alder, B. J.; Lester, W. A. 
C		Fixed-Node Quantum Monte Carlo for Molecules. J. Chem. Phys. 1982, 77, 5593−5603.
C		2.  Lee, V. G. M.; McCoy, A. B. An Efficient Approach for Studies of Water Clusters 
C		Using Diffusion Monte Carlo. J. Phys. Chem. A 2019, 123, 8063−8070
C		Steps
C		1. Read initial files (trial wave functions and generate randomly displaced 
C		geometries based on the harmonic oscillator wave function) for each individual 
C		well.
C		2. Displace walker based on random displacement and value of trial wave function.
C		3. Calculate local energy
C		4. Repeat until time=endtime
C
C		Input:
C		Input file regarding dt, endtime,random seed, where to store energy files
C		Output:
C		Wave function: 
C		number of walkers at time t,
C		ndim, nwalker coordiantes in units of bohr for each minima based on start*.dat
C		energy: energy of simulation per timestep
C		time, energy, number of walkers, acceptance value
C       initialization
        sig_stretch = 1.d0
        sig_bend = 1.d0
        pi = dacos(-1.d0)
        call initial_roh_wave(sig_stretch)
        call initial_rod_wave(sig_stretch)
        nthread = omp_get_max_threads()
        h_weight = 1.00782503223d0
        d_weight = 2.01410177812d0
        o_weight = 15.99491461957d0
        re = 1.81005d0
        rea = (104.508)*(pi/180.d0)
        c1 = 1655.4418464493720/219474.6
        c2 = 3887.5378889312051/219474.6
        g = (1/(h_weight*(re**2)))+(1/(h_weight*(re**2)))+(1/(o_weight*
     1  (re**2)))+(1/(o_weight*(re**2)))-((2*dcos(rea))/
     1  (o_weight*re*re))
        g = (1/g)*1822.88852962d0
        alpha =g*c1
        do i = 1,nsymm
            write(wavefunction(i),10) i
10          format('start',I0,'.dat')
        enddo
        call setup(dseedi,n0)
        do l = 1,nsymm
            n(l) = n0(l)
            time = 0.
            nborn = 0
            ndie = 0
            pi = dacos(-1.d0)
            call calc_r(n(l),psips(:,:,l),nmono,r)
            call calc_psi(n(l),r,sig_bend,nmono,psi_old,q_force_old,
     1      first,second)
c       calculate b matrix to get derivative in cartesian coordinates
!$omp parallel
!$omp do
            do i = 1,n(l)
                call calc_b_mat(psips(:,i,l),q_force_old(:,:,i),
     1          q_f_old(:,i),nmono,drdx(:,:,i),dr2dx2(:,:,i))
            enddo
!$omp end parallel
            drdx_old = drdx
            dr2dx2_old = dr2dx2
            veff = vavg(n0(l),first,second,psips(:,:,l))
            do i = 1,n(l)
                t_old(i) = t(i)
                v_old(i) = v(i)
            enddo
            vbar = veff
            alp = 0.5/dt
            tot_finish = 0.d0
            do while (time.lt.endtime)
                ntot(l) = n(l)
                time = time + dt
                vtot = 0.
                sig = sqrt(dt)
                call ggubs(dseedi,nthread,dseed)
!$omp parallel
!$omp do
                do i = 1,nthread
                    dseed(i) = dseed(i)*10**9.9
                    call gauss(0.d0,sig,step_chunk(:,i),
     1              (ndim*nmax)/nthread,dseed(i))
                enddo
!$omp end parallel
                it = 0
                do i = 1,nthread
                    do j = 1,(nmax*ndim)/nthread
                        it = it + 1
                        step(it) = step_chunk(j,i)
                        if (it.gt.nmax*ndim) then
                            exit
                        endif
                    enddo
                    if (it.gt.nmax*ndim) then
                        exit
                    endif
                enddo
                naccept = 0
                ip = 0
c           Calculate drift term and move the walker
                do i = 1,n(l)
                    do j = 1,ndim
                        ip = ip + 1
                        step(ip) = step(ip)/sqrt(weight(j))
                        old_psips(j,i,l) = psips(j,i,l)
                        drift(ip) = (q_f_old(j,i)*dt)/(2*weight(j))
                        step(ip) = step(ip) + drift(ip)
                        psips(j,i,l) = psips(j,i,l) + step(ip)
                    enddo
                enddo
                call ggubs(dseedi,nmax,acceptance)
C               Calculate new Psi_t term to determine if walker can make the
C               move.
                do i = 1,n(l)
                    psi_new(i) = 1.d0
                    g_top(i) = 1.d0
                    g_bottom(i) = 1.d0
                enddo
                call calc_r(n(l),psips(:,:,l),nmono,r)
                call calc_psi(n(l),r,sig_bend,nmono,psi_new,q_force_new,
     1          first_new,second_new)
!$omp parallel
!$omp do
                do i = 1,n(l)
                    call calc_b_mat(psips(:,i,l),q_force_new(:,:,i),
     1              q_f_new(:,i),nmono,drdx(:,:,i),dr2dx2(:,:,i))
                enddo
!$omp end parallel
                do i = 1,n(l)
                    psi_sq_old(i) = psi_old(i)**2
                    psi_sq_new(i) = psi_new(i)**2
                enddo
C           Calculating the green's function value to determine if
C           movement will be made
!$omp parallel
!$omp do
                do i = 1,n(l)
                    do j = 1,ndim
                        g_top(i) = g_top(i)*exp(-((old_psips(j,i,l)-
     1                  psips(j,i,l)-(q_f_new(j,i)*dt)/
     1                  (2*weight(j)))**2)/((2*dt)/weight(j))) 
                        g_bottom(i) = g_bottom(i)*exp(-((psips(j,i,l)-
     1                  old_psips(j,i,l)-(q_f_old(j,i)*dt)/
     1                  (2*weight(j)))**2)/((2*dt)/weight(j)))
                    enddo
                enddo
!$omp end parallel
                do i = 1,n(l)
                    w(i)=(psi_sq_new(i)*g_top(i))/(psi_sq_old(i)*
     1              g_bottom(i))
                    acceptance_ratio(i) = min(w(i),1.d0)
C                   If movement is accepted then copy everything the quantum
C                   force, and value of psi along with the walker, if not,
C                   then change the walker back to its previous position.
                    if (acceptance_ratio(i).gt.acceptance(i)) then
                        naccept = naccept + 1
                        first(:,i) = first_new(:,i)
                        second(:,i) = second_new(:,i)
                        q_f_old(:,i) = q_f_new(:,i)
                        psi_old(i) = psi_new(i)
                    else 
                        drdx(:,:,i) = drdx_old(:,:,i)
                        dr2dx2(:,:,i) = dr2dx2_old(:,:,i)
                        psips(:,i,l) = old_psips(:,i,l)
                        first_new(:,i) = first(:,i)
                        second_new(:,i) = second(:,i)
                    endif
                enddo
!$omp parallel
!$omp do
                do i = 1,n(l)
                    drdx_old(:,:,i) = drdx(:,:,i)
                    dr2dx2_old(:,:,i) = dr2dx2(:,:,i)
                enddo
!$omp end parallel
                neff = n(l)
C           Calculate the new dt effective to use for branching.
                dt_eff = (dfloat(naccept)/dfloat(neff))*dt
                nborn = 0
                ndie = 0
                call ggubs(dseedi,nmax,fate)
                n(l) = ntot(l)
                call calc_t(n(l),psips(:,:,l),psi_old,sig_bend,nmono,
     1          first,second,t)
                do i = 1,n(l)
                    do j = 1,ndim
                        psips(j,i,l) = psips(j,i,l)*0.52917721067d0
                    enddo
                enddo
!$omp parallel
!$omp do
                do i = 1,n(l)
                    call calcpot(nmono,v(i),psips(:,i,l))
                    v(i) = v(i)/627.509474d0
                enddo
!$omp end parallel
                do i = 1,n(l)
                    do j = 1,ndim
                        psips(j,i,l) = psips(j,i,l)/0.52917721067d0
                    enddo
                enddo
                do i = 1,n(l)
                    etot = t(i) + v(i)
                    if (v(i).lt.vmin) then
                        etot = 40.d0
                    endif
                    do j = 1,nmono
                        if ((r(1,j,i).ge.roh(nparam)).or.(r(2,j,i).ge.
     1                  roh(nparam))) then
                            etot = 40.d0
                        endif
                    enddo
                    dv = etot - veff
                    if (dv.lt.0) then
                        Pb = exp(-dv*dt) -1.
                        nb = int(Pb)
                        vtot = vtot + etot
                        if (fate(i).le.Pb-float(nb)) then
                            nb = nb + 1
                            nborn = nborn + 1
                        endif
                        do k = 1,nb
                            ntot(l) = ntot(l) + 1
                            psi_sq_old(ntot(l)) = psi_sq_old(i)
                            psi_old(ntot(l)) = psi_old(i)
                            t_old(ntot(l)) = t_old(i)
                            t(ntot(l)) = t(i)
                            v_old(ntot(l)) = v_old(i)
                            v(ntot(l)) = v(i)
                            drdx(:,:,ntot(l)) = drdx(:,:,i)
                            dr2dx2(:,:,ntot(l)) = dr2dx2(:,:,i)
                            do j = 1,ndim
                                psips(j,ntot(l),l) = psips(j,i,l)
                                q_f_old(j,ntot(l)) = q_f_old(j,i)
                            enddo
                            do j = 1,nvib
                                first(j,ntot(l)) = first(j,i)
                                second(j,ntot(l)) = second(j,i)
                            enddo
                            do j = 1,nmono
                                do m = 1,nvibw
                                    q_force_old(m,j,ntot(l))=
     1                              q_force_old(m,j,i)
                                enddo
                            enddo
                            vtot = vtot + etot
                        enddo
                    else
                        if (dv.gt.0) then
                            Pd = 1. - exp(-dv*dt)
                            if (fate(i).ge.Pd) then
                                vtot = vtot + etot
                            else
                                psips(1,i,l) = -999999.
                                ndie = ndie + 1
                            endif
                        endif
                    endif
                enddo
                n(l) = 0
                do i = 1,ntot(l)
                    if (psips(1,i,l).gt.-99998.) then
                        n(l) = n(l) + 1
                        do j = 1,ndim 
                            psips(j,n(l),l) = psips(j,i,l)
                            q_f_old(j,n(l)) = q_f_old(j,i)
                        enddo
                        do j = 1,nmono
                            do k = 1,nvibw
                                q_force_old(k,j,n(l))=q_force_old(k,j,i)
                            enddo
                        enddo
                        do j = 1,nvib
                            first(j,n(l)) = first(j,i)
                            second(j,n(l)) = second(j,i)
                        enddo
                        drdx(:,:,n(l)) = drdx(:,:,i)
                        dr2dx2(:,:,n(l)) = dr2dx2(:,:,i)
                        psi_sq_old(n(l)) = psi_sq_old(i)
                        psi_old(n(l)) = psi_old(i)
                        t_old(n(l)) = t_old(i)
                        t(n(l)) = t(i)
                        v_old(n(l)) = v_old(i)
                        v(n(l)) = v(i)
                    endif
                enddo
                veff_old = veff
                ntot(l) = n(l)
                vtot = vtot /dfloat(n(l))
                vbar = vtot
                veff = vbar + alp*log(dfloat(n0(l))/dfloat(n(l)))
                write(12,*)time,veff*219474.6,n(l),dfloat(naccept)/
     1          dfloat(neff)
c            print*, time,veff*219474.6,n(l),dfloat(naccept)/dfloat(neff)
                if (time.eq.endtime) then
                    open(unit=20,file=wavefunction(l),
     1              access='sequential',status='unknown')
                    write(20,*) n(l),dseedi,endtime
                    do i =1,n(l)
                        write(20,*) (psips(j,i,l),j=1,ndim)
                    enddo
                endif
            enddo
        enddo
        stop 
        end
C
C
        subroutine setup(dseed,n0)
        implicit real*8(a-h,o-z)
        parameter (ndim=54)
        parameter (nmax=70000)
        parameter (natoms=18)
        parameter (nsymm=8)
        common/pts/psips(ndim,nmax,nsymm),t(nmax),v(nmax)
        common/kin/weight(ndim),d(ndim)
        common/coef/s2pi
        common/pcoef/alpha,c2,re,rea,c1
        dimension n0(nsymm)
C
C		Creates an initial randomized configuration based on one
C		minima.
C		Inputs: 
C		dseed = random number
C		n0 = initial number of walkers (disregarding the number of 
C		minima (i.e. (50 000 walkers, not 50 000/number of minima 
C		walkers))
C       initialization
        call calc_structures(dseed,n0,psips)
        pi = dacos(-1.d0)
C       read the d's
C       m represents weight
        do j = 1,ndim
            if (j.lt.4) then
                weight(j) = 15.99491461957d0
            else if ((j.gt.9).and.(j.lt.13)) then
                weight(j) = 15.99491461957d0
            else if ((j.gt.18).and.(j.lt.22)) then
                weight(j) = 15.99491461957d0
            else if ((j.gt.27).and.(j.lt.31)) then
                weight(j) = 15.99491461957d0
            else if ((j.gt.36).and.(j.lt.40)) then
                weight(j) = 15.99491461957d0
            else if ((j.gt.45).and.(j.lt.49)) then
                weight(j) = 15.99491461957d0
            else 
                weight(j) = 1.00782503223d0
            endif
            weight(j) = weight(j)*1822.88852962d0
        enddo
C
C       calculate the constants
        pi = dacos(-1.d0)
        s2pi = 1./sqrt(2.*pi)
        return
        end
C
C    
        function vavg(n,first,second,psips)
        use omp_lib
        implicit real*8 (a-h,o-z)
        parameter (ndim=54)
        parameter (nvib=18)
        parameter (nvibw=3)
        parameter (nmax=70000)
        parameter (nparam=500)
        parameter (nmono=6)
        common/wav/psi_old(nmax),r(nvibw,nmono,nmax),vmin,sig_bend
        common /wa1/roh(nparam),roh1(nparam),wave_oh(nparam),
     1  feed_oh(nparam)
        common /wa2/rod(nparam),rod1(nparam),wave_od(nparam),
     1  feed_od(nparam)
        common/pcoef/alpha,c2,re,rea,c1
        common/der/drdx(nvib,ndim,nmax),dr2dx2(nvib,ndim,nmax)
        dimension eloc(nmax),first(nvib,nmax),second(nvib,nmax),
     1  psips(ndim,nmax),t(nmax),v(nmax)
C
        v = 0.
        call calc_t(n,psips,psi_old,sig_bend,nmono,first,second,t)
        do i = 1,n
            do j = 1,ndim
                psips(j,i) = psips(j,i)*0.52917721067d0
            enddo
        enddo
!$omp parallel
!$omp do
        do i = 1,n
            call calcpot(nmono,v(i),psips(:,i))
            v(i) = v(i)/627.509474d0
        enddo
!$omp end parallel
        do i = 1,n
            do j = 1,ndim
                psips(j,i) = psips(j,i)/0.52917721067d0
            enddo
        enddo
        vmin = -45.94d0/627.509474d0
        vavg = 0.d0
        do i = 1,n
            eloc(i) = t(i)+v(i)
            vavg = vavg + eloc(i)
        enddo
        vavg = vavg/dfloat(n)
        print *, vavg*219474.6
        return
        end
C
C
