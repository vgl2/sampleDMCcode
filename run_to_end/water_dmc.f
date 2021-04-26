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
        common/pts/psips(ndim,nmax),sig_bend,t(nmax),v(nmax)
        common/kin/weight(ndim),d(ndim)
        common/pcoef/alpha,c2,re,rea,c1
        common/wa1/roh(nparam),roh1(nparam),wave_oh(nparam),
     1  feed_oh(nparam)
        common/wa2/rod(nparam),rod1(nparam),wave_od(nparam),
     1  feed_od(nparam)
        common/wav/psi_old(nmax),r(nvibw,nmono,nmax),vmin
        common/der/drdx(nvib,ndim,nmax),dr2dx2(nvib,ndim,nmax)
        dimension step(nmax*ndim),fate(nmax),
     1  q_force_new(nvibw,nmono,nmax),psi_new(nmax),drift(nmax*ndim),
     1  g_top(nmax),old_psips(ndim,nmax),psi_sq_new(nmax),
     1  psi_sq_old(nmax),g_bottom(nmax),w(nmax),acceptance_ratio(nmax),
     1  acceptance(nmax),psi_old2(nmax),q_f_old(ndim,nmax),
     1  q_f_new(ndim,nmax),t_old(nmax),angle(nmax),v_old(nmax),
     1  r_old(nvibw,nmono,nmax),drdx_old(nvib,ndim,nmax),
     1  dr2dx2_old(nvib,ndim,nmax),q_force_old(nvibw,nmono,nmax),
     1  first(nvib,nmax),second(nvib,nmax),first_new(nvib,nmax),
     1  second_new(nvib,nmax),dseed(ntmax),step_chunk(ndim*nmax,ntmax)
        character (len=1024) energy,wavefunction,wf
C
        read(*,*) dt
        read(*,*) endtime
        read(*,*) wavefunction
        read(*,*) energy
        open (unit=12,file=energy,status='unknown')
        open (unit=20,file=wavefunction,status='unknown',
     1  form='unformatted')
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
C		1. Read initial files (trial wave functions and walker coordinates)
C		2. Displace walker based on random displacement and value of trial wave function.
C		3. Calculate local energy
C		4. Repeat until time=endtime
C
C		Input:
C		Input file regarding dt, endtime, where to store walkers and where to store energy files
C		Output:
C		Wave function: number of wave functions
C		number of walkers at time t, number of walkers at time 0, time taken
C		ndim, nwalker coordiantes in units of bohr
C		energy: energy of simulation per timestep
C		time, energy, number of walkers, acceptance value
C       initialization
        sig_stretch = 1.d0
        sig_bend = 1.d0
        pi = dacos(-1.d0)
c		read in the trial wave functions
        call initial_roh_wave(sig_stretch)
        call initial_rod_wave(sig_stretch)
        nthread = omp_get_max_threads()
        h_weight = 1.00782503223d0
        d_weight = 2.01410177812d0
        o_weight = 15.99491461957d0
        re = 1.81005d0
        rea = (104.508)*(pi/180.d0)
C		define the frequnecy for the vibrations (c1 denotes the HOH bend in atomic units
C		c2 denotes the OH stretch in atomic units, c2 is depricated because using wave functions
C		generated from DVR calculations)
C		g is the g matrix elements that the hoh bend motion in electron mass units.
        c1 = 1655.4418464493720/219474.6
        c2 = 3887.5378889312051/219474.6
        g = (1/(h_weight*(re**2)))+(1/(h_weight*(re**2)))+(1/(o_weight*
     1  (re**2)))+(1/(o_weight*(re**2)))-((2*dcos(rea))/
     1  (o_weight*re*re))
        g = (1/g)*1822.88852962d0
        alpha =g*c1
        call setup(dseedi,n0,n,time)
        nborn = 0
        ndie = 0
        pi = dacos(-1.d0)
C		Takes cartesian coordinates and transforms into internal coordinates for trial wave function
        call calc_r(n,psips,nmono,r)
        call calc_psi(n,r,sig_bend,nmono,psi_old,q_force_old,first,
     1  second)
c       calculate b matrix to get derivative in cartesian coordinates
!$omp parallel
!$omp do
        do i = 1,n
            call calc_b_mat(psips(:,i),q_force_old(:,:,i),q_f_old(:,i),
     1      nmono,drdx(:,:,i),dr2dx2(:,:,i))
        enddo
!$omp end parallel
        drdx_old = drdx
        dr2dx2_old = dr2dx2
        veff = vavg(n0,n,dt,first,second)
        do i = 1,n
            t_old(i) = t(i)
            v_old(i) = v(i)
        enddo
        vbar = veff
        alp = 0.5/dt
        write(20) 30
        tot_finish = 0.d0
        do while (time.lt.endtime)
            ntot = n
            time = time + dt
            vtot = 0.
            sig = sqrt(dt)
            call ggubs(dseedi,nthread,dseed)
!$omp parallel
!$omp do
            do i = 1,nthread
                dseed(i) = dseed(i)*10**9.9
                call gauss(0.d0,sig,step_chunk(:,i),(ndim*nmax)/nthread,
     1          dseed(i))
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
            do i = 1,n
                do j = 1,ndim
                    ip = ip + 1
                    step(ip) = step(ip)/sqrt(weight(j))
                    old_psips(j,i) = psips(j,i)
                    drift(ip) = (q_f_old(j,i)*dt)/(2*weight(j))
                    step(ip) = step(ip) + drift(ip)
                    psips(j,i) = psips(j,i) + step(ip)
                enddo
            enddo
            call ggubs(dseedi,nmax,acceptance)
C           Calculate new Psi_t term to determine if walker can make the
C           move.
            do i = 1,n
                psi_new(i) = 1.d0
                g_top(i) = 1.d0
                g_bottom(i) = 1.d0
            enddo
            call calc_r(n,psips,nmono,r)
            call calc_psi(n,r,sig_bend,nmono,psi_new,q_force_new,
     1      first_new,second_new)
!$omp parallel
!$omp do
            do i = 1,n
                call calc_b_mat(psips(:,i),q_force_new(:,:,i),
     1          q_f_new(:,i),nmono,drdx(:,:,i),dr2dx2(:,:,i))
            enddo
!$omp end parallel
            do i = 1,n
                psi_sq_old(i) = psi_old(i)**2
                psi_sq_new(i) = psi_new(i)**2
            enddo
C           Calculating the green's function value to determine if
C           movement will be made
!$omp parallel
!$omp do
            do i = 1,n
                do j = 1,ndim
                    g_top(i) = g_top(i)*exp(-((old_psips(j,i)-psips(j,i)
     1              -(q_f_new(j,i)*dt)/(2*weight(j)))**2)/((2*dt)/
     1              weight(j))) 
                    g_bottom(i) = g_bottom(i)*exp(-((psips(j,i)-
     1              old_psips(j,i)-(q_f_old(j,i)*dt)/(2*weight(j)))**2)/
     1              ((2*dt)/weight(j)))
                enddo
            enddo
!$omp end parallel
            do i = 1,n
               w(i)=(psi_sq_new(i)*g_top(i))/(psi_sq_old(i)*g_bottom(i))
                acceptance_ratio(i) = min(w(i),1.d0)
C               If movement is accepted then copy everything the quantum
C               force, and value of psi along with the walker, if not,
C               then change the walker back to its previous position.
                if (acceptance_ratio(i).gt.acceptance(i)) then
                    naccept = naccept + 1
                    first(:,i) = first_new(:,i)
                    second(:,i) = second_new(:,i)
                    q_f_old(:,i) = q_f_new(:,i)
                    psi_old(i) = psi_new(i)
                else 
                    drdx(:,:,i) = drdx_old(:,:,i)
                    dr2dx2(:,:,i) = dr2dx2_old(:,:,i)
                    psips(:,i) = old_psips(:,i)
                    first_new(:,i) = first(:,i)
                    second_new(:,i) = second(:,i)
                endif
            enddo
!$omp parallel
!$omp do
            do i = 1,n
                drdx_old(:,:,i) = drdx(:,:,i)
                dr2dx2_old(:,:,i) = dr2dx2(:,:,i)
            enddo
!$omp end parallel
C           Calculate the new dt effective to use for branching.
            neff = n
            dt_eff = (dfloat(naccept)/dfloat(neff))*dt
            nborn = 0
            ndie = 0
            call ggubs(dseedi,nmax,fate)
            n = ntot
            call calc_t(n,psips,psi_old,sig_bend,nmono,first,second,t)
            do i = 1,n
                do j = 1,ndim
                    psips(j,i) = psips(j,i)*0.52917721067d0
                enddo
            enddo
C			Calculate the potential Note: code for potential can be swapped out
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
C			Branching step (determine which walkers are duplicated or removed)
            it = 0
            ncount = 0
            do i = 1,n
                etot = t(i) + v(i)
                if (v(i).lt.vmin) then
                    ncount =ncount + 1
                    etot = 40.d0
                endif
                do j = 1,nmono
                    if ((r(1,j,i).ge.roh(nparam)).or.(r(2,j,i).ge.
     1              roh(nparam))) then
                        it = it + 1
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
                        ntot = ntot + 1
                        psi_sq_old(ntot) = psi_sq_old(i)
                        psi_old(ntot) = psi_old(i)
                        t(ntot) = t(i)
                        v(ntot) = v(i)
                        drdx(:,:,ntot) = drdx(:,:,i)
                        dr2dx2(:,:,ntot) = dr2dx2(:,:,i)
                        do j = 1,ndim
                            psips(j,ntot) = psips(j,i)
                            q_f_old(j,ntot) = q_f_old(j,i)
                        enddo
                        do j = 1,nvib
                            first(j,ntot) = first(j,i)
                            second(j,ntot) = second(j,i)
                        enddo
                        do j = 1,nmono
                            do l = 1,nvibw
                                q_force_old(l,j,ntot)=q_force_old(l,j,i)
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
                            psips(1,i) = -999999.
                            ndie = ndie + 1
                        endif
                    endif
                endif
            enddo
            nt = n
            birth_rate = dfloat(nborn)/dfloat(nt)
            death_rate = dfloat(ndie)/dfloat(nt)
            n = 0
            do i = 1,ntot
                if (psips(1,i).gt.-99998.) then
                    n = n + 1
                    do j = 1,ndim 
                        psips(j,n) = psips(j,i)
                        q_f_old(j,n) = q_f_old(j,i)
                    enddo
                    do j = 1,nmono
                        do k = 1,nvibw
                            q_force_old(k,j,n) = q_force_old(k,j,i)
                        enddo
                    enddo
                    do j = 1,nvib
                        first(j,n) = first(j,i)
                        second(j,n) = second(j,i)
                    enddo
                    drdx(:,:,n) = drdx(:,:,i)
                    dr2dx2(:,:,n) = dr2dx2(:,:,i)
                    psi_sq_old(n) = psi_sq_old(i)
                    psi_old(n) = psi_old(i)
                    t(n) = t(i)
                    v(n) = v(i)
                endif
            enddo
C			Save walkers in this step here
            if ((time.gt.25000).and.mod(time,2500.).eq.0) then
                write(20) n,n0,time
                do i = 1,n
                    write(20) (psips(j,i),j=1,ndim)
                enddo
            endif
C			Calculate new veff which is Eref in publications
            veff_old = veff
            ntot = n
            vtot = vtot /dfloat(n)
            vbar = vtot
            veff = vbar + alp*log(dfloat(n0)/dfloat(n))
            write(12,*)time,veff*219474.6,n,dfloat(naccept)/dfloat(neff)
c            print*, time,veff*219474.6,n,dfloat(naccept)/dfloat(neff)
        enddo
        stop 
        end
C
C
        subroutine setup(dseed,n0,n,tot_time)
        implicit real*8(a-h,o-z)
        parameter (ndim=54)
        parameter (nmax=70000)
        parameter (natoms=18)
        parameter (nsymm=8)
        common/pts/psips(ndim,nmax),sig_bend,t(nmax),v(nmax)
        common/kin/weight(ndim),d(ndim)
        common/coef/s2pi
        common/pcoef/alpha,c2,re,rea,c1
        character (len=1024) wf,wavefunction(nsymm)
        dimension coord(ndim,nmax,nsymm),nsim(nsymm),dseedi(nsymm),
     1  endtime(nsymm)
C
C		Reads in from a previous run that separated all the minima and 
C		combines all of the results from one minima together into one
C		large simulation.
C		Inputs:
C		dseed = random number
C		n0 = number of walkers at t=0
C		n = number of walkers at current timestep
C		tot_time = last recorded timestep
C		Outputs:
C		weights = the masses of each atom 
C		psips = ndim,nwalker array combining all the start*.dat files.
C       initialization by reading in old wave functions
        do i = 1,nsymm
            write(wavefunction(i),10) i
10          format('start',I0,'.dat')
            open(unit=2,access='sequential',file=wavefunction(i),
     1      status='old')
            read(2,*) nsim(i),dseedi(i),endtime(i)
            do j = 1,nsim(i)
                read(2,*) (coord(k,j,i),k=1,ndim)
            enddo
        enddo
        tot_time = endtime(1)
        call cpu_time(time)
        time = time*10**9.9
        call ggubs(time,1,pick)
        ipick = int(nsymm*pick)
        if (ipick.eq.0) ipick=1
        dseed = dseedi(int(ipick))
        n0 = 50000
        n = 0
        do i = 1,nsymm
            do j = 1,nsim(i)
                n = n + 1
                do k = 1,ndim
                    psips(k,n) = coord(k,j,i)
                enddo
            enddo
        enddo
        pi = dacos(-1.d0)
C       Describe the masses of each of atoms in electron mass units
C		(in cartesian coordinates)
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
        function vavg(n0,n,dt,first,second)
        use omp_lib
        implicit real*8 (a-h,o-z)
        parameter (ndim=54)
        parameter (nvib=18)
        parameter (nvibw=3)
        parameter (nmax=70000)
        parameter (nparam=500)
        parameter (nmono=6)
        common/wav/psi_old(nmax),r(nvibw,nmono,nmax),vmin
        common/pts/psips(ndim,nmax),sig_bend,t(nmax),v(nmax)
        common /wa1/roh(nparam),roh1(nparam),wave_oh(nparam),
     1  feed_oh(nparam)
        common /wa2/rod(nparam),rod1(nparam),wave_od(nparam),
     1  feed_od(nparam)
        common/pcoef/alpha,c2,re,rea,c1
        common/der/drdx(nvib,ndim,nmax),dr2dx2(nvib,ndim,nmax)
        dimension eloc(nmax),first(nvib,nmax),second(nvib,nmax)
C
        alp = 0.5/dt
        v = 0.
C		Calculate (T\Psi_T)/(\Psi_T) here
        call calc_t(n,psips,psi_old,sig_bend,nmono,first,second,t)
        do i = 1,n
            do j = 1,ndim
                psips(j,i) = psips(j,i)*0.52917721067d0
            enddo
        enddo
C		Calculate (V\Psi_T)/(\Psi_T) here
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
C		Minimum of potential energy surface
        vmin = -45.94d0/627.509474d0
        vavg = 0.d0
        do i = 1,n
            eloc(i) = t(i)+v(i)
            vavg = vavg + eloc(i)
        enddo
        vavg = vavg/dfloat(n)+ alp*log(dfloat(n0)/dfloat(n))
        print *, vavg*219474.6
        return
        end
C
C
