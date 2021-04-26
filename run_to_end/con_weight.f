        subroutine cont_weight(psips,dseed,n,n0,weight_tot,avg_etrial,dt
     1  ,endtime,holes)
        use omp_lib
C       Propogate my walkers to get weights based on continuous weighting 
C		used for purposes of descendant weighting. These weights will be used to 
C		calculate expectation values of multiplicative operators.
C		Inputs:
C		psips = ndim*nwalker list of coordinates of walker
C		dseed = random number for seed
C		n = number of walkers
C		n0 = initial number of walkers
C		dt = timestep of calculation
C		endtime = end of calculation in au.
C		Outputs:
C		weight_tot = total weight of all the walkers
C		avg_etrial = avg eref of the entire simulation
C		holes = coordinates of walkers that entered into holes of the potential 
        implicit real*8 (a-h,o-z)
        parameter (natoms=18)
        parameter (nmax=70000)
        parameter (nvib=18)
        parameter (ndimw=9)
        parameter (nvibw=3)
        parameter (nparam=500)
        parameter (nmono=6)
        parameter (ndim=54)
        common/kin/weight(ndim),d(ndim)
        common /wa1/roh1(nparam),wave_oh(nparam),feed_oh(nparam)
        common /wa2/rod1(nparam),wave_od(nparam),feed_od(nparam)
        common/der/drdx(nvib,ndim,nmax),dr2dx2(nvib,ndim,nmax)
        common/pcoef/alpha,c2,re,rea,c1
        dimension psips(ndim,nmax),etot_old(nmax),n_coord(nmax),
     1  step(nmax*ndim),con_weight(nmax),v(nmax),psi_old(nmax),
     1  v_initial(nmax),accept(nmax),psi_new(nmax),psi_sq_old(nmax),
     1  psi_sq_new(nmax),g_top(nmax),g_bottom(nmax),w(nmax),
     1  drift(nmax*ndim),weight_tot(nmax),q_f_old(ndim,nmax),
     1  g_diffusion(ndim),q_f_new(ndim,nmax),g_dim_top(ndim,nmax),
     1  g_top_total(nmax),g_dim_bottom(ndim,nmax),g_bottom_total(nmax),
     1  g_total(nmax),stepsize(ndim,nmax),etot(nmax),first(nvib,nmax),
     1  old_psips(ndim,nmax),dv(nmax),t(nmax),q_f_rold(nvibw,nmono,nmax)
     1  ,r_old(nvibw,nmono,nmax),r_new(nvibw,nmono,nmax),
     1  q_f_rnew(nvibw,nmono,nmax),t_old(nmax),v_old(nmax),
     1  second(nvib,nmax),first_new(nvib,nmax),second_new(nvib,nmax),
     1  drdx_new(nvib,ndim,nmax),dr2dx2_new(nvib,ndim,nmax)
        character(len=1024) holes
        open(unit=25,file=holes,status='unknown')
        pi = dacos(-1.d0)
        dseed = dseed*10**9
        tot = n
        first = n0
        do i = 1,n
            n_coord(i) = i
        enddo
        sig_bend = 1.d0
        call calc_r(n,psips,nmono,r_old)
C       Calculate Mass!!!
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
        do j = 1,ndim
            g_diffusion(j) =(2*dt)/(weight(j))
        enddo
        tote = 0.d0
C       Find total weight
        tot_descend = 0.d0
        time = 0
        do i = 1,n
            weight_tot(i) = 0.d0
            con_weight(i) = 1.d0
            tot_descend = tot_descend + con_weight(i)
        enddo
        do i = 1,nmax
            if (i.gt.n) then
                con_weight(i) = 0.d0
            endif
        enddo
        tot_run = endtime/dt
        alp = 0.5/dt
        call calc_psi(n,r_old,sig_bend,nmono,psi_old,q_f_rold,first,
     1  second)
        do i = 1,n
            call calc_b_mat(psips(:,i),q_f_rold(:,:,i),q_f_old(:,i),
     1      nmono,drdx(:,:,i),dr2dx2(:,:,i))
        enddo
        call calc_t(n,psips,psi_old,sig_bend,nmono,first,second,t)
        do i = 1,n
            do j= 1,ndim
                psips(j,i) = 0.52917721067*psips(j,i)
            enddo
        enddo
!$omp parallel
!$omp do
        do i = 1,n
            call calcpot(nmono,v_initial(i),psips(:,i))
            v_initial(i) = v_initial(i)/627.509474
        enddo
!$omp end parallel
        do i = 1,n
            do j= 1,ndim
                psips(j,i) = psips(j,i)/0.52917721067
            enddo
        enddo
        ip = 0
        do i = 1,n
            etot(i) = t(i)+v_initial(i)
        enddo
        do i = 1,n
            tote = tote + etot(i)
            etot_old(i) = etot(i)
        enddo
        vavg = tote/tot +alp*(log(dfloat(n0)/tot_descend))
        etrial = vavg
        vt =0
        sig = sqrt(dt)
        etrial_tot = 0.d0
        do while (time.lt.endtime)
            ntot = n
            naccept = 0
            time = time + dt
            call gauss(0.d0,sig,step,ndim*nmax,dseed)
            call ggubs(dseed,nmax,accept) 
            vt = 0.d0
            ip = 0
            naccept = 0
            do i = 1,n
                do j = 1,ndim
                    ip = ip + 1
                    drift(ip) = (q_f_old(j,i)*dt)/(2*weight(j))
                    step(ip) = step(ip)/sqrt(weight(j))
                    stepsize(j,i) = step(ip)
                    step(ip) = step(ip) + drift(ip)
                    old_psips(j,i) = psips(j,i)
                    psips(j,i) = psips(j,i)+ step(ip)
                enddo
            enddo
            call calc_r(n,psips,nmono,r_new)
            call calc_psi(n,r_new,sig_bend,nmono,psi_new,q_f_rnew,
     1      first_new,second_new)
            do i = 1,n
                call calc_b_mat(psips(:,i),q_f_rnew(:,:,i),q_f_new(:,i),
     1          nmono,drdx_new(:,:,i),dr2dx2_new(:,:,i))
            enddo
            do i = 1,n
                g_top_total(i) = 0.d0
                g_bottom_total(i) = 0.d0
            enddo
            do i = 1,n
                psi_sq_old(i) = psi_old(i)**2
                psi_sq_new(i) = psi_new(i)**2
                do j= 1,ndim
                    g_dim_top(j,i) = (old_psips(j,i)-psips(j,i)-
     1              ((dt*q_f_new(j,i))/(2*weight(j))))**2
                    g_top_total(i) = g_top_total(i) + 
     1              (g_dim_top(j,i)/g_diffusion(j))
                    g_dim_bottom(j,i) = (psips(j,i)-old_psips(j,i)-
     1              ((dt*q_f_old(j,i))/(2*weight(j))))**2
                    g_bottom_total(i) = g_bottom_total(i) + 
     1              (g_dim_bottom(j,i)/g_diffusion(j))
                enddo
            enddo
            do i = 1,n
                g_top(i) = exp(-g_top_total(i))
                g_bottom(i) = exp(-g_bottom_total(i))
                g_total(i) = g_top(i)/g_bottom(i)
                w(i) = (psi_sq_new(i)/psi_sq_old(i))*g_total(i)
                if (min(w(i),1.d0).lt.accept(i)) then
                    etot_old(i) = etot_old(i)
                    do j = 1,ndim
                        psips(j,i) = old_psips(j,i)
                        q_f_new(j,i) = q_f_old(j,i)
                    enddo
                    do j = 1,nmono
                        do k = 1,nvibw
                            r_new(k,j,i) = r_old(k,j,i)
                            q_f_rnew(k,j,i) = q_f_rold(k,j,i)
                        enddo
                    enddo
                    do j = 1,nvib
                        first_new(j,i) = first(j,i)
                        second_new(j,i) = second(j,i)
                    enddo
                else
                    do j = 1,nvib
                        first(j,i) = first_new(j,i)
                        second(j,i) = second_new(j,i)
                    enddo
                    naccept = naccept + 1
                    etot_old(i) = etot(i)
                    psi_old(i) = psi_new(i)
                    do j = 1,ndim
                        old_psips(j,i) = psips(j,i)
                        q_f_old(j,i) = q_f_new(j,i)
                    enddo
                    do j = 1,nmono
                        do k = 1,nvibw
                            r_old(k,j,i) = r_new(k,j,i)
                            q_f_rold(k,j,i) = q_f_rnew(k,j,i)
                        enddo
                    enddo
                    drdx(:,:,i) = drdx_new(:,:,i)
                    dr2dx2(:,:,i) = dr2dx2_new(:,:,i)
                endif
            enddo
            dt_eff = (dfloat(naccept)/dfloat(n))*dt
            call calc_t(n,psips,psi_old,sig_bend,nmono,first,second,t)
            do i = 1,n
                do j= 1,ndim
                    psips(j,i) = 0.52917721067*psips(j,i)
                enddo
            enddo
!$omp parallel
!$omp do
            do i = 1,n
                call calcpot(nmono,v(i),psips(:,i))
                v(i) = v(i)/627.509474
            enddo
!$omp end parallel
            do i = 1,n
                do j= 1,ndim
                   psips(j,i) = psips(j,i)/0.52917721067
                enddo
            enddo
            do i = 1,n
                etot(i) = t(i) + v(i)
c                etot_old(i) = t_old(i) + v_old(i)
            enddo
            do i = 1,n
                dv(i) = etot(i) - etrial
                con_weight(i) = con_weight(i)*exp(-dv(i)*dt)
            enddo
            if (time.lt.endtime) then
                ncount = 0
                do i =1,n
                    if (con_weight(i).lt.(1/tot)) then
                        ncount = ncount + 1
                        nummax = maxloc(con_weight,1)
                        n_coord(i) = n_coord(nummax)
                        con_weight(i) = con_weight(nummax)/2
                        con_weight(nummax) = con_weight(nummax)/2
                        psi_old(i) = psi_old(nummax)
                        psi_new(i) = psi_new(nummax)
                        etot(i) = etot(nummax)
                        etot_old(i) = etot_old(nummax)
                        do j = 1,ndim
                            psips(j,i) = psips(j,nummax)
                            old_psips(j,i) = old_psips(j,nummax)
                            q_f_old(j,i) = q_f_old(j,nummax)
                        enddo
                        drdx(:,:,i) = drdx(:,:,nummax)
                        dr2dx2(:,:,i) = dr2dx2(:,:,nummax)
                        do j = 1,nmono
                            do k = 1,nvibw
                                r_new(k,j,i) = r_new(k,j,nummax)
                                r_old(k,j,i) = r_old(k,j,nummax)
                                q_f_rnew(k,j,i) = q_f_rnew(k,j,nummax)
                                q_f_rold(k,j,i) = q_f_rold(k,j,nummax)  
                            enddo
                        enddo
                        do j = 1,nvib
                            first(j,i) = first(j,nummax)
                            second(j,i) = second(j,nummax)
                        enddo
                    endif
                enddo
            endif
            tot_weight = 0.d0
            do i = 1,n
                tot_weight = tot_weight + con_weight(i)
            enddo       
            vavg = 0.d0
            do i = 1,n
                vavg = vavg + con_weight(i)*etot(i)
            enddo
            veff = vavg/tot_weight
            etrial_old = etrial 
            etrial = veff + alp*(log(tot_descend/tot_weight))
c            print *, time,etrial*219474.6,tot_weight,dfloat(naccept)/
c     1      dfloat(n)
            etrial_tot = etrial_tot + etrial
        enddo
        avg_etrial = etrial_tot/tot_run
        do i = 1,n
            do j = 1,n
                if (n_coord(j).eq.i) then
                    weight_tot(i) = weight_tot(i) + con_weight(j)
                endif
            enddo
        enddo
        return
        end subroutine

