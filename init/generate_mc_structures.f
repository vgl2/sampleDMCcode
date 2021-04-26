        subroutine calc_structures(dseed,ngood,psips)
        implicit real*8(a-h,o-z)
        parameter (nmax=70000)
        parameter (ndim=54)
        parameter (nvib=48)
        parameter (nsymm=8)
        parameter (natoms=18)
        dimension freq(nvib,nsymm),evec(ndim,nvib,nsymm),
     1  coord(ndim,nmax,nsymm),x(ndim,nsymm),ngood(nsymm),
     1  psips(ndim,nmax,nsymm),coord2(3,natoms,nmax),y(3,natoms,nsymm)
C		Tables the minimum energy structures of the water clusters and randomly 
C		sample a wave function that descibes the vibrations using the 
C		harmonic oscillator to have a starting point for these calculations.
C		Inputs:
C		dseed = random number
C		ngood = total number of walkers within a particular set of 
C		configurations (nwalkers/nsymm)
C		Outputs:
C		psips = 3*natoms, nwalker/nsymm array that describes the 
C		random set of geometries sampled by a harmonic oscillator
C		for a given configuration.
        nfinal = 50000
        open(unit=499,file='hexamer-initial.dat',status='old')
        open(unit=599,file='all-h_start.dat',status='old')
        do i = 1,nsymm
            do j = 1,natoms
                read(499,*) (y(k,j,i),k=1,3)
            enddo
        enddo
        close (499)
        do i = 1,nsymm
            do j = 1,nvib
                read(599,*) freq(j,i)
                read(599,*) (evec(k,j,i),k=1,ndim)
            enddo
        enddo
        do i = 1,nsymm
            l = 0
            do j = 1,natoms
                do k = 1,3
                    l = l + 1
                    x(l,i) = y(k,j,i)/0.52917721067d0
                enddo
            enddo
        enddo
        do i = 1,nsymm
            call mc_sample(dseed,x(:,i),nsymm,nfinal,freq(:,i),
     1      evec(:,:,i),ngood(i),coord(:,:,i))
        enddo
        ip = 0
        do i = 1,nsymm
            do j = 1,ngood(i)
                ip = ip + 1
                do k = 1,ndim
                    psips(k,j,i) = coord(k,j,i)
                enddo
            enddo
        enddo
        return
        end subroutine
