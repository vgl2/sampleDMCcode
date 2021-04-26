        subroutine calc_structures(dseed,nfinal,psips)
        implicit real*8(a-h,o-z)
        parameter (nmax=70000)
        parameter (ndim=54)
        parameter (nvib=48)
        parameter (nsymm=7)
        parameter (natoms=18)
        dimension freq(nvib,nsymm),evec(ndim,nvib,nsymm),
     1  coord(ndim,nmax,nsymm),x(ndim,nsymm),ngood(nsymm),
     1  psips(ndim,nmax),coord2(3,natoms,nmax),y(3,natoms,nsymm)
C		This code reads in the minimum energy structure of all of the hexamers 
C		and their associated harmonic frequencies and eigenvectors and generates
C		a nfinal number of walkers in random configurations. 
C		Inputs:
C		dseed = random number for seed
C		Outputs:
C		nfinal = number of final number of walkers to start calculation
C		psips = number of dimensions, number of nfinal walkers randomly configured
C		to begin calculation. Cartesian coordinates in units of bohr.
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
                    psips(k,ip) = coord(k,j,i)
                enddo
            enddo
        enddo
        do while (ip.lt.nfinal)
            do i = 1,nsymm
                do j = 1,ngood(i)
                    ip = ip + 1
                    do k = 1,ndim
                        psips(k,ip) = coord(k,j,i)
                    enddo
                    if (ip.eq.nfinal) then
                        exit
                    endif
                enddo
                if (ip.eq.nfinal) then
                    exit
                endif
            enddo
        enddo
        return
        end subroutine
