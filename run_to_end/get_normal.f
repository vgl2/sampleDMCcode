        program get_normal
        implicit real*8 (a-h,o-z)
        parameter (nsymm=7)
        parameter (ndim=54)
        parameter (nvib=ndim-6)
        dimension x(ndim,nsymm),freq(nvib,nsymm),evec(ndim,nvib,nsymm)
C		Reads in minimum energy structures of water hexamer and calculates the 
C		harmonic frequencies and eigenvectors associated to frequency of 
C		water cluster.
C		Outputs:
C		freq: 3*natoms - 6 positive frequencies of the molecular vibrations.
C		evec: 3*natoms - 6 eigenvectors concerning the molecular vibrations
C		for the associated frequencies.
        open(unit=8,file='hexamer_start.dat',status='unknown')
        call calc_nm(x,freq,evec)
        write(8,*) nsymm
        do i = 1,nsymm
            do j = 1,nvib
                write(8,*) freq(j,i)
                write(8,*) (evec(k,j,i),k=1,ndim)
            enddo
        enddo
        end program
