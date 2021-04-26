        subroutine calc_nm(y,freq,evec)
        use omp_lib
        implicit real*8(a-h,o-z)
        parameter (nsymm=7)
        parameter (natoms=18)
        parameter (ndim=54)
        parameter (nvib=(3*natoms)-6)
        dimension x(3,natoms,nsymm),y(ndim,nsymm),freq(nvib,nsymm),
     1  evec(ndim,nvib,nsymm)
C		Calculate the harmonic frequencies of the water hexamer based on the initial geometries 
C		Assumes structures are minima structures 
C		(frequncies will not be useful if structures are not minima)
C		Outputs: 
C		y = nsymm geometries of water hexamer minimum energy structures 
C		in bohr.
C		freq: 3*natoms - 6 positive frequencies of the molecular vibrations.
C		evec: 3*natoms - 6 eigenvectors concerning the molecular vibrations
C		for the associated frequencies.
        open(unit=499,file='hexamer-initial.dat',status='old')
        do i = 1,nsymm
            do j = 1,natoms
                read(499,*) (x(k,j,i),k=1,3)
            enddo
        enddo
        close (499)
        do i = 1,nsymm
            l = 0
            do j = 1,natoms
                do k = 1,3
                    l = l + 1
                    y(l,i) = x(k,j,i)/0.52917721067d0
                enddo
            enddo
        enddo
!$omp parallel
!$omp do
        do i = 1,nsymm
            call calc_harmonic_frequencies(y(:,i),ndim,freq(:,i),
     1      evec(:,:,i))
        enddo
!$omp end parallel
        end subroutine
