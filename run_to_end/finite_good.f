        subroutine calc_harmonic_frequencies(coord,ip,freq_keep,
     1  evec_keep)
        implicit real*8(a-h,o-z)
        parameter (natoms=18)
        parameter (nmono=6)
        parameter (nmax=60)
        dimension x(3,natoms),coord(natoms*3),test(5,5),
     1  force(natoms*3,natoms*3),weight(natoms*3),
     1  eigenvector(natoms*3,natoms*3),scratch(natoms*3,natoms*3),
     1  freq(natoms*3),eigenvalue(natoms*3),freq_keep(natoms*3),
     1  eigenvalue_keep(natoms*3),evec(3,natoms,nmax),
     1  evec_keep(3*natoms,3*natoms)
C		calculate a normal mode calculation based on 5 point finite
C		difference approach. 
C		Inputs:
C		coord: natoms*3 array in units of bohr of the Cartesian coordinates
C		of the molecule of choice (in this case water clusters are listed
C		in the following order: (O,H,H)*natoms
C		ip = total number of dimensions (3*natoms)
C		Outputs:
C		freq_keep: 3*natoms - 6 positive frequencies of the molecular vibrations.
C		evec_keep: 3*natoms - 6 eigenvectors concerning the molecular vibrations
C		for the associated frequencies.
        pi = dacos(-1.d0)
        dx = 0.001
        do i = 1,5
            do j = 1,5
                test(i,j) = 0.d0
            enddo
        enddo
        coord(:) = coord(:)*0.52917721067d0
        call calcpot(nmono,vmin,coord(:))
        vmin = vmin/627.509474d0
        coord(:) = coord(:)/0.52917721067d0
        do k = 1,ip
            do l = k,ip
                do i = 1,5
                    if (i.eq.1) then
                        coord(k) = coord(k) - 2*dx
                    else if (i.eq.2) then
                        coord(k) = coord(k) - dx
                    else if (i.eq.3) then
                        coord(k) = coord(k) 
                    else if (i.eq.4) then  
                        coord(k) = coord(k) + dx
                    else if (i.eq.5) then  
                        coord(k) = coord(k) + 2*dx
                    endif
                        do j = 1,5
                            if (j.eq.1) then
                                coord(l) = coord(l) + 2*dx
                                coord(:) = coord(:)*0.52917721067d0
                                call calcpot(nmono,v2,coord(:))
                                v2 = v2/627.509474d0
                                coord(:) = coord(:)/0.52917721067d0
                                test(1,i) = v2
                                coord(l) = coord(l) - 2*dx
                            else if (j.eq.2) then
                                coord(l) = coord(l) + dx
                                coord(:) = coord(:)*0.52917721067d0
                                call calcpot(nmono,v2,coord(:))
                                v2 = v2/627.509474d0
                                coord(:) = coord(:)/0.52917721067d0
                                test(2,i) = v2
                                coord(l) = coord(l) - dx
                            else if (j.eq.3) then
                                coord(l) = coord(l) 
                                coord(:) = coord(:)*0.52917721067d0
                                call calcpot(nmono,v2,coord(:))
                                v2 = v2/627.509474d0
                                coord(:) = coord(:)/0.52917721067d0
                                test(3,i) = v2
                                coord(l) = coord(l)
                            else if (j.eq.4) then  
                                coord(l) = coord(l) - dx
                                coord(:) = coord(:)*0.52917721067d0
                                call calcpot(nmono,v2,coord(:))
                                v2 = v2/627.509474d0
                                coord(:) = coord(:)/0.52917721067d0
                                test(4,i) = v2
                                coord(l) = coord(l) + dx
                            else if (j.eq.5) then  
                                coord(l) = coord(l) - 2*dx
                                coord(:) = coord(:)*0.52917721067d0
                                call calcpot(nmono,v2,coord(:))
                                v2 = v2/627.509474d0
                                coord(:) = coord(:)/0.52917721067d0
                                test(5,i) = v2 
                                coord(l) = coord(l) + 2*dx
                            endif
                        enddo
                    if (i.eq.1) then
                        coord(k) = coord(k) + 2*dx
                    else if (i.eq.2) then
                        coord(k) = coord(k) + dx
                    else if (i.eq.3) then
                        coord(k) = coord(k) 
                    else if (i.eq.4) then  
                        coord(k) = coord(k) - dx
                    else if (i.eq.5) then  
                        coord(k) = coord(k) - 2*dx
                    endif
                enddo
                if (k.eq.l) then
                    force(k,l) = (-test(3,5) + (16*test(3,4))
     1              -(30*(test(3,3))) + (16*(test(3,2))) -test(3,1))/
     1              (12*dx**2)
                else 
                    force(k,l) = ((8*(test(4,5) + test(5,4) + test(1,2) 
     1              +test(2,1))) - (8*(test(2,5) + test(1,4) + test(4,1)
     1              + test(5,2))) - (test(5,5) + test(1,1) - test(5,1) -
     1              test(1,5)) + (64*(test(2,4) + test(4,2) - test(4,4)
     1              - test(2,2))))/(144.d0*dx**2)
                endif
                do i = 1,5
                    do j = 1,5
                        test(i,j) = 0.d0
                    enddo
                enddo
            enddo
        enddo
        do k =1,ip
            do l = 1,ip
                force(l,k) = force(k,l)
            enddo
        enddo
        do i = 1,ip
            if (i.lt.4) then
                weight(i) = 15.99491461957d0
            else if ((i.gt.9).and.(i.lt.13)) then
                weight(i) = 15.99491461957d0
            else if ((i.gt.18).and.(i.lt.22)) then
                weight(i) = 15.99491461957d0
            else if ((i.gt.27).and.(i.lt.31)) then
                weight(i) = 15.99491461957d0
            else if ((i.gt.36).and.(i.lt.40)) then
                weight(i) = 15.99491461957d0
            else if ((i.gt.45).and.(i.lt.49)) then
                weight(i) = 15.99491461957d0
            else 
                weight(i) = 1.00782503223d0
            endif
            weight(i) = weight(i)*1822.88852962d0
        enddo
        do i = 1,ip
            do j = 1,ip
                force(i,j) =force(i,j)/(sqrt(weight(i))*sqrt(weight(j)))
                eigenvector(i,j) = force(i,j)
            enddo
        enddo
        call house(eigenvector,ip,ip,eigenvalue,scratch)
        do i = 1,ip
            if (eigenvalue(i).gt.0.d0) then
                freq(i) = sqrt(eigenvalue(i))
                freq(i) = freq(i)/((2*pi)*137.035999139
     1          *5.2917721067e-9)
            else
                eigenvalue(i) = -eigenvalue(i)
                freq(i) = sqrt(eigenvalue(i))
                freq(i) = -freq(i)/((2*pi)*137.035999139*
     1          5.2917721067e-9)
            endif
        enddo
        it = 0
        tot_zpe = 0.d0
        do i = 1,ip
            if (i.gt.6) then
                it = it + 1
                eigenvalue_keep(it) = eigenvalue(i)
                evec_keep(:,it) = eigenvector(:,i)
                freq_keep(it) = freq(i)/219474.6
                tot_zpe = tot_zpe + freq_keep(it)
            endif
        enddo
        return
        end subroutine

