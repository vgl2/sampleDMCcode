        subroutine calc_b_mat(coord,q_f_r,q_f_cart,nmono,first_trans,
     1  second_trans)
        implicit real*8(a-h,o-z)
        parameter (nvib=18)
        parameter (nvibw=3)
        parameter (ndim=54)
        parameter (ndimw=9)
        dimension r(nvibw,nmono),coord(ndim),test_r1(5,nmono),
     1  test_r2(5,nmono),test_r3(5,nmono),first(ndimw,nvibw,nmono),
     1  first_trans(nvib,ndim),q_f_r(nvibw,nmono),q_f_cart(ndim),
     1  coord_mono(ndimw,nmono),first_mat(ndim,nvib),q_f_rfull(nvib),
     1  second(ndimw,nvibw,nmono),second_mat(ndim,nvib),
     1  second_trans(nvib,ndim)
        pi = dacos(-1.d0)
        ip = 0
        dx = 0.0001
C       Calculate internals
        ip = 0
        do i = 1,nmono
            do j = 1,ndimw
                ip = ip + 1
                coord_mono(j,i) = coord(ip)
            enddo
        enddo
        do k = 1,nmono
            r(1,k)=sqrt((coord_mono(4,k)-coord_mono(1,k))**2+(
     1      coord_mono(5,k)-coord_mono(2,k))**2+(coord_mono(6,k)-
     1      coord_mono(3,k))**2)
            r(2,k)=sqrt((coord_mono(7,k)-coord_mono(1,k))**2+(
     1      coord_mono(8,k)-coord_mono(2,k))**2+(coord_mono(9,k)-
     1      coord_mono(3,k))**2)
            r(3,k)=dacos((((coord_mono(4,k)-coord_mono(1,k))*(
     1      coord_mono(7,k)-coord_mono(1,k)))+((coord_mono(5,k)-
     1      coord_mono(2,k))*(coord_mono(8,k)-coord_mono(2,k)))+((
     1      coord_mono(6,k)-coord_mono(3,k))*(coord_mono(9,k)-
     1      coord_mono(3,k))))/(r(1,k)*r(2,k)))
        enddo
        do i = 1,nmono
            do j = 1,5
                test_r1(j,i) = 0.d0
                test_r2(j,i) = 0.d0
                test_r3(j,i) = 0.d0
            enddo
        enddo
        do j = 1,nmono
c       Do 5 point finite differnce 
            do k = 1,ndimw
                do i = 1,5
                    if (i.eq.1) then
                        coord_mono(k,j) = coord_mono(k,j) + 2*dx
                        r(1,j)=sqrt((coord_mono(4,j)-coord_mono(1,j))**2
     1                  +(coord_mono(5,j)-coord_mono(2,j))**2+(
     1                  coord_mono(6,j)-coord_mono(3,j))**2)
                        r(2,j)=sqrt((coord_mono(7,j)-coord_mono(1,j))**2
     1                  +(coord_mono(8,j)-coord_mono(2,j))**2+(
     1                  coord_mono(9,j)-coord_mono(3,j))**2)
                        r(3,j)=dacos((((coord_mono(4,j)-coord_mono(1,j))
     1                  *(coord_mono(7,j)-coord_mono(1,j)))+((
     1                  coord_mono(5,j)-coord_mono(2,j))*(
     1                  coord_mono(8,j)-coord_mono(2,j)))+((
     1                  coord_mono(6,j)-coord_mono(3,j))*
     1                  (coord_mono(9,j)-coord_mono(3,j))))/(r(1,j)*
     1                  r(2,j)))
                        test_r1(i,j) = r(1,j)
                        test_r2(i,j) = r(2,j)
                        test_r3(i,j) = r(3,j)
                        coord_mono(k,j) = coord_mono(k,j) - 2*dx
                    endif
                    if (i.eq.2) then
                        coord_mono(k,j) = coord_mono(k,j) + dx
                        r(1,j)=sqrt((coord_mono(4,j)-coord_mono(1,j))**2
     1                  +(coord_mono(5,j)-coord_mono(2,j))**2+(
     1                  coord_mono(6,j)-coord_mono(3,j))**2)
                        r(2,j)=sqrt((coord_mono(7,j)-coord_mono(1,j))**2
     1                  +(coord_mono(8,j)-coord_mono(2,j))**2+(
     1                  coord_mono(9,j)-coord_mono(3,j))**2)
                        r(3,j)=dacos((((coord_mono(4,j)-coord_mono(1,j))
     1                  *(coord_mono(7,j)-coord_mono(1,j)))+((
     1                  coord_mono(5,j)-coord_mono(2,j))*(
     1                  coord_mono(8,j)-coord_mono(2,j)))+((
     1                  coord_mono(6,j)-coord_mono(3,j))*
     1                  (coord_mono(9,j)-coord_mono(3,j))))/(r(1,j)*
     1                  r(2,j)))
                        test_r1(i,j) = r(1,j)
                        test_r2(i,j) = r(2,j)
                        test_r3(i,j) = r(3,j)
                        coord_mono(k,j) = coord_mono(k,j) - dx
                    endif
                    if (i.eq.3) then
                        coord_mono(k,j) = coord_mono(k,j) 
                        r(1,j)=sqrt((coord_mono(4,j)-coord_mono(1,j))**2
     1                  +(coord_mono(5,j)-coord_mono(2,j))**2+(
     1                  coord_mono(6,j)-coord_mono(3,j))**2)
                        r(2,j)=sqrt((coord_mono(7,j)-coord_mono(1,j))**2
     1                  +(coord_mono(8,j)-coord_mono(2,j))**2+(
     1                  coord_mono(9,j)-coord_mono(3,j))**2)
                        r(3,j)=dacos((((coord_mono(4,j)-coord_mono(1,j))
     1                  *(coord_mono(7,j)-coord_mono(1,j)))+((
     1                  coord_mono(5,j)-coord_mono(2,j))*(
     1                  coord_mono(8,j)-coord_mono(2,j)))+((
     1                  coord_mono(6,j)-coord_mono(3,j))*
     1                  (coord_mono(9,j)-coord_mono(3,j))))/(r(1,j)*
     1                  r(2,j)))
                        test_r1(i,j) = r(1,j)
                        test_r2(i,j) = r(2,j)
                        test_r3(i,j) = r(3,j)
                        coord_mono(k,j) = coord_mono(k,j) 
                    endif
                    if (i.eq.4) then
                        coord_mono(k,j) = coord_mono(k,j) - dx
                        r(1,j)=sqrt((coord_mono(4,j)-coord_mono(1,j))**2
     1                  +(coord_mono(5,j)-coord_mono(2,j))**2+(
     1                  coord_mono(6,j)-coord_mono(3,j))**2)
                        r(2,j)=sqrt((coord_mono(7,j)-coord_mono(1,j))**2
     1                  +(coord_mono(8,j)-coord_mono(2,j))**2+(
     1                  coord_mono(9,j)-coord_mono(3,j))**2)
                        r(3,j)=dacos((((coord_mono(4,j)-coord_mono(1,j))
     1                  *(coord_mono(7,j)-coord_mono(1,j)))+((
     1                  coord_mono(5,j)-coord_mono(2,j))*(
     1                  coord_mono(8,j)-coord_mono(2,j)))+((
     1                  coord_mono(6,j)-coord_mono(3,j))*
     1                  (coord_mono(9,j)-coord_mono(3,j))))/(r(1,j)*
     1                  r(2,j)))
                        test_r1(i,j) = r(1,j)
                        test_r2(i,j) = r(2,j)
                        test_r3(i,j) = r(3,j)
                        coord_mono(k,j) = coord_mono(k,j) + dx
                    endif
                    if (i.eq.5) then
                        coord_mono(k,j) = coord_mono(k,j) - 2*dx
                        r(1,j)=sqrt((coord_mono(4,j)-coord_mono(1,j))**2
     1                  +(coord_mono(5,j)-coord_mono(2,j))**2+(
     1                  coord_mono(6,j)-coord_mono(3,j))**2)
                        r(2,j)=sqrt((coord_mono(7,j)-coord_mono(1,j))**2
     1                  +(coord_mono(8,j)-coord_mono(2,j))**2+(
     1                  coord_mono(9,j)-coord_mono(3,j))**2)
                        r(3,j)=dacos((((coord_mono(4,j)-coord_mono(1,j))
     1                  *(coord_mono(7,j)-coord_mono(1,j)))+((
     1                  coord_mono(5,j)-coord_mono(2,j))*(
     1                  coord_mono(8,j)-coord_mono(2,j)))+((
     1                  coord_mono(6,j)-coord_mono(3,j))*
     1                  (coord_mono(9,j)-coord_mono(3,j))))/(r(1,j)*
     1                  r(2,j)))
                        test_r1(i,j) = r(1,j)
                        test_r2(i,j) = r(2,j)
                        test_r3(i,j) = r(3,j)
                        coord_mono(k,j) = coord_mono(k,j) + 2*dx
                    endif
                enddo
                do i = 1,nvib
                    if (i.eq.1) then
                        first(k,i,j)=(((-1.d0/12.d0)*(test_r1(1,j)))+((
     1                  2.d0/3.d0)*(test_r1(2,j)))+((-2.d0/3.d0)*(
     1                  test_r1(4,j)))+((1.d0/12.d0)*test_r1(5,j)))/(dx)
                        second(k,i,j)=(((-1.d0/12.d0)*(test_r1(1,j)))+((
     1                  4.d0/3.d0)*(test_r1(2,j)))+((-5.d0/2.d0)*
     1                  test_r1(3,j))+((4.d0/3.d0)*(test_r1(4,j)))+
     1                  ((-1.d0/12.d0)*test_r1(5,j)))/(dx**2)
                    else if (i.eq.2) then
                        first(k,i,j)=(((-1.d0/12.d0)*(test_r2(1,j)))+((
     1                  2.d0/3.d0)*(test_r2(2,j)))+((-2.d0/3.d0)*(
     1                  test_r2(4,j)))+((1.d0/12.d0)*test_r2(5,j)))/(dx)
                        second(k,i,j)=(((-1.d0/12.d0)*(test_r2(1,j)))+((
     1                  4.d0/3.d0)*(test_r2(2,j)))+((-5.d0/2.d0)*
     1                  test_r2(3,j))+((4.d0/3.d0)*(test_r2(4,j)))+
     1                  ((-1.d0/12.d0)*test_r2(5,j)))/(dx**2)
                    else if (i.eq.3) then
                        first(k,i,j)=(((-1.d0/12.d0)*(test_r3(1,j)))+(
     1                  (2.d0/3.d0)*(test_r3(2,j)))+((-2.d0/3.d0)*(
     1                  test_r3(4,j)))+((1.d0/12.d0)*test_r3(5,j)))/(dx)
                        second(k,i,j)=(((-1.d0/12.d0)*(test_r3(1,j)))+((
     1                  4.d0/3.d0)*(test_r3(2,j)))+((-5.d0/2.d0)*
     1                  test_r3(3,j))+((4.d0/3.d0)*(test_r3(4,j)))+
     1                  ((-1.d0/12.d0)*test_r3(5,j)))/(dx**2)
                    endif
                enddo
                do i = 1,5
                    test_r1(i,j) = 0.d0
                    test_r2(i,j) = 0.d0
                    test_r3(i,j) = 0.d0
                enddo
            enddo
        enddo
        ip = 0
        do j =1 ,nmono  
            do i = 1,nvibw
                ip = ip + 1
                q_f_rfull(ip) = q_f_r(i,j)
            enddo
        enddo
        do i = 1,ndim
            do j = 1,nvib
                first_mat(i,j) = 0.d0
                second_mat(i,j) = 0.d0
            enddo
        enddo
        do j = 1,nmono
            if (j.eq.1) then
                do i = 1,ndimw
                    do k = 1,nvibw
                        first_mat(i,k) = first(i,k,j)
                        second_mat(i,k) = second(i,k,j)
                    enddo
                enddo
            else if (j.eq.2) then
                do i = 1,ndimw
                    do k = 1,nvibw
                        first_mat(i+ndimw,k+nvibw) = first(i,k,j)
                        second_mat(i+ndimw,k+nvibw) = second(i,k,j)
                    enddo
                enddo
            else if (j.eq.3) then
                do i = 1,ndimw
                    do k = 1,nvibw
                        first_mat(i+2*ndimw,k+2*nvibw) = first(i,k,j)
                        second_mat(i+2*ndimw,k+2*nvibw) = second(i,k,j)
                    enddo
                enddo
            else if (j.eq.4) then
                do i = 1,ndimw
                    do k = 1,nvibw
                        first_mat(i+3*ndimw,k+3*nvibw) = first(i,k,j)
                        second_mat(i+3*ndimw,k+3*nvibw) = second(i,k,j)
                    enddo
                enddo
            else if (j.eq.5) then
                do i = 1,ndimw
                    do k = 1,nvibw
                        first_mat(i+4*ndimw,k+4*nvibw) = first(i,k,j)
                        second_mat(i+4*ndimw,k+4*nvibw) = second(i,k,j)
                    enddo
                enddo
            else if (j.eq.6) then
                do i = 1,ndimw
                    do k = 1,nvibw
                        first_mat(i+5*ndimw,k+5*nvibw) = first(i,k,j)
                        second_mat(i+5*ndimw,k+5*nvibw) = second(i,k,j)
                    enddo
                enddo
            endif
        enddo
        first_trans = transpose(first_mat)
        second_trans = transpose(second_mat)
        q_f_cart = matmul(q_f_rfull,first_trans)
        return
        end subroutine
        


