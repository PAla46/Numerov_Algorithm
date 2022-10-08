Program qm1
    IMPLICIT NONE
    Integer,parameter :: dp=selected_real_kind(14,200)

    integer :: mesh,i,icl 
    integer :: nodes,hnodes,ncross,kkk,n_iter 
    real(dp) :: xmax,dx,ddx12,norm,arg 
    real(dp) :: eup,elw,e,m,n,q 
    real(dp),allocatable :: x(:),y(:),p(:),vpot(:),f(:)
    character (len=80) :: fileout

    print '(a,$)','Max value for x (typical value:10)?'
    read(*,*) xmax 
    print '(a,$)','Number of grid points'
    read(*,*) mesh 

    allocate(x(0:mesh),y(0:mesh),p(0:mesh),vpot(0:mesh),f(0:mesh))
    dx=xmax/mesh 
    ddx12=dx*dx/12.0_dp 

    do i=0,mesh
        x(i) = float(i)*dx 
        vpot(i) = (x(i)**4) - (6*x(i)*x(i))
    end do 

    print '(a,$)','output file name ='
    read(*,'(a)') fileout 
    open(7,file=fileout,status='unknown',form='formatted')

    search_loop:do  ! start search loop
            print '(a,$)','nodes (type -1 to stop) >>'
            read(*,*) nodes 
            if (nodes<0) then
                close(7)
                deallocate(f,vpot,p,y,x)
                stop 
            end if 

            eup=maxval(vpot(:))
            elw=minval(vpot(:))
            print*,eup,elw 

            print '(a,$)','Trial energy (0=search with bisection ?)'
            read(*,*) e 
            if (e == 0.0_dp) then
                e=0.5_dp*(elw+eup)
                n_iter=1000
            else
                n_iter=1
            end if 

    iterate: do kkk=1,n_iter  ! start iterate loop
            f(0) = ddx12*(2.0_dp*(vpot(0)-e))
            icl=1
            do i=1,mesh 
                f(i)=ddx12*2.0_dp*(vpot(i)-e)
                if(f(i) == 0.0_dp) f(i)=1.d-20
                if(f(i) /= sign(f(i),f(i-1))) icl=i 
            end do 

            if(icl >= mesh-2) then 
                deallocate(f,vpot,p,y,x)
                stop 'last change of sign too far'
            else if (icl<1) then 
                deallocate(f,vpot,p,y,x)
                stop 'no classical turning point?'
            end if 

            f=1.0_dp - f 
            y=0.0_dp 
            hnodes = nodes/2
            if(2*hnodes==nodes) then
                y(0)=1.0_dp 
                m=0.5_dp*(12.0_dp-10.0_dp*f(0))
                y(1)=m*y(0)/f(1)
            else 
                y(0)=0.0_dp 
                y(1)=dx 
            end if 
            ncross=0 
            do i=1,mesh-1  !start mesh loop 
                n=(12.0_dp-10.0_dp*f(i))*y(i)
                y(i+1) = (n-f(i-1)*y(i-1))/f(i+1)
                if (y(i) /= sign(y(i),y(i+1))) ncross=ncross+1 
            end do 

            print*, kkk,e,ncross,hnodes 
    ! print*,e

            if(n_iter>1) then

                if(ncross>hnodes) then
                    eup=e 
                else 
                    elw=e 
                end if 
                e=0.5_dp*(elw+eup) 
                if(eup-elw<1.d-10) then 
                    print*,'ground state energy is:',e 
                end if 
                if (eup-elw<1.d-10) exit iterate 
            end if 

            end do iterate !end iterate loop 

        norm=0.0_dp 
        p(icl:)=0.0_dp 
        do i=0,icl 
            arg=(e-x(i)**2/2.0_dp)
            if (arg>0.0_dp) then 
                p(i) = 1.0_dp/sqrt(arg)
            else
                p(i)=0.0_dp
            end if 
            norm=norm+2.0_dp*dx*p(i)
        end do 

        norm=norm-dx*p(0)
        p(:icl-1)=p(:icl-1)/norm 
        write(7,'("# x y(x) y(x)^2 classical P(x) V")')
        ! x<0 region:
        do i=mesh,1,-1
            write(7,*) -x(i), (-1)**nodes*y(i),y(i)*y(i),p(i),vpot(i)
        end do 
        !x>0 region
        do i=0,mesh
            write(7,*) x(i),y(i),y(i)*y(i),p(i),vpot(i)
        end do 
        write(7,'(/)')

    end do search_loop  !end search loop

end program qm1 