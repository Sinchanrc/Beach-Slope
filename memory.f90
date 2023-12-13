module memory

    use initialize

    implicit none
    

    contains

    subroutine memalloc()

        implicit none

        integer :: i,j,k,m

        !Allocating arrays in the cells
        !$omp parallel do private(i,j,k) schedule (runtime) collapse(2) default(shared)
            do j=sx,ex
            do i=sy,ey

            if (dpcell(i,j)%btot/=0) then
                allocate(dpcell(i,j)%list(fac*fplistmax),dpcell(i,j)%pplist(fac*fplistmax)) 
                ! allocate(dpcell(i,j)%list2(dpcell(i,j)%btot))
                do k=1,fac*fplistmax
                    allocate(dpcell(i,j)%list(k)%interlist(3,fac*fplistmax), &
                    dpcell(i,j)%list(k)%dist(fac*fplistmax))
                    dpcell(i,j)%list(k)%interlist=0
                    dpcell(i,j)%list(k)%dist=0.0_dp

                    allocate(dpcell(i,j)%pplist(k)%meanstrn(4),dpcell(i,j)%pplist(k)%tau(4),&
                    dpcell(i,j)%pplist(k)%coff(4))

                end do

                ! do k=1,dpcell(i,j)%btot 

                !     allocate(dpcell(i,j)%list2(k)%interlist(3,fplist), &
                !     dpcell(i,j)%list2(k)%dist(fplist))
                !     dpcell(i,j)%list2(k)%interlist=0
                !     dpcell(i,j)%list2(k)%dist=0.0_dp

                ! end do


            else
                allocate(dpcell(i,j)%list(fac*fplistmax),dpcell(i,j)%pplist(fac*fplistmax))  
                do k=1,fplistmax
                    allocate(dpcell(i,j)%list(k)%interlist(3,fplistmax), &
                    dpcell(i,j)%list(k)%dist(fplistmax))
                    dpcell(i,j)%list(k)%interlist=0
                    dpcell(i,j)%list(k)%dist=0.0_dp
                    allocate(dpcell(i,j)%pplist(k)%meanstrn(4),dpcell(i,j)%pplist(k)%tau(4),&
                    dpcell(i,j)%pplist(k)%coff(4))

                end do  

            end if

                allocate(dpcell(i,j)%tn%list(ceiling(0.25*fplistmax)),&
                dpcell(i,j)%ts%list(ceiling(0.25*fplistmax)),&
                dpcell(i,j)%te%list(ceiling(0.25*fplistmax)), &
                dpcell(i,j)%tw%list(ceiling(0.25*fplistmax)),dpcell(i,j)%tne%list(ceiling(0.25*fplistmax)),&
                dpcell(i,j)%tnw%list(ceiling(0.25*fplistmax)), &
                dpcell(i,j)%tse%list(ceiling(0.25*fplistmax)),dpcell(i,j)%tsw%list(ceiling(0.25*fplistmax)))
            end do
            end do
        !$omp end parallel do

        ! allocate(probe(pbno),ploc(2,pbno),probedata(2,160))
        ! do k=1,pbno
        !     allocate(probe(k)%part%interlist(3,fplistmax),probe(k)%part%dist(fplistmax))
        ! end do
        ! probedata(:,:)=0.0_dp

    
    end subroutine memalloc

    subroutine matrixid()
        implicit none

        integer :: i,j,k,mo

        count=0

        do j=sx,ex 
            do i=sy,ey
            if (dpcell(i,j)%ptot/=0) then
                do cout=1,dpcell(i,j)%ptot

                    ! if ((dpcell(i,j)%plist(cout)%tid==3)) then
                    count=count+1
                    dpcell(i,j)%plist(cout)%pid=count
                    dpcell(i,j)%plist(cout)%matid=dpcell(i,j)%plist(cout)%pid

                end do
            end if
            end do
        end do        

        finmax=count+reserve_par

        do i=1,reserve_par
            count=count+1
            reserve%tank(i)=count
        end do

        allocate(frow(finmax+1),fvec(finmax),fsol(finmax),fguess(finmax))
        allocate(fval(finmax*ceiling(fac2*fplistmax)),fcol(finmax*ceiling(fac2*fplistmax)))
        allocate(fmatrix(finmax),perm(finmax),pguess(finmax))

        fguess=1000.0_dp
        fsol=0.0_dp
        pguess=0.0_dp

        do j=sx,ex 
            do i=sy,ey
            if (dpcell(i,j)%ptot/=0) then
                do cout=1,dpcell(i,j)%ptot

                    if ((dpcell(i,j)%plist(cout)%tid/=4)) then

                    fguess(dpcell(i,j)%plist(cout)%matid)=dpcell(i,j)%plist(cout)%pressure

                    end if


                end do
            end if
            end do
        end do

        ! ! Allocating fmatrix

            do i=1,finmax

            allocate(fmatrix(i)%val(ceiling(fac2*fplistmax)))
            allocate(fmatrix(i)%col(ceiling(fac2*fplistmax)))

            end do

            allocate(dpar(128),tmp(finmax*(2*finmax+1)+(finmax*(finmax+9))/2+1))
            fsol=fguess
            
        
    end subroutine matrixid


    
end module memory