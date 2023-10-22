module search
    use interactions
    use domain
    use functions
    use initialize
    use kernel


    implicit none
    
    contains

    ! Searching neighbourhood for interacting  particles
    subroutine nsearch(list,n,cell1,cell2,cell3,cell4,smln,ms)
        implicit none
        ! integer, parameter :: dpint = selected_real_kind(15, 100) 
        type(verlet),intent(inout) :: list
        type(cell),intent(inout) :: cell1,cell2,cell3,cell4
        type(cell) :: blockset(4)
        real(dp),intent(in) :: smln
        integer,intent(in) :: n
        integer ::counter1,counter2,intct
        real(dp) :: radial
        integer,intent(in) :: ms

        blockset(1)=cell1
        blockset(2)=cell2
        blockset(3)=cell3
        blockset(4)=cell4

        intct=0
        list%count=0
        cell1%plist(n)%mirror=.false.
        if((ms==0)) then
            do counter1=1,4
            if (blockset(counter1)%ptot/=0) then
                do counter2=1,blockset(counter1)%ptot
                radial=dist(blockset(1)%plist(n),blockset(counter1)%plist(counter2))
                if((radial<(2*smln)))then
                    if ((cell1%plist(n)%tid/=blockset(counter1)%plist(counter2)%tid) &
                    .or.((cell1%plist(n)%tid==blockset(counter1)%plist(counter2)%tid) &
                    .and.(cell1%plist(n)%pid/=blockset(counter1)%plist(counter2)%pid))) then

                    list%count=list%count+1
                    list%interlist(1,list%count)=blockset(counter1)%cellid(1)
                    list%interlist(2,list%count)=blockset(counter1)%cellid(2)
                    list%interlist(3,list%count)=counter2
                    list%dist(list%count)=radial
                    if (blockset(counter1)%plist(counter2)%tid<=2) then
                        blockset(1)%plist(n)%near_boun=.true.
                    end if
                    end if

                end if
                end do
            end if
            end do

            ! if (intct>0) then
            ! cell1%plist(n)%mirror=.true.
            ! end if

        else

            do counter1=1,1
                if (blockset(counter1)%ptot/=0) then
                    do counter2=1,blockset(counter1)%ptot
                    radial=dist(blockset(1)%plist(n),blockset(counter1)%plist(counter2))
                    if((radial<(2*smln)))then
                        if ((cell1%plist(n)%tid/=blockset(counter1)%plist(counter2)%tid) &
                        .or.((cell1%plist(n)%tid==blockset(counter1)%plist(counter2)%tid) &
                        .and.(cell1%plist(n)%pid/=blockset(counter1)%plist(counter2)%pid))) then
    
                        list%count=list%count+1
                        list%interlist(1,list%count)=blockset(counter1)%cellid(1)
                        list%interlist(2,list%count)=blockset(counter1)%cellid(2)
                        list%interlist(3,list%count)=counter2
                        list%dist(list%count)=radial
                        if (blockset(counter1)%plist(counter2)%tid<=2) then
                            blockset(1)%plist(n)%near_boun=.true.
                        end if
                        end if
    
                    end if
                    end do
                end if
            end do

            ! if (intct>0) then
            ! cell1%plist(n)%mirror=.true.
            ! end if

        end if

    end subroutine nsearch

    ! Searching neighbourhood for ghost particles
    subroutine nsearchghost(list,n,cell1,cell2,cell3,cell4 &
        ,cell5,cell6,cell7,cell8,cell9,smln) 
        implicit none
        type(verlet),intent(inout) :: list
        type(cell),intent(in) :: cell1,cell2,cell3,cell4,cell5,cell6,cell7,cell8,cell9
        type(cell) :: blockset(9)
        real(dp),intent(in) :: smln
        integer,intent(in) :: n
        integer::counter1,counter2
        real(dp) :: radial
    

        blockset(1)=cell1
        blockset(2)=cell2
        blockset(3)=cell3
        blockset(4)=cell4
        blockset(5)=cell5
        blockset(6)=cell6
        blockset(7)=cell7
        blockset(8)=cell8
        blockset(9)=cell9


        list%count=0

            do counter1=1,9
            if (blockset(counter1)%ptot/=0) then
                do counter2=1,blockset(counter1)%ptot
                radial=bpdist(blockset(1)%plist(n),blockset(counter1)%plist(counter2))
                if((radial<(2*blen*smln)).and.((blockset(counter1)%plist(counter2)%tid==3)))then
                    
                    list%count=list%count+1
                    list%interlist(1,list%count)=blockset(counter1)%cellid(1)
                    list%interlist(2,list%count)=blockset(counter1)%cellid(2)
                    list%interlist(3,list%count)=counter2
                    list%dist(list%count)=radial
                end if
                end do
            end if
            end do



    end subroutine nsearchghost

    subroutine neighbour()

        implicit none

        !$omp do private(k,i,j) schedule(runtime) collapse(2)
        do j=sx,ex
            do i=sy,ey
            if (dpcell(i,j)%ptot/=0) then
    
                do k=1,dpcell(i,j)%ptot

                ! if (dpcell(i,j)%plist(k)%tid/=4) then

                if(dpcell(i,j)%plist(k)%x<=(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y<=(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    call nsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i,j-1), &
                    dpcell(i+1,j),dpcell(i+1,j-1),h1,0)

                else if(dpcell(i,j)%plist(k)%x>(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y<=(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    call nsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i,j+1), &
                    dpcell(i+1,j),dpcell(i+1,j+1),h1,0)

                else if(dpcell(i,j)%plist(k)%x>(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y>(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    call nsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i-1,j+1), &
                    dpcell(i-1,j),dpcell(i,j+1),h1,0)

                else if(dpcell(i,j)%plist(k)%x<=(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y>(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    call nsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i-1,j), &
                    dpcell(i-1,j-1),dpcell(i,j-1),h1,0)

                else if(dpcell(i,j)%plist(k)%x==(dpcell(i,j)%xleft+ &
                    abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                    dpcell(i,j)%plist(k)%y==(dpcell(i,j)%ybot+ &
                    abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    call nsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i-1,j), &
                    dpcell(i-1,j-1),dpcell(i,j-1),h1,1)
                end if

                ! if (dpcell(i,j)%plist(k)%tid<=2) then

                

                ! call nsearchghost(dpcell(i,j)%list2(k),k,dpcell(i,j),dpcell(i,j-1), &
                ! dpcell(i,j+1),dpcell(i+1,j),dpcell(i-1,j),dpcell(i+1,j+1),dpcell(i+1,j-1)&
                ! ,dpcell(i-1,j-1),dpcell(i-1,j+1),h1)

                ! end if


                ! end if


                end do
            end if
            end do               
        end do
        !$omp end do

        return
    end subroutine neighbour

    subroutine neighbour2()

        implicit none

        type(cell) :: blockset(4)
        integer :: counter1,counter2
        real(dp) :: rad,res(4),res1,res2,res3,res4

        !$omp do private(k,i,j,rad,counter1,counter2,blockset,res,res1,res2,res3,res4) &
        !$omp schedule(runtime) collapse(2)
        do j=sx,ex
            do i=sy,ey
            if (dpcell(i,j)%ptot/=0) then
    
                do k=1,dpcell(i,j)%ptot

                res=0.0_dp
                res2=0.0_dp
                res1=0.0_dp
                res3=0.0_dp
                res4=0.0_dp

                if(dpcell(i,j)%plist(k)%x<=(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y<=(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    ! call nsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i,j-1), &
                    ! dpcell(i+1,j),dpcell(i+1,j-1),h1,0)

                    blockset(1)=dpcell(i,j)
                    blockset(2)=dpcell(i+1,j)
                    blockset(3)=dpcell(i,j-1)
                    blockset(4)=dpcell(i+1,j-1)

                    dpcell(i,j)%list(k)%count=0
                    dpcell(i,j)%plist(k)%free=.false.
                    dpcell(i,j)%pplist(k)%gradvx=0.0_dp

                    do counter1=1,4
                        ! if (blockset(counter1)%ptot/=0) then
                        do counter2=1,blockset(counter1)%ptot
                        rad=dist(dpcell(i,j)%plist(k),blockset(counter1)%plist(counter2))
                        if((rad<(2*h1)))then
                            if ((dpcell(i,j)%plist(k)%tid/=blockset(counter1)%plist(counter2)%tid) &
                            .or.((dpcell(i,j)%plist(k)%tid==blockset(counter1)%plist(counter2)%tid) &
                        .and.(dpcell(i,j)%plist(k)%pid/=blockset(counter1)%plist(counter2)%pid))) then

                            dpcell(i,j)%list(k)%count=dpcell(i,j)%list(k)%count+1
                            dpcell(i,j)%list(k)%interlist(1,dpcell(i,j)%list(k)%count)= &
                            blockset(counter1)%cellid(1)
                            dpcell(i,j)%list(k)%interlist(2,dpcell(i,j)%list(k)%count)=&
                            blockset(counter1)%cellid(2)
                            dpcell(i,j)%list(k)%interlist(3,dpcell(i,j)%list(k)%count)=counter2
                            dpcell(i,j)%list(k)%dist(dpcell(i,j)%list(k)%count)=rad

                            if (dpcell(i,j)%plist(k)%tid==3) then

                            res2=Wabx(blockset(counter1)%plist(counter2),dpcell(i,j)%plist(k),rad,h1)
                            res4=Waby(blockset(counter1)%plist(counter2),dpcell(i,j)%plist(k),rad,h1)

                            res1=(blockset(counter1)%plist(counter2)%mass)*&
                            (blockset(counter1)%plist(counter2)%x &
                            -dpcell(i,j)%plist(k)%x)/(blockset(counter1)%plist(counter2)%density)
                            res3=(blockset(counter1)%plist(counter2)%mass)*&
                            (blockset(counter1)%plist(counter2)%y &
                            -dpcell(i,j)%plist(k)%y)/(blockset(counter1)%plist(counter2)%density)

                            res(1)=res(1)+res2*res1
                            res(2)=res(2)+res4*res1
                            res(3)=res(3)+res2*res3
                            res(4)=res(4)+res4*res3

                            dpcell(i,j)%pplist(k)%gradvx=dpcell(i,j)%pplist(k)%gradvx+ &
                            ((blockset(counter1)%plist(counter2)%mass* &
                            (blockset(counter1)%plist(counter2)%x- &
                            dpcell(i,j)%plist(k)%x)*res2)/&
                            blockset(counter1)%plist(counter2)%density) + &
                            ((blockset(counter1)%plist(counter2)%mass* &
                            (blockset(counter1)%plist(counter2)%y- &
                            dpcell(i,j)%plist(k)%y)*res4)/&
                            blockset(counter1)%plist(counter2)%density)

                            end if

                            end if

                        end if
                        end do
                        ! end if
                    end do

                    if (dpcell(i,j)%plist(k)%tid==3) then
                    call invertmat2D(res)
                    if (dpcell(i,j)%pplist(k)%gradvx>=0.80_dp) then
                    dpcell(i,j)%pplist(k)%coff=res
                    else 
                    dpcell(i,j)%pplist(k)%coff=(/1.0_dp,0.0_dp,0.0_dp,1.0_dp/)
                    end if

                    if (dpcell(i,j)%pplist(k)%gradvx<=(lamfs*maxdivr)) then
                    dpcell(i,j)%plist(k)%free=.true.

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        ! if ((dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then
    
                        dpcell(y,x)%plist(pp)%vicinity=.true.

                        ! end if
    
                        end associate
                    end do

                    end if

                    end if



                else if(dpcell(i,j)%plist(k)%x>(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y<=(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    ! call nsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i,j+1), &
                    ! dpcell(i+1,j),dpcell(i+1,j+1),h1,0)

                    blockset(1)=dpcell(i,j)
                    blockset(2)=dpcell(i+1,j)
                    blockset(3)=dpcell(i,j+1)
                    blockset(4)=dpcell(i+1,j+1)

                    dpcell(i,j)%list(k)%count=0
                    dpcell(i,j)%plist(k)%free=.false.
                    dpcell(i,j)%pplist(k)%gradvx=0.0_dp

                    do counter1=1,4
                        ! if (blockset(counter1)%ptot/=0) then
                        do counter2=1,blockset(counter1)%ptot
                        rad=dist(dpcell(i,j)%plist(k),blockset(counter1)%plist(counter2))
                        if((rad<(2*h1)))then
                            if ((dpcell(i,j)%plist(k)%tid/=blockset(counter1)%plist(counter2)%tid) &
                            .or.((dpcell(i,j)%plist(k)%tid==blockset(counter1)%plist(counter2)%tid) &
                        .and.(dpcell(i,j)%plist(k)%pid/=blockset(counter1)%plist(counter2)%pid))) then

                            dpcell(i,j)%list(k)%count=dpcell(i,j)%list(k)%count+1
                            dpcell(i,j)%list(k)%interlist(1,dpcell(i,j)%list(k)%count)= &
                            blockset(counter1)%cellid(1)
                            dpcell(i,j)%list(k)%interlist(2,dpcell(i,j)%list(k)%count)=&
                            blockset(counter1)%cellid(2)
                            dpcell(i,j)%list(k)%interlist(3,dpcell(i,j)%list(k)%count)=counter2
                            dpcell(i,j)%list(k)%dist(dpcell(i,j)%list(k)%count)=rad

                            if (dpcell(i,j)%plist(k)%tid==3) then

                            res2=Wabx(blockset(counter1)%plist(counter2),dpcell(i,j)%plist(k),rad,h1)
                            res4=Waby(blockset(counter1)%plist(counter2),dpcell(i,j)%plist(k),rad,h1)

                            res1=(blockset(counter1)%plist(counter2)%mass)*&
                            (blockset(counter1)%plist(counter2)%x &
                            -dpcell(i,j)%plist(k)%x)/(blockset(counter1)%plist(counter2)%density)
                            res3=(blockset(counter1)%plist(counter2)%mass)*&
                            (blockset(counter1)%plist(counter2)%y &
                            -dpcell(i,j)%plist(k)%y)/(blockset(counter1)%plist(counter2)%density)

                            res(1)=res(1)+res2*res1
                            res(2)=res(2)+res4*res1
                            res(3)=res(3)+res2*res3
                            res(4)=res(4)+res4*res3

                            dpcell(i,j)%pplist(k)%gradvx=dpcell(i,j)%pplist(k)%gradvx+ &
                            ((blockset(counter1)%plist(counter2)%mass* &
                            (blockset(counter1)%plist(counter2)%x- &
                            dpcell(i,j)%plist(k)%x)*res2)/&
                            blockset(counter1)%plist(counter2)%density) + &
                            ((blockset(counter1)%plist(counter2)%mass* &
                            (blockset(counter1)%plist(counter2)%y- &
                            dpcell(i,j)%plist(k)%y)*res4)/&
                            blockset(counter1)%plist(counter2)%density)

                            end if

                            end if

                        end if
                        end do
                        ! end if
                    end do

                    if (dpcell(i,j)%plist(k)%tid==3) then
                    call invertmat2D(res)
                    if (dpcell(i,j)%pplist(k)%gradvx>=0.80_dp) then
                    dpcell(i,j)%pplist(k)%coff=res
                    else 
                    dpcell(i,j)%pplist(k)%coff=(/1.0_dp,0.0_dp,0.0_dp,1.0_dp/)
                    end if

                    if (dpcell(i,j)%pplist(k)%gradvx<=(lamfs*maxdivr)) then
                    dpcell(i,j)%plist(k)%free=.true.

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        ! if ((dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then
    
                        dpcell(y,x)%plist(pp)%vicinity=.true.

                        ! end if
    
                        end associate
                    end do

                    end if

                    end if

                else if(dpcell(i,j)%plist(k)%x>(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y>(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    ! call nsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i-1,j+1), &
                    ! dpcell(i-1,j),dpcell(i,j+1),h1,0)

                    blockset(1)=dpcell(i,j)
                    blockset(2)=dpcell(i-1,j)
                    blockset(3)=dpcell(i-1,j+1)
                    blockset(4)=dpcell(i,j+1)

                    dpcell(i,j)%list(k)%count=0
                    dpcell(i,j)%plist(k)%free=.false.
                    dpcell(i,j)%pplist(k)%gradvx=0.0_dp

                    do counter1=1,4
                        ! if (blockset(counter1)%ptot/=0) then
                        do counter2=1,blockset(counter1)%ptot
                        rad=dist(dpcell(i,j)%plist(k),blockset(counter1)%plist(counter2))
                        if((rad<(2*h1)))then
                            if ((dpcell(i,j)%plist(k)%tid/=blockset(counter1)%plist(counter2)%tid) &
                            .or.((dpcell(i,j)%plist(k)%tid==blockset(counter1)%plist(counter2)%tid) &
                        .and.(dpcell(i,j)%plist(k)%pid/=blockset(counter1)%plist(counter2)%pid))) then

                            dpcell(i,j)%list(k)%count=dpcell(i,j)%list(k)%count+1
                            dpcell(i,j)%list(k)%interlist(1,dpcell(i,j)%list(k)%count)= &
                            blockset(counter1)%cellid(1)
                            dpcell(i,j)%list(k)%interlist(2,dpcell(i,j)%list(k)%count)=&
                            blockset(counter1)%cellid(2)
                            dpcell(i,j)%list(k)%interlist(3,dpcell(i,j)%list(k)%count)=counter2
                            dpcell(i,j)%list(k)%dist(dpcell(i,j)%list(k)%count)=rad

                            if (dpcell(i,j)%plist(k)%tid==3) then

                            res2=Wabx(blockset(counter1)%plist(counter2),dpcell(i,j)%plist(k),rad,h1)
                            res4=Waby(blockset(counter1)%plist(counter2),dpcell(i,j)%plist(k),rad,h1)

                            res1=(blockset(counter1)%plist(counter2)%mass)*&
                            (blockset(counter1)%plist(counter2)%x &
                            -dpcell(i,j)%plist(k)%x)/(blockset(counter1)%plist(counter2)%density)
                            res3=(blockset(counter1)%plist(counter2)%mass)*&
                            (blockset(counter1)%plist(counter2)%y &
                            -dpcell(i,j)%plist(k)%y)/(blockset(counter1)%plist(counter2)%density)

                            res(1)=res(1)+res2*res1
                            res(2)=res(2)+res4*res1
                            res(3)=res(3)+res2*res3
                            res(4)=res(4)+res4*res3

                            dpcell(i,j)%pplist(k)%gradvx=dpcell(i,j)%pplist(k)%gradvx+ &
                            ((blockset(counter1)%plist(counter2)%mass* &
                            (blockset(counter1)%plist(counter2)%x- &
                            dpcell(i,j)%plist(k)%x)*res2)/&
                            blockset(counter1)%plist(counter2)%density) + &
                            ((blockset(counter1)%plist(counter2)%mass* &
                            (blockset(counter1)%plist(counter2)%y- &
                            dpcell(i,j)%plist(k)%y)*res4)/&
                            blockset(counter1)%plist(counter2)%density)

                            end if

                            end if

                        end if
                        end do
                        ! end if
                    end do

                    if (dpcell(i,j)%plist(k)%tid==3) then
                    call invertmat2D(res)
                    if (dpcell(i,j)%pplist(k)%gradvx>=0.80_dp) then
                    dpcell(i,j)%pplist(k)%coff=res
                    else 
                    dpcell(i,j)%pplist(k)%coff=(/1.0_dp,0.0_dp,0.0_dp,1.0_dp/)
                    end if

                    if (dpcell(i,j)%pplist(k)%gradvx<=(lamfs*maxdivr)) then
                    dpcell(i,j)%plist(k)%free=.true.

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        ! if ((dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then
    
                        dpcell(y,x)%plist(pp)%vicinity=.true.

                        ! end if
    
                        end associate
                    end do

                    end if

                    end if

                else if(dpcell(i,j)%plist(k)%x<=(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y>(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    ! call nsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i-1,j), &
                    ! dpcell(i-1,j-1),dpcell(i,j-1),h1,0)

                    blockset(1)=dpcell(i,j)
                    blockset(2)=dpcell(i-1,j)
                    blockset(3)=dpcell(i-1,j-1)
                    blockset(4)=dpcell(i,j-1)

                    dpcell(i,j)%list(k)%count=0
                    dpcell(i,j)%plist(k)%free=.false.
                    dpcell(i,j)%pplist(k)%gradvx=0.0_dp

                    do counter1=1,4
                        ! if (blockset(counter1)%ptot/=0) then
                        do counter2=1,blockset(counter1)%ptot
                        rad=dist(dpcell(i,j)%plist(k),blockset(counter1)%plist(counter2))
                        if((rad<(2*h1)))then
                            if ((dpcell(i,j)%plist(k)%tid/=blockset(counter1)%plist(counter2)%tid) &
                            .or.((dpcell(i,j)%plist(k)%tid==blockset(counter1)%plist(counter2)%tid) &
                        .and.(dpcell(i,j)%plist(k)%pid/=blockset(counter1)%plist(counter2)%pid))) then

                            dpcell(i,j)%list(k)%count=dpcell(i,j)%list(k)%count+1
                            dpcell(i,j)%list(k)%interlist(1,dpcell(i,j)%list(k)%count)= &
                            blockset(counter1)%cellid(1)
                            dpcell(i,j)%list(k)%interlist(2,dpcell(i,j)%list(k)%count)=&
                            blockset(counter1)%cellid(2)
                            dpcell(i,j)%list(k)%interlist(3,dpcell(i,j)%list(k)%count)=counter2
                            dpcell(i,j)%list(k)%dist(dpcell(i,j)%list(k)%count)=rad

                            if (dpcell(i,j)%plist(k)%tid==3) then

                            res2=Wabx(blockset(counter1)%plist(counter2),dpcell(i,j)%plist(k),rad,h1)
                            res4=Waby(blockset(counter1)%plist(counter2),dpcell(i,j)%plist(k),rad,h1)

                            res1=(blockset(counter1)%plist(counter2)%mass)*&
                            (blockset(counter1)%plist(counter2)%x &
                            -dpcell(i,j)%plist(k)%x)/(blockset(counter1)%plist(counter2)%density)
                            res3=(blockset(counter1)%plist(counter2)%mass)*&
                            (blockset(counter1)%plist(counter2)%y &
                            -dpcell(i,j)%plist(k)%y)/(blockset(counter1)%plist(counter2)%density)

                            res(1)=res(1)+res2*res1
                            res(2)=res(2)+res4*res1
                            res(3)=res(3)+res2*res3
                            res(4)=res(4)+res4*res3

                            dpcell(i,j)%pplist(k)%gradvx=dpcell(i,j)%pplist(k)%gradvx+ &
                            ((blockset(counter1)%plist(counter2)%mass* &
                            (blockset(counter1)%plist(counter2)%x- &
                            dpcell(i,j)%plist(k)%x)*res2)/&
                            blockset(counter1)%plist(counter2)%density) + &
                            ((blockset(counter1)%plist(counter2)%mass* &
                            (blockset(counter1)%plist(counter2)%y- &
                            dpcell(i,j)%plist(k)%y)*res4)/&
                            blockset(counter1)%plist(counter2)%density)

                            end if

                            end if

                        end if
                        end do
                        ! end if
                    end do

                    if (dpcell(i,j)%plist(k)%tid==3) then
                    call invertmat2D(res)
                    if (dpcell(i,j)%pplist(k)%gradvx>=0.80_dp) then
                    dpcell(i,j)%pplist(k)%coff=res
                    else 
                    dpcell(i,j)%pplist(k)%coff=(/1.0_dp,0.0_dp,0.0_dp,1.0_dp/)
                    end if

                    if (dpcell(i,j)%pplist(k)%gradvx<=(lamfs*maxdivr)) then
                    dpcell(i,j)%plist(k)%free=.true.

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        ! if ((dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then
    
                        dpcell(y,x)%plist(pp)%vicinity=.true.

                        ! end if
    
                        end associate
                    end do

                    end if

                    end if

                else if(dpcell(i,j)%plist(k)%x==(dpcell(i,j)%xleft+ &
                    abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                    dpcell(i,j)%plist(k)%y==(dpcell(i,j)%ybot+ &
                    abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    ! call nsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i-1,j), &
                    ! dpcell(i-1,j-1),dpcell(i,j-1),h1,1)

                    blockset(1)=dpcell(i,j)
                    blockset(2)=dpcell(i+1,j)
                    blockset(3)=dpcell(i,j-1)
                    blockset(4)=dpcell(i+1,j-1)

                    dpcell(i,j)%list(k)%count=0
                    dpcell(i,j)%plist(k)%free=.false.
                    dpcell(i,j)%pplist(k)%gradvx=0.0_dp

                    do counter1=1,1
                        ! if (blockset(counter1)%ptot/=0) then
                        do counter2=1,blockset(counter1)%ptot
                        rad=dist(dpcell(i,j)%plist(k),blockset(counter1)%plist(counter2))
                        if((rad<(2*h1)))then
                            if ((dpcell(i,j)%plist(k)%tid/=blockset(counter1)%plist(counter2)%tid) &
                            .or.((dpcell(i,j)%plist(k)%tid==blockset(counter1)%plist(counter2)%tid) &
                        .and.(dpcell(i,j)%plist(k)%pid/=blockset(counter1)%plist(counter2)%pid))) then

                            dpcell(i,j)%list(k)%count=dpcell(i,j)%list(k)%count+1
                            dpcell(i,j)%list(k)%interlist(1,dpcell(i,j)%list(k)%count)= &
                            blockset(counter1)%cellid(1)
                            dpcell(i,j)%list(k)%interlist(2,dpcell(i,j)%list(k)%count)=&
                            blockset(counter1)%cellid(2)
                            dpcell(i,j)%list(k)%interlist(3,dpcell(i,j)%list(k)%count)=counter2
                            dpcell(i,j)%list(k)%dist(dpcell(i,j)%list(k)%count)=rad

                            if (dpcell(i,j)%plist(k)%tid==3) then

                            res2=Wabx(blockset(counter1)%plist(counter2),dpcell(i,j)%plist(k),rad,h1)
                            res4=Waby(blockset(counter1)%plist(counter2),dpcell(i,j)%plist(k),rad,h1)

                            res1=(blockset(counter1)%plist(counter2)%mass)*&
                            (blockset(counter1)%plist(counter2)%x &
                            -dpcell(i,j)%plist(k)%x)/(blockset(counter1)%plist(counter2)%density)
                            res3=(blockset(counter1)%plist(counter2)%mass)*&
                            (blockset(counter1)%plist(counter2)%y &
                            -dpcell(i,j)%plist(k)%y)/(blockset(counter1)%plist(counter2)%density)

                            res(1)=res(1)+res2*res1
                            res(2)=res(2)+res4*res1
                            res(3)=res(3)+res2*res3
                            res(4)=res(4)+res4*res3

                            dpcell(i,j)%pplist(k)%gradvx=dpcell(i,j)%pplist(k)%gradvx+ &
                            ((blockset(counter1)%plist(counter2)%mass* &
                            (blockset(counter1)%plist(counter2)%x- &
                            dpcell(i,j)%plist(k)%x)*res2)/&
                            blockset(counter1)%plist(counter2)%density) + &
                            ((blockset(counter1)%plist(counter2)%mass* &
                            (blockset(counter1)%plist(counter2)%y- &
                            dpcell(i,j)%plist(k)%y)*res4)/&
                            blockset(counter1)%plist(counter2)%density)

                            end if

                            end if

                        end if
                        end do
                        ! end if
                    end do

                    if (dpcell(i,j)%plist(k)%tid==3) then
                    call invertmat2D(res)
                    if (dpcell(i,j)%pplist(k)%gradvx>=0.80_dp) then
                    dpcell(i,j)%pplist(k)%coff=res
                    else 
                    dpcell(i,j)%pplist(k)%coff=(/1.0_dp,0.0_dp,0.0_dp,1.0_dp/)
                    end if

                    if (dpcell(i,j)%pplist(k)%gradvx<=(lamfs*maxdivr)) then
                    dpcell(i,j)%plist(k)%free=.true.

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        ! if ((dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then
    
                        dpcell(y,x)%plist(pp)%vicinity=.true.

                        ! end if
    
                        end associate
                    end do

                    end if

                    end if


                end if

                end do
            end if
            end do               
        end do
        !$omp end do

        return
    end subroutine 
    
end module search