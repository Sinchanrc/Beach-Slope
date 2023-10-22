module part_shift

    use particle
    use domain
    use kernel
    use initialize

    implicit none

    contains

    subroutine exp_shift
        implicit none    
        ! Explicit particle shifting
        !$omp parallel do private(i,j,k,m,term1,term2) schedule (runtime) collapse(2) default(shared)
            do j=sx,ex
                do i=sy,ey
                    if(dpcell(i,j)%ptot/=0) then
                    
                    do k=1,dpcell(i,j)%ptot
                        dpcell(i,j)%plist(k)%ys=0.0_dp
                        dpcell(i,j)%plist(k)%xs=0.0_dp
                        dpcell(i,j)%plist(k)%vys=0.0_dp
                        dpcell(i,j)%plist(k)%vxs=0.0_dp
                        if ((dpcell(i,j)%list(k)%count/=0) &
                        ! .and.(dpcell(i,j)%plist(k)%mirror.eqv..true.) &
                        .and.(dpcell(i,j)%pplist(k)%gradvx>=0.70_dp) &
                        .and.(dpcell(i,j)%plist(k)%tid>2)) then
                        dpcell(i,j)%plist(k)%ys=0.0_dp
                        dpcell(i,j)%plist(k)%xs=0.0_dp
                        term1=0.0_dp
                        term2=0.0_dp

                        do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                            y=>dpcell(i,j)%list(k)%interlist(2,m), &
                            pp=>dpcell(i,j)%list(k)%interlist(3,m))

                            if (dpcell(i,j)%plist(k)%tid/=4) then

                            term1=term1+dpcell(i,j)%list(k)%dist(m)/ &
                            (dpcell(i,j)%list(k)%count)
                        
                            term2=term2+(dpcell(i,j)%plist(k)%x- &
                            dpcell(y,x)%plist(pp)%x)/ &
                            dpcell(i,j)%list(k)%dist(m)**3

                            end if

                            end associate
                        end do

                        dpcell(i,j)%plist(k)%xs=beta*umax*dt*(term1**2)*term2

                        term1=0.0_dp
                        term2=0.0_dp
                        if (dpcell(i,j)%list(k)%count/=0) then
                        do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                            y=>dpcell(i,j)%list(k)%interlist(2,m), &
                            pp=>dpcell(i,j)%list(k)%interlist(3,m))

                            if (dpcell(i,j)%plist(k)%tid/=4) then

                            term1=term1+dpcell(i,j)%list(k)%dist(m)/ &
                            (dpcell(i,j)%list(k)%count)
                        
                            term2=term2+(dpcell(i,j)%plist(k)%y- &
                            dpcell(y,x)%plist(pp)%y)/ &
                            dpcell(i,j)%list(k)%dist(m)**3

                            end if

                            end associate
                        end do
                        end if

                        dpcell(i,j)%plist(k)%ys=beta*umax*dt*(term1**2)*term2

                        term1=0.0_dp
                        term2=0.0_dp  

                        do m=1,dpcell(i,j)%list(k)%count
                            associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                            y=>dpcell(i,j)%list(k)%interlist(2,m), &
                            pp=>dpcell(i,j)%list(k)%interlist(3,m))

                            if (dpcell(i,j)%plist(k)%tid/=4) then

                            ! term1=term1+(Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                            ! dpcell(i,j)%list(k)%dist(m),h1)*(dpcell(y,x)%plist(pp)%vx- &
                            ! dpcell(i,j)%plist(k)%vx))&
                            ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                            term1=term1+((dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                            dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                            dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                            (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                            *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                            ! term2=term2+(Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                            ! dpcell(i,j)%list(k)%dist(m),h1)*(dpcell(y,x)%plist(pp)%vx- &
                            ! dpcell(i,j)%plist(k)%vx))&
                            ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                            term2=term2+((dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                            dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                            dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                            (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                            *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                            end if

                            end associate
                        end do

                        dpcell(i,j)%plist(k)%vxs=term1*dpcell(i,j)%plist(k)%xs &
                                                        +term2*dpcell(i,j)%plist(k)%ys


                        term1=0.0_dp
                        term2=0.0_dp

                        do m=1,dpcell(i,j)%list(k)%count
                            associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                            y=>dpcell(i,j)%list(k)%interlist(2,m), &
                            pp=>dpcell(i,j)%list(k)%interlist(3,m))

                            if (dpcell(i,j)%plist(k)%tid/=4) then

                            ! term1=term1+(Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                            ! dpcell(i,j)%list(k)%dist(m),h1)*(dpcell(y,x)%plist(pp)%vy- &
                            ! dpcell(i,j)%plist(k)%vy))&
                            ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                            term1=term1+((dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                            dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                            dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                            (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                            *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                            ! term2=term2+(Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                            ! dpcell(i,j)%list(k)%dist(m),h1)*(dpcell(y,x)%plist(pp)%vy- &
                            ! dpcell(i,j)%plist(k)%vy))&
                            ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                            term2=term2+((dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                            dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                            dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                            (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                            *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                            end if

                            end associate
                        end do

                        dpcell(i,j)%plist(k)%vys=term1*dpcell(i,j)%plist(k)%xs &
                                                        +term2*dpcell(i,j)%plist(k)%ys


                        ! dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+rx
                        ! dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+ry

                        end if
                    end do  
                                    
                    end if
                end do
            end do
        !$omp end parallel do

        ! Shifting hydrodynamic values
        !$omp parallel do schedule(dynamic) private(i,j,k) collapse(2) default(shared)
            do j=sx,ex
            do i=sy,ey
                if (dpcell(i,j)%ptot/=0) then
                do k=1,dpcell(i,j)%ptot
                    if (dpcell(i,j)%plist(k)%tid>2) then
                    dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+dpcell(i,j)%plist(k)%xs
                    dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+dpcell(i,j)%plist(k)%ys
                    dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx+dpcell(i,j)%plist(k)%vxs
                    dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy+dpcell(i,j)%plist(k)%vys
                    ! dpcell(i,j)%plist(k)%pressure=dpcell(i,j)%plist(k)%gradvx
                    end if
                end do
                end if
            end do
            end do
        !$omp end parallel do   
        
    end subroutine exp_shift

    subroutine opt_shift

        implicit none

        integer :: pint=0
        
        ! Optimized Particle Shifting
        !$omp do private(i,k,m,term1,term2,normx,normy) schedule (runtime) collapse(2)
            do j=sx,ex
            do i=sy,ey
            if(dpcell(i,j)%ptot/=0) then

            do k=1,dpcell(i,j)%ptot

                if (dpcell(i,j)%plist(k)%tid==3) then

                dpcell(i,j)%plist(k)%xs=0.0_dp
                dpcell(i,j)%plist(k)%ys=0.0_dp
                dpcell(i,j)%plist(k)%vxs=0.0_dp
                dpcell(i,j)%plist(k)%vys=0.0_dp

                if ((dpcell(i,j)%list(k)%count/=0) &
                .and.(dpcell(i,j)%pplist(k)%gradvx>=0.69_dp)) then
                    if ((dpcell(i,j)%plist(k)%free).or.(dpcell(i,j)%plist(k)%vicinity)) then
                    term1=0.0_dp
                    term2=0.0_dp
                    normx=0.0_dp
                    normy=0.0_dp
                    dpcell(i,j)%plist(k)%xs=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4) then

                        normx=normx-(dpcell(i,j)%pplist(k)%coff(1)*&
                        Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        normy=normy-(dpcell(i,j)%pplist(k)%coff(3)*&
                        Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if
            
                        end associate
                    end do

                    term1=sqrt(normx**2+normy**2)
                    normx=normx/term1
                    normy=normy/term1

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4) then

                        dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-csh*(h1**2)* &
                        ((1-normx**2)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)-normx*normy* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%xs=min(abs(dpcell(i,j)%plist(k)%xs),maxshift*dl)* &
                    dpcell(i,j)%plist(k)%xs/abs(dpcell(i,j)%plist(k)%xs)

                    term1=0.0_dp
                    term2=0.0_dp
                    dpcell(i,j)%plist(k)%ys=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4) then

                        dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%ys-csh*(h1**2)*&
                        (-normx*normy*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+(1-normy**2)* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%ys=min(abs(dpcell(i,j)%plist(k)%ys),maxshift*dl)* &
                    dpcell(i,j)%plist(k)%ys/abs(dpcell(i,j)%plist(k)%ys)

                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vxs=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys



                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1)&
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vys=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys

                    dpcell(i,j)%plist(k)%vicinity=.false.

                    else


                    dpcell(i,j)%plist(k)%xs=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4) then

                        dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-csh*(h1**2)*&
                        (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)* &
                        (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),h1)/ &
                        Wab(2*dl,h1))**4)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%xs=min(abs(dpcell(i,j)%plist(k)%xs),maxshift*dl)*&
                    dpcell(i,j)%plist(k)%xs/abs(dpcell(i,j)%plist(k)%xs)

                    dpcell(i,j)%plist(k)%ys=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4) then

                        dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%ys-csh*(h1**2)*&
                        (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)* &
                        (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),h1)/ &
                        Wab(2*dl,h1))**4)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%ys=min(abs(dpcell(i,j)%plist(k)%ys),maxshift*dl)*&
                    dpcell(i,j)%plist(k)%ys/abs(dpcell(i,j)%plist(k)%ys)


                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vxs=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys


                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vys=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys

                    end if
                end if

                end if

            end do
            end if
            end do
            end do
        !$omp end do 

        ! Shifting hydrodynamic values
        !$omp do schedule(dynamic) private(i,j,k) collapse(2)  
            do j=sx,ex
            do i=sy,ey
                if (dpcell(i,j)%ptot/=0) then
                do k=1,dpcell(i,j)%ptot
                    if (dpcell(i,j)%plist(k)%tid==3) then
                    dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+dpcell(i,j)%plist(k)%xs
                    dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+dpcell(i,j)%plist(k)%ys
                    dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx+dpcell(i,j)%plist(k)%vxs
                    dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy+dpcell(i,j)%plist(k)%vys

                    end if
                end do
                end if
            end do
            end do
        !$omp end do   
        
    end subroutine opt_shift

    subroutine opt2_shift

        implicit none

        integer :: pint=0
        real(dp) :: heff,frac=0.030_dp!0.0420_dp
        
        ! Optimized Particle Shifting
        !$omp do private(i,k,m,term1,term2,normx,normy,heff) schedule (runtime) collapse(2)
            do j=sx,ex
            do i=sy,ey
            if(dpcell(i,j)%ptot/=0) then

            do k=1,dpcell(i,j)%ptot

                if (dpcell(i,j)%plist(k)%tid==3) then

                dpcell(i,j)%plist(k)%xs=0.0_dp
                dpcell(i,j)%plist(k)%ys=0.0_dp
                dpcell(i,j)%plist(k)%vxs=0.0_dp
                dpcell(i,j)%plist(k)%vys=0.0_dp

                ! if (dpcell(i,j)%pplist(k)%inpore) then

                heff=h1!(max(prrealx,prrealy)*hfac)/(2.0_dp*sqrt(dpcell(i,j)%pplist(k)%porosity))!sqrt(por)

                ! else 

                ! heff=(max(prrealx,prrealy)*4.80_dp)/(2.0_dp)

                ! end if


                if ((dpcell(i,j)%list(k)%count/=0)) then ! &.and.(dpcell(i,j)%pplist(k)%gradvx>=0.80_dp)
                    if ((dpcell(i,j)%plist(k)%free)) then !.or.(dpcell(i,j)%plist(k)%vicinity)
                    term1=0.0_dp
                    term2=0.0_dp
                    normx=0.0_dp
                    normy=0.0_dp
                    dpcell(i,j)%plist(k)%xs=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        normx=normx-(dpcell(i,j)%pplist(k)%coff(1)*&
                        Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        normy=normy-(dpcell(i,j)%pplist(k)%coff(3)*&
                        Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if
            
                        end associate
                    end do

                    term1=sqrt(normx**2+normy**2)
                    normx=normx/term1
                    normy=normy/term1

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        ! dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-csh*(heff**2)* &
                        ! ((1-normx**2)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),heff)-normx*normy* &
                        ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
                        ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)* &
                        ! (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
                        ! Wab(dl,heff))**4)

                        dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-frac*csh*(heff**2)*&
                        (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)!* &
                        ! (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
                        ! Wab(dl1,heff))**4)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%xs=min(abs(dpcell(i,j)%plist(k)%xs),maxshift*dl)* &
                    dpcell(i,j)%plist(k)%xs/abs(dpcell(i,j)%plist(k)%xs)

                    term1=0.0_dp
                    term2=0.0_dp
                    dpcell(i,j)%plist(k)%ys=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        ! dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%ys-csh*(heff**2)*&
                        ! (-normx*normy*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),heff)+(1-normy**2)* &
                        ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
                        ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)* &
                        ! (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
                        ! Wab(dl,heff))**4)

                        dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%ys-frac*csh*(heff**2)*&
                        (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)!* &
                        ! (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
                        ! Wab(dl1,heff))**4)


                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%ys=min(abs(dpcell(i,j)%plist(k)%ys),maxshift*dl)* &
                    dpcell(i,j)%plist(k)%ys/abs(dpcell(i,j)%plist(k)%ys)

                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        ! if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        ! (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        ! term1=term1+((Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),h1))*&
                        ! (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        ! term2=term2+((Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),h1))*&
                        ! (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        ! end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vxs=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys



                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        ! if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        ! (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1)&
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        ! term1=term1+((Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),h1))*&
                        ! (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        ! term2=term2+((Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),h1))*&
                        ! (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        ! end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vys=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys

                    dpcell(i,j)%plist(k)%vicinity=.false.

                    else


                    dpcell(i,j)%plist(k)%xs=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-csh*(heff**2)*&
                        (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)* &
                        (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
                        Wab(dl1,heff))**4)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%xs=min(abs(dpcell(i,j)%plist(k)%xs),maxshift*dl)*&
                    dpcell(i,j)%plist(k)%xs/abs(dpcell(i,j)%plist(k)%xs)

                    dpcell(i,j)%plist(k)%ys=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%ys-csh*(heff**2)*&
                        (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)* &
                        (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
                        Wab(dl1,heff))**4)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%ys=min(abs(dpcell(i,j)%plist(k)%ys),maxshift*dl)*&
                    dpcell(i,j)%plist(k)%ys/abs(dpcell(i,j)%plist(k)%ys)


                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        ! if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        ! (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        ! term1=term1+((Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),h1))*&
                        ! (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        ! term2=term2+((Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),h1))*&
                        ! (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        ! end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vxs=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys


                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        ! if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        ! (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        ! term1=term1+((Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),h1))*&
                        ! (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        ! term2=term2+((Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),h1))*&
                        ! (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        ! end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vys=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys

                    end if
                end if

                end if

            end do
            end if
            end do
            end do
        !$omp end do 

        ! Shifting hydrodynamic values
        !$omp do schedule(dynamic) private(i,j,k) collapse(2)  
            do j=sx,ex
            do i=sy,ey
                if (dpcell(i,j)%ptot/=0) then
                do k=1,dpcell(i,j)%ptot
                    if (dpcell(i,j)%plist(k)%tid==3) then
                    dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+dpcell(i,j)%plist(k)%xs
                    dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+dpcell(i,j)%plist(k)%ys
                    dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx+dpcell(i,j)%plist(k)%vxs
                    dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy+dpcell(i,j)%plist(k)%vys

                    end if
                end do
                end if
            end do
            end do
        !$omp end do   
        
    end subroutine opt2_shift

    subroutine compden

        implicit none

        real(dp) :: heff

        heff=h1

            !$omp do schedule(runtime) private(k,m,i,j,term1) collapse(2)
            do j=sx,ex
                do i=sy,ey
                if (dpcell(i,j)%ptot/=0) then                
                    do k=1,dpcell(i,j)%ptot
                    if ((dpcell(i,j)%plist(k)%tid/=2)) then
                        if(dpcell(i,j)%list(k)%count/=0) then
                            dpcell(i,j)%plist(k)%tden=0.0_dp
                            term1=0.0_dp
                            do m=1,dpcell(i,j)%list(k)%count                           
                                associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                                    y=>dpcell(i,j)%list(k)%interlist(2,m), &
                                    pp=>dpcell(i,j)%list(k)%interlist(3,m))
                                
                                    term1=term1+ &
                                    Wab(dpcell(i,j)%list(k)%dist(m),h1)*dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%pplist(pp)%porosity

                                end associate
                            end do
                        end if

                        dpcell(i,j)%plist(k)%tden=term1+&
                        dpcell(i,j)%plist(k)%mass*Wo(h1)/dpcell(i,j)%pplist(k)%porosity

                    else if ((dpcell(i,j)%plist(k)%tid==2)) then

                    dpcell(i,j)%plist(k)%tden=rho

                    end if
                    end do
                    
                end if
                end do
            end do
            !$omp end do

            ! !$omp do schedule(runtime) private(k,m,i,j,term1,term2) collapse(2)
            ! do j=sx,ex
            !     do i=sy,ey
            !     if (dpcell(i,j)%ptot/=0) then                
            !         do k=1,dpcell(i,j)%ptot
            !         if ((dpcell(i,j)%plist(k)%tid==3)) then
            !             dpcell(i,j)%pplist(k)%Kstar=0.0_dp
            !             if(dpcell(i,j)%list(k)%count/=0) then
            !                 term1=0.0_dp
            !                 term2=0.0_dp                        
            !                 do m=1,dpcell(i,j)%list(k)%count                           
            !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
            !                         y=>dpcell(i,j)%list(k)%interlist(2,m), &
            !                         pp=>dpcell(i,j)%list(k)%interlist(3,m))

            !                         term1=term1+ &
            !                         ((4*(dpcell(y,x)%plist(pp)%mass/(dpcell(y,x)%plist(pp)%density+dpcell(i,j)%plist(k)%density))&
            !                         *(dpcell(i,j)%pplist(k)%tden- &
            !                 dpcell(y,x)%pplist(pp)%tden)*(dpcell(i,j)%plist(k)%x- &
            !                 dpcell(y,x)%plist(pp)%x)* &
            !                 (dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
            !                 dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)* &
            !                 Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))/ &
            !                 ((dpcell(i,j)%list(k)%dist(m)**2+lam)))

            !                 term2=term2+ &
            !                 ((4*(dpcell(y,x)%plist(pp)%mass/(dpcell(y,x)%plist(pp)%density+dpcell(i,j)%plist(k)%density))&
            !                         *(dpcell(i,j)%pplist(k)%tden- &
            !                 dpcell(y,x)%pplist(pp)%tden)*(dpcell(i,j)%plist(k)%y- &
            !                 dpcell(y,x)%plist(pp)%y)* &
            !                 (dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
            !                 dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)* &
            !                 Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))/ &
            !                 ((dpcell(i,j)%list(k)%dist(m)**2+lam)))

            !                     end associate
            !                 end do

            !                 dpcell(i,j)%pplist(k)%Kstar=term1+term2
            !             end if

            !         end if
            !         end do
                    
            !     end if
            !     end do
            ! end do
            ! !$omp end do 

            ! !$omp do schedule(runtime) private(k,m,i,j) collapse(2)
            ! do j=sx,ex
            !     do i=sy,ey
            !     if (dpcell(i,j)%ptot/=0) then                
            !         do k=1,dpcell(i,j)%ptot
            !         if ((dpcell(i,j)%plist(k)%tid==3)) then

            !         dpcell(i,j)%pplist(k)%Kstar=((dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)- &
            !         dpcell(i,j)%pplist(k)%tden)/dpcell(i,j)%pplist(k)%Kstar

            !         if (dpcell(i,j)%pplist(k)%Kstar>(0.50_dp*(heff**2))) then

            !         dpcell(i,j)%pplist(k)%Kstar=0.50_dp*(heff**2)


            !         end if

            !         end if
            !         end do
                    
            !     end if
            !     end do
            ! end do
            ! !$omp end do 

    end subroutine

    subroutine opt3_shift

        implicit none

        integer :: pint=0
        real(dp) :: heff,pscshift

        pscshift=0.0_dp
        
        ! Optimized Particle Shifting
        !$omp do private(i,k,m,term1,term2,normx,normy,heff) schedule (runtime) collapse(2)
            do j=sx,ex
            do i=sy,ey
            if(dpcell(i,j)%ptot/=0) then

            do k=1,dpcell(i,j)%ptot

                if (dpcell(i,j)%plist(k)%tid==3) then

                dpcell(i,j)%plist(k)%xs=0.0_dp
                dpcell(i,j)%plist(k)%ys=0.0_dp
                dpcell(i,j)%plist(k)%vxs=0.0_dp
                dpcell(i,j)%plist(k)%vys=0.0_dp

                heff=h1

                if ((dpcell(i,j)%list(k)%count/=0) &
                .and.(dpcell(i,j)%pplist(k)%gradvx>=0.80_dp)) then
                    if ((dpcell(i,j)%plist(k)%free).or.(dpcell(i,j)%plist(k)%vicinity)) then
                    term1=0.0_dp
                    term2=0.0_dp
                    normx=0.0_dp
                    normy=0.0_dp
                    dpcell(i,j)%plist(k)%xs=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        normx=normx-(dpcell(i,j)%pplist(k)%coff(1)*&
                        Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        normy=normy-(dpcell(i,j)%pplist(k)%coff(3)*&
                        Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if
            
                        end associate
                    end do

                    term1=sqrt(normx**2+normy**2)
                    normx=normx/term1
                    normy=normy/term1

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-dpcell(i,j)%pplist(k)%Kstar* &
                        ((1.0_dp-normx**2)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)-normx*normy* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%xs=min(abs(dpcell(i,j)%plist(k)%xs),maxshift*dl)* &
                    dpcell(i,j)%plist(k)%xs/abs(dpcell(i,j)%plist(k)%xs)

                    term1=0.0_dp
                    term2=0.0_dp
                    dpcell(i,j)%plist(k)%ys=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%ys-dpcell(i,j)%pplist(k)%Kstar*&
                        (-normx*normy*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+(1.0_dp-normy**2)* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%ys=min(abs(dpcell(i,j)%plist(k)%ys),maxshift*dl)* &
                    dpcell(i,j)%plist(k)%ys/abs(dpcell(i,j)%plist(k)%ys)

                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vxs=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys



                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1)&
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vys=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys

                    dpcell(i,j)%plist(k)%vicinity=.false.

                    else


                    dpcell(i,j)%plist(k)%xs=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-dpcell(i,j)%pplist(k)%Kstar*&
                        (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)* &
                        (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
                        Wab(2*dl,heff))**4)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%xs=min(abs(dpcell(i,j)%plist(k)%xs),maxshift*dl)*&
                    dpcell(i,j)%plist(k)%xs/abs(dpcell(i,j)%plist(k)%xs)

                    dpcell(i,j)%plist(k)%ys=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%ys-dpcell(i,j)%pplist(k)%Kstar*&
                        (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)* &
                        (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
                        Wab(2*dl,heff))**4)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%ys=min(abs(dpcell(i,j)%plist(k)%ys),maxshift*dl)*&
                    dpcell(i,j)%plist(k)%ys/abs(dpcell(i,j)%plist(k)%ys)


                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vxs=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys


                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vys=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys

                    end if
                end if

                end if

            end do
            end if
            end do
            end do
        !$omp end do 

        ! Shifting hydrodynamic values
        !$omp do schedule(dynamic) private(i,j,k) collapse(2)  
            do j=sx,ex
            do i=sy,ey
                if (dpcell(i,j)%ptot/=0) then
                do k=1,dpcell(i,j)%ptot
                    if (dpcell(i,j)%plist(k)%tid==3) then
                    dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+dpcell(i,j)%plist(k)%xs
                    dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+dpcell(i,j)%plist(k)%ys
                    dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx+dpcell(i,j)%plist(k)%vxs
                    dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy+dpcell(i,j)%plist(k)%vys

                    end if
                end do
                end if
            end do
            end do
        !$omp end do   
        
    end subroutine opt3_shift

    subroutine exp1_shift

        implicit none

        ! integer :: pint=0

        real(dp) :: intx,inty,heff
        
        ! Optimized Particle Shifting
        !$omp do private(i,k,m,term1,term2,normx,normy,intx,inty,heff) schedule (runtime) collapse(2)
            do j=sx,ex
            do i=sy,ey
            if(dpcell(i,j)%ptot/=0) then

            do k=1,dpcell(i,j)%ptot

                if (dpcell(i,j)%plist(k)%tid==3) then

                dpcell(i,j)%plist(k)%xs=0.0_dp
                dpcell(i,j)%plist(k)%ys=0.0_dp
                dpcell(i,j)%plist(k)%vxs=0.0_dp
                dpcell(i,j)%plist(k)%vys=0.0_dp
                intx=0.0_dp
                inty=0.0_dp

                ! if (dpcell(i,j)%pplist(k)%inpore) then

                heff=(max(prrealx,prrealy)*4.80_dp)/(2.0_dp*sqrt(dpcell(i,j)%pplist(k)%porosity))

                ! else 

                ! heff=(max(prrealx,prrealy)*4.80_dp)/(2.0_dp)

                ! end if

                if ((dpcell(i,j)%list(k)%count/=0) &
                .and.(dpcell(i,j)%pplist(k)%gradvx>=0.70_dp)) then
                    if ((dpcell(i,j)%plist(k)%free).or.(dpcell(i,j)%plist(k)%vicinity)) then
                    term1=0.0_dp
                    term2=0.0_dp
                    normx=0.0_dp
                    normy=0.0_dp
                    dpcell(i,j)%plist(k)%xs=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4.and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        normx=normx-(dpcell(i,j)%pplist(k)%coff(1)*&
                        Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        normy=normy-(dpcell(i,j)%pplist(k)%coff(3)*&
                        Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if
            
                        end associate
                    end do

                    term1=sqrt(normx**2+normy**2)
                    normx=normx/term1
                    normy=normy/term1

                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4.and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                            term1=term1+dpcell(i,j)%list(k)%dist(m)/ &
                            (dpcell(i,j)%list(k)%count)
                        
                            term2=term2+(dpcell(i,j)%plist(k)%x- &
                            dpcell(y,x)%plist(pp)%x)/ &
                            dpcell(i,j)%list(k)%dist(m)**3

                        end if

                        end associate
                    end do

                    intx=beta*umax*dt*(term1**2)*term2

                    term1=0.0_dp
                    term2=0.0_dp
                    if (dpcell(i,j)%list(k)%count/=0) then
                    do m=1,dpcell(i,j)%list(k)%count
                    associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4.and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        term1=term1+dpcell(i,j)%list(k)%dist(m)/ &
                        (dpcell(i,j)%list(k)%count)
                    
                        term2=term2+(dpcell(i,j)%plist(k)%y- &
                        dpcell(y,x)%plist(pp)%y)/ &
                        dpcell(i,j)%list(k)%dist(m)**3

                        end if

                        end associate
                    end do
                    end if

                    inty=beta*umax*dt*(term1**2)*term2

                    dpcell(i,j)%plist(k)%xs=(1-normx**2)*intx-(normx*normy)*inty
                    dpcell(i,j)%plist(k)%ys=(1-normy**2)*inty-(normx*normy)*intx


                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4.and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vxs=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys

                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4.and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1)&
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vys=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys

                    dpcell(i,j)%plist(k)%vicinity=.false.

                    else


                    dpcell(i,j)%plist(k)%xs=0.0_dp

                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4.and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                            term1=term1+dpcell(i,j)%list(k)%dist(m)/ &
                            (dpcell(i,j)%list(k)%count)
                        
                            term2=term2+(dpcell(i,j)%plist(k)%x- &
                            dpcell(y,x)%plist(pp)%x)/ &
                            dpcell(i,j)%list(k)%dist(m)**3

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%xs=beta*umax*dt*(term1**2)*term2

                    dpcell(i,j)%plist(k)%ys=0.0_dp

                    term1=0.0_dp
                    term2=0.0_dp
                    if (dpcell(i,j)%list(k)%count/=0) then
                    do m=1,dpcell(i,j)%list(k)%count
                    associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4.and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        term1=term1+dpcell(i,j)%list(k)%dist(m)/ &
                        (dpcell(i,j)%list(k)%count)
                    
                        term2=term2+(dpcell(i,j)%plist(k)%y- &
                        dpcell(y,x)%plist(pp)%y)/ &
                        dpcell(i,j)%list(k)%dist(m)**3

                        end if

                        end associate
                    end do
                    end if

                    dpcell(i,j)%plist(k)%ys=beta*umax*dt*(term1**2)*term2


                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4.and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vxs=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys


                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4.and. &
                        (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vys=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys

                    end if
                end if

                end if

            end do
            end if
            end do
            end do
        !$omp end do 

        ! Shifting hydrodynamic values
        !$omp do schedule(dynamic) private(i,j,k) collapse(2)  
            do j=sx,ex
            do i=sy,ey
                if (dpcell(i,j)%ptot/=0) then
                do k=1,dpcell(i,j)%ptot
                    if (dpcell(i,j)%plist(k)%tid==3) then
                    dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+dpcell(i,j)%plist(k)%xs
                    dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+dpcell(i,j)%plist(k)%ys
                    ! dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx+dpcell(i,j)%plist(k)%vxs
                    ! dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy+dpcell(i,j)%plist(k)%vys

                    end if
                end do
                end if
            end do
            end do
        !$omp end do   
        
    end subroutine exp1_shift

    subroutine opt1_shift

        implicit none

        ! integer :: pint=0
        real(dp) :: A=2.0_dp,heff,frac=1.0_dp
        
        ! Optimized Particle Shifting
        !$omp do private(i,k,m,term1,term2,normx,normy,heff) schedule (runtime) collapse(2)
            do j=sx,ex
            do i=sy,ey
            if(dpcell(i,j)%ptot/=0) then

            do k=1,dpcell(i,j)%ptot

                if (dpcell(i,j)%plist(k)%tid>2) then

                dpcell(i,j)%plist(k)%xs=0.0_dp
                dpcell(i,j)%plist(k)%ys=0.0_dp
                dpcell(i,j)%plist(k)%vxs=0.0_dp
                dpcell(i,j)%plist(k)%vys=0.0_dp

                ! if (dpcell(i,j)%pplist(k)%inpore) then

                heff=h1

                ! else 

                ! heff=(max(prrealx,prrealy)*4.80_dp)/(2.0_dp)

                ! end if

                if ((dpcell(i,j)%list(k)%count/=0) &
                .and.(dpcell(i,j)%pplist(k)%gradvx>=0.70_dp)) then
                    if ((dpcell(i,j)%plist(k)%free).or.(dpcell(i,j)%plist(k)%vicinity)) then
                    term1=0.0_dp
                    term2=0.0_dp
                    normx=0.0_dp
                    normy=0.0_dp
                    dpcell(i,j)%plist(k)%xs=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        normx=normx-(dpcell(i,j)%pplist(k)%coff(1)*&
                        Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        normy=normy-(dpcell(i,j)%pplist(k)%coff(3)*&
                        Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)
            
                        end associate
                    end do

                    term1=sqrt(normx**2+normy**2)
                    normx=normx/term1
                    normy=normy/term1

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        ! dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-csh*(heff)*sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)*dt* &                        
                        ! ((1-normx**2)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),heff)-normx*normy* &
                        ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
                        ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-A*frac*(heff)*&
                        sqrt((dpcell(i,j)%plist(k)%vx)**2+(dpcell(i,j)%plist(k)%vy)**2)*dt* &                        
                        (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%xs=min(abs(dpcell(i,j)%plist(k)%xs),maxshift*dl)* &
                    dpcell(i,j)%plist(k)%xs/abs(dpcell(i,j)%plist(k)%xs)

                    term1=0.0_dp
                    term2=0.0_dp
                    dpcell(i,j)%plist(k)%ys=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        ! dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%ys-csh*(heff)*sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)*dt* &                    
                        ! (-normx*normy*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),heff)+(1-normy**2)* &
                        ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
                        ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-A*frac*(heff)*&
                        sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)*dt* &                        
                        (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)


                        end associate
                    end do

                    dpcell(i,j)%plist(k)%ys=min(abs(dpcell(i,j)%plist(k)%ys),maxshift*dl)* &
                    dpcell(i,j)%plist(k)%ys/abs(dpcell(i,j)%plist(k)%ys)

                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vxs=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys


                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))


                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vys=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys

                    dpcell(i,j)%plist(k)%vicinity=.false.

                    else


                    dpcell(i,j)%plist(k)%xs=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))


                        ! dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-csh*(heff)*sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)*dt* &
                        ! (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),heff))&
                        ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)* &
                        ! (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
                        ! Wab(2*dl,heff))**4)

                        dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-A*(heff)*&
                        sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)*dt* &                        
                        (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)


                        end associate
                    end do

                    dpcell(i,j)%plist(k)%xs=min(abs(dpcell(i,j)%plist(k)%xs),maxshift*dl)*&
                    dpcell(i,j)%plist(k)%xs/abs(dpcell(i,j)%plist(k)%xs)

                    dpcell(i,j)%plist(k)%ys=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))


                        ! dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%ys-csh*(heff)*sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)*dt* &
                        ! (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),heff))&
                        ! *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)* &
                        ! (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
                        ! Wab(2*dl,heff))**4)

                        dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-A*(heff)*&
                        sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)*dt* &                        
                        (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)


                        end associate
                    end do

                    dpcell(i,j)%plist(k)%ys=min(abs(dpcell(i,j)%plist(k)%ys),maxshift*dl)*&
                    dpcell(i,j)%plist(k)%ys/abs(dpcell(i,j)%plist(k)%ys)

                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vxs=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys


                    term1=0.0_dp
                    term2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vys=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys

                    end if
                end if

                end if

            end do
            end if
            end do
            end do
        !$omp end do 

        ! Shifting hydrodynamic values
        !$omp do schedule(dynamic) private(i,j,k) collapse(2)  
            do j=sx,ex
            do i=sy,ey
                if (dpcell(i,j)%ptot/=0) then
                do k=1,dpcell(i,j)%ptot
                    if (dpcell(i,j)%plist(k)%tid>2) then
                    dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+dpcell(i,j)%plist(k)%xs
                    dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+dpcell(i,j)%plist(k)%ys
                    dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx+dpcell(i,j)%plist(k)%vxs/dpcell(i,j)%pplist(k)%porosity
                    dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy+dpcell(i,j)%plist(k)%vys/dpcell(i,j)%pplist(k)%porosity
                    end if
                end do
                end if
            end do
            end do
        !$omp end do   
        
    end subroutine opt1_shift

    subroutine implicit_shift

        use functions
        use solver

        implicit none

        ! real(dp),dimension(:) :: pguess(finmax)
        real(dp) :: heff,frac=1.0_dp

        heff=h1

        ! pguess=0.0_dp


        !Density for internal particles
        !$omp parallel default(shared)
        !$omp do schedule(runtime) private(k,m,term1,i,j) collapse(2)
        do j=sx,ex
            do i=sy,ey
            if (dpcell(i,j)%ptot/=0) then
                
                do k=1,dpcell(i,j)%ptot
                if ((dpcell(i,j)%plist(k)%tid>2)) then
                    if(dpcell(i,j)%list(k)%count/=0) then
                        dpcell(i,j)%pplist(k)%tden=0.0_dp
                        ! dpcell(i,j)%pplist(k)%Kstar=0.0_dp
                        do m=1,dpcell(i,j)%list(k)%count                           
                            associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                                y=>dpcell(i,j)%list(k)%interlist(2,m), &
                                pp=>dpcell(i,j)%list(k)%interlist(3,m))

                                if (dpcell(y,x)%plist(pp)%tid/=4) then
                            
                                dpcell(i,j)%pplist(k)%tden=dpcell(i,j)%pplist(k)%tden+ &
                                Wab(dpcell(i,j)%list(k)%dist(m),h1)*dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%pplist(pp)%porosity

                                end if

                            end associate
                        end do
                    end if

                    dpcell(i,j)%pplist(k)%tden=dpcell(i,j)%pplist(k)%tden+&
                    dpcell(i,j)%plist(k)%mass*Wo(h1)/dpcell(i,j)%pplist(k)%porosity

                elseif ((dpcell(i,j)%plist(k)%tid<=2)) then

                dpcell(i,j)%pplist(k)%tden=rho

                end if
                end do
                
            end if
            end do
        end do
        !$omp end do 

        !$omp do schedule(runtime) private(k,m,term1,i,j) collapse(2)
        do j=sx,ex
            do i=sy,ey
            if (dpcell(i,j)%ptot/=0) then
                
                do k=1,dpcell(i,j)%ptot
                if ((dpcell(i,j)%plist(k)%tid>2)) then
                    if(dpcell(i,j)%list(k)%count/=0) then
                        ! dpcell(i,j)%pplist(k)%tden=0.0_dp
                        dpcell(i,j)%pplist(k)%Kstar=0.0_dp
                        do m=1,dpcell(i,j)%list(k)%count                           
                            associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                                y=>dpcell(i,j)%list(k)%interlist(2,m), &
                                pp=>dpcell(i,j)%list(k)%interlist(3,m))

                                if (dpcell(y,x)%plist(pp)%tid/=4) then

                                dpcell(i,j)%pplist(k)%Kstar=dpcell(i,j)%pplist(k)%Kstar+ &
                                Wab(dpcell(i,j)%list(k)%dist(m),h1)* &
                                dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%pplist(pp)%tden

                                end if

                            end associate
                        end do
                    end if

                    dpcell(i,j)%pplist(k)%Kstar=dpcell(i,j)%pplist(k)%Kstar+ &
                    dpcell(i,j)%plist(k)%mass*Wo(h1)/dpcell(i,j)%pplist(k)%tden

                elseif ((dpcell(i,j)%plist(k)%tid<=2)) then

                dpcell(i,j)%pplist(k)%Kstar=1.0_dp

                end if
                end do
                
            end if
            end do
        end do
        !$omp end do 

        ! Preparing the coff matrix for fluid part in CSR
        !$omp do schedule(runtime) &
        !$omp private(m,term1,term2,k,i,j,lamp) collapse(2)
        do j=sx,ex
            do i=sy,ey
            if (dpcell(i,j)%ptot/=0) then

                do k=1,dpcell(i,j)%ptot
                        
                if (dpcell(i,j)%plist(k)%tid/=4) then
                associate(pos=>dpcell(i,j)%plist(k)%matid,num=>dpcell(i,j)%list(k)%count)
                    term1=0.0_dp
                    term2=0.0_dp
                    fmatrix(pos)%sz=num+1
                    fmatrix(pos)%val(:)=0.0_dp
                    fmatrix(pos)%col(:)=pos
                    fvec(pos)=0.0_dp
                    num2=dpcell(i,j)%list(k)%count
                    if((num/=0)) then   !.or.(dpcell(i,j)%plist(k)%free/=1)
                    do m=1,num2
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                            if (dpcell(y,x)%plist(pp)%tid/=4) then

                                if (dpcell(i,j)%plist(k)%free) then

                                    lamp=0.0_dp!(1.0_dp+dpcell(y,x)%plist(pp)%density/(dpcell(y,x)%pplist(k)%porosity*rhomax))
    
                                elseif((pll<dpcell(i,j)%pplist(k)%gradvx).and.(dpcell(i,j)%pplist(k)%gradvx<pul)) then
    
                                lamp=0.50_dp*(1.0_dp-cos(22.0_dp*(dpcell(i,j)%pplist(k)%gradvx-pll)/ &
                                                ((pul-pll)*7.0_dp)))
                                else 
    
                                    lamp =1.0_dp
    
                                end if

                            if (dpcell(i,j)%plist(k)%tid>2) then

                                term1=2*dpcell(y,x)%plist(pp)%mass*(dpcell(i,j)%pplist(k)%porosity**2)*(&
                                Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)/ &
                                (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                (dpcell(y,x)%pplist(pp)%tden*dpcell(i,j)%pplist(k)%tden))

                                term2=2*dpcell(y,x)%plist(pp)%mass*(dpcell(i,j)%pplist(k)%porosity**2)*(&
                                Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)/ &
                                (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                (dpcell(y,x)%pplist(pp)%tden*dpcell(i,j)%pplist(k)%tden))                               

                                fmatrix(pos)%val(m)=(-(term1+term2))*lamp
                                fmatrix(pos)%col(m)=dpcell(y,x)%plist(pp)%matid

                                fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(term1+term2)


                                ! fvec(pos)=(1.0_dp-dpcell(i,j)%pplist(k)%Kstar)*lamp/real(dt**2,dp) !lamp

                                if ((dpcell(i,j)%plist(k)%free).and.(dpcell(i,j)%plist(k)%near_boun)) then


                                    fmatrix(pos)%val(m)=(-(term1+term2))
                                    fvec(pos)=fvec(pos)+frac*csh*(heff**2)*&
                                    (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                                    dpcell(i,j)%list(k)%dist(m),heff))&
                                    *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)/&
                                    ((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k)))*dt+lam*dt) + &
                                    frac*csh*(heff**2)*&
                                    (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                                    dpcell(i,j)%list(k)%dist(m),heff))&
                                    *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)/&
                                    ((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k)))*(dt**2)+lam*(dt**2))

                                else

                                    fvec(pos)=(1.0_dp-dpcell(i,j)%pplist(k)%Kstar)*lamp/real(dt**2,dp)

                                end if



                            elseif (dpcell(i,j)%plist(k)%tid==2) then

                                term1=2*dpcell(y,x)%plist(pp)%mass*(&
                                Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)/ &
                                (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                (dpcell(y,x)%plist(pp)%density*dpcell(i,j)%pplist(k)%tden))

                                term2=2*dpcell(y,x)%plist(pp)%mass*(&
                                Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)/ &
                                (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                (dpcell(y,x)%plist(pp)%density*dpcell(i,j)%pplist(k)%tden))

                                

                                if (dpcell(y,x)%plist(pp)%tid<=2) then

                                fmatrix(pos)%val(m)=-(term1+term2)!*dpcell(y,x)%pplist(pp)%lamp
                                    fmatrix(pos)%col(m)=dpcell(y,x)%plist(pp)%matid
                                    fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(term1+term2)!*dpcell(y,x)%pplist(pp)%lamp

                                    fvec(pos)=(1.0_dp-dpcell(i,j)%pplist(k)%Kstar)/real(dt**2,dp)



                                end if



                                ! fmatrix(pos)%val(1)=-1.0_dp!(term1+term2)
                                ! fmatrix(pos)%col(1)=dpcell(i,j)%plist(k)%wall(1)!dpcell(y,x)%plist(pp)%matid
                                ! fmatrix(pos)%val(num+1)=1.0_dp!fmatrix(pos)%val(num+1)+(term1+term2)


                            elseif (dpcell(i,j)%plist(k)%tid==1) then

                                if (dpcell(y,x)%plist(pp)%tid==3) then

                                term1=2*dpcell(i,j)%plist(k)%mass*(dpcell(y,x)%pplist(pp)%porosity**2)*(&
                                Wabx(dpcell(i,j)%plist(k),dpcell(y,x)%plist(pp),dpcell(i,j)%list(k)%dist(m),h1)) &
                                *(-dpcell(i,j)%plist(k)%x+dpcell(y,x)%plist(pp)%x)/ &
                                (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                (dpcell(i,j)%pplist(k)%tden*dpcell(y,x)%pplist(pp)%tden))

                                term2=2*dpcell(i,j)%plist(k)%mass*(dpcell(y,x)%pplist(pp)%porosity**2)*(&
                                Waby(dpcell(i,j)%plist(k),dpcell(y,x)%plist(pp),dpcell(i,j)%list(k)%dist(m),h1)) &
                                *(-dpcell(i,j)%plist(k)%y+dpcell(y,x)%plist(pp)%y)/ &
                                (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                (dpcell(i,j)%pplist(k)%tden*dpcell(y,x)%pplist(pp)%tden))

                                ! term1=2*dpcell(y,x)%plist(pp)%mass*(dpcell(i,j)%pplist(k)%porosity**2)*(&
                                ! Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                ! *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)/ &
                                ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                ! (dpcell(y,x)%plist(pp)%density*dpcell(i,j)%pplist(k)%tden))

                                ! term2=2*dpcell(y,x)%plist(pp)%mass*(dpcell(i,j)%pplist(k)%porosity**2)*(&
                                ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                ! *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)/ &
                                ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                ! (dpcell(y,x)%plist(pp)%density*dpcell(i,j)%pplist(k)%tden))

                                    fmatrix(pos)%val(m)=-(term1+term2)!*dpcell(y,x)%pplist(pp)%lamp
                                fmatrix(pos)%col(m)=dpcell(y,x)%plist(pp)%matid
                                fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(term1+term2)!*dpcell(y,x)%pplist(pp)%lamp

                                    fvec(pos)=(1.0_dp-dpcell(i,j)%pplist(k)%Kstar)/real(dt**2,dp)

                                
                                elseif (dpcell(y,x)%plist(pp)%tid<=2) then

                                    term1=2*dpcell(y,x)%plist(pp)%mass*(&
                                    Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                    *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)/ &
                                    (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                    (dpcell(y,x)%pplist(pp)%tden*dpcell(i,j)%pplist(k)%tden))

                                    term2=2*dpcell(y,x)%plist(pp)%mass*(&
                                    Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                    *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)/ &
                                    (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                    (dpcell(y,x)%pplist(pp)%tden*dpcell(i,j)%pplist(k)%tden))

                                    fmatrix(pos)%val(m)=-(term1+term2)!*dpcell(y,x)%pplist(pp)%lamp
                                    fmatrix(pos)%col(m)=dpcell(y,x)%plist(pp)%matid
                                    fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(term1+term2)!*dpcell(y,x)%pplist(pp)%lamp

                                    fvec(pos)=(1.0_dp-dpcell(i,j)%pplist(k)%Kstar)/real(dt**2,dp)
                                    ! fvec(pos)=fvec(pos)+(term1+term2)/real(dt,dp)  !lamp
                                    

                                end if
                            
                            end if

                        end if


                        end associate
                    end do
                    end if

                    fmatrix(pos)%col(num+1)=pos
                    if ((num==0)) then  !.or.( dpcell(i,j)%plist(k)%free==1)
                    fmatrix(pos)%val(num+1)=1.0_dp
                    end if

                    term1=fmatrix(pos)%val(num+1)

                    ! if (dpcell(i,j)%plist(k)%tid/=2) then

                    fvec(pos)=fvec(pos)/real(term1,dp)   ! Jacobi Preconditioning
                    
                    do m=1,ceiling(fac2*fplistmax)
                    fmatrix(pos)%val(m)=fmatrix(pos)%val(m)/real(term1,dp) !Jacobi Preconditioning
                    end do 

                    ! end if


                end associate

                end if

                end do

            end if            
            end do
        end do
        !$omp end do
        !$omp end parallel

        call pformat(fmatrix,fval,frow,fcol,finmax)
        call bicgstab(tl,pguess,finmax,fval,frow,fcol,fvec,fsol)

        ! Assigning pressures
        !$omp parallel default(shared)
        !$omp do schedule(runtime) private(i,k) collapse(2)       
            do j=sx,ex
                do i=sy,ey
                    if (dpcell(i,j)%ptot/=0) then

                        do k=1,dpcell(i,j)%ptot

                        if (dpcell(i,j)%plist(k)%tid/=4) then

                        dpcell(i,j)%pplist(k)%pshift=fsol(dpcell(i,j)%plist(k)%matid)

                        if (dpcell(i,j)%plist(k)%tid<=2) then

                            dpcell(i,j)%pplist(k)%pshift=dpcell(i,j)%pplist(k)%pshift&
                            +2*rho*abs(g)*min(prrealx,prrealy) 

                        end if

                        end if

                        end do

                    end if
                end do
            end do
        !$omp end do 

        ! New velocity calculations for particles 
        !$omp do private(m,term1,i,k) schedule (runtime) collapse(2)  
        do j=sx,ex
            do i=sy,ey
                if(dpcell(i,j)%ptot/=0) then

                do k=1,dpcell(i,j)%ptot
                    if (dpcell(i,j)%plist(k)%tid==3) then
                    !New fluid particle velocities
                    dpcell(i,j)%plist(k)%vxs=0.0d0
                    dpcell(i,j)%plist(k)%vys=0.0d0
                    if (dpcell(i,j)%list(k)%count/=0) then
                    do m=1,dpcell(i,j)%list(k)%count
                    associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4) then

                        term1=((dpcell(y,x)%plist(pp)%mass*dpcell(i,j)%pplist(k)%porosity/dpcell(y,x)%pplist(pp)%tden)*&
                        ((dpcell(y,x)%pplist(pp)%pshift)+(dpcell(i,j)%pplist(k)%pshift))* &
                        (dpcell(i,j)%pplist(k)%coff(1) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))/((dpcell(i,j)%pplist(k)%tden))

                        ! term1=((dpcell(y,x)%plist(pp)%mass*dpcell(i,j)%pplist(k)%porosity/dpcell(y,x)%pplist(pp)%tden)*&
                        ! ((dpcell(y,x)%pplist(pp)%pshift)+(dpcell(i,j)%pplist(k)%pshift))* &
                        ! (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),h1)))/((dpcell(i,j)%pplist(k)%tden))


                        dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vxs-term1*dt


                        term1=((dpcell(y,x)%plist(pp)%mass*dpcell(i,j)%pplist(k)%porosity/dpcell(y,x)%pplist(pp)%tden)*&
                        ((dpcell(y,x)%pplist(pp)%pshift)+(dpcell(i,j)%pplist(k)%pshift))* &
                        (dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))/((dpcell(i,j)%pplist(k)%tden))

                        ! term1=((dpcell(y,x)%plist(pp)%mass*dpcell(i,j)%pplist(k)%porosity/dpcell(y,x)%pplist(pp)%tden)*&
                        ! ((dpcell(y,x)%pplist(pp)%pshift)+(dpcell(i,j)%pplist(k)%pshift))* &
                        ! (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        ! dpcell(i,j)%list(k)%dist(m),h1)))/((dpcell(i,j)%pplist(k)%tden))

                        dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vys-term1*dt

                        end if

                        end associate
                    end do
                    end if

                    end if

                end do

                end if
            end do
        end do
        !$omp end do  

        ! New position calculations for fluid particles
        !$omp do schedule (runtime) private(i,k,m,term1,term2) collapse(2)     
        do j=sx,ex
            do i=sy,ey
                if(dpcell(i,j)%ptot/=0) then

                do k=1,dpcell(i,j)%ptot

                    if (dpcell(i,j)%plist(k)%tid==3) then


                    dpcell(i,j)%plist(k)%xs=(dt*dpcell(i,j)%plist(k)%vxs)

                    dpcell(i,j)%plist(k)%ys=(dt*dpcell(i,j)%plist(k)%vys)
                    dpcell(i,j)%plist(k)%vxs=0.0_dp
                    dpcell(i,j)%plist(k)%vys=0.0_dp

                    term1=0.0d0
                    term2=0.0d0  

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if
                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vxs=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys


                    term1=0.0d0
                    term2=0.0d0

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4) then

                        term1=term1+((dpcell(i,j)%pplist(k)%coff(1)&
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        term2=term2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),&
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(dpcell(y,x)%plist(pp),&
                        dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
                        (dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy))&
                        *(dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                        end if
                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vys=term1*dpcell(i,j)%plist(k)%xs &
                                                    +term2*dpcell(i,j)%plist(k)%ys

                    end if

                end do

                end if
            end do
        end do
        !$omp end do

        ! Shifting hydrodynamic values
        !$omp do schedule(dynamic) private(i,j,k) collapse(2)  
            do j=sx,ex
            do i=sy,ey
                if (dpcell(i,j)%ptot/=0) then
                do k=1,dpcell(i,j)%ptot
                    if (dpcell(i,j)%plist(k)%tid==3) then
                    dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+dpcell(i,j)%plist(k)%xs
                    dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+dpcell(i,j)%plist(k)%ys
                    dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx+dpcell(i,j)%plist(k)%vxs
                    dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy+dpcell(i,j)%plist(k)%vys

                    end if
                end do
                end if
            end do
            end do
        !$omp end do
        !$omp end parallel         
        
    end subroutine implicit_shift
    

end module part_shift
