module integrator

    use initialize
    use particle
    use domain
    use kernel

    implicit none

    contains

    subroutine cellshift
        implicit none 
        
        integer :: i,j,k,m
        
        ! Preparing particle transfers to surrounding cells    
        !$omp do private(i,j,k) schedule (runtime) collapse(2)
            do j=sx,ex
                do i=sy,ey
                dpcell(i,j)%temfct=0
                dpcell(i,j)%tn%count=0
                dpcell(i,j)%ts%count=0
                dpcell(i,j)%te%count=0
                dpcell(i,j)%tw%count=0
                dpcell(i,j)%tne%count=0
                dpcell(i,j)%tnw%count=0
                dpcell(i,j)%tse%count=0
                dpcell(i,j)%tsw%count=0

                if (dpcell(i,j)%ptot/=0) then

                    do k=1,dpcell(i,j)%ptot ! must be carried out serially
                    if (dpcell(i,j)%plist(k)%tid==3) then
                    associate(x1=>dpcell(i,j)%plist(k)%x,y1=>dpcell(i,j)%plist(k)%y)
                        if((x1>=dpcell(i,j)%xleft).and.(x1<dpcell(i,j)%xright)) then
                        if((y1>=dpcell(i,j)%ybot).and.(y1<dpcell(i,j)%ytop))then
                            dpcell(i,j)%temfct=dpcell(i,j)%temfct+1
                            dpcell(i,j)%ftemp(dpcell(i,j)%temfct)=dpcell(i,j)%plist(k)
                        elseif(y1<dpcell(i,j)%ybot) then
                            dpcell(i,j)%ts%count=dpcell(i,j)%ts%count+1
                            dpcell(i,j)%ts%list(dpcell(i,j)%ts%count)=dpcell(i,j)%plist(k)
                        elseif(y1>=dpcell(i,j)%ytop) then
                            dpcell(i,j)%tn%count=dpcell(i,j)%tn%count+1
                            dpcell(i,j)%tn%list(dpcell(i,j)%tn%count)=dpcell(i,j)%plist(k)
                        end if
                        elseif(x1>=dpcell(i,j)%xright) then
                        if(y1>=dpcell(i,j)%ytop)then
                            dpcell(i,j)%tne%count=dpcell(i,j)%tne%count+1
                            dpcell(i,j)%tne%list(dpcell(i,j)%tne%count)=dpcell(i,j)%plist(k)
                        elseif((y1>=dpcell(i,j)%ybot).and.(y1<dpcell(i,j)%ytop)) then
                            dpcell(i,j)%te%count=dpcell(i,j)%te%count+1
                            dpcell(i,j)%te%list(dpcell(i,j)%te%count)=dpcell(i,j)%plist(k)
                        elseif(y1<dpcell(i,j)%ybot) then
                            dpcell(i,j)%tse%count=dpcell(i,j)%tse%count+1
                            dpcell(i,j)%tse%list(dpcell(i,j)%tse%count)=dpcell(i,j)%plist(k)
                        end if
                        elseif (x1<dpcell(i,j)%xleft) then
                        if (y1>=dpcell(i,j)%ytop) then
                            dpcell(i,j)%tnw%count=dpcell(i,j)%tnw%count+1
                            dpcell(i,j)%tnw%list(dpcell(i,j)%tnw%count)=dpcell(i,j)%plist(k)
                        elseif((y1>=dpcell(i,j)%ybot).and.(y1<dpcell(i,j)%ytop)) then
                            dpcell(i,j)%tw%count=dpcell(i,j)%tw%count+1
                            dpcell(i,j)%tw%list(dpcell(i,j)%tw%count)=dpcell(i,j)%plist(k)
                        elseif (y1<dpcell(i,j)%ybot) then
                            dpcell(i,j)%tsw%count=dpcell(i,j)%tsw%count+1
                            dpcell(i,j)%tsw%list(dpcell(i,j)%tsw%count)=dpcell(i,j)%plist(k)
                        end if
                        end if
                    end associate
                    end if
                    end do

                end if
                end do
            end do
        !$omp end do

        ! Particle transfers    
        !$omp do private(i,j,k) schedule (runtime) collapse(2)
            do j=sx,ex
                do i=sy,ey
                ! From here everything must be done sequentially

                dpcell(i,j)%ptot=dpcell(i,j)%btot+dpcell(i,j)%temfct

                if (dpcell(i,j)%temfct/=0) then
                    do k=1,dpcell(i,j)%temfct
                    dpcell(i,j)%plist(k+dpcell(i,j)%btot)=dpcell(i,j)%ftemp(k)
                    end do
                end if

                if ((dpcell(i+1,j)%tn%count/=0)) then
                    do k=1,dpcell(i+1,j)%tn%count
                    dpcell(i,j)%plist(k+dpcell(i,j)%ptot)=dpcell(i+1,j)%tn%list(k)
                    end do
                end if
                dpcell(i,j)%ptot=dpcell(i,j)%ptot+(dpcell(i+1,j)%tn%count)

                if (dpcell(i-1,j)%ts%count/=0) then
                    do k=1,dpcell(i-1,j)%ts%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i-1,j)%ts%list(k)
                    end do
                end if

                if (dpcell(i,j+1)%tw%count/=0) then
                    do k=1,dpcell(i,j+1)%tw%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i,j+1)%tw%list(k)
                    end do
                end if

                if (dpcell(i,j-1)%te%count/=0) then
                    do k=1,dpcell(i,j-1)%te%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i,j-1)%te%list(k)
                    end do
                end if

                if (dpcell(i-1,j-1)%tse%count/=0) then
                    do k=1,dpcell(i-1,j-1)%tse%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i-1,j-1)%tse%list(k)
                    end do
                end if

                if (dpcell(i+1,j-1)%tne%count/=0) then
                    do k=1,dpcell(i+1,j-1)%tne%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i+1,j-1)%tne%list(k)
                    end do
                end if

                if (dpcell(i-1,j+1)%tsw%count/=0) then
                    do k=1,dpcell(i-1,j+1)%tsw%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i-1,j+1)%tsw%list(k)
                    end do
                end if

                if (dpcell(i+1,j+1)%tnw%count/=0) then
                    do k=1,dpcell(i+1,j+1)%tnw%count
                    dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                    dpcell(i,j)%plist(dpcell(i,j)%ptot)=dpcell(i+1,j+1)%tnw%list(k)
                    end do
                end if

                end do
            end do
        !$omp end do
        
    end subroutine cellshift

    subroutine boun_vel()
        implicit none

        logical :: pcounter

        real(dp) :: t1

        integer :: i,j,k,m

        !$omp do schedule (runtime) private(m,t1,k,i,j,pcounter) collapse(2)
        do j=sx,ex
            do i=sy,ey 

            do k=1,dpcell(i,j)%ptot

            if (dpcell(i,j)%plist(k)%tid<=2) then
            pcounter=.false.
            t1=0.0_dp
            dpcell(i,j)%plist(k)%vx=0.0_dp
            dpcell(i,j)%plist(k)%vy=0.0_dp
            dpcell(i,j)%plist(k)%vxs=0.0_dp
            dpcell(i,j)%plist(k)%vys=0.0_dp

            do m=1,dpcell(i,j)%list(k)%count
                associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                y=>dpcell(i,j)%list(k)%interlist(2,m), &
                pp=>dpcell(i,j)%list(k)%interlist(3,m))

                if (dpcell(y,x)%plist(pp)%tid==3) then

                pcounter=.true.

                dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx- &
                    dpcell(y,x)%plist(pp)%mass*Wab(dpcell(i,j)%list(k)%dist(m),h1)* &
                    dpcell(y,x)%plist(pp)%vx/dpcell(y,x)%plist(pp)%density

                dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy- &
                    dpcell(y,x)%plist(pp)%mass*Wab(dpcell(i,j)%list(k)%dist(m),h1)* &
                    dpcell(y,x)%plist(pp)%vy/dpcell(y,x)%plist(pp)%density

                t1=t1+dpcell(y,x)%plist(pp)%mass*Wab(dpcell(i,j)%list(k)%dist(m),h1)/&
                        dpcell(y,x)%plist(pp)%density

                end if

                end associate

            end do
            if (t1>=1e-6) then
            dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx/t1
            ! dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vx
            dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy/t1            
            ! dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vy

            end if

            end if

            end do

            end do
        end do
        !$omp end do
        
    end subroutine boun_vel

    subroutine int_vel
        implicit none

        real(dp) :: vart,ps=1.0_dp,t1,t2,pvol,p_dist

        integer :: i,j,k,m

        vart=0.0_dp

        !$omp do schedule (runtime) private(m,t1,t2,k,i,j) collapse(2)
        do j=sx,ex
            do i=sy,ey 

            do k=1,dpcell(i,j)%ptot

            t1=0.0_dp
            t2=0.0_dp
            if(dpcell(i,j)%plist(k)%tid==3) then
            dpcell(i,j)%pplist(k)%varts=0.0_dp
            if (dpcell(i,j)%list(k)%count/=0) then

            if (dpcell(i,j)%plist(k)%free) then

            do m=1,dpcell(i,j)%list(k)%count
                associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                y=>dpcell(i,j)%list(k)%interlist(2,m), &
                pp=>dpcell(i,j)%list(k)%interlist(3,m))

                t1=t1+(dpcell(y,x)%plist(pp)%vx-dpcell(i,j)%plist(k)%vx)* &
                (dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)* &
                Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))* &
                (dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                t2=t2+(dpcell(y,x)%plist(pp)%vy-dpcell(i,j)%plist(k)%vy)* &
                (dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)* &
                Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))* &
                (dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)

                end associate
            end do

            dpcell(i,j)%pplist(k)%varts=ps*(dl**2)*sqrt(t1**2+t2**2)*&
            sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)/ &
                (umax+1e-6)

            end if


            end if


            end if

            end do

            end do
        end do
        !$omp end do

        ! Non-pressure velocity calculation for fluid particles
        !$omp do schedule (runtime) private(m,t1,t2,k,i,j,vart,pvol,p_dist) collapse(2)
        do j=sx,ex
            do i=sy,ey            
            
            do k=1,dpcell(i,j)%ptot

            t1=0.0_dp
            t2=0.0_dp
            if(dpcell(i,j)%plist(k)%tid==3) then
            dpcell(i,j)%plist(k)%vxs=0.0_dp
            dpcell(i,j)%plist(k)%vys=0.0_dp
            if (dpcell(i,j)%list(k)%count/=0) then

            do m=1,dpcell(i,j)%list(k)%count
                associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                y=>dpcell(i,j)%list(k)%interlist(2,m), &
                pp=>dpcell(i,j)%list(k)%interlist(3,m))

                pvol=dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density
                p_dist=(dpcell(i,j)%list(k)%dist(m)**2+lam)**(-1)

                ! if (dpcell(y,x)%plist(pp)%tid/=4) then

                if (dpcell(i,j)%plist(k)%free) then
                vart=dpcell(i,j)%pplist(k)%varts
                elseif ((.not.(dpcell(i,j)%plist(k)%free)).and. &
                (dpcell(y,x)%plist(pp)%free)) then
                vart=Wab(dpcell(i,j)%list(k)%dist(m),h1)*dpcell(y,x)%pplist(pp)%varts &
                        /Wo(h1)

                else
                vart=0.0_dp
                end if

                    t1=(2*pvol*(dpcell(i,j)%plist(k)%vx- &
                    dpcell(y,x)%plist(pp)%vx)*(dpcell(i,j)%plist(k)%x- &
                    dpcell(y,x)%plist(pp)%x)* &
                    (dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                    dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)* &
                    Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))*p_dist


                    t2=(2*pvol*(dpcell(i,j)%plist(k)%vx- &
                    dpcell(y,x)%plist(pp)%vx)*(dpcell(i,j)%plist(k)%y- &
                    dpcell(y,x)%plist(pp)%y)* &
                    (dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                    dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)* &
                    Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))*p_dist

                    dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vxs+ &
                    (t1+t2)*((2*(mu/dpcell(i,j)%plist(k)%density)* &
                    (mu/dpcell(y,x)%plist(pp)%density)/ &
                    ((mu/dpcell(i,j)%plist(k)%density)+(mu/dpcell(y,x)%plist(pp)%density)))+vart)

                    if (dpcell(y,x)%plist(pp)%tid>2) then

                    if ((dpcell(i,j)%pplist(k)%nut<1e-6).and.(dpcell(y,x)%pplist(pp)%nut<1e-6)) then
                        
                    dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vxs+ &
                    (t1+t2)*(dpcell(i,j)%pplist(k)%nut+dpcell(y,x)%pplist(pp)%nut)*0.50_dp
                    else 
                        dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vxs+ &
                    (t1+t2)*((2*dpcell(i,j)%pplist(k)%nut*dpcell(y,x)%pplist(pp)%nut)/ &
                    (dpcell(i,j)%pplist(k)%nut+dpcell(y,x)%pplist(pp)%nut))

                    end if

                    ! else

                    !     dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vxs+ &
                    ! (t1+t2)*(dpcell(i,j)%pplist(k)%nut)

                    end if


                    t1=(2*pvol*(dpcell(i,j)%plist(k)%vy- &
                    dpcell(y,x)%plist(pp)%vy)*(dpcell(i,j)%plist(k)%x- &
                    dpcell(y,x)%plist(pp)%x)* &
                    (dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                    dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)* &
                    Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))*p_dist


                    t2=(2*pvol*(dpcell(i,j)%plist(k)%vy- &
                    dpcell(y,x)%plist(pp)%vy)*(dpcell(i,j)%plist(k)%y- &
                    dpcell(y,x)%plist(pp)%y)* &
                    (dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                    dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)* &
                    Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))*p_dist

                    dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vys+ &
                    (t1+t2)*((2*(mu/dpcell(i,j)%plist(k)%density)* &
                    (mu/dpcell(y,x)%plist(pp)%density)/ &
                    ((mu/dpcell(i,j)%plist(k)%density)+(mu/dpcell(y,x)%plist(pp)%density)))+vart)

                    if (dpcell(y,x)%plist(pp)%tid>2) then

                    if ((dpcell(i,j)%pplist(k)%nut<1e-6).and.(dpcell(y,x)%pplist(pp)%nut<1e-6)) then
                        
                        dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vys+ &
                        (t1+t2)*(dpcell(i,j)%pplist(k)%nut+dpcell(y,x)%pplist(pp)%nut)*0.50_dp
                    else 
                            dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vys+ &
                        (t1+t2)*((2*dpcell(i,j)%pplist(k)%nut*dpcell(y,x)%pplist(pp)%nut)/ &
                        (dpcell(i,j)%pplist(k)%nut+dpcell(y,x)%pplist(pp)%nut))
    
                    end if

                    ! else

                    !     dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vys+ &
                    !     (t1+t2)*(dpcell(i,j)%pplist(k)%nut)

                    end if

                    ! Turbulence Calculations

                    if (dpcell(i,j)%plist(k)%tid==3) then

                    ! if (dpcell(y,x)%plist(pp)%tid==3) then

                    ! t1=(dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%pplist(pp)%tke+ &
                    ! dpcell(i,j)%pplist(k)%tke)* &
                    ! (dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                    ! dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)* &
                    ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))) &
                        ! /(dpcell(y,x)%plist(pp)%density)

                    t1=(pvol*(dpcell(y,x)%pplist(pp)%tke- &
                    dpcell(i,j)%pplist(k)%tke)* &
                    (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                    dpcell(i,j)%list(k)%dist(m),h1)))

                    ! dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vxs+(t1+t2)*dpcell(i,j)%pplist(k)%porosity
                    dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vxs- &
                    (t1)*dpcell(i,j)%pplist(k)%porosity*2.0_dp/3

                    ! t1=(dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%pplist(pp)%tke+ &
                    ! dpcell(i,j)%pplist(k)%tke)* &
                    ! (dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                    ! dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)* &
                    ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))) &
                    ! /(dpcell(y,x)%plist(pp)%density)

                    t1=(pvol*(dpcell(y,x)%pplist(pp)%tke- &
                    dpcell(i,j)%pplist(k)%tke)* &
                    (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                    dpcell(i,j)%list(k)%dist(m),h1)))

                    ! dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vys+(t1+t2)*dpcell(i,j)%pplist(k)%porosity
                    dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vys-&
                    (t1)*dpcell(i,j)%pplist(k)%porosity*2.0_dp/3

                    end if

                    ! end if

                end associate
            end do
            end if
                dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vx+dt* &
                dpcell(i,j)%plist(k)%vxs+dt*dpcell(i,j)%pplist(k)%resistx*dpcell(i,j)%pplist(k)%porosity

                dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vy+dt* &
                (dpcell(i,j)%plist(k)%vys+g*dpcell(i,j)%pplist(k)%porosity)+ &
                dt*dpcell(i,j)%pplist(k)%resisty*dpcell(i,j)%pplist(k)%porosity          
            end if

                end do
            end do
        end do
        !$omp end do
        
    end subroutine int_vel

    subroutine comp_vel
        implicit none

        integer :: i,j,k,m
        
        real(dp) :: t1 

        ! New velocity calculations for particles 
        !$omp do private(m,t1,i,k,j) schedule (runtime) collapse(2)
        do j=sx,ex
            do i=sy,ey
                if(dpcell(i,j)%ptot/=0) then
                
                do k=1,dpcell(i,j)%ptot
                    if (dpcell(i,j)%plist(k)%tid>2) then
                    !New fluid particle velocities
                    dpcell(i,j)%plist(k)%vx=0.0_dp
                    dpcell(i,j)%plist(k)%vy=0.0_dp
                    if (dpcell(i,j)%list(k)%count/=0) then

                    do m=1,dpcell(i,j)%list(k)%count
                    associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        ! t1=(dpcell(y,x)%pplist(pp)%porosity**2)*(dpcell(y,x)%plist(pp)%mass/(dpcell(y,x)%plist(pp)%density**2))* &
                        !     (dpcell(y,x)%plist(pp)%pressure-dpcell(i,j)%plist(k)%pressure)* &
                        !     (dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        ! dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)* &
                        ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))



                        ! t1=2*((dpcell(y,x)%plist(pp)%mass*dpcell(i,j)%pplist(k)%porosity/dpcell(y,x)%plist(pp)%density)*&
                        ! ((dpcell(y,x)%plist(pp)%pressure*dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)+&
                        ! (dpcell(i,j)%plist(k)%pressure*dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(pp)%porosity))* &
                        ! (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        ! dpcell(i,j)%list(k)%dist(m),h1))) &
                        ! /((dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(pp)%porosity)* &
                        ! ((dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(pp)%porosity) &
                        ! +(dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)))

                        t1=2*((dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)*(dpcell(i,j)%pplist(k)%porosity**2)*&
                        ((dpcell(y,x)%plist(pp)%pressure*dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)+&
                        (dpcell(i,j)%plist(k)%pressure*dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(pp)%porosity))* &
                        (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        dpcell(i,j)%list(k)%dist(m),h1)))/((dpcell(i,j)%plist(k)%density)* &
                        ((dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)+ &
                        (dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(pp)%porosity)))


                        dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx-t1*dt

                        ! t1=(dpcell(y,x)%pplist(pp)%porosity**2)*(dpcell(y,x)%plist(pp)%mass/(dpcell(y,x)%plist(pp)%density**2))* &
                        !     (dpcell(y,x)%plist(pp)%pressure-dpcell(i,j)%plist(k)%pressure)* &
                        !     (dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        ! dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)* &
                        ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))


                        ! t1=2*((dpcell(y,x)%plist(pp)%mass*dpcell(i,j)%pplist(k)%porosity/dpcell(y,x)%plist(pp)%density)*&
                        ! ((dpcell(y,x)%plist(pp)%pressure*dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)+&
                        ! (dpcell(i,j)%plist(k)%pressure*dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(pp)%porosity))* &
                        ! (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        ! dpcell(i,j)%list(k)%dist(m),h1))) &
                        ! /((dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(pp)%porosity)* &
                        ! ((dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(pp)%porosity) &
                        ! +(dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)))

                        t1=2*((dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)*(dpcell(i,j)%pplist(k)%porosity**2)*&
                        ((dpcell(y,x)%plist(pp)%pressure*dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)+&
                        (dpcell(i,j)%plist(k)%pressure*dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(pp)%porosity))* &
                        (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        dpcell(i,j)%list(k)%dist(m),h1)))/((dpcell(i,j)%plist(k)%density)* &
                        ((dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)+ &
                        (dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(pp)%porosity)))


                        dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy-t1*dt

                        end associate
                    end do
                    end if

                    dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vxs + &
                    dpcell(i,j)%plist(k)%vx

                    dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vys + &
                    dpcell(i,j)%plist(k)%vy

                    dpcell(i,j)%pplist(k)%vel=sqrt(dpcell(i,j)%plist(k)%vx**2+ &
                    dpcell(i,j)%plist(k)%vy**2)/dpcell(i,j)%pplist(k)%porosity

                    end if

                end do

                end if
            end do
        end do
        !$omp end do        
        
    end subroutine comp_vel

    subroutine comp_pos
        
        implicit none   
        
        integer :: i,j,k,m

        ! New position calculations for fluid particles
        !$omp do schedule (runtime) private(i,k,j) collapse(2)   
        do j=sx,ex
            do i=sy,ey
                if(dpcell(i,j)%ptot/=0) then

                do k=1,dpcell(i,j)%ptot

                    if (dpcell(i,j)%plist(k)%tid==3) then

                    dpcell(i,j)%plist(k)%x=(dpcell(i,j)%plist(k)%x+ &
                    dpcell(i,j)%plist(k)%xs+dt*dpcell(i,j)%plist(k)%vx/dpcell(i,j)%pplist(k)%porosity)/2.0_dp

                    dpcell(i,j)%plist(k)%y=(dpcell(i,j)%plist(k)%y + &
                    dpcell(i,j)%plist(k)%ys+dt*dpcell(i,j)%plist(k)%vy/dpcell(i,j)%pplist(k)%porosity)/2.0_dp


                    end if

                end do

                end if
            end do
        end do
        !$omp end do        
        
    end subroutine comp_pos

    subroutine timestep

        use,intrinsic :: ieee_arithmetic

        implicit none

        integer :: i,j,k,m
        
        ! Finding max velocity
        !$omp single
        umax=0.0_dp
        numax=0.0_dp
        !$omp end single

        !$omp do schedule(runtime) private(i,j) collapse(2)
            do j=sx,ex 
                do i=sy,ey
                    if (dpcell(i,j)%ptot/=0) then
                    dpcell(i,j)%maxvel=maxval(abs(dpcell(i,j)%pplist(1:dpcell(i,j)%ptot)%vel))
                    dpcell(i,j)%maxeddy=maxval(abs(dpcell(i,j)%pplist(1:dpcell(i,j)%ptot)%nut))
                    end if
                    ! umax=max(umax,dpcell(i,j)%maxvel)
                    ! numax=max(numax,dpcell(i,j)%maxeddy)
                end do
            end do
            
        !$omp end do

        !$omp do schedule(runtime) private(i,j) reduction(max:umax,numax) collapse(2)
            do j=sx,ex 
                do i=sy,ey

                    umax=max(umax,dpcell(i,j)%maxvel)
                    numax=max(numax,dpcell(i,j)%maxeddy)

                end do
            end do
            
        !$omp end do

        !$omp single
        if (t<=1.0_dp) then
        dt=min(sig1*(dl)/umax,sig1*((dl)**2)/((mu/rho)+numax),0.0010_dp)
        elseif ((t>1.0_dp).and.(t<=2.0_dp)) then
        dt=real(min(real(sig1*(dl)/umax),real(sig1*((dl)**2)*rho/mu),0.0015))
        ! elseif ((t>2.0_dp).and.(t<=3.0_dp)) then
        ! dt=real(min(real(sig1*(dl)/umax),real(sig1*((dl)**2)*rho/mu),0.002))
        ! elseif ((t>3.0_dp).and.(t<=4.0_dp)) then
        ! dt=real(min(real(sig1*(dl)/umax),real(sig1*((dl)**2)*rho/mu),0.004))
        else
        dt=min(sig1*(dl)/umax,sig1*((dl)**2)/((mu/rho)+numax),0.00150_dp)
        end if
        if((.not.(ieee_is_finite(numax))).or.(ieee_is_nan(numax))) then
        numax=0.0_dp
        end if
        dtsol=sig1*(dl**2)/(Dm+numax/tschmidt)
        if (dtsol<dt) then
        solsteps=ceiling(dt/dtsol)
        dtsol=dt/solsteps
        end if

        !$omp end single 
        
    end subroutine timestep
    

end module integrator
