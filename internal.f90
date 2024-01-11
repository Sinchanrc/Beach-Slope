module internal

    use initialize
    use particle

    implicit none

    integer :: i,j,k,m
    
    contains

    subroutine gen_fluid()

        implicit none

        integer :: step=0

        real(dp) :: bounlen,bounlen2

        bounlen=fpy*2*prrealy/sqrt(por)
        bounlen2=open_lhs*2*prrealy

        allocate(input(fpy))

        do i=fpy,1,-1 
            input(i)%y=((brrealy)*((2*bl)-1))+(fpy-i)*2*prrealy/sqrt(por)+prrealy/sqrt(por)
            input(i)%x=((brrealx)*((2*bl)-1))+(fpx-1)*2*prrealx/sqrt(por)+prrealx/sqrt(por)+domain_shift
        end do

        !$omp parallel default(shared)
        !Setting up fluid particle positions and pressure
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j =1,fpx+5
            do i=fpy,1,-1
                flist(i,j)%y=((brrealy)*((2*bl)-1))+(fpy-i)*2*prrealy/sqrt(por)+prrealy/sqrt(por)
                ! if (mod(i,2)==0) then
                ! flist(i,j)%x=((brrealx)*((2*bl)-1))+(j-1)*2*prrealx/sqrt(por)+2*prrealx/sqrt(por)
                ! else
                flist(i,j)%x=((brrealx)*((2*bl)-1))+(j-1)*2*prrealx/sqrt(por)+prrealx/sqrt(por)+domain_shift
                ! end if
            flist(i,j)%vx=entry_vel*3.0_dp*(((flist(i,j)%y-(brrealy)*(2*bl-1))/bounlen) &
                            -0.5_dp*((flist(i,j)%y-(brrealy)*(2*bl-1))/bounlen)**2)
            flist(i,j)%vy=0.0_dp

            if (j>fpx) then
                flist(i,j)%buffer=.true.
            end if

            ! if (i==1) then
            !     flist(i,j)%pressure=0.0_dp

            ! else
                flist(i,j)%pressure=(-flist(i,j)%y+((brrealy*distfac)*((2*bl)-1))+wc+prrealy/sqrt(por))*rho*abs(g)
            
            ! end if

            end do
            end do
        !$omp end do

        !Setting particle id for fluid particles
        !$omp single 
            count=0
            finmin=count+1
            do j =1,fpx+5
                do i=fpy,1,-1
    
                ! if ((mod(i,2)/=0)) then
                    count=count+1
                    flist(i,j)%pid=count
                    flist(i,j)%tid=3
                ! else 
                !     if (j/=fpx) then
                !         count=count+1
                !     flist(i,j)%pid=count
                !     flist(i,j)%tid=3
                !     else 

                !     flist(i,j)%tid=0

                !     end if

                ! end if

                if(((flist(i,j)%y-yl)-line_grad*(flist(i,j)%x-xl))>0.0) then
                    flist(i,j)%tid=0
                end if
                
                end do
            end do
        
        !$omp end single

        !Distributing fluid particles to cells
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j=2,cellx-1
                do i=2,celly-1
                do l1=1,fpx+5
                    do k=1,fpy
                    if ((flist(k,l1)%x>=dpcell(i,j)%xleft) .and. &
                        (flist(k,l1)%x<dpcell(i,j)%xright).and. &
                        (flist(k,l1)%y>=dpcell(i,j)%ybot).and. &
                        (flist(k,l1)%y<dpcell(i,j)%ytop).and.(flist(k,l1)%tid/=0))then
                        dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)=flist(k,l1)

                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%mass=fmass
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%density=rho

                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%oden=dpcell(i,j)%plist&
                        (dpcell(i,j)%ptot)%density
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%ovol=fmass/rho

                    end if
                    end do
                end do
                end do
            end do
        !$omp end do

        ! Assigning locations to entry buffers in RHS
        !$omp do schedule(runtime) private(i,j) collapse(2)
            do j=2,cellx-1
                do i=2,celly-1

                    do k=1,fpy
                    if ((input(k)%x>=dpcell(i,j)%xleft) .and. &
                        (input(k)%x<dpcell(i,j)%xright).and. &
                        (input(k)%y>=dpcell(i,j)%ybot).and. &
                        (input(k)%y<dpcell(i,j)%ytop))then

                            if (.not.(allocated(dpcell(i,j)%ebuffpt))) then
                                allocate(dpcell(i,j)%ebuffpt(7))
                            end if
                            dpcell(i,j)%entrybuff=.true.
                            dpcell(i,j)%ebuff=dpcell(i,j)%ebuff+1
                            dpcell(i,j)%ebuffpt(dpcell(i,j)%ebuff)=input(k)
                            

                    end if
                    end do

                end do
            end do
        !$omp end do

        !Counting number of cells on RHS
        !$omp single

            do j=2,cellx-1
                do i=2,celly-1

                    if (dpcell(i,j)%entrybuff) then

                        entrycounter1=entrycounter1+1

                    end if

                    if (dpcell(i,j)%exitbuff) then

                        exitcounter=exitcounter+1

                    end if

                end do
            end do

        !$omp end single 

        !$omp single

            allocate(entrycell1(entrycounter1),exitcell(exitcounter))
            step=1

            ! Pointing to cells containing entry points on RHS
            do j=2,cellx-1
                do i=2,celly-1

                    if (dpcell(i,j)%entrybuff) then

                        entrycell1(step)%bcell=>dpcell(i,j)
                        step=step+1

                    end if

                end do
            end do

            step=1

            ! Pointing to cells containing exitcells on LHS
            do j=2,cellx-1
                do i=2,celly-1

                    if (dpcell(i,j)%exitbuff) then

                        exitcell(step)%bcell=>dpcell(i,j)
                        step=step+1

                    end if

                end do
            end do
        
            deallocate(flist,input)
            fpy=floor(real(coastal_ht,dp)/(2*real(prrealy,dp)))!+1
            fpx=floor(real(2.0_dp,dp)/(2*real(prrealx,dp)))!+1
            allocate(flist(fpy,fpx))
        !$omp end single

!!!!!!!!!!!!!!!!!!! Setting up porous media(non-deform)!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Setting up fluid particle positions and pressure
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j =1,fpx
                do i=fpy,1,-1
    
    
                flist(i,j)%y=((brrealy)*((2*bl)-1))+(fpy-i)*2*prrealy+prrealy+set_ht
                ! if (mod(i,2)==0) then
                ! flist(i,j)%x=((brrealx)*((2*bl)-1))+(j-1)*2*prrealx+2*prrealx
                ! else
                flist(i,j)%x=((brrealx)*((2*bl)-1))+(j-1)*2*prrealx+prrealx+domain_shift
                ! end if
                flist(i,j)%vx=0.0_dp!-2.5_dp/(3600*24*por)
                flist(i,j)%vy=0.0_dp
                if (i==1) then
                    flist(i,j)%pressure=0.0_dp
    
                else
                    flist(i,j)%pressure=(-flist(i,j)%y+((brrealy*distfac)*((2*bl)-1))+coastal_ht &
                    +prrealy+set_ht)*rho*abs(g)
                
                end if

                end do
                end do
        !$omp end do

        !Setting particle id for fluid particles
        !$omp single 
            do j =1,fpx
                do i=fpy,1,-1

                    ! if ((mod(i,2)/=0)) then
                        count=count+1
                        flist(i,j)%pid=count
                        flist(i,j)%tid=3
                    ! else 
                    !     if (j/=fpx) then
                    !         count=count+1
                    !     flist(i,j)%pid=count 
                    !     flist(i,j)%tid=3
                    !     else 
    
                    !     flist(i,j)%tid=0
    
                    !     end if
    
                    ! end if
    
                ! if ((mod(i,2)/=0)) then
                    ! count=count+1
                    ! flist(i,j)%pid=count
                    ! flist(i,j)%tid=3
                ! else 

                        if(((flist(i,j)%y-yl+prrealy)-line_grad*(flist(i,j)%x-xl))<0.0) then
                            flist(i,j)%tid=0

                        else
                            flist(i,j)%y=flist(i,j)%y+2*dl

                        end if


                ! end if
                
                end do
            end do
        
        !$omp end single

        !Distributing fluid particles to cells
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j=2,cellx-1
                do i=2,celly-1
                do l1=1,fpx
                    do k=1,fpy

                    if ((flist(k,l1)%x>=dpcell(i,j)%xleft) .and. &
                        (flist(k,l1)%x<dpcell(i,j)%xright).and. &
                        (flist(k,l1)%y>=dpcell(i,j)%ybot).and. &
                        (flist(k,l1)%y<dpcell(i,j)%ytop).and.(flist(k,l1)%tid/=0))then
                        dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)=flist(k,l1)
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%mass=fmass*1.0
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%density=rho*1.0
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%oden=dpcell(i,j)%plist&
                        (dpcell(i,j)%ptot)%density

                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%ovol=fmass/rho

                    end if


                    end do
                end do
                end do
            end do
        !$omp end do

        !$omp single

            deallocate(flist)
            allocate(flist(open_lhs,5))

        !$omp end single

        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j =1,5
                do i=open_lhs,1,-1
                    flist(i,j)%y=2*bny*brrealy-brrealy+prrealy+(open_lhs-i)*2*prrealy

                    flist(i,j)%x=(j-1)*2*prrealx+brrealx

                flist(i,j)%vx=(-bounlen/bounlen2)*entry_vel*3.0_dp*(((flist(i,j)%y-(brrealy)*(2*bny-1))/bounlen2) &
                                -0.5_dp*((flist(i,j)%y-(brrealy)*(2*bny-1))/bounlen2)**2)
                flist(i,j)%vy=0.0_dp

                    flist(i,j)%buffer=.true.

                    flist(i,j)%pressure=(-flist(i,j)%y+((brrealy*distfac)*((2*bl)-1))+wc+prrealy)*rho*abs(g)

                end do
                end do
            !$omp end do

        !$omp single 
            do j =1,5
                do i=1,open_lhs

                        count=count+1
                        flist(i,j)%pid=count
                        flist(i,j)%tid=3

                
                end do
            end do
        
        !$omp end single


        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j=2,cellx-1
                do i=2,celly-1
                do l1=1,5
                    do k=1,open_lhs

                    if ((flist(k,l1)%x>=dpcell(i,j)%xleft) .and. &
                        (flist(k,l1)%x<dpcell(i,j)%xright).and. &
                        (flist(k,l1)%y>=dpcell(i,j)%ybot).and. &
                        (flist(k,l1)%y<dpcell(i,j)%ytop).and.(flist(k,l1)%tid/=0))then
                        dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)=flist(k,l1)
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%mass=fmass*1.0
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%density=rho*1.0
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%oden=dpcell(i,j)%plist&
                        (dpcell(i,j)%ptot)%density

                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%ovol=fmass/rho

                    end if


                    end do
                end do
                end do
            end do
        !$omp end do

        !$omp single
        
            deallocate(flist)

        !$omp end single





            



        !$omp end parallel
    
    end subroutine gen_fluid
    
end module internal