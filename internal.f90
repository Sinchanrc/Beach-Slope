module internal

    use initialize
    use particle

    implicit none

    integer :: i,j,k,m
    
    contains

    subroutine gen_fluid()

        implicit none

        !$omp parallel default(shared)
        !Setting up fluid particle positions and pressure
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j =1,fpx
            do i=fpy,1,-1
                flist(i,j)%y=((brrealy)*((2*bl)-1))+(fpy-i)*2*prrealy/sqrt(por)+prrealy/sqrt(por)
                if (mod(i,2)==0) then
                flist(i,j)%x=((brrealx)*((2*bl)-1))+(j-1)*2*prrealx/sqrt(por)+2*prrealx/sqrt(por)
                else
                flist(i,j)%x=((brrealx)*((2*bl)-1))+(j-1)*2*prrealx/sqrt(por)+prrealx/sqrt(por)
                end if
            flist(i,j)%vx=0.0_dp
            flist(i,j)%vy=0.0_dp

            if (i==1) then
                flist(i,j)%pressure=0.0_dp

            else
                flist(i,j)%pressure=(-flist(i,j)%y+((brrealy*distfac)*((2*bl)-1))+wc+prrealy/sqrt(por))*rho*abs(g)
            
            end if

            end do
            end do
        !$omp end do

        !Setting particle id for fluid particles
        !$omp single 
            count=0
            finmin=count+1
            do j =1,fpx
                do i=fpy,1,-1
    
                if ((mod(i,2)/=0)) then
                    count=count+1
                    flist(i,j)%pid=count
                    flist(i,j)%tid=3
                else 
                    if (j/=fpx) then
                        count=count+1
                    flist(i,j)%pid=count
                    flist(i,j)%tid=3
                    else 

                    flist(i,j)%tid=0

                    end if

                end if

                if(((flist(i,j)%y-yl)-line_grad*(flist(i,j)%x-xl))>0.0) then
                    flist(i,j)%tid=0
                end if
                
                end do
            end do
            ! finmax=count

            ! fmass=fmass*fpx*fpy/finmax
        
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
                        ! if (dpcell(i,j)%plist(dpcell(i,j)%ptot)%x<=((wl/10.0_dp)+(2*bl)*brrealx)) then
                        ! dpcell(i,j)%plist(dpcell(i,j)%ptot)%mass=fmass*rhomax/rhomin
                        ! dpcell(i,j)%plist(dpcell(i,j)%ptot)%density=rhomax
                        ! dpcell(i,j)%plist(dpcell(i,j)%ptot)%ovol=dpcell(i,j)%plist(dpcell(i,j)%ptot)%mass/&
                        !                                             dpcell(i,j)%plist(dpcell(i,j)%ptot)%density
                        ! dpcell(i,j)%plist(dpcell(i,j)%ptot)%con=0.50_dp
                        ! else
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%mass=fmass
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%density=rho

                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%oden=dpcell(i,j)%plist&
                        (dpcell(i,j)%ptot)%density
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%ovol=fmass/rho
                        ! dpcell(i,j)%plist(dpcell(i,j)%ptot)%con=-0.50_dp
                        ! end if
                    end if
                    end do
                end do
                end do
            end do
        !$omp end do

        !$omp single
        
            deallocate(flist)
            fpy=floor(real(coastal_ht,dp)/(2*real(prrealy,dp)))+1
            fpx=floor(real(2.0_dp,dp)/(2*real(prrealx,dp)))!+1
            allocate(flist(fpy,fpx))
        !$omp end single

!!!!!!!!!!!!!!!!!!! Setting up porous media(non-deform)!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !Setting up fluid particle positions and pressure
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j =1,fpx
                do i=fpy,1,-1
    
    
                flist(i,j)%y=((brrealy)*((2*bl)-1))+(fpy-i)*2*prrealy+prrealy+set_ht
                if (mod(i,2)==0) then
                flist(i,j)%x=((brrealx)*((2*bl)-1))+(j-1)*2*prrealx+2*prrealx
                else
                flist(i,j)%x=((brrealx)*((2*bl)-1))+(j-1)*2*prrealx+prrealx
                end if
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

                    if ((mod(i,2)/=0)) then
                        count=count+1
                        flist(i,j)%pid=count
                        flist(i,j)%tid=3
                    else 
                        if (j/=fpx) then
                            count=count+1
                        flist(i,j)%pid=count 
                        flist(i,j)%tid=3
                        else 
    
                        flist(i,j)%tid=0
    
                        end if
    
                    end if
    
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
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%mass=fmass*1.05
                        dpcell(i,j)%plist(dpcell(i,j)%ptot)%density=rho*1.05
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