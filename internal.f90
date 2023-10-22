module internal

    use initialize
    use particle

    implicit none
    
    contains

    subroutine gen_fluid()

        implicit none

        !$omp parallel default(shared)
        !Setting up fluid particle positions and pressure
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j =1,fpx
            do i=fpy,1,-1
            flist(i,j)%y=((brrealy)*((2*bl)-1))+(fpy-i)*2*prrealy+prrealy
            ! if (mod(i,2)==0) then
            ! flist(i,j)%x=((brrealx)*((2*bl)-1))+(j-1)*2*prrealx+2*prrealx
            ! else
            flist(i,j)%x=((brrealx)*((2*bl)-1))+(j-1)*2*prrealx+prrealx
            ! end if
            flist(i,j)%vx=0.0_dp
            flist(i,j)%vy=0.0_dp

            end do
            end do
        !$omp end do

        !Setting particle id for fluid particles
        !$omp single 
            count=0
            finmin=count+1
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
                        ! dpcell(i,j)%plist(dpcell(i,j)%ptot)%ovol=dpcell(i,j)%plist(dpcell(i,j)%ptot)%mass/&
                        !                                             dpcell(i,j)%plist(dpcell(i,j)%ptot)%density
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
        !$omp end single
        !$omp end parallel
    
    end subroutine gen_fluid
    
end module internal