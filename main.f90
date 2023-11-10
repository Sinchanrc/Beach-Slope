program dam_break

    use initialize
    use setup
    use boundary
    use solver
    use output
    use internal
    use memory
    use search
    use isph
    use integrator
    use turbulence
    use part_shift
    use gradient
    use scalar
    use porous

    implicit none

    integer :: j1,i1

    call init
    call gen_boundary
    call gen_fluid    
    call memalloc

    call print_fluid
    call print_porous
    call print_fixbd
    call print_ghostbd
    call matrixid
    iter=iter+1

    !     do while(iter<51)

    !         told=t
    !         t=t+dt
            
    !         ! if (dtsol>dt) then
    !         ! !$omp parallel default(shared)
    !         ! call scalart
    !         ! call scalarupdate(dt)
    !         ! !$omp end parallel
    !         ! else
    !         ! do i=1,solsteps
    !         ! !$omp parallel default(shared)
    !         ! call scalart
    !         ! call scalarupdate(dtsol)
    !         ! !$omp end parallel
    !         ! end do
    !         ! end if

    !         !$omp parallel default(shared)

    !         call projection 
    !         call cellshift
    !         call neighbour
    !         call effpor
    !         call freesurf
    !         ! call compcorr(3,1)
    !         call eddyvis
    !         ! call ghost_en
    !         call int_vel
    !         !$omp end parallel
    !         call resetid
    !         call ppesolve
        
    !         !$omp parallel default(shared)
    !         call comp_vel        
    !         call comp_pos
    !         call cellshift
    !         call neighbour
    !         call effpor
    !         call freesurf
    !         ! call compcorr(3,1)
    !         ! call comp_ghost
    !         ! call boun_vel 
    !         call opt2_shift
    !         ! call massupdate
    !         call timestep
    !         ! call eddyvis 
    !         !$omp end parallel

    !         ! call implicit_shift()

    !         ! !$omp parallel default(shared)
    !         ! ! call massupdate
    !         ! call timestep
    !         ! call eddyvis
    !         ! !$omp end parallel

    !         ! if (((told*sqrt(abs(g)/wc))<iter*displaytime).and. &
    !         ! ((t*sqrt(abs(g)/wc))>=iter*displaytime)) then
            
    !         ! call probevalue
    !         ! call print_fluid
    !         iter=iter+1
    !         ! end if

    !     end do

    !     iter=1

    !     do j1=sx,ex 
    !         do i1=sy,ey
    !         if (dpcell(i1,j1)%ptot/=0) then
    !             do cout=1,dpcell(i1,j1)%ptot

    !                 if ((dpcell(i1,j1)%plist(cout)%tid==3)) then

    !                     dpcell(i1,j1)%plist(cout)%vx=0.0_dp
    !                     dpcell(i1,j1)%plist(cout)%vy=-2.5_dp/(3600*24)

    !                     if(((dpcell(i1,j1)%plist(cout)%y-yl-prrealy)-line_grad* &
    !                     (dpcell(i1,j1)%plist(cout)%x-xl))>0.0) then
    !                         dpcell(i1,j1)%plist(cout)%vx=-2.5_dp/(3600*24*por)
    !                     end if

    !                 end if


    !             end do
    !         end if
    !         end do
    !     end do

    ! do while(iter<4001)

    !     told=t
    !     t=t+dt
        
    !     ! if (dtsol>dt) then
    !     ! !$omp parallel default(shared)
    !     ! call scalart
    !     ! call scalarupdate(dt)
    !     ! !$omp end parallel
    !     ! else
    !     ! do i=1,solsteps
    !     ! !$omp parallel default(shared)
    !     ! call scalart
    !     ! call scalarupdate(dtsol)
    !     ! !$omp end parallel
    !     ! end do
    !     ! end if

    !     !$omp parallel default(shared)

    !     call projection 
    !     call cellshift
    !     call neighbour
    !     call effpor
    !     call freesurf
    !     ! call compcorr(3,1)
    !     call eddyvis
    !     ! call ghost_en
    !     call int_vel
    !     !$omp end parallel
    !     call resetid
    !     call ppesolve
    
    !     !$omp parallel default(shared)
    !     call comp_vel        
    !     call comp_pos
    !     call cellshift
    !     call neighbour
    !     call effpor
    !     call freesurf
    !     ! call compcorr(3,1)
    !     ! call comp_ghost
    !     ! call boun_vel 
    !     call opt2_shift
    !     ! call massupdate
    !     call timestep
    !     ! call eddyvis 
    !     !$omp end parallel

    !     ! call implicit_shift()

    !     ! !$omp parallel default(shared)
    !     ! ! call massupdate
    !     ! call timestep
    !     ! call eddyvis
    !     ! !$omp end parallel

    !     ! if (((told*sqrt(abs(g)/wc))<iter*displaytime).and. &
    !     ! ((t*sqrt(abs(g)/wc))>=iter*displaytime)) then
        
    !     ! call probevalue
    !     call print_fluid
    !     iter=iter+1
    !     ! end if

    ! end do

    
end program dam_break

