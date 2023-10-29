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
    use mls
    use pprobe

    implicit none

    call init
    call gen_boundary
    call gen_fluid    
    call memalloc

    ! ploc(1,1)=5.366*wc+(2*bl-1)*brrealx!+2*prrealx
    ! ploc(2,1)=0.19*wc+(2*bl-1)*brrealy
    ! call probesetup(ploc)
    ! call print_probe

    call print_fluid
    call print_fixbd
    call print_ghostbd
    call matrixid
    iter=iter+1
    ! write(*,*)brrealx

    ! do while((t*sqrt(abs(g)/wc))<time)
    do while(iter<411)

        told=t
        t=t+dt
        
        ! if (dtsol>dt) then
        ! !$omp parallel default(shared)
        ! call scalart
        ! call scalarupdate(dt)
        ! !$omp end parallel
        ! else
        ! do i=1,solsteps
        ! !$omp parallel default(shared)
        ! call scalart
        ! call scalarupdate(dtsol)
        ! !$omp end parallel
        ! end do
        ! end if

        !$omp parallel default(shared)

        call projection 
        call cellshift
        call neighbour
        ! call effpor
        call freesurf
        ! call compcorr(3,1)
        call eddyvis
        ! call ghost_en
        call int_vel
        !$omp end parallel

        call ppesolve
    
        !$omp parallel default(shared)
        call comp_vel        
        call comp_pos
        call cellshift
        call neighbour
        ! call effpor
        call freesurf
        ! call compcorr(3,1)
        ! call comp_ghost
        ! call boun_vel 
        call opt2_shift
        ! call massupdate
        call timestep
        ! call eddyvis 
        !$omp end parallel

        ! call implicit_shift()

        ! !$omp parallel default(shared)
        ! ! call massupdate
        ! call timestep
        ! call eddyvis
        ! !$omp end parallel

        ! if (((told*sqrt(abs(g)/wc))<iter*displaytime).and. &
        ! ((t*sqrt(abs(g)/wc))>=iter*displaytime)) then
        
        ! call probevalue
        call print_fluid
        iter=iter+1
        ! end if

    end do
    call print_free
    
end program dam_break

