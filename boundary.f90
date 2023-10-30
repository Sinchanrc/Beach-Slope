module boundary

    use initialize
    use setup

    implicit none

    public
    
    contains

    subroutine gen_boundary()

        implicit none

        integer :: i,j,k,m

        ! Setting left hand vertical boundary positions
        !$omp parallel default(shared)

        !$omp do schedule(runtime) private(i,j) collapse(2) 

            do j =1,bl
                do i=bny,1,-1            

                    blist(i,j)%y=((bny-i)*brrealy*2.0)

                    blist(i,j)%x=(j-1)*brrealx*2

                    allocate(blist(i,j)%posshift(2))

                end do
            end do

        !$omp end do

        !Allocating index number to left hand vertical boundary particles 
        !$omp single
            icount=0
            count=0
            binmin=icount+1

            do j =1,bl
                do i=bny,1,-1
            
                if (j==bl) then

                if ((i<=(bny-bl+1))) then
                

                icount=icount+1
                blist(i,j)%tid=1
                blist(i,j)%pid=icount
                ! blist(i,j)%xnorm=.true.

                end if

                if ((i<bl)) then

                blist(i,j)%tid=2


                end if
                

                end if

                if((j/=bl).or.((j==bl).and.(i>(bny-bl+1))))then

                    icount=icount+1 !count=count+1
                    blist(i,j)%tid=2
                    blist(i,j)%pid=icount !icount
                    ! blist(i,j)%xnorm=.true.
                    ! if (i>=(bny-bl+1)) then
                    !     blist(i,j)%xnorm=.true.
                    !     blist(i,j)%ynorm=.true.
                    ! end if
                    
                end if

                end do
            end do

            do j =bl,1,-1
                do i=1,bny

                    matidct=matidct+1

                    if ((i<(bny-bl+1))) then

                    blist(i,j)%posshift(1)= 2*(((2*bl-1)*brrealx)-blist(i,j)%x)
                    blist(i,j)%posshift(2)= 0.0_dp

                    else 

                    blist(i,j)%posshift(1)= 2*(((2*bl-1)*brrealx)-blist(i,j)%x)
                    blist(i,j)%posshift(2)= 2*(((2*bl-1)*brrealy)-blist(i,j)%y)

                    end if

                    if (((j==bl).and.(i<=(bny-bl+1)))) then

                        
                        blist(i,j)%matid=matidct

                    elseif (((j==bl).and.(i>(bny-bl+1)))) then

                        allocate(blist(i,j)%wall(1))
                        blist(i,j)%wall(1)=blist((bny-bl+1),bl)%matid

                        blist(i,j)%matid=matidct

                        
                    end if

                    if (j/=bl) then

                        allocate(blist(i,j)%wall(1))
                        blist(i,j)%wall(1)=blist(i,(j+1))%matid

                        blist(i,j)%matid=matidct

                    end if

                end do 
            end do

            do i=1,(bl-1)
            if ((i<bl)) then

            allocate(blist(i,bl)%wall(1))
            blist(i,bl)%wall(1)=blist(bl,bl)%matid

            end if


            end do

            

        !$omp end single

        ! Moving boundary particles to cell boundary list
        !$omp do schedule(runtime) private(i,m,k) collapse(2) 
            do j=2,cellx-1
                do i=2,celly-1
                do m=1,bl
                    do k=1,bny

                    if ((blist(k,m)%x>=dpcell(i,j)%xleft) .and. &
                        (blist(k,m)%x<dpcell(i,j)%xright).and. &
                        (blist(k,m)%y>=dpcell(i,j)%ybot).and. &
                        (blist(k,m)%y<dpcell(i,j)%ytop)) then

                        if ((blist(k,m)%tid==2)) then
                            dpcell(i,j)%btot=dpcell(i,j)%btot+1
                            dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                            dpcell(i,j)%gcount=dpcell(i,j)%gcount+1
                            dpcell(i,j)%plist(dpcell(i,j)%btot)=blist(k,m)
                            if (.not.(blist(k,m)%ynorm)) then
                            dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=fmass
                            ! dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=4*brrealx*brrealy*rho!fmass
                            else
                                dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=fmass
                                ! dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=4*brrealx*brrealy*rho!fmass
                            end if
                            dpcell(i,j)%plist(dpcell(i,j)%btot)%density=rho

                        end if

                        if ((blist(k,m)%tid==1)) then

                            dpcell(i,j)%btot=dpcell(i,j)%btot+1
                            dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                            dpcell(i,j)%plist(dpcell(i,j)%btot)=blist(k,m)

                            dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=fmass
                            ! dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=4*brrealx*brrealy*rho!fmass

                            dpcell(i,j)%plist(dpcell(i,j)%btot)%density=rho


                        end if


                    end if

                    end do
                end do
                end do
            end do
        !$omp end do

        !Setting right hand vertical boundary positions
        !$omp do schedule(runtime) private(i,j) collapse(2) 

            do j =1,bl
                do i=bny,1,-1            
    
                blist(i,j)%x=blist(i,j)%x+(2*bl-1)*brrealx+brrealx+(bnx-2*bl-1)*brrealx*2+2*brrealx
                ! blist(i,j)%x=blist(i,j)%x+L+2*brrealx+brrealx+brrealx
                blist(i,j)%pid=0
                end do
            end do
    
        !$omp end do

        !Allocating index number to right hand vertical boundary particles 
        !$omp single 
            do j =1,bl
                do i=bny,1,-1
    
                if (((j==1))) then

                if (i<=(bny-bl+1)) then

                    icount=icount+1
                    blist(i,j)%tid=1
                    blist(i,j)%pid=icount

                    ! blist(i,j)%xnorm=.true.

                    end if

                if ((i<bl)) then

                blist(i,j)%tid=2


                end if

                end if
    
                if((j/=1).or.((j==1).and.(i>(bny-bl+1))))then
    
                    icount=icount+1 !count=count+1
                    blist(i,j)%tid=2
                    blist(i,j)%pid=icount !count

                    ! blist(i,j)%xnorm=.true.
                    ! if (i>=(bny-bl+1)) then
                    !     blist(i,j)%ynorm=.true.
                    ! end if
                    
                end if
    
                end do
            end do

            do j =1,bl
                do i=1,bny

                    if ((i<(bny-bl+1))) then

                    blist(i,j)%posshift(1)= 2*(-blist(i,j)%x+(blist(1,1)%x-brrealx))
                    blist(i,j)%posshift(2)= 0.0_dp

                    else 

                    blist(i,j)%posshift(1)= 2*(-blist(i,j)%x+(blist(1,1)%x-brrealx))
                    blist(i,j)%posshift(2)= 2*(((2*bl-1)*brrealy)-blist(i,j)%y)

                    end if

                    if (((j==1).and.(i<=(bny-bl+1)))) then

                        matidct=matidct+1
                        blist(i,j)%matid=matidct

                    elseif (((j==1).and.(i>(bny-bl+1)))) then

                        if (.not.(allocated(blist(i,j)%wall))) then
                        allocate(blist(i,j)%wall(1))
                        end if
                        blist(i,j)%wall(1)=blist((bny-bl+1),1)%matid
                        matidct=matidct+1
                        blist(i,j)%matid=matidct
                        
                    end if

                    if (j/=1) then

                        if (.not.(allocated(blist(i,j)%wall))) then
                        allocate(blist(i,j)%wall(1))
                        end if
                        blist(i,j)%wall(1)=blist(i,(j-1))%matid
                        matidct=matidct+1
                        blist(i,j)%matid=matidct

                    end if

                end do 
            end do

            do i=1,(bl-1)
            if ((i<bl)) then

            if (.not.(allocated(blist(i,1)%wall))) then
                        allocate(blist(i,1)%wall(1))
            end if
            blist(i,1)%wall(1)=blist(bl,1)%matid

            end if


            end do

            
        
        !$omp end single

        ! Moving boundary particles to cell boundary list
        !$omp do schedule(runtime) private(i,m,k) collapse(2) 
            do j=2,cellx-1
                do i=2,celly-1
                do m=1,bl
                    do k=1,bny

                    if ((blist(k,m)%x>=dpcell(i,j)%xleft) .and. &
                        (blist(k,m)%x<dpcell(i,j)%xright).and. &
                        (blist(k,m)%y>=dpcell(i,j)%ybot).and. &
                        (blist(k,m)%y<dpcell(i,j)%ytop)) then

                        if ((blist(k,m)%tid==2)) then
                            dpcell(i,j)%btot=dpcell(i,j)%btot+1
                            dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                            dpcell(i,j)%gcount=dpcell(i,j)%gcount+1
                            dpcell(i,j)%plist(dpcell(i,j)%btot)=blist(k,m)
                            if (.not.(blist(k,m)%ynorm)) then
                                dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=fmass
                            ! dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=4*brrealx*brrealy*rho!fmass
                            else
                                dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=fmass
                                ! dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=4*brrealx*brrealy*rho!fmass
                            end if
                            dpcell(i,j)%plist(dpcell(i,j)%btot)%density=rho


                        end if

                        if ((blist(k,m)%tid==1)) then
                            dpcell(i,j)%btot=dpcell(i,j)%btot+1
                            dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                            dpcell(i,j)%plist(dpcell(i,j)%btot)=blist(k,m)

                            dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=fmass
                            ! dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=4*brrealx*brrealy*rho!fmass

                            dpcell(i,j)%plist(dpcell(i,j)%btot)%density=rho

                        end if


                    end if

                    end do
                end do
                end do
            end do
        !$omp end do

        !$omp single
            deallocate(blist)
            allocate(blist(bl,(bnx-2*bl)))
        !$omp end single

        !Setting bottom boundary positions
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j =1,bnx-2*bl
                do i=1,bl           
    
                blist(i,j)%y=(bl-i)*2*brrealy
                blist(i,j)%x=(2*bl-1)*brrealx+brrealx+(j-1)*brrealx*2
                allocate(blist(i,j)%posshift(2))
    
                end do
            end do
        !$omp end do
    
        !Allocating index number to bottom boundary particles
        !$omp single 
    
            do j =1,bnx-2*bl
                do i=1,bl
                    if (i==1) then

                        icount=icount+1
                        blist(i,j)%tid=1
                        blist(i,j)%pid=icount

                    end if
                    ! blist(i,j)%ynorm=.true.
                    ! blist(i,j)%xnorm=.false.
                    if((i/=1))then

                        icount=icount+1 !count=count+1
                        blist(i,j)%tid=2
                        blist(i,j)%pid=icount !count

                    end if
                end do
            end do
            binmax=icount

            do j =1,bnx-2*bl
                do i=1,bl

                    blist(i,j)%posshift(1)= 0.0_dp

                    blist(i,j)%posshift(2)= 2*(((2*bl-1)*brrealy)-blist(i,j)%y)

                    if (i==1) then

                        matidct=matidct+1
                        blist(i,j)%matid=matidct

                    else

                        allocate(blist(i,j)%wall(1))
                        blist(i,j)%wall(1)=blist(1,j)%matid
                        matidct=matidct+1
                        blist(i,j)%matid=matidct
                        
                    end if

                end do 
            end do
    
        !$omp end single
    
        ! Moving boundary particles to cell boundary list
        !$omp do schedule(runtime) private(i,m,k) collapse(2) 
            do j=2,cellx-1
                do i=2,celly-1
                do m=1,bnx-2*bl
                    do k=1,bl

                    if ((blist(k,m)%x>=dpcell(i,j)%xleft) .and. &
                        (blist(k,m)%x<dpcell(i,j)%xright).and. &
                        (blist(k,m)%y>=dpcell(i,j)%ybot).and. &
                        (blist(k,m)%y<dpcell(i,j)%ytop)) then

                        if ((blist(k,m)%tid==2)) then
                            dpcell(i,j)%btot=dpcell(i,j)%btot+1
                            dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                            dpcell(i,j)%gcount=dpcell(i,j)%gcount+1
                            dpcell(i,j)%plist(dpcell(i,j)%btot)=blist(k,m)
                            if (.not.(blist(k,m)%ynorm)) then
                                dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=fmass
                            ! dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=4*brrealx*brrealy*rho!fmass
                            else
                                dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=fmass
                                ! dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=4*brrealx*brrealy*rho!fmass
                            end if
                            dpcell(i,j)%plist(dpcell(i,j)%btot)%density=rho

                        end if

                        if ((blist(k,m)%tid==1)) then
                            dpcell(i,j)%btot=dpcell(i,j)%btot+1
                            dpcell(i,j)%ptot=dpcell(i,j)%ptot+1
                            dpcell(i,j)%plist(dpcell(i,j)%btot)=blist(k,m)

                            dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=fmass
                            ! dpcell(i,j)%plist(dpcell(i,j)%btot)%mass=4*brrealx*brrealy*rho!fmass

                            dpcell(i,j)%plist(dpcell(i,j)%btot)%density=rho

                        end if
                        
                    end if

                    end do
                end do
                end do
            end do
        !$omp end do
        
        !$omp single
            deallocate(blist)
        !$omp end single

                    !$omp single
        
            allocate(blist(spy,spx))
        !$omp end single

        !Setting up soil particle positions and pressure
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j =1,spx
            do i=spy,1,-1


            blist(i,j)%y=(spy-i)*2*solidy+solidy+((brrealy)*((2*bl)-1))-6*solidy
            blist(i,j)%x=((brrealx)*((2*bl)-1))+(j-1)*2*solidx+brrealx-6*solidx

            end do
            end do
        !$omp end do

        !Setting particle id for fluid particles
        !$omp single 
            count=0

            do j =1,spx
                do i=spy,1,-1
    
                    count=count+1
                    blist(i,j)%pid=count
                    blist(i,j)%tid=4
                    ! if(((xl-blist(i,j)%x)*(yu-blist(i,j)%y)-&
                    ! (yl-blist(i,j)%y)*(xu-blist(i,j)%x))<0.0) then
                    !     blist(i,j)%pid=0
                    ! end if

                    if(((blist(i,j)%y-yl)-line_grad*(blist(i,j)%x-xl))>0.0) then
                        blist(i,j)%pid=0
                    end if
                
                
                end do
            end do
        
        !$omp end single

        !Distributing porous particles to cells
        !$omp do schedule(runtime) private(i,j) collapse(2) 
            do j=2,cellx-1
                do i=2,celly-1
                do l1=1,spx
                    do k=1,spy

                    if ((blist(k,l1)%x>=dpcell(i,j)%xleft) .and. &
                        (blist(k,l1)%x<dpcell(i,j)%xright).and. &
                        (blist(k,l1)%y>=dpcell(i,j)%ybot).and. &
                        (blist(k,l1)%y<dpcell(i,j)%ytop).and.(blist(k,l1)%pid/=0))then
                        
                        

                        if (.not.(allocated(dpcell(i,j)%porlist))) then
                            allocate(dpcell(i,j)%porlist(fplistmax))
                        end if

                        dpcell(i,j)%porct=dpcell(i,j)%porct+1
                        dpcell(i,j)%porlist(dpcell(i,j)%porct)=blist(k,l1)
                        dpcell(i,j)%porlist(dpcell(i,j)%porct)%mass=soill*soilh*bulkden/(spx*spy)
                        dpcell(i,j)%porlist(dpcell(i,j)%porct)%density=bulkden

                    end if

                    end do
                end do
                end do
            end do
        !$omp end do

        !$omp single
        
            deallocate(blist)
        !$omp end single



        !$omp end parallel

        return
    end subroutine gen_boundary
    
end module boundary