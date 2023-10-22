module porous

    use particle
    use kernel
    use domain
    use initialize

    implicit none

    contains

    subroutine porsearch(list,n,cell1,cell2,cell3,cell4,smln,solidfrac)
        implicit none
        type(verlet),intent(inout) :: list
        type(cell),intent(inout) :: cell1,cell2,cell3,cell4
        type(cell) :: blockset(4)
        real(dp),intent(in) :: smln
        integer,intent(in) :: n
        integer ::counter1,counter2,intct
        real(dp) :: radial
        real(dp),intent(out) :: solidfrac


        blockset(1)=cell1
        blockset(2)=cell2
        blockset(3)=cell3
        blockset(4)=cell4

        solidfrac=0.0_dp

        ! intct=0
        ! list%count=0
        ! cell1%plist(n)%mirror=.false.

            do counter1=1,4
            if (blockset(counter1)%porct/=0) then
                do counter2=1,blockset(counter1)%porct
                radial=dist(blockset(1)%plist(n),blockset(counter1)%porlist(counter2))
                if((radial<(2*smln)))then
                    ! if ((cell1%plist(n)%tid/=blockset(counter1)%plist(counter2)%tid) &
                    ! .or.((cell1%plist(n)%tid==blockset(counter1)%plist(counter2)%tid) &
                    ! .and.(cell1%plist(n)%pid/=blockset(counter1)%plist(counter2)%pid))) then

                    solidfrac=solidfrac+(1.0_dp-por)*Wab(radial,h1)* &
                                        blockset(counter1)%porlist(counter2)%mass &
                                        /blockset(counter1)%porlist(counter2)%density
                    ! end if

                end if
                end do
            end if
            end do

    end subroutine porsearch

    subroutine effpor
    
        implicit none

        real(dp) :: solidfrac,lamk=150.0_dp,d50=0.01590_dp,Fch,kper

        !$omp do schedule (runtime) private(m,term1,term2,k,i,j,solidfrac,Fch,kper) collapse(2)
        do j=sx,ex
            do i=sy,ey 

            do k=1,dpcell(i,j)%ptot

                solidfrac=0.0_dp
                dpcell(i,j)%pplist(k)%resistx=0.0_dp
                dpcell(i,j)%pplist(k)%resisty=0.0_dp
                dpcell(i,j)%pplist(k)%transporx=0.0_dp
                dpcell(i,j)%pplist(k)%transpory=0.0_dp

                if (dpcell(i,j)%plist(k)%tid==3) then

                dpcell(i,j)%pplist(k)%porosity=1.0_dp
                

                if ((dpcell(i,j)%plist(k)%x<=(((brrealx)*((2*bl)-1))+0.30_dp+brrealx+2*h1)).or. &
                (dpcell(i,j)%plist(k)%x>=(((brrealx)*((2*bl)-1)) &
                +0.30_dp+brrealx+(spx-1)*2*solidx))) then

                dpcell(i,j)%pplist(k)%inpore=.false.

                else 

                dpcell(i,j)%pplist(k)%inpore=.true.

                end if

                if(dpcell(i,j)%plist(k)%x<=(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y<=(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    call porsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i,j-1), &
                    dpcell(i+1,j),dpcell(i+1,j-1),h1,solidfrac)

                else if(dpcell(i,j)%plist(k)%x>(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y<=(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    call porsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i,j+1), &
                    dpcell(i+1,j),dpcell(i+1,j+1),h1,solidfrac)

                else if(dpcell(i,j)%plist(k)%x>(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y>(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    call porsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i-1,j+1), &
                    dpcell(i-1,j),dpcell(i,j+1),h1,solidfrac)

                else if(dpcell(i,j)%plist(k)%x<=(dpcell(i,j)%xleft+ &
                abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                dpcell(i,j)%plist(k)%y>(dpcell(i,j)%ybot+ &
                abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    call porsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i-1,j), &
                    dpcell(i-1,j-1),dpcell(i,j-1),h1,solidfrac)

                else if(dpcell(i,j)%plist(k)%x==(dpcell(i,j)%xleft+ &
                    abs((dpcell(i,j)%xleft-dpcell(i,j)%xright)/2)).and. &
                    dpcell(i,j)%plist(k)%y==(dpcell(i,j)%ybot+ &
                    abs((dpcell(i,j)%ytop-dpcell(i,j)%ybot)/2))) then
                    call porsearch(dpcell(i,j)%list(k),k,dpcell(i,j),dpcell(i-1,j), &
                    dpcell(i-1,j-1),dpcell(i,j-1),h1,solidfrac)
                end if


                ! if (dpcell(i,j)%list(k)%count/=0) then

                ! do m=1,dpcell(i,j)%list(k)%count
                !     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                !     y=>dpcell(i,j)%list(k)%interlist(2,m), &
                !     pp=>dpcell(i,j)%list(k)%interlist(3,m))

                !     if (dpcell(y,x)%plist(pp)%tid==4) then

                !         solidfrac=solidfrac+(1.0_dp-por)*Wab(dpcell(i,j)%list(k)%dist(m),h1)* &
                !                         dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density

                !         dpcell(i,j)%pplist(k)%transporx=dpcell(i,j)%pplist(k)%transporx+((1.0_dp-por)* &
                !                         Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)* &
                !                         (dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density))

                !         dpcell(i,j)%pplist(k)%transpory=dpcell(i,j)%pplist(k)%transpory+((1.0_dp-por)* &
                !                         Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)* &
                !                         (dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density))

                        

                !     end if


                !     end associate
                ! end do

                ! end if

                if (solidfrac>(1.0_dp-por)) then

                        solidfrac=1.0_dp-por

                end if

                dpcell(i,j)%pplist(k)%porosity=1.0_dp-solidfrac
                dpcell(i,j)%plist(k)%density=dpcell(i,j)%pplist(k)%porosity*rho!fmass/(4*prrealx*prrealy)

                kper=(dpcell(i,j)%pplist(k)%porosity**3)*(d50**2)/(lamk*(1.0_dp-dpcell(i,j)%pplist(k)%porosity)**2)

                Fch=1.750_dp/(sqrt(150*dpcell(i,j)%pplist(k)%porosity**3))

                dpcell(i,j)%pplist(k)%resistx=-((mu/(fmass/(4*prrealx*prrealy)))*dpcell(i,j)%plist(k)%vx/kper) &
                -(Fch*dpcell(i,j)%plist(k)%vx*sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2))/sqrt(kper)

                dpcell(i,j)%pplist(k)%resisty=-((mu/(fmass/(4*prrealx*prrealy)))*dpcell(i,j)%plist(k)%vy/kper) &
                -(Fch*dpcell(i,j)%plist(k)%vy*sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2))/sqrt(kper)



                ! dpcell(i,j)%pplist(k)%resistx=-((mu/rho)* &
                ! ((1.0_dp-dpcell(i,j)%pplist(k)%porosity)**2)*lamk &
                ! *dpcell(i,j)%plist(k)%vx/((dpcell(i,j)%pplist(k)%porosity**2)*(d50**2))) &
                ! -(1.750_dp*dpcell(i,j)%plist(k)%vx*sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2) &
                ! *(1.0_dp-dpcell(i,j)%pplist(k)%porosity)*sqrt(lamk)/(sqrt(150*dpcell(i,j)%pplist(k)%porosity**2) &
                ! *d50))

                ! dpcell(i,j)%pplist(k)%resisty=-((mu/rho)* &
                ! ((1.0_dp-dpcell(i,j)%pplist(k)%porosity)**2)*lamk &
                ! *dpcell(i,j)%plist(k)%vy/((dpcell(i,j)%pplist(k)%porosity**2)*(d50**2))) &
                ! -(1.750_dp*dpcell(i,j)%plist(k)%vy*sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2) &
                ! *(1.0_dp-dpcell(i,j)%pplist(k)%porosity)*sqrt(lamk)/(sqrt(150*dpcell(i,j)%pplist(k)%porosity**2) &
                ! *d50))

                ! dpcell(i,j)%pplist(k)%resistx=-((mu/rho)* &
                ! ((1.0_dp-dpcell(i,j)%pplist(k)%porosity)**2)*lamk &
                ! *dpcell(i,j)%plist(k)%vx/((dpcell(i,j)%pplist(k)%porosity)*(d50**2))) &
                ! -(1.750_dp*dpcell(i,j)%plist(k)%vx*sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2) &
                ! *(1.0_dp-dpcell(i,j)%pplist(k)%porosity)*sqrt(lamk)/(sqrt(150*dpcell(i,j)%pplist(k)%porosity) &
                ! *d50))

                ! dpcell(i,j)%pplist(k)%resisty=-((mu/rho)* &
                ! ((1.0_dp-dpcell(i,j)%pplist(k)%porosity)**2)*lamk &
                ! *dpcell(i,j)%plist(k)%vy/((dpcell(i,j)%pplist(k)%porosity)*(d50**2))) &
                ! -(1.750_dp*dpcell(i,j)%plist(k)%vy*sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2) &
                ! *(1.0_dp-dpcell(i,j)%pplist(k)%porosity)*sqrt(lamk)/(sqrt(150*dpcell(i,j)%pplist(k)%porosity) &
                ! *d50))

                end if

            end do

            end do
        end do
        !$omp end do
        
    end subroutine effpor
    

end module porous
