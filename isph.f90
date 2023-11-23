module isph

    use particle
    use initialize
    use domain
    use kernel
    use functions
    use solver
    use output


    implicit none

    contains

    subroutine projection
        implicit none

        integer :: i,j,k,m
        
        ! Intermediate fluid position calc
        !$omp do schedule(runtime) collapse(2) private(i,j,k)
            do j=sx,ex
                do i=sy,ey
                if (dpcell(i,j)%ptot/=0) then
                
                    do k=1,dpcell(i,j)%ptot

                    if (dpcell(i,j)%plist(k)%tid==3) then
                    dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%x
                    dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+(dt* &
                    dpcell(i,j)%plist(k)%vx/dpcell(i,j)%pplist(k)%porosity)
                    dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%y
                    dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+(dt* &
                    dpcell(i,j)%plist(k)%vy/dpcell(i,j)%pplist(k)%porosity)
                    end if
                    end do

                end if
                end do
            end do
        !$omp end do        
        

        
    end subroutine projection
    
    subroutine freesurf
        implicit none 

        real(dp) :: heff,t1,t2
        integer :: i,j,k,m


        ! Free surface(dpcell(y,x)%plist(pp)%tid==3) then
        !$omp do private(m,t1,t2,i,k,j,heff) schedule (runtime) collapse(2)    
        do j=sx,ex
            do i=sy,ey
            if(dpcell(i,j)%ptot/=0) then
                do k=1,dpcell(i,j)%ptot
                ! Free Surface Identification
                if ((dpcell(i,j)%plist(k)%tid==3).and. &
                (.not.(dpcell(i,j)%plist(k)%buffer))) then
                ! dpcell(i,j)%plist(k)%free=.false.
                ! dpcell(i,j)%pplist(k)%gradvx=0.0_dp
                ! t1=0.0_dp
                ! t2=0.0_dp

                ! if (dpcell(i,j)%pplist(k)%inpore) then

                heff=h1!(max(prrealx,prrealy)*hfac)/(2.0_dp*sqrt(dpcell(i,j)%pplist(k)%porosity))  !*sqrt(por)

                ! else 

                ! heff=(max(prrealx,prrealy)*hfac)/(2.0_dp)

                ! end if

                ! do m=1,dpcell(i,j)%list(k)%count
                !     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                !     y=>dpcell(i,j)%list(k)%interlist(2,m), &
                !     pp=>dpcell(i,j)%list(k)%interlist(3,m))

                !     if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                !     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                !     t1=t1+((dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%plist(pp)%x- &
                !     dpcell(i,j)%plist(k)%x)*Wabx(dpcell(y,x)%plist(pp),&
                !     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))/&
                !     dpcell(y,x)%plist(pp)%density)

                !     t2=t2+((dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%plist(pp)%y- &
                !     dpcell(i,j)%plist(k)%y)*Waby(dpcell(y,x)%plist(pp),&
                !     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))/&
                !     dpcell(y,x)%plist(pp)%density)

                !     end if

                !     end associate
                ! end do

                ! dpcell(i,j)%pplist(k)%gradvx=t1+t2
                if(dpcell(i,j)%pplist(k)%gradvx<=(lamfs*maxdivr)) then
                    ! dpcell(i,j)%plist(k)%free=.true.
                
                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if ((dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then
    
                        dpcell(y,x)%plist(pp)%vicinity=.true.

                        end if
    
                        end associate
                    end do

                end if

                ! if (dpcell(i,j)%pplist(k)%gradvx<=pll) then

                !     dpcell(i,j)%pplist(k)%lamp=0.50_dp

                ! elseif((pll<dpcell(i,j)%pplist(k)%gradvx).and.(dpcell(i,j)%pplist(k)%gradvx<pul)) then

                ! dpcell(i,j)%pplist(k)%lamp=0.50_dp*(1.0_dp-cos(22.0_dp*(dpcell(i,j)%pplist(k)%gradvx-pll)/ &
                !                 ((pul-pll)*7.0_dp)))
                ! else 

                !     dpcell(i,j)%pplist(k)%lamp =1.0_dp

                ! end if


                end if

                end do
            end if
            end do
        end do
        !$omp end do
        
    end subroutine freesurf

    subroutine resetid
        implicit none 

        integer :: mat_count,t1,i,j

        mat_count=0

        do j=1,cellx 
            do i=1,celly
            if (dpcell(i,j)%ptot/=0) then
                do t1=1,dpcell(i,j)%ptot

                    mat_count=mat_count+1
                    dpcell(i,j)%plist(t1)%matid=mat_count

                end do
            end if
            end do
        end do
    
        
    end subroutine resetid

    subroutine ppesolve

        implicit none

        integer :: i,j,k,m

        real(dp) :: lamp,lamp2,t1,t2,rho_i,rho_j,pvol,p_dist,rho_ij
        !$omp parallel do schedule(runtime) default(shared)
        do i=1,finmax 
            fmatrix(i)%sz=0
            fmatrix(i)%val(:)=0.0_dp
            fmatrix(i)%col(:)=i
            fvec(i)=0.0_dp
        end do
        !$omp end parallel do
        
        ! Preparing the coff matrix for fluid part in CSR
        !$omp parallel do schedule(runtime) default(shared) &
        !$omp private(m,t1,t2,k,i,j,lamp,lamp2,rho_i,rho_j,pvol,p_dist,rho_ij) collapse(2)
        do j=sx,ex
            do i=sy,ey
            ! if (dpcell(i,j)%ptot/=0) then

                do k=1,dpcell(i,j)%ptot
                        
                ! if (dpcell(i,j)%plist(k)%tid/=4) then
                associate(pos=>dpcell(i,j)%plist(k)%matid,num=>dpcell(i,j)%list(k)%count)
                    t1=0.0_dp
                    t2=0.0_dp
                    fmatrix(pos)%sz=num+1
                    fmatrix(pos)%val(:)=0.0_dp
                    fmatrix(pos)%col(:)=pos
                    fvec(pos)=0.0_dp
                    num2=dpcell(i,j)%list(k)%count
                    rho_i=dpcell(i,j)%plist(k)%density*dpcell(i,j)%pplist(k)%porosity**(-1)
                    if((num/=0)) then   !.or.(dpcell(i,j)%plist(k)%free/=1)
                    do m=1,num2
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                            ! if (dpcell(y,x)%plist(pp)%tid/=4) then

                            if (dpcell(i,j)%plist(k)%free) then

                                lamp=(1.0_dp+dpcell(y,x)%plist(pp)%density/(dpcell(y,x)%pplist(pp)%porosity*rho*1.0))**(-1)

                            ! elseif((pll<dpcell(i,j)%pplist(k)%gradvx).and.(dpcell(i,j)%pplist(k)%gradvx<pul)) then

                            ! lamp=((1.0_dp+dpcell(y,x)%plist(pp)%density/ &
                            ! (dpcell(y,x)%pplist(pp)%porosity*rho*1.0))**(-1)) &
                            ! *(1.0_dp-cos(22.0_dp*(dpcell(i,j)%pplist(k)%gradvx-pll)/ &
                            !                 ((pul-pll)*7.0_dp)))
                            else 

                                lamp =1.0_dp

                            end if

                            ! rho_i=dpcell(i,j)%plist(k)%density*dpcell(i,j)%pplist(k)%porosity**(-2)
                            rho_j=dpcell(y,x)%plist(pp)%density*dpcell(y,x)%pplist(pp)%porosity**(-1)
                            pvol=dpcell(y,x)%plist(pp)%mass*dpcell(y,x)%plist(pp)%density**(-1)
                            p_dist=((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)**(-1)

                            rho_ij=2*dpcell(y,x)%plist(pp)%density*dpcell(i,j)%plist(k)%density* &
                                (((dpcell(y,x)%pplist(pp)%porosity)*(dpcell(i,j)%pplist(k)%porosity))**(-1))* &
                                ((dpcell(y,x)%plist(pp)%density*dpcell(y,x)%pplist(pp)%porosity**(-1))+ &
                                (dpcell(i,j)%plist(k)%density*dpcell(i,j)%pplist(k)%porosity**(-1)))**(-1)

                            if (dpcell(i,j)%plist(k)%tid==3) then

                                t1=2*pvol*(&
                                (dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)* &
                                Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))) &
                                *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)*(p_dist)

                                t2=2*pvol*(&
                                (dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)* &
                                Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))) &
                                *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)*(p_dist)                              

                                if (.not.(dpcell(y,x)%plist(pp)%buffer)) then
                                fmatrix(pos)%val(m)=(-(t1+t2))*lamp
                                fmatrix(pos)%col(m)=dpcell(y,x)%plist(pp)%matid

                                else 
                                    fvec(pos)=fvec(pos)+(t1+t2)*lamp*dpcell(y,x)%plist(pp)%pressure

                                end if

                                fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(t1+t2)


                                if (dpcell(y,x)%plist(pp)%tid==3) then


                                t1=(pvol*(dpcell(y,x)%plist(pp)%vxs*rho_j/dpcell(i,j)%pplist(k)%porosity- &
                                dpcell(i,j)%plist(k)%vxs*rho_i/dpcell(i,j)%pplist(k)%porosity)* &
                                (dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)* &
                                Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))


                                t2=(pvol*(dpcell(y,x)%plist(pp)%vys*rho_j/dpcell(i,j)%pplist(k)%porosity- &
                                dpcell(i,j)%plist(k)%vys*rho_i/dpcell(i,j)%pplist(k)%porosity)* &
                                (dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)* &
                                Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))

                                else 

                                    t1=(pvol*(dpcell(y,x)%plist(pp)%vxs*rho_j/dpcell(i,j)%pplist(k)%porosity- &
                                    dpcell(i,j)%plist(k)%vxs*rho_i/dpcell(i,j)%pplist(k)%porosity)* &
                                    (dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                    dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)* &
                                    Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))


                                    t2=(pvol*(dpcell(y,x)%plist(pp)%vys*rho_j/dpcell(i,j)%pplist(k)%porosity- &
                                    dpcell(i,j)%plist(k)%vys*rho_i/dpcell(i,j)%pplist(k)%porosity)* &
                                    (dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                    dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)* &
                                    Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))

                                end if

                                fvec(pos)=fvec(pos)+(t1+t2)*lamp/real(dt,dp) !lamp



                            elseif (dpcell(i,j)%plist(k)%tid==2) then


                                t1=2*pvol*(&
                                Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)*(p_dist)

                                t2=2*pvol*(&
                                Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)*(p_dist)
                                
                                if (dpcell(y,x)%plist(pp)%tid<=2) then

                                    fmatrix(pos)%val(m)=-(t1+t2)!*dpcell(y,x)%pplist(pp)%lamp
                                        fmatrix(pos)%col(m)=dpcell(y,x)%plist(pp)%matid
                                        fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(t1+t2)!*dpcell(y,x)%pplist(pp)%lamp
    
                                end if


                                ! fmatrix(pos)%val(1)=-1.0_dp!(t1+t2)
                                ! fmatrix(pos)%col(1)=dpcell(i,j)%plist(k)%wall(1)!dpcell(y,x)%plist(pp)%matid
                                ! fmatrix(pos)%val(num+1)=1.0_dp!fmatrix(pos)%val(num+1)+(t1+t2)

                                ! end if

                            elseif (dpcell(i,j)%plist(k)%tid==1) then

                                if (dpcell(y,x)%plist(pp)%tid==3) then

                                    if (dpcell(y,x)%plist(pp)%free) then

                                        lamp2=(1.0_dp+dpcell(i,j)%plist(k)%density/(rho*1.0))**(-1)
        
                                    ! elseif((pll<dpcell(y,x)%pplist(pp)%gradvx).and.(dpcell(y,x)%pplist(pp)%gradvx<pul)) then
        
                                    ! lamp2=((1.0_dp+dpcell(i,j)%plist(k)%density/(rho*1.0))**(-1)) &
                                    ! *(1.0_dp-cos(22.0_dp*(dpcell(y,x)%pplist(pp)%gradvx-pll)/ &
                                    !                 ((pul-pll)*7.0_dp)))
                                    else 
        
                                        lamp2 =1.0_dp
        
                                    end if

                                    if (dpcell(y,x)%plist(pp)%buffer) then

                                        lamp2=1.0_dp

                                    end if

                                t1=2*dpcell(i,j)%plist(k)%mass*(&
                                (dpcell(y,x)%pplist(pp)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                dpcell(i,j)%list(k)%dist(m),h1)+dpcell(y,x)%pplist(pp)%coff(2)* &
                                Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))) &
                                *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)*(p_dist)/ &
                                ((dpcell(i,j)%plist(k)%density))

                                t2=2*dpcell(i,j)%plist(k)%mass*(&
                                (dpcell(y,x)%pplist(pp)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                dpcell(i,j)%list(k)%dist(m),h1)+dpcell(y,x)%pplist(pp)%coff(4)* &
                                Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))) &
                                *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)*(p_dist)/ &
                                ((dpcell(i,j)%plist(k)%density))

                                ! t1=2*pvol*(&
                                ! Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                ! *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)*(p_dist)

                                ! t2=2*pvol*(&
                                ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                ! *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)*(p_dist)

                                if (.not.(dpcell(y,x)%plist(pp)%buffer)) then
                                fmatrix(pos)%val(m)=-(t1+t2)*lamp2!*dpcell(y,x)%pplist(pp)%lamp
                                fmatrix(pos)%col(m)=dpcell(y,x)%plist(pp)%matid
                                else 
                                    fvec(pos)=fvec(pos)+(t1+t2)*lamp2*dpcell(y,x)%plist(pp)%pressure
                                end if

                                fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(t1+t2)!*dpcell(y,x)%pplist(pp)%lamp


                                    t1=(pvol*(dpcell(y,x)%plist(pp)%vxs*rho_j/dpcell(i,j)%pplist(k)%porosity- &
                                    dpcell(i,j)%plist(k)%vxs*rho_i/dpcell(i,j)%pplist(k)%porosity)* &
                                    (dpcell(y,x)%pplist(pp)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                    dpcell(i,j)%list(k)%dist(m),h1)+dpcell(y,x)%pplist(pp)%coff(2)* &
                                    Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))


                                    t2=(pvol*(dpcell(y,x)%plist(pp)%vys*rho_j/dpcell(i,j)%pplist(k)%porosity- &
                                    dpcell(i,j)%plist(k)%vys*rho_i/dpcell(i,j)%pplist(k)%porosity)* &
                                    (dpcell(y,x)%pplist(pp)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                    dpcell(i,j)%list(k)%dist(m),h1)+dpcell(y,x)%pplist(pp)%coff(4)* &
                                    Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))


                                    fvec(pos)=fvec(pos)+(t1+t2)*lamp2/real(dt,dp)  !lamp

                                
                                elseif (dpcell(y,x)%plist(pp)%tid<=2) then

                                    ! t1=2*dpcell(y,x)%plist(pp)%mass*(&
                                    ! Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                    ! *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)/ &
                                    ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                    ! (dpcell(y,x)%plist(pp)%density*dpcell(i,j)%plist(k)%density))

                                    ! t2=2*dpcell(y,x)%plist(pp)%mass*(&
                                    ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                    ! *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)/ &
                                    ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                    ! (dpcell(y,x)%plist(pp)%density*dpcell(i,j)%plist(k)%density))

                                    t1=2*pvol*(&
                                    Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                    *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)*(p_dist)

                                    t2=2*pvol*(&
                                    Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                    *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)*(p_dist)

                                    fmatrix(pos)%val(m)=-(t1+t2)!*dpcell(y,x)%pplist(pp)%lamp
                                    fmatrix(pos)%col(m)=dpcell(y,x)%plist(pp)%matid
                                    fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(t1+t2)!*dpcell(y,x)%pplist(pp)%lamp
                                    

                                end if
                            
                            end if

                        ! end if


                        end associate
                    end do
                    end if

                    ! if (.not.(dpcell(i,j)%plist(k)%free)) then

                    ! fvec(pos)=fvec(pos)+alpha*(1.0_dp-dpcell(i,j)%plist(k)%tden/dpcell(i,j)%plist(k)%density)/(dt**2)
                    ! end if

                    fmatrix(pos)%col(num+1)=pos
                    if ((num==0)) then  !.or.( dpcell(i,j)%plist(k)%free==1)
                    fmatrix(pos)%val(num+1)=1.0_dp
                    end if

                    t1=fmatrix(pos)%val(num+1)

                    ! if (dpcell(i,j)%plist(k)%tid/=2) then

                    fvec(pos)=fvec(pos)/real(t1,dp)   ! Jacobi Preconditioning
                    
                    do m=1,ceiling(fac2*fplistmax)
                    fmatrix(pos)%val(m)=fmatrix(pos)%val(m)/real(t1,dp) !Jacobi Preconditioning
                    end do 

                    ! end if


                end associate

                ! end if

                end do

            ! end if            
            end do
        end do
        !$omp end parallel do

        ! call format(fmatrix,fval,frow,fcol,finmax)
        ! call bicgstab(tl,fguess,finmax,fval,frow,fcol,fvec,fsol)
        call fgmres
        
        ! Assigning pressures
        !$omp parallel do schedule(runtime) default(shared) private(i,k,j) collapse(2)       
            do j=sx,ex
                do i=sy,ey
                    if (dpcell(i,j)%ptot/=0) then

                        do k=1,dpcell(i,j)%ptot

                        if (.not.(dpcell(i,j)%plist(k)%buffer)) then
                        
                        dpcell(i,j)%plist(k)%pressure=fsol(dpcell(i,j)%plist(k)%matid)

                        ! if (dpcell(i,j)%plist(k)%tid<=2) then

                        !     dpcell(i,j)%plist(k)%pressure=dpcell(i,j)%plist(k)%pressure &
                        !     +2*rho*abs(g)*min(prrealx,prrealy) 

                        end if

                        ! end if

                        end do

                    end if
                end do
            end do
        !$omp end parallel do 

    end subroutine ppesolve

end module isph



