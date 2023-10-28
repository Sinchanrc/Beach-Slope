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

        real(dp) :: heff


        ! Free surface(dpcell(y,x)%plist(pp)%tid==3) then
        !$omp do private(m,term1,term2,i,k,j,heff) schedule (runtime) collapse(2)    
        do j=sx,ex
            do i=sy,ey
            if(dpcell(i,j)%ptot/=0) then
                do k=1,dpcell(i,j)%ptot
                ! Free Surface Identification
                if (dpcell(i,j)%plist(k)%tid==3) then
                dpcell(i,j)%plist(k)%free=.false.
                dpcell(i,j)%pplist(k)%gradvx=0.0_dp
                term1=0.0_dp
                term2=0.0_dp

                ! if (dpcell(i,j)%pplist(k)%inpore) then

                heff=h1!(max(prrealx,prrealy)*hfac)/(2.0_dp*sqrt(dpcell(i,j)%pplist(k)%porosity))  !*sqrt(por)

                ! else 

                ! heff=(max(prrealx,prrealy)*hfac)/(2.0_dp)

                ! end if

                do m=1,dpcell(i,j)%list(k)%count
                    associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                    y=>dpcell(i,j)%list(k)%interlist(2,m), &
                    pp=>dpcell(i,j)%list(k)%interlist(3,m))

                    if ((dpcell(y,x)%plist(pp)%tid/=4).and. &
                    (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

                    term1=term1+((dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%plist(pp)%x- &
                    dpcell(i,j)%plist(k)%x)*Wabx(dpcell(y,x)%plist(pp),&
                    dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))/&
                    dpcell(y,x)%plist(pp)%density)

                    term2=term2+((dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%plist(pp)%y- &
                    dpcell(i,j)%plist(k)%y)*Waby(dpcell(y,x)%plist(pp),&
                    dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))/&
                    dpcell(y,x)%plist(pp)%density)

                    end if

                    end associate
                end do

                dpcell(i,j)%pplist(k)%gradvx=term1+term2
                if(dpcell(i,j)%pplist(k)%gradvx<=(lamfs*maxdivr)) then
                    dpcell(i,j)%plist(k)%free=.true.
                
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

    subroutine ppesolve

        implicit none

        ! real(dp) :: lamp
        
        ! Preparing the coff matrix for fluid part in CSR
        !$omp parallel do schedule(runtime) default(shared) &
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

                                lamp=0.50_dp!(1.0_dp+dpcell(y,x)%plist(pp)%density/(dpcell(y,x)%pplist(k)%porosity*rhomax))

                            elseif((pll<dpcell(i,j)%pplist(k)%gradvx).and.(dpcell(i,j)%pplist(k)%gradvx<pul)) then

                            lamp=0.50_dp*(1.0_dp-cos(22.0_dp*(dpcell(i,j)%pplist(k)%gradvx-pll)/ &
                                            ((pul-pll)*7.0_dp)))
                            else 

                                lamp =1.0_dp

                            end if

                            if (dpcell(i,j)%plist(k)%tid>2) then

                                ! term1=4*dpcell(y,x)%plist(pp)%mass*dpcell(i,j)%pplist(k)%porosity*(&
                                ! Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                ! *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)/ &
                                ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                ! (dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(k)%porosity &
                                ! +dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)*dpcell(y,x)%plist(pp)%density)

                                ! term2=4*dpcell(y,x)%plist(pp)%mass*dpcell(i,j)%pplist(k)%porosity*(&
                                ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                ! *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)/ &
                                ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                ! (dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(k)%porosity &
                                ! +dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)*dpcell(y,x)%plist(pp)%density) 

                                term1=2*dpcell(y,x)%plist(pp)%mass*(dpcell(i,j)%pplist(k)%porosity**2)*(&
                                Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)/ &
                                (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                (dpcell(y,x)%plist(pp)%density))

                                term2=2*dpcell(y,x)%plist(pp)%mass*(dpcell(i,j)%pplist(k)%porosity**2)*(&
                                Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)/ &
                                (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                (dpcell(y,x)%plist(pp)%density))                               

                                fmatrix(pos)%val(m)=(-(term1+term2))*lamp
                                fmatrix(pos)%col(m)=dpcell(y,x)%plist(pp)%matid

                                fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(term1+term2)


                                term1=(dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%plist(pp)%vxs*dpcell(y,x)%plist(pp)%density- &
                                dpcell(i,j)%plist(k)%vxs*dpcell(i,j)%plist(k)%density)* &
                                (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                dpcell(i,j)%list(k)%dist(m),h1)))/ &
                                (dpcell(y,x)%plist(pp)%density)


                                term2=(dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%plist(pp)%vys*dpcell(y,x)%plist(pp)%density- &
                                dpcell(i,j)%plist(k)%vys*dpcell(i,j)%plist(k)%density)* &
                                (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                dpcell(i,j)%list(k)%dist(m),h1)))/ &
                                (dpcell(y,x)%plist(pp)%density)

                                fvec(pos)=fvec(pos)+(term1+term2)*lamp/real(dt,dp) !lamp



                            elseif (dpcell(i,j)%plist(k)%tid==2) then

                                term1=4*dpcell(y,x)%plist(pp)%mass*dpcell(i,j)%pplist(k)%porosity*(&
                                Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)/ &
                                (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                (dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(k)%porosity &
                                +dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)*dpcell(y,x)%plist(pp)%density)

                                term2=4*dpcell(y,x)%plist(pp)%mass*dpcell(i,j)%pplist(k)%porosity*(&
                                Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)/ &
                                (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                (dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(k)%porosity &
                                +dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)*dpcell(y,x)%plist(pp)%density) 

                            ! term1=2*dpcell(y,x)%plist(pp)%mass*(dpcell(i,j)%pplist(k)%porosity**2)*(&
                            ! Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                            ! *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)/ &
                            ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                            ! (dpcell(y,x)%plist(pp)%density*dpcell(i,j)%plist(k)%density))

                            ! term2=2*dpcell(y,x)%plist(pp)%mass*(dpcell(i,j)%pplist(k)%porosity**2)*(&
                            ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                            ! *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)/ &
                            ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                            ! (dpcell(y,x)%plist(pp)%density*dpcell(i,j)%plist(k)%density)) 
                                
                                if (dpcell(y,x)%plist(pp)%tid<=2) then

                                    fmatrix(pos)%val(m)=-(term1+term2)!*dpcell(y,x)%pplist(pp)%lamp
                                        fmatrix(pos)%col(m)=dpcell(y,x)%plist(pp)%matid
                                        fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(term1+term2)!*dpcell(y,x)%pplist(pp)%lamp
    
                                end if


                                ! fmatrix(pos)%val(1)=-1.0_dp!(term1+term2)
                                ! fmatrix(pos)%col(1)=dpcell(i,j)%plist(k)%wall(1)!dpcell(y,x)%plist(pp)%matid
                                ! fmatrix(pos)%val(num+1)=1.0_dp!fmatrix(pos)%val(num+1)+(term1+term2)

                                ! end if

                            elseif (dpcell(i,j)%plist(k)%tid==1) then

                                if (dpcell(y,x)%plist(pp)%tid==3) then

                                !     term1=4*dpcell(y,x)%plist(pp)%mass*dpcell(i,j)%pplist(k)%porosity*(&
                                ! Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                ! *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)/ &
                                ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                ! (dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(k)%porosity &
                                ! +dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)*dpcell(y,x)%plist(pp)%density)

                                ! term2=4*dpcell(y,x)%plist(pp)%mass*dpcell(i,j)%pplist(k)%porosity*(&
                                ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                ! *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)/ &
                                ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                ! (dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(k)%porosity &
                                ! +dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)*dpcell(y,x)%plist(pp)%density) 

                                !     term1=4*dpcell(i,j)%plist(k)%mass*dpcell(y,x)%pplist(k)%porosity*(&
                                ! Wabx(dpcell(i,j)%plist(k),dpcell(y,x)%plist(pp),dpcell(i,j)%list(k)%dist(m),h1)) &
                                ! *(-dpcell(i,j)%plist(k)%x+dpcell(y,x)%plist(pp)%x)/ &
                                ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                ! (dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(k)%porosity &
                                ! +dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)*dpcell(i,j)%plist(k)%density)

                                ! term2=4*dpcell(i,j)%plist(k)%mass*dpcell(y,x)%pplist(k)%porosity*(&
                                ! Waby(dpcell(i,j)%plist(k),dpcell(y,x)%plist(pp),dpcell(i,j)%list(k)%dist(m),h1)) &
                                ! *(-dpcell(i,j)%plist(k)%y+dpcell(y,x)%plist(pp)%y)/ &
                                ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                ! (dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(k)%porosity &
                                ! +dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)*dpcell(i,j)%plist(k)%density)

                                term1=2*dpcell(i,j)%plist(k)%mass*(dpcell(y,x)%pplist(pp)%porosity**2)*(&
                                Wabx(dpcell(i,j)%plist(k),dpcell(y,x)%plist(pp),dpcell(i,j)%list(k)%dist(m),h1)) &
                                *(-dpcell(i,j)%plist(k)%x+dpcell(y,x)%plist(pp)%x)/ &
                                (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                (dpcell(i,j)%plist(k)%density))

                                term2=2*dpcell(i,j)%plist(k)%mass*(dpcell(y,x)%pplist(pp)%porosity**2)*(&
                                Waby(dpcell(i,j)%plist(k),dpcell(y,x)%plist(pp),dpcell(i,j)%list(k)%dist(m),h1)) &
                                *(-dpcell(i,j)%plist(k)%y+dpcell(y,x)%plist(pp)%y)/ &
                                (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                (dpcell(i,j)%plist(k)%density))

                                ! term1=2*dpcell(y,x)%plist(pp)%mass*(&
                                ! Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                ! *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)/ &
                                ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                ! (dpcell(y,x)%plist(pp)%density))

                                ! term2=2*dpcell(y,x)%plist(pp)%mass*(&
                                ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                ! *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)/ &
                                ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                ! (dpcell(y,x)%plist(pp)%density))

                                fmatrix(pos)%val(m)=-(term1+term2)!*dpcell(y,x)%pplist(pp)%lamp
                                fmatrix(pos)%col(m)=dpcell(y,x)%plist(pp)%matid
                                fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(term1+term2)!*dpcell(y,x)%pplist(pp)%lamp


                                    term1=(dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%plist(pp)%vxs*dpcell(y,x)%plist(pp)%density- &
                                    dpcell(i,j)%plist(k)%vxs*dpcell(i,j)%plist(k)%density)* &
                                    (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                    dpcell(i,j)%list(k)%dist(m),h1)))/ &
                                    (dpcell(y,x)%plist(pp)%density)


                                    term2=(dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%plist(pp)%vys*dpcell(y,x)%plist(pp)%density- &
                                    dpcell(i,j)%plist(k)%vys*dpcell(i,j)%plist(k)%density)* &
                                    (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                    dpcell(i,j)%list(k)%dist(m),h1)))/ &
                                    (dpcell(y,x)%plist(pp)%density)

                                    fvec(pos)=fvec(pos)+(term1+term2)/real(dt,dp)  !lamp

                                
                                elseif (dpcell(y,x)%plist(pp)%tid<=2) then

                                    term1=4*dpcell(y,x)%plist(pp)%mass*dpcell(i,j)%pplist(k)%porosity*(&
                                    Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                    *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)/ &
                                    (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                    (dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(k)%porosity &
                                    +dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)*dpcell(y,x)%plist(pp)%density)

                                    term2=4*dpcell(y,x)%plist(pp)%mass*dpcell(i,j)%pplist(k)%porosity*(&
                                    Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                    *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)/ &
                                    (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                    (dpcell(y,x)%plist(pp)%density/dpcell(y,x)%pplist(k)%porosity &
                                    +dpcell(i,j)%plist(k)%density/dpcell(i,j)%pplist(k)%porosity)*dpcell(y,x)%plist(pp)%density)

                                    ! term1=2*dpcell(y,x)%plist(pp)%mass*(&
                                    ! Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                    ! *(dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)/ &
                                    ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                    ! (dpcell(y,x)%plist(pp)%density*dpcell(i,j)%plist(k)%density))

                                    ! term2=2*dpcell(y,x)%plist(pp)%mass*(&
                                    ! Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
                                    ! *(dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)/ &
                                    ! (((dist(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k))**2)+lam)* &
                                    ! (dpcell(y,x)%plist(pp)%density*dpcell(i,j)%plist(k)%density))

                                    fmatrix(pos)%val(m)=-(term1+term2)!*dpcell(y,x)%pplist(pp)%lamp
                                    fmatrix(pos)%col(m)=dpcell(y,x)%plist(pp)%matid
                                    fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(term1+term2)!*dpcell(y,x)%pplist(pp)%lamp

                                    ! term1=(dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%plist(pp)%vxs- &
                                    ! dpcell(i,j)%plist(k)%vxs)* &
                                    ! (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                    ! dpcell(i,j)%list(k)%dist(m),h1)))/ &
                                    ! (dpcell(y,x)%plist(pp)%density)

                                    ! term2=(dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%plist(pp)%vys- &
                                    ! dpcell(i,j)%plist(k)%vys)* &
                                    ! (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                                    ! dpcell(i,j)%list(k)%dist(m),h1)))/ &
                                    ! (dpcell(y,x)%plist(pp)%density)

                                    ! fvec(pos)=fvec(pos)+(term1+term2)/real(dt,dp)  !lamp
                                    

                                end if
                            
                            end if

                        end if


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
        !$omp end parallel do

        call pformat(fmatrix,fval,frow,fcol,finmax)
        call bicgstab(tl,fguess,finmax,fval,frow,fcol,fvec,fsol)
        
        ! Assigning pressures
        !$omp parallel do schedule(runtime) default(shared) private(i,k,j) collapse(2)       
            do j=sx,ex
                do i=sy,ey
                    if (dpcell(i,j)%ptot/=0) then

                        do k=1,dpcell(i,j)%ptot
                        
                        dpcell(i,j)%plist(k)%pressure=fsol(dpcell(i,j)%plist(k)%matid)

                        ! if (dpcell(i,j)%plist(k)%tid<=2) then

                        !     dpcell(i,j)%plist(k)%pressure=dpcell(i,j)%plist(k)%pressure &
                        !     +2*rho*abs(g)*min(prrealx,prrealy) 

                        ! end if

                        end do

                    end if
                end do
            end do
        !$omp end parallel do 

    end subroutine ppesolve

end module isph



