module turbulence

    use particle
    use initialize
    use domain
    use kernel

    implicit none

    contains

    subroutine eddyvis

        implicit none
        integer :: i,j,k,m
        
        !Mean strain and SGS tensor calculation for fluid particles
        !$omp parallel do schedule (runtime) private(m,i,k,j) collapse(2) default(shared)
        do j=sx,ex
            do i=sy,ey
                if(dpcell(i,j)%ptot/=0) then
                
                do k=1,dpcell(i,j)%ptot
                
                if ((dpcell(i,j)%plist(k)%tid==3).and.(dpcell(i,j)%pplist(k)%gradvx>=0.70_dp)) then
                    dpcell(i,j)%pplist(k)%meanstrn=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count ! Calculating mean strain tensor
                        associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                        y=>dpcell(i,j)%list(k)%interlist(2,m), &
                        pp=>dpcell(i,j)%list(k)%interlist(3,m))

                        if (dpcell(y,x)%plist(pp)%tid/=4) then

                        dpcell(i,j)%pplist(k)%meanstrn(1)=dpcell(i,j)%pplist(k)%meanstrn(1)+ &
                            (dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%plist(pp)%vx/dpcell(y,x)%pplist(pp)%porosity- &
                            dpcell(i,j)%plist(k)%vx/dpcell(i,j)%pplist(k)%porosity)* &
                            (dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                            dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)* &
                            Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))) &
                            /dpcell(y,x)%plist(pp)%density

                        dpcell(i,j)%pplist(k)%meanstrn(4)=dpcell(i,j)%pplist(k)%meanstrn(4)+ &
                            (dpcell(y,x)%plist(pp)%mass*(dpcell(y,x)%plist(pp)%vy/dpcell(y,x)%pplist(pp)%porosity- &
                            dpcell(i,j)%plist(k)%vy/dpcell(i,j)%pplist(k)%porosity)* &
                            (dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                            dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)* &
                            Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))) &
                            /dpcell(y,x)%plist(pp)%density

                        dpcell(i,j)%pplist(k)%meanstrn(2)=dpcell(i,j)%pplist(k)%meanstrn(2)+ &
                        (0.50_dp*dpcell(y,x)%plist(pp)%mass*((dpcell(y,x)%plist(pp)%vx/dpcell(y,x)%pplist(pp)%porosity- &
                        dpcell(i,j)%plist(k)%vx/dpcell(i,j)%pplist(k)%porosity)* &
                        (dpcell(i,j)%pplist(k)%coff(3)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))+ &
                        (dpcell(y,x)%plist(pp)%vy/dpcell(y,x)%pplist(pp)%porosity- &
                        dpcell(i,j)%plist(k)%vy/dpcell(i,j)%pplist(k)%porosity)* &
                        (dpcell(i,j)%pplist(k)%coff(1)*Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                        dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)* &
                        Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))) &
                        /dpcell(y,x)%plist(pp)%density

                        end if

                        end associate
                    end do
                    
                    dpcell(i,j)%pplist(k)%meanstrn(3)=dpcell(i,j)%pplist(k)%meanstrn(2)


                    ! Calculating local strain

                    ! dpcell(i,j)%pplist(k)%S=sqrt(2*(dpcell(i,j)%pplist(k)%meanstrn(1)**2+ &
                    ! dpcell(i,j)%pplist(k)%meanstrn(4)**2))

                    dpcell(i,j)%pplist(k)%S=sqrt(2*(dpcell(i,j)%pplist(k)%meanstrn(1)**2+ &
                    dpcell(i,j)%pplist(k)%meanstrn(4)**2)+(dpcell(i,j)%pplist(k)%meanstrn(2) &
                    +dpcell(i,j)%pplist(k)%meanstrn(3))**2)

                    ! Calculating Eddy viscosity
                    dpcell(i,j)%pplist(k)%nut=((cs*dl)**2)*dpcell(i,j)%pplist(k)%S

                    ! Calculating tke
                    dpcell(i,j)%pplist(k)%tke=(cv/ce)*(dl*dpcell(i,j)%pplist(k)%S)**2

                    ! Calculating SGS stress tensor
                    ! dpcell(i,j)%pplist(k)%tau(1)=2*dpcell(i,j)%pplist(k)%meanstrn(1)*dpcell(i,j)%pplist(k)%nut &
                    ! -2*dpcell(i,j)%pplist(k)%tke/3.0_dp
                    
                    ! dpcell(i,j)%pplist(k)%tau(2)=2*dpcell(i,j)%pplist(k)%meanstrn(2)*dpcell(i,j)%pplist(k)%nut

                    ! dpcell(i,j)%pplist(k)%tau(3)=2*dpcell(i,j)%pplist(k)%meanstrn(3)*dpcell(i,j)%pplist(k)%nut

                    ! dpcell(i,j)%pplist(k)%tau(4)=2*dpcell(i,j)%pplist(k)%meanstrn(4)*dpcell(i,j)%pplist(k)%nut &
                    ! -2*dpcell(i,j)%pplist(k)%tke/3.0_dp
                else
                    dpcell(i,j)%pplist(k)%tau(:)=0.0_dp
                    dpcell(i,j)%pplist(k)%tke=0.0_dp
                    dpcell(i,j)%pplist(k)%nut=0.0_dp
                end if
                end do
                end if         
            end do
        end do
        !$omp end parallel do
        
        
    end subroutine eddyvis
    

end module turbulence
