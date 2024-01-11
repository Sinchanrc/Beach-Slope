module part_shift

    use particle
    use domain
    use kernel
    use initialize

    implicit none

    contains

    ! subroutine exp_shift
    !     implicit none
        
    !     real(dp) :: t1,t2
    !     integer :: i,j,k,m
    !     ! Explicit particle shifting
    !     !$omp parallel do private(i,j,k,m,t1,t2) schedule (runtime) collapse(2) default(shared)
    !         do j=sx,ex
    !             do i=sy,ey
    !                 if(dpcell(i,j)%ptot/=0) then
                    
    !                 do k=1,dpcell(i,j)%ptot
    !                     dpcell(i,j)%plist(k)%ys=0.0_dp
    !                     dpcell(i,j)%plist(k)%xs=0.0_dp
    !                     dpcell(i,j)%plist(k)%vys=0.0_dp
    !                     dpcell(i,j)%plist(k)%vxs=0.0_dp
    !                     if ((dpcell(i,j)%list(k)%count/=0) &
    !                     ! .and.(dpcell(i,j)%plist(k)%mirror.eqv..true.) &
    !                     .and.(dpcell(i,j)%pplist(k)%gradvx>=0.70_dp) &
    !                     .and.(dpcell(i,j)%plist(k)%tid>2)) then
    !                     dpcell(i,j)%plist(k)%ys=0.0_dp
    !                     dpcell(i,j)%plist(k)%xs=0.0_dp
    !                     t1=0.0_dp
    !                     t2=0.0_dp

    !                     do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                         y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                         pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                         if (dpcell(i,j)%plist(k)%tid/=4) then

    !                         t1=t1+dpcell(i,j)%list(k)%dist(m)/ &
    !                         (dpcell(i,j)%list(k)%count)
                        
    !                         t2=t2+(dpcell(i,j)%plist(k)%x- &
    !                         x%x)/ &
    !                         dpcell(i,j)%list(k)%dist(m)**3

    !                         end if

    !                         end associate
    !                     end do

    !                     dpcell(i,j)%plist(k)%xs=beta*umax*dt*(t1**2)*t2

    !                     t1=0.0_dp
    !                     t2=0.0_dp
    !                     if (dpcell(i,j)%list(k)%count/=0) then
    !                     do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                         y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                         pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                         if (dpcell(i,j)%plist(k)%tid/=4) then

    !                         t1=t1+dpcell(i,j)%list(k)%dist(m)/ &
    !                         (dpcell(i,j)%list(k)%count)
                        
    !                         t2=t2+(dpcell(i,j)%plist(k)%y- &
    !                         x%y)/ &
    !                         dpcell(i,j)%list(k)%dist(m)**3

    !                         end if

    !                         end associate
    !                     end do
    !                     end if

    !                     dpcell(i,j)%plist(k)%ys=beta*umax*dt*(t1**2)*t2

    !                     t1=0.0_dp
    !                     t2=0.0_dp  

    !                     do m=1,dpcell(i,j)%list(k)%count
    !                         associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                         y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                         pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                         if (dpcell(i,j)%plist(k)%tid/=4) then

    !                         ! t1=t1+(Wabx(x,dpcell(i,j)%plist(k),&
    !                         ! dpcell(i,j)%list(k)%dist(m),h1)*(x%vx- &
    !                         ! dpcell(i,j)%plist(k)%vx))&
    !                         ! *(x%mass/x%density)

    !                         t1=t1+((dpcell(i,j)%pplist(k)%coff(1)*Wabx(x,dpcell(i,j)%plist(k),&
    !                         dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                         dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                         (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                         *(x%mass/x%density)

    !                         ! t2=t2+(Waby(x,dpcell(i,j)%plist(k),&
    !                         ! dpcell(i,j)%list(k)%dist(m),h1)*(x%vx- &
    !                         ! dpcell(i,j)%plist(k)%vx))&
    !                         ! *(x%mass/x%density)

    !                         t2=t2+((dpcell(i,j)%pplist(k)%coff(3)*Wabx(x,dpcell(i,j)%plist(k),&
    !                         dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                         dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                         (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                         *(x%mass/x%density)

    !                         end if

    !                         end associate
    !                     end do

    !                     dpcell(i,j)%plist(k)%vxs=t1*dpcell(i,j)%plist(k)%xs &
    !                                                     +t2*dpcell(i,j)%plist(k)%ys


    !                     t1=0.0_dp
    !                     t2=0.0_dp

    !                     do m=1,dpcell(i,j)%list(k)%count
    !                         associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                         y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                         pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                         if (dpcell(i,j)%plist(k)%tid/=4) then

    !                         ! t1=t1+(Wabx(x,dpcell(i,j)%plist(k),&
    !                         ! dpcell(i,j)%list(k)%dist(m),h1)*(x%vy- &
    !                         ! dpcell(i,j)%plist(k)%vy))&
    !                         ! *(x%mass/x%density)

    !                         t1=t1+((dpcell(i,j)%pplist(k)%coff(1)*Wabx(x,dpcell(i,j)%plist(k),&
    !                         dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                         dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                         (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                         *(x%mass/x%density)

    !                         ! t2=t2+(Waby(x,dpcell(i,j)%plist(k),&
    !                         ! dpcell(i,j)%list(k)%dist(m),h1)*(x%vy- &
    !                         ! dpcell(i,j)%plist(k)%vy))&
    !                         ! *(x%mass/x%density)

    !                         t2=t2+((dpcell(i,j)%pplist(k)%coff(3)*Wabx(x,dpcell(i,j)%plist(k),&
    !                         dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                         dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                         (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                         *(x%mass/x%density)

    !                         end if

    !                         end associate
    !                     end do

    !                     dpcell(i,j)%plist(k)%vys=t1*dpcell(i,j)%plist(k)%xs &
    !                                                     +t2*dpcell(i,j)%plist(k)%ys


    !                     ! dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+rx
    !                     ! dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+ry

    !                     end if
    !                 end do  
                                    
    !                 end if
    !             end do
    !         end do
    !     !$omp end parallel do

    !     ! Shifting hydrodynamic values
    !     !$omp parallel do schedule(dynamic) private(i,j,k) collapse(2) default(shared)
    !         do j=sx,ex
    !         do i=sy,ey
    !             if (dpcell(i,j)%ptot/=0) then
    !             do k=1,dpcell(i,j)%ptot
    !                 if (dpcell(i,j)%plist(k)%tid>2) then
    !                 dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+dpcell(i,j)%plist(k)%xs
    !                 dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+dpcell(i,j)%plist(k)%ys
    !                 dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx+dpcell(i,j)%plist(k)%vxs
    !                 dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy+dpcell(i,j)%plist(k)%vys
    !                 ! dpcell(i,j)%plist(k)%pressure=dpcell(i,j)%plist(k)%gradvx
    !                 end if
    !             end do
    !             end if
    !         end do
    !         end do
    !     !$omp end parallel do   
        
    ! end subroutine exp_shift

    ! subroutine opt_shift

    !     implicit none

    !     integer :: pint=0
    !     integer :: i,j,k,m
    !     real(dp) :: t1,t2
        
    !     ! Optimized Particle Shifting
    !     !$omp do private(i,k,m,t1,t2,normx,normy) schedule (runtime) collapse(2)
    !         do j=sx,ex
    !         do i=sy,ey
    !         if(dpcell(i,j)%ptot/=0) then

    !         do k=1,dpcell(i,j)%ptot

    !             if (dpcell(i,j)%plist(k)%tid==3) then

    !             dpcell(i,j)%plist(k)%xs=0.0_dp
    !             dpcell(i,j)%plist(k)%ys=0.0_dp
    !             dpcell(i,j)%plist(k)%vxs=0.0_dp
    !             dpcell(i,j)%plist(k)%vys=0.0_dp

    !             if ((dpcell(i,j)%list(k)%count/=0) &
    !             .and.(dpcell(i,j)%pplist(k)%gradvx>=0.69_dp)) then
    !                 if ((dpcell(i,j)%plist(k)%free).or.(dpcell(i,j)%plist(k)%vicinity)) then
    !                 t1=0.0_dp
    !                 t2=0.0_dp
    !                 normx=0.0_dp
    !                 normy=0.0_dp
    !                 dpcell(i,j)%plist(k)%xs=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4) then

    !                     normx=normx-(dpcell(i,j)%pplist(k)%coff(1)*&
    !                     Wabx(x,dpcell(i,j)%plist(k), &
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)* &
    !                     Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))&
    !                     *(x%mass/x%density)

    !                     normy=normy-(dpcell(i,j)%pplist(k)%coff(3)*&
    !                     Wabx(x,dpcell(i,j)%plist(k), &
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)* &
    !                     Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))&
    !                     *(x%mass/x%density)

    !                     end if
            
    !                     end associate
    !                 end do

    !                 t1=sqrt(normx**2+normy**2)
    !                 normx=normx/t1
    !                 normy=normy/t1

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4) then

    !                     dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-csh*(h1**2)* &
    !                     ((1-normx**2)*Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)-normx*normy* &
    !                     Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%xs=min(abs(dpcell(i,j)%plist(k)%xs),maxshift*dl)* &
    !                 dpcell(i,j)%plist(k)%xs/abs(dpcell(i,j)%plist(k)%xs)

    !                 t1=0.0_dp
    !                 t2=0.0_dp
    !                 dpcell(i,j)%plist(k)%ys=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4) then

    !                     dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%ys-csh*(h1**2)*&
    !                     (-normx*normy*Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)+(1-normy**2)* &
    !                     Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%ys=min(abs(dpcell(i,j)%plist(k)%ys),maxshift*dl)* &
    !                 dpcell(i,j)%plist(k)%ys/abs(dpcell(i,j)%plist(k)%ys)

    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4) then

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vxs=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys



    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4) then

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1)&
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vys=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys

    !                 dpcell(i,j)%plist(k)%vicinity=.false.

    !                 else


    !                 dpcell(i,j)%plist(k)%xs=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4) then

    !                     dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-csh*(h1**2)*&
    !                     (Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1))&
    !                     *(x%mass/x%density)* &
    !                     (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),h1)/ &
    !                     Wab(2*dl,h1))**4)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%xs=min(abs(dpcell(i,j)%plist(k)%xs),maxshift*dl)*&
    !                 dpcell(i,j)%plist(k)%xs/abs(dpcell(i,j)%plist(k)%xs)

    !                 dpcell(i,j)%plist(k)%ys=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4) then

    !                     dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%ys-csh*(h1**2)*&
    !                     (Waby(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1))&
    !                     *(x%mass/x%density)* &
    !                     (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),h1)/ &
    !                     Wab(2*dl,h1))**4)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%ys=min(abs(dpcell(i,j)%plist(k)%ys),maxshift*dl)*&
    !                 dpcell(i,j)%plist(k)%ys/abs(dpcell(i,j)%plist(k)%ys)


    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4) then

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vxs=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys


    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4) then

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2) &
    !                     *Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4) &
    !                     *Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vys=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys

    !                 end if
    !             end if

    !             end if

    !         end do
    !         end if
    !         end do
    !         end do
    !     !$omp end do 

    !     ! Shifting hydrodynamic values
    !     !$omp do schedule(dynamic) private(i,j,k) collapse(2)  
    !         do j=sx,ex
    !         do i=sy,ey
    !             if (dpcell(i,j)%ptot/=0) then
    !             do k=1,dpcell(i,j)%ptot
    !                 if (dpcell(i,j)%plist(k)%tid==3) then
    !                 dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+dpcell(i,j)%plist(k)%xs
    !                 dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+dpcell(i,j)%plist(k)%ys
    !                 dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx+dpcell(i,j)%plist(k)%vxs
    !                 dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy+dpcell(i,j)%plist(k)%vys

    !                 end if
    !             end do
    !             end if
    !         end do
    !         end do
    !     !$omp end do   
        
    ! end subroutine opt_shift

    subroutine opt2_shift

        implicit none

        integer :: pint
        integer :: i,j,k,m
        real(dp) :: heff,frac=0.0010_dp,t1,t2,xs,ys,vxs,vys!0.0420_dp
        
        ! Optimized Particle Shifting
        !$omp do private(i,k,m,t1,t2,normx,normy,heff,pint,xs,ys,vxs,vys) &
        !$omp schedule (runtime) collapse(2)
            do j=sx,ex
            do i=sy,ey
            if(dpcell(i,j)%ptot/=0) then

            do k=1,dpcell(i,j)%ptot

                if ((dpcell(i,j)%plist(k)%tid==3).and. &
                (.not.(dpcell(i,j)%plist(k)%buffer))) then

                xs=0.0_dp
                ys=0.0_dp
                vxs=0.0_dp
                vys=0.0_dp

                heff=h1

                if ((dpcell(i,j)%list(k)%count/=0)) then
                    t1=0.0_dp
                    t2=0.0_dp
                    normx=0.0_dp
                    normy=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%nh(m)%part, &
                            y=>dpcell(i,j)%list(k)%pnh(m)%ppart, &
                            z=>dpcell(i,j)%list(k)%klt)

                        normx=normx-(dpcell(i,j)%pplist(k)%coff(1)*&
                        z(1,m)+dpcell(i,j)%pplist(k)%coff(2)* &
                        z(2,m))&
                        *(x%mass/x%density)

                        normy=normy-(dpcell(i,j)%pplist(k)%coff(3)*&
                        z(1,m)+dpcell(i,j)%pplist(k)%coff(4)* &
                        z(2,m))&
                        *(x%mass/x%density)

                        end associate
                    end do

                    t1=sqrt(normx**2+normy**2)
                    normx=normx/t1
                    normy=normy/t1

                    normx=merge(normx,0.0_dp,dpcell(i,j)%plist(k)%free)
                    normy=merge(normy,0.0_dp,dpcell(i,j)%plist(k)%free)
                    pint=merge(1,1,(.not.(dpcell(i,j)%plist(k)%free)))

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%nh(m)%part, &
                            y=>dpcell(i,j)%list(k)%pnh(m)%ppart, &
                            z=>dpcell(i,j)%list(k)%klt)

                        xs=xs-csh*(heff**2)* &
                        ((1-normx**2)*z(1,m)-normx*normy* &
                        z(2,m))&
                        *(x%mass/x%density)* &
                        (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
                        Wab(dl,heff))**4)

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%xs=xs

                    t1=0.0_dp
                    t2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%nh(m)%part, &
                            y=>dpcell(i,j)%list(k)%pnh(m)%ppart, &
                            z=>dpcell(i,j)%list(k)%klt)

                        ys=ys-csh*(heff**2)*&
                        (-normx*normy*z(1,m)+(1-normy**2)* &
                        z(2,m))&
                        *(x%mass/x%density)* &
                        (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
                        Wab(dl,heff))**4)

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%ys=ys

                    t1=0.0_dp
                    t2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%nh(m)%part, &
                            y=>dpcell(i,j)%list(k)%pnh(m)%ppart, &
                            z=>dpcell(i,j)%list(k)%klt)

                        t1=t1+((dpcell(i,j)%pplist(k)%coff(1) &
                        *z(1,m)+dpcell(i,j)%pplist(k)%coff(2)*z(2,m))*&
                        (x%vx-dpcell(i,j)%plist(k)%vx))&
                        *(x%mass/x%density)

                        t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *z(1,m)+dpcell(i,j)%pplist(k)%coff(4)*z(2,m))*&
                        (x%vx-dpcell(i,j)%plist(k)%vx))&
                        *(x%mass/x%density)

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vxs=t1*xs+t2*ys

                    t1=0.0_dp
                    t2=0.0_dp

                    do m=1,dpcell(i,j)%list(k)%count
                        associate(x=>dpcell(i,j)%list(k)%nh(m)%part, &
                            y=>dpcell(i,j)%list(k)%pnh(m)%ppart, &
                            z=>dpcell(i,j)%list(k)%klt)

                        t1=t1+((dpcell(i,j)%pplist(k)%coff(1)&
                        *z(1,m)+dpcell(i,j)%pplist(k)%coff(2)*z(2,m))*&
                        (x%vy-dpcell(i,j)%plist(k)%vy))&
                        *(x%mass/x%density)

                        t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
                        *z(1,m)+dpcell(i,j)%pplist(k)%coff(4)*z(2,m))*&
                        (x%vy-dpcell(i,j)%plist(k)%vy))&
                        *(x%mass/x%density)

                        end associate
                    end do

                    dpcell(i,j)%plist(k)%vys=t1*xs+t2*ys

                end if

                end if

            end do
            end if
            end do
            end do
        !$omp end do 

        ! Shifting hydrodynamic values
        !$omp do schedule(dynamic) private(i,j,k) collapse(2)  
            do j=sx,ex
            do i=sy,ey
                if (dpcell(i,j)%ptot/=0) then
                do k=1,dpcell(i,j)%ptot
                    if ((dpcell(i,j)%plist(k)%tid==3).and. &
                    (.not.(dpcell(i,j)%plist(k)%buffer))) then
                    dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+dpcell(i,j)%plist(k)%xs
                    dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+dpcell(i,j)%plist(k)%ys
                    dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx+dpcell(i,j)%plist(k)%vxs
                    dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy+dpcell(i,j)%plist(k)%vys

                    end if
                end do
                end if
            end do
            end do
        !$omp end do   
        
    end subroutine opt2_shift

    ! subroutine opt3_shift

    !     implicit none

    !     integer :: pint=0
    !     integer :: i,j,k,m
    !     real(dp) :: heff,pscshift,t1,t2

    !     pscshift=0.0_dp
        
    !     ! Optimized Particle Shifting
    !     !$omp do private(i,k,m,t1,t2,normx,normy,heff) schedule (runtime) collapse(2)
    !         do j=sx,ex
    !         do i=sy,ey
    !         if(dpcell(i,j)%ptot/=0) then

    !         do k=1,dpcell(i,j)%ptot

    !             if (dpcell(i,j)%plist(k)%tid==3) then

    !             dpcell(i,j)%plist(k)%xs=0.0_dp
    !             dpcell(i,j)%plist(k)%ys=0.0_dp
    !             dpcell(i,j)%plist(k)%vxs=0.0_dp
    !             dpcell(i,j)%plist(k)%vys=0.0_dp

    !             heff=h1

    !             if ((dpcell(i,j)%list(k)%count/=0) &
    !             .and.(dpcell(i,j)%pplist(k)%gradvx>=0.80_dp)) then
    !                 if ((dpcell(i,j)%plist(k)%free).or.(dpcell(i,j)%plist(k)%vicinity)) then
    !                 t1=0.0_dp
    !                 t2=0.0_dp
    !                 normx=0.0_dp
    !                 normy=0.0_dp
    !                 dpcell(i,j)%plist(k)%xs=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if ((x%tid/=4).and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     normx=normx-(dpcell(i,j)%pplist(k)%coff(1)*&
    !                     Wabx(x,dpcell(i,j)%plist(k), &
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)* &
    !                     Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
    !                     *(x%mass/x%density)

    !                     normy=normy-(dpcell(i,j)%pplist(k)%coff(3)*&
    !                     Wabx(x,dpcell(i,j)%plist(k), &
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)* &
    !                     Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
    !                     *(x%mass/x%density)

    !                     end if
            
    !                     end associate
    !                 end do

    !                 t1=sqrt(normx**2+normy**2)
    !                 normx=normx/t1
    !                 normy=normy/t1

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if ((x%tid/=4).and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-dpcell(i,j)%pplist(k)%Kstar* &
    !                     ((1.0_dp-normx**2)*Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)-normx*normy* &
    !                     Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%xs=min(abs(dpcell(i,j)%plist(k)%xs),maxshift*dl)* &
    !                 dpcell(i,j)%plist(k)%xs/abs(dpcell(i,j)%plist(k)%xs)

    !                 t1=0.0_dp
    !                 t2=0.0_dp
    !                 dpcell(i,j)%plist(k)%ys=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if ((x%tid/=4).and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%ys-dpcell(i,j)%pplist(k)%Kstar*&
    !                     (-normx*normy*Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+(1.0_dp-normy**2)* &
    !                     Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%ys=min(abs(dpcell(i,j)%plist(k)%ys),maxshift*dl)* &
    !                 dpcell(i,j)%plist(k)%ys/abs(dpcell(i,j)%plist(k)%ys)

    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if ((x%tid/=4).and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vxs=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys



    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if ((x%tid/=4).and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1)&
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vys=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys

    !                 dpcell(i,j)%plist(k)%vicinity=.false.

    !                 else


    !                 dpcell(i,j)%plist(k)%xs=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if ((x%tid/=4).and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-dpcell(i,j)%pplist(k)%Kstar*&
    !                     (Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff))&
    !                     *(x%mass/x%density)* &
    !                     (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
    !                     Wab(2*dl,heff))**4)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%xs=min(abs(dpcell(i,j)%plist(k)%xs),maxshift*dl)*&
    !                 dpcell(i,j)%plist(k)%xs/abs(dpcell(i,j)%plist(k)%xs)

    !                 dpcell(i,j)%plist(k)%ys=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if ((x%tid/=4).and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%ys-dpcell(i,j)%pplist(k)%Kstar*&
    !                     (Waby(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff))&
    !                     *(x%mass/x%density)* &
    !                     (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
    !                     Wab(2*dl,heff))**4)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%ys=min(abs(dpcell(i,j)%plist(k)%ys),maxshift*dl)*&
    !                 dpcell(i,j)%plist(k)%ys/abs(dpcell(i,j)%plist(k)%ys)


    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if ((x%tid/=4).and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vxs=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys


    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if ((x%tid/=4).and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vys=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys

    !                 end if
    !             end if

    !             end if

    !         end do
    !         end if
    !         end do
    !         end do
    !     !$omp end do 

    !     ! Shifting hydrodynamic values
    !     !$omp do schedule(dynamic) private(i,j,k) collapse(2)  
    !         do j=sx,ex
    !         do i=sy,ey
    !             if (dpcell(i,j)%ptot/=0) then
    !             do k=1,dpcell(i,j)%ptot
    !                 if (dpcell(i,j)%plist(k)%tid==3) then
    !                 dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+dpcell(i,j)%plist(k)%xs
    !                 dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+dpcell(i,j)%plist(k)%ys
    !                 dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx+dpcell(i,j)%plist(k)%vxs
    !                 dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy+dpcell(i,j)%plist(k)%vys

    !                 end if
    !             end do
    !             end if
    !         end do
    !         end do
    !     !$omp end do   
        
    ! end subroutine opt3_shift

    ! subroutine exp1_shift

    !     implicit none

    !     ! integer :: pint=0

    !     real(dp) :: intx,inty,heff,t1,t2
    !     integer :: i,j,k,m
        
    !     ! Optimized Particle Shifting
    !     !$omp do private(i,k,m,t1,t2,normx,normy,intx,inty,heff) schedule (runtime) collapse(2)
    !         do j=sx,ex
    !         do i=sy,ey
    !         if(dpcell(i,j)%ptot/=0) then

    !         do k=1,dpcell(i,j)%ptot

    !             if (dpcell(i,j)%plist(k)%tid==3) then

    !             dpcell(i,j)%plist(k)%xs=0.0_dp
    !             dpcell(i,j)%plist(k)%ys=0.0_dp
    !             dpcell(i,j)%plist(k)%vxs=0.0_dp
    !             dpcell(i,j)%plist(k)%vys=0.0_dp
    !             intx=0.0_dp
    !             inty=0.0_dp

    !             ! if (dpcell(i,j)%pplist(k)%inpore) then

    !             heff=(max(prrealx,prrealy)*4.80_dp)/(2.0_dp*sqrt(dpcell(i,j)%pplist(k)%porosity))

    !             ! else 

    !             ! heff=(max(prrealx,prrealy)*4.80_dp)/(2.0_dp)

    !             ! end if

    !             if ((dpcell(i,j)%list(k)%count/=0) &
    !             .and.(dpcell(i,j)%pplist(k)%gradvx>=0.70_dp)) then
    !                 if ((dpcell(i,j)%plist(k)%free).or.(dpcell(i,j)%plist(k)%vicinity)) then
    !                 t1=0.0_dp
    !                 t2=0.0_dp
    !                 normx=0.0_dp
    !                 normy=0.0_dp
    !                 dpcell(i,j)%plist(k)%xs=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4.and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     normx=normx-(dpcell(i,j)%pplist(k)%coff(1)*&
    !                     Wabx(x,dpcell(i,j)%plist(k), &
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)* &
    !                     Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
    !                     *(x%mass/x%density)

    !                     normy=normy-(dpcell(i,j)%pplist(k)%coff(3)*&
    !                     Wabx(x,dpcell(i,j)%plist(k), &
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)* &
    !                     Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
    !                     *(x%mass/x%density)

    !                     end if
            
    !                     end associate
    !                 end do

    !                 t1=sqrt(normx**2+normy**2)
    !                 normx=normx/t1
    !                 normy=normy/t1

    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4.and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                         t1=t1+dpcell(i,j)%list(k)%dist(m)/ &
    !                         (dpcell(i,j)%list(k)%count)
                        
    !                         t2=t2+(dpcell(i,j)%plist(k)%x- &
    !                         x%x)/ &
    !                         dpcell(i,j)%list(k)%dist(m)**3

    !                     end if

    !                     end associate
    !                 end do

    !                 intx=beta*umax*dt*(t1**2)*t2

    !                 t1=0.0_dp
    !                 t2=0.0_dp
    !                 if (dpcell(i,j)%list(k)%count/=0) then
    !                 do m=1,dpcell(i,j)%list(k)%count
    !                 associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4.and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     t1=t1+dpcell(i,j)%list(k)%dist(m)/ &
    !                     (dpcell(i,j)%list(k)%count)
                    
    !                     t2=t2+(dpcell(i,j)%plist(k)%y- &
    !                     x%y)/ &
    !                     dpcell(i,j)%list(k)%dist(m)**3

    !                     end if

    !                     end associate
    !                 end do
    !                 end if

    !                 inty=beta*umax*dt*(t1**2)*t2

    !                 dpcell(i,j)%plist(k)%xs=(1-normx**2)*intx-(normx*normy)*inty
    !                 dpcell(i,j)%plist(k)%ys=(1-normy**2)*inty-(normx*normy)*intx


    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4.and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vxs=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys

    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4.and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1)&
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vys=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys

    !                 dpcell(i,j)%plist(k)%vicinity=.false.

    !                 else


    !                 dpcell(i,j)%plist(k)%xs=0.0_dp

    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4.and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                         t1=t1+dpcell(i,j)%list(k)%dist(m)/ &
    !                         (dpcell(i,j)%list(k)%count)
                        
    !                         t2=t2+(dpcell(i,j)%plist(k)%x- &
    !                         x%x)/ &
    !                         dpcell(i,j)%list(k)%dist(m)**3

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%xs=beta*umax*dt*(t1**2)*t2

    !                 dpcell(i,j)%plist(k)%ys=0.0_dp

    !                 t1=0.0_dp
    !                 t2=0.0_dp
    !                 if (dpcell(i,j)%list(k)%count/=0) then
    !                 do m=1,dpcell(i,j)%list(k)%count
    !                 associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4.and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     t1=t1+dpcell(i,j)%list(k)%dist(m)/ &
    !                     (dpcell(i,j)%list(k)%count)
                    
    !                     t2=t2+(dpcell(i,j)%plist(k)%y- &
    !                     x%y)/ &
    !                     dpcell(i,j)%list(k)%dist(m)**3

    !                     end if

    !                     end associate
    !                 end do
    !                 end if

    !                 dpcell(i,j)%plist(k)%ys=beta*umax*dt*(t1**2)*t2


    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4.and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vxs=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys


    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4.and. &
    !                     (dpcell(i,j)%list(k)%dist(m)<=(2*heff))) then

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     end if

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vys=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys

    !                 end if
    !             end if

    !             end if

    !         end do
    !         end if
    !         end do
    !         end do
    !     !$omp end do 

    !     ! Shifting hydrodynamic values
    !     !$omp do schedule(dynamic) private(i,j,k) collapse(2)  
    !         do j=sx,ex
    !         do i=sy,ey
    !             if (dpcell(i,j)%ptot/=0) then
    !             do k=1,dpcell(i,j)%ptot
    !                 if (dpcell(i,j)%plist(k)%tid==3) then
    !                 dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+dpcell(i,j)%plist(k)%xs
    !                 dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+dpcell(i,j)%plist(k)%ys
    !                 ! dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx+dpcell(i,j)%plist(k)%vxs
    !                 ! dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy+dpcell(i,j)%plist(k)%vys

    !                 end if
    !             end do
    !             end if
    !         end do
    !         end do
    !     !$omp end do   
        
    ! end subroutine exp1_shift

    ! subroutine opt1_shift

    !     implicit none

    !     ! integer :: pint=0
    !     integer :: i,j,k,m
    !     real(dp) :: A=2.0_dp,heff,frac=1.0_dp,t1,t2
        
    !     ! Optimized Particle Shifting
    !     !$omp do private(i,k,m,t1,t2,normx,normy,heff) schedule (runtime) collapse(2)
    !         do j=sx,ex
    !         do i=sy,ey
    !         if(dpcell(i,j)%ptot/=0) then

    !         do k=1,dpcell(i,j)%ptot

    !             if (dpcell(i,j)%plist(k)%tid>2) then

    !             dpcell(i,j)%plist(k)%xs=0.0_dp
    !             dpcell(i,j)%plist(k)%ys=0.0_dp
    !             dpcell(i,j)%plist(k)%vxs=0.0_dp
    !             dpcell(i,j)%plist(k)%vys=0.0_dp

    !             ! if (dpcell(i,j)%pplist(k)%inpore) then

    !             heff=h1

    !             ! else 

    !             ! heff=(max(prrealx,prrealy)*4.80_dp)/(2.0_dp)

    !             ! end if

    !             if ((dpcell(i,j)%list(k)%count/=0) &
    !             .and.(dpcell(i,j)%pplist(k)%gradvx>=0.70_dp)) then
    !                 if ((dpcell(i,j)%plist(k)%free).or.(dpcell(i,j)%plist(k)%vicinity)) then
    !                 t1=0.0_dp
    !                 t2=0.0_dp
    !                 normx=0.0_dp
    !                 normy=0.0_dp
    !                 dpcell(i,j)%plist(k)%xs=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     normx=normx-(dpcell(i,j)%pplist(k)%coff(1)*&
    !                     Wabx(x,dpcell(i,j)%plist(k), &
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)* &
    !                     Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
    !                     *(x%mass/x%density)

    !                     normy=normy-(dpcell(i,j)%pplist(k)%coff(3)*&
    !                     Wabx(x,dpcell(i,j)%plist(k), &
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)* &
    !                     Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
    !                     *(x%mass/x%density)
            
    !                     end associate
    !                 end do

    !                 t1=sqrt(normx**2+normy**2)
    !                 normx=normx/t1
    !                 normy=normy/t1

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     ! dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-csh*(heff)* &
    !                     ! sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)*dt* &                        
    !                     ! ((1-normx**2)*Wabx(x,dpcell(i,j)%plist(k),&
    !                     ! dpcell(i,j)%list(k)%dist(m),heff)-normx*normy* &
    !                     ! Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
    !                     ! *(x%mass/x%density)

    !                     dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-A*frac*(heff)*&
    !                     sqrt((dpcell(i,j)%plist(k)%vx)**2+(dpcell(i,j)%plist(k)%vy)**2)*dt* &                        
    !                     (Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff))&
    !                     *(x%mass/x%density)

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%xs=min(abs(dpcell(i,j)%plist(k)%xs),maxshift*dl)* &
    !                 dpcell(i,j)%plist(k)%xs/abs(dpcell(i,j)%plist(k)%xs)

    !                 t1=0.0_dp
    !                 t2=0.0_dp
    !                 dpcell(i,j)%plist(k)%ys=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     ! dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%ys-csh*(heff)* &
    !                     ! sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)*dt* &                    
    !                     ! (-normx*normy*Wabx(x,dpcell(i,j)%plist(k),&
    !                     ! dpcell(i,j)%list(k)%dist(m),heff)+(1-normy**2)* &
    !                     ! Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))&
    !                     ! *(x%mass/x%density)

    !                     dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-A*frac*(heff)*&
    !                     sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)*dt* &                        
    !                     (Waby(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff))&
    !                     *(x%mass/x%density)


    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%ys=min(abs(dpcell(i,j)%plist(k)%ys),maxshift*dl)* &
    !                 dpcell(i,j)%plist(k)%ys/abs(dpcell(i,j)%plist(k)%ys)

    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1)*Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3)*Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vxs=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys


    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))


    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1)*Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3)*Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vys=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys

    !                 dpcell(i,j)%plist(k)%vicinity=.false.

    !                 else


    !                 dpcell(i,j)%plist(k)%xs=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))


    !                     ! dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-csh*(heff)* &
    !                     ! sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)*dt* &
    !                     ! (Wabx(x,dpcell(i,j)%plist(k),&
    !                     ! dpcell(i,j)%list(k)%dist(m),heff))&
    !                     ! *(x%mass/x%density)* &
    !                     ! (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
    !                     ! Wab(2*dl,heff))**4)

    !                     dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-A*(heff)*&
    !                     sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)*dt* &                        
    !                     (Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff))&
    !                     *(x%mass/x%density)


    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%xs=min(abs(dpcell(i,j)%plist(k)%xs),maxshift*dl)*&
    !                 dpcell(i,j)%plist(k)%xs/abs(dpcell(i,j)%plist(k)%xs)

    !                 dpcell(i,j)%plist(k)%ys=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))


    !                     ! dpcell(i,j)%plist(k)%ys=dpcell(i,j)%plist(k)%ys-csh*(heff)* &
    !                     ! sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)*dt* &
    !                     ! (Waby(x,dpcell(i,j)%plist(k),&
    !                     ! dpcell(i,j)%list(k)%dist(m),heff))&
    !                     ! *(x%mass/x%density)* &
    !                     ! (1.0_dp+0.20_dp*pint*(Wab(dpcell(i,j)%list(k)%dist(m),heff)/ &
    !                     ! Wab(2*dl,heff))**4)

    !                     dpcell(i,j)%plist(k)%xs=dpcell(i,j)%plist(k)%xs-A*(heff)*&
    !                     sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)*dt* &                        
    !                     (Waby(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff))&
    !                     *(x%mass/x%density)


    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%ys=min(abs(dpcell(i,j)%plist(k)%ys),maxshift*dl)*&
    !                 dpcell(i,j)%plist(k)%ys/abs(dpcell(i,j)%plist(k)%ys)

    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1)*Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3)*Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vxs=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys


    !                 t1=0.0_dp
    !                 t2=0.0_dp

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1)*Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3)*Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),heff)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),heff))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vys=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys

    !                 end if
    !             end if

    !             end if

    !         end do
    !         end if
    !         end do
    !         end do
    !     !$omp end do 

    !     ! Shifting hydrodynamic values
    !     !$omp do schedule(dynamic) private(i,j,k) collapse(2)  
    !         do j=sx,ex
    !         do i=sy,ey
    !             if (dpcell(i,j)%ptot/=0) then
    !             do k=1,dpcell(i,j)%ptot
    !                 if (dpcell(i,j)%plist(k)%tid>2) then
    !                 dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+dpcell(i,j)%plist(k)%xs
    !                 dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+dpcell(i,j)%plist(k)%ys
    !                 dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx+dpcell(i,j)%plist(k)%vxs/dpcell(i,j)%pplist(k)%porosity
    !                 dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy+dpcell(i,j)%plist(k)%vys/dpcell(i,j)%pplist(k)%porosity
    !                 end if
    !             end do
    !             end if
    !         end do
    !         end do
    !     !$omp end do   
        
    ! end subroutine opt1_shift

    ! subroutine coll_pos_shift
    !     implicit none
        
    !     real(dp) :: t1,t2
    !     integer :: i,j,k,m
    !     real(dp),parameter :: cr=0.0_dp,disfac=1.0_dp
    !     real(dp) :: dcol
    !     ! Explicit particle shifting
    !     !$omp do private(i,j,k,m,t1,t2,dcol) schedule (runtime) collapse(2)
    !         do j=sx,ex
    !             do i=sy,ey
    !                 if(dpcell(i,j)%ptot/=0) then
                    
    !                 do k=1,dpcell(i,j)%ptot
    !                     dpcell(i,j)%plist(k)%ys=0.0_dp
    !                     dpcell(i,j)%plist(k)%xs=0.0_dp

    !                     if ((dpcell(i,j)%list(k)%count/=0) &
    !                     .and.(dpcell(i,j)%pplist(k)%gradvx>=0.70_dp) &
    !                     .and.(dpcell(i,j)%plist(k)%tid==3).and. &
    !                     (.not.(dpcell(i,j)%plist(k)%buffer))) then
    !                     dpcell(i,j)%plist(k)%ys=0.0_dp
    !                     dpcell(i,j)%plist(k)%xs=0.0_dp
    !                     t1=0.0_dp
    !                     t2=0.0_dp

    !                     do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                         y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                         pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                         dcol=0.5_dp*2*prrealx*disfac* &
    !                         ((sqrt(dpcell(i,j)%pplist(k)%porosity)**(-1))+ &
    !                         (sqrt(y%porosity)**(-1)))

    !                         t1=t1+(x%mass)
                        
    !                         t2=t2+merge(-1.0_dp*(x%mass*(1.0_dp+cr)* &
    !                         dpcell(i,j)%list(k)%dist(m)-dcol)*((dpcell(i,j)%plist(k)%x- &
    !                         x%x)/dpcell(i,j)%list(k)%dist(m)),0.0_dp, &
    !                         (dpcell(i,j)%list(k)%dist(m)<dcol))

    !                         end associate
    !                     end do

    !                     dpcell(i,j)%plist(k)%xs=t2/(dpcell(i,j)%plist(k)%mass+t1)

    !                     t1=0.0_dp
    !                     t2=0.0_dp
    !                     if (dpcell(i,j)%list(k)%count/=0) then
    !                     do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                         y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                         pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                         dcol=0.5_dp*2*prrealy*disfac* &
    !                         ((sqrt(dpcell(i,j)%pplist(k)%porosity)**(-1))+ &
    !                         (sqrt(y%porosity)**(-1)))

    !                         t1=t1+(x%mass)
                        
    !                         t2=t2+merge(-1.0_dp*(x%mass*(1.0_dp+cr)* &
    !                         dpcell(i,j)%list(k)%dist(m)-dcol)*((dpcell(i,j)%plist(k)%y- &
    !                         x%y)/dpcell(i,j)%list(k)%dist(m)),0.0_dp, &
    !                         (dpcell(i,j)%list(k)%dist(m)<dcol))

    !                         end associate
    !                     end do
    !                     end if

    !                     dpcell(i,j)%plist(k)%ys=t2/(dpcell(i,j)%plist(k)%mass+t1)

    !                     end if
    !                 end do  
                                    
    !                 end if
    !             end do
    !         end do
    !     !$omp end do

    !     ! Shifting hydrodynamic values
    !     !$omp do schedule(dynamic) private(i,j,k) collapse(2)
    !         do j=sx,ex
    !         do i=sy,ey
    !             if (dpcell(i,j)%ptot/=0) then
    !             do k=1,dpcell(i,j)%ptot
    !                 if ((dpcell(i,j)%pplist(k)%gradvx>=0.70_dp) &
    !                 .and.(dpcell(i,j)%plist(k)%tid==3).and. &
    !                 (.not.(dpcell(i,j)%plist(k)%buffer))) then
    !                 dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+dpcell(i,j)%plist(k)%xs
    !                 dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+dpcell(i,j)%plist(k)%ys

    !                 end if
    !             end do
    !             end if
    !         end do
    !         end do
    !     !$omp end do   
        
    ! end subroutine coll_pos_shift

    ! subroutine implicit_shift

    !     use functions
    !     use solver

    !     implicit none

    !     ! real(dp),dimension(:) :: pguess(finmax)
    !     real(dp) :: heff,frac=1.0_dp,lamp,t1,t2
    !     integer :: i,j,k,m

    !     heff=h1

    !     ! pguess=0.0_dp


    !     !Density for internal particles
    !     !$omp parallel default(shared)
    !     !$omp do schedule(runtime) private(k,m,t1,i,j) collapse(2)
    !     do j=sx,ex
    !         do i=sy,ey
    !         if (dpcell(i,j)%ptot/=0) then
                
    !             do k=1,dpcell(i,j)%ptot
    !             if ((dpcell(i,j)%plist(k)%tid>2)) then
    !                 if(dpcell(i,j)%list(k)%count/=0) then
    !                     dpcell(i,j)%pplist(k)%tden=0.0_dp
    !                     ! dpcell(i,j)%pplist(k)%Kstar=0.0_dp
    !                     do m=1,dpcell(i,j)%list(k)%count                           
    !                         associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                             y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                             pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                             if (x%tid/=4) then
                            
    !                             dpcell(i,j)%pplist(k)%tden=dpcell(i,j)%pplist(k)%tden+ &
    !                             Wab(dpcell(i,j)%list(k)%dist(m),h1)*x%mass/y%porosity

    !                             end if

    !                         end associate
    !                     end do
    !                 end if

    !                 dpcell(i,j)%pplist(k)%tden=dpcell(i,j)%pplist(k)%tden+&
    !                 dpcell(i,j)%plist(k)%mass*Wo(h1)/dpcell(i,j)%pplist(k)%porosity

    !             elseif ((dpcell(i,j)%plist(k)%tid<=2)) then

    !             dpcell(i,j)%pplist(k)%tden=rho

    !             end if
    !             end do
                
    !         end if
    !         end do
    !     end do
    !     !$omp end do 

    !     !$omp do schedule(runtime) private(k,m,t1,i,j) collapse(2)
    !     do j=sx,ex
    !         do i=sy,ey
    !         if (dpcell(i,j)%ptot/=0) then
                
    !             do k=1,dpcell(i,j)%ptot
    !             if ((dpcell(i,j)%plist(k)%tid>2)) then
    !                 if(dpcell(i,j)%list(k)%count/=0) then
    !                     ! dpcell(i,j)%pplist(k)%tden=0.0_dp
    !                     dpcell(i,j)%pplist(k)%Kstar=0.0_dp
    !                     do m=1,dpcell(i,j)%list(k)%count                           
    !                         associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                             y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                             pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                             if (x%tid/=4) then

    !                             dpcell(i,j)%pplist(k)%Kstar=dpcell(i,j)%pplist(k)%Kstar+ &
    !                             Wab(dpcell(i,j)%list(k)%dist(m),h1)* &
    !                             x%mass/y%tden

    !                             end if

    !                         end associate
    !                     end do
    !                 end if

    !                 dpcell(i,j)%pplist(k)%Kstar=dpcell(i,j)%pplist(k)%Kstar+ &
    !                 dpcell(i,j)%plist(k)%mass*Wo(h1)/dpcell(i,j)%pplist(k)%tden

    !             elseif ((dpcell(i,j)%plist(k)%tid<=2)) then

    !             dpcell(i,j)%pplist(k)%Kstar=1.0_dp

    !             end if
    !             end do
                
    !         end if
    !         end do
    !     end do
    !     !$omp end do 

    !     ! Preparing the coff matrix for fluid part in CSR
    !     !$omp do schedule(runtime) &
    !     !$omp private(m,t1,t2,k,i,j,lamp) collapse(2)
    !     do j=sx,ex
    !         do i=sy,ey
    !         if (dpcell(i,j)%ptot/=0) then

    !             do k=1,dpcell(i,j)%ptot
                        
    !             if (dpcell(i,j)%plist(k)%tid/=4) then
    !             associate(pos=>dpcell(i,j)%plist(k)%matid,num=>dpcell(i,j)%list(k)%count)
    !                 t1=0.0_dp
    !                 t2=0.0_dp
    !                 fmatrix(pos)%sz=num+1
    !                 fmatrix(pos)%val(:)=0.0_dp
    !                 fmatrix(pos)%col(:)=pos
    !                 fvec(pos)=0.0_dp
    !                 num2=dpcell(i,j)%list(k)%count
    !                 if((num/=0)) then   !.or.(dpcell(i,j)%plist(k)%free/=1)
    !                 do m=1,num2
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                         if (x%tid/=4) then

    !                             if (dpcell(i,j)%plist(k)%free) then

    !                                 lamp=0.0_dp!(1.0_dp+x%density/(dpcell(y,x)%pplist(k)%porosity*rhomax))
    
    !                             elseif((pll<dpcell(i,j)%pplist(k)%gradvx).and.(dpcell(i,j)%pplist(k)%gradvx<pul)) then
    
    !                             lamp=0.50_dp*(1.0_dp-cos(22.0_dp*(dpcell(i,j)%pplist(k)%gradvx-pll)/ &
    !                                             ((pul-pll)*7.0_dp)))
    !                             else 
    
    !                                 lamp =1.0_dp
    
    !                             end if

    !                         if (dpcell(i,j)%plist(k)%tid>2) then

    !                             t1=2*x%mass*(dpcell(i,j)%pplist(k)%porosity**2)*(&
    !                             Wabx(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
    !                             *(dpcell(i,j)%plist(k)%x-x%x)/ &
    !                             (((dist(x,dpcell(i,j)%plist(k))**2)+lam)* &
    !                             (y%tden*dpcell(i,j)%pplist(k)%tden))

    !                             t2=2*x%mass*(dpcell(i,j)%pplist(k)%porosity**2)*(&
    !                             Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
    !                             *(dpcell(i,j)%plist(k)%y-x%y)/ &
    !                             (((dist(x,dpcell(i,j)%plist(k))**2)+lam)* &
    !                             (y%tden*dpcell(i,j)%pplist(k)%tden))                               

    !                             fmatrix(pos)%val(m)=(-(t1+t2))*lamp
    !                             fmatrix(pos)%col(m)=x%matid

    !                             fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(t1+t2)


    !                             ! fvec(pos)=(1.0_dp-dpcell(i,j)%pplist(k)%Kstar)*lamp/real(dt**2,dp) !lamp

    !                             if ((dpcell(i,j)%plist(k)%free).and.(dpcell(i,j)%plist(k)%near_boun)) then


    !                                 fmatrix(pos)%val(m)=(-(t1+t2))
    !                                 fvec(pos)=fvec(pos)+frac*csh*(heff**2)*&
    !                                 (Wabx(x,dpcell(i,j)%plist(k),&
    !                                 dpcell(i,j)%list(k)%dist(m),heff))&
    !                                 *(x%mass/x%density)/&
    !                                 ((dist(x,dpcell(i,j)%plist(k)))*dt+lam*dt) + &
    !                                 frac*csh*(heff**2)*&
    !                                 (Waby(x,dpcell(i,j)%plist(k),&
    !                                 dpcell(i,j)%list(k)%dist(m),heff))&
    !                                 *(x%mass/x%density)/&
    !                                 ((dist(x,dpcell(i,j)%plist(k)))*(dt**2)+lam*(dt**2))

    !                             else

    !                                 fvec(pos)=(1.0_dp-dpcell(i,j)%pplist(k)%Kstar)*lamp/real(dt**2,dp)

    !                             end if



    !                         elseif (dpcell(i,j)%plist(k)%tid==2) then

    !                             t1=2*x%mass*(&
    !                             Wabx(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
    !                             *(dpcell(i,j)%plist(k)%x-x%x)/ &
    !                             (((dist(x,dpcell(i,j)%plist(k))**2)+lam)* &
    !                             (x%density*dpcell(i,j)%pplist(k)%tden))

    !                             t2=2*x%mass*(&
    !                             Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
    !                             *(dpcell(i,j)%plist(k)%y-x%y)/ &
    !                             (((dist(x,dpcell(i,j)%plist(k))**2)+lam)* &
    !                             (x%density*dpcell(i,j)%pplist(k)%tden))

                                

    !                             if (x%tid<=2) then

    !                             fmatrix(pos)%val(m)=-(t1+t2)!*y%lamp
    !                                 fmatrix(pos)%col(m)=x%matid
    !                                 fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(t1+t2)!*y%lamp

    !                                 fvec(pos)=(1.0_dp-dpcell(i,j)%pplist(k)%Kstar)/real(dt**2,dp)



    !                             end if



    !                             ! fmatrix(pos)%val(1)=-1.0_dp!(t1+t2)
    !                             ! fmatrix(pos)%col(1)=dpcell(i,j)%plist(k)%wall(1)!x%matid
    !                             ! fmatrix(pos)%val(num+1)=1.0_dp!fmatrix(pos)%val(num+1)+(t1+t2)


    !                         elseif (dpcell(i,j)%plist(k)%tid==1) then

    !                             if (x%tid==3) then

    !                             t1=2*dpcell(i,j)%plist(k)%mass*(y%porosity**2)*(&
    !                             Wabx(dpcell(i,j)%plist(k),x,dpcell(i,j)%list(k)%dist(m),h1)) &
    !                             *(-dpcell(i,j)%plist(k)%x+x%x)/ &
    !                             (((dist(x,dpcell(i,j)%plist(k))**2)+lam)* &
    !                             (dpcell(i,j)%pplist(k)%tden*y%tden))

    !                             t2=2*dpcell(i,j)%plist(k)%mass*(y%porosity**2)*(&
    !                             Waby(dpcell(i,j)%plist(k),x,dpcell(i,j)%list(k)%dist(m),h1)) &
    !                             *(-dpcell(i,j)%plist(k)%y+x%y)/ &
    !                             (((dist(x,dpcell(i,j)%plist(k))**2)+lam)* &
    !                             (dpcell(i,j)%pplist(k)%tden*y%tden))

    !                             ! t1=2*x%mass*(dpcell(i,j)%pplist(k)%porosity**2)*(&
    !                             ! Wabx(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
    !                             ! *(dpcell(i,j)%plist(k)%x-x%x)/ &
    !                             ! (((dist(x,dpcell(i,j)%plist(k))**2)+lam)* &
    !                             ! (x%density*dpcell(i,j)%pplist(k)%tden))

    !                             ! t2=2*x%mass*(dpcell(i,j)%pplist(k)%porosity**2)*(&
    !                             ! Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
    !                             ! *(dpcell(i,j)%plist(k)%y-x%y)/ &
    !                             ! (((dist(x,dpcell(i,j)%plist(k))**2)+lam)* &
    !                             ! (x%density*dpcell(i,j)%pplist(k)%tden))

    !                                 fmatrix(pos)%val(m)=-(t1+t2)!*y%lamp
    !                             fmatrix(pos)%col(m)=x%matid
    !                             fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(t1+t2)!*y%lamp

    !                                 fvec(pos)=(1.0_dp-dpcell(i,j)%pplist(k)%Kstar)/real(dt**2,dp)

                                
    !                             elseif (x%tid<=2) then

    !                                 t1=2*x%mass*(&
    !                                 Wabx(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
    !                                 *(dpcell(i,j)%plist(k)%x-x%x)/ &
    !                                 (((dist(x,dpcell(i,j)%plist(k))**2)+lam)* &
    !                                 (y%tden*dpcell(i,j)%pplist(k)%tden))

    !                                 t2=2*x%mass*(&
    !                                 Waby(x,dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)) &
    !                                 *(dpcell(i,j)%plist(k)%y-x%y)/ &
    !                                 (((dist(x,dpcell(i,j)%plist(k))**2)+lam)* &
    !                                 (y%tden*dpcell(i,j)%pplist(k)%tden))

    !                                 fmatrix(pos)%val(m)=-(t1+t2)!*y%lamp
    !                                 fmatrix(pos)%col(m)=x%matid
    !                                 fmatrix(pos)%val(num+1)=fmatrix(pos)%val(num+1)+(t1+t2)!*y%lamp

    !                                 fvec(pos)=(1.0_dp-dpcell(i,j)%pplist(k)%Kstar)/real(dt**2,dp)
    !                                 ! fvec(pos)=fvec(pos)+(t1+t2)/real(dt,dp)  !lamp
                                    

    !                             end if
                            
    !                         end if

    !                     end if


    !                     end associate
    !                 end do
    !                 end if

    !                 fmatrix(pos)%col(num+1)=pos
    !                 if ((num==0)) then  !.or.( dpcell(i,j)%plist(k)%free==1)
    !                 fmatrix(pos)%val(num+1)=1.0_dp
    !                 end if

    !                 t1=fmatrix(pos)%val(num+1)

    !                 ! if (dpcell(i,j)%plist(k)%tid/=2) then

    !                 fvec(pos)=fvec(pos)/real(t1,dp)   ! Jacobi Preconditioning
                    
    !                 do m=1,ceiling(fac2*fplistmax)
    !                 fmatrix(pos)%val(m)=fmatrix(pos)%val(m)/real(t1,dp) !Jacobi Preconditioning
    !                 end do 

    !                 ! end if


    !             end associate

    !             end if

    !             end do

    !         end if            
    !         end do
    !     end do
    !     !$omp end do
    !     !$omp end parallel

    !     call pformat(fmatrix,fval,frow,fcol,finmax)
    !     call bicgstab(tl,pguess,finmax,fval,frow,fcol,fvec,fsol)

    !     ! Assigning pressures
    !     !$omp parallel default(shared)
    !     !$omp do schedule(runtime) private(i,k) collapse(2)       
    !         do j=sx,ex
    !             do i=sy,ey
    !                 if (dpcell(i,j)%ptot/=0) then

    !                     do k=1,dpcell(i,j)%ptot

    !                     if (dpcell(i,j)%plist(k)%tid/=4) then

    !                     dpcell(i,j)%pplist(k)%pshift=fsol(dpcell(i,j)%plist(k)%matid)

    !                     if (dpcell(i,j)%plist(k)%tid<=2) then

    !                         dpcell(i,j)%pplist(k)%pshift=dpcell(i,j)%pplist(k)%pshift&
    !                         +2*rho*abs(g)*min(prrealx,prrealy) 

    !                     end if

    !                     end if

    !                     end do

    !                 end if
    !             end do
    !         end do
    !     !$omp end do 

    !     ! New velocity calculations for particles 
    !     !$omp do private(m,t1,i,k) schedule (runtime) collapse(2)  
    !     do j=sx,ex
    !         do i=sy,ey
    !             if(dpcell(i,j)%ptot/=0) then

    !             do k=1,dpcell(i,j)%ptot
    !                 if (dpcell(i,j)%plist(k)%tid==3) then
    !                 !New fluid particle velocities
    !                 dpcell(i,j)%plist(k)%vxs=0.0d0
    !                 dpcell(i,j)%plist(k)%vys=0.0d0
    !                 if (dpcell(i,j)%list(k)%count/=0) then
    !                 do m=1,dpcell(i,j)%list(k)%count
    !                 associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4) then

    !                     t1=((x%mass*dpcell(i,j)%pplist(k)%porosity/y%tden)*&
    !                     ((y%pshift)+(dpcell(i,j)%pplist(k)%pshift))* &
    !                     (dpcell(i,j)%pplist(k)%coff(1) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))/((dpcell(i,j)%pplist(k)%tden))

    !                     ! t1=((x%mass*dpcell(i,j)%pplist(k)%porosity/y%tden)*&
    !                     ! ((y%pshift)+(dpcell(i,j)%pplist(k)%pshift))* &
    !                     ! (Wabx(x,dpcell(i,j)%plist(k),&
    !                     ! dpcell(i,j)%list(k)%dist(m),h1)))/((dpcell(i,j)%pplist(k)%tden))


    !                     dpcell(i,j)%plist(k)%vxs=dpcell(i,j)%plist(k)%vxs-t1*dt


    !                     t1=((x%mass*dpcell(i,j)%pplist(k)%porosity/y%tden)*&
    !                     ((y%pshift)+(dpcell(i,j)%pplist(k)%pshift))* &
    !                     (dpcell(i,j)%pplist(k)%coff(3) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1)))/((dpcell(i,j)%pplist(k)%tden))

    !                     ! t1=((x%mass*dpcell(i,j)%pplist(k)%porosity/y%tden)*&
    !                     ! ((y%pshift)+(dpcell(i,j)%pplist(k)%pshift))* &
    !                     ! (Waby(x,dpcell(i,j)%plist(k),&
    !                     ! dpcell(i,j)%list(k)%dist(m),h1)))/((dpcell(i,j)%pplist(k)%tden))

    !                     dpcell(i,j)%plist(k)%vys=dpcell(i,j)%plist(k)%vys-t1*dt

    !                     end if

    !                     end associate
    !                 end do
    !                 end if

    !                 end if

    !             end do

    !             end if
    !         end do
    !     end do
    !     !$omp end do  

    !     ! New position calculations for fluid particles
    !     !$omp do schedule (runtime) private(i,k,m,t1,t2) collapse(2)     
    !     do j=sx,ex
    !         do i=sy,ey
    !             if(dpcell(i,j)%ptot/=0) then

    !             do k=1,dpcell(i,j)%ptot

    !                 if (dpcell(i,j)%plist(k)%tid==3) then


    !                 dpcell(i,j)%plist(k)%xs=(dt*dpcell(i,j)%plist(k)%vxs)

    !                 dpcell(i,j)%plist(k)%ys=(dt*dpcell(i,j)%plist(k)%vys)
    !                 dpcell(i,j)%plist(k)%vxs=0.0_dp
    !                 dpcell(i,j)%plist(k)%vys=0.0_dp

    !                 t1=0.0d0
    !                 t2=0.0d0  

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4) then

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                     (x%vx-dpcell(i,j)%plist(k)%vx))&
    !                     *(x%mass/x%density)

    !                     end if
    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vxs=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys


    !                 t1=0.0d0
    !                 t2=0.0d0

    !                 do m=1,dpcell(i,j)%list(k)%count
    !                     associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
    !                     y=>dpcell(i,j)%list(k)%interlist(2,m), &
    !                     pp=>dpcell(i,j)%list(k)%interlist(3,m))

    !                     if (x%tid/=4) then

    !                     t1=t1+((dpcell(i,j)%pplist(k)%coff(1)&
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(2)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     t2=t2+((dpcell(i,j)%pplist(k)%coff(3) &
    !                     *Wabx(x,dpcell(i,j)%plist(k),&
    !                     dpcell(i,j)%list(k)%dist(m),h1)+dpcell(i,j)%pplist(k)%coff(4)*Waby(x,&
    !                     dpcell(i,j)%plist(k),dpcell(i,j)%list(k)%dist(m),h1))*&
    !                     (x%vy-dpcell(i,j)%plist(k)%vy))&
    !                     *(x%mass/x%density)

    !                     end if
    !                     end associate
    !                 end do

    !                 dpcell(i,j)%plist(k)%vys=t1*dpcell(i,j)%plist(k)%xs &
    !                                                 +t2*dpcell(i,j)%plist(k)%ys

    !                 end if

    !             end do

    !             end if
    !         end do
    !     end do
    !     !$omp end do

    !     ! Shifting hydrodynamic values
    !     !$omp do schedule(dynamic) private(i,j,k) collapse(2)  
    !         do j=sx,ex
    !         do i=sy,ey
    !             if (dpcell(i,j)%ptot/=0) then
    !             do k=1,dpcell(i,j)%ptot
    !                 if (dpcell(i,j)%plist(k)%tid==3) then
    !                 dpcell(i,j)%plist(k)%x=dpcell(i,j)%plist(k)%x+dpcell(i,j)%plist(k)%xs
    !                 dpcell(i,j)%plist(k)%y=dpcell(i,j)%plist(k)%y+dpcell(i,j)%plist(k)%ys
    !                 dpcell(i,j)%plist(k)%vx=dpcell(i,j)%plist(k)%vx+dpcell(i,j)%plist(k)%vxs
    !                 dpcell(i,j)%plist(k)%vy=dpcell(i,j)%plist(k)%vy+dpcell(i,j)%plist(k)%vys

    !                 end if
    !             end do
    !             end if
    !         end do
    !         end do
    !     !$omp end do
    !     !$omp end parallel         
        
    ! end subroutine implicit_shift
    

end module part_shift
