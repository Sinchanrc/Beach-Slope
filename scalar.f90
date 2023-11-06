module scalar

    use particle
    use initialize
    use domain
    use kernel

    implicit none


    contains


! subroutine scalart

!     real(dp) :: t1

!     integer :: i,j,k,m

!     !$omp do schedule (runtime) private(m,k,i,j,t1) collapse(2)
!         do j=sx,ex
!             do i=sy,ey            
        
!             do k=1,dpcell(i,j)%ptot

!             if(dpcell(i,j)%plist(k)%tid==3) then

!             dpcell(i,j)%pplist(k)%cdiff=0.0_dp

!             if (dpcell(i,j)%list(k)%count/=0) then

!             do m=1,dpcell(i,j)%list(k)%count
!                 associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
!                 y=>dpcell(i,j)%list(k)%interlist(2,m), &
!                 pp=>dpcell(i,j)%list(k)%interlist(3,m))

!                 if (dpcell(y,x)%plist(pp)%tid==3) then

!                 if ((dpcell(i,j)%pplist(k)%nut<1e-6).and.(dpcell(y,x)%pplist(pp)%nut<1e-6)) then

!                     t1=Dm+(dpcell(i,j)%pplist(k)%nut+dpcell(y,x)%pplist(pp)%nut)*0.50_dp/tschmidt
!                 else

!                 t1=2*(((Dm+dpcell(i,j)%pplist(k)%nut/tschmidt)*(Dm+dpcell(y,x)%pplist(pp)%nut/tschmidt))/&
!                 ((Dm+dpcell(i,j)%pplist(k)%nut/tschmidt)+(Dm+dpcell(y,x)%pplist(pp)%nut/tschmidt)))

!                 end if

!                 dpcell(i,j)%pplist(k)%cdiff=dpcell(i,j)%pplist(k)%cdiff + &
!                 (2*t1*(dpcell(i,j)%plist(k)%con-dpcell(y,x)%plist(pp)%con)* &
!                 (dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)* &
!                 (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
!                 dpcell(i,j)%list(k)%dist(m),h1))*(dpcell(y,x)%plist(pp)%mass)/&
!                 ((dpcell(i,j)%list(k)%dist(m)**2+lam)*dpcell(y,x)%plist(pp)%density))+ &
!                 (2*t1*(dpcell(i,j)%plist(k)%con-dpcell(y,x)%plist(pp)%con)* &
!                 (dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)* &
!                 (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
!                 dpcell(i,j)%list(k)%dist(m),h1))*(dpcell(y,x)%plist(pp)%mass)/&
!                 ((dpcell(i,j)%list(k)%dist(m)**2+lam)*dpcell(y,x)%plist(pp)%density))

!                 end if
!                 end associate

!             end do

!             end if

!             end if

!             end do

!             end do
!         end do
!     !$omp end do

! end subroutine


! subroutine scalarupdate(stime)

!     use particle

!     implicit none 
!     real(dp),intent(in) :: stime
!     integer :: i,j,k,m

!     !$omp do schedule (runtime) private(k,i,j) collapse(2)
!         do j=sx,ex
!             do i=sy,ey            
!                 do k=1,dpcell(i,j)%ptot

!                 if(dpcell(i,j)%plist(k)%tid==3) then

!                     dpcell(i,j)%plist(k)%con=dpcell(i,j)%plist(k)%con+dpcell(i,j)%pplist(k)%cdiff*stime


!                 end if

!                 end do
!             end do
!         end do
!     !$omp end do

! end subroutine

! subroutine massupdate

!     implicit none 

!     integer :: i,j,k,m

!     !$omp do schedule (runtime) private(k,i,j) collapse(2)
!         do j=sx,ex
!             do i=sy,ey            
!                 do k=1,dpcell(i,j)%ptot

!                 if(dpcell(i,j)%plist(k)%tid==3) then

!                     dpcell(i,j)%plist(k)%density=0.5*(rhomax+rhomin)+dpcell(i,j)%plist(k)%con*(rhomax-rhomin)
!                     dpcell(i,j)%plist(k)%mass=dpcell(i,j)%plist(k)%density*dpcell(i,j)%plist(k)%ovol


!                 end if

!                 end do
!             end do
!         end do
!     !$omp end do

! end subroutine

subroutine scalart

    integer :: i,j,k,m
    real(dp) :: t1,t2,t3

    !$omp do schedule (runtime) private(m,k,i,j,t1,t2,t3) collapse(2)
        do j=sx,ex
            do i=sy,ey            
        
            do k=1,dpcell(i,j)%ptot

            if(dpcell(i,j)%plist(k)%tid==3) then

            dpcell(i,j)%pplist(k)%cdiff=0.0_dp

            if (dpcell(i,j)%list(k)%count/=0) then

                dpcell(i,j)%plist(k)%con=(dpcell(i,j)%plist(k)%mass* &
                dpcell(i,j)%pplist(k)%porosity)/dpcell(i,j)%plist(k)%ovol

            do m=1,dpcell(i,j)%list(k)%count
                associate(x=>dpcell(i,j)%list(k)%interlist(1,m), &
                y=>dpcell(i,j)%list(k)%interlist(2,m), &
                pp=>dpcell(i,j)%list(k)%interlist(3,m))

                if (dpcell(y,x)%plist(pp)%tid==3) then

                    dpcell(y,x)%plist(k)%con=(dpcell(y,x)%plist(pp)%mass* &
                dpcell(y,x)%pplist(pp)%porosity)/dpcell(y,x)%plist(pp)%ovol

                if ((dpcell(i,j)%pplist(k)%nut<1e-6).and.(dpcell(y,x)%pplist(pp)%nut<1e-6)) then

                    t1=Dm+(dpcell(i,j)%pplist(k)%nut+dpcell(y,x)%pplist(pp)%nut)*0.50_dp/tschmidt
                else

                t1=2*(((Dm+dpcell(i,j)%pplist(k)%nut/tschmidt)*(Dm+dpcell(y,x)%pplist(pp)%nut/tschmidt))/&
                ((Dm+dpcell(i,j)%pplist(k)%nut/tschmidt)+(Dm+dpcell(y,x)%pplist(pp)%nut/tschmidt)))

                end if

                t2=2*dpcell(i,j)%pplist(k)%porosity*dpcell(y,x)%pplist(pp)%porosity/&
                    (dpcell(i,j)%pplist(k)%porosity+dpcell(y,x)%pplist(pp)%porosity)
                
                t3=(dpcell(i,j)%plist(k)%mass/dpcell(i,j)%plist(k)%density)**2 + &
                (dpcell(y,x)%plist(pp)%mass/dpcell(y,x)%plist(pp)%density)**2

                dpcell(i,j)%pplist(k)%cdiff=dpcell(i,j)%pplist(k)%cdiff + &
                (2*t1*t2*t3*(dpcell(i,j)%plist(k)%con-dpcell(y,x)%plist(pp)%con)* &
                (dpcell(i,j)%plist(k)%x-dpcell(y,x)%plist(pp)%x)* &
                (Wabx(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                dpcell(i,j)%list(k)%dist(m),h1))/(dpcell(i,j)%list(k)%dist(m)**2+lam))+ &
                (2*t1*t2*t3*(dpcell(i,j)%plist(k)%con-dpcell(y,x)%plist(pp)%con)* &
                (dpcell(i,j)%plist(k)%y-dpcell(y,x)%plist(pp)%y)* &
                (Waby(dpcell(y,x)%plist(pp),dpcell(i,j)%plist(k), &
                dpcell(i,j)%list(k)%dist(m),h1))/(dpcell(i,j)%list(k)%dist(m)**2+lam))

                end if
                end associate

            end do

            end if

            end if

            end do

            end do
        end do
    !$omp end do

end subroutine


subroutine scalarupdate(stime)

    use particle

    implicit none 
    integer :: i,j,k,m
    real(dp) :: t1,t2,t3
    real(dp),intent(in) :: stime

    !$omp do schedule (runtime) private(k,i,j) collapse(2)
        do j=sx,ex
            do i=sy,ey            
                do k=1,dpcell(i,j)%ptot

                if(dpcell(i,j)%plist(k)%tid==3) then

                    ! dpcell(i,j)%plist(k)%con=dpcell(i,j)%plist(k)%con+dpcell(i,j)%pplist(k)%cdiff*stime

                    dpcell(i,j)%plist(k)%mass=dpcell(i,j)%plist(k)%mass+dpcell(i,j)%pplist(k)%cdiff*stime


                end if

                end do
            end do
        end do
    !$omp end do

end subroutine

subroutine massupdate

    implicit none 

    integer :: i,j,k,m
    real(dp) :: t1,t2,t3

    !$omp do schedule (runtime) private(k,i,j) collapse(2)
        do j=sx,ex
            do i=sy,ey            
                do k=1,dpcell(i,j)%ptot

                if(dpcell(i,j)%plist(k)%tid==3) then

                    ! dpcell(i,j)%plist(k)%density=0.5*(rhomax+rhomin)+dpcell(i,j)%plist(k)%con*(rhomax-rhomin)
                    ! dpcell(i,j)%plist(k)%mass=dpcell(i,j)%plist(k)%density*dpcell(i,j)%plist(k)%ovol

                    dpcell(i,j)%plist(k)%oden=dpcell(i,j)%plist(k)%mass/dpcell(i,j)%plist(k)%ovol


                end if

                end do
            end do
        end do
    !$omp end do

end subroutine
    

end module
