module scalar

    use particle
    use initialize
    use domain
    use kernel

    implicit none


    contains

    subroutine conupdate()

        use particle
    
        implicit none 
        integer :: i,j,k
    
        !$omp do schedule (runtime) private(k,i,j) collapse(2)
            do j=sx,ex
                do i=sy,ey            
                    do k=1,dpcell(i,j)%ptot
    
                    if(dpcell(i,j)%plist(k)%tid==3) then
    
    
                        dpcell(i,j)%plist(k)%con=dpcell(i,j)%plist(k)%mass/dpcell(i,j)%plist(k)%ovol
    
    
                    end if
    
                    end do
                end do
            end do
        !$omp end do
    
    end subroutine



subroutine scalart

    integer :: i,j,k,m
    real(dp) :: t1,t2,t3

    !$omp do schedule (runtime) private(m,k,i,j,t1,t2,t3) collapse(2)
        do j=sx,ex
            do i=sy,ey            
        
            do k=1,dpcell(i,j)%ptot

            if(dpcell(i,j)%plist(k)%tid==3) then

            dpcell(i,j)%plist(k)%con=0.0_dp

            dpcell(i,j)%plist(k)%con=dpcell(i,j)%plist(k)%mass/dpcell(i,j)%plist(k)%ovol

            do m=1,dpcell(i,j)%list(k)%count
                associate(x=>dpcell(i,j)%list(k)%nh(m)%part, &
                    y=>dpcell(i,j)%list(k)%pnh(m)%ppart, &
                    z=>dpcell(i,j)%list(k)%klt)

                if (x%tid==3) then

                    x%con=x%mass/x%ovol

                if ((dpcell(i,j)%pplist(k)%nut<1e-6).and.(y%nut<1e-6)) then

                    t1=Dm+(dpcell(i,j)%pplist(k)%nut+y%nut)*0.50_dp/tschmidt
                else

                t1=2*(((Dm+dpcell(i,j)%pplist(k)%nut/tschmidt)*(Dm+y%nut/tschmidt))/&
                ((Dm+dpcell(i,j)%pplist(k)%nut/tschmidt)+(Dm+y%nut/tschmidt)))

                end if

                t2=2*dpcell(i,j)%pplist(k)%porosity*y%porosity/&
                    (dpcell(i,j)%pplist(k)%porosity+y%porosity)
                
                t3=(dpcell(i,j)%plist(k)%ovol**2)*(dpcell(i,j)%pplist(k)%porosity+y%porosity)/&
                (dpcell(i,j)%pplist(k)%porosity*y%porosity)

                dpcell(i,j)%plist(k)%con=dpcell(i,j)%plist(k)%con + &
                (t1*t2*t3*(dpcell(i,j)%plist(k)%con-x%con)* &
                (dpcell(i,j)%plist(k)%x-x%x)*z(1,m)/(dpcell(i,j)%list(k)%dist(m)**2+lam))+ &
                (t1*t2*t3*(dpcell(i,j)%plist(k)%con-x%con)* &
                (dpcell(i,j)%plist(k)%y-x%y)*z(2,m)/(dpcell(i,j)%list(k)%dist(m)**2+lam))

                end if
                end associate

            end do

            end if

            end do

            end do
        end do
    !$omp end do

end subroutine


subroutine scalarupdate(stime)

    use particle

    implicit none 
    integer :: i,j,k
    real(dp),intent(in) :: stime

    !$omp do schedule (runtime) private(k,i,j) collapse(2)
        do j=sx,ex
            do i=sy,ey            
                do k=1,dpcell(i,j)%ptot

                if(dpcell(i,j)%plist(k)%tid==3) then

                    ! dpcell(i,j)%plist(k)%con=dpcell(i,j)%plist(k)%con+dpcell(i,j)%pplist(k)%cdiff*stime

                    dpcell(i,j)%plist(k)%mass=dpcell(i,j)%plist(k)%mass+dpcell(i,j)%plist(k)%con*stime


                end if

                end do
            end do
        end do
    !$omp end do

end subroutine

subroutine massupdate

    implicit none 

    integer :: i,j,k

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
