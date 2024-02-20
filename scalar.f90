module scalar

    use particle
    use initialize
    use domain
    use kernel

    implicit none


    contains

    subroutine massupdate(stime)

        use particle
    
        implicit none 
        integer :: i,j,k
        real(dp),intent(in) :: stime
    
        !$omp do schedule (runtime) private(k,i,j) collapse(2)
            do j=sx,ex
                do i=sy,ey            
                    do k=1,dpcell(i,j)%ptot
    
                    if((dpcell(i,j)%plist(k)%tid==3).and. &
                    (.not.(dpcell(i,j)%plist(k)%buffer))) then
    
    
                        dpcell(i,j)%plist(k)%con=dpcell(i,j)%plist(k)%con+ &
                        dpcell(i,j)%pplist(k)%cdiff*stime
    
    
                    end if
    
                    end do
                end do
            end do
        !$omp end do
    
    end subroutine

    subroutine scalart

        integer :: i,j,k,m
        real(dp) :: t1x,t1y,t2,t3,con1,con2,Dxi,Dyi,Dxj,Dyj,vi,vj
        real(dp),parameter :: al=0.001_dp,at=0.0001_dp,Dm=1e-9

        !$omp do schedule (runtime) private(m,k,i,j,t1x,t1y,t2,t3,con1,con2,Dxi,Dyi,Dxj,Dyj,vi,vj) collapse(2)
            do j=sx,ex
                do i=sy,ey            
            
                do k=1,dpcell(i,j)%ptot

                if((dpcell(i,j)%plist(k)%tid==3).and. &
                (.not.(dpcell(i,j)%plist(k)%buffer))) then

                dpcell(i,j)%pplist(k)%cdiff=0.0_dp

                con1=dpcell(i,j)%plist(k)%con

                vi=sqrt(dpcell(i,j)%plist(k)%vx**2+dpcell(i,j)%plist(k)%vy**2)

                Dxi=Dm+((al*dpcell(i,j)%plist(k)%vx**2)+(at*dpcell(i,j)%plist(k)%vy**2)+ &
                ((al-at)*dpcell(i,j)%plist(k)%vx*dpcell(i,j)%plist(k)%vy))/(vi*dpcell(i,j)%pplist(k)%porosity)

                Dyi=Dm+((at*dpcell(i,j)%plist(k)%vx**2)+(al*dpcell(i,j)%plist(k)%vy**2)+ &
                ((al-at)*dpcell(i,j)%plist(k)%vx*dpcell(i,j)%plist(k)%vy))/(vi*dpcell(i,j)%pplist(k)%porosity)

                do m=1,dpcell(i,j)%list(k)%count
                    associate(x=>dpcell(i,j)%list(k)%nh(m)%part, &
                        y=>dpcell(i,j)%list(k)%pnh(m)%ppart, &
                        z=>dpcell(i,j)%list(k)%klt)

                    if ((x%tid==3).and.(.not.(x%buffer))) then

                        con2=x%con

                        vj=sqrt(x%vx**2+x%vy**2)

                        Dxj=Dm+((al*x%vx**2)+(at*x%vy**2)+((al-at)*x%vx*x%vy))/(vj*y%porosity)

                        Dyj=Dm+((at*x%vx**2)+(al*x%vy**2)+((al-at)*x%vx*x%vy))/(vj*y%porosity)

                        if ((dpcell(i,j)%pplist(k)%porosity>=0.65_dp).or.(y%porosity>=0.65_dp)) then

                        if ((dpcell(i,j)%pplist(k)%nut<1e-6).and.(y%nut<1e-6)) then

                            t1x=Dm+(dpcell(i,j)%pplist(k)%nut+y%nut)*0.50_dp/tschmidt
                        else

                        t1x=2*(((Dm+dpcell(i,j)%pplist(k)%nut/tschmidt)*(Dm+y%nut/tschmidt))/&
                        ((Dm+dpcell(i,j)%pplist(k)%nut/tschmidt)+(Dm+y%nut/tschmidt)))
                        

                        end if
                        t1y=t1x

                        else

                            t1x=2.0_dp*Dxi*Dxj/(Dxi+Dxj)
                            t1y=2.0_dp*Dyi*Dyj/(Dyi+Dyj)


                        end if

                        ! t1=Dm

                        t2=2*dpcell(i,j)%pplist(k)%porosity*y%porosity/&
                            (dpcell(i,j)%pplist(k)%porosity+y%porosity)
                        
                        t3=(dpcell(i,j)%plist(k)%ovol)*(dpcell(i,j)%pplist(k)%porosity+y%porosity)/&
                        (dpcell(i,j)%pplist(k)%porosity*y%porosity)

                        dpcell(i,j)%pplist(k)%cdiff=dpcell(i,j)%pplist(k)%cdiff + &
                        (t1x*t2*t3*(con1-con2)* &
                        (dpcell(i,j)%plist(k)%x-x%x)*z(1,m)/(dpcell(i,j)%list(k)%dist(m)**2+lam))+ &
                        (t1y*t2*t3*(con1-con2)* &
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

                    if((dpcell(i,j)%plist(k)%tid==3).and. &
                    (.not.(dpcell(i,j)%plist(k)%buffer))) then

                        dpcell(i,j)%plist(k)%con=dpcell(i,j)%plist(k)%con+ &
                        dpcell(i,j)%pplist(k)%cdiff*stime
                        


                    end if

                    end do
                end do
            end do
        !$omp end do

    end subroutine

    subroutine densityupdate

        implicit none 

        integer :: i,j,k

        !$omp do schedule (runtime) private(k,i,j) collapse(2)
            do j=sx,ex
                do i=sy,ey            
                    do k=1,dpcell(i,j)%ptot

                    if((dpcell(i,j)%plist(k)%tid==3).and. &
                    (.not.(dpcell(i,j)%plist(k)%buffer))) then

                        ! dpcell(i,j)%plist(k)%density=0.5*(rhomax+rhomin)+dpcell(i,j)%plist(k)%con*(rhomax-rhomin)
                        ! dpcell(i,j)%plist(k)%mass=dpcell(i,j)%plist(k)%density*dpcell(i,j)%plist(k)%ovol
                        dpcell(i,j)%plist(k)%oden=0.5*(rhomax+rhomin)+dpcell(i,j)%plist(k)%con*(rhomax-rhomin)
                        dpcell(i,j)%plist(k)%mass=dpcell(i,j)%plist(k)%oden*dpcell(i,j)%plist(k)%ovol
                        dpcell(i,j)%plist(k)%density=dpcell(i,j)%plist(k)%oden*dpcell(i,j)%pplist(k)%porosity


                    end if

                    end do
                end do
            end do
        !$omp end do

    end subroutine
    

end module
