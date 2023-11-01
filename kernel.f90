!------------------------------------------------------------------------------
! IITKGP, PhD
!------------------------------------------------------------------------------
!
! MODULE: Kernel and Kernel Derivatives
!
!> @author
!> Sinchan Roy Chowdhury}
!
! DESCRIPTION: 
!> Functions to calculate kernel values and their derivatives for particle pairs
!  and Moving least squares values for ghost boundary
!
! REVISION HISTORY:
! 
! 
!------------------------------------------------------------------------------




module kernel

    use particle
    use domain
    ! use functions

    implicit none
    
    contains
    
    function Wab(radial,smln) result(res2)
        implicit none
        real(dp) :: res2
        real(dp),intent(in) :: smln
        real(dp),intent(in) :: radial 

        res2=0.0_dp

        ! res2=(7.0_dp/(4*3.140_dp*(smln**2)))*(1+(2*radial/smln))* & !Wendland c2
        !         (1.0_dp-(radial/(2*smln)))**4

        ! res2=(9.0_dp/(4*3.140_dp*(smln**2)))*(1+(3*radial/smln)+ &
        ! (35.0_dp/12)*((radial/smln)**2))* & !Wendland c4
        ! (1.0_dp-(radial/(2*smln)))**6

        res2=(2/(3.14*smln**2))*((3*(radial/(4*smln))**2)- &
        (3*radial/(4*smln))+3/4.0d0)          !Quadratic

    end function Wab

    function Wo(smln) result(res2)
        implicit none
        real(dp) :: res2
        real(dp),intent(in) ::smln 

        res2=0.0_dp

        ! res2=(7.0_dp/(4.0_dp*3.140_dp*(smln**2))) !Wendland c2
        ! res2=(9.0_dp/(4.0_dp*3.140_dp*(smln**2))) !Wendland c4

        res2=3/(2*3.14*smln**2) !Quadratic


    end function Wo

    function Wabx(p1,p2,radial,smln) result(res)
        implicit none
        type(particles),intent(in) :: p1,p2
        real(dp) :: res
        real(dp),intent(in) :: smln
        real(dp) ,intent(in):: radial

        res=0.0_dp
        ! res=(35.0_dp*(p1%x-p2%x)*((1.0_dp-radial/(2*smln)))**3)/ & !Wendland c2
        !     (4.0_dp*3.140_dp*smln**4)

        ! res=((p1%x-p2%x)*((1.0_dp-radial/(2*smln)))**5)* & !Wendland c4
        ! ((35.0_dp*radial/(3*(smln**3)))+(14.0_dp/(3*(smln**2))))* &
        ! (9.0_dp/(4.0_dp*3.140_dp*(smln**2)))

        res=(2/(3.14*smln**2))*((3/(8*smln**2))-(3/(4*smln*radial))) &
                *(p2%x-p1%x)                         !Quadratic


    end function Wabx

    function Waby(p1,p2,radial,smln) result(res)
        implicit none
        type(particles),intent(in) :: p1,p2
        real(dp) :: res
        real(dp),intent(in) :: smln
        real(dp) ,intent(in):: radial

        res=0.0_dp
        ! res=(35.0_dp*(p1%y-p2%y)*((1.0_dp-radial/(2*smln)))**3)/ & !Wendland
        !     (4.0_dp*3.140_dp*smln**4)

        ! res=((p1%y-p2%y)*((1.0_dp-radial/(2*smln)))**5)* & !Wendland c4
        ! ((35.0_dp*radial/(3*(smln**3)))+(14.0_dp/(3*(smln**2))))* &
        ! (9.0_dp/(4.0_dp*3.140_dp*(smln**2)))

        res=(2/(3.14*smln**2))*((3/(8*smln**2))-(3/(4*smln*radial))) &
        *(p2%y-p1%y)                         !Quadratic


    end function Waby

end module kernel