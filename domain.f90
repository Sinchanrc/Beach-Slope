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
!> Defining the cell type the spatial domain is comprised of.
!
! REVISION HISTORY:
! 
! 
!------------------------------------------------------------------------------

module domain

    use particle
    use interactions
    use transfer

    implicit none

    public
    
    type cell

        !plist=> list of all particles in cell
        type(particles),allocatable,dimension(:),public :: plist,ftemp,porlist
        !Vectors used for particles transfer from 1 cell to the other
        type(trfvec),public ::  tw,te,tn,ts,tne,tse,tnw,tsw
        !list=>list storing particle interactions in the cell for fluid/boundary
        type(verlet),allocatable,public :: list(:),list2(:)
        type(properties),allocatable :: pplist(:)
        integer,dimension(:),allocatable,public :: cellid,exitlist ! Cell index
        !xleft,xright,ytop,ybot=values defining the extent of the cell
        type(buffer),allocatable :: ebuffpt(:)
        real(dp)  :: xleft=0.0_dp,xright=0.0_dp,ytop=0.0_dp,ybot=0.0_dp,maxvel=0.0_dp,maxeddy=0.0_dp
        integer,public :: ptot=0,btot=0,temfct=0,gcount=0,porct=0,elist=0
        logical :: entrybuff=.false.,exitbuff=.false.,cutoff=.false.
        

    end type cell

    
end module domain