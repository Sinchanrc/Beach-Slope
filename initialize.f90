include 'mkl_rci.f90'
include 'mkl_blas.f90'

module initialize

use particle
use interactions
use domain

implicit none

public

! L=length of domain, H=height of domain,pr=half fluid particle spacing,br=half boundary particle spacing
! h1=half smoothing length,wc=water column height,wl=length of water column,rho=density
! ,tl=solver tolerance
real(dp),public :: wc=0.650_dp,L=3.10_dp ,H=1.0_dp,&
                    h1,pr=0.0050_dp*0.90_dp,br=0.00450_dp*0.90_dp, & !5.366*0.3 1.2
                    wl=3.10_dp,rho=1000.0_dp,g=-9.810_dp,tl=1e-3,Dm=1e-9,&
                    csh=0.50_dp,rx,ry,lid_driven=0.0010_dp,rhomin=1000.0_dp, &
                    rhomax,atwood=0.25870_dp,tschmidt=1.0_dp,soill=3.150_dp,soilh=0.750_dp,set_ht=0.20_dp,&
                    outlet_ht=0.7_dp,coastal_ht=0.450_dp,por=0.380_dp,bulkden=1800_dp

! dt=time step,t=simulation time
real(dp),public :: dt=0.001,t=0.0_dp,told=0.0_dp,time=7.60_dp,displaytime=0.050_dp,dtsol=0.010_dp !0.00001

! fmass=fluid particle mass,prrealx=half particle spacing x dir,prrealy=half particle spacing y direction
! brrealx,brrealy=half boundary spacing ,lam=coefff to prevent singularity,mu=viscosity,beta=coff to control particle shifting
! sig1,sig2=time step control coefficient
real(dp),public :: fmass=0.0_dp,bfdist=0.0_dp,prrealx=0.0_dp,prrealy=0.0_dp,brrealx=0.0_dp,blen=1.0_dp,distfac=1.0_dp, &
                        brrealy=0.0_dp,lam=0.0_dp,mu=0.0010_dp,beta=0.010_dp,sig1=0.20_dp,&
                        sig2=0.20_dp,delt=0.10_dp,maxshift=0.10_dp,alpha=0.010_dp,dl1=0.0_dp,&
                        solidx=0.0_dp,solidy=0.0_dp,line_grad=15.0_dp,xl,yl,xu,yu

! r=particle shifting value,lamfs=surface tracking coeffincient,umax=max velocity
real(dp),public :: r=0.0_dp,cs=0.150_dp,cv=0.080_dp,maxdivr=2.0_dp, &
                        ce=1.0_dp,dl,lamfs=0.80_dp,umax=0.0_dp,ker=0.0_dp,normx,normy,co,t_gam=7, &
                        pll=1.60_dp,pul=1.80_dp,hfac=4.80_dp,fac2=0.80_dp,numax=0.0_dp !term1=0.0_dp,term2=0.0_dp,

! bcor1,bcor2=domain corner points,ref1,ref2=reference points for mirroring interpolation nodes
! type(corner) :: bcor1,bcor2,ref1,ref2
! real(dp),allocatable,public :: bval(:),bvec(:),bsol(:),bguess(:)

!fval=array for storing nonzeros,fvec=rhs vector for solver,fsol=solution vector,fguess=initial guess vector
real(dp),allocatable,public :: fval(:),fvec(:),fsol(:),fguess(:),ploc(:,:),probedata(:,:),pguess(:)

!cellx=domain cells in x dir,celly=domain cells in y dir,bl=number of boundary layers
!maxdiv= maximun values of surface divergence
integer,public :: cellx=0,celly=0,fplist=0,fplistmax=0,l1=0,bp=0,bl=4,num2,spx,spy, &
                bnx=0,bny=0,count=0,fpx=0,fpy=0,totc=0,incr=4,cout=0,tempct=0, &
                binmax=0,binmin=0,finmax=0,finmin=0,totalct=0,pbno=1,matidct=0,solsteps=0 !,i=0,j=0,k=0,m=0

!frow=vector for carrying row values,fcol=vector for carrying correspomding column values
integer,allocatable,public :: brow(:),bcol(:),frow(:),fcol(:)
type(particles),allocatable,public :: blist(:,:),flist(:,:)
type(matrixsetup),allocatable,public :: bmatrix(:),fmatrix(:)
integer,public ::sx,ex,sy,ey,threads,modifier=1000,fac=1,iter=0,icount=0
integer,parameter,public :: ghost=1 ! 0=> No ghost boundary,1=>ghost + fixed,2=> only ghost boundary
character(len=70),public:: result
type(cell),allocatable,dimension(:,:),public,target :: dpcell
type(particles),save ,public:: bcor1,bcor2,ref1,ref2
logical ::swcsph=.false.,lid=.true.

integer :: pt(64),maxfct=1,mnum=1,mtype=11,phase=33, &
                nrhs=1,iparm(64),msglvl=0,err1
integer,allocatable :: perm(:)

type(pprob),allocatable :: probe(:)

integer :: ipar(128),rci_request,gmresct=0 
real(dp),allocatable :: dpar(:),tmp(:)

end module initialize
