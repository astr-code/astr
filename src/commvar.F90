!+---------------------------------------------------------------------+
!| This module is to define common variables.                          |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module commvar
  !
  use commtype
  !
  implicit none
  !
  integer :: ia,ja,ka,im,jm,km,is,ie,js,je,ks,ke
  integer :: hm
  integer :: numq,num_species,ndims,navg,ninit
  character(len=4) :: conschm,difschm
  character(len=64) :: gridfile,solidfile
  !+---------------------+---------------------------------------------+
  !|            ia,ja,ka | the dimension of the entire domain          | 
  !|            im,jm,km | the dimension of the local domain           | 
  !|   is,ie,js,je,ks,ke | start and end nodes number.                 |
  !|               ndims | The dimension of problem not equations.     |
  !|                  nh | level of halo nodes.                        |
  !|                numq | number of independent variables.            |
  !|     conschm,difschm | the scheme of solving convectional and      |
  !|                     | diffusional terms.                          |
  !|         num_species | number of species                           |
  !|                navg | frequency of averaging.                     |
  !|               ninit | initialisation method.                      |
  !|            gridfile | the gridfile.                               |
  !|       solidbodyfile | the file contains solid body geometry .     |
  !+---------------------+---------------------------------------------+
  logical :: lihomo,ljhomo,lkhomo
  logical :: nondimen,diffterm,lfilter,lreadgrid,lwrite,lavg,lfftk,    &
             limmbou
  character(len=3) :: rkscheme
  !+---------------------+---------------------------------------------+
  !| lihomo,ljhomo,lkhomo| to define homogeneous direction.            |
  !|             nondimen| if the equation is non-dimensional          |
  !|             diffterm| if the diffusion terms is solved.           |
  !|                     | .f. means Euler equations are solved.       |
  !|             lfilter | to activate filer flag                      |
  !|              lwrite | write samples or not.                       |
  !|                lavg | average the flow field or not .             |
  !|               lfftk | to use fft in the k direction.              |
  !|             limmbou | to use immersed boundary method             |
  !|            rkscheme | which rk method to use.                     |
  !+---------------------+---------------------------------------------+
  !
  integer :: npdci,npdcj,npdck
  !+---------------------+---------------------------------------------+
  !|               npdci | to control scheme at boundary.              |
  !+---------------------+---------------------------------------------+
  real(8) :: xmax,xmin,ymax,ymin,zmax,zmin,voldom,dxyzmax,dxyzmin
  real(8) :: alfa_filter,bfacmpld
  integer :: kcutoff
  !+---------------------+---------------------------------------------+
  !|                *mix | min coordinates                             |
  !|                *max | max coordinates                             |
  !|         alfa_filter | the parameter to control width of filter.   |
  !|            bfacmpld | b factor for mp-ld scheme.                  |
  !|             kcutoff | cutoff wavenumber when fft used.            |
  !|              voldom | total volume of the domain.                 |
  !|     dxyzmax,dxyzmin | characteristic grid spacing.                |
  !+---------------------+---------------------------------------------+
  integer :: nstep,maxstep,nwrite,nlstep,filenumb,nmonitor
  integer,allocatable :: imon(:,:)
  real(8) :: time,deltat
  real(8) :: uinf,vinf,winf,roinf,pinf,tinf
  real(8) :: ref_t,reynolds,mach,rgas,gamma,prandtl
  real(8) :: const1,const2,const3,const4,const5,const6,const7
  real(8) :: tempconst,tempconst1
  !+---------------------+---------------------------------------------+
  !|               nstep | the total time step number.                 |
  !|             maxstep | the max step to run.                        |
  !|              nwrite | frequency of output                         |
  !|              nlstep | frequency of list flowstate.                |
  !|            nmonitor | number of montors                           |
  !|                imon | monitor coordinates                         |
  !|            filenumb | filenumber                                  |
  !|                time | total time of computation.                  |
  !|              deltat | time step.                                  |
  !|               ref_t | reference temperature.                      |
  !|            reynolds | Reynolds number.                            |
  !|                mach | Mach number.                                |
  !|                rgas | gas constant, p=ro*rgas*T.                  |
  !|               gamma | specific heat ratio.                        |
  !|             prandtl | Prandtl number.                             |
  !|                uinf | infinite velocity u                         |     
  !|                vinf | infinite velocity v                         |     
  !|                winf | infinite velocity w                         |     
  !|               roinf | infinite density                            |  
  !|                tinf | infinite temperature                        |      
  !|                pinf | infinite pressure                           |   
  !+---------------------+---------------------------------------------+
  character(len=8) :: flowtype
  !+---------------------+---------------------------------------------+
  !|            flowtype | to define the type of flow.                 |
  !+---------------------+---------------------------------------------+
  real(8) :: ctime(10)=0.d0
  !+---------------------+---------------------------------------------+
  !|             cputime | computational time statistics.              |
  !+---------------------+---------------------------------------------+
  integer :: bctype(6)
  real(8) :: twall(6)
  real(8) :: force(3)
  !+---------------------+---------------------------------------------+
  !|              bctype | define type of boundary condition.          |
  !|               twall | wall temperature.                           |
  !|               force | body force.                                 |
  !+---------------------+---------------------------------------------+
  logical :: lisponge,ljsponge,lksponge
  integer :: spg_imin,spg_imax,spg_jmin,spg_jmax,spg_kmin,spg_kmax
  !+---------------------+---------------------------------------------+
  !|   spg_imin,spg_imax |                                             |
  !|   spg_jmin,spg_jmax | number of nodes in the sponge layer near    |
  !|   spg_kmin,spg_kmax | each boundary                               |
  !|               twall | wall temperature.                           |
  !|               force | body force.                                 |
  !+---------------------+---------------------------------------------+
  logical :: lchardecomp
  integer :: recon_schem
  !+---------------------+---------------------------------------------+
  !|         lchardecomp | local character decomposition used to not   |
  !|         recon_schem | scheme of reconstruction method.            |
  !+---------------------+---------------------------------------------+
  logical :: lrestart
  !+---------------------+---------------------------------------------+
  !|            lrestart | to assign the start mode. t=restart, f=new  |
  !+---------------------+---------------------------------------------+
  integer :: nsolid
  type(solid),allocatable,target :: immbody(:)
  type(sboun),allocatable,target :: immbnod(:)
  !+---------------------+---------------------------------------------+
  !|           num_solid | number of solid bodies                      |
  !|             immbody | the immersed body.                          |
  !|             immbnod | the boundary nodes of immersed body.        |
  !+---------------------+---------------------------------------------+
  !
  parameter(hm=5)
  !
end module commvar
!+---------------------------------------------------------------------+
!| The end of the module commvar.                                      |
!+---------------------------------------------------------------------+