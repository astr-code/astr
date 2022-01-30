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
#ifdef COMB
  use cantera
#endif
  !
  implicit none
  !
  integer :: ia,ja,ka,im,jm,km,is,ie,js,je,ks,ke
  integer :: hm
  integer :: numq,num_species,ndims,ninit,num_modequ
  character(len=10) :: turbmode
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
  !|               ninit | initialisation method.                      |
  !|            gridfile | the gridfile.                               |
  !|       solidbodyfile | the file contains solid body geometry .     |
  !+---------------------+---------------------------------------------+
  logical :: lihomo,ljhomo,lkhomo
  logical :: nondimen,diffterm,lfilter,lreadgrid,lfftk,limmbou
  character(len=3) :: rkscheme
  !+---------------------+---------------------------------------------+
  !| lihomo,ljhomo,lkhomo| to define homogeneous direction.            |
  !|             nondimen| if the equation is non-dimensional          |
  !|             diffterm| if the diffusion terms is solved.           |
  !|                     | .f. means Euler equations are solved.       |
  !|             lfilter | to activate filer flag                      |
  !|               lfftk | to use fft in the k direction.              |
  !|             limmbou | to use immersed boundary method             |
  !|            rkscheme | which rk method to use.                     |
  !+---------------------+---------------------------------------------+
  !
  logical :: lwsequ,lwslic,lavg
  !+---------------------+---------------------------------------------+
  !|              lwsequ | write flowfield sequence or not.            |
  !|              lwslic | write slices or not.                        |
  !|                lavg | average the flow field or not .             |
  !+---------------------+---------------------------------------------+
  !
  integer :: feqchkpt,feqwsequ,feqslice,feqlist,feqavg
  !+---------------------+---------------------------------------------+
  !|            feqchkpt | frequency of writing checkpoint             |
  !|            feqwsequ | frequency of writing flowfield sequence     |
  !|            feqslice | frequency of writing slices.                |
  !|             feqlist | frequency of listing log.                   |
  !|              feqavg | frequency of averaging flowfield.           |
  !+---------------------+---------------------------------------------+
  !
  integer :: npdci,npdcj,npdck
  !+---------------------+---------------------------------------------+
  !|               npdci | to control scheme at boundary.              |
  !+---------------------+---------------------------------------------+
  real(8) :: xmax,xmin,ymax,ymin,zmax,zmin,voldom,dxyzmax,dxyzmin
  real(8) :: alfa_filter,bfacmpld,shkcrt
  integer :: kcutoff
  !+---------------------+---------------------------------------------+
  !|                *mix | min coordinates                             |
  !|                *max | max coordinates                             |
  !|         alfa_filter | the parameter to control width of filter.   |
  !|            bfacmpld | b factor for MP-LD scheme.                  |
  !|              shkcrt | criteria to identify shock.                 |
  !|             kcutoff | cutoff wavenumber when fft used.            |
  !|              voldom | total volume of the domain.                 |
  !|     dxyzmax,dxyzmin | characteristic grid spacing.                |
  !+---------------------+---------------------------------------------+
  integer :: nstep,nsrpt,maxstep,filenumb,nmonitor,fnumslic
  integer :: islice,jslice,kslice
  integer,allocatable :: imon(:,:)
  real(8) :: time,deltat
  real(8) :: uinf,vinf,winf,roinf,pinf,tinf
  real(8) :: ref_t,reynolds,mach,rgas,gamma,prandtl
  real(8),allocatable :: schmidt(:)
  real(8) :: const1,const2,const3,const4,const5,const6,const7
  real(8) :: tempconst,tempconst1
  real(8),allocatable :: spcinf(:)
  !+---------------------+---------------------------------------------+
  !|               nstep | the total time step number.                 |
  !|             maxstep | the max step to run.                        |
  !|islice,jslice,kslice | position of slices.                         |
  !|            nmonitor | number of montors                           |
  !|                imon | monitor coordinates                         |
  !|            filenumb | filenumber                                  |
  !|            fnumslic | file number for slices                      |
  !|                time | total time of computation.                  |
  !|              deltat | time step.                                  |
  !|               ref_t | reference temperature.                      |
  !|            reynolds | Reynolds number.                            |
  !|                mach | Mach number.                                |
  !|                rgas | gas constant, p=ro*rgas*T.                  |
  !|               gamma | specific heat ratio.                        |
  !|             prandtl | Prandtl number.                             |
  !|             schmidt | the Schmidt number array for each scalar.   |
  !|                uinf | infinite velocity u                         |     
  !|                vinf | infinite velocity v                         |     
  !|                winf | infinite velocity w                         |     
  !|               roinf | infinite density                            |  
  !|                tinf | infinite temperature                        |      
  !|                pinf | infinite pressure                           |   
  !+---------------------+---------------------------------------------+
  character(len=16) :: flowtype
  !+---------------------+---------------------------------------------+
  !|            flowtype | to define the type of flow.                 |
  !+---------------------+---------------------------------------------+
  real(8) :: preptime=0.d0
  real(8) :: ctime(32)=0.d0
  !+---------------------+---------------------------------------------+
  !|             cputime | computational time statistics.              |
  !+---------------------+---------------------------------------------+
  real(8) :: force(3)
  !+---------------------+---------------------------------------------+
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
  integer :: nsolid,imbroot
  type(solid),allocatable,target :: immbody(:)
  type(sboun),allocatable,target :: immbond(:)
  integer,allocatable :: imb_node_have(:),imb_node_need(:),        &
                         num_icell_rank(:),num_ighost_rank(:)
  !+---------------------+---------------------------------------------+
  !|           num_solid | number of solid bodies                      |
  !|             immbody | the immersed body.                          |
  !|             immbond | the boundary nodes of immersed body.        |
  !+---------------------+---------------------------------------------+
  logical :: lreport
  !+---------------------+---------------------------------------------+
  !|             lreport | to control report of subroutines            |
  !+---------------------+---------------------------------------------+
  !
#ifdef COMB
  type(phase_t) :: mixture
  logical :: lcomb
  character(len=255) :: chemfile
  character(len=3) :: odetype
#endif 
  !
  parameter(hm=5)
  !
end module commvar
!+---------------------------------------------------------------------+
!| The end of the module commvar.                                      |
!+---------------------------------------------------------------------+