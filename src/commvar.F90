!+---------------------------------------------------------------------+
!| This module is to define common variables.                          |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module commvar
  !
  implicit none
  !
  integer :: ia,ja,ka,im,jm,km,is,ie,js,je,ks,ke
  integer :: hm
  integer :: numq,num_species,ndims
  character(len=4) :: conschm,difschm
  character(len=64) :: gridfile 
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
  !|            gridfile | the gridfile.                               |
  !+---------------------+---------------------------------------------+
  logical :: lihomo,ljhomo,lkhomo
  logical :: nondimen,diffterm,lfilter,lreadgrid
  !+---------------------+---------------------------------------------+
  !| lihomo,ljhomo,lkhomo| to define homogeneous direction.            |
  !|             nondimen| if the equation is non-dimensional          |
  !|             diffterm| if the diffusion terms is solved.           |
  !|                     | .f. means Euler equations are solved.       |
  !|             lfilter | to activate filer flag                      |
  !+---------------------+---------------------------------------------+
  !
  integer :: npdci,npdcj,npdck
  !+---------------------+---------------------------------------------+
  !|               npdci | to control scheme at boundary.              |
  !+---------------------+---------------------------------------------+
  real(8) :: xmax,xmin,ymax,ymin,zmax,zmin,voldom
  real(8) :: alfa_filter
  !+---------------------+---------------------------------------------+
  !|                *mix | min coordinates                             |
  !|                *max | max coordinates                             |
  !|         alfa_filter | the parameter to control width of filter.   |
  !|              voldom | total volume of the domain.                 |
  !+---------------------+---------------------------------------------+
  integer :: nstep,maxstep,nwrite,nlstep,filenumb
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

  !
  parameter(hm=5)
  !
end module commvar
!+---------------------------------------------------------------------+
!| The end of the module commvar.                                      |
!+---------------------------------------------------------------------+