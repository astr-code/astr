!+---------------------------------------------------------------------+
!| This module contains subroutines/functions for thermo-chemistry     |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 13-Aug-2020  | Created by Z.X. Chen @ Cambridge                     |
!+---------------------------------------------------------------------+
module thermchem
  !
  use commvar, only: num_species
#ifdef COMB
  use cantera
#endif 
  !
  implicit none
  !
  Interface rgcmix
  module procedure rgcmix_scar
  module procedure rgcmix_1d
  module procedure rgcmix_3d
  end Interface rgcmix
  !
#ifdef COMB
  !
  integer :: ncstep,nbody=0,ngibb=0,nlind=0,ntroe=0,nsrif=0
  real(8),parameter :: rguniv=8.3142d3
  real(8),parameter :: alamdc=2.58d-5,rlamda=7.0d-1,tlamda=2.98d2
  real(8) :: prefgb,alamda
  !
  character(len=10),allocatable :: spcsym(:),bdysym(:)
  character(len=256) :: chemxmlfile
  integer,allocatable :: &
    ntint(:),ncofcp(:,:),nsslen(:),nsspec(:,:),nrslen(:),nrspec(:,:),npslen(:) &
    ,npspec(:,:),nrclen(:),nrcpec(:,:),npclen(:),npcpec(:,:),mblist(:)         &
    ,mglist(:),mllist(:),mtlist(:),mslist(:),ncpoly(:,:),ncpom1(:,:)           &
    ,ncenth(:,:),ncenpy(:,:)
  real(8),allocatable :: &
    wmolar(:),clewis(:),tintlo(:,:),tinthi(:,:),amolcp(:,:,:),Arrhenius(:,:)   &
    ,crspec(:,:),cpspec(:,:),diffmu(:,:),effy3b(:,:),rclind(:,:),rctroe(:,:)   &
    ,rcsrif(:,:),diffmw(:,:),ovwmol(:),rgspec(:),amascp(:,:,:),amascv(:,:,:)   &
    ,amasct(:,:,:),amasch(:,:,:),amasce(:,:,:),amascs(:,:,:),amolgb(:,:,:)     &
    ,olewis(:),wirate(:)
    !
  type(phase_t) :: mixture
#endif
  !
  character(len=5) :: tranmod='mixav'
  logical :: ctrflag=.true.
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine reads the chemistry data from cantera format file |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-Feb-2021  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  subroutine chemread(cheminfile)
    !
    ! arguments
    character(len=*),intent(in) :: cheminfile
    !
#ifdef COMB
    ! local data
    integer :: &
      js,is,it,icp,jr,ir,j3b,i3b,ilind,jlind,itroe,jtroe,isrif,jsrif,  &
      ireac,iprod,icol,ispa,jreac,jprod
    real(8), allocatable :: &
      rcoeff(:,:),pcoeff(:,:),Arrhenius0(:,:),troeparams(:,:), &
      sriparams(:,:)
    real(8) :: dum
    character(len=10) :: spcnm,eunit,reactype,reverse,fallofftype,phase_id
    character(len=100) :: stringline
    !
    phase_id='gas'
    !---CANTERA---
    mixture=importPhase(cheminfile,trim(phase_id))
    if(speciesIndex(mixture,'O')<0 .and. &
    speciesIndex(mixture,'o')<0) &
      stop '!! Species "O" not exist, check phase_id !!'
    ncstep=nReactions(mixture)
    num_species=nSpecies(mixture)
    prefgb=refPressure(mixture)
    call setPressure(mixture,prefgb)
    chemxmlfile=cheminfile
    !
    !===============ALLOCATE=============================
    allocate(spcsym(num_species),wmolar(num_species),                  &
             clewis(num_species),ntint(num_species),                   &
             ovwmol(num_species),rgspec(num_species))
    !
    !---CANTERA---
    call getMolecularWeights(mixture,wmolar)
    !
    do js=1,num_species
      call getSpeciesName(mixture,js,spcsym(js))
    enddo
    !
    clewis(:)=1.d0
    !
    ! open(unit=11,file='datin/thermchem.txt',status='old',form='formatted')
    ! call skipline(11,3)
    ! do js=1,num_species
    !   read(11,*) spcsym(js),clewis(js)
    !   call getSpeciesName(mixture,js,spcnm)
    !   if(trim(spcnm)/=trim(spcsym(js))) &
    !     stop '!! Species name inconsistent ctr vs thermtxt !!'
    !   ! print*,js,spcsym(js),clewis(js),wmolar(js)
    ! enddo
    ! !
    ! !By default two-level NASA Polynomials
    ! ntint(:)=2
    ! species: do js=1,num_species
    !   if(.not.allocated(tintlo)) &
    !     allocate(tintlo(ntint(js),num_species),                     &
    !              tinthi(ntint(js),num_species),                     &
    !              ncofcp(ntint(js),num_species))
    !   tmp_levels: do it=1,ntint(js)
    !     read(11,*) spcnm,dum,tintlo(it,js), &
    !       tinthi(it,js),ncofcp(it,js)
    !     ! print*,spcnm,tintlo(it,js),tinthi(it,js),ncofcp(it,js)
    !     if(trim(spcnm)/=trim(spcsym(js))) &
    !       stop '!! Species name inconsistent at NASA Polynomials !!'
    !     !===============ALLOCATE=============================
    !     if(.not.allocated(amolcp)) allocate( &
    !       amolcp(ncofcp(it,js),ntint(js),num_species))
    !     !
    !     read(11,*)(amolcp(icp,it,js),icp=1,ncofcp(it,js))
    !     ! write(*,*)(amolcp(icp,it,js),icp=1,ncofcp(it,js))
    !     !
    !   enddo tmp_levels
    !       !
    ! enddo species
    ! !
    ! !===============ALLOCATE=============================
    ! allocate(Arrhenius(3,ncstep),nsslen(ncstep),nrslen(ncstep), &
    !          npslen(ncstep),nrclen(ncstep),npclen(ncstep), &
    !          rcoeff(ncstep,num_species),pcoeff(ncstep,num_species), &
    !          Arrhenius0(3,ncstep),troeparams(4,ncstep),sriparams(5,ncstep))
    ! !===============ALLOCATE=============================
    ! allocate(bdysym(ncstep),mblist(ncstep),mglist(ncstep),             &
    !          effy3b(num_species,ncstep),mllist(ncstep),                &
    !          mtlist(ncstep),mslist(ncstep))
    ! mblist(:)=0
    ! mglist(:)=0
    ! mllist(:)=0
    ! mtlist(:)=0
    ! mslist(:)=0
    ! effy3b(:,:)=1.d0
    ! troeparams(:,:)=0.d0
    ! icol=0
    ! !STEP RATE DATA
    ! call skipline(11,3)
    ! reactsteps: do jr=1,ncstep
    !   !
    !   read(11,*)ir,reactype,reverse
    !   ! write(*,'(I3,2A)') ir,' ',reactype
    !   !
    !   if(ir/=jr) stop '!! Reaction number inconsistent !!'
    !   !
    !   if(trim(reverse)=="yes") then
    !     ngibb=ngibb+1
    !     mglist(jr)=ir
    !   endif
    !   !
    !   call skipline(11,1)
    !   !
    !   read(11,*)(Arrhenius(icp,ir),icp=1,3),eunit
    !   !CONVERT TO J/KMOL
    !   if(eunit=="cal") Arrhenius(3,ir)=Arrhenius(3,ir)*4.186798d3
    !   ! write(*,'(I3,3(1PE12.4),2A)')ir,(Arrhenius(icp,ir),icp=1,3),' ',eunit
    !   !
    !   if(trim(reactype)/="elementary") then
    !     !
    !     read(11,'(A)')stringline
    !     ! write(*,'(A)')trim(stringline)
    !     !
    !     if(trim(stringline)/='') then
    !       !
    !       nbody=nbody+1
    !       mblist(ir)=nbody
    !       ! GET THREEBODY EFFICIENCIES
    !       do while(len(trim(stringline))>0)
    !         !
    !         icol=index(stringline,':')
    !         ispa=index(stringline,' ')
    !         ! print*,icol,ispa
    !         spcnm=stringline(1:icol-1)
    !         !---CANTERA---
    !         is=speciesIndex(mixture,trim(spcnm))
    !         !
    !         if(is==0) then
    !           stringline=''
    !           backspace(11)
    !         else
    !           ! print*,ir,is,spcnm,stringline(icol+1:ispa-1)
    !           read(stringline(icol+1:ispa-1),'(E5.3)')effy3b(is,nbody)
    !           stringline=stringline(ispa+1:)
    !           do while(index(stringline,' ')==1 &
    !                   .and.len(trim(stringline))>0)
    !             stringline=stringline(2:)
    !           enddo
    !         endif
    !         !
    !       enddo
    !       !
    !     endif
    !     !
    !     if(trim(reactype)=="falloff") then
    !       !
    !       read(11,*)fallofftype
    !       !
    !       read(11,*)(Arrhenius0(icp,ir),icp=1,3),eunit
    !       !CONVERT TO J/KMOL
    !       if(eunit=="cal") Arrhenius0(3,ir)=Arrhenius0(3,ir)*4.186798d3
    !       ! write(*,'(I3,3(1PE12.4),2A)')ir,(Arrhenius0(icp,ir),icp=1,3),' ',eunit
    !       !
    !       if(trim(fallofftype)=="Lindemann") then
    !         !
    !         nlind=nlind+1
    !         mllist(ir)=nlind
    !         !
    !       elseif(trim(fallofftype)=="Troe") then
    !         !
    !         ntroe=ntroe+1
    !         mtlist(ir)=ntroe

    !         read(11,'(A)')stringline
    !         stringline=trim(stringline)
    !         icp=0
    !         do while(len(trim(stringline))>0)
    !           icp=icp+1
    !           ispa=index(stringline,' ')
    !           read(stringline(1:ispa-1),*)troeparams(icp,ir)
    !           stringline=stringline(ispa+1:)
    !           do while(index(stringline,' ')==1 &
    !                   .and.len(trim(stringline))>0)
    !             stringline=stringline(2:)
    !           enddo
    !         enddo
    !         ! write(*,*)(troeparams(icp,jr),icp=1,4)
    !         !
    !       elseif(trim(fallofftype)=="SRI") then
    !         !
    !         nsrif=nsrif+1
    !         mslist(ir)=nsrif
    !         print*,' !!Warning - SRI reactions not validated!!'
    !         read(11,*)(sriparams(icp,ir),icp=1,5)
    !         !
    !       else
    !         !
    !         stop '!!Error - fallofftype not recognised!!'
    !         !
    !       endif
    !       !
    !     endif
    !     !
    !   endif
    !   !
    !   call skipline(11,1)

    ! enddo reactsteps
    ! !
    ! close(11) !thermchem.txt
    ! !
    ! !STEP SPECIES-LIST
    ! !===============ALLOCATE=============================
    ! ! MAX 10 SPECIES IN ONE STEP
    ! allocate(nsspec(10,ncstep),nrspec(10,ncstep),npspec(10,ncstep), &
    !          nrcpec(10,ncstep),crspec(10,ncstep), &
    !          npcpec(10,ncstep),cpspec(10,ncstep), &
    !          diffmu(10,ncstep),diffmw(10,ncstep))
    ! rcoeff(:,:)=0.d0
    ! pcoeff(:,:)=0.d0
    ! !
    ! do jr=1,ncstep
    !   is=0
    !   ireac=0
    !   iprod=0
    !   do js=1,num_species
    !     rcoeff(jr,js)=reactantStoichCoeff(mixture,js,jr)
    !     pcoeff(jr,js)=productStoichCoeff(mixture,js,jr)
    !     if(rcoeff(jr,js)>0.d0) then
    !       is=is+1
    !       nsspec(is,jr)=js
    !       jreac=nint(rcoeff(jr,js))
    !       do while(jreac>0)
    !         ireac=ireac+1
    !         nrspec(ireac,jr)=js
    !         jreac=jreac-1
    !       enddo
    !     endif
    !     if(pcoeff(jr,js)>0.d0) then
    !       is=is+1
    !       nsspec(is,jr)=js
    !       jprod=nint(pcoeff(jr,js))
    !       do while(jprod>0)
    !         iprod=iprod+1
    !         npspec(iprod,jr)=js
    !         jprod=jprod-1
    !       enddo
    !     endif
    !   enddo
    !   nsslen(jr)=is
    !   nrslen(jr)=ireac
    !   npslen(jr)=iprod
    !   ! if(jr==1)print*,'Step species-list:'
    !   ! print*,jr,nsslen(jr)
    !   ! do js=1,nsslen(jr)
    !   !   print*,jr,js,nsspec(js,jr)
    !   ! enddo
    ! enddo
    ! !
    ! !STEP REACTANT-LIST
    ! ! print*,'Step reactant-list:'
    ! ! do jr=1,ncstep
    ! !   print*,jr,nrslen(jr)
    ! !   do js=1,nrslen(jr)
    ! !     print*,jr,js,nrspec(js,jr)
    ! !   enddo
    ! ! enddo
    ! !
    ! !STEP PRODUCT-LIST
    ! ! print*,'Step product-list:'
    ! ! do jr=1,ncstep
    ! !   print*,jr,npslen(jr)
    ! !   do js=1,npslen(jr)
    ! !     print*,jr,js,npspec(js,jr)
    ! !   enddo
    ! ! enddo
    ! !
    ! !STEP REACTANT non-int COEFFICIENT-LIST (not used)
    ! ! print*,'Step reactant coefficient-list:'
    ! nrclen(:)=0
    ! nrcpec(:,:)=0
    ! crspec(:,:)=0.d0
    ! !
    ! !STEP PRODUCT non-int COEFFICIENT-LIST (not used)
    ! ! print*,'Step product coefficient-list:'
    ! npclen(:)=0
    ! npcpec(:,:)=0
    ! cpspec(:,:)=0.d0
    ! !
    ! !SPECIES DELTA-LIST
    ! ! print*,'Species delta-list:'
    ! do jr=1,ncstep
    !   do js=1,nsslen(jr)
    !     diffmu(js,jr)=pcoeff(jr,nsspec(js,jr)) &
    !                         -rcoeff(jr,nsspec(js,jr))
    !     ! print*,jr,js,diffmu(js,jr)
    !   enddo
    ! enddo
    ! !
    ! !THIRD-BODY LIST
    ! ! print*,'Third-body list:'
    ! ! print*,nbody
    ! do j3b=1,nbody
    !   write(bdysym(j3b),'(I4.4)')j3b
    !   bdysym(j3b)='M'//bdysym(j3b)
    !   ! print*,j3b,bdysym(j3b)
    ! enddo
    ! !
    ! !THIRD-BODY STEP-LIST
    ! if(nbody>0) then
    !   do jr=1,ncstep
    !     ! print*,jr,mblist(jr)
    !   enddo
    ! endif
    ! !
    ! !THIRD-BODY EFFICIENCIES
    ! ! do j3b=1,nbody
    ! !   do js=1,num_species
    ! !     print*,j3b,js,effy3b(js,j3b)
    ! !   enddo
    ! ! enddo
    ! !
    ! !GIBBS STEP-LIST
    ! ! write(*,'(I5)')ngibb
    ! ! if(ngibb>0) then
    ! !   do jr=1,ncstep
    ! !     print*,jr,mglist(jr)
    ! !   enddo
    ! ! endif
    ! !
    ! !LINDEMANN STEPS
    ! ! write(*,'(I5)')nlind
    ! allocate(rclind(4,nlind))
    ! if(nlind>0) then
    !   do jr=1,ncstep
    !     ! print*,jr,mllist(jr)
    !     ilind=mllist(jr)
    !     if(ilind>0) then
    !       do icp=1,3
    !         rclind(icp,ilind)=Arrhenius0(icp,jr)
    !       enddo
    !       rclind(4,ilind)=1.d0
    !       ! if(ilind>0)print*,ilind,(rclind(icp,ilind),icp=1,4)
    !     endif
    !   enddo
    ! endif
    ! !
    ! !
    ! !TROE STEPS
    ! ! write(*,'(I5)')ntroe
    ! allocate(rctroe(12,ntroe))
    ! if(ntroe>0) then
    !   do jr=1,ncstep
    !     ! print*,jr,mtlist(jr)
    !     itroe=mtlist(jr)
    !     if(itroe>0) then
    !       do icp=1,3
    !         rctroe(icp,itroe)=Arrhenius0(icp,jr)
    !       enddo
    !       rctroe(4,itroe)=troeparams(1,jr) !A
    !       rctroe(5,itroe)=troeparams(3,jr) !T1
    !       rctroe(6,itroe)=troeparams(4,jr) !T2
    !       rctroe(7,itroe)=troeparams(2,jr) !T3
    !       rctroe(8,itroe)=-0.4d0
    !       rctroe(9,itroe)=-0.67d0
    !       rctroe(10,itroe)=0.75d0
    !       rctroe(11,itroe)=-1.27d0
    !       rctroe(12,itroe)=0.14d0
    !       ! if(itroe>0)print*,itroe,(rctroe(icp,itroe),icp=1,12)
    !     endif
    !   enddo
    ! endif
    ! !
    ! !SRI STEPS
    ! ! write(*,'(I5)')nsrif
    ! allocate(rcsrif(8,nsrif))
    ! if(nsrif>0) then
    !   do jr=1,ncstep
    !     ! print*,jr,mslist(jr)
    !     isrif=mslist(jr)
    !     do icp=1,3
    !       rcsrif(icp,isrif)=Arrhenius0(icp,jr)
    !     enddo
    !     rcsrif(4,itroe)=sriparams(1,jr) !a
    !     rcsrif(5,itroe)=sriparams(2,jr) !b
    !     rcsrif(6,itroe)=sriparams(3,jr) !c
    !     rcsrif(7,itroe)=sriparams(4,jr) !d
    !     rcsrif(8,itroe)=sriparams(5,jr) !e
    !     ! if(isrif>0)print*,isrif,(rcsrif(icp,isrif),icp=1,8)
    !   enddo
    ! endif
    ! !
    ! !=============================================================
    ! !===============EVALUATE DERIVED QUANTITIES===================
    ! !=============================================================
    ! !CONVERT RATE PARAMETERS
    ! do ir=1,ncstep
    !   Arrhenius(1,ir)=log(Arrhenius(1,ir))
    !   Arrhenius(3,ir)=Arrhenius(3,ir)/rguniv
    ! enddo
    ! !
    ! !LINDEMANN STEP RATE DATA
    ! do ilind=1,nlind
    !   rclind(1,ilind)=log(rclind(1,ilind))
    !   rclind(3,ilind)=rclind(3,ilind)/rguniv
    ! enddo
    ! !
    ! !TROE FORM STEP RATE DATA
    ! do itroe=1,ntroe
    !   rctroe(1,itroe)=log(rctroe(1,itroe))
    !   rctroe(3,itroe)=rctroe(3,itroe)/rguniv
    !   rctroe(5,itroe)=-1.0d0/rctroe(5,itroe)!T1
    !   rctroe(7,itroe)=-1.0d0/rctroe(7,itroe)!T3
    ! enddo
    ! !
    ! !SRI FORM STEP RATE DATA
    ! do isrif=1,nsrif
    !   rcsrif(1,isrif)=log(rcsrif(1,isrif))
    !   rcsrif(3,isrif)=rcsrif(3,isrif)/rguniv
    !   rcsrif(5,isrif)=-rcsrif(5,isrif) !b
    !   rcsrif(6,isrif)=-1.0d0/rcsrif(6,isrif) !c
    ! enddo
    ! !
    ! !STOICHIOMETRIC COEFFICIENTS TIMES MOLAR MASS
    ! do ir=1,ncstep
    !   do is=1,nsslen(ir)
    !     js=nsspec(is,ir)
    !     diffmw(is,ir)=diffmu(is,ir)*wmolar(js)
    !   enddo
    ! enddo
    !
    !RECIPROCAL OF MOLAR MASS
    ovwmol(:)=1.0d0/wmolar(:)
    !SPECIFIC mixture CONSTANT
    rgspec(:)=rguniv*ovwmol(:)
    !
#endif
    !
  end subroutine chemread
  !+-------------------------------------------------------------------+
  !| The end of the subroutine chemread_ctr                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine prints the chemistry data for display.            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Aug-2020  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  subroutine chemrep(filename)
    ! arguments
    character(len=*),intent(in) :: filename
    !
#ifdef COMB
    ! local data
    integer :: js,icp,jcp,kcp,jr,istr1,istr2,ks
    ! real(8) :: ttemp(5),fornow,ttold(5)
    ! character*132 char132
    ! character*10 char10
    ! character*5 char5
    ! character*4 char4
    ! character*1 char1
    ! !
    open(16,file=filename)
    !
    ! SPECIES LIST, 43 char length per line
    write(16,'(A)')' +------------ Chemical Data -------------+'
    write(16,'(A,I5,A18)')' Number of species:',num_species,''
    write(16,'(A7,3X,A7,3X,A8,4X,A9,A2)') &
    ' Index','Species','Mol.Mass','Lewis No.',''
    do js=1,num_species
      write(16,'(A2,I5,3X,A7,3X,1PE9.3,3X,1PE9.3,A2)') &
      ' ',js,spcsym(js),wmolar(js),clewis(js),''
    enddo
    !
    ! write(16,*)
    ! !THERMODYNAMIC DATA
    ! write(16,'(A,A14)')' Species thermodynamic data:',''
    ! write(16,'(A,1PE12.4,A10)')' Reference pressure:',prefgb,''
    ! write(16,*)'Spec.  No of T intervals', &
    !   '  Interval   T low       T high       No of coeffs'
    ! do js=1,num_species
    !   icp=1
    !   write(16,'(I5,6X,I5,9X,I5,8X,2(1PE12.4),I8)') &
    !     js,ntint(js),icp,tintlo(icp,js),  &
    !     tinthi(icp,js),ncofcp(icp,js)
    !   do icp=2,ntint(js)
    !     write(16,'(25X,I5,8X,2(1PE12.4),I8)') &
    !       icp,tintlo(icp,js),tinthi(icp,js),ncofcp(icp,js)
    !   enddo
    ! enddo
    ! write(16,*)'Cp coeffs by mass'
    ! write(16,*)' Spec.  T int.  Coeff no.  Coeff.'
    ! do js=1,num_species
    !   do icp=1,ntint(js)
    !     jcp=1
    !     write(16,'(I5,2X,I5,5X,I5,4X,1PE15.7)') &
    !       js,icp,jcp,amascp(jcp,icp,js)
    !     do jcp=2,ncofcp(icp,js)
    !       write(16,'(17X,I5,4X,1PE15.7)')jcp,amascp(jcp,icp,js)
    !     enddo
    !   enddo
    ! enddo
    ! !
    ! write(16,*)
    ! write(16,*)'Mass-specific Cp, Enthalpy, Entropy;  Molar Gibbs fn.:'
    ! write(16,*)' Spec.  T int.  Temp.       Cp',  &
    !                   '          Enthalpy    Entropy     Molar Gibbs'
    ! species: do js=1,num_species
    !   tmp_levels: do icp=1,ntint(js)
    !     tmp_limits: do kcp=1,2
    !       !Temperature
    !       if(kcp==1)ttemp(1)=tintlo(icp,js)
    !       if(kcp==2)ttemp(1)=tinthi(icp,js)
    !       !MASS-SPECifIC CP
    !       fornow=amascp(ncpoly(icp,js),icp,js)
    !       do jcp=ncpom1(icp,js),1,-1
    !         fornow=fornow*ttemp(1)+amascp(jcp,icp,js)
    !       enddo
    !       ttemp(2)=fornow
    !       !MASS-SPECifIC ENTHALPY
    !       fornow=amasch(ncpoly(icp,js),icp,js)
    !       do jcp=ncpom1(icp,js),1,-1
    !         fornow=fornow*ttemp(1)+amasch(jcp,icp,js)
    !       enddo
    !       fornow=amasch(ncenth(icp,js),icp,js)+fornow*ttemp(1)
    !       ttemp(3)=fornow
    !       !MASS-SPECifIC ENTROPY
    !       fornow=amascs(ncpoly(icp,js),icp,js)
    !       do jcp=ncpom1(icp,js),2,-1
    !         fornow=fornow*ttemp(1)+amascs(jcp,icp,js)
    !       enddo
    !       fornow=amascs(ncenpy(icp,js),icp,js)+fornow*ttemp(1)  &
    !              +amascs(1,icp,js)*log(ttemp(1))
    !       ttemp(4)=fornow

    !       !MOLAR GIBBS FUNCTION: GIBBS/(R^0 T) WITH PRESSURE TERM
    !       fornow=amolgb(ncpoly(icp,js),icp,js)
    !       do jcp=ncpom1(icp,js),1,-1
    !         fornow=amolgb(jcp,icp,js)+fornow*ttemp(1)
    !       enddo
    !       fornow=amolgb(ncenth(icp,js),icp,js)/ttemp(1) &
    !               - amolgb(ncenpy(icp,js),icp,js)*log(ttemp(1)) &
    !               - fornow
    !       ttemp(5)=fornow

    !       if(kcp==1)then
    !         !
    !         if(icp==1)then
    !           write(16,'(I5,2X,I5,X,"l",X,5(1PE13.5))') &
    !                js,icp,(ttemp(jcp),jcp=1,5)
    !         else
    !           write(16,'(7X,I5,X,"l",X,5(1PE13.5))')  &
    !                      icp,(ttemp(jcp),jcp=1,5)
    !           do jcp=1,5
    !             if(abs(ttold(jcp)-ttemp(jcp))>abs(1.0d-3*ttemp(jcp)))then
    !               write(16,*)'Warning: INDATA: Mismatched thermo data'
    !               write(16,'(I7,1PE12.4)')jcp,ttemp(jcp)
    !             endif
    !           enddo
    !         endif
    !         !
    !       else
    !         write(16,'(7X,I5,X,"h",X,5(1PE13.5))')icp,(ttemp(jcp),jcp=1,5)
    !       endif
    !       !
    !       if(kcp==2)then
    !         do jcp=1,5
    !           ttold(jcp)=ttemp(jcp)
    !         enddo
    !       endif
    !       !
    !     enddo tmp_limits
    !   enddo tmp_levels
    ! enddo species
    ! !
    ! ! REACTION DATA
    ! write(16,*)'Reaction mechanism:'
    ! write(16,*)'  Number of steps:'
    ! write(16,'(I5)')ncstep
    ! do jr=1,ncstep
    !   write(char5,'(I5)')jr
    !   istr1=1
    !   istr2=len(char5)+3
    !   char132(istr1:istr2)=char5//'  '
    !   do js=1,nrslen(jr)
    !     char10=spcsym(nrspec(js,jr))
    !     istr1=istr2+1
    !     istr2=istr1+len(char10)
    !     char132(istr1:istr2)=char10
    !     istr1=istr2+1
    !     istr2=istr1
    !     char132(istr1:istr2)=' + '
    !   enddo
    !   if(mblist(jr)>0)then
    !     char10=bdysym(mblist(jr))
    !     istr1=istr2+1
    !     istr2=istr1+len(char10)
    !     char132(istr1:istr2)=char10
    !   else
    !     istr2=istr2 - 1
    !   endif
    !   char4=' => '
    !   if(mglist(jr)>0)char4=' == '
    !   istr1=istr2+1
    !   istr2=istr1+len(char4)
    !   char132(istr1:istr2)=char4
    !   !
    !   do js=1,nsslen(jr)
    !     icp=0
    !     do ks=1,nrslen(jr)
    !       if(nrspec(ks,jr)==nsspec(js,jr))then
    !         icp=icp+1
    !       endif
    !     enddo
    !     icp=nint(diffmu(js,jr))+icp
    !     if(icp>0)then
    !       if(icp>1)then
    !         write(char1,'(I1)')icp
    !         istr1=istr2+1
    !         istr2=istr1+len(char1)
    !         char132(istr1:istr2)=char1
    !       endif
    !       char10=spcsym(nsspec(js,jr))
    !       istr1=istr2+1
    !       istr2=istr1+len(char10)
    !       char132(istr1:istr2)=char10
    !       istr1=istr2+1
    !       istr2=istr1
    !       char132(istr1:istr2)=' + '
    !     endif
    !   enddo
    !   !
    !   if(mblist(js)>0)then
    !     char10=bdysym(mblist(js))
    !     istr1=istr2+1
    !     istr2=istr1+len(char10)
    !     char132(istr1:istr2)=char10
    !   else
    !     istr2=istr2 - 1
    !   endif
    !   !
    !   write(16,'(A)')char132(1:istr2)
    !   !
    !   write(16,'(A)')'Step species-list:'
    !   write(16,'(2I5)')jr,nsslen(jr)
    !   do js=1,nsslen(jr)
    !     write(16,'(3I5)')jr,js,nsspec(js,jr)
    !   enddo
    !   !
    !   write(16,'(A)')'Step reactant-list:'
    !   write(16,'(2I5)')jr,nrslen(jr)
    !   do js=1,nrslen(jr)
    !     write(16,'(3I5)')jr,js,nrspec(js,jr)
    !   enddo
    !   !
    !   write(16,'(A)')'Step product-list:'
    !   write(16,'(2I5)')jr,npslen(jr)
    !   do js=1,npslen(jr)
    !     write(16,'(3I5)')jr,js,npspec(js,jr)
    !   enddo
    !   !
    !   write(16,'(A)')'Step delta-list:'
    !   do js=1,nsslen(jr)
    !     write(16,'(2I5,1X,1PE12.5)')jr,js,diffmu(js,jr)
    !   enddo
    !   !
    ! enddo
    ! !
    ! write(16,*)
    ! !
    ! write(16,*)'Reaction parameters A, n, E:'
    ! do jr=1,ncstep
    !   write(16,'(I5,3(2X,1PE12.4))') &
    !        jr,(Arrhenius(icp,jr),icp=1,3)
    ! enddo
    ! !
    ! if(nlind>0)then
    ! write(16,*)
    ! write(16,*)'Lindemann parameters:'
    ! do jr=1,ncstep
    !   if(mllist(jr)/=0)then
    !     write(16,'(I5,3(2X,1PE12.4))') &
    !          jr,(rclind(icp,mllist(jr)),icp=1,4)
    !   endif
    ! enddo
    ! endif
    ! !
    ! if(ntroe>0)then
    ! write(16,*)
    ! write(16,*)'Troe parameters:'
    ! do jr=1,ncstep
    !   if(mtlist(jr)/=0)then
    !     write(16,'(I5,6(2X,1PE12.4))')jr,(rctroe(icp,mtlist(jr)),icp=1,6)
    !     write(16,'(I5,6(2X,1PE12.4))')jr,(rctroe(icp,mtlist(jr)),icp=7,12)
    !   endif
    ! enddo
    ! endif
    ! !
    ! !
    ! if(nsrif>0)then
    ! write(16,*)
    ! write(16,*)'SRI parameters:'
    ! do jr=1,ncstep
    !   if(mslist(jr)/=0)then
    !     write(16,'(I5,8(2X,1PE12.4))')jr,(rcsrif(icp,mslist(jr)),icp=1,8)
    !   endif
    ! enddo
    ! endif
    ! !
    ! write(16,*)
    ! write(16,*)'Third body efficiencies:'
    ! do ks=1,nbody
    !   write(16,'(I5,3X,A10)')ks,bdysym(ks)
    !   do js=1,num_species
    !     write(16,'(I5,3X,A10,1PE12.4)') &
    !          js,spcsym(js),effy3b(js,ks)
    !   enddo
    ! enddo
    ! !
    ! write(16,'(A)')
    ! write(16,'(A)')' End of chemical data'
    ! write(16,'(A)')' +----------------------------------------+'
    ! !
    close(16)
    !
    print*,' << ',filename
    !
#endif
    !
  end subroutine chemrep
  !+-------------------------------------------------------------------+
  !| The end of the function chemrep.                                  |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to compute thermodynamic quantities       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-Aug-2020: Created by Z.X. Chen @ Cambridge                     |
  !+-------------------------------------------------------------------+
  subroutine thermdyn
    !
#ifdef COMB
    use commvar,  only : num_species
    use constdef, only : pi
    !
    ! local data
    integer :: is,it,itnm1,jc,icp
    !
    real(8) :: &
      tbreak,twidth,ovtwid,ovtwrp,ttemp1,cpmol1,ttemp2,cpmol2,deltcp   &
      ,giblet
    !
    !===============ALLOCATE=============================
    allocate( &
      !  ncpoly(ntint(1),num_species),ncpom1(ntint(1),num_species)       &
      ! ,ncenth(ntint(1),num_species),ncenpy(ntint(1),num_species)       &
      ! ,amascp(size(amolcp,1),size(amolcp,2),size(amolcp,3))            &
      ! ,amascv(size(amolcp,1),size(amolcp,2),size(amolcp,3))            &
      ! ,amasct(size(amolcp,1),size(amolcp,2),size(amolcp,3))            &
      ! ,amasch(size(amolcp,1),size(amolcp,2),size(amolcp,3))            &
      ! ,amasce(size(amolcp,1),size(amolcp,2),size(amolcp,3))            &
      ! ,amascs(size(amolcp,1),size(amolcp,2),size(amolcp,3))            &
      ! ,amolgb(size(amolcp,1),size(amolcp,2),size(amolcp,3)),           &
      olewis(num_species),wirate(num_species)   &
      )

    ! !NUMBER OF CP POLYNOMIAL COEFFICIENTS AND DITTO MINUS ONE
    ! !INDEX NUMBERS OF ENTHALPY AND ENTROPY COEFFICIENTS
    ! do is=1,num_species
    !   do it=1,ntint(is)
    !     ncpoly(it,is)=ncofcp(it,is)-2
    !     ncpom1(it,is)=ncpoly(it,is)-1
    !     ncenth(it,is)=ncofcp(it,is)-1
    !     ncenpy(it,is)=ncofcp(it,is)
    !   enddo
    ! enddo
    ! !
    ! !BLENDER DATA
    ! tbreak=tinthi(1,1)
    ! twidth=2.0d1
    ! ovtwid=1.0d0/twidth
    ! ovtwrp=1.0d0/(sqrt(pi)*twidth)
    ! !
    ! !CHECK CP COEFFICIENTS FOR CONTINUITY AT INTERVAL BREAKPOINTS
    ! do is=1,num_species
    !   do it=2,ntint(is)
    !     itnm1=it-1
    !     ttemp1=tinthi(itnm1,is)
    !     cpmol1=amolcp(ncpoly(itnm1,is),itnm1,is)
    !     do jc=ncpom1(itnm1,is),1,-1
    !       cpmol1=cpmol1*ttemp1+amolcp(jc,itnm1,is)
    !     enddo
    !     ttemp2=tintlo(it,is)
    !     cpmol2=amolcp(ncpoly(it,is),it,is)
    !     do jc=ncpom1(it,is),1,-1
    !       cpmol2=cpmol2*ttemp2+amolcp(jc,it,is)
    !     enddo
    !     deltcp=cpmol2-cpmol1
    !     ! write(*,'(2I3,5E12.4)')is,it,ttemp1,ttemp2,cpmol1,       &
    !     !                        cpmol2,deltcp
    !     amolcp(1,it,is)=amolcp(1,it,is)-deltcp
    !   enddo !it
    ! enddo !is
    ! !
    ! !RECHECK CP COEFFICIENTS FOR CONTINUITY AT INTERVAL BREAKPOINTS
    ! do is=1,num_species
    !   do it=2,ntint(is)
    !     itnm1=it-1
    !     ttemp1=tinthi(itnm1,is)
    !     cpmol1=amolcp(ncpoly(itnm1,is),itnm1,is)
    !     do jc=ncpom1(itnm1,is),1,-1
    !       cpmol1=cpmol1*ttemp1+amolcp(jc,itnm1,is)
    !     enddo
    !     ttemp2=tintlo(it,is)
    !     cpmol2=amolcp(ncpoly(it,is),it,is)
    !     do jc=ncpom1(it,is),1,-1
    !       cpmol2=cpmol2*ttemp2+amolcp(jc,it,is)
    !     enddo
    !     deltcp=cpmol2-cpmol1
    !     ! write(*,'(2I3,5E12.4)')is,it,ttemp1,ttemp2,cpmol1,       &
    !     !                        cpmol2,deltcp
    !   enddo !it
    ! enddo !is
    ! !
    ! !SPECIFIC HEAT CAPACITY CP PER UNIT MASS
    ! do is=1,num_species
    !   do it=1,ntint(is)
    !     do icp=1,ncofcp(it,is)
    !       amascp(icp,it,is)  &
    !       = amolcp(icp,it,is)*rgspec(is)
    !       ! write(*,'(3I3,2E12.4)')is,it,icp,                      &
    !       !          amolcp(icp,it,is),amascp(icp,it,is)
    !     enddo
    !   enddo
    ! enddo
    ! !
    ! !CHECK CP COEFFICIENTS FOR CONTINUITY AT INTERVAL BREAKPOINTS
    ! do is=1,num_species
    !   do it=2,ntint(is)
    !   itnm1=it-1
    !   ttemp1=tinthi(itnm1,is)
    !   cpmol1=amascp(ncpoly(itnm1,is),itnm1,is)
    !   do jc=ncpom1(itnm1,is),1,-1
    !     cpmol1=cpmol1*ttemp1+amascp(jc,itnm1,is)
    !   enddo
    !   ttemp2=tintlo(it,is)
    !   cpmol2=amascp(ncpoly(it,is),it,is)
    !   do jc=ncpom1(it,is),1,-1
    !     cpmol2=cpmol2*ttemp2+amascp(jc,it,is)
    !   enddo
    !   deltcp=cpmol2-cpmol1
    !   ! write(*,'(2I3,5E12.4)')is,it,ttemp1,ttemp2,cpmol1,      &
    !   !                        cpmol2,deltcp
    !   amascp(1,it,is)=amascp(1,it,is)-deltcp
    !   enddo !it
    ! enddo !is
    ! !
    ! !RECHECK CP COEFFICIENTS FOR CONTINUITY AT INTERVAL BREAKPOINTS
    ! do is=1,num_species
    !   do it=2,ntint(is)
    !     itnm1=it-1
    !     ttemp1=tinthi(itnm1,is)
    !     cpmol1=amascp(ncpoly(itnm1,is),itnm1,is)
    !     do jc=ncpom1(itnm1,is),1,-1
    !       cpmol1=cpmol1*ttemp1+amascp(jc,itnm1,is)
    !     enddo
    !     ttemp2=tintlo(it,is)
    !     cpmol2=amascp(ncpoly(it,is),it,is)
    !     do jc=ncpom1(it,is),1,-1
    !       cpmol2=cpmol2*ttemp2+amascp(jc,it,is)
    !     enddo
    !     deltcp=cpmol2-cpmol1
    !     ! write(*,'(2I3,5E12.4)')is,it,ttemp1,ttemp2,cpmol1,      &
    !     !                        cpmol2,deltcp
    !   enddo !it
    ! enddo !is
    ! !
    ! !SPECIFIC HEAT CAPACITY CV PER UNIT MASS
    ! do is=1,num_species
    !   do it=1,ntint(is)
    !     do icp=1,ncofcp(it,is)
    !       amascv(icp,it,is)=amascp(icp,it,is)
    !     enddo
    !     amascv(1,it,is)=amascp(1,it,is)-rgspec(is)
    !   enddo
    ! enddo
    ! !CHECK CV COEFFICIENTS FOR CONTINUITY AT INTERVAL BREAKPOINTS
    ! do is=1,num_species
    !   do it=2,ntint(is)
    !     itnm1=it-1
    !     ttemp1=tinthi(itnm1,is)
    !     cpmol1=amascv(ncpoly(itnm1,is),itnm1,is)
    !     do jc=ncpom1(itnm1,is),1,-1
    !       cpmol1=cpmol1*ttemp1+amascv(jc,itnm1,is)
    !     enddo
    !     ttemp2=tintlo(it,is)
    !     cpmol2=amascv(ncpoly(it,is),it,is)
    !     do jc=ncpom1(it,is),1,-1
    !       cpmol2=cpmol2*ttemp2+amascv(jc,it,is)
    !     enddo
    !     deltcp=cpmol2-cpmol1
    !     ! write(*,'(2I3,5E12.4)')is,it,ttemp1,ttemp2,cpmol1,     &
    !     !                        cpmol2,deltcp
    !   enddo !it
    ! enddo !is
    ! !
    ! !COEFFICIENTS FOR TEMPERATURE
    ! do is=1,num_species
    !   do it=1,ntint(is)
    !     amasct(1,it,is)=amascp(1,it,is)-rgspec(is)
    !       do icp=2,ncpoly(it,is)
    !         amasct(icp,it,is)=amascp(icp,it,is)/real(icp)
    !       enddo
    !     amasct(ncenth(it,is),it,is)=0.d0
    !     amasct(ncenpy(it,is),it,is)=0.d0
    !   enddo
    ! enddo
    ! !
    ! !COEFFICIENTS FOR ENTHALPY PER UNIT MASS
    ! do is=1,num_species
    !   do it=1,ntint(is)
    !     amasch(1,it,is)=amascp(1,it,is)
    !       do icp=2,ncpoly(it,is)
    !         amasch(icp,it,is)=amascp(icp,it,is)/real(icp)
    !       enddo
    !     amasch(ncenth(it,is),it,is)  &
    !     = amascp(ncenth(it,is),it,is)
    !     amasch(ncenpy(it,is),it,is)=0.d0
    !   enddo
    ! enddo
    ! !
    ! !CHECK ENTHALPY COEFFICIENTS FOR CONTINUITY AT INTERVAL BREAKPOINTS
    ! do is=1,num_species
    !   do it=2,ntint(is)
    !     itnm1=it-1
    !     ttemp1=tinthi(itnm1,is)
    !     cpmol1=amasch(ncpoly(itnm1,is),itnm1,is)
    !     do jc=ncpom1(itnm1,is),1,-1
    !       cpmol1=cpmol1*ttemp1+amasch(jc,itnm1,is)
    !     enddo
    !     cpmol1=cpmol1*ttemp1  &
    !            + amasch(ncenth(itnm1,is),itnm1,is)
    !     ttemp2=tintlo(it,is)
    !     cpmol2=amasch(ncpoly(it,is),it,is)
    !     do jc=ncpom1(it,is),1,-1
    !       cpmol2=cpmol2*ttemp2+amasch(jc,it,is)
    !     enddo
    !     cpmol2=cpmol2*ttemp2  &
    !            + amasch(ncenth(it,is),it,is)
    !     deltcp=cpmol2-cpmol1
    !     ! write(*,'(2I3,5E12.4)')is,it,ttemp1,ttemp2,cpmol1,     &
    !     !                        cpmol2,deltcp
    !     amasch(ncenth(it,is),it,is)  &
    !     = amasch(ncenth(it,is),it,is)-deltcp
    !     amascp(ncenth(it,is),it,is)  &
    !     = amasch(ncenth(it,is),it,is)
    !   enddo !it
    ! enddo !is
    ! !
    ! !RECHECK ENTHALPY COEFFICIENTS FOR CONTINUITY AT INTERVAL BREAKPOINTS
    ! do is=1,num_species
    !   do it=2,ntint(is)
    !     itnm1=it-1
    !     ttemp1=tinthi(itnm1,is)
    !     cpmol1=amasch(ncpoly(itnm1,is),itnm1,is)
    !     do jc=ncpom1(itnm1,is),1,-1
    !       cpmol1=cpmol1*ttemp1+amasch(jc,itnm1,is)
    !     enddo
    !     cpmol1=cpmol1*ttemp1  &
    !            + amasch(ncenth(itnm1,is),itnm1,is)
    !     ttemp2=tintlo(it,is)
    !     cpmol2=amasch(ncpoly(it,is),it,is)
    !     do jc=ncpom1(it,is),1,-1
    !       cpmol2=cpmol2*ttemp2+amasch(jc,it,is)
    !     enddo
    !     cpmol2=cpmol2*ttemp2  &
    !            + amasch(ncenth(it,is),it,is)
    !     deltcp=cpmol2-cpmol1
    !     ! write(*,'(2I3,5E12.4)')is,it,ttemp1,ttemp2,cpmol1,     &
    !     !                        cpmol2,deltcp
    !   enddo !it
    ! enddo !is
    ! !
    ! !COEFFICIENTS FOR INTERNAL ENERGY PER UNIT MASS
    ! do is=1,num_species
    !   do it=1,ntint(is)
    !     amasce(1,it,is)=amasch(1,it,is)-rgspec(is)
    !       do icp=2,ncpoly(it,is)
    !         amasce(icp,it,is)=amasch(icp,it,is)
    !       enddo
    !     amasce(ncenth(it,is),it,is)  &
    !     = amasch(ncenth(it,is),it,is)
    !     amasce(ncenpy(it,is),it,is)=0.d0
    !   enddo
    ! enddo
    ! !
    ! !CHECK INTERNAL ENERGY COEFFICIENTS FOR CONTINUITY AT BREAKPOINTS
    ! do is=1,num_species
    !   do it=2,ntint(is)
    !     itnm1=it-1
    !     ttemp1=tinthi(itnm1,is)
    !     cpmol1=amasce(ncpoly(itnm1,is),itnm1,is)
    !     do jc=ncpom1(itnm1,is),1,-1
    !       cpmol1=cpmol1*ttemp1+amasce(jc,itnm1,is)
    !     enddo
    !     cpmol1=cpmol1*ttemp1  &
    !            + amasce(ncenth(itnm1,is),itnm1,is)
    !     ttemp2=tintlo(it,is)
    !     cpmol2=amasch(ncpoly(it,is),it,is)
    !     do jc=ncpom1(it,is),1,-1
    !       cpmol2=cpmol2*ttemp2+amasce(jc,it,is)
    !     enddo
    !     cpmol2=cpmol2*ttemp2  &
    !            + amasce(ncenth(it,is),it,is)
    !     deltcp=cpmol2-cpmol1
    !     ! write(*,'(2I3,5E12.4)')is,it,ttemp1,ttemp2,cpmol1,     &
    !     !                        cpmol2,deltcp
    !   enddo !it
    ! enddo !is
    ! !
    ! !COEFFICIENTS FOR ENTROPY PER UNIT MASS
    ! do is=1,num_species
    !   do it=1,ntint(is)
    !     amascs(1,it,is)=amascp(1,it,is)
    !       do icp=2,ncpoly(it,is)
    !         amascs(icp,it,is)=amascp(icp,it,is)          &
    !                                 /real(icp-1)
    !       enddo
    !     amascs(ncenth(it,is),it,is)=0.d0
    !     amascs(ncenpy(it,is),it,is)  &
    !     = amascp(ncenpy(it,is),it,is)
    !   enddo
    ! enddo
    ! !
    ! !CHECK ENTROPY COEFFICIENTS FOR CONTINUITY AT INTERVAL BREAKPOINTS
    ! do is=1,num_species
    !   do it=2,ntint(is)
    !     itnm1=it-1
    !     ttemp1=tinthi(itnm1,is)
    !     cpmol1=amascs(ncpoly(itnm1,is),itnm1,is)
    !     do jc=ncpom1(itnm1,is),2,-1
    !       cpmol1=cpmol1*ttemp1+amascs(jc,itnm1,is)
    !     enddo
    !     cpmol1=cpmol1*ttemp1  &
    !            + amascs(ncenpy(itnm1,is),itnm1,is)  &
    !            + amascs(1,itnm1,is)*log(ttemp1)
    !     ttemp2=tintlo(it,is)
    !     cpmol2=amascs(ncpoly(it,is),it,is)
    !     do jc=ncpom1(it,is),2,-1
    !       cpmol2=cpmol2*ttemp2+amascs(jc,it,is)
    !     enddo
    !     cpmol2=cpmol2*ttemp2  &
    !            + amascs(ncenpy(it,is),it,is)  &
    !            + amascs(1,it,is)*log(ttemp2)
    !     deltcp=cpmol2-cpmol1
    !     ! write(*,'(2I3,5E12.4)')is,it,ttemp1,ttemp2,cpmol1,     &
    !     !                        cpmol2,deltcp
    !     amascs(ncenpy(it,is),it,is)  &
    !     = amascs(ncenpy(it,is),it,is)-deltcp
    !   enddo !it
    ! enddo !is
    ! !
    ! !RECHECK ENTROPY COEFFICIENTS FOR CONTINUITY AT INTERVAL BREAKPOINTS
    ! do is=1,num_species
    !   do it=2,ntint(is)
    !     itnm1=it-1
    !     ttemp1=tinthi(itnm1,is)
    !     cpmol1=amascs(ncpoly(itnm1,is),itnm1,is)
    !     do jc=ncpom1(itnm1,is),2,-1
    !       cpmol1=cpmol1*ttemp1+amascs(jc,itnm1,is)
    !     enddo
    !     cpmol1=cpmol1*ttemp1  &
    !            + amascs(ncenpy(itnm1,is),itnm1,is)  &
    !            + amascs(1,itnm1,is)*log(ttemp1)
    !     ttemp2=tintlo(it,is)
    !     cpmol2=amascs(ncpoly(it,is),it,is)
    !     do jc=ncpom1(it,is),2,-1
    !       cpmol2=cpmol2*ttemp2+amascs(jc,it,is)
    !     enddo
    !     cpmol2=cpmol2*ttemp2  &
    !            + amascs(ncenpy(it,is),it,is)  &
    !            + amascs(1,it,is)*log(ttemp2)
    !     deltcp=cpmol2-cpmol1
    !     ! write(*,'(2I3,5E12.4)')is,it,ttemp1,ttemp2,cpmol1,     &
    !     !                        cpmol2,deltcp
    !   enddo !it
    ! enddo !is
    ! !
    ! !COEFFICIENTS FOR GIBBS FUNCTION PER MOLE
    ! !ACTUALLY GIBBS/(R^0 T) WITH PRESSURE TERM
    ! giblet=log(prefgb/rguniv)
    ! do is=1,num_species
    !   do it=1,ntint(is)
    !     amolgb(1,it,is)  &
    !     = amolcp(ncenpy(it,is),it,is)  &
    !       - amolcp(1,it,is) + giblet
    !     do icp=2,ncpoly(it,is)
    !       amolgb(icp,it,is)=amolcp(icp,it,is)            &
    !                               /real(icp*(icp-1))
    !     enddo
    !     amolgb(ncenth(it,is),it,is)  &
    !     = amolcp(ncenth(it,is),it,is)
    !     amolgb(ncenpy(it,is),it,is)  &
    !     = amolcp(1,it,is) - 1.0d0
    !   enddo !it
    ! enddo !is
    !
    !RECIPROCAL OF LEWIS NUMBER
    do is=1,num_species
      olewis(is)=1.0d0/clewis(is)
    enddo
    !
    !CONDUCTIVITY COEFFICIENT
    alamda=alamdc*exp(-rlamda*log(tlamda))
    !
#endif
    !
  end subroutine thermdyn
  !+-------------------------------------------------------------------+
  !| The end of the subroutine cheminit.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This funcion computes speed of sound.                             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-Aug-2020  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  subroutine aceval(tmp,spc,css)
    !
    ! arguments
    real(8),intent(in) :: tmp,spc(:)
    real(8),intent(out) :: css
    !
#ifdef COMB
    ! local data
    real(8) :: cpcmix,gamrgc
    !
    call setState_TPY(mixture,tmp,prefgb,spc(:))
    cpcmix=cp_mass(mixture)
    !
    gamrgc=rgcmix(spc)*cpcmix/(cpcmix-rgcmix(spc))
    !
    css=sqrt(tmp*gamrgc)
    !
#endif 
    !
  end subroutine aceval
  !+-------------------------------------------------------------------+
  !| The end of the subroutine aceval.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This funcion computes Specific Heats.                             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 07-Jan-2021  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  real(8) function gammarmix(tmp,spc)
    !
    ! arguments
    real(8),intent(in),optional :: tmp,spc(:)
    !
#ifdef COMB
    ! local data
    real(8) :: cpcmix,gamrgc
    !
    call cpeval(tmp=tmp,spc=spc,cp=cpcmix)
    !
    gammarmix=cpcmix/(cpcmix-rgcmix(spc))
    !
    return
    !
#endif 
    !
  end function gammarmix
  !+-------------------------------------------------------------------+
  !| The end of the subroutine gammarmix.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This funcion computes pressure using thermal EoS.                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-Aug-2020  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
!   function thermeos(den,tmp,spc,prs)
!     !
!     ! arguments
!     real(8) :: thermeos
!     real(8),intent(in) :: spc(:)
!     real(8),intent(in),optional :: den,tmp,prs
!     !
! #ifdef COMB
!     if(present(den).and.present(tmp)) then
!       thermeos = den*tmp*rgcmix(spc)
!     elseif(present(den).and.present(prs)) then
!       thermeos = prs/den/rgcmix(spc)
!     elseif(present(tmp).and.present(prs)) then
!       thermeos = prs/tmp/rgcmix(spc)
!     else
!       stop ' !! unable to use thermal EoS !!'
!     endif
!     !
! #endif 
!     !
!   end function thermeos
  !+-------------------------------------------------------------------+
  !| The end of the subroutine thermeos.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This funcion computes mixture mixture constant.                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Aug-2020  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  function rgcmix_scar(spc) result(vout)
    !
    real(8),intent(in) :: spc(:)
    real(8) :: vout
    !
#ifdef COMB
    vout = sum(spc(:)*rgspec(:))
#endif
    !
  end function rgcmix_scar
  !
  function rgcmix_1d(spc,dim) result(vout)
    !
    integer,intent(in) :: dim
    real(8),intent(in) :: spc(:,:)
    real(8) :: vout(dim)
    !local
    integer :: j
    !
#ifdef COMB
    do j=1,dim
      vout(j) = sum(spc(j,:)*rgspec(:))
    enddo 
#endif
    !
  end function rgcmix_1d
  !
  function rgcmix_3d(spc,dim) result(vout)
    ! 
    integer,intent(in) :: dim(3)
    real(8),intent(in) :: spc(:,:,:,:)
    real(8) :: vout(dim(1),dim(2),dim(3))
    !
    !local
    integer :: i,j,k
    !
#ifdef COMB
    !
    do i=1,dim(1)
      do j=1,dim(2)
        do k=1,dim(3)
          vout(i,j,k) = sum(spc(i,j,k,:)*rgspec(:))
        enddo
      enddo 
    enddo 
    !
#endif
    !
  end function rgcmix_3d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine rgcmix.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine computes mixture specific heat capacity           |
  !| and mixture gas constant, or internal energy depending on input   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-Aug-2020  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  subroutine cpeval(tmp,spc,eng,cp,ke)
    !
    ! arguments
    real(8),intent(in) :: spc(:),tmp
    real(8),intent(in),optional :: ke
    real(8),intent(out),optional :: cp,eng
    !
#ifdef COMB
    ! local data
    integer :: is,it,icp,jt
    real(8) :: fornow
    logical :: lctr_cp=.true.
    !
    it=1
    !
    species: do is=1,num_species
      !
      if(present(cp)) then
        !
        if(present(eng)) &
            stop ' !! Conflict - both cp and e are given!!'
            !
        if(lctr_cp .and. ctrflag) then
          !
          call setState_TPY(mixture,tmp,prefgb,spc(:))
          cp=cp_mass(mixture)
          exit species
          !
        else
          !
          ! if(is==1)cp=0.0d0
          ! fornow=amascp(ncpoly(it,is),it,is)
          ! do icp=ncpom1(it,is),1,-1
          !   fornow=amascp(icp,it,is)+fornow*tmp
          ! enddo
          ! cp=cp+spc(is)*fornow
          !
        endif
        !
      elseif(present(eng)) then
        !
        if(lctr_cp .and. ctrflag) then
          !
          call setState_TPY(mixture,tmp,prefgb,spc(:))
          eng=intEnergy_mass(mixture)
          exit species
          !
        else
          !
          ! if(is==1)eng=0.0d0
          ! fornow=amasch(ncpoly(it,is),it,is)
          ! do icp=ncpom1(it,is),1,-1
          !   fornow=amasch(icp,it,is)+fornow*tmp
          ! enddo
          ! fornow=amasch(ncenth(it,is),it,is)+fornow*tmp
          ! !
          ! eng=(eng+spc(is)*fornow)
          !
        endif
        !
      else
        !
        stop ' !! Error - neither cp nor e is given !!'
        !
      endif !present(cp)
      !
    enddo species
    !
    if(present(eng)) then
      if(present(ke)) then
        !
        if(lctr_cp .and. ctrflag) then
          continue
        else
          ! eng=eng-rgcmix(spc)*tmp
        endif
        !
        eng=eng+ke
        !
      else
        stop ' !! kinetic energy not given for total energy !!'
      endif
    endif
    !
#endif
    !
  end subroutine cpeval
  !+-------------------------------------------------------------------+
  !| The end of the subroutine cpeval.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine computes local mixture viscosity, thermal         |
  !| conductivity, and species mass diffusivity.                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-Oct-2020  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  subroutine tranco(den,tmp,cp,mu,lam,rhodi,spc,rhodij)
    !
    use commvar, only: prandtl
    !
    ! arguments
    real(8),intent(in) :: tmp,spc(:),den
    real(8),intent(out),optional :: mu,lam,rhodi(:),cp,rhodij(:,:)
    !
#ifdef COMB
    ! local data
    integer :: js
    real(8) :: lamocp
    !
    call setState_TRY(mixture,tmp,den,spc(:))
    !
    if(present(lam)) lam=thermalConductivity(mixture)
    !
    if(present(mu)) mu=viscosity(mixture)
    !
    if(present(cp)) cp=cp_mass(mixture)
    !
    if(present(rhodi) .or. present(rhodij)) then
      !
      select case(tranmod)
        !
        case('mixav')
          call getMixDiffCoeffs(mixture,rhodi(:))
          rhodi(:)=den*rhodi(:)
          !
        case('multi')
          do js=1,num_species
            call getMultiDiffCoeffs(mixture,js,rhodij(js,:))
            rhodij(js,:)=den*rhodij(js,:)
          enddo
          !
        case default
          ! 
          rhodi(:)=thermalConductivity(mixture) &
                    /cp_mass(mixture)*olewis(:)
      end select
      !
    endif
    !
#endif
    !
  end subroutine tranco
  !+-------------------------------------------------------------------+
  !| The end of the subroutine tranco.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine computes local species enthalpy.                  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-Oct-2020  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  subroutine enthpy(tmp,hi)
    !
    ! arguments
    real(8),intent(in) :: tmp
    real(8),intent(out) :: hi(num_species)
    !
#ifdef COMB
    ! local data
    integer :: is,it,jt,icp
    real(8) :: fornow
    !
    if(.true. .and. ctrflag) then
      !
      call setTemperature(mixture,tmp)
      call getEnthalpies_RT(mixture,hi(:))
      hi(:)=hi(:)*rgspec(:)*tmp
      !
    else
      !
      it=1
      do is=1,num_species
        !
        if(is==1) then
          do jt=1,ntint(is)
            if(tmp>tinthi(ntint(is),is)) then
              print*,' !! ENTH Error - Temperature out of bound!!, T =',tmp
            endif
            if(tmp>tinthi(jt,is)) it=it+1
          enddo
        endif
        !
        fornow=amasch(ncpoly(it,is),it,is)
        do icp=ncpom1(it,is),1,-1
          fornow=amasch(icp,it,is)+fornow*tmp
        enddo
        hi(is)=amasch(ncenth(it,is),it,is)+fornow*tmp
        !
      enddo
      !
    endif
    !
#endif
    !
  end subroutine enthpy
  !+-------------------------------------------------------------------+
  !| The end of the subroutine enthpy.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine computes chemcal reaction rates                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 1-Nov-2020  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  subroutine chemrate(den,tmp,spc,wi)
    !
    ! arguments
    real(8),intent(in) :: den,tmp,spc(:)
    real(8),optional,intent(out) :: wi(:)
    !
#ifdef COMB
    ! local data
    integer :: js,jr,j3b,jsspec,it,jt,icp
    real(8) :: &
      fornow,tbconc,rfwd,rfld,rfln,rbln,gibbs,rbwd,stoi,p_rdc,ftc, &
      const1,const2,rfsr,rftr,wi_molar(num_species),prs
    logical :: flag3by
    !
    if(ctrflag) then
      call setState_TRY(mixture,tmp,den,spc(:))
      call getNetProductionRates(mixture,wi_molar(:))
      wirate(:)=wi_molar(:)*wmolar(:)
      if(present(wi)) wi(:)=wirate(:)
    endif
#endif
    !
  end subroutine chemrate
  !+-------------------------------------------------------------------+
  !| The end of the subroutine chemrate.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine computes heat release rate.                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 19-Nov-2020  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  real(8) function heatrate(den,tmp,spc)
    !
    use commvar, only: num_species
    !
    ! arguments
    real(8),intent(in) :: den,tmp,spc(:)
    !
#ifdef COMB
    ! local data
    real(8) :: hi(num_species)
    !
    call enthpy(tmp,hi)
    call chemrate(den,tmp,spc(:))
    !
    heatrate=-1.d0*sum(wirate(:)*hi(:))
    !
#endif
    !
  end function heatrate
  !+-------------------------------------------------------------------+
  !| The end of the subroutine aceval.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to skip a given number of lines while read.    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-Feb-2021  | Created by Z.X. Chen @ Cambridge                   |
  !+-------------------------------------------------------------------+
  subroutine skipline(fileunit,nlines)
    !
    integer, intent(in) :: fileunit,nlines
    !
    integer :: i,n
    !
    do i=1,nlines
      read(fileunit,*)
    enddo
    !
  end subroutine skipline
  !+-------------------------------------------------------------------+
  !| The end of the subroutine skipline.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is convert between mass and mole fractions.       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 3-Mar-2021  | Created by Z.X. Chen @ Cambridge                    |
  !+-------------------------------------------------------------------+
  subroutine convertxiyi(fracin,fracout,mode)
    !
    character(len=*), intent(in) :: mode
    real(8), intent(in) :: fracin(:)
    real(8), intent(out) :: fracout(:)
    !
#ifdef COMB
    !
    if(mode=='X2Y') then
      fracout(:)=(fracin(:)*wmolar(:))/sum(fracin(:)*wmolar(:))
    elseif(mode=='Y2X') then
      fracout(:)=(fracin(:)/wmolar(:))/sum(fracin(:)/wmolar(:))
    else
      stop ' !!Error - wrong mode given in convertxiyi!!'
    endif
    !
#endif
    !
  end subroutine convertxiyi
  !+-------------------------------------------------------------------+
  !| The end of the subroutine skipline.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine finds a species index using its name string.      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 9-Mar-2021  | Created by Z.X. Chen @ Cambridge                    |
  !+-------------------------------------------------------------------+
  integer function spcindex(spcname)
    !
    character(len=*), intent(in) :: spcname
    !
#ifdef COMB
    spcindex=speciesIndex(mixture,spcname)
#endif
    !
  end function spcindex
  !+-------------------------------------------------------------------+
  !| The end of the subroutine spcindex.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is find a species name string using its index.    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 9-Mar-2021  | Created by Z.X. Chen @ Cambridge                    |
  !+-------------------------------------------------------------------+
  character(len=10) function spcname(spcindex)
    !
    integer, intent(in) :: spcindex
    !
#ifdef COMB
    call getSpeciesName(mixture,spcindex,spcname)
#endif
    !
  end function spcname
  !+-------------------------------------------------------------------+
  !| The end of the subroutine spcindex.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is implicit Euler ODE solver.                     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-Nov-2021  | Created by Z.X. Chen @ Peking University           |
  !+-------------------------------------------------------------------+
  subroutine imp_euler_ode(den,tmp,spc,dt)
    !
    ! arguments
    real(8), intent(in) :: den,tmp,dt
    real(8), intent(inout) :: spc(:)
    !
#ifdef COMB
    ! local data
    integer :: is,iter
    real(8) :: differ,sumy,cmolrates(num_species),dmolrates(num_species) &
              ,spc1(num_species),spc2(num_species)  
    !
    differ=1.d0
    iter=0
    spc1(:)=spc(:)
    spc2(:)=0.d0
    !
    do while(differ>1.d-6)
      !
      differ=0.d0
      sumy=0.d0
      call setState_TRY(mixture,tmp,den,spc1(:))
      call getCreationRates(mixture,cmolrates(:))
      call getDestructionRates(mixture,dmolrates(:))
      do is=1,num_species
        if(spc1(is)<1.d-15 .and. dmolrates(is)<1.d-15) then 
          spc2(is)=(spc(is)+dt*cmolrates(is)*wmolar(is)/den)
        else
          spc2(is)=spc1(is) &
                   *(spc(is)+dt*cmolrates(is)*wmolar(is)/den) &
                   /(spc1(is)+dt*dmolrates(is)*wmolar(is)/den)
        endif
        if(spc1(is)>1.d-9) &
          differ=max(differ,abs(log10(spc2(is))/log10(spc1(is))-1.d0))
        sumy=sumy+max(spc2(is),0.d0)
      enddo
      !
      ! sumy=1.d0
      spc1(1:num_species)=max(spc2(1:num_species)/sumy,0.d0)
      !
      if(iter<1000) then
        iter=iter+1        
      else
        print*,tmp,spc,differ
        print*,' !!Error - implicit Euler ODE failed!!'
        stop
      endif
      !
    enddo
    !
    spc(:)=spc1(:)
    ! 
#endif
    !
  end subroutine imp_euler_ode
  !+-------------------------------------------------------------------+
  !| The end of the subroutine spcindex.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine calculates termpature from energy.                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-Jan-2022  | Created by Z.X. Chen @ Peking University           |
  !+-------------------------------------------------------------------+
  subroutine temperature_calc(tmp,den,spc,eint)
    !
    real(8), intent(in) :: den,spc(:),eint
    real(8), intent(inout) :: tmp
    !
#ifdef COMB
    call setState_TRY(mixture,tmp,den,spc(:))
    call setState_UV(mixture,eint,1.d0/den)
    tmp=temperature(mixture)
#endif
    !
  end subroutine temperature_calc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine spcindex.                               |
  !+-------------------------------------------------------------------+
  !
end module thermchem
!+---------------------------------------------------------------------+
!| The end of the module thermchem.                                    |
!+---------------------------------------------------------------------+
