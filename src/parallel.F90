!+---------------------------------------------------------------------+
!| This module contains subroutines relating to parallelisation.       |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module parallel
  !
  use mpi
  use commvar,   only : im,jm,km,hm,ia,ja,ka,lihomo,ljhomo,lkhomo,     &
                        npdci,npdcj,npdck,lfftk
  use stlaio,  only: get_unit
  !
  implicit none
  !
  interface bcast
    module procedure bcast_cha
    module procedure bcast_chafix_ary
    module procedure bcast_int
    module procedure bcast_r8
    module procedure bcast_log
    module procedure bcast_int_ary
    module procedure bcast_int_ary2
    module procedure bcast_int_ary3
    module procedure bcast_r8_ary
    module procedure bcast_r8_ary2
    module procedure bcast_r8_ary3
    !
    module procedure bcast_solid_ary
  end interface
  !
  interface dataswap
    module procedure array3d_sendrecv_log
    module procedure array3d_sendrecv_int
    module procedure array3d_sendrecv
    module procedure array4d_sendrecv
    module procedure array5d_sendrecv
  end interface
  !
  interface datasync
    module procedure array3d_sync
  end interface
  !
  interface pmerg
    module procedure pmerg_sbou
  end interface
  !
  interface psum
    module procedure psum_int
    module procedure psum_int_ary
    module procedure psum_r8
    module procedure psum_r8_ary
    module procedure psum_r8_ary_2d
  end interface
  !
  interface pmax
    module procedure pmax_int
    module procedure pmax_r8
  end interface
  !
  interface pmin
    module procedure pmin_int
    module procedure pmin_r8
  end interface
  !
  interface por
    module procedure por_log
    module procedure por_log_array1d
  end interface
  !
  interface ptabupd
    module procedure updatable_int
    module procedure updatable_int_a2a_v
    module procedure updatable_rel_a2a_v
    module procedure updatable_rel2d_a2a_v
  end interface
  !
  interface pgather
    module procedure pgather_rel2d_array_rank
    module procedure pgather_rel2d_array
    module procedure pgather_int2d_array
    module procedure pgather_int1d_array
    module procedure pgather_cha1_array
  end interface
  !
  interface pscatter
    module procedure pscatter_rel2d
  end interface
  !
  interface rmdup
     module procedure rmdup_i4
     module procedure rmdup_i4_2d
  end interface
  !
  integer :: mpirank,mpisize,mpirankmax
  integer :: isize,jsize,ksize,irkm,jrkm,krkm,irk,jrk,krk,ig0,jg0,kg0, &
             irk_islice
  integer :: mpileft,mpiright,mpidown,mpiup,mpifront,mpiback,mpitag
  character(len=8) :: mpirankname
  logical :: lio
  integer :: status(mpi_status_size)
  integer :: mpi_imin,mpi_jmin,group_imin,group_jmin,mpi_group_world,  &
             mpi_islice,group_islice
  character(mpi_max_processor_name) :: processor_name
  !
  contains
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this subroutine is used to initial the mpi environment.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! writen by fang jian, 2010-08-26.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpiinitial
    !
    integer :: ierr,namelen
    !
    call mpi_init(ierr)
    !
    call mpi_comm_rank(mpi_comm_world,mpirank,ierr)
    ! get rank number
    call mpi_comm_size(mpi_comm_world,mpisize,ierr)
    ! get total number of nodes
    call mpi_get_processor_name(processor_name,namelen,ierr)
    ! ger processor names
    mpirankmax=mpisize-1
    !
    if(mpisize==1) then
      print*,' ** the program is working in serial environment!'
    else
      if(mpirank==0)  then
        write(*,'(A,I0,A)')' ** the program is working with ',mpisize, &
                                                                ' ranks'
      endif
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
    end if
    !
    mpitag=100
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    ! print*,' ** rank',mpirank,'is on processor',processor_name
    !
  end subroutine mpiinitial
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The end of the subroutine mpiinitial
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to assign the distributions size on the ranks
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2013-7-19.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpisizedis
    !
    logical :: lallo
    integer :: nfactor(1:30)
    integer :: i,j,k,n,n1,n2,n3,nsize,kaa
    integer(8) :: nvar1,nvar2,nvar3
    !
    integer :: ndims
    !
    if(mpisize==1) then
       !
      isize=1
      jsize=1
      ksize=1
      !
      return
      !
    endif
    !
    if(ja==0 .and. ka==0) then
      ndims=1
    elseif(ka==0) then
      ndims=2
    else
      ndims=3
    endif
    !
    if(ka==0) then 
      kaa=ka+1
    else
      kaa=ka
    endif
    !
    if(mpirank==0) then
      !
      write(*,'(A,I0)')'  ** dimension of domain: ',ndims
      !
      lallo=.false.
      !
      nvar2=2**30
      !
      if(mpisize>2**30) then
        print*,' !! Number of processors is too large.'
        stop
      end if
      !
      n=0
      do i=1,mpisize
        !
        if(mod(mpisize,i)==0) then
            n=n+1
            nfactor(n)=mpisize/i
        end if
        !
      end do
      print*,' ** Number of factor is ',n
      !
      do n1=1,n
      do n2=1,n
      do n3=1,n
        !
        !
        if(ndims==1) then
          !
          nsize=nfactor(n1)*nfactor(n2)*nfactor(n3)
          if(nfactor(n2) .ne. 1) cycle
          if(nfactor(n3) .ne. 1) cycle
          !
          ! print*,n1,n2,n3,'|',nfactor(n1),nfactor(n2),nfactor(n3)
        elseif(ndims==2 .or. lfftk) then
          !
          nsize=nfactor(n1)*nfactor(n2)*nfactor(n3)
          if(nfactor(n3) .ne. 1) cycle
          if(nfactor(n2) == 1) cycle
          if(nfactor(n1) == 1) cycle
          !
        else
          !
          if(nfactor(n3)>1 .and. nfactor(n2)>1 .and. nfactor(n1)>1) then
            nsize=nfactor(n1)*nfactor(n2)*nfactor(n3)
          else
           ! for 3D  calculation, the mpisize at k direction should not
           ! be 1
           cycle
          end if
          !
        end if
        !
        if(nsize .eq. mpisize) then
          !
          nvar1=ja*kaa*nfactor(n1)
          !
          nvar1=nvar1+ia*kaa*nfactor(n2)
          !
          nvar1=nvar1+ia*ja*nfactor(n3)
          !
          ! print*,nvar1,nvar2,'-',nfactor(n1),nfactor(n2),nfactor(n3)
          !
          if(nvar1<nvar2) then
            !
            nvar2=nvar1
            !
            isize=nfactor(n1)
            jsize=nfactor(n2)
            ksize=nfactor(n3)
            !
            print*,' ** isze,jsize,ksize= ',isize,jsize,ksize
            !
            lallo=.true.
            !
          end if
          !
        end if
        !
      end do
      end do
      end do
      !
      write(*,'(3(A,I0))')'  ** mpi size= ',isize,' x ',jsize,' x ',ksize
      !
      if(.not. lallo) then
        !
        print*,' !! Size of ranks can not be allocated !!'
        stop
        !
      end if
      !
    end if
    !
  end subroutine mpisizedis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The end of the subroutine mpisizedis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to preprocess the input file to generate
  ! files for each process.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2010-08-26.
  ! Improved by Fang Jian, 2013-07-24.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine parapp
    !
    integer, allocatable, dimension(:,:,:) :: imp,jmp,kmp
    integer :: ierr,nproc
    integer :: n,n1,n2
    integer :: nsi,nsj,nsk,fh
    !
    ! integer :: n,n1,n2,ns,i2,j2,k2,nar1,nsi,nsj,nsk,ina,jna,kna
    ! integer :: ispa,jspa,kspa
    ! real(8) :: var1
    ! character(len=32) :: cmdcptemp
    ! !
    ! !
    if(mpirank==0) then
      !
      nproc=isize*jsize*ksize
      !
      irkm=isize-1
      jrkm=jsize-1
      krkm=ksize-1
      !
      if(nproc .ne. mpisize) then
        print*,' !! error!, nproc /= mpisize !!'
        print*,' ** nproc=',nproc
        print*,' ** mpisize=',mpisize
        stop
      end if
      !
      ! allocate block size for each rank
      allocate( imp(0:isize-1,0:jsize-1,0:ksize-1),                  &
                jmp(0:isize-1,0:jsize-1,0:ksize-1),                  &
                kmp(0:isize-1,0:jsize-1,0:ksize-1)                   )
      !
      n1=ia/isize
      n2=mod(ia,isize)
      !
      do n=0,isize-1
        imp(n,:,:)=n1
      end do
      !
      if(n2 .ne. 0) then
        !
        do n=isize-1,isize-n2,-1
          imp(n,:,:)=imp(n,:,:)+1
        end do
        !
      end if
      !
      n1=ja/jsize
      n2=mod(ja,jsize)
      !
      do n=0,jsize-1
        jmp(:,n,:)=n1
      end do
      !
      if(n2 .ne. 0) then
        !
        do n=jsize-1,jsize-n2,-1
          jmp(:,n,:)=jmp(:,n,:)+1
        end do
        !
      end if
      !
      n1=ka/ksize
      n2=mod(ka,ksize)
      !
      do n=0,ksize-1
        kmp(:,:,n)=n1
      end do
      !
      if(n2 .ne. 0) then
        !
        do n=ksize-1,ksize-n2,-1
          kmp(:,:,n)=kmp(:,:,n)+1
        end do
        !
      end if
      !
      fh=get_unit()
      !
      open(fh,file='datin/parallel.info',form='formatted')
      write(fh,"(3(A9,1x))")'isize','jsize','ksize'
      write(fh,"(3(I9,1x))")isize,jsize,ksize
      write(fh,"(10(A9,1x))")'Rank','Irk','Jrk','Krk','IM','JM','KM',  &
                                                      'I0','J0','K0'
      nsi=0
      nsj=0
      nsk=0
      do krk=0,ksize-1
      do jrk=0,jsize-1
      do irk=0,isize-1
        !
        n=krk*(isize*jsize)+jrk*isize+irk
        !
        ! Output file of rank information.
        write(fh,"(10(I9,1x))")n,irk,jrk,krk,imp(irk,jrk,krk),       &
                            jmp(irk,jrk,krk),kmp(irk,jrk,krk),       &
                            nsi,nsj,nsk
        !
        nsi=nsi+imp(irk,0,0)
      enddo
        nsi=0
        nsj=nsj+jmp(0,jrk,0)
      enddo
        nsi=0
        nsj=0
        nsk=nsk+kmp(0,0,krk)
      enddo
      !
      close(fh)
      print*,' << parallel.info ... done !'
      !
    endif
    !
    call bcast(irkm)
    call bcast(jrkm)
    call bcast(krkm)
    !
    call bcast(isize)
    call bcast(jsize)
    call bcast(ksize)
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    if(mpirank==0) print*,' ** parallel processing ... done.'
    !
    return
    !
  end subroutine parapp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! the end of the subroutine parapp.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to initilize parallel parameters
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2010-08-27.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine parallelini
    !
    use commvar,   only : is,ie,js,je,ks,ke,islice,jslice,kslice
    !
    ! local data
    integer :: n,nn,nsize1,nsize,ierr,fh,ni,nj,nk,newsize,nrk
    integer,allocatable:: ntemp(:,:),nrank(:,:,:)
    integer,allocatable :: rank_use(:)
    integer,allocatable,dimension(:) :: nrect,irkg,jrkg,krkg,          &
                                        img,jmg,kmg,i0g,j0g,k0g
    !
    !
    if(mpirank>=0 .and. mpirank<=9)then
      write(mpirankname,'(7h.000000,i1)') mpirank
    elseif(mpirank>=10 .and. mpirank<=99)then
      write(mpirankname,'(6h.00000,i2)') mpirank
    elseif(mpirank>=100 .and. mpirank<=999)then
      write(mpirankname,'(5h.0000,i3)') mpirank
    elseif(mpirank>=1000 .and. mpirank<=9999)then
      write(mpirankname,'(4h.000,i4)') mpirank
    elseif(mpirank>=10000 .and. mpirank<=99999)then
      write(mpirankname,'(3h.00,i5)') mpirank
    elseif(mpirank>=100000 .and. mpirank<=999999)then
      write(mpirankname,'(2h.0,i6)') mpirank
    elseif(mpirank>=1000000 .and. mpirank<=9999999)then
      write(mpirankname,'(1h.,i7)') mpirank
    else
      print *, ' !! Error: rank number not in the range [0,9999999]'
      stop
    end if
    !
    lio=.false.
    !
    allocate( nrank(0:irkm,0:jrkm,0:krkm) )
    !
    if(mpirank==0) then
      !
      allocate( ntemp(1:11,0:mpirankmax)                               )
      allocate( irkg(0:mpirankmax),jrkg(0:mpirankmax),                 &
                krkg(0:mpirankmax),img(0:mpirankmax),                  &
                jmg(0:mpirankmax),kmg(0:mpirankmax),                   &
                i0g(0:mpirankmax),j0g(0:mpirankmax),                   &
                k0g(0:mpirankmax)                                     )
      !
      fh=get_unit()
      open(fh,file='datin/parallel.info',form='formatted')
      read(fh,"(//)")
      do n=0,mpirankmax
        !
        read(fh,*)ntemp(1,n),ntemp(2,n),ntemp(3,n),ntemp(4,n),         &
                  ntemp(5,n),ntemp(6,n),ntemp(7,n),ntemp(8,n),         &
                  ntemp(9,n),ntemp(10,n)
        !
        irkg(n)=ntemp(2,n)
        jrkg(n)=ntemp(3,n)
        krkg(n)=ntemp(4,n)
        img(n)=ntemp(5,n)
        jmg(n)=ntemp(6,n)
        kmg(n)=ntemp(7,n)
        i0g(n)=ntemp(8,n)
        j0g(n)=ntemp(9,n)
        k0g(n)=ntemp(10,n)
        !
        nrank(irkg(n),jrkg(n),krkg(n))=n
        !
      end do
      close(fh)
      print*,' >> parallel.info ...done.'
      !
      ntemp(11,:)=0
      !
      nsize=(img(0)+1)*(jmg(0)+1)*(kmg(0)+1)
      ntemp(11,0)=1
      !
      do n=1,mpirankmax
        !
        nsize1=(img(n)+1)*(jmg(n)+1)*(kmg(n)+1)
        if(nsize1<nsize) then
          nsize=nsize1
          ntemp(11,n)=1
          ntemp(11,n-1)=0
        end if
        !
      end do
      !
      ! do n=1,mpirankmax
      !   call MPI_SEND(ntemp(1,n),11,MPI_INTEGER,n,mpitag,              &
      !                                               mpi_comm_world,ierr)
      ! end do
      ! mpitag=mpitag+1
      ! !
      ! deallocate( ntemp )
    endif
    !
    allocate( nrect(1:11) )
    !
    call mpi_scatter(ntemp(1,0),11,mpi_integer,                        &
                     nrect(1),  11,mpi_integer,                        &
                     0,mpi_comm_world,ierr)
      !
      !
    !   call MPI_RECV(nrect(1),11,MPI_INTEGER,0,mpitag,                  &
    !                                          mpi_comm_world,status,ierr)
    !   mpitag=mpitag+1
    !   !
    ! print*,mpirank,'|',nrect(1)
    !
    irk=nrect(2)
    jrk=nrect(3)
    krk=nrect(4)
    im=nrect(5)
    jm=nrect(6)
    km=nrect(7)
    ig0=nrect(8)
    jg0=nrect(9)
    kg0=nrect(10)
    !
    !
    if(nrect(11)==1) lio=.true.
    !
    if(mpirank==0) deallocate( ntemp )
    !
    deallocate( nrect )
    !
    !
    if(lihomo) then
      !
      is=0
      ie=im
      !
      if(irk==0) then
        mpileft=krk*(isize*jsize)+jrk*isize+irkm
        mpiright=krk*(isize*jsize)+jrk*isize+irk+1
      elseif(irk==irkm) then
        mpileft=krk*(isize*jsize)+jrk*isize+irk-1
        mpiright=krk*(isize*jsize)+jrk*isize+0
      else
        mpileft=krk*(isize*jsize)+jrk*isize+irk-1
        mpiright=krk*(isize*jsize)+jrk*isize+irk+1
      end if
      !
      npdci=3
      !
    else
      !
      if(irk==0) then
        is=1
        ie=im
        !
        npdci=1
        !
        mpileft=MPI_PROC_NULL
        mpiright=krk*(isize*jsize)+jrk*isize+irk+1
      elseif(irk==irkm) then
        is=0
        ie=im-1
        !
        npdci=2
        mpileft=krk*(isize*jsize)+jrk*isize+irk-1
        mpiright=MPI_PROC_NULL
      else
        is=0
        ie=im
        !
        npdci=3
        !
        mpileft =krk*(isize*jsize)+jrk*isize+irk-1
        mpiright=krk*(isize*jsize)+jrk*isize+irk+1
      end if
      !
    end if
    !
    if(ljhomo) then
      !
      js=0
      je=jm
      !
      npdcj=3
      !
      if(jsize==1) then
        mpidown=MPI_PROC_NULL
        mpiup  =MPI_PROC_NULL
      else
        if(jrk==0) then
          mpidown=krk*(isize*jsize)+   jrkm*isize+irk
          mpiup  =krk*(isize*jsize)+(jrk+1)*isize+irk
        elseif(jrk==jrkm) then
          mpidown=krk*(isize*jsize)+(jrk-1)*isize+irk
          mpiup  =krk*(isize*jsize)+      0*isize+irk
        else
          mpidown=krk*(isize*jsize)+(jrk-1)*isize+irk
          mpiup  =krk*(isize*jsize)+(jrk+1)*isize+irk
        end if
      endif
      !
    else
      !
      if(jsize==1) then
        mpidown=MPI_PROC_NULL
        mpiup  =MPI_PROC_NULL
      else
        if(jrk==0) then
          js=1
          je=jm
          !
          npdcj=1
          !
          mpidown=MPI_PROC_NULL
          mpiup  =krk*(isize*jsize)+(jrk+1)*isize+irk
        elseif(jrk==jrkm) then
          js=0
          je=jm-1
          !
          npdcj=2
          !
          mpidown=krk*(isize*jsize)+(jrk-1)*isize+irk
          mpiup=MPI_PROC_NULL
        else
          js=0
          je=jm
          !
          npdcj=3
          !
          mpidown=krk*(isize*jsize)+(jrk-1)*isize+irk
          mpiup  =krk*(isize*jsize)+(jrk+1)*isize+irk
        end if
      endif
      !
    end if
    !
    if(lkhomo) then
      !
      ks=0
      ke=km
      !
      npdck=3
      !
      if(ksize==1) then
        mpiback=MPI_PROC_NULL
        mpifront=MPI_PROC_NULL
      else
        if(krk==0) then
          mpiback=    krkm*(isize*jsize)+jrk*isize+irk
          mpifront=(krk+1)*(isize*jsize)+jrk*isize+irk
        elseif(krk==krkm) then
          mpiback=(krk-1)*(isize*jsize)+jrk*isize+irk
          mpifront=  0.d0*(isize*jsize)+jrk*isize+irk
        else
          mpiback= (krk-1)*(isize*jsize)+jrk*isize+irk
          mpifront=(krk+1)*(isize*jsize)+jrk*isize+irk
        end if
      end if
      !
    else
      !
      if(ksize==1) then
        mpiback=MPI_PROC_NULL
        mpifront=MPI_PROC_NULL
      else
        !
        if(krk==0) then
          ks=1
          ke=km
          !
          npdck=1
          !
          mpiback=MPI_PROC_NULL
          mpifront=(krk+1)*(isize*jsize)+jrk*isize+irk
        elseif(krk==krkm) then
          ks=0
          ke=km-1
          !
          npdck=2
          !
          mpiback=(krk-1)*(isize*jsize)+jrk*isize+irk
          mpifront=MPI_PROC_NULL
        else
          ks=0
          ke=km
          !
          npdck=3
          !
          mpiback=(krk-1)*(isize*jsize)+jrk*isize+irk
          mpifront=(krk+1)*(isize*jsize)+jrk*isize+irk
        end if
        !
      end if
      !
    end if 
    !
    if(lio) write(*,'(A,I0,A)')'  ** rank ',mpirank,' is the I/O rank'
    !
    call bcast(nrank)
    !
    allocate(rank_use(jsize*ksize))
    rank_use=-1
    n=0
    do nk=0,ksize-1
    do nj=0,jsize-1
      n=n+1
      rank_use(n)=nrank(0,nj,nk)
    end do
    end do
    !
    call mpi_comm_group(mpi_comm_world,mpi_group_world,ierr)
    call mpi_group_incl(mpi_group_world,size(rank_use),rank_use,group_imin,ierr)
    call mpi_comm_create(mpi_comm_world,group_imin,mpi_imin,ierr)
    !
    if(irk==0) call mpi_comm_size(mpi_imin,newsize,ierr)
    if(mpirank==0) write(*,'(A,I0)') &
      '  ** new communicator: mpi_imin  ... created, size: ',newsize
    !
    ! set sub communicator for islice
    irk_islice=-1
    if(mpirank==0) then
      !
      do nrk=0,mpirankmax
        !
        if(islice>=i0g(nrk) .and. islice<i0g(nrk)+img(nrk)) then
          irk_islice=irkg(nrk)
          print*,' ** islice is at irk=',irk_islice
          exit
        endif
        !
      enddo
      !
    endif
    !
    call bcast(irk_islice)
    !
    if(irk_islice>0) then
      !
      rank_use=-1
      n=0
      do nk=0,ksize-1
      do nj=0,jsize-1
        n=n+1
        rank_use(n)=nrank(irk_islice,nj,nk)
      end do
      end do
      !
      call mpi_group_incl(mpi_group_world,size(rank_use),rank_use,group_islice,ierr)
      call mpi_comm_create(mpi_comm_world,group_islice,mpi_islice,ierr)
      if(irk==irk_islice) then
        call mpi_comm_size(mpi_islice,newsize,ierr)
        if(jrk==0 .and. krk==0) then
          print*,' ** new communicator: mpi_islice  ... created, size: ',newsize
        endif
      endif
      !
    endif
    !
    deallocate(rank_use)
    ! end of set sub communicator for islice
    !
    allocate(rank_use(isize*ksize))
    n=0
    do nk=0,ksize-1
    do ni=0,isize-1
      n=n+1
      rank_use(n)=nrank(ni,0,nk)
    end do
    end do
    !
    call mpi_group_incl(mpi_group_world,size(rank_use),rank_use,group_jmin,ierr)
    call mpi_comm_create(mpi_comm_world,group_jmin,mpi_jmin,ierr)
    !
    if(jrk==0) call mpi_comm_size(mpi_jmin,newsize,ierr)
    if(mpirank==0) write(*,'(A,I0)') &
      '  ** new communicator: mpi_jmin  ... created, size: ',newsize
    !
    deallocate(rank_use)
    !
    deallocate(nrank)
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
  end subroutine parallelini
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The end of the subroutine parallelini
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to barrier all ranks                      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 19-July-2019: Created by J. Fang @ STFC Daresbury Laboratory      |
  !+-------------------------------------------------------------------+
  subroutine mpibar
    !
    integer :: ierr
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
  end subroutine mpibar
  !+-------------------------------------------------------------------+
  !| The end of the subroutine mpibar.                                 |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to finalise mpi and stop the program.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 19-July-2019: Created by J. Fang @ STFC Daresbury Laboratory      |
  !+-------------------------------------------------------------------+
  subroutine mpistop
    !
    integer :: ierr
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    call mpi_finalize(ierr)
    !
    if(lio) print*,' ** The job is done!'
    !
    stop
    !
  end subroutine mpistop
  !+-------------------------------------------------------------------+
  !| The end of the subroutine mpistop.                                |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to broadcase variables and arraies to all |
  !|  ranks.                                                           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 19-July-2019: Created by J. Fang @ STFC Daresbury Laboratory      |
  !+-------------------------------------------------------------------+
  subroutine bcast_cha(vario)
    !
    ! arguments
    character(len=*),intent(inout) :: vario
    !
    ! local data
    integer :: len,ierr
    !
    len=len_trim(vario)
    call MPI_BCAST(len,1,mpi_integer,0,mpi_comm_world,ierr)
    !
    call MPI_BCAST(vario,len,mpi_character,0,mpi_comm_world,ierr)
    !
    vario=vario(1:len)
    !
  end subroutine bcast_cha
  !
  subroutine bcast_chafix_ary(vario)
    !
    ! arguments
    character(len=*),intent(inout) :: vario(:)
    !
    ! local data
    integer :: nsize,ierr
    !
    nsize=size(vario)
    !
    call MPI_BCAST(vario,nsize,mpi_character,0,mpi_comm_world,ierr)
    !
  end subroutine bcast_chafix_ary
  !
  subroutine bcast_int(vario,comm)
    !
    ! arguments
    integer,intent(inout) :: vario
    integer,intent(in),optional :: comm
    !
    ! local data
    integer :: ierr
    !
    if(present(comm)) then
      call mpi_bcast(vario,1,mpi_integer,0,comm,ierr)
    else
      call mpi_bcast(vario,1,mpi_integer,0,mpi_comm_world,ierr)
    endif
    !
    !
  end subroutine bcast_int
  !
  subroutine bcast_int_ary(vario)
    !
    ! arguments
    integer,intent(inout) :: vario(:)
    !
    ! local data
    integer :: nsize,ierr
    !
    nsize=size(vario)
    !
    call MPI_BCAST(vario,nsize,mpi_integer,0,mpi_comm_world,ierr)
    !
  end subroutine bcast_int_ary
  !
  subroutine bcast_int_ary2(vario)
    !
    ! arguments
    integer,intent(inout) :: vario(:,:)
    !
    ! local data
    integer :: nsize,ierr
    !
    nsize=size(vario,1)*size(vario,2)
    !
    call MPI_BCAST(vario,nsize,mpi_integer,0,mpi_comm_world,ierr)
    !
  end subroutine bcast_int_ary2
  !
  subroutine bcast_int_ary3(vario)
    !
    ! arguments
    integer,intent(inout) :: vario(:,:,:)
    !
    ! local data
    integer :: nsize,ierr
    !
    nsize=size(vario,1)*size(vario,2)*size(vario,3)
    !
    call MPI_BCAST(vario,nsize,mpi_integer,0,mpi_comm_world,ierr)
    !
  end subroutine bcast_int_ary3
  !
  subroutine bcast_r8(vario,comm)
    !
    ! arguments
    real(8),intent(inout) :: vario
    integer,intent(in),optional :: comm
    !
    ! local data
    integer :: ierr
    !
    if(present(comm)) then
      call mpi_bcast(vario,1,mpi_real8,0,comm,ierr)
    else
      call mpi_bcast(vario,1,mpi_real8,0,mpi_comm_world,ierr)
    endif
    !
    !
  end subroutine bcast_r8
  !
  subroutine bcast_r8_ary(vario)
    !
    ! arguments
    real(8),intent(inout) :: vario(:)
    !
    ! local data
    integer :: nsize,ierr
    !
    nsize=size(vario)
    !
    call MPI_BCAST(vario,nsize,mpi_real8,0,mpi_comm_world,ierr)
    !
  end subroutine bcast_r8_ary
  !
  subroutine bcast_r8_ary2(vario,rank)
    !
    ! arguments
    real(8),intent(inout) :: vario(:,:)
    integer,intent(in),optional :: rank
    !
    ! local data
    integer :: nsize,ierr,root
    !
    nsize=size(vario,1)*size(vario,2)
    !
    if(present(rank)) then
      root=rank
    else
      root=0
    endif
    !
    call MPI_BCAST(vario,nsize,mpi_real8,root,mpi_comm_world,ierr)
    !
  end subroutine bcast_r8_ary2
  !
  subroutine bcast_r8_ary3(vario)
    !
    ! arguments
    real(8),intent(inout) :: vario(:,:,:)
    !
    ! local data
    integer :: nsize,ierr
    !
    nsize=size(vario,1)*size(vario,2)*size(vario,3)
    !
    call MPI_BCAST(vario,nsize,mpi_real8,0,mpi_comm_world,ierr)
    !
  end subroutine bcast_r8_ary3
  !
  subroutine bcast_log(vario)
    !
    ! arguments
    logical,intent(inout) :: vario
    !
    ! local data
    integer :: ierr
    !
    call MPI_BCAST(vario,1,mpi_logical,0,mpi_comm_world,ierr)
    !
  end subroutine bcast_log
  !+-------------------------------------------------------------------+
  !| The end of the subroutine bcast.                                  |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to broadcase solid array.                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 02-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine bcast_solid_ary(vario)
    !
    use commtype,  only : solid,triangle
    !
    ! arguments
    type(solid),allocatable,intent(inout) :: vario(:)
    !
    ! local data
    integer :: nsize,ierr,js,jf
    !
    nsize=size(vario)
    !
    call bcast_int(nsize)
    !
    if(mpirank.ne.0) then
      if(allocated(vario)) deallocate(vario)
      allocate(vario(1:nsize))
    endif
    !
    do js=1,nsize
      !
      call bcast_cha(vario(js)%name)
      call bcast_r8_ary(vario(js)%xmin)
      call bcast_r8_ary(vario(js)%xmax)
      call bcast_r8_ary(vario(js)%xref)
      call bcast_r8_ary(vario(js)%xcen)
      call bcast_int(vario(js)%num_face)
      call bcast_int(vario(js)%num_edge)
      !
      if(mpirank.ne.0) then
        call vario(js)%alloface()
        call vario(js)%alloedge()
      endif
      !
      do jf=1,vario(js)%num_face
        call bcast_r8_ary(vario(js)%face(jf)%a)
        call bcast_r8_ary(vario(js)%face(jf)%b)
        call bcast_r8_ary(vario(js)%face(jf)%c)
        call bcast_r8_ary(vario(js)%face(jf)%normdir)
        call bcast_r8(vario(js)%face(jf)%area)
      enddo
      !
      do jf=1,vario(js)%num_edge
        call bcast_r8_ary(vario(js)%edge(jf)%a)
        call bcast_r8_ary(vario(js)%edge(jf)%b)
        call bcast_r8_ary(vario(js)%edge(jf)%normdir)
      enddo
      !
    enddo
    !
    !
  end subroutine bcast_solid_ary
  !+-------------------------------------------------------------------+
  !| The end of the subroutine bcast_solid_ary.                        |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to update and synchronise the table.      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 30-September-2019: Created by J. Fang @ STFC Daresbury Laboratory |
  !+-------------------------------------------------------------------+
  subroutine updatable_int(var,table,offset,debug)
    !
    ! arguments
    integer,intent(in) :: var
    integer,optional,intent(out) :: table(0:mpirankmax)
    integer,optional,intent(out) :: offset
    logical,intent(in),optional :: debug
    !
    ! local data
    integer :: ierr,i
    integer :: vta(0:mpirankmax)
    logical :: ldebug
    !
    if(present(debug)) then
      ldebug=debug
    else
      ldebug=.false.
    endif
    !
    call mpi_allgather(var,1,mpi_integer,                              &
                       vta,1,mpi_integer,mpi_comm_world,ierr)
    !
    if(present(table)) then
      table=vta
    endif
    !
    if(present(offset)) then
      !
      if(mpirank==0) then
        offset=0
      else
        !
        offset=0
        do i=0,mpirank-1
          offset=offset+vta(i)
        enddo
        !
      endif
      !
    endif
    !
  end subroutine updatable_int
  !
  subroutine updatable_int_a2a_v(datasend,datarecv,                    &
                                 sendtabl,recvtabl)
    !
    ! arguments
    integer,intent(in) :: datasend(:)
    integer,allocatable,intent(out) :: datarecv(:)
    integer,intent(in) :: sendtabl(0:),recvtabl(0:)
    !
    ! local data
    integer :: nvar,ierr,incode,jrank,recvsize
    integer,allocatable :: senddispls(:),recvdispls(:)
    !
    ! start
    !
    nvar=1
    !
    allocate(senddispls(0:mpirankmax),recvdispls(0:mpirankmax))
    !
    senddispls=0
    recvdispls=0
    do jrank=1,mpirankmax
      senddispls(jrank)=senddispls(jrank-1)+sendtabl(jrank-1)
      recvdispls(jrank)=recvdispls(jrank-1)+recvtabl(jrank-1)
    enddo
    recvsize=recvdispls(mpirankmax)+recvtabl(mpirankmax)
    allocate(datarecv(recvsize))
    !
    call mpi_alltoallv(datasend, sendtabl, senddispls, mpi_integer, &
                       datarecv, recvtabl, recvdispls, mpi_integer, &
                       mpi_comm_world, ierr)
    !
    deallocate(senddispls,recvdispls)
    !
  end subroutine updatable_int_a2a_v
  !
  subroutine updatable_rel_a2a_v(datasend,datarecv,                    &
                                 sendtabl,recvtabl)
    !
    ! arguments
    real(8),intent(in) :: datasend(:)
    real(8),allocatable,intent(out) :: datarecv(:)
    integer,intent(in) :: sendtabl(0:),recvtabl(0:)
    !
    ! local data
    integer :: nvar,ierr,incode,jrank,recvsize
    integer,allocatable :: senddispls(:),recvdispls(:)
    !
    ! start
    !
    nvar=1
    !
    allocate(senddispls(0:mpirankmax),recvdispls(0:mpirankmax))
    !
    senddispls=0
    recvdispls=0
    do jrank=1,mpirankmax
      senddispls(jrank)=senddispls(jrank-1)+sendtabl(jrank-1)
      recvdispls(jrank)=recvdispls(jrank-1)+recvtabl(jrank-1)
    enddo
    recvsize=recvdispls(mpirankmax)+recvtabl(mpirankmax)
    allocate(datarecv(recvsize))
    !
    call mpi_alltoallv(datasend, sendtabl, senddispls, mpi_real8, &
                       datarecv, recvtabl, recvdispls, mpi_real8, &
                       mpi_comm_world, ierr)
    !
    deallocate(senddispls,recvdispls)
    !
  end subroutine updatable_rel_a2a_v
  !
  subroutine updatable_rel2d_a2a_v(datasend,datarecv,                 &
                                   sendtabl,recvtabl)
    !
    ! arguments
    real(8),intent(in) :: datasend(:,:)
    real(8),allocatable,intent(out) :: datarecv(:,:)
    integer,intent(in) :: sendtabl(0:),recvtabl(0:)
    !
    ! local data
    integer :: nvar,ierr,incode,jrank,recvsize,n2size
    integer,allocatable :: senddispls(:),recvdispls(:),                &
                           sendtabl2(:),recvtabl2(:)
    !
    ! start
    !
    nvar=1
    !
    n2size=size(datasend,2)
    !
    allocate(senddispls(0:mpirankmax),recvdispls(0:mpirankmax),        &
             sendtabl2(0:mpirankmax),recvtabl2(0:mpirankmax))
    !
    sendtabl2=sendtabl*n2size
    recvtabl2=recvtabl*n2size
    senddispls=0
    recvdispls=0
    do jrank=1,mpirankmax
      senddispls(jrank)=senddispls(jrank-1)+sendtabl2(jrank-1)
      recvdispls(jrank)=recvdispls(jrank-1)+recvtabl2(jrank-1)
    enddo
    recvsize=recvdispls(mpirankmax)+recvtabl(mpirankmax)
    allocate(datarecv(recvsize/n2size,n2size))
    !
    call mpi_alltoallv(datasend, sendtabl2, senddispls, mpi_real8, &
                       datarecv, recvtabl2, recvdispls, mpi_real8, &
                       mpi_comm_world, ierr)
    !
    deallocate(senddispls,recvdispls,sendtabl2,recvtabl2)
    !
  end subroutine updatable_rel2d_a2a_v
  !+-------------------------------------------------------------------+
  !| The end of the subroutine updatable_int.                          |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to gather arraies.                        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-Jul-2021: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine pgather_rel2d_array_rank(array,data,rank)
    !
    ! arguments
    real(8),intent(in) :: array(:,:)
    integer,intent(in) :: rank
    real(8),intent(out),allocatable :: data(:,:)
    !
    !
    ! local data
    integer :: displs(0:mpirankmax),counts(0:mpirankmax)
    integer :: ierr,jrank,nrecv(1),nsize,dim1,dim2
    !
    ! get the size of data to gather
    dim1=size(array,1)
    dim2=size(array,2)
    !
    dim1=pmax(dim1)
    !
    nrecv(1)=size(array)
    !
    call mpi_gather(nrecv, 1, mpi_integer, counts, 1, mpi_integer,  &
                    rank,mpi_comm_world, ierr)
    ! 
    if(mpirank==rank) then
      !
      displs(0)=0
      do jrank=1,mpirankmax
        displs(jrank)=displs(jrank-1)+counts(jrank-1)
      enddo
      !
      nsize=displs(mpirankmax)+counts(mpirankmax)
      !
      allocate(data(dim1,nsize/dim1))
      !
    endif
    !
    call mpi_gatherv(array, nrecv(1), mpi_real8,              &
                      data, counts, displs, mpi_real8, rank,  &
                      mpi_comm_world, ierr)
    return
    !
  end subroutine pgather_rel2d_array_rank
  !
  subroutine pgather_int1d_array(array,data)
    !
    ! arguments
    integer,intent(in) :: array(:)
    integer,intent(out),allocatable :: data(:)
    !
    !
    ! local data
    integer :: displs(0:mpirankmax),counts(0:mpirankmax)
    integer :: ierr,jrank,nrecv(1),nsize
    !
    ! get the size of data to gather
    nrecv(1)=size(array)
    !
    call mpi_allgather(nrecv, 1, mpi_integer, counts, 1, mpi_integer,  &
                       mpi_comm_world, ierr)
    !
    displs(0)=0
    do jrank=1,mpirankmax
      displs(jrank)=displs(jrank-1)+counts(jrank-1)
    enddo
    !
    nsize=displs(mpirankmax)+counts(mpirankmax)
    allocate(data(1:nsize))
    !
    call mpi_allgatherv(array, nrecv(1), mpi_integer,       &
                        data,  counts, displs, mpi_integer, &
                        mpi_comm_world, ierr)
    !
  end subroutine pgather_int1d_array
  !
  subroutine pgather_rel2d_array(array,data)
    !
    ! arguments
    real(8),intent(in) :: array(:,:)
    real(8),intent(out),allocatable :: data(:,:)
    !
    !
    ! local data
    integer :: displs(0:mpirankmax),counts(0:mpirankmax)
    integer :: ierr,jrank,nrecv(1),nsize,dim1,dim2
    !
    ! get the size of data to scatter
    dim1=size(array,1)
    dim2=size(array,2)
    !
    nrecv(1)=size(array)
    !
    dim1=pmax(dim1)
    !
    call mpi_allgather(nrecv, 1, mpi_integer, counts, 1, mpi_integer,  &
                    mpi_comm_world, ierr)
    !
    displs(0)=0
    do jrank=1,mpirankmax
      displs(jrank)=displs(jrank-1)+counts(jrank-1)
    enddo
    !
    nsize=displs(mpirankmax)+counts(mpirankmax)
    !
    allocate(data(dim1,nsize/dim1))
    !
    call mpi_allgatherv(array, nrecv(1), mpi_real8,       &
                        data,  counts, displs, mpi_real8, &
                        mpi_comm_world, ierr)
    !
  end subroutine pgather_rel2d_array
  !
  subroutine pgather_int2d_array(array,data)
    !
    ! arguments
    integer,intent(in) :: array(:,:)
    integer,intent(out),allocatable :: data(:,:)
    !
    !
    ! local data
    integer :: displs(0:mpirankmax),counts(0:mpirankmax)
    integer :: ierr,jrank,nrecv(1),nsize,dim1,dim2
    !
    ! get the size of data to scatter
    dim1=size(array,1)
    dim2=size(array,2)
    !
    nrecv(1)=size(array)
    !
    dim1=pmax(dim1)
    !
    call mpi_allgather(nrecv, 1, mpi_integer, counts, 1, mpi_integer,  &
                    mpi_comm_world, ierr)
    !
    displs(0)=0
    do jrank=1,mpirankmax
      displs(jrank)=displs(jrank-1)+counts(jrank-1)
    enddo
    !
    nsize=displs(mpirankmax)+counts(mpirankmax)
    !
    allocate(data(dim1,nsize/dim1))
    !
    call mpi_allgatherv(array, nrecv(1), mpi_integer,       &
                        data,  counts, displs, mpi_integer, &
                        mpi_comm_world, ierr)
    !
  end subroutine pgather_int2d_array
  !
  subroutine pgather_cha1_array(array,data)
    !
    ! arguments
    character(len=1),intent(in) :: array(:)
    character(len=1),intent(out),allocatable :: data(:)
    !
    !
    ! local data
    integer :: displs(0:mpirankmax),counts(0:mpirankmax)
    integer :: ierr,jrank,nrecv(1),nsize,dim1
    !
    ! get the size of data to scatter
    dim1=size(array)
    !
    nrecv(1)=dim1
    !
    call mpi_allgather(nrecv, 1, mpi_integer, counts, 1, mpi_integer,  &
                    mpi_comm_world, ierr)
    !
    displs(0)=0
    do jrank=1,mpirankmax
      displs(jrank)=displs(jrank-1)+counts(jrank-1)
    enddo
    !
    nsize=displs(mpirankmax)+counts(mpirankmax)
    !
    allocate(data(nsize))
    !
    call mpi_allgatherv(array, nrecv(1), mpi_character,       &
                        data,  counts, displs, mpi_character, &
                        mpi_comm_world, ierr)
    !
  end subroutine pgather_cha1_array
  !+-------------------------------------------------------------------+
  !| The end of the subroutine pgather.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to scatter data according to the table.   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 15-Jul-2021: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine pscatter_rel2d(array,data,table,offset,rank)
    !
    ! arguments
    real(8),intent(in) :: array(:,:)
    integer,intent(in) :: table(0:),rank
    integer,intent(in),optional :: offset(0:)
    real(8),intent(out),allocatable :: data(:,:)
    !
    ! local data
    integer :: displs(0:mpirankmax),counts(0:mpirankmax)
    integer :: ierr,jrank,dim1
    !
    dim1=size(array,1)
    !
    counts=dim1*table
    !
    if(present(offset)) then
      displs=offset
    else
      displs(0)=0
      do jrank=1,mpirankmax
        displs(jrank)=displs(jrank-1)+counts(jrank-1)
      enddo
    endif
    !
    allocate(data(dim1,table(mpirank)))
    !
    call mpi_scatterv(array, counts, displs,  mpi_real8,       &
                      data,  counts(mpirank), mpi_real8, rank, &
                      mpi_comm_world,ierr)
    !
    return
    !
  end subroutine pscatter_rel2d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine pscatter.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to merge arraies from differnet ranks and |
  !| broadcast it.                                                     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 15-Jul-2021: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine pmerg_sbou(var,nvar,vmerg)
    !
    use commtype,  only : sboun
    !
    ! arguments
    type(sboun),allocatable,intent(in) :: var(:)
    type(sboun),allocatable,intent(out) :: vmerg(:)
    integer,intent(in),optional :: nvar
    !
    ! local data
    integer :: nsize,ntotle,jb
    real(8),allocatable :: sb_rea(:,:),sbt_rea(:,:)
    integer,allocatable :: sb_int(:,:),sbt_int(:,:)
    character(len=1),allocatable :: sb_cha(:),sbt_cha(:)
    !
    integer :: ns_tab(0:mpirankmax)
    !
    if(present(nvar)) then
      nsize=nvar
    else
      nsize=size(var)
    endif
    !
    call ptabupd(var=nsize,table=ns_tab)
    !
    allocate(sb_rea(3,nsize),sb_int(3,nsize),sb_cha(nsize))
    !
    ! merge x 
    do jb=1,nsize
      sb_rea(:,jb)=var(jb)%x(:)
    enddo
    !
    call pgather(sb_rea,sbt_rea)
    !
    ntotle=psum(nsize)
    !
    allocate(vmerg(1:ntotle))
    !
    do jb=1,ntotle
      vmerg(jb)%x(:)=sbt_rea(:,jb)
    enddo
    !
    ! end of merge x 
    !
    ! merge normdir 
    do jb=1,nsize
      sb_rea(:,jb)=var(jb)%normdir(:)
    enddo
    !
    call pgather(sb_rea,sbt_rea)
    !
    do jb=1,ntotle
      vmerg(jb)%normdir(:)=sbt_rea(:,jb)
    enddo
    !
    ! end of merge normdir
    !
    ! merge ximag 
    do jb=1,nsize
      sb_rea(:,jb)=var(jb)%ximag(:)
    enddo
    !
    call pgather(sb_rea,sbt_rea)
    !
    do jb=1,ntotle
      vmerg(jb)%ximag(:)=sbt_rea(:,jb)
    enddo
    !
    ! end of merge ximag
    !
    ! merge dis2ghost and dist2image 
    do jb=1,nsize
      sb_rea(1,jb)=var(jb)%dis2image
      sb_rea(2,jb)=var(jb)%dis2ghost
      sb_rea(3,jb)=0.d0
    enddo
    !
    call pgather(sb_rea,sbt_rea)
    !
    do jb=1,ntotle
      vmerg(jb)%dis2image=sbt_rea(1,jb)
      vmerg(jb)%dis2ghost=sbt_rea(2,jb)
    enddo
    !
    ! end of merge ximag
    !
    ! merge igh 
    do jb=1,nsize
      sb_int(:,jb)=var(jb)%igh(:)
    enddo
    !
    call pgather(sb_int,sbt_int)
    !
    do jb=1,ntotle
      vmerg(jb)%igh(:)=sbt_int(:,jb)
    enddo
    ! end of merge igh
    !
    ! merge nodetype 
    do jb=1,nsize
      sb_cha(jb)=var(jb)%nodetype
    enddo
    !
    call pgather(sb_cha,sbt_cha)
    !
    do jb=1,ntotle
      vmerg(jb)%nodetype=sbt_cha(jb)
    enddo
    ! end of merge igh
    !
    ! print*,mpirank,'|',nsize
    ! !
    ! if(mpirank==4) then
    !   print*,mpirank,'|',sb_int(:,1)
    ! endif
    ! !
    ! if(mpirank==0) then
    !   print*,mpirank,'|',sbt_int(:,717)
    ! endif
    !
    return
    !
  end subroutine pmerg_sbou
  !+-------------------------------------------------------------------+
  !| The end of the subroutine pmerg_sbou.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is used to sum a number across ranks                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 01-October-2019: Created by J. Fang @ STFC Daresbury Laboratory   |
  !+-------------------------------------------------------------------+
  integer function  psum_int(var)
    !
    ! arguments
    integer,intent(in) :: var
    !
    ! local data
    integer :: ierr
    !
    call mpi_allreduce(var,psum_int,1,mpi_integer,mpi_sum,             &
                                                    mpi_comm_world,ierr)
    !
  end function psum_int
  !
  function psum_int_ary(var) result(varsum)
    !
    ! arguments
    integer,intent(in) :: var(:)
    integer,allocatable :: varsum(:)
    !
    ! local data
    integer :: ierr,nsize
    !
    nsize=size(var)
    !
    allocate(varsum(nsize))
    !
    call mpi_allreduce(var,varsum,nsize,mpi_integer,mpi_sum,           &
                                                    mpi_comm_world,ierr)
    !
    return
    !
  end function psum_int_ary
  !
  real(8) function  psum_r8(var)
    !
    ! arguments
    real(8),intent(in) :: var
    !
    ! local data
    integer :: ierr
    !
    call mpi_allreduce(var,psum_r8,1,mpi_real8,mpi_sum,               &
                                                    mpi_comm_world,ierr)
    !
  end function psum_r8
  !
  function psum_r8_ary(var) result(varsum)
    !
    ! arguments
    real(8),intent(in) :: var(:)
    real(8),allocatable :: varsum(:)
    !
    ! local data
    integer :: ierr,nsize
    !
    nsize=size(var)
    !
    allocate(varsum(nsize))
    !
    call mpi_allreduce(var,varsum,nsize,mpi_real8,mpi_sum,             &
                                                    mpi_comm_world,ierr)
    !
    return
    !
  end function psum_r8_ary
  !
  function psum_r8_ary_2d(var) result(varsum)
    !
    ! arguments
    real(8),intent(in) :: var(:,:)
    real(8),allocatable :: varsum(:,:)
    !
    ! local data
    integer :: ierr,nsize1,nsize2
    !
    nsize1=size(var,1)
    nsize2=size(var,2)
    !
    allocate(varsum(nsize1,nsize2))
    !
    call mpi_allreduce(var,varsum,nsize1*nsize2,mpi_real8,mpi_sum,     &
                                                    mpi_comm_world,ierr)
    !
    return
    !
  end function psum_r8_ary_2d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine psum.                                   |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is used to get the max number of all ranks          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 07-October-2019: Created by J. Fang @ STFC Daresbury Laboratory   |
  !| 09-January-2020: add a real8 interface, by J. Fang @ STFC         |
  !|                   Daresbury Laboratory                            |
  !+-------------------------------------------------------------------+
  integer function  pmax_int(var)
    !
    ! arguments
    integer,intent(in) :: var
    !
    ! local data
    integer :: ierr
    !
    call mpi_allreduce(var,pmax_int,1,mpi_integer,mpi_max,             &
                                                    mpi_comm_world,ierr)
    !
  end function pmax_int
  !
  real(8) function  pmax_r8(var)
    !
    ! arguments
    real(8),intent(in) :: var
    !
    ! local data
    integer :: ierr
    !
    call mpi_allreduce(var,pmax_r8,1,mpi_real8,mpi_max,               &
                                                    mpi_comm_world,ierr)
    !
  end function pmax_r8
  !+-------------------------------------------------------------------+
  !| The end of the subroutine pmax.                                   |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| Thes functions are used to get the min number of all ranks        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-January-2020: Created by J. Fang @ STFC Daresbury Laboratory   |
  !+-------------------------------------------------------------------+
  integer function  pmin_int(var)
    !
    ! arguments
    integer,intent(in) :: var
    !
    ! local data
    integer :: ierr
    !
    call mpi_allreduce(var,pmin_int,1,mpi_integer,mpi_min,             &
                                                    mpi_comm_world,ierr)
    !
  end function pmin_int
  !
  real(8) function  pmin_r8(var)
    !
    ! arguments
    real(8),intent(in) :: var
    !
    ! local data
    integer :: ierr
    !
    call mpi_allreduce(var,pmin_r8,1,mpi_real8,mpi_min,               &
                                                    mpi_comm_world,ierr)
    !
  end function pmin_r8
  !+-------------------------------------------------------------------+
  !| The end of the subroutine pmin.                                   |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this subroutine is used to message passing grids coordinates.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! writen by fang jian, 2010-09-02.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gridsendrecv
    !
    use commarray, only : x
    !
    ! local data
    integer :: ierr
    integer :: i,j,k,m,n
    integer :: ncou
    real(8), allocatable, dimension(:,:,:,:) :: sendbuf1,sendbuf2,     &
                                                recvbuf1,recvbuf2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! sendbuf1: send buffer
    ! sendbuf2: send buffer
    ! recvbuf1: received buffer
    ! recvbuf2: received buffer
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! print*,mpirank,'|',mpileft,mpiright,mpidown,mpiup,mpiback,mpifront
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ncou=(jm+1)*(km+1)*hm*3
    !
    allocate( sendbuf1(1:hm,0:jm,0:km,1:3),                            &
              sendbuf2(1:hm,0:jm,0:km,1:3),                            &
              recvbuf1(1:hm,0:jm,0:km,1:3),                            &
              recvbuf2(1:hm,0:jm,0:km,1:3)                             )
    !
    ! pack the left/right send buffer
    do n=1,hm
      sendbuf1(n,0:jm,0:km,1:3)=x(n,0:jm,0:km,1:3)   -x(0,0:jm,0:km,1:3)
      sendbuf2(n,0:jm,0:km,1:3)=x(im-n,0:jm,0:km,1:3)-x(im,0:jm,0:km,1:3)
    enddo
    !
    ! message passing
    call mpi_sendrecv(sendbuf1, ncou, mpi_real8, mpileft,  mpitag,     &
                      recvbuf1, ncou, mpi_real8, mpiright, mpitag,     &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    call mpi_sendrecv(sendbuf2, ncou, mpi_real8, mpiright, mpitag,    &
                      recvbuf2, ncou, mpi_real8, mpileft,  mpitag,     &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    !! unpack the received left the packet
    if(mpileft==MPI_PROC_NULL) then
      !
      do n=1,hm
        x(-n,0:jm,0:km,1:3)=2.d0*x(0,0:jm,0:km,1:3)-x(n,0:jm,0:km,1:3) ! even
      enddo
      !
    else
      !
      do n=1,hm
        x(-n,0:jm,0:km,1:3)=recvbuf2(n,0:jm,0:km,1:3)+x(0,0:jm,0:km,1:3)
      enddo
      !
    end if
    !
    ! unpack the received right the packet
    if(mpiright==MPI_PROC_NULL) then
      !
      do n=1,hm
        x(im+n,0:jm,0:km,1:3)=2.d0*x(im,0:jm,0:km,1:3)-x(im-n,0:jm,0:km,1:3)
      enddo
      !
    else
      do n=1,hm
        x(im+n,0:jm,0:km,1:3)=recvbuf1(n,0:jm,0:km,1:3)+x(im,0:jm,0:km,1:3)
      enddo
      !
    end if
    !
    deallocate( sendbuf1,sendbuf2,recvbuf1,recvbuf2 )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Finish message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(jsize>1) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(km+1)*hm*3
      !
      allocate( sendbuf1(0:im,1:hm,0:km,1:3),                          &
                sendbuf2(0:im,1:hm,0:km,1:3),                          &
                recvbuf1(0:im,1:hm,0:km,1:3),                          &
                recvbuf2(0:im,1:hm,0:km,1:3)                           )
      !
      ! pack the down send/up buffer
      do n=1,hm
        sendbuf1(0:im,n,0:km,1:3)=x(0:im,   n,0:km,1:3)-x(0:im, 0,0:km,1:3)
        sendbuf2(0:im,n,0:km,1:3)=x(0:im,jm-n,0:km,1:3)-x(0:im,jm,0:km,1:3)
      enddo
      !
      ! Message passing
      call mpi_sendrecv(sendbuf1,ncou,mpi_real8,mpidown,mpitag,        &
                        recvbuf1,ncou,mpi_real8,mpiup,mpitag,          &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sendbuf2,ncou,mpi_real8,mpiup,mpitag,          &
                        recvbuf2,ncou,mpi_real8,mpidown,mpitag,        &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      ! unpack the received up the packet
      if(mpidown==MPI_PROC_NULL) then
        !
        do n=1,hm
          x(0:im,-n,0:km,1:3)=2.d0*x(0:im,0,0:km,1:3)-x(0:im,n,0:km,1:3)
        enddo
        !
      else
        !
        do n=1,hm
          x(0:im,-n,0:km,1:3)=recvbuf2(0:im,n,0:km,1:3)+x(0:im,0,0:km,1:3)
        enddo
        !
      end if
      !
      ! unpack the received down the packet
      if(mpiup==MPI_PROC_NULL) then
        do n=1,hm
          x(0:im,jm+n,0:km,1:3)=2.d0*x(0:im,jm,0:km,1:3)-x(0:im,jm-n,0:km,1:3)
        enddo
      else
        !
        do n=1,hm
          x(0:im,jm+n,0:km,1:3)=recvbuf1(0:im,n,0:km,1:3)+x(0:im,jm,0:km,1:3)
        enddo
        !
      endif
      !
      deallocate( sendbuf1,sendbuf2,recvbuf1,recvbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
      !
      if(ja==0) then
        !
        do j=-hm,hm
          x(0:im,j,0:km,1)=x(0:im,0,0:km,1)
          x(0:im,j,0:km,2)=x(0:im,0,0:km,2)+real(j,8)
          x(0:im,j,0:km,3)=x(0:im,0,0:km,3)
        enddo
        !
      endif
      !
    endif
    !
    if(ksize>1) then
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(jm+1)*hm*3
      !
      allocate( sendbuf1(0:im,0:jm,1:hm,1:3),                          &
                sendbuf2(0:im,0:jm,1:hm,1:3),                          &
                recvbuf1(0:im,0:jm,1:hm,1:3),                          &
                recvbuf2(0:im,0:jm,1:hm,1:3)                           )
      !
      ! pack the back/front send buffer
      do n=1,hm
        sendbuf1(0:im,0:jm,n,1:3)=x(0:im,0:jm,   n,1:3)-x(0:im,0:jm, 0,1:3)
        sendbuf2(0:im,0:jm,n,1:3)=x(0:im,0:jm,km-n,1:3)-x(0:im,0:jm,km,1:3)
      enddo
      !
      ! Message passing
      call MPI_SENDRECV(sendbuf1,ncou,MPI_REAL8,mpiback,mpitag,        &
                        recvbuf1,ncou,MPI_REAL8,mpifront,mpitag,       &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call MPI_SENDRECV(sendbuf2,ncou,MPI_REAL8,mpifront,mpitag,       &
                        recvbuf2,ncou,MPI_REAL8,mpiback,mpitag,        &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      ! unpack the received back the packet
      if(mpiback==MPI_PROC_NULL) then
        !
        ! do n=1,hm
        !   x(0:im,0:jm,-n,1:3)=2.d0*x(0:im,0:jm,0,1:3)-x(0:im,0:jm,n,1:3)
        ! end do
        !
      else
        !
        do n=1,hm
          x(0:im,0:jm,-n,1:3)=recvbuf2(0:im,0:jm,n,1:3)+x(0:im,0:jm,0,1:3)
        end do
        !
      end if
      !
      ! unpack the received front the packet
      if(mpifront==MPI_PROC_NULL) then
        !
        ! do n=1,hm
        !   x(0:im,0:jm,km+n,1:3)=2.d0*x(0:im,0:jm,km,1:3)-x(0:im,0:jm,km-n,1:3)
        ! end do
        !
      else
        !
        do n=1,hm
          x(0:im,0:jm,km+n,1:3)=recvbuf1(0:im,0:jm,n,1:3)+x(0:im,0:jm,km,1:3)
        end do
        !
      end if
      !
      deallocate( sendbuf1,sendbuf2,recvbuf1,recvbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else
      !
      if(ka==0) then
        !
        do k=-hm,hm
          x(0:im,0:jm,k,1:2)=x(0:im,0:jm,0,1:2)
          x(0:im,0:jm,k,3)=x(0:im,0:jm,0,3)+real(k,8)
        enddo
        !
      endif
      !
    end if
    !
    if(lio) print*,' ** grid coordinates swapped'
    !
    return
    !
  end subroutine gridsendrecv
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! the end of the subroutine gridsendrecv.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to swap a 3-D logical tensor.             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| -Jul-2021 | Created by J. Fang @ Warringon.                       |
  !+-------------------------------------------------------------------+
  subroutine array3d_sendrecv_log(var,subtime)
    !
    ! arguments
    logical,intent(inout) :: var(-hm:im+hm,-hm:jm+hm,-hm:km+hm)
    real(8),intent(inout),optional :: subtime
    !
    ! logical data
    integer :: ncou,j,k
    integer :: ierr
    logical,allocatable,dimension(:,:,:) :: sbuf1,sbuf2,rbuf1,rbuf2
    real(8) :: time_beg
    !
    if(present(subtime)) time_beg=ptime()
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! buf1: send buffer
    ! buf2: send buffer
    ! buf1: redevice buffer
    ! buf2: redevice buffer
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ncou=(jm+1)*(km+1)*hm
    !
    allocate( sbuf1(1:hm, 0:jm,0:km),sbuf2(1:hm, 0:jm,0:km),           &
              rbuf1(1:hm, 0:jm,0:km),rbuf2(-hm:-1,0:jm,0:km) )
    !
    if(mpileft .ne. MPI_PROC_NULL) then
      ! pack the left send buffer
      sbuf1(1:hm,0:jm,0:km)=var(1:hm,0:jm,0:km)
    endif
    if(mpiright .ne. MPI_PROC_NULL) then
      ! pack the right send buffer
      sbuf2(1:hm,0:jm,0:km)=var(im-hm:im-1,0:jm,0:km)
    endif
    !
    ! Message passing
    call mpi_sendrecv(sbuf1,ncou,mpi_logical,mpileft, mpitag,            &
                      rbuf1,ncou,mpi_logical,mpiright,mpitag,            &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    call mpi_sendrecv(sbuf2,ncou,mpi_logical,mpiright,mpitag,            &
                      rbuf2,ncou,mpi_logical,mpileft, mpitag,            &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    if(mpiright .ne. MPI_PROC_NULL) then
      !
      ! unpack the received the packet from right
      var(im+1:im+hm,0:jm,0:km)=rbuf1(1:hm,0:jm,0:km)
      !
    end if
      !
    if(mpileft .ne. MPI_PROC_NULL) then
      !
      ! unpack the received the packet from left
      var(-hm:-1,0:jm,0:km)=rbuf2(-hm:-1,0:jm,0:km)
      !
      ! var(0,0:jm,0:km)=0.5d0*( var(0,0:jm,0:km) +                      &
      !                        rbuf2(0,0:jm,0:km) )
    end if
    !
    deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Finish message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(jsize==1 .and. ljhomo) then
      !
      if(ja==0) then 
        do j=-hm,hm
          var(0:im,j,0:km)    =var(0:im,0,0:km)
        enddo
      else
        var(0:im,-hm:-1,0:km)    =var(0:im,jm-hm:jm-1,0:km)
        var(0:im,jm+1:jm+hm,0:km)=var(0:im,1:hm,0:km)
      endif
      !
    else
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(km+1)*hm
      !
      allocate( sbuf1(0:im,1:hm, 0:km),sbuf2(0:im,1:hm, 0:km),         &
                rbuf1(0:im,1:hm, 0:km),rbuf2(0:im,-hm:-1,0:km) )
      !
      if(mpidown .ne. MPI_PROC_NULL) then
        ! pack the upper send buffer
        sbuf1(0:im,1:hm,0:km)=var(0:im,1:hm,0:km)
      endif
      if(mpiup .ne. MPI_PROC_NULL) then
        ! pack the down send buffer
        sbuf2(0:im,1:hm,0:km)=var(0:im,jm-hm:jm-1,0:km)
      end if
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_logical,mpidown,mpitag,             &
                        rbuf1,ncou,mpi_logical,mpiup,  mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_logical,mpiup,  mpitag,             &
                        rbuf2,ncou,mpi_logical,mpidown,mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpiup .ne. MPI_PROC_NULL) then
        ! unpack the received the packet from up
        var(0:im,jm+1:jm+hm,0:km)=rbuf1(0:im,1:hm,0:km)
        !
      endif
      !
      if(mpidown .ne. MPI_PROC_NULL) then
        ! unpack the received the packet from down
        var(0:im,-hm:-1,0:km)=rbuf2(0:im,-hm:-1,0:km) 
        !
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
    endif
    !
    if(ksize==1 .and. lkhomo) then
      !
      if(ka==0) then
        do k=-hm,hm
          var(0:im,0:jm,k)    =var(0:im,0:jm,0)
        enddo
      else
        var(0:im,0:jm,-hm:-1)    =var(0:im,0:jm,km-hm:km-1)
        var(0:im,0:jm,km+1:km+hm)=var(0:im,0:jm,1:hm)
      endif
      !             
    else
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(jm+1)*hm
      !
      allocate( sbuf1(0:im,0:jm, 1:hm),                      &
                sbuf2(0:im,0:jm, 1:hm),                      &
                rbuf1(0:im,0:jm, 1:hm),                      &
                rbuf2(0:im,0:jm,-hm:-1) )
      !
      if(mpiback .ne. MPI_PROC_NULL) then
        ! pack the back send buffer
        sbuf1(0:im,0:jm,1:hm)=var(0:im,0:jm,1:hm)
      endif
      if(mpifront .ne. MPI_PROC_NULL) then
        ! pack the front send buffer
        sbuf2(0:im,0:jm,1:hm)=var(0:im,0:jm,km-hm:km-1)
      endif
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_logical,mpiback, mpitag,          &
                        rbuf1,ncou,mpi_logical,mpifront,mpitag,          &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_logical,mpifront,mpitag,          &
                        rbuf2,ncou,mpi_logical,mpiback, mpitag,          &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpifront .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet from front
        var(0:im,0:jm,km+1:km+hm)=rbuf1(0:im,0:jm,1:hm)
        !
      end if
      !
      if(mpiback .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet back
        var(0:im,0:jm,-hm:-1)=rbuf2(0:im,0:jm,-hm:-1)
        !
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    if(present(subtime)) subtime=subtime+ptime()-time_beg
    !
    return
    !
  end subroutine array3d_sendrecv_log
  !+-------------------------------------------------------------------+
  !| The end of the subroutine array3d_sendrecv_log.                   |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to swap a 3-D integer tensor.             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| -Jul-2021 | Created by J. Fang @ Warringon.                       |
  !+-------------------------------------------------------------------+
  subroutine array3d_sendrecv_int(var,subtime)
    !
    ! arguments
    integer,intent(inout) :: var(-hm:im+hm,-hm:jm+hm,-hm:km+hm)
    real(8),intent(inout),optional :: subtime
    !
    ! logical data
    integer :: ncou,j,k
    integer :: ierr
    integer,allocatable,dimension(:,:,:) :: sbuf1,sbuf2,rbuf1,rbuf2
    real(8) :: time_beg
    !
    if(present(subtime)) time_beg=ptime()
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! buf1: send buffer
    ! buf2: send buffer
    ! buf1: redevice buffer
    ! buf2: redevice buffer
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ncou=(jm+1)*(km+1)*hm
    !
    allocate( sbuf1(1:hm, 0:jm,0:km),sbuf2(1:hm, 0:jm,0:km),           &
              rbuf1(1:hm, 0:jm,0:km),rbuf2(-hm:-1,0:jm,0:km) )
    !
    if(mpileft .ne. MPI_PROC_NULL) then
      ! pack the left send buffer
      sbuf1(1:hm,0:jm,0:km)=var(1:hm,0:jm,0:km)
    endif
    if(mpiright .ne. MPI_PROC_NULL) then
      ! pack the right send buffer
      sbuf2(1:hm,0:jm,0:km)=var(im-hm:im-1,0:jm,0:km)
    endif
    !
    ! Message passing
    call mpi_sendrecv(sbuf1,ncou,mpi_integer,mpileft, mpitag,            &
                      rbuf1,ncou,mpi_integer,mpiright,mpitag,            &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    call mpi_sendrecv(sbuf2,ncou,mpi_integer,mpiright,mpitag,            &
                      rbuf2,ncou,mpi_integer,mpileft, mpitag,            &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    if(mpiright .ne. MPI_PROC_NULL) then
      !
      ! unpack the received the packet from right
      var(im+1:im+hm,0:jm,0:km)=rbuf1(1:hm,0:jm,0:km)
      !
    end if
      !
    if(mpileft .ne. MPI_PROC_NULL) then
      !
      ! unpack the received the packet from left
      var(-hm:-1,0:jm,0:km)=rbuf2(-hm:-1,0:jm,0:km)
      !
      ! var(0,0:jm,0:km)=0.5d0*( var(0,0:jm,0:km) +                      &
      !                        rbuf2(0,0:jm,0:km) )
    end if
    !
    deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Finish message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(jsize==1 .and. ljhomo) then
      !
      if(ja==0) then 
        do j=-hm,hm
          var(0:im,j,0:km)    =var(0:im,0,0:km)
        enddo
      else
        var(0:im,-hm:-1,0:km)    =var(0:im,jm-hm:jm-1,0:km)
        var(0:im,jm+1:jm+hm,0:km)=var(0:im,1:hm,0:km)
      endif
      !
    else
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(km+1)*hm
      !
      allocate( sbuf1(0:im,1:hm, 0:km),sbuf2(0:im,1:hm, 0:km),         &
                rbuf1(0:im,1:hm, 0:km),rbuf2(0:im,-hm:-1,0:km) )
      !
      if(mpidown .ne. MPI_PROC_NULL) then
        ! pack the upper send buffer
        sbuf1(0:im,1:hm,0:km)=var(0:im,1:hm,0:km)
      endif
      if(mpiup .ne. MPI_PROC_NULL) then
        ! pack the down send buffer
        sbuf2(0:im,1:hm,0:km)=var(0:im,jm-hm:jm-1,0:km)
      end if
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_integer,mpidown,mpitag,             &
                        rbuf1,ncou,mpi_integer,mpiup,  mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_integer,mpiup,  mpitag,             &
                        rbuf2,ncou,mpi_integer,mpidown,mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpiup .ne. MPI_PROC_NULL) then
        ! unpack the received the packet from up
        var(0:im,jm+1:jm+hm,0:km)=rbuf1(0:im,1:hm,0:km)
        !
      endif
      !
      if(mpidown .ne. MPI_PROC_NULL) then
        ! unpack the received the packet from down
        var(0:im,-hm:-1,0:km)=rbuf2(0:im,-hm:-1,0:km) 
        !
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
    endif
    !
    if(ksize==1 .and. lkhomo) then
      !
      if(ka==0) then
        do k=-hm,hm
          var(0:im,0:jm,k)    =var(0:im,0:jm,0)
        enddo
      else
        var(0:im,0:jm,-hm:-1)    =var(0:im,0:jm,km-hm:km-1)
        var(0:im,0:jm,km+1:km+hm)=var(0:im,0:jm,1:hm)
      endif
      !             
    else
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(jm+1)*hm
      !
      allocate( sbuf1(0:im,0:jm, 1:hm),                      &
                sbuf2(0:im,0:jm, 1:hm),                      &
                rbuf1(0:im,0:jm, 1:hm),                      &
                rbuf2(0:im,0:jm,-hm:-1) )
      !
      if(mpiback .ne. MPI_PROC_NULL) then
        ! pack the back send buffer
        sbuf1(0:im,0:jm,1:hm)=var(0:im,0:jm,1:hm)
      endif
      if(mpifront .ne. MPI_PROC_NULL) then
        ! pack the front send buffer
        sbuf2(0:im,0:jm,1:hm)=var(0:im,0:jm,km-hm:km-1)
      endif
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_integer,mpiback, mpitag,          &
                        rbuf1,ncou,mpi_integer,mpifront,mpitag,          &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_integer,mpifront,mpitag,          &
                        rbuf2,ncou,mpi_integer,mpiback, mpitag,          &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpifront .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet from front
        var(0:im,0:jm,km+1:km+hm)=rbuf1(0:im,0:jm,1:hm)
        !
      end if
      !
      if(mpiback .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet back
        var(0:im,0:jm,-hm:-1)=rbuf2(0:im,0:jm,-hm:-1)
        !
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    if(present(subtime)) subtime=subtime+ptime()-time_beg
    !
    return
    !
  end subroutine array3d_sendrecv_int
  !+-------------------------------------------------------------------+
  !| The end of the subroutine array3d_sendrecv_int.                   |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to swap a 3-D tensor.                     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-Feb-2021 | Created by J. Fang @ Warringon.                     |
  !+-------------------------------------------------------------------+
  subroutine array3d_sendrecv(var,subtime)
    !
    ! arguments
    real(8),intent(inout) :: var(-hm:im+hm,-hm:jm+hm,-hm:km+hm)
    real(8),intent(inout),optional :: subtime
    !
    ! logical data
    integer :: ncou,j,k
    integer :: ierr
    real(8),allocatable,dimension(:,:,:) :: sbuf1,sbuf2,rbuf1,rbuf2
    real(8) :: time_beg
    !
    if(present(subtime)) time_beg=ptime()
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! buf1: send buffer
    ! buf2: send buffer
    ! buf1: redevice buffer
    ! buf2: redevice buffer
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ncou=(jm+1)*(km+1)*hm
    !
    allocate( sbuf1(1:hm, 0:jm,0:km),sbuf2(1:hm, 0:jm,0:km),           &
              rbuf1(1:hm, 0:jm,0:km),rbuf2(-hm:-1,0:jm,0:km) )
    !
    if(mpileft .ne. MPI_PROC_NULL) then
      ! pack the left send buffer
      sbuf1(1:hm,0:jm,0:km)=var(1:hm,0:jm,0:km)
    endif
    if(mpiright .ne. MPI_PROC_NULL) then
      ! pack the right send buffer
      sbuf2(1:hm,0:jm,0:km)=var(im-hm:im-1,0:jm,0:km)
    endif
    !
    ! Message passing
    call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpileft, mpitag,            &
                      rbuf1,ncou,mpi_real8,mpiright,mpitag,            &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpiright,mpitag,            &
                      rbuf2,ncou,mpi_real8,mpileft, mpitag,            &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    if(mpiright .ne. MPI_PROC_NULL) then
      !
      ! unpack the received the packet from right
      var(im+1:im+hm,0:jm,0:km)=rbuf1(1:hm,0:jm,0:km)
      !
    end if
      !
    if(mpileft .ne. MPI_PROC_NULL) then
      !
      ! unpack the received the packet from left
      var(-hm:-1,0:jm,0:km)=rbuf2(-hm:-1,0:jm,0:km)
      !
      ! var(0,0:jm,0:km)=0.5d0*( var(0,0:jm,0:km) +                      &
      !                        rbuf2(0,0:jm,0:km) )
    end if
    !
    deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Finish message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(jsize==1 .and. ljhomo) then
      !
      if(ja==0) then 
        do j=-hm,hm
          var(0:im,j,0:km)    =var(0:im,0,0:km)
        enddo
      else
        var(0:im,-hm:-1,0:km)    =var(0:im,jm-hm:jm-1,0:km)
        var(0:im,jm+1:jm+hm,0:km)=var(0:im,1:hm,0:km)
      endif
      !
    else
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(km+1)*hm
      !
      allocate( sbuf1(0:im,1:hm, 0:km),sbuf2(0:im,1:hm, 0:km),         &
                rbuf1(0:im,1:hm, 0:km),rbuf2(0:im,-hm:-1,0:km) )
      !
      if(mpidown .ne. MPI_PROC_NULL) then
        ! pack the upper send buffer
        sbuf1(0:im,1:hm,0:km)=var(0:im,1:hm,0:km)
      endif
      if(mpiup .ne. MPI_PROC_NULL) then
        ! pack the down send buffer
        sbuf2(0:im,1:hm,0:km)=var(0:im,jm-hm:jm-1,0:km)
      end if
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpidown,mpitag,             &
                        rbuf1,ncou,mpi_real8,mpiup,  mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpiup,  mpitag,             &
                        rbuf2,ncou,mpi_real8,mpidown,mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpiup .ne. MPI_PROC_NULL) then
        ! unpack the received the packet from up
        var(0:im,jm+1:jm+hm,0:km)=rbuf1(0:im,1:hm,0:km)
        !
      endif
      !
      if(mpidown .ne. MPI_PROC_NULL) then
        ! unpack the received the packet from down
        var(0:im,-hm:-1,0:km)=rbuf2(0:im,-hm:-1,0:km) 
        !
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
    endif
    !
    if(ksize==1 .and. lkhomo) then
      !
      if(ka==0) then
        do k=-hm,hm
          var(0:im,0:jm,k)    =var(0:im,0:jm,0)
        enddo
      else
        var(0:im,0:jm,-hm:-1)    =var(0:im,0:jm,km-hm:km-1)
        var(0:im,0:jm,km+1:km+hm)=var(0:im,0:jm,1:hm)
      endif
      !             
    else
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(jm+1)*hm
      !
      allocate( sbuf1(0:im,0:jm, 1:hm),                      &
                sbuf2(0:im,0:jm, 1:hm),                      &
                rbuf1(0:im,0:jm, 1:hm),                      &
                rbuf2(0:im,0:jm,-hm:-1) )
      !
      if(mpiback .ne. MPI_PROC_NULL) then
        ! pack the back send buffer
        sbuf1(0:im,0:jm,1:hm)=var(0:im,0:jm,1:hm)
      endif
      if(mpifront .ne. MPI_PROC_NULL) then
        ! pack the front send buffer
        sbuf2(0:im,0:jm,1:hm)=var(0:im,0:jm,km-hm:km-1)
      endif
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpiback, mpitag,          &
                        rbuf1,ncou,mpi_real8,mpifront,mpitag,          &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpifront,mpitag,          &
                        rbuf2,ncou,mpi_real8,mpiback, mpitag,          &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpifront .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet from front
        var(0:im,0:jm,km+1:km+hm)=rbuf1(0:im,0:jm,1:hm)
        !
      end if
      !
      if(mpiback .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet back
        var(0:im,0:jm,-hm:-1)=rbuf2(0:im,0:jm,-hm:-1)
        !
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    if(present(subtime)) subtime=subtime+ptime()-time_beg
    !
    return
    !
  end subroutine array3d_sendrecv
  !+-------------------------------------------------------------------+
  !| The end of the subroutine array3d_sendrecv.                       |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to synconize the interface variables by   |
  !| swap and averageing                                               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-Jun-2021 | Created by J. Fang @ Warringon.                     |
  !+-------------------------------------------------------------------+
  subroutine array3d_sync(var,subtime)
    !
    ! arguments
    real(8),intent(inout) :: var(0:im,0:jm,0:km)
    real(8),intent(inout),optional :: subtime
    !
    ! logical data
    integer :: ncou,j,k
    integer :: ierr
    real(8),allocatable,dimension(:,:) :: sbuf1,sbuf2,rbuf1,rbuf2
    real(8) :: time_beg
    !
    if(present(subtime)) time_beg=ptime()
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! buf1: send buffer
    ! buf2: send buffer
    ! buf1: redevice buffer
    ! buf2: redevice buffer
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ncou=(jm+1)*(km+1)
    !
    allocate( sbuf1(0:jm,0:km),sbuf2(0:jm,0:km),           &
              rbuf1(0:jm,0:km),rbuf2(0:jm,0:km) )
    !
    if(mpileft .ne. MPI_PROC_NULL) then
      ! pack the left send buffer
      sbuf1(0:jm,0:km)=var(0,0:jm,0:km)
    endif
    if(mpiright .ne. MPI_PROC_NULL) then
      ! pack the right send buffer
      sbuf2(0:jm,0:km)=var(im,0:jm,0:km)
    endif
    !
    ! Message passing
    call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpileft, mpitag,            &
                      rbuf1,ncou,mpi_real8,mpiright,mpitag,            &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpiright,mpitag,            &
                      rbuf2,ncou,mpi_real8,mpileft, mpitag,            &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    if(mpiright .ne. MPI_PROC_NULL) then
      !
      ! unpack the received the packet from right
      var(im,0:jm,0:km)=0.5d0*(var(im,0:jm,0:km)+rbuf1(0:jm,0:km))
      !
    end if
      !
    if(mpileft .ne. MPI_PROC_NULL) then
      !
      ! unpack the received the packet from left
      var(0,0:jm,0:km)=0.5d0*(var(0,0:jm,0:km)+rbuf2(0:jm,0:km) )
    end if
    !
    deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Finish message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(jsize==1 .and. ljhomo) then
      !
      var(0:im,0,0:km)=0.5d0*(var(0:im,0,0:km)+var(0:im,jm,0:km))
      var(0:im,jm,0:km)=var(0:im,0,0:km)
      !
    else
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(km+1)
      !
      allocate( sbuf1(0:im,0:km),sbuf2(0:im,0:km),                    &
                rbuf1(0:im,0:km),rbuf2(0:im,0:km) )
      !
      if(mpidown .ne. MPI_PROC_NULL) then
        ! pack the upper send buffer
        sbuf1(0:im,0:km)=var(0:im,0,0:km)
      endif
      if(mpiup .ne. MPI_PROC_NULL) then
        ! pack the down send buffer
        sbuf2(0:im,0:km)=var(0:im,jm,0:km)
      end if
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpidown,mpitag,             &
                        rbuf1,ncou,mpi_real8,mpiup,  mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpiup,  mpitag,             &
                        rbuf2,ncou,mpi_real8,mpidown,mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpiup .ne. MPI_PROC_NULL) then
        ! unpack the received the packet from up
        var(0:im,jm,0:km)=0.5d0*(var(0:im,jm,0:km)+rbuf1(0:im,0:km))
        !
      endif
      !
      if(mpidown .ne. MPI_PROC_NULL) then
        ! unpack the received the packet from down
        var(0:im,0,0:km)=0.5d0*(var(0:im,0,0:km)+rbuf2(0:im,0:km))
        !
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
    endif
    !
    if(ksize==1 .and. lkhomo) then
      !
      var(0:im,0:jm,0)    =0.5d0*(var(0:im,0:jm,0)+var(0:im,0:jm,km))
      var(0:im,0:jm,km)   =var(0:im,0:jm,0)
      !             
    else
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(jm+1)
      !
      allocate( sbuf1(0:im,0:jm),                      &
                sbuf2(0:im,0:jm),                      &
                rbuf1(0:im,0:jm),                      &
                rbuf2(0:im,0:jm) )
      !
      if(mpiback .ne. MPI_PROC_NULL) then
        ! pack the back send buffer
        sbuf1(0:im,0:jm)=var(0:im,0:jm,0)
      endif
      if(mpifront .ne. MPI_PROC_NULL) then
        ! pack the front send buffer
        sbuf2(0:im,0:jm)=var(0:im,0:jm,km)
      endif
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpiback, mpitag,          &
                        rbuf1,ncou,mpi_real8,mpifront,mpitag,          &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpifront,mpitag,          &
                        rbuf2,ncou,mpi_real8,mpiback, mpitag,          &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpifront .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet from front
        var(0:im,0:jm,km)=0.5d0*(var(0:im,0:jm,km)+rbuf1(0:im,0:jm))
        !
      end if
      !
      if(mpiback .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet back
        var(0:im,0:jm,0)=0.5d0*(var(0:im,0:jm,0)+rbuf2(0:im,0:jm))
        !
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    if(present(subtime)) subtime=subtime+ptime()-time_beg
    !
    return
    !
  end subroutine array3d_sync
  !+-------------------------------------------------------------------+
  !| The end of the subroutine array3d_sync.                           |
  !+-------------------------------------------------------------------+
  !!
  !!+------------------------------------------------------------------+
  !| This subroutine is used to swap a 4-D tensor.                     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-Feb-2021 | Created by J. Fang @ Warringon.                     |
  !+-------------------------------------------------------------------+
  subroutine array4d_sendrecv(var,direction,debug,subtime)
    !
    ! arguments
    real(8),intent(inout) :: var(-hm:,-hm:,-hm:,1:)
    integer,intent(in),optional :: direction
    logical,intent(in),optional :: debug
    real(8),intent(inout),optional :: subtime
    !
    ! logical data
    integer :: ncou,nx,dir
    integer :: ierr,j,k
    real(8),allocatable,dimension(:,:,:,:) :: sbuf1,sbuf2,rbuf1,rbuf2
    real(8) :: time_beg
    !
    if(present(subtime)) time_beg=ptime()
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! buf1: send buffer
    ! buf2: send buffer
    ! buf1: redevice buffer
    ! buf2: redevice buffer
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(present(direction)) then
      dir=direction
    else
      dir=0
    endif
    !
    nx=size(var,4)
    !
    if(dir==1 .or. dir==0) then
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in i direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(jm+1)*(km+1)*nx*hm
      !
      allocate( sbuf1(1:hm, 0:jm,0:km,1:nx),                        &
                sbuf2(1:hm, 0:jm,0:km,1:nx),                        &
                rbuf1(1:hm, 0:jm,0:km,1:nx),                        &
                rbuf2(-hm:-1,0:jm,0:km,1:nx) )
      !
      if(mpileft .ne. MPI_PROC_NULL) then
        ! pack the left send buffer
        sbuf1(1:hm,0:jm,0:km,:)=var(1:hm,0:jm,0:km,:)
      endif
      if(mpiright .ne. MPI_PROC_NULL) then
        ! pack the right send buffer
        sbuf2(1:hm,0:jm,0:km,:)=var(im-hm:im-1,0:jm,0:km,:)
      endif
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpileft, mpitag,            &
                        rbuf1,ncou,mpi_real8,mpiright,mpitag,            &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpiright,mpitag,            &
                        rbuf2,ncou,mpi_real8,mpileft, mpitag,            &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpiright .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet from right
        var(im+1:im+hm,0:jm,0:km,:)=rbuf1(1:hm,0:jm,0:km,:)
        !
      end if
        !
      if(mpileft .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet from left
        var(-hm:-1,0:jm,0:km,:)=rbuf2(-hm:-1,0:jm,0:km,:)
        !
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in i direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    if(dir==2 .or. dir==0) then
      !
      if(jsize==1 .and. ljhomo) then
        !
        if(ja==0) then
          do j=-hm,hm
            var(0:im,j,0:km,:)=var(0:im,0,0:km,:)
          enddo
        else
          var(0:im,-hm:-1,0:km,:)    =var(0:im,jm-hm:jm-1,0:km,:)
          var(0:im,jm+1:jm+hm,0:km,:)=var(0:im,1:hm,0:km,:)
        endif
        !
      else
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Message pass in j direction.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ncou=(im+1)*(km+1)*nx*hm
        !
        allocate( sbuf1(0:im,1:hm, 0:km,1:nx),                        &
                  sbuf2(0:im,1:hm, 0:km,1:nx),                        &
                  rbuf1(0:im,1:hm, 0:km,1:nx),                        &
                  rbuf2(0:im,-hm:-1,0:km,1:nx) )
        !
        if(mpidown .ne. MPI_PROC_NULL) then
          ! pack the upper send buffer
          sbuf1(0:im,1:hm,0:km,:)=var(0:im,1:hm,0:km,:)
        endif
        if(mpiup .ne. MPI_PROC_NULL) then
          ! pack the down send buffer
          sbuf2(0:im,1:hm,0:km,:)=var(0:im,jm-hm:jm-1,0:km,:)
        end if
        !
        ! Message passing
        call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpidown,mpitag,             &
                          rbuf1,ncou,mpi_real8,mpiup,  mpitag,             &
                                                 mpi_comm_world,status,ierr)
        mpitag=mpitag+1
        call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpiup,  mpitag,             &
                          rbuf2,ncou,mpi_real8,mpidown,mpitag,             &
                                                 mpi_comm_world,status,ierr)
        mpitag=mpitag+1
        !
        if(mpiup .ne. MPI_PROC_NULL) then
          ! unpack the received the packet from up
          var(0:im,jm+1:jm+hm,0:km,:)=rbuf1(0:im,1:hm,0:km,:)
          !
        endif
        !
        if(mpidown .ne. MPI_PROC_NULL) then
          ! unpack the received the packet from down
          var(0:im,-hm:-1,0:km,:)=rbuf2(0:im,-hm:-1,0:km,:) 
          !
        end if
        !
        deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Finish message pass in j direction.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
      endif
      !
    endif

    if(dir==3 .or. dir==0) then
      !
      !
      if(ksize==1 .and. lkhomo) then
        !
        if(ka==0) then
          do k=-hm,hm
            var(0:im,0:jm,k,:)=var(0:im,0:jm,0,:)
          enddo
        else
          var(0:im,0:jm, -hm:-1   ,:)=var(0:im,0:jm,km-hm:km-1,:)
          var(0:im,0:jm,km+1:km+hm,:)=var(0:im,0:jm,    1:hm,  :)
        endif
        !             
      else
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Message pass in k direction.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ncou=(im+1)*(jm+1)*nx*hm
        !
        allocate( sbuf1(0:im,0:jm, 1:hm,1:nx),                      &
                  sbuf2(0:im,0:jm, 1:hm,1:nx),                      &
                  rbuf1(0:im,0:jm, 1:hm,1:nx),                      &
                  rbuf2(0:im,0:jm,-hm:-1,1:nx) )
        !
        if(mpiback .ne. MPI_PROC_NULL) then
          ! pack the back send buffer
          sbuf1(0:im,0:jm,1:hm,:)=var(0:im,0:jm,1:hm,:)
        endif
        if(mpifront .ne. MPI_PROC_NULL) then
          ! pack the front send buffer
          sbuf2(0:im,0:jm,1:hm,:)=var(0:im,0:jm,km-hm:km-1,:)
        endif
        !
        ! Message passing
        call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpiback, mpitag,          &
                          rbuf1,ncou,mpi_real8,mpifront,mpitag,          &
                                               mpi_comm_world,status,ierr)
        mpitag=mpitag+1
        call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpifront,mpitag,          &
                          rbuf2,ncou,mpi_real8,mpiback, mpitag,          &
                                               mpi_comm_world,status,ierr)
        mpitag=mpitag+1
        !
        if(mpifront .ne. MPI_PROC_NULL) then
          !
          ! unpack the received the packet from front
          var(0:im,0:jm,km+1:km+hm,:)=rbuf1(0:im,0:jm,1:hm,:)
          !
          !
        end if
        !
        if(mpiback .ne. MPI_PROC_NULL) then
          !
          ! unpack the received the packet back
          var(0:im,0:jm,-hm:-1,:)=rbuf2(0:im,0:jm,-hm:-1,:)
          !
        end if
        !
        deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Finish message pass in k direction.
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      endif
      !
    endif
    !
    if(present(subtime)) subtime=subtime+ptime()-time_beg
    !
    return
    !
  end subroutine array4d_sendrecv
  !+-------------------------------------------------------------------+
  !| The end of the subroutine array4d_sendrecv.                       |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to swap a 5-D tensor.                     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-Feb-2021 | Created by J. Fang @ Warringon.                     |
  !+-------------------------------------------------------------------+
  subroutine array5d_sendrecv(var,subtime)
    !
    ! arguments
    real(8),intent(inout) :: var(-hm:,-hm:,-hm:,1:,1:)
    real(8),intent(inout),optional :: subtime
    !
    ! logical data
    integer :: ncou,nx,mx
    integer :: ierr,j,k
    real(8),allocatable,dimension(:,:,:,:,:) :: sbuf1,sbuf2,rbuf1,rbuf2
    real(8) :: time_beg
    !
    if(present(subtime)) time_beg=ptime()
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! buf1: send buffer
    ! buf2: send buffer
    ! buf1: redevice buffer
    ! buf2: redevice buffer
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    nx=size(var,4)
    mx=size(var,5)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ncou=(jm+1)*(km+1)*nx*mx*hm
    !
    allocate( sbuf1(1:hm, 0:jm,0:km,1:nx,1:mx),                        &
              sbuf2(1:hm, 0:jm,0:km,1:nx,1:mx),                        &
              rbuf1(1:hm, 0:jm,0:km,1:nx,1:mx),                        &
              rbuf2(-hm:-1,0:jm,0:km,1:nx,1:mx) )
    !
    if(mpileft .ne. MPI_PROC_NULL) then
      ! pack the left send buffer
      sbuf1(1:hm,0:jm,0:km,:,:)=var(1:hm,0:jm,0:km,:,:)
    endif
    if(mpiright .ne. MPI_PROC_NULL) then
      ! pack the right send buffer
      sbuf2(1:hm,0:jm,0:km,:,:)=var(im-hm:im-1,0:jm,0:km,:,:)
    endif
    !
    ! Message passing
    call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpileft, mpitag,            &
                      rbuf1,ncou,mpi_real8,mpiright,mpitag,            &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpiright,mpitag,            &
                      rbuf2,ncou,mpi_real8,mpileft, mpitag,            &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    if(mpiright .ne. MPI_PROC_NULL) then
      !
      ! unpack the received the packet from right
      var(im+1:im+hm,0:jm,0:km,:,:)=rbuf1(1:hm,0:jm,0:km,:,:)
      !
    end if
      !
    if(mpileft .ne. MPI_PROC_NULL) then
      !
      ! unpack the received the packet from left
      var(-hm:-1,0:jm,0:km,:,:)=rbuf2(-hm:-1,0:jm,0:km,:,:)
      !
    end if
    !
    deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Finish message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(jsize==1 .and. ljhomo) then
      !
      if(ja==0) then 
        do j=-hm,hm
          var(0:im,j,0:km,:,:)=var(0:im,0,0:km,:,:)
        enddo
      else
        var(0:im, -hm:-1   ,0:km,:,:)=var(0:im,jm-hm:jm-1,0:km,:,:)
        var(0:im,jm+1:jm+hm,0:km,:,:)=var(0:im,    1:hm,  0:km,:,:)
      endif
      !  
    else
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(km+1)*nx*mx*hm
      !
      allocate( sbuf1(0:im,1:hm, 0:km,1:nx,1:mx),                        &
                sbuf2(0:im,1:hm, 0:km,1:nx,1:mx),                        &
                rbuf1(0:im,1:hm, 0:km,1:nx,1:mx),                        &
                rbuf2(0:im,-hm:-1,0:km,1:nx,1:mx) )
      !
      if(mpidown .ne. MPI_PROC_NULL) then
        ! pack the upper send buffer
        sbuf1(0:im,1:hm,0:km,:,:)=var(0:im,1:hm,0:km,:,:)
      endif
      if(mpiup .ne. MPI_PROC_NULL) then
        ! pack the down send buffer
        sbuf2(0:im,1:hm,0:km,:,:)=var(0:im,jm-hm:jm-1,0:km,:,:)
      end if
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpidown,mpitag,             &
                        rbuf1,ncou,mpi_real8,mpiup,  mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpiup,  mpitag,             &
                        rbuf2,ncou,mpi_real8,mpidown,mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpiup .ne. MPI_PROC_NULL) then
        ! unpack the received the packet from up
        var(0:im,jm+1:jm+hm,0:km,:,:)=rbuf1(0:im,1:hm,0:km,:,:)
        !
      endif
      !
      if(mpidown .ne. MPI_PROC_NULL) then
        ! unpack the received the packet from down
        var(0:im,-hm:-1,0:km,:,:)=rbuf2(0:im,-hm:-1,0:km,:,:) 
        !
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    if(ksize==1 .and. lkhomo) then
      !
      if(ka==0) then
        do k=-hm,hm
          var(0:im,0:jm,k,:,:)=var(0:im,0:jm,0,:,:)
        enddo
      else
        var(0:im,0:jm, -hm:-1   ,:,:)=var(0:im,0:jm,km-hm:km-1,:,:)
        var(0:im,0:jm,km+1:km+hm,:,:)=var(0:im,0:jm,    1:hm,  :,:)
      endif
      !             
    else
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(jm+1)*nx*mx*hm
      !
      allocate( sbuf1(0:im,0:jm, 1:hm,1:nx,1:mx),                      &
                sbuf2(0:im,0:jm, 1:hm,1:nx,1:mx),                      &
                rbuf1(0:im,0:jm, 1:hm,1:nx,1:mx),                      &
                rbuf2(0:im,0:jm,-hm:-1,1:nx,1:mx) )
      !
      if(mpiback .ne. MPI_PROC_NULL) then
        ! pack the back send buffer
        sbuf1(0:im,0:jm,1:hm,:,:)=var(0:im,0:jm,1:hm,:,:)
      endif
      if(mpifront .ne. MPI_PROC_NULL) then
        ! pack the front send buffer
        sbuf2(0:im,0:jm,1:hm,:,:)=var(0:im,0:jm,km-hm:km-1,:,:)
      endif
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpiback, mpitag,          &
                        rbuf1,ncou,mpi_real8,mpifront,mpitag,          &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpifront,mpitag,          &
                        rbuf2,ncou,mpi_real8,mpiback, mpitag,          &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpifront .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet from front
        var(0:im,0:jm,km+1:km+hm,:,:)=rbuf1(0:im,0:jm,1:hm,:,:)
        !
      end if
      !
      if(mpiback .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet back
        var(0:im,0:jm,-hm:-1,:,:)=rbuf2(0:im,0:jm,-hm:-1,:,:)
        !
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    if(present(subtime)) subtime=subtime+ptime()-time_beg
    !
    return
    !
  end subroutine array5d_sendrecv
  !+-------------------------------------------------------------------+
  !| The end of the subroutine array5d_sendrecv.                       |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to swap q and update the flow variables.  |
  !| The flow variables at interfaces are also synchronized.           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-Feb-2021 | Created by J. Fang @ Warringon.                     |
  !+-------------------------------------------------------------------+
  subroutine qswap(subtime)
    !
    use commvar,   only: numq,turbmode
    use commarray, only: q,rho,vel,prs,tmp,spc,tke,omg
    use fludyna,   only: q2fvar
    !
    ! argument
    real(8),intent(inout),optional :: subtime
    !
    ! local data
    integer :: ncou
    integer :: ierr,j,k
    real(8),allocatable,dimension(:,:,:,:) :: sbuf1,sbuf2,rbuf1,rbuf2
    real(8) :: time_beg
    !
    if(present(subtime)) time_beg=ptime()
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! buf1: send buffer
    ! buf2: send buffer
    ! buf1: redevice buffer
    ! buf2: redevice buffer
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ncou=(jm+1)*(km+1)*numq*(hm+1)
    !
    allocate( sbuf1(0:hm, 0:jm,0:km,1:numq),                        &
              sbuf2(0:hm, 0:jm,0:km,1:numq),                        &
              rbuf1(0:hm, 0:jm,0:km,1:numq),                        &
              rbuf2(-hm:0,0:jm,0:km,1:numq) )
    !
    if(mpileft .ne. MPI_PROC_NULL) then
      ! pack the left send buffer
      sbuf1(0:hm,0:jm,0:km,:)=q(0:hm,0:jm,0:km,:)
    endif
    if(mpiright .ne. MPI_PROC_NULL) then
      ! pack the right send buffer
      sbuf2(0:hm,0:jm,0:km,:)=q(im-hm:im,0:jm,0:km,:)
    endif
    !
    ! Message passing
    call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpileft, mpitag,            &
                      rbuf1,ncou,mpi_real8,mpiright,mpitag,            &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpiright,mpitag,            &
                      rbuf2,ncou,mpi_real8,mpileft, mpitag,            &
                                             mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    if(mpiright .ne. MPI_PROC_NULL) then
      !
      ! unpack the received the packet from right
      q(im+1:im+hm,0:jm,0:km,:)=rbuf1(1:hm,0:jm,0:km,:)
      !
      q(im,0:jm,0:km,:)=0.5d0*( q(im,0:jm,0:km,:) +            &
                                rbuf1(0,0:jm,0:km,:) )
      !
      if(trim(turbmode)=='k-omega') then
        call q2fvar(q=q(im:im+hm,0:jm,0:km,:),                         &
                                     density=rho(im:im+hm,0:jm,0:km),  &
                                    velocity=vel(im:im+hm,0:jm,0:km,:),&
                                    pressure=prs(im:im+hm,0:jm,0:km),  &
                                 temperature=tmp(im:im+hm,0:jm,0:km),  &
                                     species=spc(im:im+hm,0:jm,0:km,:),&
                                         tke=tke(im:im+hm,0:jm,0:km),  &
                                       omega=omg(im:im+hm,0:jm,0:km) )
      elseif(trim(turbmode)=='none' .or. trim(turbmode)=='udf1') then
        call q2fvar(q=q(im:im+hm,0:jm,0:km,:),                         &
                                     density=rho(im:im+hm,0:jm,0:km),  &
                                    velocity=vel(im:im+hm,0:jm,0:km,:),&
                                    pressure=prs(im:im+hm,0:jm,0:km),  &
                                 temperature=tmp(im:im+hm,0:jm,0:km),  &
                                     species=spc(im:im+hm,0:jm,0:km,:) )
      else
        stop ' !! ERROR 1 turbmode @ qswap'
      endif
      ! 
    end if
    !
    if(mpileft .ne. MPI_PROC_NULL) then
      !
      ! unpack the received the packet from left
      q(-hm:-1,0:jm,0:km,:)=rbuf2(-hm:-1,0:jm,0:km,:)
      !
      q(0,0:jm,0:km,:)=0.5d0*( q(0,0:jm,0:km,:) +              &
                               rbuf2(0,0:jm,0:km,:) )
      !
      if(trim(turbmode)=='k-omega') then
        call q2fvar(q=q(-hm:0,0:jm,0:km,:),                            &
                                       density=rho(-hm:0,0:jm,0:km),   &
                                      velocity=vel(-hm:0,0:jm,0:km,:), &
                                      pressure=prs(-hm:0,0:jm,0:km),   &
                                   temperature=tmp(-hm:0,0:jm,0:km),   &
                                       species=spc(-hm:0,0:jm,0:km,:), &
                                           tke=tke(-hm:0,0:jm,0:km),   &
                                         omega=omg(-hm:0,0:jm,0:km)    )
      elseif(trim(turbmode)=='none' .or. trim(turbmode)=='udf1') then
        call q2fvar(q=q(-hm:0,0:jm,0:km,:),                            &
                                       density=rho(-hm:0,0:jm,0:km),   &
                                      velocity=vel(-hm:0,0:jm,0:km,:), &
                                      pressure=prs(-hm:0,0:jm,0:km),   &
                                   temperature=tmp(-hm:0,0:jm,0:km),   &
                                       species=spc(-hm:0,0:jm,0:km,:)  )
      else
        stop ' !! ERROR 2 turbmode @ qswap'
      endif
      ! 
    end if
    !
    deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Finish message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(jsize==1 .and. ljhomo) then
      !
      if(jm==0) then
        do j=-hm,hm
          q(0:im,j,0:km,:)=q(0:im,0,0:km,:)
        enddo
      else
        q(0:im, -hm:-1   ,0:km,:)=q(0:im,jm-hm:jm-1,0:km,:)
        q(0:im,jm+1:jm+hm,0:km,:)=q(0:im,    1:hm,  0:km,:)
        !
        q(0:im, 0,0:km,:)=0.5d0*(q(0:im,0,0:km,:)+q(0:im,jm,0:km,:))
        q(0:im,jm,0:km,:)=q(0:im,0,0:km,:)
      endif
      !
      if(trim(turbmode)=='k-omega') then
        call q2fvar(q=q(0:im,-hm:0,0:km,:),                            &
                                       density=rho(0:im,-hm:0,0:km),   &
                                      velocity=vel(0:im,-hm:0,0:km,:), &
                                      pressure=prs(0:im,-hm:0,0:km),   &
                                   temperature=tmp(0:im,-hm:0,0:km),   &
                                       species=spc(0:im,-hm:0,0:km,:), &
                                           tke=tke(0:im,-hm:0,0:km),   &
                                         omega=omg(0:im,-hm:0,0:km)    )
        call q2fvar(q=q(0:im,jm:jm+hm,0:km,:),                         &
                                  density=rho(0:im,jm:jm+hm,0:km),     &
                                 velocity=vel(0:im,jm:jm+hm,0:km,:),   &
                                 pressure=prs(0:im,jm:jm+hm,0:km),     &
                              temperature=tmp(0:im,jm:jm+hm,0:km),     &
                                  species=spc(0:im,jm:jm+hm,0:km,:),   &
                                      tke=tke(0:im,jm:jm+hm,0:km),     &
                                    omega=omg(0:im,jm:jm+hm,0:km)      )
      elseif(trim(turbmode)=='none' .or. trim(turbmode)=='udf1') then
        call q2fvar(q=q(0:im,-hm:0,0:km,:),                            &
                                       density=rho(0:im,-hm:0,0:km),   &
                                      velocity=vel(0:im,-hm:0,0:km,:), &
                                      pressure=prs(0:im,-hm:0,0:km),   &
                                   temperature=tmp(0:im,-hm:0,0:km),   &
                                       species=spc(0:im,-hm:0,0:km,:)  )
        call q2fvar(q=q(0:im,jm:jm+hm,0:km,:),                         &
                                  density=rho(0:im,jm:jm+hm,0:km),     &
                                 velocity=vel(0:im,jm:jm+hm,0:km,:),   &
                                 pressure=prs(0:im,jm:jm+hm,0:km),     &
                              temperature=tmp(0:im,jm:jm+hm,0:km),     &
                                  species=spc(0:im,jm:jm+hm,0:km,:)    )
      else
        stop ' !! ERROR 3 turbmode @ qswap'
      endif
      !
    else
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(km+1)*numq*(hm+1)
      !
      allocate( sbuf1(0:im,0:hm, 0:km,1:numq),                        &
                sbuf2(0:im,0:hm, 0:km,1:numq),                        &
                rbuf1(0:im,0:hm, 0:km,1:numq),                        &
                rbuf2(0:im,-hm:0,0:km,1:numq) )
      !
      if(mpidown .ne. MPI_PROC_NULL) then
        ! pack the upper send buffer
        sbuf1(0:im,0:hm,0:km,:)=q(0:im,0:hm,0:km,:)
      endif
      if(mpiup .ne. MPI_PROC_NULL) then
        ! pack the down send buffer
        sbuf2(0:im,0:hm,0:km,:)=q(0:im,jm-hm:jm,0:km,:)
      end if
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpidown,mpitag,             &
                        rbuf1,ncou,mpi_real8,mpiup,  mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpiup,  mpitag,             &
                        rbuf2,ncou,mpi_real8,mpidown,mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpiup .ne. MPI_PROC_NULL) then
        ! unpack the received the packet from up
        q(0:im,jm+1:jm+hm,0:km,:)=rbuf1(0:im,1:hm,0:km,:)
        !
        q(0:im,jm,0:km,:)=0.5d0*( q(0:im,jm,0:km,:) +            &
                                  rbuf1(0:im, 0,0:km,:) )
        !
        if(trim(turbmode)=='k-omega') then
          call q2fvar(q=q(0:im,jm:jm+hm,0:km,:),                           &
                                         density=rho(0:im,jm:jm+hm,0:km),  &
                                        velocity=vel(0:im,jm:jm+hm,0:km,:),&
                                        pressure=prs(0:im,jm:jm+hm,0:km),  &
                                     temperature=tmp(0:im,jm:jm+hm,0:km),  &
                                         species=spc(0:im,jm:jm+hm,0:km,:),&
                                             tke=tke(0:im,jm:jm+hm,0:km),  &
                                           omega=omg(0:im,jm:jm+hm,0:km)   )
        elseif(trim(turbmode)=='none' .or. trim(turbmode)=='udf1') then
          call q2fvar(q=q(0:im,jm:jm+hm,0:km,:),                           &
                                         density=rho(0:im,jm:jm+hm,0:km),  &
                                        velocity=vel(0:im,jm:jm+hm,0:km,:),&
                                        pressure=prs(0:im,jm:jm+hm,0:km),  &
                                     temperature=tmp(0:im,jm:jm+hm,0:km),  &
                                         species=spc(0:im,jm:jm+hm,0:km,:) )
        else
          stop ' !! ERROR 4 turbmode @ qswap'
        endif
      endif
      !
      if(mpidown .ne. MPI_PROC_NULL) then
        ! unpack the received the packet from down
        q(0:im,-hm:-1,0:km,:)=rbuf2(0:im,-hm:-1,0:km,:) 
        !
        q(0:im,0,0:km,:)=0.5d0*( q(0:im, 0,0:km,:) +            &
                                 rbuf2(0:im, 0,0:km,:) )
        !
        if(trim(turbmode)=='k-omega') then
          call q2fvar(q=q(0:im,-hm:0,0:km,:),                              &
                                         density=rho(0:im,-hm:0,0:km),     &
                                        velocity=vel(0:im,-hm:0,0:km,:),   &
                                        pressure=prs(0:im,-hm:0,0:km),     &
                                     temperature=tmp(0:im,-hm:0,0:km),     &
                                         species=spc(0:im,-hm:0,0:km,:),   &
                                             tke=tke(0:im,-hm:0,0:km),     &
                                           omega=omg(0:im,-hm:0,0:km)      )
        elseif(trim(turbmode)=='none' .or. trim(turbmode)=='udf1') then
          call q2fvar(q=q(0:im,-hm:0,0:km,:),                              &
                                         density=rho(0:im,-hm:0,0:km),     &
                                        velocity=vel(0:im,-hm:0,0:km,:),   &
                                        pressure=prs(0:im,-hm:0,0:km),     &
                                     temperature=tmp(0:im,-hm:0,0:km),     &
                                         species=spc(0:im,-hm:0,0:km,:)    )
        else
          stop ' !! ERROR 5 turbmode @ qswap'
        endif
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    if(ksize==1 .and. lkhomo) then
      !
      if(ka==0) then
        do k=-hm,hm
          q(0:im,0:jm,k,:)=q(0:im,0:jm,0,:)
        enddo
      else
        q(0:im,0:jm, -hm:-1   ,:)=q(0:im,0:jm,km-hm:km-1,:)
        q(0:im,0:jm,km+1:km+hm,:)=q(0:im,0:jm,    1:hm,  :)
        !
        q(0:im,0:jm,0,:)=0.5d0*(q(0:im,0:jm,0,:)+q(0:im,0:jm,km,:))
        q(0:im,0:jm,km,:)=q(0:im,0:jm,0,:)
      endif
      !
      if(trim(turbmode)=='k-omega') then
        call q2fvar(q=q(0:im,0:jm,-hm:0,:),                              &
                                       density=rho(0:im,0:jm,-hm:0),     &
                                      velocity=vel(0:im,0:jm,-hm:0,:),   &
                                      pressure=prs(0:im,0:jm,-hm:0),     &
                                   temperature=tmp(0:im,0:jm,-hm:0),     &
                                       species=spc(0:im,0:jm,-hm:0,:),   &
                                           tke=tke(0:im,0:jm,-hm:0),     &
                                         omega=omg(0:im,0:jm,-hm:0)    )
        call q2fvar(q=q(0:im,0:jm,km:km+hm,:),                           &
                                  density=rho(0:im,0:jm,km:km+hm),       &
                                 velocity=vel(0:im,0:jm,km:km+hm,:),     &
                                 pressure=prs(0:im,0:jm,km:km+hm),       &
                              temperature=tmp(0:im,0:jm,km:km+hm),       &
                                  species=spc(0:im,0:jm,km:km+hm,:),     &
                                      tke=tke(0:im,0:jm,km:km+hm),       &
                                    omega=omg(0:im,0:jm,km:km+hm)      )
      elseif(trim(turbmode)=='none' .or. trim(turbmode)=='udf1') then
        call q2fvar(q=q(0:im,0:jm,-hm:0,:),                              &
                                       density=rho(0:im,0:jm,-hm:0),     &
                                      velocity=vel(0:im,0:jm,-hm:0,:),   &
                                      pressure=prs(0:im,0:jm,-hm:0),     &
                                   temperature=tmp(0:im,0:jm,-hm:0),     &
                                       species=spc(0:im,0:jm,-hm:0,:)    )
        call q2fvar(q=q(0:im,0:jm,km:km+hm,:),                           &
                                  density=rho(0:im,0:jm,km:km+hm),       &
                                 velocity=vel(0:im,0:jm,km:km+hm,:),     &
                                 pressure=prs(0:im,0:jm,km:km+hm),       &
                              temperature=tmp(0:im,0:jm,km:km+hm),       &
                                  species=spc(0:im,0:jm,km:km+hm,:)      )
      else
        stop ' !! ERROR 6 turbmode @ qswap'
      endif
      !
    else
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(jm+1)*numq*(hm+1)
      !
      allocate( sbuf1(0:im,0:jm, 0:hm,1:numq),                      &
                sbuf2(0:im,0:jm, 0:hm,1:numq),                      &
                rbuf1(0:im,0:jm, 0:hm,1:numq),                      &
                rbuf2(0:im,0:jm,-hm:0,1:numq) )
      !
      if(mpiback .ne. MPI_PROC_NULL) then
        ! pack the back send buffer
        sbuf1(0:im,0:jm,0:hm,:)=q(0:im,0:jm,0:hm,:)
      endif
      if(mpifront .ne. MPI_PROC_NULL) then
        ! pack the front send buffer
        sbuf2(0:im,0:jm,0:hm,:)=q(0:im,0:jm,km-hm:km,:)
      endif
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpiback, mpitag,          &
                        rbuf1,ncou,mpi_real8,mpifront,mpitag,          &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpifront,mpitag,          &
                        rbuf2,ncou,mpi_real8,mpiback, mpitag,          &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpifront .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet from front
        q(0:im,0:jm,km+1:km+hm,:)=rbuf1(0:im,0:jm,1:hm,:)
        !
        q(0:im,0:jm,km,:)=0.5d0*( q(0:im,0:jm,km,:) +              &
                                  rbuf1(0:im,0:jm, 0,:) )
        !
        if(trim(turbmode)=='k-omega') then
          call q2fvar(q=q(0:im,0:jm,km:km+hm,:),                         &
                                    density=rho(0:im,0:jm,km:km+hm),     &
                                   velocity=vel(0:im,0:jm,km:km+hm,:),   &
                                   pressure=prs(0:im,0:jm,km:km+hm),     &
                                temperature=tmp(0:im,0:jm,km:km+hm),     &
                                    species=spc(0:im,0:jm,km:km+hm,:),   &
                                        tke=tke(0:im,0:jm,km:km+hm),     &
                                      omega=omg(0:im,0:jm,km:km+hm)    )
        elseif(trim(turbmode)=='none' .or. trim(turbmode)=='udf1') then
          call q2fvar(q=q(0:im,0:jm,km:km+hm,:),                         &
                                    density=rho(0:im,0:jm,km:km+hm),     &
                                   velocity=vel(0:im,0:jm,km:km+hm,:),   &
                                   pressure=prs(0:im,0:jm,km:km+hm),     &
                                temperature=tmp(0:im,0:jm,km:km+hm),     &
                                    species=spc(0:im,0:jm,km:km+hm,:)    )
        else
          stop ' !! ERROR 7 turbmode @ qswap'
        endif
        !
      end if
      !
      if(mpiback .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet back
        q(0:im,0:jm,-hm:-1,:)=rbuf2(0:im,0:jm,-hm:-1,:)
        !
        q(0:im,0:jm,0,:)=0.5d0*( q(0:im,0:jm,0,:) +                &
                                 rbuf2(0:im,0:jm,0,:)  )
        !
        if(trim(turbmode)=='k-omega') then
          call q2fvar(q=q(0:im,0:jm,-hm:0,:),                             &
                                         density=rho(0:im,0:jm,-hm:0),    &
                                        velocity=vel(0:im,0:jm,-hm:0,:),  &
                                        pressure=prs(0:im,0:jm,-hm:0),    &
                                     temperature=tmp(0:im,0:jm,-hm:0),    &
                                         species=spc(0:im,0:jm,-hm:0,:),  &
                                             tke=tke(0:im,0:jm,-hm:0),    &
                                           omega=omg(0:im,0:jm,-hm:0)   )
        elseif(trim(turbmode)=='none' .or. trim(turbmode)=='udf1') then
          call q2fvar(q=q(0:im,0:jm,-hm:0,:),                             &
                                         density=rho(0:im,0:jm,-hm:0),    &
                                        velocity=vel(0:im,0:jm,-hm:0,:),  &
                                        pressure=prs(0:im,0:jm,-hm:0),    &
                                     temperature=tmp(0:im,0:jm,-hm:0),    &
                                         species=spc(0:im,0:jm,-hm:0,:)   )
        else
          stop ' !! ERROR 7 turbmode @ qswap'
        endif
        !
        !
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    if(mpitag>10000) mpitag=100
    ! reset mpitag
    !
    if(present(subtime)) subtime=subtime+ptime()-time_beg
    !
    return
    !
  end subroutine qswap
  !+-------------------------------------------------------------------+
  !| The end of the subroutine qswap.                                  |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| The wraper of MPI_Wtime                                           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-November-2019: Created by J. Fang @ STFC Daresbury Laboratory  |
  !+-------------------------------------------------------------------+
  real(8) function ptime()
    !
    ptime=MPI_Wtime()
    !
    return
    !
  end function ptime
  !+-------------------------------------------------------------------+
  !| The end of the function ptime.                                    |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to synconize immbnod(jb)%inmg             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-Jul-2021: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine syncinmg(asboun)
    !
    use commtype,  only : sboun
    !
    ! arguments
    type(sboun),intent(inout) :: asboun
    !
    ! local data
    integer :: ierr,jrank
    real(8) :: vta(0:mpirankmax),vtmin
    integer :: ijk(3,0:mpirankmax)
    !
    call mpi_allgather(asboun%dis_imga_inmg,1,mpi_real8,              &
                       vta,1,mpi_real8,mpi_comm_world,ierr)
    call mpi_allgather(asboun%inmg,3,mpi_integer,              &
                       ijk,3,mpi_integer,mpi_comm_world,ierr)
    !
    ! print*,asboun%dis_imga_inmg,asboun%inmg
    vtmin=1.d15
    do jrank=0,mpirankmax
      !
      if(vta(jrank)<vtmin) then
        vtmin=vta(jrank)
        asboun%dis_imga_inmg=vtmin
        asboun%inmg=ijk(:,jrank)
      endif
      !
    enddo
    !
  end subroutine syncinmg
  !+-------------------------------------------------------------------+
  !| The end of the subroutine syncinmg.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to synconize immbnod(jb)%isup             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-Jul-2021: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine syncisup(asboun,snodes,nsnode)
    !
    use commtype,  only : sboun
    !
    ! arguments
    type(sboun),intent(inout) :: asboun
    integer,intent(in) :: snodes(:,:)
    integer,intent(in) :: nsnode
    !
    ! local data
    integer :: jb
    integer :: ns_tab(0:mpirankmax)
    integer,allocatable :: sb_int(:,:),sbt_int(:,:)
    !
    call ptabupd(var=nsnode,table=ns_tab)
    !
    allocate(sb_int(3,nsnode))
    !
    do jb=1,nsnode
      sb_int(:,jb)=snodes(jb,:)
    enddo
    !
    call pgather(sb_int,sbt_int)
    !
    call rmdup(sbt_int)
    !
    if(size(sbt_int,2)>0) then
      !
      allocate(asboun%isup(size(sbt_int,2),3))
      !
      do jb=1,size(sbt_int,2)
        asboun%isup(jb,:)=sbt_int(:,jb)
      enddo
      !
    else
      print*,' !! WARNING !! no supporting nodes can be found @ syncisup'
      stop
    endif
    !
    ! if(size(sbt_int,2)>9) then
    !   print*,mpirank,'|',nsnode,size(sbt_int,2)
    !   print*,sbt_int(:,:)
    !   print*,'---------------------------'
    !   print*,sbt_int(:,:)
    !   stop
    ! endif
    !
    !
    ! do jb=1,nsnode
    !   sb_rea(:,jb)=var(jb)%x(:)
    ! enddo
    ! call mpi_allgather(snodes,1,mpi_integer,              &
    !                    vta,   1,mpi_integer, mpi_comm_world,ierr)
    
    ! print*,mpirank,'|',sbt_int
    !
  end subroutine syncisup
  !+-------------------------------------------------------------------+
  !| The end of the subroutine syncisup.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to synconize immbnod(jb)%weig             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-Jul-2021: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine syncweig(asboun)
    !
    use commtype,  only : sboun
    !
    ! arguments
    type(sboun),intent(inout) :: asboun
    !
    ! local data
    integer :: jb,nweig
    real(8),allocatable :: sb_rea(:)
    real(8) :: epslion=1.d-16
    real(8) :: varcon
    !
    nweig=size(asboun%weig)
    !
    allocate(sb_rea(nweig))
    !
    sb_rea=asboun%weig
    !
    sb_rea=psum(sb_rea)
    !
    varcon=0.d0
    do jb=1,nweig
      !
      if(sb_rea(jb)<=epslion) then
        sb_rea=0.d0
        sb_rea(jb)=1.d0
        varcon=1.d0
        exit
      else
        sb_rea(jb)=1.d0/sb_rea(jb)
        varcon=varcon+sb_rea(jb)
      endif
      !
    enddo
    sb_rea=sb_rea/varcon
    !
    asboun%weig=sb_rea
    !
    deallocate(sb_rea)
    !
    return
    !
  end subroutine syncweig
  !+-------------------------------------------------------------------+
  !| The end of the subroutine syncisup.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| this subroutine removes duplicate entries from a sorted array     | 
  !+-------------------------------------------------------------------+
  !| ======                                                            | 
  !| AUTHOR                                                            |
  !| ------                                                            |
  !| 15-Oct-2018:  CREATED by J. Fang @ STFC Daresbury Laboratory      |
  !+-------------------------------------------------------------------+
  subroutine rmdup_i4(ary)
    !
    ! arguments
    integer,allocatable,intent(inout) :: ary(:)
    !
    ! local variable
    integer :: i,nsz,nh,is
    logical,allocatable :: mask(:)
    integer,allocatable :: buffer(:)
    !
    if(.not. allocated(ary)) return
    !
    nsz=size(ary)
    !
    allocate(mask(nsz))
    !
    mask=.true.
    nh=nsz
    do i=nsz,2,-1
      if(any(ary(1:i-1)==ary(i))) then
        mask(i)=.false.
        nh=nh-1
      endif
    enddo
    !
    allocate(buffer(1:nh))
    is=0
    do i=1,nsz
      if(mask(i)) then
        is=is+1
        buffer(is)=ary(i)
      endif
    enddo
    !
    deallocate(ary)
    call move_alloc(buffer,ary)
    !
    deallocate(mask)
    !
    !print*,ary
    !print*,'============='
  end subroutine rmdup_i4
  !
  subroutine rmdup_i4_2d(ary)
    !
    ! arguments
    integer,allocatable,intent(inout) :: ary(:,:)
    !
    ! local variable
    integer :: i,j,nsz2,nsz1,nh,js
    logical,allocatable :: mask(:)
    integer,allocatable :: buffer(:,:)
    !
    if(.not. allocated(ary)) return
    !
    nsz1=size(ary,1)
    nsz2=size(ary,2)
    !
    allocate(mask(nsz2))
    !
    mask=.true.
    nh=nsz2
    do j=nsz2,2,-1
      !
      do i=1,j-1
        !
        if(ary(1,i)==ary(1,j) .and. &
           ary(2,i)==ary(2,j) .and. &
           ary(3,i)==ary(3,j) ) then
          !
          mask(j)=.false.
          nh=nh-1
          !
          exit
          !
        endif
        !
      enddo
      !
    enddo
    !
    allocate(buffer(1:nsz1,1:nh))
    js=0
    do j=1,nsz2
      if(mask(j)) then
        js=js+1
        buffer(:,js)=ary(:,j)
      endif
    enddo
    !
    deallocate(ary)
    call move_alloc(buffer,ary)
    !
    deallocate(mask)
    !
  end subroutine rmdup_i4_2d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine rmdup_i4.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| Thes functions are used to get the or of logical variable of all  |
  !|  ranks                                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-January-2020: Created by J. Fang @ STFC Daresbury Laboratory   |
  !+-------------------------------------------------------------------+
  logical function por_log(var)
    !
    ! arguments
    logical,intent(in) :: var
    !
    ! local data
    integer :: ierr
    !
    call mpi_allreduce(var,por_log,1,mpi_logical,MPI_LOR,              &
                                                    mpi_comm_world,ierr)
    !
    return
    !
  end function por_log
  !
  function por_log_array1d(var) result(vout)
    !
    ! arguments
    logical,intent(in) :: var(:)
    logical :: vout(size(var))
    !
    ! local data
    integer :: ierr
    integer :: nsize
    !
    nsize=size(var)
    !
    call mpi_allreduce(var,vout,nsize,mpi_logical,MPI_LOR,             &
                                                    mpi_comm_world,ierr)
    !
    return
    !
  end function por_log_array1d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine por_log.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| The  function is to determin icell from different processor.      |
  !| the fisrt non-zero icell will be broadcasted to all processors    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-January-2020: Created by J. Fang @ STFC Daresbury Laboratory   |
  !+-------------------------------------------------------------------+
  function pgeticell(icellin) result(icellout)
    !
    ! arguments
    integer,intent(in) :: icellin(3)
    integer :: icellout(3)
    !
    ! local data
    integer :: ierr,jrank
    integer :: recvbuf(1:3,0:mpirankmax)
    !
    call mpi_gather(icellin,3,mpi_integer,                             &
                    recvbuf,3,mpi_integer, 0 , mpi_comm_world,ierr)
    !
    if(mpirank==0) then
      do jrank=0,mpirankmax
        if(recvbuf(1,jrank)>0 .and. recvbuf(2,jrank)>0 ) then
          icellout(:)=recvbuf(:,jrank)
          exit
        endif
      enddo
    endif
    !
    call bcast(icellout)
    !
    return
    !
  end function pgeticell
  !
  subroutine pcollecicell(icellin,abound)
    !
    use commtype, only : sboun
    !
    ! arguments
    integer,intent(in) :: icellin(:,:)
    type(sboun),intent(inout) :: abound(:)
    !
    ! local data
    integer :: ierr,jrank,csize,jb
    integer,allocatable :: recvbuf(:,:,:),icellout(:,:)
    !
    csize=size(icellin,1)
    allocate(recvbuf(1:csize,1:3,0:mpirankmax),icellout(1:csize,1:3))
    !
    call mpi_gather(icellin,csize*3,mpi_integer,                       &
                    recvbuf,csize*3,mpi_integer, 0 , mpi_comm_world,ierr)
    !
    if(mpirank==0) then
      !
      do jb=1,csize
        !
        do jrank=0,mpirankmax
          if(recvbuf(jb,1,jrank)>0) then
            ! the first non-zero element
            icellout(jb,:)=recvbuf(jb,:,jrank)
            exit
          endif
        enddo
        !
        if(jrank==mpirankmax+1) then
          stop ' ERROR 1 @ pcollecicell: all icell is 0'
          print*,' ** jrank=',jrank
          print*,' ** jb=',jb
          print*,' ** recvbuf(jb,1,:)',recvbuf(jb,1,:)
        endif
        !
      enddo
      !
    endif
    !
    call bcast(icellout)
    !
    do jb=1,csize
      abound(jb)%icell=icellout(jb,:)
    enddo
    !
    deallocate(recvbuf,icellout)
    !
    return
    !
  end subroutine pcollecicell
  !+-------------------------------------------------------------------+
  !| The end of the subroutine pcollecicell.                           |
  !+-------------------------------------------------------------------+
  !
  !
end module parallel
!+---------------------------------------------------------------------+
!| The end of the module parallel.                                     |
!+---------------------------------------------------------------------+