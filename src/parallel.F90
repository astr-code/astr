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
    module procedure bcast_r8_ary
    module procedure bcast_r8_ary2
    module procedure bcast_r8_ary3
  end interface
  !
  integer :: mpirank,mpisize,mpirankmax
  integer :: isize,jsize,ksize,irkm,jrkm,krkm,irk,jrk,krk,ig0,jg0,kg0
  integer :: mpileft,mpiright,mpidown,mpiup,mpifront,mpiback
  character(len=8) :: mpirankname
  logical :: lio
  integer :: status(mpi_status_size)
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
    use commvar, only : ia,ja,ka
    !
    logical :: lallo
    integer :: nfactor(1:30)
    integer :: i,j,k,n,n1,n2,n3,nsize
    integer(8) :: nvar1,nvar2,nvar3
    !
    logical :: l2dcomp
    !
    l2dcomp=.false.
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
    if(mpirank==0) then
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
      n=1
      do i=1,mpisize
        !
        if(mod(mpisize,i)==0) then
            nfactor(n)=mpisize/i
            n=n+1
        end if
        !
      end do
      print*,' ** Number of factor is ',n
      !
      do n1=1,n
      do n2=1,n
      do n3=1,n
        !
        if(l2dcomp) then
          !
          nsize=nfactor(n1)*nfactor(n2)*nfactor(n3)
          if(nfactor(n3) .ne. 1) cycle
          !
        else
          if(nfactor(n3)>1) then
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
          nvar3=ja*ka
          nvar3=nvar3*nfactor(n1)
          nvar1=nvar3
          !
          nvar3=ia*ka
          nvar3=nvar3*nfactor(n2)
          nvar1=nvar1+nvar3
          !
          nvar3=ia*ja
          nvar3=nvar3*nfactor(n3)
          nvar1=nvar1+nvar3
          !
          if(nvar1<nvar2 .and. nfactor(n1)>1 .and. nfactor(n2)>1) then
            !
            nvar2=nvar1
            !
            isize=nfactor(n1)
            jsize=nfactor(n2)
            ksize=nfactor(n3)
            !!
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
        call mpistop
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
    use commvar, only : ia,ja,ka,lihomo,ljhomo,lkhomo
    !
    integer, allocatable, dimension(:,:,:) :: imp,jmp,kmp
    integer :: ierr,nproc
    integer :: n,n1,n2
    integer :: nsi,nsj,nsk
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
      !
      open(13,file='datin/parallel.info',form='formatted')
      write(13,"(3(A9,1x))")'isize','jsize','ksize'
      write(13,"(3(I9,1x))")isize,jsize,ksize
      write(13,"(10(A9,1x))")'Rank','Irk','Jrk','Krk','IM','JM','KM',  &
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
        write(13,"(10(I9,1x))")n,irk,jrk,krk,imp(irk,jrk,krk),       &
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
      close(13)
      print*,' << parallel.info ... done !'
      !
    endif
    !
    call bcast(ia)
    call bcast(ja)
    call bcast(ka)
    !
    call bcast(irkm)
    call bcast(jrkm)
    call bcast(krkm)
    !
    call bcast(isize)
    call bcast(jsize)
    call bcast(ksize)
    !
    call bcast(lihomo)
    call bcast(ljhomo)
    call bcast(lkhomo)
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
    use commvar, only: im,jm,km,is,ie,js,je,ks,ke,lihomo,ljhomo,lkhomo,&
                       npdci,npdcj,npdck
    !
    ! local data
    integer :: n,nn,nsize1,nsize,ierr
    integer,allocatable:: ntemp(:,:),nrank(:,:,:)
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
      open(12,file='datin/parallel.info',form='formatted')
      read(12,"(//)")
      do n=0,mpirankmax
        !
        read(12,*)ntemp(1,n),ntemp(2,n),ntemp(3,n),ntemp(4,n),         &
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
      close(12)
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
      irk=irkg(0)
      jrk=jrkg(0)
      krk=krkg(0)
      im=img(0)
      jm=jmg(0)
      km=kmg(0)
      ig0=ntemp(8,0)
      jg0=ntemp(9,0)
      kg0=ntemp(10,0)
      !
      if(ntemp(11,0)==1) lio=.true.
      !
      do n=1,mpirankmax
        call MPI_SEND(ntemp(1,n),11,MPI_INTEGER,n,99,                  &
                                                    mpi_comm_world,ierr)
      end do
      !
      deallocate( ntemp )
    else
      !
      allocate( nrect(1:11) )
      !
      call MPI_RECV(nrect(1),11,MPI_INTEGER,0,99,                      &
                                             mpi_comm_world,status,ierr)
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
      if(nrect(11)==1) lio=.true.
      !
      deallocate( nrect )
      !
    end if
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
      !
      npdcj=3
    else
      !
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
    ! n=(irkm+1)*(jrkm+1)*(krkm+1)
    ! call MPI_BCAST(nrank,n,mpi_integer,0,mpi_comm_world,ierr)
    !
    if(lio) write(*,'(A,I0,A)')'  ** rank ',mpirank,' is the I/O rank'
    !
    ! print*,mpirank,'|',im,jm,km
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
    call mpi_finalize(ierr)
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
  subroutine bcast_int(vario)
    !
    ! arguments
    integer,intent(inout) :: vario
    !
    ! local data
    integer :: ierr
    !
    call MPI_BCAST(vario,1,mpi_integer,0,mpi_comm_world,ierr)
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
  subroutine bcast_r8(vario)
    !
    ! arguments
    real(8),intent(inout) :: vario
    !
    ! local data
    integer :: ierr
    !
    call MPI_BCAST(vario,1,mpi_real8,0,mpi_comm_world,ierr)
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
  subroutine bcast_r8_ary2(vario)
    !
    ! arguments
    real(8),intent(inout) :: vario(:,:)
    !
    ! local data
    integer :: nsize,ierr
    !
    nsize=size(vario,1)*size(vario,2)
    !
    call MPI_BCAST(vario,nsize,mpi_real8,0,mpi_comm_world,ierr)
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
end module parallel
!+---------------------------------------------------------------------+
!| The end of the module parallel.                                     |
!+---------------------------------------------------------------------+