module pastr_multiblock_type

    use iso_fortran_env, only: wp => real64
    use pastr_constdef
    use pastr_commvar,  only: nhalo,lihomo,ljhomo,lkhomo
    use pastr_fdm
    use pastr_commtype
    use module_tpatch

    implicit none


  type :: link_type

      logical :: active=.false.
      integer :: id
      integer :: block_id(2),patch_id(2)

      contains

      procedure :: info => print_link

    end type link_type

    type :: patch_type
      
      integer :: id,blk_id
      integer :: imin,imax,jmin,jmax,kmin,kmax
      integer :: i_send_min,i_send_max,j_send_min,j_send_max,k_send_min,k_send_max
      integer :: i_recv_min,i_recv_max,j_recv_min,j_recv_max,k_recv_min,k_recv_max
      integer :: nvar,nvarmax

      character(len=1) :: ntyp ! internal (i) or bc (n)
      character(len=2) :: cdir
      character(len=5) :: btyp ! internal (i) or bc (n)
      
      integer,allocatable :: isplit(:),jsplit(:),ksplit(:)

      type(link_type) :: linker
      logical :: shifted=.false.

      character(len=6),allocatable :: varname(:)

      real(wp) :: xrange(2,3)
      real(wp),allocatable :: x(:,:,:,:),xc(:,:,:,:),var(:,:,:,:)
      real(wp),allocatable :: dsend(:,:,:,:),drecv(:,:,:,:)
      character(len=1),allocatable :: csend(:,:,:,:),crecv(:,:,:,:)

      contains

      procedure :: init => alloc_patch
      procedure :: destruct => patch_destruct
      procedure :: copy_head => patch_copy_head
      procedure :: info => print_patch
      procedure :: shift => patch_perid_shift
      procedure :: alloc_data => patch_ijk_range

    end type patch_type

    type :: block_type

      integer :: id
      integer :: im,jm,km
      integer :: nvar,nhalo,npatch

      character(len=6) :: name

      character(len=6),allocatable :: varname(:)

      character(len=1),allocatable :: node_state(:,:,:)
      real(wp),allocatable :: x(:,:,:,:),var(:,:,:,:)
      real(wp),allocatable :: dxi(:,:,:,:,:)
      real(wp),allocatable :: dvel(:,:,:,:,:)
      real(wp),allocatable :: qcriterion(:,:,:)

      type(patch_type),allocatable :: patch(:)

      type(fdm_type),allocatable :: fdi(:,:),fdj(:,:),fdk(:,:)

      contains

      procedure :: init => alloc_block
      procedure :: patch_init => block_patch_init
      procedure :: patch_destruct => block_patch_destruct
      procedure :: patch_data_give => block_patch_data
      procedure :: solver_init => block_solver_init

    end type block_type


contains
  subroutine alloc_block(ablock)

    class(block_type),target :: ablock
    integer :: n

    nhalo=4

    allocate(ablock%varname(ablock%nvar))
    allocate(ablock%x( -nhalo:ablock%im+nhalo, &
                       -nhalo:ablock%jm+nhalo, &
                       -nhalo:ablock%km+nhalo,3) )
    allocate(ablock%var(-nhalo:ablock%im+nhalo, &
                        -nhalo:ablock%jm+nhalo, &
                        -nhalo:ablock%km+nhalo, ablock%nvar))
    allocate( ablock%node_state(-nhalo:ablock%im+nhalo, &
                                -nhalo:ablock%jm+nhalo, &
                                -nhalo:ablock%km+nhalo) )

    print*,' ** block ',ablock%name,' initiated.'

  end subroutine alloc_block

  subroutine block_solver_init(ablock)

    class(block_type),target :: ablock

    integer :: i,j,k

    allocate(ablock%fdi(0:ablock%jm,0:ablock%km))
    allocate(ablock%fdj(0:ablock%im,0:ablock%km))
    allocate(ablock%fdk(0:ablock%im,0:ablock%jm))

    do k=0,ablock%km
    do j=0,ablock%jm
      call ablock%fdi(j,k)%init(ablock%im,nhalo,ablock%node_state(:,j,k))
    enddo
    enddo

    do k=0,ablock%km
    do i=0,ablock%im
      call ablock%fdj(i,k)%init(ablock%jm,nhalo,ablock%node_state(i,:,k))
    enddo
    enddo

    do j=0,ablock%jm
    do i=0,ablock%im
      call ablock%fdk(i,j)%init(ablock%km,nhalo,ablock%node_state(i,j,:))
    enddo
    enddo

  end subroutine block_solver_init

  subroutine alloc_patch(apatch)

    class(patch_type),target :: apatch

    apatch%nvar = 4

    if(.not.allocated(apatch%varname)) allocate(apatch%varname(apatch%nvar))
    if(.not.allocated(apatch%x)) allocate(apatch%x( apatch%imin:apatch%imax, &
                                                    apatch%jmin:apatch%jmax, &
                                                    apatch%kmin:apatch%kmax,3) )
    if(.not.allocated(apatch%var)) allocate(apatch%var(apatch%imin:apatch%imax, &
                                                       apatch%jmin:apatch%jmax, &
                                                       apatch%kmin:apatch%kmax, apatch%nvar))

    write(*,'(2(A,I0),A)')'  ** patch ',apatch%id,' of block ',apatch%blk_id,' initiated.'

  end subroutine alloc_patch

  subroutine patch_destruct(apatch)

    class(patch_type),target :: apatch

    if(allocated(apatch%isplit))  deallocate(apatch%isplit)
    if(allocated(apatch%jsplit))  deallocate(apatch%jsplit)
    if(allocated(apatch%ksplit))  deallocate(apatch%ksplit)
    if(allocated(apatch%varname)) deallocate(apatch%varname)
    if(allocated(apatch%x))       deallocate(apatch%x)
    if(allocated(apatch%xc))      deallocate(apatch%xc)
    if(allocated(apatch%var))     deallocate(apatch%var)
    if(allocated(apatch%dsend))   deallocate(apatch%dsend)
    if(allocated(apatch%drecv))   deallocate(apatch%drecv)
    if(allocated(apatch%csend))   deallocate(apatch%csend)
    if(allocated(apatch%crecv))   deallocate(apatch%crecv)

  end subroutine patch_destruct

  subroutine patch_ijk_range(apatch)

    class(patch_type),target :: apatch

    apatch%i_send_min=apatch%imin;  apatch%i_send_max=apatch%imax
    apatch%j_send_min=apatch%jmin;  apatch%j_send_max=apatch%jmax
    apatch%k_send_min=apatch%kmin;  apatch%k_send_max=apatch%kmax

    apatch%i_recv_min=apatch%imin;  apatch%i_recv_max=apatch%imax
    apatch%j_recv_min=apatch%jmin;  apatch%j_recv_max=apatch%jmax
    apatch%k_recv_min=apatch%kmin;  apatch%k_recv_max=apatch%kmax

    if(apatch%cdir=='i-') then
      apatch%i_send_min=1
      apatch%i_send_max=nhalo

      apatch%i_recv_min=-nhalo
      apatch%i_recv_max=-1
    elseif(apatch%cdir=='i+') then
      apatch%i_send_min=apatch%imax-nhalo
      apatch%i_send_max=apatch%imax-1

      apatch%i_recv_min=apatch%imax+1
      apatch%i_recv_max=apatch%imax+nhalo
    elseif(apatch%cdir=='j-') then
      apatch%j_send_min=1
      apatch%j_send_max=nhalo

      apatch%j_recv_min=-nhalo
      apatch%j_recv_max=-1
    elseif(apatch%cdir=='j+') then
      apatch%j_send_min=apatch%jmax-nhalo
      apatch%j_send_max=apatch%jmax-1

      apatch%j_recv_min=apatch%jmax+1
      apatch%j_recv_max=apatch%jmax+nhalo
    elseif(apatch%cdir=='k-') then
      apatch%k_send_min=1
      apatch%k_send_max=nhalo

      apatch%k_recv_min=-nhalo
      apatch%k_recv_max=-1
    elseif(apatch%cdir=='k+') then
      apatch%k_send_min=apatch%kmax-nhalo
      apatch%k_send_max=apatch%kmax-1

      apatch%k_recv_min=apatch%kmax+1
      apatch%k_recv_max=apatch%kmax+nhalo
    endif

    apatch%nvarmax=4
    allocate(apatch%dsend(apatch%i_send_min:apatch%i_send_max, &
                          apatch%j_send_min:apatch%j_send_max, &
                          apatch%k_send_min:apatch%k_send_max, 1:apatch%nvarmax))
    allocate(apatch%drecv(apatch%i_recv_min:apatch%i_recv_max, &
                          apatch%j_recv_min:apatch%j_recv_max, &
                          apatch%k_recv_min:apatch%k_recv_max, 1:apatch%nvarmax))

    allocate(apatch%csend(apatch%i_send_min:apatch%i_send_max, &
                          apatch%j_send_min:apatch%j_send_max, &
                          apatch%k_send_min:apatch%k_send_max,1))
    allocate(apatch%crecv(apatch%i_recv_min:apatch%i_recv_max, &
                          apatch%j_recv_min:apatch%j_recv_max, &
                          apatch%k_recv_min:apatch%k_recv_max,1))

    
  end subroutine patch_ijk_range

  function patch_copy_head(apatch) result(des_patch)

    class(patch_type),target :: apatch
    type(patch_type) :: des_patch

    des_patch%id     =apatch%id      
    des_patch%blk_id =apatch%blk_id 

    des_patch%imin   =apatch%imin
    des_patch%imax   =apatch%imax
    des_patch%jmin   =apatch%jmin
    des_patch%jmax   =apatch%jmax
    des_patch%kmin   =apatch%kmin
    des_patch%kmax   =apatch%kmax

    des_patch%nvar   =apatch%nvar

    des_patch%ntyp   =apatch%ntyp
    des_patch%cdir   =apatch%cdir
    des_patch%btyp   =apatch%btyp

  end function patch_copy_head

  subroutine print_patch(apatch)

    class(patch_type),target :: apatch

    write(*,'(2(A,I0),A,A)')'  ** patch id:',apatch%id,     &
                                ' block id:',apatch%blk_id, &
                               ' direction:',apatch%cdir
    write(*,'(6(A,I0),A)')'  ** patch range:(',apatch%imin,'-',apatch%imax, &
                                         ')x(',apatch%jmin,'-',apatch%jmax, &
                                         ')x(',apatch%kmin,'-',apatch%kmax,')'
    ! if(allocated(apatch%isplit)) then
    !   write(*,*)' ** patch to be split at i=',apatch%isplit
    ! endif
    ! if(allocated(apatch%jsplit)) then
    !   write(*,*)' ** patch to be split at j=',apatch%jsplit
    ! endif
    ! if(allocated(apatch%ksplit)) then
    !   write(*,*)' ** patch to be split at k=',apatch%ksplit
    ! endif

  end subroutine print_patch

  subroutine print_link(alink)

    class(link_type),target :: alink

    write(*,'(5(A,I0))')'  ** link',alink%id,              &
                 ' linking block: ',alink%block_id(1),     &
                         ' patch: ',alink%patch_id(1),     &
                     ' <-> block: ',alink%block_id(2),     &
                         ' patch: ',alink%patch_id(2)

  end subroutine print_link

  subroutine block_patch_destruct(ablock)

    class(block_type),target :: ablock
    
    integer :: i

    if( allocated(ablock%patch)) then
      do i=1,ablock%npatch
        if(allocated(ablock%patch(i)%isplit))  deallocate(ablock%patch(i)%isplit)
        if(allocated(ablock%patch(i)%jsplit))  deallocate(ablock%patch(i)%jsplit)
        if(allocated(ablock%patch(i)%ksplit))  deallocate(ablock%patch(i)%ksplit)
        if(allocated(ablock%patch(i)%varname)) deallocate(ablock%patch(i)%varname)
        if(allocated(ablock%patch(i)%x))       deallocate(ablock%patch(i)%x)
        if(allocated(ablock%patch(i)%xc))      deallocate(ablock%patch(i)%xc)
        if(allocated(ablock%patch(i)%var))     deallocate(ablock%patch(i)%var)
        if(allocated(ablock%patch(i)%dsend))   deallocate(ablock%patch(i)%dsend)
        if(allocated(ablock%patch(i)%drecv))   deallocate(ablock%patch(i)%drecv)
        if(allocated(ablock%patch(i)%csend))   deallocate(ablock%patch(i)%csend)
        if(allocated(ablock%patch(i)%crecv))   deallocate(ablock%patch(i)%crecv)
      enddo
      deallocate(ablock%patch)
    endif

  end subroutine block_patch_destruct

  subroutine block_patch_init(ablock)

    class(block_type),target :: ablock

    integer :: n,i,j,k
    integer :: ioff, joff, koff
    integer :: ilo, ihi, jlo, jhi, klo, khi

    ablock%npatch=6

    allocate(ablock%patch(ablock%npatch))
    do n=1,ablock%npatch
      ablock%patch(n)%blk_id=ablock%id
      ablock%patch(n)%id=n
    enddo

    ablock%patch(1)%cdir='i-'
    ablock%patch(1)%imin=0
    ablock%patch(1)%imax=0
    ablock%patch(1)%jmin=0
    ablock%patch(1)%jmax=ablock%jm
    ablock%patch(1)%kmin=0
    ablock%patch(1)%kmax=ablock%km

    ablock%patch(2)%cdir='i+'
    ablock%patch(2)%imin=ablock%im
    ablock%patch(2)%imax=ablock%im
    ablock%patch(2)%jmin=0
    ablock%patch(2)%jmax=ablock%jm
    ablock%patch(2)%kmin=0
    ablock%patch(2)%kmax=ablock%km

    ablock%patch(3)%cdir='j-'
    ablock%patch(3)%imin=0
    ablock%patch(3)%imax=ablock%im
    ablock%patch(3)%jmin=0
    ablock%patch(3)%jmax=0
    ablock%patch(3)%kmin=0
    ablock%patch(3)%kmax=ablock%km

    ablock%patch(4)%cdir='j+'
    ablock%patch(4)%imin=0
    ablock%patch(4)%imax=ablock%im
    ablock%patch(4)%jmin=ablock%jm
    ablock%patch(4)%jmax=ablock%jm
    ablock%patch(4)%kmin=0
    ablock%patch(4)%kmax=ablock%km

    ablock%patch(5)%cdir='k-'
    ablock%patch(5)%imin=0
    ablock%patch(5)%imax=ablock%im
    ablock%patch(5)%jmin=0
    ablock%patch(5)%jmax=ablock%jm
    ablock%patch(5)%kmin=0
    ablock%patch(5)%kmax=0

    ablock%patch(6)%cdir='k+'
    ablock%patch(6)%imin=0
    ablock%patch(6)%imax=ablock%im
    ablock%patch(6)%jmin=0
    ablock%patch(6)%jmax=ablock%jm
    ablock%patch(6)%kmin=ablock%km
    ablock%patch(6)%kmax=ablock%km

    if(lihomo) then
      ablock%patch(1)%btyp='perid'
      ablock%patch(2)%btyp='perid'
    endif
    if(ljhomo) then
      ablock%patch(3)%btyp='perid'
      ablock%patch(4)%btyp='perid'
    endif
    if(lkhomo) then
      ablock%patch(5)%btyp='perid'
      ablock%patch(6)%btyp='perid'
    endif

    do n=1,ablock%npatch

      call ablock%patch(n)%init()
      
      ablock%patch(n)%varname(:)=ablock%varname(:)
      ilo=ablock%patch(n)%imin
      ihi=ablock%patch(n)%imax
      jlo=ablock%patch(n)%jmin
      jhi=ablock%patch(n)%jmax
      klo=ablock%patch(n)%kmin
      khi=ablock%patch(n)%kmax

      ablock%patch(n)%x(ilo:ihi,jlo:jhi,klo:khi,:)  =ablock%x(ilo:ihi,jlo:jhi,klo:khi,:)
      ablock%patch(n)%var(ilo:ihi,jlo:jhi,klo:khi,:)=ablock%var(ilo:ihi,jlo:jhi,klo:khi,:)

      !-----------------------------------------
      ! Determine collapsed direction (only once)
      !-----------------------------------------
      if (ilo == ihi) then
        ioff = 0; joff = 1; koff = 1
      elseif (jlo == jhi) then
        ioff = 1; joff = 0; koff = 1
      elseif (klo == khi) then
        ioff = 1; joff = 1; koff = 0
      else
        print*, "!! error @ block_patch_init",ilo,ihi,jlo,jhi,klo,khi
        stop 1
      endif

      !-----------------------------------------
      ! Compute generic index ranges
      !-----------------------------------------
      ilo = ilo + ioff
      jlo = jlo + joff
      klo = klo + koff

      !-----------------------------------------
      ! Allocate unified xc
      !-----------------------------------------
      allocate(ablock%patch(n)%xc(ilo:ihi, jlo:jhi, klo:khi, 1:3))

      !-----------------------------------------
      !  loop 
      !-----------------------------------------
      do k = klo, khi
        do j = jlo, jhi
          do i = ilo, ihi
            ablock%patch(n)%xc(i,j,k,:) = 0.25d0 * ( ablock%patch(n)%x(i-ioff, j-joff, k-koff, :) &
                                                   + ablock%patch(n)%x(i-ioff, j-joff, k,      :) &
                                                   + ablock%patch(n)%x(i,      j-joff, k-koff, :) &
                                                   + ablock%patch(n)%x(i,      j,      k,      :) )
          end do
        end do
      end do

      ablock%patch(n)%xrange(1,1)=minval(ablock%patch(n)%xc(:,:,:,1))
      ablock%patch(n)%xrange(2,1)=maxval(ablock%patch(n)%xc(:,:,:,1))
      ablock%patch(n)%xrange(1,2)=minval(ablock%patch(n)%xc(:,:,:,2))
      ablock%patch(n)%xrange(2,2)=maxval(ablock%patch(n)%xc(:,:,:,2))
      ablock%patch(n)%xrange(1,3)=minval(ablock%patch(n)%xc(:,:,:,3))
      ablock%patch(n)%xrange(2,3)=maxval(ablock%patch(n)%xc(:,:,:,3))

    enddo

  end subroutine block_patch_init

  subroutine patch_perid_shift(apatch,dir)

    use pastr_commvar, only: lx,ly,lz

    class(patch_type),target :: apatch
    character(len=1),intent(in) :: dir

    if(.not.(apatch%shifted) .and. dir=='+') then

      if(apatch%btyp=='perid') then

        if(apatch%cdir=='i-') then
          apatch%x(:,:,:,1)=apatch%x(:,:,:,1) + lx 
        endif
        if(apatch%cdir=='j-') then
          apatch%x(:,:,:,2)=apatch%x(:,:,:,2) + ly 
        endif
        if(apatch%cdir=='k-') then
          apatch%x(:,:,:,3)=apatch%x(:,:,:,3) + lz 
        endif

        if(apatch%cdir=='i-') then
          apatch%xrange(:,1)=apatch%xrange(:,1) + lx 
        endif
        if(apatch%cdir=='j-') then
          apatch%xrange(:,2)=apatch%xrange(:,2) + ly 
        endif
        if(apatch%cdir=='k-') then
          apatch%xrange(:,3)=apatch%xrange(:,3) + lz 
        endif

        apatch%shifted=.true.

      endif

    endif

    if(apatch%shifted .and. dir=='-') then

      if(apatch%btyp=='perid') then

        if(apatch%cdir=='i-') then
          apatch%x(:,:,:,1)=apatch%x(:,:,:,1) - lx 
        endif
        if(apatch%cdir=='j-') then
          apatch%x(:,:,:,2)=apatch%x(:,:,:,2) - ly 
        endif
        if(apatch%cdir=='k-') then
          apatch%x(:,:,:,3)=apatch%x(:,:,:,3) - lz 
        endif

        if(apatch%cdir=='i-') then
          apatch%xrange(:,1)=apatch%xrange(:,1) - lx 
        endif
        if(apatch%cdir=='j-') then
          apatch%xrange(:,2)=apatch%xrange(:,2) - ly 
        endif
        if(apatch%cdir=='k-') then
          apatch%xrange(:,3)=apatch%xrange(:,3) - lz 
        endif

        apatch%shifted=.false.

      endif

    endif

  end subroutine patch_perid_shift

  subroutine block_patch_data(ablock)

    class(block_type),target :: ablock

    integer :: n,i,j,k
    integer :: ioff, joff, koff
    integer :: ilo, ihi, jlo, jhi, klo, khi

    do n=1,ablock%npatch

      call ablock%patch(n)%destruct()

      call ablock%patch(n)%init()
      
      ablock%patch(n)%varname(:)=ablock%varname(:)
      ilo=ablock%patch(n)%imin
      ihi=ablock%patch(n)%imax
      jlo=ablock%patch(n)%jmin
      jhi=ablock%patch(n)%jmax
      klo=ablock%patch(n)%kmin
      khi=ablock%patch(n)%kmax

      ablock%patch(n)%x(ilo:ihi,jlo:jhi,klo:khi,:)  =ablock%x(ilo:ihi,jlo:jhi,klo:khi,:)
      ablock%patch(n)%var(ilo:ihi,jlo:jhi,klo:khi,:)=ablock%var(ilo:ihi,jlo:jhi,klo:khi,:)

      !-----------------------------------------
      ! Determine collapsed direction (only once)
      !-----------------------------------------
      if (ilo == ihi) then
        ioff = 0; joff = 1; koff = 1
      elseif (jlo == jhi) then
        ioff = 1; joff = 0; koff = 1
      elseif (klo == khi) then
        ioff = 1; joff = 1; koff = 0
      else
        print*, "!! error @ block_patch_init",ilo,ihi,jlo,jhi,klo,khi
        stop 1
      endif

      !-----------------------------------------
      ! Compute generic index ranges
      !-----------------------------------------
      ilo = ilo + ioff
      jlo = jlo + joff
      klo = klo + koff

      !-----------------------------------------
      ! Allocate unified xc
      !-----------------------------------------
      allocate(ablock%patch(n)%xc(ilo:ihi, jlo:jhi, klo:khi, 1:3))

      !-----------------------------------------
      !  loop 
      !-----------------------------------------
      do k = klo, khi
        do j = jlo, jhi
          do i = ilo, ihi
            ablock%patch(n)%xc(i,j,k,:) = 0.25d0 * ( ablock%patch(n)%x(i-ioff, j-joff, k-koff, :) &
                                                   + ablock%patch(n)%x(i-ioff, j-joff, k,      :) &
                                                   + ablock%patch(n)%x(i,      j-joff, k-koff, :) &
                                                   + ablock%patch(n)%x(i,      j,      k,      :) )
          end do
        end do
      end do

      ablock%patch(n)%xrange(1,1)=minval(ablock%patch(n)%xc(:,:,:,1))
      ablock%patch(n)%xrange(2,1)=maxval(ablock%patch(n)%xc(:,:,:,1))
      ablock%patch(n)%xrange(1,2)=minval(ablock%patch(n)%xc(:,:,:,2))
      ablock%patch(n)%xrange(2,2)=maxval(ablock%patch(n)%xc(:,:,:,2))
      ablock%patch(n)%xrange(1,3)=minval(ablock%patch(n)%xc(:,:,:,3))
      ablock%patch(n)%xrange(2,3)=maxval(ablock%patch(n)%xc(:,:,:,3))

    enddo

  end subroutine block_patch_data


    subroutine block_define(multi_block,nblocks,pblocks)

      use pastr_commtype
      use pastr_commvar, only : im,jm,km

      logical,intent(out) :: multi_block
      integer,intent(out) :: nblocks
      type(tblock),intent(out),allocatable :: pblocks(:)

      integer :: i
      logical :: lfex

      inquire(file='blockdef.txt',exist=lfex)

      multi_block=.true.

      if(lfex) then

        open(12,file='blockdef.txt')
        read(12,*)nblocks
        allocate(pblocks(nblocks))
        do i=1,nblocks
          read(12,*)pblocks(i)%ilo,pblocks(i)%ihi, &
                    pblocks(i)%jlo,pblocks(i)%jhi, &
                    pblocks(i)%klo,pblocks(i)%khi
          pblocks(i)%id=i
          pblocks(i)%im=pblocks(i)%ihi-pblocks(i)%ilo
          pblocks(i)%jm=pblocks(i)%jhi-pblocks(i)%jlo
          pblocks(i)%km=pblocks(i)%khi-pblocks(i)%klo
          write(pblocks(i)%name,'(A,I5.5)')'b',i
        enddo
        close(12)
        print*,' >> blockdef.txt'

      else

        nblocks=1
        allocate(pblocks(nblocks))
        i=1
        pblocks(i)%ilo=0
        pblocks(i)%ihi=im
        pblocks(i)%jlo=0
        pblocks(i)%jhi=jm
        pblocks(i)%klo=0
        pblocks(i)%khi=km

        pblocks(i)%id=i

        pblocks(i)%im=pblocks(i)%ihi-pblocks(i)%ilo
        pblocks(i)%jm=pblocks(i)%jhi-pblocks(i)%jlo
        pblocks(i)%km=pblocks(i)%khi-pblocks(i)%klo
        write(pblocks(i)%name,'(A,I5.5)')'b',i

      endif

      print*,' ** nblocks: ',nblocks

      call block_tpatch_init(pblocks)

      call block_connection_build(pblocks)

      call block_patch_ijk(pblocks)

      call block_node_stat_setup(pblocks)

      call tblock_solver_init(pblocks)

      return

    end subroutine block_define

    subroutine insert_unique_string(s, list, found, pos)
       character(len=*), intent(in)              :: s
       character(len=*), allocatable, intent(inout) :: list(:)
       logical,intent(out) :: found
       integer,intent(out) :: pos
    
       integer :: i, n

       found=.false.
       pos=-1
    
       if (.not. allocated(list)) then
          allocate(list(1))
          list(1) = s
          pos=1
          return
       end if
    
       n = size(list)
       do i = 1, n
          if (trim(list(i)) == trim(s)) then
            found=.true.
            pos=i
            return
          endif
       end do
    
       list = [list, s]
       pos  = n+1

    end subroutine insert_unique_string

    subroutine fdm_init(ablock,dim,dir)

      use pastr_fdm
      use pastr_commtype

      type(tblock),intent(inout) :: ablock
      integer,intent(in) :: dim
      character(len=1),intent(in) :: dir

      integer :: count,i,j,k
      character(len=dim+1+nhalo+nhalo) :: stemp
      character(len=dim+1+nhalo+nhalo),allocatable :: store(:)
      character(len=1),allocatable :: ns(:)

      logical :: found
      integer :: pos

      allocate(ns(-nhalo:dim+nhalo))

      count=0

      if(dir=='i') then
        allocate(ablock%fdi_map(0:ablock%jm,0:ablock%km))

        do k=0,ablock%km
        do j=0,ablock%jm
  
          do i=-nhalo,ablock%im+nhalo
            stemp(i+nhalo+1:i+nhalo+1)=ablock%nds(i,j,k)
          enddo
  
          call insert_unique_string(stemp,store,found,pos)
  
          ablock%fdi_map(j,k)=pos
  
          if(.not. found) then
  
            count=count+1
  
          endif
  
        enddo
        enddo

        allocate(ablock%fdi(count))

        do i=1,count
          do j=-nhalo,ablock%im+nhalo
            ns(j:j)=store(i)(j+nhalo+1:j+nhalo+1)
          enddo
          call ablock%fdi(i)%init(ablock%im,nhalo,ns)
        enddo

      elseif(dir=='j') then
        allocate(ablock%fdj_map(0:ablock%im,0:ablock%km))

        do k=0,ablock%km
        do i=0,ablock%im
  
          do j=-nhalo,ablock%jm+nhalo
            stemp(j+nhalo+1:j+nhalo+1)=ablock%nds(i,j,k)
          enddo
  
          call insert_unique_string(stemp,store,found,pos)
  
          ablock%fdj_map(i,k)=pos
  
          if(.not. found) then
  
            count=count+1
  
          endif
        enddo
        enddo

        allocate(ablock%fdj(count))
        do i=1,count
          do j=-nhalo,ablock%jm+nhalo
            ns(j:j)=store(i)(j+nhalo+1:j+nhalo+1)
          enddo
          call ablock%fdj(i)%init(ablock%jm,nhalo,ns)
        enddo
        
      elseif(dir=='k') then
        allocate(ablock%fdk_map(0:ablock%im,0:ablock%jm))
        do j=0,ablock%jm
        do i=0,ablock%im
  
          do k=-nhalo,ablock%km+nhalo
            stemp(k+nhalo+1:k+nhalo+1)=ablock%nds(i,j,k)
          enddo
  
          call insert_unique_string(stemp,store,found,pos)
  
          ablock%fdk_map(i,j)=pos
  
          if(.not. found) then
  
            count=count+1
  
          endif
  
        enddo
        enddo

        allocate(ablock%fdk(count))

        do i=1,count
          do j=-nhalo,ablock%km+nhalo
            ns(j:j)=store(i)(j+nhalo+1:j+nhalo+1)
          enddo
          call ablock%fdk(i)%init(ablock%km,nhalo,ns)
        enddo

      endif

    end subroutine fdm_init

    subroutine tblock_solver_init(blocks)

      use pastr_commtype

      type(tblock),intent(inout),target :: blocks(:)

      integer :: i,j,k,n
      type(tblock), pointer :: b

      do n=1,size(blocks)

        b=>blocks(n)

        call fdm_init(b,b%im,'i')
        call fdm_init(b,b%jm,'j')
        call fdm_init(b,b%km,'k')

      enddo

      ! stop

    end subroutine tblock_solver_init

    subroutine block_tpatch_init(blocks)

      use pastr_commtype

      type(tblock),intent(inout),target :: blocks(:)

      integer :: i
      type(tpatch) :: patmp
      type(tblock), pointer :: b

      do i=1,size(blocks)

        b=>blocks(i)

        if(b%ilo==b%ihi) then
          call block_patch_create(b,'j-')
          call block_patch_create(b,'j+')
          call block_patch_create(b,'k-')
          call block_patch_create(b,'k+')
        elseif(b%jlo==b%jhi) then
          call block_patch_create(b,'i-')
          call block_patch_create(b,'i+')
          call block_patch_create(b,'k-')
          call block_patch_create(b,'k+')
        elseif(b%klo==b%khi) then
          call block_patch_create(b,'i-')
          call block_patch_create(b,'i+')
          call block_patch_create(b,'j-')
          call block_patch_create(b,'j+')
        else
          call block_patch_create(b,'i-')
          call block_patch_create(b,'i+')
          call block_patch_create(b,'j-')
          call block_patch_create(b,'j+')
          call block_patch_create(b,'k-')
          call block_patch_create(b,'k+')
        endif

        ! do j=1,b%npatch
      enddo

    end subroutine block_tpatch_init

    subroutine block_patch_create(ablock,dir)
      
      use pastr_commtype

      type(tblock),intent(inout) :: ablock
      character(len=*),intent(in) :: dir

      type(tpatch) :: patmp

      patmp%dir=dir

      patmp%bct='bou'

      patmp%ilo=ablock%ilo
      patmp%ihi=ablock%ihi
      patmp%jlo=ablock%jlo
      patmp%jhi=ablock%jhi
      patmp%klo=ablock%klo
      patmp%khi=ablock%khi

      select case (dir)
        case ('i-')
          patmp%ilo=ablock%ilo
          patmp%ihi=ablock%ilo
        case ('i+')
          patmp%ilo=ablock%ihi
          patmp%ihi=ablock%ihi
        case ('j-')
          patmp%jlo=ablock%jlo
          patmp%jhi=ablock%jlo
        case ('j+')
          patmp%jlo=ablock%jhi
          patmp%jhi=ablock%jhi
        case ('k-')
          patmp%klo=ablock%klo
          patmp%khi=ablock%klo
        case ('k+')
          patmp%klo=ablock%khi
          patmp%khi=ablock%khi
        case default
           stop ' 1 @ patch_create'
      end select

      call ablock%add_patch(patmp)

    end subroutine block_patch_create

    subroutine patch_match_split(pach_a,pach_b,newp)

      use pastr_commtype
      use pastr_commvar, only : im,jm,km

      type(tpatch),intent(in)  :: pach_a,pach_b
      type(tpatch),allocatable,intent(out) :: newp(:)

      integer :: ilo_A,jlo_A,klo_A,ilo_B,jlo_B,klo_B
      integer :: nnew

      if (face_overlap(pach_a,pach_b)) then

        if(pach_a%dir(1:1)=='i') then
          call split_iface(pach_a,pach_b,newp)
        endif

      endif
      
    end subroutine patch_match_split

    logical function overlap_1d(a0,a1,b0,b1)
      integer, intent(in) :: a0,a1,b0,b1
      overlap_1d = max(a0,b0) < min(a1,b1)
    end function overlap_1d

    logical function face_overlap(p1,p2)

      use pastr_commtype

      type(tpatch), intent(in) :: p1,p2

      face_overlap = .false.

      if(p1%dir(1:1) .ne. p2%dir(1:1)) return
      !not the same direction

      select case(p1%dir(1:1))
      case('i') ! I-face
        if (p1%ilo /= p2%ilo) return
        face_overlap = overlap_1d(p1%jlo,p1%jhi,p2%jlo,p2%jhi) .and. &
                       overlap_1d(p1%klo,p1%khi,p2%klo,p2%khi)
    
      case('j') ! J-face
        if (p1%jlo /= p2%jlo) return
        face_overlap = overlap_1d(p1%ilo,p1%ihi,p2%ilo,p2%ihi) .and. &
                       overlap_1d(p1%klo,p1%khi,p2%klo,p2%khi)
    
      case('k') ! K-face
        if (p1%klo /= p2%klo) return
        face_overlap = overlap_1d(p1%ilo,p1%ihi,p2%ilo,p2%ihi) .and. &
                       overlap_1d(p1%jlo,p1%jhi,p2%jlo,p2%jhi)
      end select

    end function face_overlap

    subroutine split_iface(pA,pB,out)

      use pastr_commtype

      type(tpatch), intent(in)  :: pA,pB
      type(tpatch), allocatable, intent(out) :: out(:)
    
      integer :: j0,j1,k0,k1,j,k

      integer,allocatable :: knode(:),jnode(:)
      integer :: nout,pos

      nout = 0
    
      j0 = max(pA%jlo,pB%jlo)
      j1 = min(pA%jhi,pB%jhi)
      k0 = max(pA%klo,pB%klo)
      k1 = min(pA%khi,pB%khi)

      call insert_element(a=knode, pos=1, value=pA%klo)
      call insert_element(a=knode, pos=2, value=pA%khi)
      call insert_element(a=jnode, pos=1, value=pA%jlo)
      call insert_element(a=jnode, pos=2, value=pA%jhi)
    
      pos=1
      ! bottom strip
      if (pA%klo < k0) then
        pos=pos+1
        call insert_element(a=knode, pos=pos, value=k0)
      endif
    
      ! top strip
      if (k1 < pA%khi) then
        pos=pos+1
        call insert_element(a=knode, pos=pos, value=k1)
      endif
    
      ! left strip
      pos=1
      if (pA%jlo < j0) then
        pos=pos+1
        call insert_element(a=jnode, pos=pos, value=j0)
      endif
    
      ! right strip
      if (j1 < pA%jhi) then
        pos=pos+1
        call insert_element(a=jnode, pos=pos, value=j1)
      endif

      nout=(size(knode)-1)*(size(jnode)-1)
      allocate(out(nout))

      nout=0
      do k=1,size(knode)-1
      do j=1,size(jnode)-1
        nout=nout+1
        out(nout)=pA
        out(nout)%jlo=jnode(j)
        out(nout)%jhi=jnode(j+1)
        out(nout)%klo=knode(k)
        out(nout)%khi=knode(k+1)
      enddo
      enddo

    end subroutine split_iface

    subroutine insert_element(a, pos, value)

      integer, allocatable, intent(inout) :: a(:)
      integer, intent(in) :: pos, value
    
      integer, allocatable :: tmp(:)
      integer :: n
    
      if(.not. allocated(a)) then
        allocate(a(1))
        a(1)=value
      else
        n = size(a)
        allocate(tmp(n+1))
      
        tmp(1:pos-1) = a(1:pos-1)
        tmp(pos)     = value
        tmp(pos+1:)  = a(pos:)
        call move_alloc(tmp, a)
      endif

      return

    end subroutine

    subroutine find_shared_face(A, B)

      use pastr_commtype
      use pastr_commvar, only : im,jm,km

      type(tpatch)  :: A, B
   
      integer :: ilo, ihi, jlo, jhi, klo, khi
      integer :: ilo_A,ihi_A,jlo_A,jhi_A,klo_A,khi_A, &
                 ilo_B,ihi_B,jlo_B,jhi_B,klo_B,khi_B
      integer :: npatch_A

      type(tpatch) :: patmp
      type(tpatch) :: patch
   
      ilo_A=A%ilo
      jlo_A=A%jlo
      klo_A=A%klo

      ilo_B=B%ilo
      jlo_B=B%jlo
      klo_B=B%klo

      if(lihomo .and. ilo_A==0) ilo_A=ilo_A+im
      if(lihomo .and. ilo_B==0) ilo_B=ilo_B+im
      if(ljhomo .and. jlo_A==0) jlo_A=jlo_A+jm
      if(ljhomo .and. jlo_B==0) jlo_B=jlo_B+jm
      if(lkhomo .and. klo_A==0) klo_A=klo_A+km
      if(lkhomo .and. klo_B==0) klo_B=klo_B+km

      ! ---------- i-face ----------
      if( ilo_A == ilo_B .and.                      &
          A%jlo == B%jlo .and. A%jhi == B%jhi .and. &
          A%klo == B%klo .and. A%khi == B%khi ) then

        A%nbr_block=B%blk_id
        A%nbr_patch=B%id
        A%bct='i'
        B%nbr_block=A%blk_id
        B%nbr_patch=A%id
        B%bct='i'
      endif

      ! ---------- j-face ----------
      if( jlo_A == jlo_B .and.                      &
          A%ilo == B%ilo .and. A%ihi == B%ihi .and. &
          A%klo == B%klo .and. A%khi == B%khi ) then

        A%nbr_block=B%blk_id
        A%nbr_patch=B%id
        A%bct='int'
        B%nbr_block=A%blk_id
        B%nbr_patch=A%id
        B%bct='int'
      endif

      ! ---------- k-face ----------
      if( klo_A == klo_B .and.                      &
          A%ilo == B%ilo .and. A%ihi == B%ihi .and. &
          A%jlo == B%jlo .and. A%jhi == B%jhi ) then

        A%nbr_block=B%blk_id
        A%nbr_patch=B%id
        A%bct='int'
        B%nbr_block=A%blk_id
        B%nbr_patch=A%id
        B%bct='int'
      endif

   end subroutine find_shared_face

   pure elemental subroutine reverse_dir(cdir)

     character(len=2), intent(inout) :: cdir
   
     select case (cdir(2:2))
     case ('+')
        cdir(2:2) = '-'
     case ('-')
        cdir(2:2) = '+'
     case default
        error stop 'reverse_dir: invalid direction'
     end select

   end subroutine reverse_dir

   subroutine block_patch_ijk(blocks)

      use pastr_commtype

      type(tblock) :: blocks(:)

      integer :: i,j,nblk

      do i = 1, size(blocks)

        do j=1,blocks(i)%npatches

          blocks(i)%patches(j)%is=blocks(i)%patches(j)%ilo-blocks(i)%ilo
          blocks(i)%patches(j)%ie=blocks(i)%patches(j)%ihi-blocks(i)%ilo
          blocks(i)%patches(j)%js=blocks(i)%patches(j)%jlo-blocks(i)%jlo
          blocks(i)%patches(j)%je=blocks(i)%patches(j)%jhi-blocks(i)%jlo
          blocks(i)%patches(j)%ks=blocks(i)%patches(j)%klo-blocks(i)%klo
          blocks(i)%patches(j)%ke=blocks(i)%patches(j)%khi-blocks(i)%klo
        enddo

      enddo

   end subroutine block_patch_ijk

   subroutine block_connection_build(blocks)

      use pastr_commtype

      type(tblock) :: blocks(:)

      integer :: i,j,m,n
      integer :: np1,np2
      logical :: has_face
      type(tpatch) :: tmp
      integer :: count,npatch,nblk

      type(tpatch),allocatable ::patch2a(:)

      nblk=size(blocks)

      count = 0
      do i = 1, nblk
      do m = 1,blocks(i)%npatches
        do j = 1, nblk
        do n = 1,blocks(j)%npatches

          call patch_match_split(blocks(i)%patches(m),blocks(j)%patches(n),patch2a)

          call patch_replace_insert(blocks(i)%patches,patch2a,m)

        end do
        end do
      end do
      end do

      do i = 1, nblk
        blocks(i)%npatches=size(blocks(i)%patches)
        do m = 1,blocks(i)%npatches
          blocks(i)%patches(m)%id=m
        enddo
      enddo

      do i = 1, nblk
      do m = 1,blocks(i)%npatches
        do j = 1, nblk
        do n = 1,blocks(j)%npatches
          if(i==j .and. m==n) cycle
          call find_shared_face(blocks(i)%patches(m),blocks(j)%patches(n))
        end do
        end do
      end do
      end do

      ! do i = 1, nblk
      ! do m = 1,blocks(i)%npatches
      !   if(blocks(i)%patches(m)%bct=='int') then
      !     print*,blocks(i)%patches(m)%blk_id,blocks(i)%patches(m)%id,blocks(i)%patches(m)%bct,'<-->', &
      !            blocks(blocks(i)%patches(m)%nbr_block)%patches(blocks(i)%patches(m)%nbr_patch)%blk_id, &
      !            blocks(blocks(i)%patches(m)%nbr_block)%patches(blocks(i)%patches(m)%nbr_patch)%id, &
      !            blocks(blocks(i)%patches(m)%nbr_block)%patches(blocks(i)%patches(m)%nbr_patch)%bct
      !   endif
      ! enddo
      ! enddo

   
    end subroutine block_connection_build

    subroutine patch_replace_insert(patches,patchrep,pos)
      use pastr_commtype
      
      type(tpatch),allocatable,intent(inout) :: patches(:)
      type(tpatch),allocatable,intent(inout) :: patchrep(:)
      integer,intent(in) :: pos

      integer :: total_size

      type(tpatch),allocatable :: temp(:)

      if(.not. allocated(patchrep)) return

      if(size(patchrep)==1) then
        ! replace only, actually nothing need to do
      else
        total_size=size(patches)+size(patchrep)-1

        allocate(temp(total_size))

        temp(1:size(patches))=patches(:) ! copy

        temp(pos)=patchrep(1)            ! replace

        temp(size(patches)+1:total_size)=patchrep(2:size(patchrep)) ! insert

        call move_alloc(temp,patches)

      endif

    end subroutine patch_replace_insert

    subroutine block_grid_bc(blocks)

      use pastr_commtype
      use pastr_commvar, only : nhalo

      type(tblock) :: blocks(:)

      integer :: i,j,k,m,n,is,ie,js,je,ks,ke

      do m = 1, size(blocks)
        
        do n=1,blocks(m)%npatches

          if(blocks(m)%patches(n)%bct=='bou') then

            is=blocks(m)%patches(n)%is
            ie=blocks(m)%patches(n)%ie
            js=blocks(m)%patches(n)%js
            je=blocks(m)%patches(n)%je
            ks=blocks(m)%patches(n)%ks
            ke=blocks(m)%patches(n)%ke

            select case (blocks(m)%patches(n)%dir)
               case ('i-')
                 do i=1,nhalo
                   blocks(m)%x(-i,js:je,ks:ke)=2.d0*blocks(m)%x(0,js:je,ks:ke)-blocks(m)%x(i,js:je,ks:ke)
                   blocks(m)%y(-i,js:je,ks:ke)=2.d0*blocks(m)%y(0,js:je,ks:ke)-blocks(m)%y(i,js:je,ks:ke)
                   blocks(m)%z(-i,js:je,ks:ke)=2.d0*blocks(m)%z(0,js:je,ks:ke)-blocks(m)%z(i,js:je,ks:ke)

                   blocks(m)%xt%spi0(-i,js:je,ks:ke)=2.d0*blocks(m)%xt%core(0,js:je,ks:ke)-blocks(m)%xt%core(i,js:je,ks:ke)
                   blocks(m)%yt%spi0(-i,js:je,ks:ke)=2.d0*blocks(m)%yt%core(0,js:je,ks:ke)-blocks(m)%yt%core(i,js:je,ks:ke)
                   blocks(m)%zt%spi0(-i,js:je,ks:ke)=2.d0*blocks(m)%zt%core(0,js:je,ks:ke)-blocks(m)%zt%core(i,js:je,ks:ke)
                 enddo
               case ('i+')
                 do i=1,nhalo
                   blocks(m)%x(blocks(m)%im+i,js:je,ks:ke)=2.d0*blocks(m)%x(blocks(m)%im,js:je,ks:ke)- &
                                                                blocks(m)%x(blocks(m)%im-i,js:je,ks:ke)
                   blocks(m)%y(blocks(m)%im+i,js:je,ks:ke)=2.d0*blocks(m)%y(blocks(m)%im,js:je,ks:ke)- &
                                                                blocks(m)%y(blocks(m)%im-i,js:je,ks:ke)
                   blocks(m)%z(blocks(m)%im+i,js:je,ks:ke)=2.d0*blocks(m)%z(blocks(m)%im,js:je,ks:ke)- &
                                                                blocks(m)%z(blocks(m)%im-i,js:je,ks:ke)

                   blocks(m)%xt%spim(i,js:je,ks:ke)=2.d0*blocks(m)%xt%core(blocks(m)%im,js:je,ks:ke)- &
                                                         blocks(m)%xt%core(blocks(m)%im-i,js:je,ks:ke)
                   blocks(m)%yt%spim(i,js:je,ks:ke)=2.d0*blocks(m)%yt%core(blocks(m)%im,js:je,ks:ke)- &
                                                         blocks(m)%yt%core(blocks(m)%im-i,js:je,ks:ke)
                   blocks(m)%zt%spim(i,js:je,ks:ke)=2.d0*blocks(m)%zt%core(blocks(m)%im,js:je,ks:ke)- &
                                                         blocks(m)%zt%core(blocks(m)%im-i,js:je,ks:ke)
                 enddo
               case ('j-')
                 do j=1,nhalo
                   blocks(m)%x(is:ie,-j,ks:ke)=2.d0*blocks(m)%x(is:ie,0,ks:ke)-blocks(m)%x(is:ie,j,ks:ke)
                   blocks(m)%y(is:ie,-j,ks:ke)=2.d0*blocks(m)%y(is:ie,0,ks:ke)-blocks(m)%y(is:ie,j,ks:ke)
                   blocks(m)%z(is:ie,-j,ks:ke)=2.d0*blocks(m)%z(is:ie,0,ks:ke)-blocks(m)%z(is:ie,j,ks:ke)

                   blocks(m)%xt%spj0(is:ie,-j,ks:ke)=2.d0*blocks(m)%xt%core(is:ie,0,ks:ke)-blocks(m)%xt%core(is:ie,j,ks:ke)
                   blocks(m)%yt%spj0(is:ie,-j,ks:ke)=2.d0*blocks(m)%yt%core(is:ie,0,ks:ke)-blocks(m)%yt%core(is:ie,j,ks:ke)
                   blocks(m)%zt%spj0(is:ie,-j,ks:ke)=2.d0*blocks(m)%zt%core(is:ie,0,ks:ke)-blocks(m)%zt%core(is:ie,j,ks:ke)
                 enddo
               case ('j+')
                 do j=1,nhalo
                   blocks(m)%x(is:ie,blocks(m)%jm+j,ks:ke)=2.d0*blocks(m)%x(is:ie,blocks(m)%jm,ks:ke) - &
                                                                blocks(m)%x(is:ie,blocks(m)%jm-j,ks:ke)
                   blocks(m)%y(is:ie,blocks(m)%jm+j,ks:ke)=2.d0*blocks(m)%y(is:ie,blocks(m)%jm,ks:ke) - &
                                                                blocks(m)%y(is:ie,blocks(m)%jm-j,ks:ke)
                   blocks(m)%z(is:ie,blocks(m)%jm+j,ks:ke)=2.d0*blocks(m)%z(is:ie,blocks(m)%jm,ks:ke) - &
                                                                blocks(m)%z(is:ie,blocks(m)%jm-j,ks:ke)

                   blocks(m)%xt%spjm(is:ie,j,ks:ke)=2.d0*blocks(m)%xt%core(is:ie,blocks(m)%jm,ks:ke) - &
                                                         blocks(m)%xt%core(is:ie,blocks(m)%jm-j,ks:ke)
                   blocks(m)%yt%spjm(is:ie,j,ks:ke)=2.d0*blocks(m)%yt%core(is:ie,blocks(m)%jm,ks:ke) - &
                                                         blocks(m)%yt%core(is:ie,blocks(m)%jm-j,ks:ke)
                   blocks(m)%zt%spjm(is:ie,j,ks:ke)=2.d0*blocks(m)%zt%core(is:ie,blocks(m)%jm,ks:ke) - &
                                                         blocks(m)%zt%core(is:ie,blocks(m)%jm-j,ks:ke)
                 enddo
               case ('k-')
                 do k=1,nhalo
                   blocks(m)%x(is:ie,js:je,-k)=2.d0*blocks(m)%x(is:ie,js:je,0)-blocks(m)%x(is:ie,js:je,k)
                   blocks(m)%y(is:ie,js:je,-k)=2.d0*blocks(m)%y(is:ie,js:je,0)-blocks(m)%y(is:ie,js:je,k)
                   blocks(m)%z(is:ie,js:je,-k)=2.d0*blocks(m)%z(is:ie,js:je,0)-blocks(m)%z(is:ie,js:je,k)

                   blocks(m)%xt%spk0(is:ie,js:je,-k)=2.d0*blocks(m)%xt%core(is:ie,js:je,0)-blocks(m)%xt%core(is:ie,js:je,k)
                   blocks(m)%yt%spk0(is:ie,js:je,-k)=2.d0*blocks(m)%yt%core(is:ie,js:je,0)-blocks(m)%yt%core(is:ie,js:je,k)
                   blocks(m)%zt%spk0(is:ie,js:je,-k)=2.d0*blocks(m)%zt%core(is:ie,js:je,0)-blocks(m)%zt%core(is:ie,js:je,k)
                 enddo
               case ('k+')
                 do k=1,nhalo
                   blocks(m)%x(is:ie,js:je,blocks(m)%km+k)=2.d0*blocks(m)%x(is:ie,js:je,blocks(m)%km) - &
                                                                blocks(m)%x(is:ie,js:je,blocks(m)%km-k)
                   blocks(m)%y(is:ie,js:je,blocks(m)%km+k)=2.d0*blocks(m)%y(is:ie,js:je,blocks(m)%km) - &
                                                                blocks(m)%y(is:ie,js:je,blocks(m)%km-k)
                   blocks(m)%z(is:ie,js:je,blocks(m)%km+k)=2.d0*blocks(m)%z(is:ie,js:je,blocks(m)%km) - &
                                                                blocks(m)%z(is:ie,js:je,blocks(m)%km-k)

                   blocks(m)%xt%spkm(is:ie,js:je,k)=2.d0*blocks(m)%xt%core(is:ie,js:je,blocks(m)%km) - &
                                                         blocks(m)%xt%core(is:ie,js:je,blocks(m)%km-k)
                   blocks(m)%yt%spkm(is:ie,js:je,k)=2.d0*blocks(m)%yt%core(is:ie,js:je,blocks(m)%km) - &
                                                         blocks(m)%yt%core(is:ie,js:je,blocks(m)%km-k)
                   blocks(m)%zt%spkm(is:ie,js:je,k)=2.d0*blocks(m)%zt%core(is:ie,js:je,blocks(m)%km) - &
                                                         blocks(m)%zt%core(is:ie,js:je,blocks(m)%km-k)
                 enddo
               case default
                  stop '1 @ block_grid_bc'
            end select

          endif

        enddo

      enddo

      call block_data_swap(blocks,'xyz')

    end subroutine block_grid_bc

    subroutine block_node_stat_setup(blocks)

      use pastr_commtype
      use pastr_commvar, only : nhalo

      type(tblock) :: blocks(:)

      integer :: i,j,k,n
      integer :: ilo,ihi,jlo,jhi,klo,khi

      do n = 1, size(blocks)

        allocate(blocks(n)%nds(-nhalo:blocks(n)%im+nhalo, &
                               -nhalo:blocks(n)%jm+nhalo, &
                               -nhalo:blocks(n)%km+nhalo) )
        blocks(n)%nds='u' !undefined
        blocks(n)%nds(0:blocks(n)%im, &
                      0:blocks(n)%jm, &
                      0:blocks(n)%km)='f' ! all nodes are fluid by default

        ! set interface nodes back to fluids
        do j=1,blocks(n)%npatches
          if(blocks(n)%patches(j)%bct=='bou') then
            ilo=blocks(n)%patches(j)%is
            jlo=blocks(n)%patches(j)%js
            klo=blocks(n)%patches(j)%ks
            ihi=blocks(n)%patches(j)%ie
            jhi=blocks(n)%patches(j)%je
            khi=blocks(n)%patches(j)%ke
            blocks(n)%nds(ilo:ihi,jlo:jhi,klo:khi) ='b'
          endif
        enddo

      enddo

      call block_data_swap(blocks,'nds')

    end subroutine block_node_stat_setup

    subroutine block_data_swap(blocks,dtype)

      use pastr_commtype
      use pastr_commvar, only : nhalo

      type(tblock) :: blocks(:)
      character(len=*),intent(in) :: dtype

      integer :: i,j,m,n,nswap

      select case(dtype)
      case('nds')
        nswap=1
      case('xyz')
        nswap=3
      case('rswap')
        nswap=blocks(1)%nswp
      end select
      
      do i = 1, size(blocks)
        
        do j=1,blocks(i)%npatches

          if(blocks(i)%patches(j)%bct=='int') then

            call blocks(i)%patches(j)%alloc_swap(nhalo,dtype,nswap)
            call pack_patch_data(blocks(i)%patches(j),blocks(i),nhalo,'pack',dtype)
  
            blocks(i)%patches(j)%data_not_swaped=.true.

          endif

        enddo

      enddo

      do i = 1, size(blocks)
        
        do j=1,blocks(i)%npatches

          if(blocks(i)%patches(j)%bct=='int') then
            m=blocks(i)%patches(j)%nbr_block
            n=blocks(i)%patches(j)%nbr_patch

            if(blocks(i)%patches(j)%data_not_swaped) then
              call patch_data_swap(blocks(i)%patches(j),blocks(m)%patches(n),dtype)
            endif

            blocks(i)%patches(j)%data_not_swaped=.false.
            blocks(m)%patches(n)%data_not_swaped=.false.
          endif

        enddo

      enddo


      do i = 1, size(blocks)
        
        do j=1,blocks(i)%npatches

          if(blocks(i)%patches(j)%bct=='int') then

            call pack_patch_data(blocks(i)%patches(j),blocks(i),nhalo,'unpack',dtype)

          endif

        enddo

      enddo

    end subroutine block_data_swap

    subroutine patch_data_swap(patch_a,patch_b,dtype)
      
      type(tpatch),intent(inout) :: patch_a,patch_b
      character(len=*),intent(in) :: dtype

      character(len=1),allocatable :: ctemp(:,:,:)
      real(wp),allocatable :: rtemp(:,:,:,:)

      if(dtype=='nds') then
        ctemp=patch_a%cswap

        patch_a%cswap=patch_b%cswap
  
        patch_b%cswap=ctemp
  
        deallocate(ctemp)
      elseif(dtype=='xyz' .or. dtype=='rswap') then
        rtemp=patch_a%rswap
  
        patch_a%rswap=patch_b%rswap
  
        patch_b%rswap=rtemp
        deallocate(rtemp)
      else
        stop '1 @ patch_data_swap'
      endif

    end subroutine patch_data_swap

    subroutine block_grid_jacobian(blocks)


      type(tblock),intent(inout),target :: blocks(:)

      integer :: i,j,k,n,imap
      type(tblock), pointer :: b
      real(wp),allocatable :: dx(:,:,:,:,:)
      real(wp) :: var1

      do n=1,size(blocks)

        b=>blocks(n)

        allocate( b%ddi(0:b%im,0:b%jm,0:b%km,1:3), &
                  b%ddj(0:b%im,0:b%jm,0:b%km,1:3), &
                  b%ddk(0:b%im,0:b%jm,0:b%km,1:3) )

        allocate( dx(0:b%im,0:b%jm,0:b%km,1:3,1:3) )

        do k=0,b%km
        do j=0,b%jm
          imap=b%fdi_map(j,k)
          dx(:,j,k,1,1)=b%fdi(imap)%cal(b%x(:,j,k),nhalo)
          dx(:,j,k,2,1)=b%fdi(imap)%cal(b%y(:,j,k),nhalo)
          dx(:,j,k,3,1)=b%fdi(imap)%cal(b%z(:,j,k),nhalo)
        enddo
        enddo

        do k=0,b%km
        do i=0,b%im
          imap=b%fdj_map(i,k)
          dx(i,:,k,1,2)=b%fdj(imap)%cal(b%x(i,:,k),nhalo)
          dx(i,:,k,2,2)=b%fdj(imap)%cal(b%y(i,:,k),nhalo)
          dx(i,:,k,3,2)=b%fdj(imap)%cal(b%z(i,:,k),nhalo)
        enddo
        enddo

        do j=0,b%jm
        do i=0,b%im
          imap=b%fdk_map(i,j)
          dx(i,j,:,1,3)=b%fdk(imap)%cal(b%x(i,j,:),nhalo)
          dx(i,j,:,2,3)=b%fdk(imap)%cal(b%y(i,j,:),nhalo)
          dx(i,j,:,3,3)=b%fdk(imap)%cal(b%z(i,j,:),nhalo)
        enddo
        enddo


        do k=0,b%km
        do j=0,b%jm
        do i=0,b%im
          var1= dx(i,j,k,1,1)*dx(i,j,k,2,2)*dx(i,j,k,3,3)                        &
               +dx(i,j,k,1,2)*dx(i,j,k,2,3)*dx(i,j,k,3,1)                        &
               +dx(i,j,k,1,3)*dx(i,j,k,2,1)*dx(i,j,k,3,2)                        &
               -dx(i,j,k,1,3)*dx(i,j,k,2,2)*dx(i,j,k,3,1)                        &
               -dx(i,j,k,1,2)*dx(i,j,k,2,1)*dx(i,j,k,3,3)                        &
               -dx(i,j,k,1,1)*dx(i,j,k,2,3)*dx(i,j,k,3,2)


          b%ddi(i,j,k,1)=dx(i,j,k,2,2)*dx(i,j,k,3,3)-dx(i,j,k,2,3)*dx(i,j,k,3,2)
          b%ddi(i,j,k,2)=dx(i,j,k,1,3)*dx(i,j,k,3,2)-dx(i,j,k,1,2)*dx(i,j,k,3,3)
          b%ddi(i,j,k,3)=dx(i,j,k,1,2)*dx(i,j,k,2,3)-dx(i,j,k,1,3)*dx(i,j,k,2,2)

          b%ddj(i,j,k,1)=dx(i,j,k,2,3)*dx(i,j,k,3,1)-dx(i,j,k,2,1)*dx(i,j,k,3,3)
          b%ddj(i,j,k,2)=dx(i,j,k,1,1)*dx(i,j,k,3,3)-dx(i,j,k,1,3)*dx(i,j,k,3,1)
          b%ddj(i,j,k,3)=dx(i,j,k,1,3)*dx(i,j,k,2,1)-dx(i,j,k,1,1)*dx(i,j,k,2,3)

          b%ddk(i,j,k,1)=dx(i,j,k,2,1)*dx(i,j,k,3,2)-dx(i,j,k,2,2)*dx(i,j,k,3,1)
          b%ddk(i,j,k,2)=dx(i,j,k,1,2)*dx(i,j,k,3,1)-dx(i,j,k,1,1)*dx(i,j,k,3,2)
          b%ddk(i,j,k,3)=dx(i,j,k,1,1)*dx(i,j,k,2,2)-dx(i,j,k,1,2)*dx(i,j,k,2,1)

          b%ddi(i,j,k,1)=b%ddi(i,j,k,1)/var1
          b%ddi(i,j,k,2)=b%ddi(i,j,k,2)/var1
          b%ddi(i,j,k,3)=b%ddi(i,j,k,3)/var1
          b%ddj(i,j,k,1)=b%ddj(i,j,k,1)/var1
          b%ddj(i,j,k,2)=b%ddj(i,j,k,2)/var1
          b%ddj(i,j,k,3)=b%ddj(i,j,k,3)/var1
          b%ddk(i,j,k,1)=b%ddk(i,j,k,1)/var1
          b%ddk(i,j,k,2)=b%ddk(i,j,k,2)/var1
          b%ddk(i,j,k,3)=b%ddk(i,j,k,3)/var1

        end do
        end do
        end do

        deallocate(dx)

        write(*,'(A,I0,A)')'  ** grid jacobian for block ',n,' calculated.'

      enddo


    end subroutine block_grid_jacobian


    subroutine block_res_cal(blocks)

      type(tblock),intent(inout),target :: blocks(:)

      integer :: i,j,k,n,imap
      type(tblock), pointer :: b

      logical :: block_swap=.false.

      do n=1,size(blocks)

        b=>blocks(n)

        do i=1,b%nres
          
          do j=1,b%nvar
            
            if(trim(b%resname(i))==trim(b%varname(j))) then
              b%res(:,:,:,i)=b%var(:,:,:,j)
            endif

          enddo

          if(trim(b%resname(i))=='omegax') then

            block_swap=.true.
            
            call b%init_rswap(2,nhalo)

            do j=1,b%nvar

              if(trim(b%varname(j))=='u1') then
                b%rswap(0:b%im,0:b%jm,0:b%km,1)=b%var(:,:,:,j)
              endif
              if(trim(b%varname(j))=='u2') then
                b%rswap(0:b%im,0:b%jm,0:b%km,2)=b%var(:,:,:,j)
              endif

            enddo

          endif

        enddo

      enddo

      ! call block_data_swap(blocks,'rswap')

    end subroutine block_res_cal

    subroutine pack_patch_data(apatch,ablock,nhalo,exetype,datatype)
      
      type(tpatch),intent(inout) :: apatch
      type(tblock),intent(inout) :: ablock
      character(len=*),intent(in) :: exetype,datatype
      integer,intent(in) :: nhalo

      integer :: i,j,k

      if(datatype=='nds') then
      
        if(apatch%dir(1:2)=='i-') then
          if(exetype=='pack') then
            apatch%cswap(0:nhalo,apatch%js:apatch%je,apatch%ks:apatch%ke) = &
              ablock%nds(0:nhalo,apatch%js:apatch%je,apatch%ks:apatch%ke)
          elseif(exetype=='unpack') then
              ablock%nds(-nhalo:0,apatch%js:apatch%je,apatch%ks:apatch%ke) = &
            apatch%cswap( 0:nhalo,apatch%js:apatch%je,apatch%ks:apatch%ke)
          endif
        elseif(apatch%dir(1:2)=='i+') then
          if(exetype=='pack') then
            apatch%cswap(              0:nhalo,    apatch%js:apatch%je,apatch%ks:apatch%ke) = &
              ablock%nds(ablock%im-nhalo:ablock%im,apatch%js:apatch%je,apatch%ks:apatch%ke)
          elseif(exetype=='unpack') then
              ablock%nds(ablock%im:ablock%im+nhalo,apatch%js:apatch%je,apatch%ks:apatch%ke) = &
            apatch%cswap(        0:nhalo,          apatch%js:apatch%je,apatch%ks:apatch%ke)
          endif
        elseif(apatch%dir(1:2)=='j-') then
          if(exetype=='pack') then
            apatch%cswap(apatch%is:apatch%ie,0:nhalo,apatch%ks:apatch%ke) = &
              ablock%nds(apatch%is:apatch%ie,0:nhalo,apatch%ks:apatch%ke)
          elseif(exetype=='unpack') then
              ablock%nds(apatch%is:apatch%ie,-nhalo:0,apatch%ks:apatch%ke) = &
            apatch%cswap(apatch%is:apatch%ie,     0:nhalo,apatch%ks:apatch%ke)
          endif
        elseif(apatch%dir(1:2)=='j+') then
          if(exetype=='pack') then
            apatch%cswap(apatch%is:apatch%ie,              0:nhalo,    apatch%ks:apatch%ke) = &
              ablock%nds(apatch%is:apatch%ie,ablock%jm-nhalo:ablock%jm,apatch%ks:apatch%ke)
          elseif(exetype=='unpack') then
              ablock%nds(apatch%is:apatch%ie,ablock%jm:ablock%jm+nhalo,apatch%ks:apatch%ke) = &
            apatch%cswap(apatch%is:apatch%ie,        0:nhalo,    apatch%ks:apatch%ke)
          endif
        elseif(apatch%dir(1:2)=='k-') then
          if(exetype=='pack') then
            apatch%cswap(apatch%is:apatch%ie,apatch%js:apatch%je,0:nhalo) = &
              ablock%nds(apatch%is:apatch%ie,apatch%js:apatch%je,0:nhalo)
          elseif(exetype=='unpack') then
              ablock%nds(apatch%is:apatch%ie,apatch%js:apatch%je,-nhalo:0) = &
            apatch%cswap(apatch%is:apatch%ie,apatch%js:apatch%je,     0:nhalo)
          endif
        elseif(apatch%dir(1:2)=='k+') then

          if(exetype=='pack') then
            apatch%cswap(apatch%is:apatch%ie,apatch%js:apatch%je,              0:nhalo) = &
              ablock%nds(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km-nhalo:ablock%km)
          elseif(exetype=='unpack') then
              ablock%nds(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km:ablock%km+nhalo) = &
            apatch%cswap(apatch%is:apatch%ie,apatch%js:apatch%je,        0:nhalo)
          endif
        endif

      elseif(datatype=='xyz') then

        ! print*,exetype

        call pack_data_slap(ablock%xt,apatch,1,nhalo,datatype,exetype)
        call pack_data_slap(ablock%yt,apatch,2,nhalo,datatype,exetype)
        call pack_data_slap(ablock%zt,apatch,3,nhalo,datatype,exetype)
      
        if(apatch%dir(1:2)=='i-') then
          if(exetype=='pack') then

            do i=0,nhalo
              apatch%rswap(i,apatch%js:apatch%je,apatch%ks:apatch%ke,1) = &
                  ablock%x(i,apatch%js:apatch%je,apatch%ks:apatch%ke)   - &
                  ablock%x(0,apatch%js:apatch%je,apatch%ks:apatch%ke)
              apatch%rswap(i,apatch%js:apatch%je,apatch%ks:apatch%ke,2) = &
                  ablock%y(i,apatch%js:apatch%je,apatch%ks:apatch%ke)   - &
                  ablock%y(0,apatch%js:apatch%je,apatch%ks:apatch%ke)
              apatch%rswap(i,apatch%js:apatch%je,apatch%ks:apatch%ke,3) = &
                  ablock%z(i,apatch%js:apatch%je,apatch%ks:apatch%ke)   - &
                  ablock%z(0,apatch%js:apatch%je,apatch%ks:apatch%ke)
            enddo

          elseif(exetype=='unpack') then

            do i=0,nhalo
                  ablock%x(i-nhalo,apatch%js:apatch%je,apatch%ks:apatch%ke)   = &
              apatch%rswap(      i,apatch%js:apatch%je,apatch%ks:apatch%ke,1) + &
                  ablock%x(      0,apatch%js:apatch%je,apatch%ks:apatch%ke)
                  ablock%y(i-nhalo,apatch%js:apatch%je,apatch%ks:apatch%ke)   = &
              apatch%rswap(      i,apatch%js:apatch%je,apatch%ks:apatch%ke,2) + &
                  ablock%y(      0,apatch%js:apatch%je,apatch%ks:apatch%ke)
                  ablock%z(i-nhalo,apatch%js:apatch%je,apatch%ks:apatch%ke)   = &
              apatch%rswap(      i,apatch%js:apatch%je,apatch%ks:apatch%ke,3) + &
                  ablock%z(      0,apatch%js:apatch%je,apatch%ks:apatch%ke)
            enddo

          endif
        elseif(apatch%dir(1:2)=='i+') then
          if(exetype=='pack') then

            do i=0,nhalo
              apatch%rswap(                i,apatch%js:apatch%je,apatch%ks:apatch%ke,1) = &
                  ablock%x(ablock%im-nhalo+i,apatch%js:apatch%je,apatch%ks:apatch%ke)   - &
                  ablock%x(        ablock%im,apatch%js:apatch%je,apatch%ks:apatch%ke)
              apatch%rswap(                i,apatch%js:apatch%je,apatch%ks:apatch%ke,2) = &
                  ablock%y(ablock%im-nhalo+i,apatch%js:apatch%je,apatch%ks:apatch%ke)   - &
                  ablock%y(        ablock%im,apatch%js:apatch%je,apatch%ks:apatch%ke)
              apatch%rswap(                i,apatch%js:apatch%je,apatch%ks:apatch%ke,3) = &
                  ablock%z(ablock%im-nhalo+i,apatch%js:apatch%je,apatch%ks:apatch%ke)   - &
                  ablock%z(        ablock%im,apatch%js:apatch%je,apatch%ks:apatch%ke)
            enddo

          elseif(exetype=='unpack') then

            do i=0,nhalo
                  ablock%x(ablock%im+i,apatch%js:apatch%je,apatch%ks:apatch%ke)   = &
              apatch%rswap(          i,apatch%js:apatch%je,apatch%ks:apatch%ke,1) + &
                  ablock%x(  ablock%im,apatch%js:apatch%je,apatch%ks:apatch%ke)

                  ablock%y(ablock%im+i,apatch%js:apatch%je,apatch%ks:apatch%ke)   = &
              apatch%rswap(          i,apatch%js:apatch%je,apatch%ks:apatch%ke,2) + &
                  ablock%y(  ablock%im,apatch%js:apatch%je,apatch%ks:apatch%ke)

                  ablock%z(ablock%im+i,apatch%js:apatch%je,apatch%ks:apatch%ke)   = &
              apatch%rswap(          i,apatch%js:apatch%je,apatch%ks:apatch%ke,3) + &
                  ablock%z(  ablock%im,apatch%js:apatch%je,apatch%ks:apatch%ke)
            enddo

          endif
        elseif(apatch%dir(1:2)=='j-') then
          if(exetype=='pack') then

            do j=0,nhalo
              apatch%rswap(apatch%is:apatch%ie,j,apatch%ks:apatch%ke,1) = &
                  ablock%x(apatch%is:apatch%ie,j,apatch%ks:apatch%ke)   - &
                  ablock%x(apatch%is:apatch%ie,0,apatch%ks:apatch%ke) 
              apatch%rswap(apatch%is:apatch%ie,j,apatch%ks:apatch%ke,2) = &
                  ablock%y(apatch%is:apatch%ie,j,apatch%ks:apatch%ke)   - &
                  ablock%y(apatch%is:apatch%ie,0,apatch%ks:apatch%ke)
              apatch%rswap(apatch%is:apatch%ie,j,apatch%ks:apatch%ke,3) = &
                  ablock%z(apatch%is:apatch%ie,j,apatch%ks:apatch%ke)   - &
                  ablock%z(apatch%is:apatch%ie,0,apatch%ks:apatch%ke)
            enddo

          elseif(exetype=='unpack') then

            do j=0,nhalo
                  ablock%x(apatch%is:apatch%ie,j-nhalo,apatch%ks:apatch%ke)   = &
              apatch%rswap(apatch%is:apatch%ie,      j,apatch%ks:apatch%ke,1) + &
                  ablock%x(apatch%is:apatch%ie,      0,apatch%ks:apatch%ke)

                  ablock%y(apatch%is:apatch%ie,j-nhalo,apatch%ks:apatch%ke)   = &
              apatch%rswap(apatch%is:apatch%ie,      j,apatch%ks:apatch%ke,2) + &
                  ablock%y(apatch%is:apatch%ie,      0,apatch%ks:apatch%ke)

                  ablock%z(apatch%is:apatch%ie,j-nhalo,apatch%ks:apatch%ke)   = &
              apatch%rswap(apatch%is:apatch%ie,      j,apatch%ks:apatch%ke,3) + &
                  ablock%z(apatch%is:apatch%ie,      0,apatch%ks:apatch%ke)
            enddo

          endif
        elseif(apatch%dir(1:2)=='j+') then
          if(exetype=='pack') then

            do j=0,nhalo
              apatch%rswap(apatch%is:apatch%ie,                j,apatch%ks:apatch%ke,1) = &
                  ablock%x(apatch%is:apatch%ie,ablock%jm-nhalo+j,apatch%ks:apatch%ke)   - &
                  ablock%x(apatch%is:apatch%ie,        ablock%jm,apatch%ks:apatch%ke) 
              apatch%rswap(apatch%is:apatch%ie,                j,apatch%ks:apatch%ke,2) = &
                  ablock%y(apatch%is:apatch%ie,ablock%jm-nhalo+j,apatch%ks:apatch%ke)   - &
                  ablock%y(apatch%is:apatch%ie,        ablock%jm,apatch%ks:apatch%ke)
              apatch%rswap(apatch%is:apatch%ie,                j,apatch%ks:apatch%ke,3) = &
                  ablock%z(apatch%is:apatch%ie,ablock%jm-nhalo+j,apatch%ks:apatch%ke)   - &
                  ablock%z(apatch%is:apatch%ie,        ablock%jm,apatch%ks:apatch%ke)
            enddo

          elseif(exetype=='unpack') then

            do j=0,nhalo
                  ablock%x(apatch%is:apatch%ie,ablock%jm+j,apatch%ks:apatch%ke)   = &
              apatch%rswap(apatch%is:apatch%ie,          j,apatch%ks:apatch%ke,1) + &
                  ablock%x(apatch%is:apatch%ie,  ablock%jm,apatch%ks:apatch%ke)

                  ablock%y(apatch%is:apatch%ie,ablock%jm+j,apatch%ks:apatch%ke)   = &
              apatch%rswap(apatch%is:apatch%ie,          j,apatch%ks:apatch%ke,2) + &
                  ablock%y(apatch%is:apatch%ie,  ablock%jm,apatch%ks:apatch%ke)

                  ablock%z(apatch%is:apatch%ie,ablock%jm+j,apatch%ks:apatch%ke)   = &
              apatch%rswap(apatch%is:apatch%ie,          j,apatch%ks:apatch%ke,3) + &
                  ablock%z(apatch%is:apatch%ie,  ablock%jm,apatch%ks:apatch%ke)
            enddo

          endif
        elseif(apatch%dir(1:2)=='k-') then
          if(exetype=='pack') then

            do k=0,nhalo
              apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,k,1) = &
                  ablock%x(apatch%is:apatch%ie,apatch%js:apatch%je,k)   - &
                  ablock%x(apatch%is:apatch%ie,apatch%js:apatch%je,0)
              apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,k,2) = &
                  ablock%y(apatch%is:apatch%ie,apatch%js:apatch%je,k)   - &
                  ablock%y(apatch%is:apatch%ie,apatch%js:apatch%je,0)
              apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,k,3) = &
                  ablock%z(apatch%is:apatch%ie,apatch%js:apatch%je,k)   - &
                  ablock%z(apatch%is:apatch%ie,apatch%js:apatch%je,0)
            enddo

          elseif(exetype=='unpack') then

            do k=0,nhalo
                  ablock%x(apatch%is:apatch%ie,apatch%js:apatch%je,k-nhalo)   = &
              apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,    k,1)   + &
                  ablock%x(apatch%is:apatch%ie,apatch%js:apatch%je,0)

                  ablock%y(apatch%is:apatch%ie,apatch%js:apatch%je,k-nhalo)   = &
              apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,    k,2)   + &
                  ablock%y(apatch%is:apatch%ie,apatch%js:apatch%je,0)

                  ablock%z(apatch%is:apatch%ie,apatch%js:apatch%je,k-nhalo)   = &
              apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,    k,3)   + &
                  ablock%z(apatch%is:apatch%ie,apatch%js:apatch%je,0)
            enddo

          endif
        elseif(apatch%dir(1:2)=='k+') then

          if(exetype=='pack') then

            do k=0,nhalo
              apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,k,1) = &
                  ablock%x(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km-nhalo+k)   - &
                  ablock%x(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km)
              apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,k,2) = &
                  ablock%y(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km-nhalo+k)   - &
                  ablock%y(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km)
              apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,k,3) = &
                  ablock%z(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km-nhalo+k)   - &
                  ablock%z(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km)
            enddo

          elseif(exetype=='unpack') then

            do k=0,nhalo
                  ablock%x(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km+k)   = &
              apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,        k,1)   + &
                  ablock%x(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km)

                  ablock%y(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km+k)   = &
              apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,        k,2)   + &
                  ablock%y(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km)

                  ablock%z(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km+k)  = &
              apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,        k,3)  + &
                  ablock%z(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km)
            enddo

          endif

        endif

      elseif(datatype=='rswap') then
      
        if(apatch%dir(1:2)=='i-') then
          if(exetype=='pack') then
            apatch%rswap(1:nhalo,apatch%js:apatch%je,apatch%ks:apatch%ke,:) = &
            ablock%rswap(1:nhalo,apatch%js:apatch%je,apatch%ks:apatch%ke,:)
          elseif(exetype=='unpack') then
            ablock%rswap(-nhalo:-1,   apatch%js:apatch%je,apatch%ks:apatch%ke,:)   = &
            apatch%rswap(     1:nhalo,apatch%js:apatch%je,apatch%ks:apatch%ke,:)
          endif
        elseif(apatch%dir(1:2)=='i+') then
          if(exetype=='pack') then
            apatch%rswap(              1:nhalo,      apatch%js:apatch%je,apatch%ks:apatch%ke,:) = &
            ablock%rswap(ablock%im-nhalo:ablock%im-1,apatch%js:apatch%je,apatch%ks:apatch%ke,:)
          elseif(exetype=='unpack') then
            ablock%rswap(ablock%im+1:ablock%im+nhalo,apatch%js:apatch%je,apatch%ks:apatch%ke,:)   = &
            apatch%rswap(          1:nhalo,          apatch%js:apatch%je,apatch%ks:apatch%ke,:)
          endif
        elseif(apatch%dir(1:2)=='j-') then
          if(exetype=='pack') then
            apatch%rswap(apatch%is:apatch%ie,1:nhalo,apatch%ks:apatch%ke,:) = &
            ablock%rswap(apatch%is:apatch%ie,1:nhalo,apatch%ks:apatch%ke,:)
          elseif(exetype=='unpack') then
            ablock%rswap(apatch%is:apatch%ie,-nhalo:-1,   apatch%ks:apatch%ke,:) = &
            apatch%rswap(apatch%is:apatch%ie,     1:nhalo,apatch%ks:apatch%ke,:)
          endif
        elseif(apatch%dir(1:2)=='j+') then
          if(exetype=='pack') then
            apatch%rswap(apatch%is:apatch%ie,              1:nhalo,      apatch%ks:apatch%ke,:) = &
            ablock%rswap(apatch%is:apatch%ie,ablock%jm-nhalo:ablock%jm-1,apatch%ks:apatch%ke,:)
          elseif(exetype=='unpack') then
            ablock%rswap(apatch%is:apatch%ie,ablock%jm+1:ablock%jm+nhalo,apatch%ks:apatch%ke,:)   = &
            apatch%rswap(apatch%is:apatch%ie,          1:nhalo,          apatch%ks:apatch%ke,:)
          endif
        elseif(apatch%dir(1:2)=='k-') then
        
          if(exetype=='pack') then
            apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,1:nhalo,:) = &
            ablock%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,1:nhalo,:)
          elseif(exetype=='unpack') then
            ablock%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,-nhalo:-1,:)   = &
            apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,  1:nhalo,:)
          endif

        elseif(apatch%dir(1:2)=='k+') then

          if(exetype=='pack') then
            apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,1:nhalo,:) = &
            ablock%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km-nhalo:ablock%km-1,:)
          elseif(exetype=='unpack') then
            ablock%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,ablock%km+1:ablock%km+nhalo,:)   = &
            apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,        1:nhalo,:)
          endif

        endif

      endif

    end subroutine pack_patch_data

end module pastr_multiblock_type