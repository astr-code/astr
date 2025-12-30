module pastr_multiblock_type

    use iso_fortran_env, only: wp => real64
    use pastr_constdef
    use pastr_commvar,  only: nhalo,lihomo,ljhomo,lkhomo

    implicit none

    type :: stencil_type
      real(wp),pointer :: a(:)
    end type stencil_type
    
    type :: fdm_type
      integer :: first_node,last_node,dim
      type(stencil_type),allocatable :: s(:)

      contains

      procedure :: init =>fdm_solver_init
      procedure :: cal=>fdm_solver_operator
      procedure :: info=>fdm_solver_print

    end type fdm_type

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
      call ablock%fdi(j,k)%init(ablock%im,ablock%node_state(:,j,k))
    enddo
    enddo

    do k=0,ablock%km
    do i=0,ablock%im
      call ablock%fdj(i,k)%init(ablock%jm,ablock%node_state(i,:,k))
    enddo
    enddo

    do j=0,ablock%jm
    do i=0,ablock%im
      call ablock%fdk(i,j)%init(ablock%km,ablock%node_state(i,j,:))
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

    subroutine fdm_solver_init(asolver,dim,node_state)
      use pastr_commvar, only: stencil_1
      class(fdm_type),target :: asolver
      integer,intent(in) :: dim

      character(len=1),intent(in) :: node_state(-nhalo:dim+nhalo)

      integer :: n
      character(len=9) :: cstencil
      
      asolver%dim=dim
      allocate(asolver%s(0:dim))

      do n=0,dim
        cstencil=node_state(n-3)//node_state(n-2)//node_state(n-1)// &
                             '-'//node_state(n)//'-'//               &
                 node_state(n+1)//node_state(n+2)//node_state(n+3)

        allocate(asolver%s(n)%a(-3:3))
        asolver%s(n)%a=0.d0
        
        if(cstencil(1:9)=='fff-f-fff' .or. &
           cstencil(1:9)=='bff-f-fff' .or. &
           cstencil(1:9)=='fff-f-ffb' .or. &
           cstencil(1:9)=='bff-f-ffb' ) then
          ! sixth-order central scheme
          asolver%s(n)%a(-3) =-num1d60
          asolver%s(n)%a(-2) = 0.15d0
          asolver%s(n)%a(-1) =-0.75d0
          
          asolver%s(n)%a( 1) = 0.75d0
          asolver%s(n)%a( 2) =-0.15d0
          asolver%s(n)%a( 3) = num1d60
        elseif(cstencil(2:8)=='ff-f-ff' .or. &
               cstencil(2:8)=='bf-f-ff' .or. &
               cstencil(2:8)=='ff-f-fb' .or. &
               cstencil(2:8)=='bf-f-fb' ) then
          ! fourth-order central scheme
          asolver%s(n)%a(-2) = num1d12
          asolver%s(n)%a(-1) =-num2d3

          asolver%s(n)%a( 1) = num2d3
          asolver%s(n)%a( 2) =-num1d12
        elseif(cstencil(3:7)=='f-f-f' .or. &
               cstencil(3:7)=='b-f-f' .or. &
               cstencil(3:7)=='f-f-b' .or. &
               cstencil(3:7)=='b-f-b'  ) then
         ! second-order central scheme
          asolver%s(n)%a(-1) =-0.5d0

          asolver%s(n)%a( 1) = 0.5d0
        elseif(cstencil(4:7)=='-f-f' .or. &
               cstencil(4:7)=='-b-f'  ) then
          ! second-order biased scheme, i- direction
          asolver%s(n)%a( 0) =-1.5d0
          asolver%s(n)%a( 1) = 2.d0
          asolver%s(n)%a( 2) =-0.5d0
        elseif(cstencil(3:6)=='f-f-' .or. &
               cstencil(3:6)=='f-b-' ) then
          ! second-order biased scheme, i+ direction
          asolver%s(n)%a(-2) =+0.5d0
          asolver%s(n)%a(-1) =-2.d0
          asolver%s(n)%a( 0) = 1.5d0
        elseif(cstencil(1:9)=='bbb-b-bbb') then
          asolver%s(n)%a(-3) =-num1d60
          asolver%s(n)%a(-2) = 0.15d0
          asolver%s(n)%a(-1) =-0.75d0
          
          asolver%s(n)%a( 1) = 0.75d0
          asolver%s(n)%a( 2) =-0.15d0
          asolver%s(n)%a( 3) = num1d60
        elseif(cstencil(2:8)=='bb-b-bb') then
          asolver%s(n)%a(-2) = num1d12
          asolver%s(n)%a(-1) =-num2d3

          asolver%s(n)%a( 1) = num2d3
          asolver%s(n)%a( 2) =-num1d12
        elseif(cstencil(3:7)=='b-b-b') then
          asolver%s(n)%a(-1) =-0.5d0

          asolver%s(n)%a( 1) = 0.5d0
        elseif(cstencil(4:7)=='-b-b') then
          asolver%s(n)%a( 0) =-1.5d0
          asolver%s(n)%a( 1) = 2.d0
          asolver%s(n)%a( 2) =-0.5d0
        elseif(cstencil(3:6)=='b-b-') then
          asolver%s(n)%a(-2) =+0.5d0
          asolver%s(n)%a(-1) =-2.d0
          asolver%s(n)%a( 0) = 1.5d0
        endif

      enddo

    end subroutine fdm_solver_init

    
  function fdm_solver_operator(asolver,f) result(df)
      class(fdm_type),target :: asolver
      real(wp),intent(in) :: f(-nhalo:asolver%dim+nhalo)
      real(wp) :: df(0:asolver%dim)

      integer :: n,j

      df=0.d0

      do n=0,asolver%dim
        do j=-3,3
          df(n)=df(n)+asolver%s(n)%a(j)*f(n+j)
        enddo
      enddo

    end function fdm_solver_operator

    subroutine fdm_solver_print(asolver)
      class(fdm_type),target :: asolver
      integer :: n

      do n=0,size(asolver%s)-1
        print*,asolver%s(n)%a(:)
      enddo

    end subroutine fdm_solver_print

    subroutine block_define(multi_block,nblocks,pblocks)

      use pastr_commtype

      logical,intent(out) :: multi_block
      integer,intent(out) :: nblocks
      type(tblock),intent(out),allocatable :: pblocks(:)

      integer :: i
      logical :: lfex

      inquire(file='blockdef.txt',exist=lfex)
      if(lfex) then
        multi_block=.true.
      else
        multi_block=.false.
        return
      endif

      open(12,file='blockdef.txt')
      read(12,*)nblocks
      allocate(pblocks(nblocks))
      do i=1,nblocks
        read(12,*)pblocks(i)%ilo,pblocks(i)%ihi, &
                  pblocks(i)%jlo,pblocks(i)%jhi, &
                  pblocks(i)%klo,pblocks(i)%khi
        pblocks(i)%im=pblocks(i)%ihi-pblocks(i)%ilo
        pblocks(i)%jm=pblocks(i)%jhi-pblocks(i)%jlo
        pblocks(i)%km=pblocks(i)%khi-pblocks(i)%klo
      enddo
      close(12)
      print*,' >> blockdef.txt'

      print*,' ** nblocks: ',nblocks

      return

    end subroutine block_define

end module pastr_multiblock_type