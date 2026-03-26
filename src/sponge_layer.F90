module sponge_layer

  use constdef
  
  implicit none

  logical :: lsponge,lsponge_loc,lspg_i0,lspg_im,lspg_j0,lspg_jm,      &
                                 lspg_k0,lspg_km
  integer :: spg_i0,spg_im,spg_j0,spg_jm,spg_k0,spg_km
  integer :: spg_i0_beg,spg_i0_end,spg_im_beg,spg_im_end, &
             spg_j0_beg,spg_j0_end,spg_jm_beg,spg_jm_end, &
             spg_k0_beg,spg_k0_end,spg_km_beg,spg_km_end

  real(8),allocatable,dimension(:,:,:) :: sponge_damp_coef_i0,sponge_damp_coef_im, &
                                          sponge_damp_coef_j0,sponge_damp_coef_jm, &
                                          sponge_damp_coef_k0,sponge_damp_coef_km, &
                                          sponge_damp_coef
  contains

  !+-------------------------------------------------------------------+
  !| This subroutine is to initilise sponge layer.                     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 11-03-2021  | Created by J. Fang @ Warrington                     |
  !|             | (have not consider subdomain situation)             |
  !| 22-02-2024  | Applied only for Cartesian mesh, include subdomain  |
  !| 19-09-2024  | modified for body fitted mesh.                      |
  !| 03-03-2024  | leave this as udf, based on the geom rather i,j,k.  |
  !+-------------------------------------------------------------------+
  subroutine spongelayerini
    !
    use commvar,  only : spg_def

    if(spg_def=='layer') then
      call spongelayer_define_ijk
    elseif(spg_def=='circl') then
      call spongelayer_define_circle
    else
      print*,'spg_def',spg_def
      stop ' !! sponge layer not defined !!'
    endif
    !
  end subroutine spongelayerini

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used for spatial filter the conservative variable
  ! for stabilizing the computation.
  ! 2-order Explicit filter is incorporated.
  ! Ref1: Datta V. Gaitonde and Miguel R. Visbal, AIAA JOURNAL Vol.38,
  !      No.11, November 2000. 
  ! Ref2: Xavier Gloerfelt and Philippe Lafon, Computers & Fluids, 2008,
  !       37: 388-401.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spongefilter
    
    use commvar,  only : spg_def

    if(spg_def=='layer') then
      call spongefilter_layer
    elseif(spg_def=='circl') then
      call spongefilter_global
    endif

  end subroutine spongefilter
  !
  subroutine spongefilter_layer
    !
    use commvar,  only : is,ie,js,je,ks,ke,numq
    use commarray,only : x,q
    use parallel, only : dataswap,mpirank,kg0
    !
    integer :: i,j,k,n
    real(8) :: var1
    real(8),allocatable :: qtemp(:,:,:,:)
    !
    ! sponger layer attached at the im end.
    if(lspg_i0) then
      !
      call dataswap(q,direction=1)
      !
      if(spg_i0_beg>=0) then
        !
        print*,mpirank,'|',spg_i0_beg,spg_i0_end
        !
        allocate(qtemp(spg_i0_beg:spg_i0_end,js:je,ks:ke,1:numq))
        !
        do k=ks,ke
        do j=js,je
          !
          do i=spg_i0_beg,spg_i0_end
            !
            var1=sponge_damp_coef_i0(i,j,k)
            !
            do n=1,numq
              qtemp(i,j,k,n)=(1.d0-var1)*q(i,j,k,n)+        &
                             num1d6*var1*(q(i+1,j,k,n)+     &
                                          q(i-1,j,k,n)+     &
                                          q(i,j+1,k,n)+     &
                                          q(i,j-1,k,n)+     &
                                          q(i,j,k+1,n)+     &
                                          q(i,j,k-1,n)      )
            enddo
            !
            ! if(j==0) print*,i,sponge_damp_coef_i0(i,j,k)
            !
          enddo
          !
        enddo
        enddo
        !
        do k=ks,ke
        do j=js,je
          !
          do i=spg_i0_beg,spg_i0_end
            !
            do n=1,numq
              q(i,j,k,n)=qtemp(i,j,k,n)
            enddo
            !
          enddo
          !
        enddo
        enddo
        !
        deallocate(qtemp)
        !
      endif
      !
    endif
    !
    ! sponger layer attached at the im end.
    if(lspg_im) then
      !
      call dataswap(q,direction=1)
      !
      if(spg_im_beg>=0) then
        !
        allocate(qtemp(spg_im_beg:spg_im_end,js:je,ks:ke,1:numq))
        !
        do k=ks,ke
        do j=js,je
          !
          do i=spg_im_beg,spg_im_end
            !
            var1=sponge_damp_coef_im(i,j,k)
            !
            do n=1,numq
              qtemp(i,j,k,n)=(1.d0-var1)*q(i,j,k,n)+        &
                             num1d6*var1*(q(i+1,j,k,n)+     &
                                          q(i-1,j,k,n)+     &
                                          q(i,j+1,k,n)+     &
                                          q(i,j-1,k,n)+     &
                                          q(i,j,k+1,n)+     &
                                          q(i,j,k-1,n)      )
            enddo
            !
          enddo
          !
        enddo
        enddo
        !
        do k=ks,ke
        do j=js,je
          !
          do i=spg_im_beg,spg_im_end
            !
            do n=1,numq
              q(i,j,k,n)=qtemp(i,j,k,n)
            enddo
            !
          enddo
          !
        enddo
        enddo
        !
        deallocate(qtemp)
        !
      endif
      !
    endif
    !
    ! sponger layer attached at the jm end.
    if(lspg_jm) then
      !
      call dataswap(q,direction=2)
      !
      if(spg_jm_beg>=0) then
        !
        allocate(qtemp(is:ie,spg_jm_beg:spg_jm_end,ks:ke,1:numq))
        !
        do k=ks,ke
        do i=is,ie
          !
          do j=spg_jm_beg,spg_jm_end
            !
            var1=sponge_damp_coef_jm(i,j,k)
            !
            do n=1,numq
              qtemp(i,j,k,n)=(1.d0-var1)*q(i,j,k,n)+        &
                             num1d6*var1*(q(i+1,j,k,n)+     &
                                          q(i-1,j,k,n)+     &
                                          q(i,j+1,k,n)+     &
                                          q(i,j-1,k,n)+     &
                                          q(i,j,k+1,n)+     &
                                          q(i,j,k-1,n)      )
            enddo
            !
          enddo
          !
        enddo
        enddo
        !
        do k=ks,ke
        do i=is,ie
          !
          do j=spg_jm_beg,spg_jm_end
            !
            do n=1,numq
              q(i,j,k,n)=qtemp(i,j,k,n)
            enddo
            !
          enddo
          !
        enddo
        enddo
        !
        deallocate(qtemp)
        !
      endif
      !
    endif

    ! sponger layer attached at the k=0 end.
    if(lspg_k0) then
      !
      call dataswap(q,direction=3)
      !
      if(spg_k0_beg>=0) then
        !
        allocate(qtemp(is:ie,js:je,spg_k0_beg:spg_k0_end,1:numq))

        do k=spg_k0_beg,spg_k0_end
        do j=js,je
        do i=is,ie

          var1=sponge_damp_coef_k0(i,j,k)
          do n=1,numq
            qtemp(i,j,k,n)=(1.d0-var1)*q(i,j,k,n)+        &
                           num1d6*var1*(q(i+1,j,k,n)+     &
                                        q(i-1,j,k,n)+     &
                                        q(i,j+1,k,n)+     &
                                        q(i,j-1,k,n)+     &
                                        q(i,j,k+1,n)+     &
                                        q(i,j,k-1,n)      )
          enddo

        enddo
        enddo
        enddo
        !
        do k=spg_k0_beg,spg_k0_end
        do j=js,je
        do i=is,ie
          q(i,j,k,1:numq)=qtemp(i,j,k,1:numq)
        enddo
        enddo
        enddo
        !
        deallocate(qtemp)
        !
      endif
      !
    endif

    ! sponger layer attached at the k=km end.
    if(lspg_km) then
      
      call dataswap(q,direction=3)
      !
      if(spg_km_beg>=0) then
        !
        allocate(qtemp(is:ie,js:je,spg_km_beg:spg_km_end,1:numq))

        do k=spg_km_beg,spg_km_end
        do j=js,je
        do i=is,ie

          var1=sponge_damp_coef_km(i,j,k)
          
          do n=1,numq
            qtemp(i,j,k,n)=(1.d0-var1)*q(i,j,k,n)+        &
                           num1d6*var1*(q(i+1,j,k,n)+     &
                                        q(i-1,j,k,n)+     &
                                        q(i,j+1,k,n)+     &
                                        q(i,j-1,k,n)+     &
                                        q(i,j,k+1,n)+     &
                                        q(i,j,k-1,n)      )
          enddo
          
        enddo
        enddo
        enddo
        !
        do k=spg_km_beg,spg_km_end
        do j=js,je
        do i=is,ie
          q(i,j,k,1:numq)=qtemp(i,j,k,1:numq)
        enddo
        enddo
        enddo
        !
        deallocate(qtemp)
        !
      endif
      !
    endif

  end subroutine spongefilter_layer

  subroutine spongefilter_global

    use commvar,  only : is,ie,js,je,ks,ke,numq
    use commarray,only : x,q
    use parallel, only : dataswap
    !
    integer :: i,j,k,n
    real(8) :: var1
    real(8),allocatable :: qtemp(:,:,:,:)
    !
    ! sponger layer attached at the im end.

    if(lsponge) call dataswap(q)

    if(lsponge_loc) then

      allocate(qtemp(is:ie,js:je,ks:ke,1:numq))
  
      do k=ks,ke
      do j=js,je
      do i=is,ie
        var1=sponge_damp_coef(i,j,k)

        do n=1,numq
          qtemp(i,j,k,n)=(1.d0-var1)*q(i,j,k,n)+        &
                         num1d6*var1*(q(i+1,j,k,n)+     &
                                      q(i-1,j,k,n)+     &
                                      q(i,j+1,k,n)+     &
                                      q(i,j-1,k,n)+     &
                                      q(i,j,k+1,n)+     &
                                      q(i,j,k-1,n)      )
        enddo

      enddo
      enddo
      enddo

      q(is:ie,js:je,ks:ke,1:numq)=qtemp(is:ie,js:je,ks:ke,1:numq)

      deallocate(qtemp)

    endif

  end subroutine spongefilter_global
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine spongefilter.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spongelayer_define_circle

    use commvar,  only : is,ie,js,je,ks,ke,spg_def
    use commarray,only : x
    use parallel, only : pmax,por,mpirank

    ! local data
    real(8),parameter :: dampfac=0.05d0
    integer :: i,j,k
    real(8) :: var1,var2,xc,yc,zc,range_spange,max_dis

    lsponge_loc=.false.

    allocate( sponge_damp_coef(is:ie,js:je,ks:ke) )

    xc=0.d0
    yc=0.d0
    zc=0.d0

    range_spange=0.06d0 !60mm

    max_dis=0.d0

    do k=ks,ke
    do j=js,je
    do i=is,ie

      var1=sqrt((x(i,j,k,1)-xc)**2+(x(i,j,k,2)-yc)**2+(x(i,j,k,3)-zc)**2)

      if (var1>=range_spange) then

        var2=(var1-range_spange)**2

        lsponge_loc=.true.
      else
        var2=0.d0
      endif

      sponge_damp_coef(i,j,k)=var2

      max_dis=max(max_dis,var2)

    enddo
    enddo
    enddo

    max_dis=pmax(max_dis)

    lsponge=por(lsponge_loc)

    if(mpirank==0) then
      !
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                        *** sponge layer ***'
      if(lsponge) then
        write(*,'(32X,A,A)')' sponge layer definition: ',spg_def
        write(*,'(A,F6.3)')'  based on the distance to the center, activation when distance >',range_spange
      endif
      write(*,'(2X,62A)')('-',i=1,62)
      !
    endif

    sponge_damp_coef=sponge_damp_coef/max_dis*dampfac

    ! call tecbin('testout/tecsponge'//mpirankname//'.plt',                &
    !                                   x(is:ie,js:je,ks:ke,1),'x',         &
    !                                   x(is:ie,js:je,ks:ke,2),'y',         &
    !                                   x(is:ie,js:je,ks:ke,3),'z',         &
    !                                   sponge_damp_coef(is:ie,js:je,ks:ke),'sp' )


  end subroutine spongelayer_define_circle

  subroutine spongelayer_define_ijk
    !
    use commvar, only : is,ie,js,je,ks,ke
    use commarray,only: x
    use parallel,only : ig0,jg0,kg0,bcast,psum,irk,jrk,krk,irkm,         &
                        jrkm,mpi_igroup,mpi_jgroup,mpi_kgroup,mpistop,   &
                        mpileft,mpiright,mpidown,mpiup,mpifront,mpiback, &
                        precv,psend,ia,ja,ka,mpirank
    !
    ! local data
    real(8),parameter :: dampfac=0.05d0
    real(8),allocatable :: length_sponge(:,:),length_sponge_b(:,:),    &
                           length_sponge_local(:,:)
    real(8) :: dis,var1,var2,var3
    !
    integer :: i,j,k,iglobal_spg_beg,iglobal_beg,iglobal_end,          &
                     jglobal_spg_beg,jglobal_beg,jglobal_end,          &
                     kglobal_spg_beg,kglobal_beg,kglobal_end
    !
    lspg_i0=.false.
    lspg_im=.false.
    lspg_j0=.false.
    lspg_jm=.false.
    lspg_k0=.false.
    lspg_km=.false.
    !
    if(spg_i0>0) then
      lspg_i0=.true.
    endif
    !
    if(spg_im>0) then
      lspg_im=.true.
    endif
    !
    if(spg_j0>0) then
      lspg_j0=.true.
    endif
    !
    if(spg_jm>0) then
      lspg_jm=.true.
    endif
    !
    if(spg_k0>0) then
      lspg_k0=.true.
    endif
    !
    if(spg_km>0) then
      lspg_km=.true.
    endif
    !
    if(mpirank==0) then
      !
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                        *** sponge layer ***'
      if(lspg_i0 .or. lspg_im) then
        write(*,'(22X,3(A,I4))')'i direction:   0 ~',spg_i0,         &
                                        ' ....... ',ia-spg_im,' ~ ',ia
      endif
      if(lspg_j0 .or. lspg_jm) then
        write(*,'(22X,3(A,I4))')'j direction:   0 ~',spg_j0,         &
                                        ' ....... ',ja-spg_jm,' ~ ',ja
      endif
      if(lspg_k0 .or. lspg_km) then
        write(*,'(22X,3(A,I4))')'k direction:   0 ~',spg_k0,         &
                                        ' ....... ',ka-spg_km,' ~ ',ka
      endif
      write(*,'(2X,62A)')('-',i=1,62)
      !
    endif
    !
    if(lspg_i0) then
      call layer_setup(dir='i0',spg_beg=spg_i0_beg, &
                                spg_end=spg_i0_end, &
                                sponge_damp_coef=sponge_damp_coef_i0)
    endif
    if(lspg_im) then
      call layer_setup(dir='im',spg_beg=spg_im_beg, &
                                spg_end=spg_im_end, &
                                sponge_damp_coef=sponge_damp_coef_im)
    endif
    if(lspg_jm) then
      call layer_setup(dir='jm',spg_beg=spg_jm_beg, &
                                spg_end=spg_jm_end, &
                                sponge_damp_coef=sponge_damp_coef_jm)
    endif
    if(lspg_k0) then
      call layer_setup(dir='k0',spg_beg=spg_k0_beg, &
                                spg_end=spg_k0_end, &
                                sponge_damp_coef=sponge_damp_coef_k0)
    endif

    if(lspg_km) then
      call layer_setup(dir='km',spg_beg=spg_km_beg, &
                                spg_end=spg_km_end, &
                                sponge_damp_coef=sponge_damp_coef_km)
    endif
    !
  end subroutine spongelayer_define_ijk
  !+-------------------------------------------------------------------+
  !| The end of the subroutine spongelayer_define_ijk.                 |
  !+-------------------------------------------------------------------+
  !
  subroutine layer_setup(dir,spg_beg,spg_end,sponge_damp_coef)

    use commvar, only : im,jm,km,is,ie,js,je,ks,ke
    use parallel,only : ig0,jg0,kg0,mpi_igroup,mpi_jgroup,mpi_kgroup, &
                        psum,mpileft,mpiright,mpidown,mpiup,mpifront, &
                        mpiback,precv,psend,ia,ja,ka,mpirank,         &
                        mpirankname,mpistop
    use commarray,only: x
    use tecio

    character(len=*),intent(in) :: dir
    integer,intent(out) :: spg_beg,spg_end
    real(8),allocatable,intent(out) :: sponge_damp_coef(:,:,:)

    integer :: i,j,k,global_spg_beg,global_beg,global_end
    real(8),allocatable :: length_sponge(:,:),length_sponge_b(:,:),    &
                           length_sponge_local(:,:)

    real(8),parameter :: dampfac=0.05d0
    real(8) :: var1,var2,var3,dis

    if(dir=='i0') then

      global_spg_beg=spg_i0

      global_beg=ig0
      global_end=ig0+im
      !
      if(global_end<=global_spg_beg) then
        ! the sponger layer is completely within the domain
        spg_beg=is
        spg_end=ie
      elseif(global_beg>global_spg_beg) then
        ! the sponger layer is not within the domain
        spg_beg=-1
        spg_end=-1
      elseif(global_beg<=global_spg_beg .and. global_end>=global_spg_beg) then
        ! the sponger layer is partly within the domain
        spg_beg=is
        spg_end=global_spg_beg-global_beg
      else
        stop ' error 1: local domain define error @ spongelayerini'
      endif
      !
      allocate( length_sponge(0:jm,0:km),length_sponge_b(0:jm,0:km), &
                length_sponge_local(0:jm,0:km) )
      allocate( sponge_damp_coef(spg_beg:spg_end,0:jm,0:km) )

      length_sponge      =0.d0
      length_sponge_b    =0.d0
      length_sponge_local=0.d0

      call precv(vario=length_sponge_b,recv_dir=mpiright,tag=21)

      if(spg_end>0) then
        !
        do k=0,km
        do j=0,jm
          !
          do i=spg_end-1,spg_beg,-1
            length_sponge_local(j,k)=length_sponge_local(j,k)+     &
                               sqrt((x(i+1,j,k,1)-x(i,j,k,1))**2 + &
                                    (x(i+1,j,k,2)-x(i,j,k,2))**2 + &
                                    (x(i+1,j,k,3)-x(i,j,k,3))**2)
          enddo
          !
        enddo
        enddo
        !
        length_sponge=length_sponge_b+length_sponge_local
        !
      endif
      !
      call psend(varin=length_sponge,send_dir=mpileft,tag=21)
      !
      length_sponge=psum(length_sponge_local,comm=mpi_igroup)
      !
      if(spg_end>0) then

        do k=0,km
        do j=0,jm
          !
          dis =length_sponge_b(j,k)
          var2=length_sponge(j,k)**2
          do i=spg_end,spg_beg,-1
            if(i<spg_end) then
              var3=sqrt((x(i+1,j,k,1)-x(i,j,k,1))**2 + &
                        (x(i+1,j,k,2)-x(i,j,k,2))**2 + &
                        (x(i+1,j,k,3)-x(i,j,k,3))**2)
            else
              var3=0.d0
            endif
            dis=dis + var3
            sponge_damp_coef(i,j,k)=dampfac*dis**2/var2
          enddo
          !
        enddo
        enddo

      endif

      deallocate(length_sponge,length_sponge_b,length_sponge_local)

    endif

    ! ----------------------------------------------------------------

    if(dir=='im') then

      global_spg_beg=ia-spg_im

      global_beg=ig0
      global_end=ig0+im
      !
      if(global_spg_beg<=global_beg) then
        ! the sponger layer is completely within the domain
        spg_beg=is
        spg_end=ie
      elseif(global_spg_beg>=global_beg .and. global_spg_beg<=global_end) then
        ! the sponger layer is partly within the domain
        spg_im_beg=global_spg_beg-global_beg
        spg_im_end=ie
      elseif(global_spg_beg>global_end) then
        ! the sponger layer is not within the domain
        spg_beg=-1
        spg_end=-1
      else
        stop ' error 2: local domain define error @ spongelayerini'
      endif
      !
      allocate( length_sponge(0:jm,0:km),length_sponge_b(0:jm,0:km), &
                length_sponge_local(0:jm,0:km) )
      allocate( sponge_damp_coef(spg_beg:spg_end,0:jm,0:km) )

      length_sponge      =0.d0
      length_sponge_b    =0.d0
      length_sponge_local=0.d0

      call precv(vario=length_sponge_b,recv_dir=mpileft,tag=22)

      if(spg_end>0) then
        !
        do k=0,km
        do j=0,jm
          !
          do i=spg_beg+1,spg_end
            length_sponge_local(j,k)=length_sponge_local(j,k)+     &
                               sqrt((x(i,j,k,1)-x(i-1,j,k,1))**2 + &
                                    (x(i,j,k,2)-x(i-1,j,k,2))**2 + &
                                    (x(i,j,k,3)-x(i-1,j,k,3))**2)
          enddo
          !
        enddo
        enddo
        !
        length_sponge=length_sponge_b+length_sponge_local
        !
      endif
      !
      call psend(varin=length_sponge,send_dir=mpiright,tag=22)
      !
      length_sponge=psum(length_sponge_local,comm=mpi_igroup)
      !
      if(spg_end>0) then

        do k=0,km
        do j=0,jm
          !
          dis =length_sponge_b(j,k)
          var2=length_sponge(j,k)**2
          do i=spg_beg,spg_end
            if(i>spg_beg) then
              var3=sqrt((x(i,j,k,1)-x(i-1,j,k,1))**2 + &
                        (x(i,j,k,2)-x(i-1,j,k,2))**2 + &
                        (x(i,j,k,3)-x(i-1,j,k,3))**2)
            else
              var3=0.d0
            endif
            dis=dis + var3
            sponge_damp_coef(i,j,k)=dampfac*dis**2/var2
          enddo
          !
        enddo
        enddo

      endif

      deallocate(length_sponge,length_sponge_b,length_sponge_local)

    endif

    ! ----------------------------------------------------------------

    if(dir=='jm') then

      global_spg_beg=ja-spg_jm

      global_beg=jg0
      global_end=jg0+jm
      !
      if(global_spg_beg<=global_beg) then
        ! the sponger layer is completely within the domain
        spg_beg=js
        spg_end=je
      elseif(global_spg_beg>=global_beg .and. global_spg_beg<=global_end) then
        ! the sponger layer is partly within the domain
        spg_beg=global_spg_beg-global_beg
        spg_end=je
      elseif(global_spg_beg>global_end) then
        ! the sponger layer is not within the domain
        spg_beg=-1
        spg_end=-1
      else
        stop ' error 2: local domain define error @ spongelayerini'
      endif
      !
      allocate( length_sponge(0:im,0:km),length_sponge_b(0:im,0:km), &
                length_sponge_local(0:im,0:km) )
      allocate( sponge_damp_coef(0:im,spg_beg:spg_end,0:km) )

      length_sponge      =0.d0
      length_sponge_b    =0.d0
      length_sponge_local=0.d0

      call precv(vario=length_sponge_b,recv_dir=mpidown,tag=23)
      
      if(spg_end>0) then
        !
        do k=0,km
        do i=0,im
          !
          do j=spg_beg+1,spg_end
            length_sponge_local(i,k)=length_sponge_local(i,k)+     &
                               sqrt((x(i,j,k,1)-x(i,j-1,k,1))**2 + &
                                    (x(i,j,k,2)-x(i,j-1,k,2))**2 + &
                                    (x(i,j,k,3)-x(i,j-1,k,3))**2)
          enddo
          !
        enddo
        enddo
        !
        length_sponge=length_sponge_b+length_sponge_local
        !
      endif
      !
      call psend(varin=length_sponge,send_dir=mpiup,tag=23)
      !
      length_sponge=psum(length_sponge_local,comm=mpi_jgroup)
      !
      if(spg_end>0) then

        do k=0,km
        do i=0,im
          !
          dis =length_sponge_b(i,k)
          var2=length_sponge(i,k)**2
          do j=spg_beg,spg_end
            if(j>spg_beg) then
              var3=sqrt((x(i,j,k,1)-x(i,j-1,k,1))**2 + &
                        (x(i,j,k,2)-x(i,j-1,k,2))**2 + &
                        (x(i,j,k,3)-x(i,j-1,k,3))**2)
            else
              var3=0.d0
            endif
            dis=dis + var3
            sponge_damp_coef(i,j,k)=dampfac*dis**2/var2
          enddo
          !
        enddo
        enddo

      endif

      deallocate(length_sponge,length_sponge_b,length_sponge_local)

    endif

    ! ----------------------------------------------------------------

    if(dir=='k0') then

      global_spg_beg=spg_k0

      global_beg=kg0
      global_end=kg0+km
      !
      if(global_end<=global_spg_beg) then
        ! the sponger layer is completely within the domain
        spg_beg=ks
        spg_end=ke
      elseif(global_beg>global_spg_beg) then
        ! the sponger layer is not within the domain
        spg_beg=-1
        spg_end=-1
      elseif(global_beg<=global_spg_beg .and. global_end>=global_spg_beg) then
        ! the sponger layer is partly within the domain
        spg_beg=ks
        spg_end=global_spg_beg-global_beg
      else
        stop ' error 1: local domain define error @ spongelayerini'
      endif
      
      allocate( length_sponge(0:im,0:jm),length_sponge_b(0:im,0:jm), &
                length_sponge_local(0:im,0:jm) )
      allocate( sponge_damp_coef(0:im,0:jm,spg_beg:spg_end) )

      length_sponge      =0.d0
      length_sponge_b    =0.d0
      length_sponge_local=0.d0

      call precv(vario=length_sponge_b,recv_dir=mpifront,tag=25)

      if(spg_end>0) then
        !
        do j=0,jm
        do i=0,im
          !
          do k=spg_end-1,spg_beg,-1
            length_sponge_local(i,j)=length_sponge_local(i,j)+     &
                               sqrt((x(i,j,k+1,1)-x(i,j,k,1))**2 + &
                                    (x(i,j,k+1,2)-x(i,j,k,2))**2 + &
                                    (x(i,j,k+1,3)-x(i,j,k,3))**2)
          enddo
          !
        enddo
        enddo
        !
        length_sponge=length_sponge_b+length_sponge_local
        !
      endif
      !
      call psend(varin=length_sponge,send_dir=mpiback,tag=25)
      !
      length_sponge=psum(length_sponge_local,comm=mpi_kgroup)
      !
      if(spg_end>0) then

        do j=0,jm
        do i=0,im
          !
          dis =length_sponge_b(j,k)
          var2=length_sponge(j,k)**2
          do k=spg_end,spg_beg,-1
            if(k<spg_end) then
              var3=sqrt((x(i,j,k+1,1)-x(i,j,k,1))**2 + &
                        (x(i,j,k+1,2)-x(i,j,k,2))**2 + &
                        (x(i,j,k+1,3)-x(i,j,k,3))**2)
            else
              var3=0.d0
            endif
            dis=dis + var3
            sponge_damp_coef(i,j,k)=dampfac*dis**2/var2
          enddo
          !
        enddo
        enddo

      endif

      ! if(spg_end>0) then
      ! call tecbin('testout/tecsponge_k0'//mpirankname//'.plt',                &
      !                                   x(0:im,0:jm,spg_beg:spg_end,1),'x',         &
      !                                   x(0:im,0:jm,spg_beg:spg_end,2),'y',         &
      !                                   x(0:im,0:jm,spg_beg:spg_end,3),'z',         &
      !                    sponge_damp_coef(0:im,0:jm,spg_beg:spg_end),'sp' )

      ! endif
      deallocate(length_sponge,length_sponge_b,length_sponge_local)

    endif

    ! ----------------------------------------------------------------
    
    if(dir=='km') then

      global_spg_beg=ka-spg_km

      global_beg=kg0
      global_end=kg0+km
      !
      if(global_spg_beg<=global_beg) then
        ! the sponger layer is completely within the domain
        spg_beg=ks
        spg_end=ke
      elseif(global_spg_beg>=global_beg .and. global_spg_beg<=global_end) then
        ! the sponger layer is partly within the domain
        spg_beg=global_spg_beg-global_beg
        spg_end=ke
      elseif(global_spg_beg>global_end) then
        ! the sponger layer is not within the domain
        spg_beg=-1
        spg_end=-1
      else
        stop ' error 2: local domain define error @ spongelayerini'
      endif
      !
      allocate( length_sponge(0:im,0:jm),length_sponge_b(0:im,0:jm), &
                length_sponge_local(0:im,0:jm) )
      allocate( sponge_damp_coef(0:im,0:jm,spg_beg:spg_end) )

      length_sponge      =0.d0
      length_sponge_b    =0.d0
      length_sponge_local=0.d0

      call precv(vario=length_sponge_b,recv_dir=mpiback,tag=26)
      
      if(spg_end>0) then
        !
        do j=0,jm
        do i=0,im
          !
          do k=spg_beg+1,spg_end
            length_sponge_local(i,j)=length_sponge_local(i,j)+     &
                               sqrt((x(i,j,k,1)-x(i,j,k-1,1))**2 + &
                                    (x(i,j,k,2)-x(i,j,k-1,2))**2 + &
                                    (x(i,j,k,3)-x(i,j,k-1,3))**2)
          enddo
          !
        enddo
        enddo
        !
        length_sponge=length_sponge_b+length_sponge_local
        !
      endif
      !
      call psend(varin=length_sponge,send_dir=mpifront,tag=26)
      !
      length_sponge=psum(length_sponge_local,comm=mpi_kgroup)
      !
      if(spg_end>0) then

        do j=0,jm
        do i=0,im
          !
          dis =length_sponge_b(i,j)
          var2=length_sponge(i,j)**2
          do k=spg_beg,spg_end
            if(k>spg_beg) then
              var3=sqrt((x(i,j,k,1)-x(i,j,k-1,1))**2 + &
                        (x(i,j,k,2)-x(i,j,k-1,2))**2 + &
                        (x(i,j,k,3)-x(i,j,k-1,3))**2)
            else
              var3=0.d0
            endif
            dis=dis + var3
            sponge_damp_coef(i,j,k)=dampfac*dis**2/var2
          enddo
          !
        enddo
        enddo

      endif

      deallocate(length_sponge,length_sponge_b,length_sponge_local)

      ! if(spg_end>0) then
      ! call tecbin('testout/tecsponge_km'//mpirankname//'.plt',                &
      !                                   x(0:im,0:jm,spg_beg:spg_end,1),'x',         &
      !                                   x(0:im,0:jm,spg_beg:spg_end,2),'y',         &
      !                                   x(0:im,0:jm,spg_beg:spg_end,3),'z',         &
      !                    sponge_damp_coef(0:im,0:jm,spg_beg:spg_end),'sp' )

      ! endif

    endif


  end subroutine layer_setup

end module sponge_layer