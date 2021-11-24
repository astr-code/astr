
program mpisizedis
    !
    implicit none
    !
    logical :: lallo
    integer :: nfactor(1:30)
    integer :: i,j,k,n,n1,n2,n3,nsize,kaa
    integer(8) :: nvar1,nvar2,nvar3
    !
    integer :: ndims,ia,ja,ka,isize,jsize,ksize,mpisize,mpirank
    logical :: lfftk
    !
    ndims=3
    mpisize=128*128
    ia=2240
    ja=160
    ka=768
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
          !
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
          nvar1=int(ja,8)*int(kaa,8)*int(nfactor(n1),8)
          !
          nvar1=nvar1+int(ia,8)*int(kaa,8)*int(nfactor(n2),8)
          !
          nvar1=nvar1+int(ia,8)*int(ja,8)*int(nfactor(n3),8)
          !
          print*,nvar1,nvar2,'-',nfactor(n1),nfactor(n2),nfactor(n3),'|',ia*kaa*nfactor(n2),ia*kaa,nfactor(n2)
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
  end program mpisizedis