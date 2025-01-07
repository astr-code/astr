
  subroutine convrsdflux(timerept)
    !
    use commvar,  only: im,jm,km,hm,numq,num_species,num_modequ,       &
                        conschm,npdci,npdcj,npdck,is,ie,js,je,ks,ke
    use commarray,only: q,vel,rho,prs,tmp,spc,dxi,jacob,qrhs
    use commfunc, only: ddfc,stable_euler_flux
    use comsolver, only : alfa_con,cci,ccj,cck
    !
    ! arguments
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: i,j,k,jspc,jmod,n,nterm
    real(8),allocatable :: fcs(:,:),dfcs(:,:),uu(:)
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    nterm=numq+5
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(fcs(1:nterm,-hm:im+hm),dfcs(0:im,1:numq),uu(-hm:im+hm))
    do k=ks,ke
    do j=js,je
      !
      fcs(:,1) =dxi(:,j,k,1,1)*vel(:,j,k,1)+dxi(:,j,k,1,2)*vel(:,j,k,2) +  &
                dxi(:,j,k,1,3)*vel(:,j,k,3)
      fcs(:,2) =rho(:,j,k)
      fcs(:,3) =vel(:,j,k,1)
      fcs(:,4) =vel(:,j,k,2)
      fcs(:,5) =vel(:,j,k,3)
      fcs(:,6) =(q(:,j,k,5)+prs(:,j,k))/rho(:,j,k)
      fcs(:,7) =prs(:,j,k)
      fcs(:,8) =dxi(:,j,k,1,1)
      fcs(:,9) =dxi(:,j,k,1,2)
      fcs(:,10)=dxi(:,j,k,1,3)

      if(num_species>0) then
        do jspc=1,num_species
          fcs(:,11+jspc)=q(:,j,k,numq+jspc)/rho(:,j,k)
        enddo
      endif
      
      dfcs(:,:)=stable_euler_flux(fcs,npdci,nterm,im)

      do i=is,ie
        qrhs(i,j,k,:)=qrhs(i,j,k,:)+jacob(i,j,k)*dfcs(i,:)
      enddo

    enddo
    enddo
    deallocate(fcs,dfcs,uu)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ndims>=2) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! calculating along j direction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(fcs(-hm:jm+hm,1:numq),dfcs(0:jm,1:numq),uu(-hm:jm+hm))
      do k=ks,ke
      do i=is,ie
        !
        uu(:)=dxi(i,:,k,2,1)*vel(i,:,k,1)+dxi(i,:,k,2,2)*vel(i,:,k,2) +  &
              dxi(i,:,k,2,3)*vel(i,:,k,3)
        fcs(:,1)=jacob(i,:,k)*  q(i,:,k,1)*uu
        fcs(:,2)=jacob(i,:,k)*( q(i,:,k,2)*uu+dxi(i,:,k,2,1)*prs(i,:,k) )
        fcs(:,3)=jacob(i,:,k)*( q(i,:,k,3)*uu+dxi(i,:,k,2,2)*prs(i,:,k) )
        fcs(:,4)=jacob(i,:,k)*( q(i,:,k,4)*uu+dxi(i,:,k,2,3)*prs(i,:,k) )
        fcs(:,5)=jacob(i,:,k)*( q(i,:,k,5)+prs(i,:,k) )*uu
        !
        if(num_species>0) then
          n=5
          do jspc=1,num_species
            fcs(:,n+jspc)=jacob(i,:,k)*q(i,:,k,n+jspc)*uu
          enddo
        endif
        !
        if(num_modequ>0) then
          n=5+num_species
          do jmod=1,num_modequ
            fcs(:,n+jmod)=jacob(i,:,k)*q(i,:,k,n+jmod)*uu
          enddo
        endif
        !
        do n=1,numq
          dfcs(:,n)=ddfc(fcs(:,n),conschm,npdcj,jm,alfa_con,ccj)
        enddo
        !
        qrhs(i,js:je,k,:)=qrhs(i,js:je,k,:)+dfcs(js:je,:)
        !
        !
      enddo
      enddo
      deallocate(fcs,dfcs,uu)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! end calculating along j direction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    if(ndims==3) then
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! calculating along j direction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(fcs(-hm:km+hm,1:numq),dfcs(0:km,1:numq),uu(-hm:km+hm))
      do j=js,je
      do i=is,ie
        !
        uu(:)=dxi(i,j,:,3,1)*vel(i,j,:,1)+dxi(i,j,:,3,2)*vel(i,j,:,2) +  &
              dxi(i,j,:,3,3)*vel(i,j,:,3)
        fcs(:,1)=jacob(i,j,:)*  q(i,j,:,1)*uu
        fcs(:,2)=jacob(i,j,:)*( q(i,j,:,2)*uu+dxi(i,j,:,3,1)*prs(i,j,:) )
        fcs(:,3)=jacob(i,j,:)*( q(i,j,:,3)*uu+dxi(i,j,:,3,2)*prs(i,j,:) )
        fcs(:,4)=jacob(i,j,:)*( q(i,j,:,4)*uu+dxi(i,j,:,3,3)*prs(i,j,:) )
        fcs(:,5)=jacob(i,j,:)*( q(i,j,:,5)+prs(i,j,:) )*uu
        !
        if(num_species>0) then
          n=5
          do jspc=1,num_species
            fcs(:,n+jspc)=jacob(i,j,:)*q(i,j,:,n+jspc)*uu
          enddo
        endif
        !
        if(num_modequ>0) then
          n=5+num_species
          do jmod=1,num_modequ
            fcs(:,n+jmod)=jacob(i,j,:)*q(i,j,:,n+jmod)*uu
          enddo
        endif
        !
        do n=1,numq
          dfcs(:,n)=ddfc(fcs(:,n),conschm,npdck,km,alfa_con,cck,lfft=lfftk)
        enddo
        !
        qrhs(i,j,ks:ke,:)=qrhs(i,j,ks:ke,:)+dfcs(ks:ke,:)
        !
      enddo
      enddo
      deallocate(fcs,dfcs,uu)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! end calculating along j direction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
    endif
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='convrsdcal6', &
                                             timecost=subtime, &
                                              message='diffusion term with central scheme')
    endif
    !
    return
    !
  end subroutine convrsdflux