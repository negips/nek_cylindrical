!====================================================================== 
!
!     Author: Prabal Negi
!     Description: Generic subroutines that are frequently used
!     
!     Routines:
!     ortho_subspace          : Orthogonalize vector with subspace
!     tensor3_op              : Tensor operations, opx,opyt,opzt
!     tensorx_op              : Tensor operation, opx
!     tensory_op              : Tensor operation, opyt
!     tensorz_op              : Tensor operation, opzt
!
!      
!====================================================================== 
      subroutine ortho_subspace(r,nt,h,V,ldv,k,wgt,ifwgt,ngs,wk1,wk2)

      implicit none

      integer nt              ! Length of the vector r
      integer ldv             ! Leading dimension of V
      integer k               ! No of Columns in V
      real r(nt)              ! Vector to orthogonalize
      real V(ldv,k)           ! Orthogonalizing Space
      real wgt(nt)            ! Weights
      logical ifwgt           ! If Weighted orthogonalization
      integer ngs             ! No. of Gram-Schmidt
      real h(k)               ! Projections on V

!     Work Arrays      
      real wk1(k)
      real wk2(k)

      integer igs,i

      real vlsc2,vlsc3        ! Functions

!     Zero projections      
      call rzero(h,k)

      do igs = 1,ngs
!       Gram-Schmidt:
        do i=1,k
          if (ifwgt) then
            wk1(i)=vlsc3(r,V(1,i),wgt,nt)       ! wk1 = (Bw,V )
          else
            wk1(i)=vlsc2(r,V(1,i),nt)           ! wk1 = (w,V )
          endif
        enddo                                             
        call gop(wk1,wk2,'+  ',k)               ! sum over all procs

        do i=1,k
          call add2s2(r,V(1,i),-wk1(i),nt)      ! r = r - V*wk1
          h(i) = h(i) + wk1(i)                  ! h = h + wk1 
        enddo
      enddo       ! igs 

      return
      end subroutine ortho_subspace
!---------------------------------------------------------------------- 
      subroutine tensor3_op(fldf,fld,nx,ny,nz,opx,opyt,opzt,mx,my,mz)

      implicit none

      include 'SIZE'

      integer nx,ny,nz
      integer mx  ! No. of rows of opx
      integer my  ! No. of columns of opyt
      integer mz  ! No. of columns of opzt

      real fld  (nx*ny*nz)
      real fldf (mx*my*mz)

!     x operator is applied directly 
      real opx(mx*nx)

!     y,z operators apply as transposed operators 
      real opyt(ny*my)
      real opzt(nz*mz)

      integer lxwk
      parameter (lxwk=lx1+5)        ! set arbitrarily

      real op_wk(lxwk**ldim,2)
      common /tensor_op_work/ op_wk

      if ((mx*my*mz.gt.(lxwk**ldim)) .and. 
     $    (nx*ny*nz.gt.(lxwk**ldim))) then
        if (nio.eq.0) write(6,*) 
     $    'Inadequate workspace size in tensor3_op'
        call exitt
      endif  

      call tensorx_op(op_wk(1,1),fld,nx,ny,nz,opx,mx)
      if (ndim.eq.3) then
        call tensory_op(op_wk(1,2),op_wk(1,1),mx,ny,nz,opyt,my)
        call tensorz_op(fldf,op_wk(1,2),mx,my,nz,opzt,mz)
      else
        call tensory_op(fldf,op_wk(1,1),mx,ny,nz,opyt,my)
      endif  


      return
      end subroutine tensor3_op
!---------------------------------------------------------------------- 

      subroutine tensorx_op(fldo,fld,nx,ny,nz,opx,mx)

      implicit none

      integer nx,ny,nz
      integer mx              ! no of rows

      real fld  (1)
      real fldo (1)

      real opx(mx*nx)

      call mxm  (opx,mx,fld,nx,fldo,ny*nz)

      return
      end subroutine tensorx_op
!---------------------------------------------------------------------- 

      subroutine tensory_op(fldo,fld,nx,ny,nz,opy,my)

      implicit none

      integer nx,ny,nz
      integer my              ! no of columns

      real fld  (1)
      real fldo (1)

      real opy(ny*my)
      integer i,j,k

      i = 1
      j = 1
      do k = 1,nz
        call mxm (fld(i),nx,opy,ny,fldo(j),my)
        i = i + nx*ny
        j = j + nx*my
      enddo  

      return
      end subroutine tensory_op
!---------------------------------------------------------------------- 

      subroutine tensorz_op(fldo,fld,nx,ny,nz,opz,mz)

      implicit none

      include 'SIZE'

      integer nx,ny,nz
      integer mz              ! no of columns

      real fld  (1)
      real fldo (1)

      real opz(nz*mz)

      if (ndim.eq.3) then
        call mxm  (fld,nx*ny,opz,nz,fldo,mz)
      endif  


      return
      end subroutine tensorz_op
!---------------------------------------------------------------------- 

