!====================================================================== 
!
!     Author: Prabal Negi
!     Description: Generic subroutines that are frequently used
!     
!     Routines:
!     ortho_subspace          : Orthogonalize vector with subspace
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

