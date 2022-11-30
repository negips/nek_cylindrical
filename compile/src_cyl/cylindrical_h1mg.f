!====================================================================== 
!     Author: Prabal Negi
!     Description: Coarse grid preconditioner for cylindrical solve.
!      
!      
!====================================================================== 
      subroutine get_local_crs_galerkin_cyl(a,ncl,nxc,h1,h2,w1,w2)

c     This routine generates Nelv submatrices of order ncl using
c     Galerkin projection

      implicit none

      include 'SIZE'

      integer ncl,nxc

      real    a(ncl,ncl,1),h1(1),h2(1)
      real    w1(lx1*ly1*lz1,nelv),w2(lx1*ly1*lz1,nelv)

      integer lcrd
      parameter (lcrd=lx1**ldim)

      real b
      common /ctmp1z/ b(lcrd,8)

      integer i,j,e
      integer isd,imsh,nxyz
      real vlsc2

      do j=1,ncl
         call gen_crs_basis(b(1,j),j) ! bi- or tri-linear interpolant
      enddo

      isd  = 1
      imsh = 1

      nxyz = lx1*ly1*lz1
      do j = 1,ncl
         do e = 1,nelv
            call copy(w1(1,e),b(1,j),nxyz)
         enddo

         call axhelm (w2,w1,h1,h2,imsh,isd)        ! A^e * bj

         do e = 1,nelv
         do i = 1,ncl
            a(i,j,e) = vlsc2(b(1,i),w2(1,e),nxyz)  ! bi^T * A^e * bj
         enddo
         enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine axhelm_cyl_pr (au,u,helm1,helm2,imesh,isd)

!     Compute the (Helmholtz) matrix-vector product,
!     AU = helm1*[A]u + helm2*[B]u, for NEL elements.

      implicit none

      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
      include 'PARALLEL'
      include 'CTIMER'

      real    wddx,wddyt,wddzt
      COMMON /FASTAX/ WDDX(LX1,LX1),WDDYT(LY1,LY1),WDDZT(LZ1,LZ1)
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV

      REAL           AU    (LX1,LY1,LZ1,1)
     $ ,             U     (LX1,LY1,LZ1,1)
     $ ,             HELM1 (LX1,LY1,LZ1,1)
     $ ,             HELM2 (LX1,LY1,LZ1,1)

      real dudr,duds,dudt,tmp1,tmp2,tmp3
      COMMON /CTMP1/ DUDR  (LX1,LY1,LZ1)
     $ ,             DUDS  (LX1,LY1,LZ1)
     $ ,             DUDT  (LX1,LY1,LZ1)
     $ ,             TMP1  (LX1,LY1,LZ1)
     $ ,             TMP2  (LX1,LY1,LZ1)
     $ ,             TMP3  (LX1,LY1,LZ1)

      REAL           TM1   (LX1,LY1,LZ1)
      REAL           TM2   (LX1,LY1,LZ1)
      REAL           TM3   (LX1,LY1,LZ1)
      REAL           DUAX  (LX1)
      REAL           YSM1  (LX1)
      EQUIVALENCE    (DUDR,TM1),(DUDS,TM2),(DUDT,TM3)

      integer i,j,iz,e
      integer imesh,isd
      integer nel,nxy,nyz,nxz,nxyz,ntot
      real h1,term1,term2


      naxhm = naxhm + 1
      etime1 = dnekclock()

      nel=nelt
      if (imesh.eq.1) nel=nelv

      NXY=lx1*ly1
      NYZ=ly1*lz1
      NXZ=lx1*lz1
      NXYZ=lx1*ly1*lz1
      NTOT=NXYZ*NEL

      IF (.NOT.IFSOLV) CALL SETFAST(HELM1,HELM2,IMESH)
      CALL RZERO (AU,NTOT)

      do 100 e=1,nel
C
        if (ifaxis) call setaxdy ( ifrzer(e) )
C
        IF (ldim.EQ.2) THEN
C
C          2-d case ...............
C
           if (iffast(e)) then
C
C          Fast 2-d mode: constant properties and undeformed element
C
           h1 = helm1(1,1,1,e)
           call mxm   (wddx,lx1,u(1,1,1,e),lx1,tm1,nyz)
           call mxm   (u(1,1,1,e),lx1,wddyt,ly1,tm2,ly1)
           call col2  (tm1,g4m1(1,1,1,e),nxyz)
           call col2  (tm2,g5m1(1,1,1,e),nxyz)
           call add3  (au(1,1,1,e),tm1,tm2,nxyz)
           call cmult (au(1,1,1,e),h1,nxyz)
C
           else
C
C
           call mxm  (dxm1,lx1,u(1,1,1,e),lx1,dudr,nyz)
           call mxm  (u(1,1,1,e),lx1,dytm1,ly1,duds,ly1)
           call col3 (tmp1,dudr,g1m1(1,1,1,e),nxyz)
           call col3 (tmp2,duds,g2m1(1,1,1,e),nxyz)
           if (ifdfrm(e)) then
              call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
              call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
           endif
           call col2 (tmp1,helm1(1,1,1,e),nxyz)
           call col2 (tmp2,helm1(1,1,1,e),nxyz)
           call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
           call mxm  (tmp2,lx1,dym1,ly1,tm2,ly1)
           call add2 (au(1,1,1,e),tm1,nxyz)
           call add2 (au(1,1,1,e),tm2,nxyz)

        endif
C
        else
C
C       3-d case ...............
C
           if (iffast(e)) then
C
C          Fast 3-d mode: constant properties and undeformed element
C
           h1 = helm1(1,1,1,e)
           call mxm   (wddx,lx1,u(1,1,1,e),lx1,tm1,nyz)
           do 5 iz=1,lz1
           call mxm   (u(1,1,iz,e),lx1,wddyt,ly1,tm2(1,1,iz),ly1)
 5         continue
           call mxm   (u(1,1,1,e),nxy,wddzt,lz1,tm3,lz1)
           call col2  (tm1,g4m1(1,1,1,e),nxyz)
           call col2  (tm2,g5m1(1,1,1,e),nxyz)
           call col2  (tm3,g6m1(1,1,1,e),nxyz)
           call add3  (au(1,1,1,e),tm1,tm2,nxyz)
           call add2  (au(1,1,1,e),tm3,nxyz)
           call cmult (au(1,1,1,e),h1,nxyz)
C
           else
C
C
           call mxm(dxm1,lx1,u(1,1,1,e),lx1,dudr,nyz)
           do 10 iz=1,lz1
              call mxm(u(1,1,iz,e),lx1,dytm1,ly1,duds(1,1,iz),ly1)
   10      continue
           call mxm     (u(1,1,1,e),nxy,dztm1,lz1,dudt,lz1)
           call col3    (tmp1,dudr,g1m1(1,1,1,e),nxyz)
           call col3    (tmp2,duds,g2m1(1,1,1,e),nxyz)
           call col3    (tmp3,dudt,g3m1(1,1,1,e),nxyz)
           if (ifdfrm(e)) then
              call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
              call addcol3 (tmp1,dudt,g5m1(1,1,1,e),nxyz)
              call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
              call addcol3 (tmp2,dudt,g6m1(1,1,1,e),nxyz)
              call addcol3 (tmp3,dudr,g5m1(1,1,1,e),nxyz)
              call addcol3 (tmp3,duds,g6m1(1,1,1,e),nxyz)
           endif
           call col2 (tmp1,helm1(1,1,1,e),nxyz)
           call col2 (tmp2,helm1(1,1,1,e),nxyz)
           call col2 (tmp3,helm1(1,1,1,e),nxyz)
           call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
           do 20 iz=1,lz1
              call mxm(tmp2(1,1,iz),lx1,dym1,ly1,tm2(1,1,iz),ly1)
   20      continue
           call mxm  (tmp3,nxy,dzm1,lz1,tm3,lz1)
           call add2 (au(1,1,1,e),tm1,nxyz)
           call add2 (au(1,1,1,e),tm2,nxyz)
           call add2 (au(1,1,1,e),tm3,nxyz)
C
           endif
c
        endif
C
 100  continue
C
      if (ifh2) call addcol4 (au,helm2,bm1,u,ntot)

!     Add contributions of additional terms in cylindrical formulation
!     Not done yet      

      taxhm=taxhm+(dnekclock()-etime1)
      return
      end
C
c=======================================================================






























