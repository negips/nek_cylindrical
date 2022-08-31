!======================================================================
!     Author: Prabal Negi      
!     Description: Routines for 2D polar coordinates implementation
!
!====================================================================== 

       subroutine opgradt_polar(outx,outy,outz,inpfld)
C------------------------------------------------------------------------
C
C     Compute DTx, DTy, DTz of an input field INPFLD 
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      real outx   (lx1,ly1,lz1,1)
      real outy   (lx1,ly1,lz1,1)
      real outz   (lx1,ly1,lz1,1)
      real inpfld (lx2,ly2,lz2,1)
C
      call cdtp_polar (outx,inpfld,rxm2,sxm2,txm2,1)
      call cdtp_polar (outy,inpfld,rym2,sym2,tym2,2)
C
      return
      end
c-----------------------------------------------------------------------

      subroutine cdtp_polar (dtx,x,rm2,sm2,tm2,isd)

!     Compute DT*X (entire field)
!     I assume y == Radial
!     and      x == \theta

!     I assume ifaxis == False
!     But BM2 includes a factor of R.        
!     I also assume we are not dealing with the singular case
!     R=0 (yet).        

      implicit none  
  
      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'IXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
      include 'ESOLV'
      include 'CTIMER'
C
      real dtx  (lx1*ly1*lz1,lelv)
      real x    (lx2*ly2*lz2,lelv)
      real rm2  (lx2*ly2*lz2,lelv)
      real sm2  (lx2*ly2*lz2,lelv)
      real tm2  (lx2*ly2*lz2,lelv)

      real wx,ta1,ta2,ta3
      common /ctmp1/ wx  (lx1*ly1*lz1)
     $ ,             ta1 (lx1*ly1*lz1)
     $ ,             ta2 (lx1*ly1*lz1)
     $ ,             ta3 (lx1*ly1,lz1)

      REAL           DUAX(LX1)

      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV

      integer e

      integer isd,i1,i2,ix,iy,iz
      integer n1,n2,ly12,nxy1,nxyz1,nxyz2,nyz2

#ifdef TIMER
      if (icalld.eq.0) tcdtp=0.0
      icalld=icalld+1
      ncdtp=icalld
      etime1=dnekclock()
#endif

      nxyz1 = lx1*ly1*lz1
      nxyz2 = lx2*ly2*lz2
      nyz2  = ly2*lz2
      nxy1  = lx1*ly1

      n1    = lx1*ly1
      n2    = lx1*ly2

      do e=1,nelv

C       Use the appropriate derivative- and interpolation operator in 
C       the y-direction (= radial direction if axisymmetric).
        ly12   = ly1*ly2
        if (ifrzer(e)) then
           call copy (iym12,iam12,ly12)
           call copy (dym12,dam12,ly12)
           call copy (w3m2,w2am2,nxyz2)
        else
           call copy (iym12,icm12,ly12)
           call copy (dym12,dcm12,ly12)
           call copy (w3m2,w2cm2,nxyz2)
        endif

!       Collocate with weights
        if(ifsplit) then

        else
!         if (.not.ifaxis) call col3 (wx,w3m2,x(1,e),nxyz2)

          if (ifrzer(e)) then
              call col3    (wx,x(1,e),bm2(1,1,1,e),nxyz2)
              call invcol2 (wx,jacm2(1,1,1,e),nxyz2)
          else
              call col3    (wx,w3m2,x(1,e),nxyz2)
              call col2    (wx,ym2(1,1,1,e),nxyz2)
          endif
        endif

        call col3 (ta1,wx,rm2(1,e),nxyz2)
        call mxm  (dxtm12,lx1,ta1,lx2,ta2,nyz2)
        call mxm  (ta2,lx1,iym12,ly2,dtx(1,e),ly1)

        call col3 (ta1,wx,sm2(1,e),nxyz2)
        call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)
        call mxm  (ta2,lx1,dym12,ly2,ta1,ly1)

        call add2 (dtx(1,e),ta1,nxyz1)

C     Add an extra diagonal term in the radial/theta directions
C     direction (only if solving the momentum equations and ISD=2)
C     NOTE: lz1=lz2=1

      if (isd.eq.2) then         
        call col3    (ta1,x(1,e),bm2(1,1,1,e),nxyz2)
        call invcol2 (ta1,ym2(1,1,1,e),nxyz2)
        call mxm     (ixtm12,lx1,ta1,lx2,ta2,ly2)
        call mxm     (ta2,lx1,iym12,ly2,ta1,ly1)
        call add2    (dtx(1,e),ta1,nxyz1)
      endif

      enddo
C
#ifdef TIMER
      tcdtp=tcdtp+(dnekclock()-etime1)
#endif
      return
      end
!---------------------------------------------------------------------- 

