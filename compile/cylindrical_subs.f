!======================================================================
!     Author: Prabal Negi      
!     Description: Routines for 3D cylindrical solve implementation
!
!====================================================================== 

       subroutine opgradt_cylindrical(outx,outy,outz,inpfld)
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
      call cdtp (outx,inpfld,rxm2,sxm2,txm2,1)
      call cdtp (outy,inpfld,rym2,sym2,tym2,2)
      if (ldim.eq.3) 
     $   call cdtp (outz,inpfld,rzm2,szm2,tzm2,3)
C
      return
      end
c-----------------------------------------------------------------------

      subroutine cdtp_cylindrical (dtx,x,rm2,sm2,tm2,isd)
C-------------------------------------------------------------
C
C     Compute DT*X (entire field)
C
C-------------------------------------------------------------
      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'IXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
      include 'ESOLV'
C
      real dtx  (lx1*ly1*lz1,lelv)
      real x    (lx2*ly2*lz2,lelv)
      real rm2  (lx2*ly2*lz2,lelv)
      real sm2  (lx2*ly2*lz2,lelv)
      real tm2  (lx2*ly2*lz2,lelv)
C
      common /ctmp1/ wx  (lx1*ly1*lz1)
     $ ,             ta1 (lx1*ly1*lz1)
     $ ,             ta2 (lx1*ly1*lz1)
     $ ,             ta3 (lx1*ly1,lz1)

      REAL           DUAX(LX1)
c
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
      include 'CTIMER'

      integer e
C
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
        if (ifaxis) then
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
       endif
C
C      Collocate with weights
C
       if(ifsplit) then
         call col3 (wx,bm1(1,1,1,e),x(1,e),nxyz1)
         call invcol2(wx,jacm1(1,1,1,e),nxyz1)
       else
         if (.not.ifaxis) call col3 (wx,w3m2,x(1,e),nxyz2)
C
         if (ifaxis) then
            if (ifrzer(e)) then
                call col3    (wx,x(1,e),bm2(1,1,1,e),nxyz2)
                call invcol2 (wx,jacm2(1,1,1,e),nxyz2)
            else
                call col3    (wx,w3m2,x(1,e),nxyz2)
                call col2    (wx,ym2(1,1,1,e),nxyz2)
            endif
         endif
       endif
C
       if (ldim.eq.2) then
         if (.not.ifdfrm(e) .and. ifalgn(e)) then
C
            if (      ifrsxy(e).and.isd.eq.1  .or. 
     $           .not.ifrsxy(e).and.isd.eq.2) then
C
               call col3 (ta1,wx,rm2(1,e),nxyz2)
               call mxm  (dxtm12,lx1,ta1,lx2,ta2,nyz2)
               call mxm  (ta2,lx1,iym12,ly2,dtx(1,e),ly1)
            else
               call col3 (ta1,wx,sm2(1,e),nxyz2)
               call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)
               call mxm  (ta2,lx1,dym12,ly2,dtx(1,e),ly1)
            endif
         else
            call col3 (ta1,wx,rm2(1,e),nxyz2)
            call mxm  (dxtm12,lx1,ta1,lx2,ta2,nyz2)
            call mxm  (ta2,lx1,iym12,ly2,dtx(1,e),ly1)

            call col3 (ta1,wx,sm2(1,e),nxyz2)
            call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)
            call mxm  (ta2,lx1,dym12,ly2,ta1,ly1)

            call add2 (dtx(1,e),ta1,nxyz1)
         endif

       else
         if (ifsplit) then

            call col3 (ta1,wx,rm2(1,e),nxyz2)
            call mxm  (dxtm12,lx1,ta1,lx2,dtx(1,e),nyz2)
            call col3 (ta1,wx,sm2(1,e),nxyz2)
            i1 = 1
            i2 = 1
            do iz=1,lz2
               call mxm  (ta1(i2),lx1,dym12,ly2,ta2(i1),ly1)
               i1 = i1 + n1
               i2 = i2 + n2
            enddo
            call add2 (dtx(1,e),ta2,nxyz1)
            call col3 (ta1,wx,tm2(1,e),nxyz2)
            call mxm  (ta1,nxy1,dzm12,lz2,ta2,lz1)
            call add2 (dtx(1,e),ta2,nxyz1)

         else

            call col3 (ta1,wx,rm2(1,e),nxyz2)
            call mxm  (dxtm12,lx1,ta1,lx2,ta2,nyz2)
            i1 = 1
            i2 = 1
            do iz=1,lz2
               call mxm  (ta2(i2),lx1,iym12,ly2,ta1(i1),ly1)
               i1 = i1 + n1
               i2 = i2 + n2
            enddo
            call mxm  (ta1,nxy1,izm12,lz2,dtx(1,e),lz1)

            call col3 (ta1,wx,sm2(1,e),nxyz2)
            call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)
            i1 = 1
            i2 = 1
            do iz=1,lz2
               call mxm  (ta2(i2),lx1,dym12,ly2,ta1(i1),ly1)
               i1 = i1 + n1
               i2 = i2 + n2
            enddo
            call mxm  (ta1,nxy1,izm12,lz2,ta2,lz1)
            call add2 (dtx(1,e),ta2,nxyz1)

            call col3 (ta1,wx,tm2(1,e),nxyz2)
            call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)
            i1 = 1
            i2 = 1
            do iz=1,lz2
               call mxm  (ta2(i2),lx1,iym12,ly2,ta1(i1),ly1)
               i1 = i1 + n1
               i2 = i2 + n2
            enddo
            call mxm  (ta1,nxy1,dzm12,lz2,ta2,lz1)
            call add2 (dtx(1,e),ta2,nxyz1)

         endif

       endif         
C
C     If axisymmetric, add an extra diagonal term in the radial 
C     direction (only if solving the momentum equations and ISD=2)
C     NOTE: lz1=lz2=1
C
C
      if(ifsplit) then

       if (ifaxis.and.(isd.eq.4)) then
        call copy    (ta1,x(1,e),nxyz1)
        if (ifrzer(e)) THEN
           call rzero(ta1, lx1)
           call mxm  (x  (1,e),lx1,datm1,ly1,duax,1)
           call copy (ta1,duax,lx1)
        endif
        call col2    (ta1,baxm1(1,1,1,e),nxyz1)
        call add2    (dtx(1,e),ta1,nxyz1)
       endif

      else

       if (ifaxis.and.(isd.eq.2)) then
         call col3    (ta1,x(1,e),bm2(1,1,1,e),nxyz2)
         call invcol2 (ta1,ym2(1,1,1,e),nxyz2)
         call mxm     (ixtm12,lx1,ta1,lx2,ta2,ly2)
         call mxm     (ta2,lx1,iym12,ly2,ta1,ly1)
         call add2    (dtx(1,e),ta1,nxyz1)
       endif

      endif

      enddo
C
#ifdef TIMER
      tcdtp=tcdtp+(dnekclock()-etime1)
#endif
      return
      end
!---------------------------------------------------------------------- 

