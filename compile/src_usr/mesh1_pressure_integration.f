!======================================================================
!     Author: Prabal Negi
!     Description: Integrate pressure terms using GLL grid
!     Routines: cM1dabdtp()          : E operator on M1
!              crhodabdtp()          : E operator with variable density
!               opgradtM1()          : p(D*v)  on M1  
!               opdivM1()            : q*(D*v) on M1               
!               cdtM1p()             : DT*p on M1
!               MultDM1()            : D*x  on M1 
!               rho_opdivM1()        : q*\rho*(D*v) on M1               
!               rho_MultDM1()        : \rho*D*x  on M1 
!               rho_opdiv()          : q*\rho*(D*v) on M2
!               rho_MultD()          : \rho*D*x  on M2
!               cdab_full_dtp()      : E operator with full mass matrix inversion
!
!======================================================================       
C-----------------------------------------------------------------------
      subroutine cM1dabdtp (ap,wp,h1,h2,h2inv,intype)

C     INTYPE= 0  Compute the matrix-vector product    DA(-1)DT*p
C     INTYPE= 1  Compute the matrix-vector product    D(B/DT)(-1)DT*p
C     INTYPE=-1  Compute the matrix-vector product    D(A+B/DT)(-1)DT*p

      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'         ! ifanls
!      include 'TOTAL'
      REAL           AP    (LX2,LY2,LZ2,1)
      REAL           WP    (LX2,LY2,LZ2,1)
      REAL           H1    (LX1,LY1,LZ1,1)
      REAL           H2    (LX1,LY1,LZ1,1)
      REAL           H2INV (LX1,LY1,LZ1,1)

      real ta1,ta2,ta3,tb1,tb2,tb3
      COMMON /SCRNS/ TA1 (LX1,LY1,LZ1,LELV)
     $ ,             TA2 (LX1,LY1,LZ1,LELV)
     $ ,             TA3 (LX1,LY1,LZ1,LELV)
     $ ,             TB1 (LX1,LY1,LZ1,LELV)
     $ ,             TB2 (LX1,LY1,LZ1,LELV)
     $ ,             TB3 (LX1,LY1,LZ1,LELV)

      integer intype
      real tolhin,dtbdi

      call opgradtM1(ta1,ta2,ta3,wp)
      if ((intype.eq.0).or.(intype.eq.-1)) then
         tolhin=tolhs
         call ophinv (tb1,tb2,tb3,ta1,ta2,ta3,h1,h2,tolhin,nmxv)
      else
         if (ifanls) then
            dtbdi = dt/bd(1)   ! scale by dt*backwd-diff coefficient
            CALL OPBINV1(TB1,TB2,TB3,TA1,TA2,TA3,dtbdi)
         else
            CALL OPBINV (TB1,TB2,TB3,TA1,TA2,TA3,H2INV)
         endif
      endif
      call opdivM1 (ap,tb1,tb2,tb3)

      return
      end
C
C-----------------------------------------------------------------------
      subroutine crhodabdtp (ap,wp,h1,h2,h2inv,intype)

C     INTYPE= 0  Compute the matrix-vector product    DA(-1)DT*p
C     INTYPE= 1  Compute the matrix-vector product    D(B/DT)(-1)DT*p
C     INTYPE=-1  Compute the matrix-vector product    D(A+B/DT)(-1)DT*p

      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'         ! ifanls
!      include 'TOTAL'
      REAL           AP    (LX2,LY2,LZ2,1)
      REAL           WP    (LX2,LY2,LZ2,1)
      REAL           H1    (LX1,LY1,LZ1,1)
      REAL           H2    (LX1,LY1,LZ1,1)
      REAL           H2INV (LX1,LY1,LZ1,1)

      real ta1,ta2,ta3,tb1,tb2,tb3
      COMMON /SCRNS/ TA1 (LX1,LY1,LZ1,LELV)
     $ ,             TA2 (LX1,LY1,LZ1,LELV)
     $ ,             TA3 (LX1,LY1,LZ1,LELV)
     $ ,             TB1 (LX1,LY1,LZ1,LELV)
     $ ,             TB2 (LX1,LY1,LZ1,LELV)
     $ ,             TB3 (LX1,LY1,LZ1,LELV)

      integer intype
      real tolhin,dtbdi

      call opgradt(ta1,ta2,ta3,wp)
      if ((intype.eq.0).or.(intype.eq.-1)) then
         tolhin=tolhs
         call ophinv (tb1,tb2,tb3,ta1,ta2,ta3,h1,h2,tolhin,nmxv)
      else
         if (ifanls) then
            dtbdi = dt/bd(1)   ! scale by dt*backwd-diff coefficient
            CALL OPBINV1(TB1,TB2,TB3,TA1,TA2,TA3,dtbdi)
         else
            CALL OPBINV (TB1,TB2,TB3,TA1,TA2,TA3,H2INV)
         endif
      endif
      call rho_opdiv (ap,tb1,tb2,tb3)

      return
      end
C
C-----------------------------------------------------------------------
      subroutine cdab_full_dtp (ap,wp,h1,h2,h2inv,intype)

C     INTYPE= 0  Compute the matrix-vector product    DA(-1)DT*p
C     INTYPE= 1  Compute the matrix-vector product    D(B/DT)(-1)DT*p
C     INTYPE=-1  Compute the matrix-vector product    D(A+B/DT)(-1)DT*p

      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'         ! ifanls
      include 'SOLN'          ! V?MASK

      REAL           AP    (LX2,LY2,LZ2,1)
      REAL           WP    (LX2,LY2,LZ2,1)
      REAL           H1    (LX1,LY1,LZ1,1)
      REAL           H2    (LX1,LY1,LZ1,1)
      REAL           H2INV (LX1,LY1,LZ1,1)

      real ta1,ta2,ta3,tb1,tb2,tb3
      COMMON /SCRNS/ TA1 (LX1,LY1,LZ1,LELV)
     $ ,             TA2 (LX1,LY1,LZ1,LELV)
     $ ,             TA3 (LX1,LY1,LZ1,LELV)
     $ ,             TB1 (LX1,LY1,LZ1,LELV)
     $ ,             TB2 (LX1,LY1,LZ1,LELV)
     $ ,             TB3 (LX1,LY1,LZ1,LELV)

      integer intype
      real tolhin,dtbdi

      logical iffullmass
      integer maxiter

      iffullmass = .false.

      call opgradt(ta1,ta2,ta3,wp)
      if ((intype.eq.0).or.(intype.eq.-1)) then
        tolhin=tolhs
        call ophinv (tb1,tb2,tb3,ta1,ta2,ta3,h1,h2,tolhin,nmxv)
      else
        if (ifanls) then
          dtbdi = dt/bd(1)   ! scale by dt*backwd-diff coefficient
          CALL OPBINV1(TB1,TB2,TB3,TA1,TA2,TA3,dtbdi)
        else
          if (iffullmass) then
            maxiter = 200
            call opcopy(tb1,tb2,tb3,ta1,ta2,ta3)
!            call op_gmres(tb1,tb2,tb3,h2,maxiter)
          else
            CALL OPBINV (TB1,TB2,TB3,TA1,TA2,TA3,H2INV)
          endif
        endif
      endif
      call opdiv (ap,tb1,tb2,tb3)

      return
      end
C
C-----------------------------------------------------------------------

      subroutine opdivM1(outfld,inpx,inpy,inpz)

C     Compute OUTFLD = SUMi Di*INPi, 
C     the divergence of the vector field (INPX,INPY,INPZ)


      implicit none

      include 'SIZE'
      include 'GEOM'

      real outfld (lx2,ly2,lz2,1)
      real inpx   (lx1,ly1,lz1,1)
      real inpy   (lx1,ly1,lz1,1)
      real inpz   (lx1,ly1,lz1,1)

      real work
      common /ctmp0/ work (lx2,ly2,lz2,lelv)
      
      integer iflg,ntot2

      iflg = 1

      ntot2 = lx2*ly2*lz2*nelv
      call multdM1 (work,inpx,rxm1,sxm1,txm1,1,iflg)
      call copy  (outfld,work,ntot2)
      call multdM1 (work,inpy,rym1,sym1,tym1,2,iflg)
      call add2  (outfld,work,ntot2)
      if (ldim.eq.3) then
         call multdM1 (work,inpz,rzm1,szm1,tzm1,3,iflg)
         call add2  (outfld,work,ntot2)
      endif
C
      return
      end
C
c-----------------------------------------------------------------------
      subroutine opgradtM1(outx,outy,outz,inpfld)

C     Compute DTx, DTy, DTz of an input field INPFLD
C     Evaluated on Mesh 1        

      implicit none

      include 'SIZE'
      include 'GEOM'

      real outx   (lx1,ly1,lz1,1)
      real outy   (lx1,ly1,lz1,1)
      real outz   (lx1,ly1,lz1,1)
      real inpfld (lx2,ly2,lz2,1)
C
      call cdtM1p (outx,inpfld,rxm1,sxm1,txm1,1)
      call cdtM1p (outy,inpfld,rym1,sym1,tym1,2)
      if (ldim.eq.3) 
     $   call cdtM1p (outz,inpfld,rzm1,szm1,tzm1,3)
C
      return
      end
c-----------------------------------------------------------------------

      subroutine cdtM1p (dtx,x,rm1,sm1,tm1,isd)

!     Compute DT*X (entire field)
!     Evaluated on the M1 Mesh.        

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
      real rm1  (lx1*ly1*lz1,lelv)
      real sm1  (lx1*ly1*lz1,lelv)
      real tm1  (lx1*ly1*lz1,lelv)

      real wx,ta1,ta2,ta3
      common /ctmp1/ wx  (lx1*ly1*lz1)
     $ ,             ta1 (lx1*ly1*lz1)
     $ ,             ta2 (lx1*ly1*lz1)
     $ ,             ta3 (lx1*ly1,lz1)

      REAL           DUAX(LX1)

      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV

      integer isd
      integer e
      integer nxyz1,nxyz2,nyz1,nyz2,nxy1
      integer n1,n2,i1,i2,iz
C
#ifdef TIMER
      if (icalld.eq.0) tcdtp=0.0
      icalld=icalld+1
      ncdtp=icalld
      etime1=dnekclock()
#endif

      nxyz1 = lx1*ly1*lz1
      nxyz2 = lx2*ly2*lz2
      nyz1  = ly1*lz1
      nyz2  = ly2*lz2
      nxy1  = lx1*ly1

      n1    = lx1*ly1
      n2    = lx1*ly2

      do e=1,nelv

C       Use the appropriate derivative- and interpolation operator in 
C       the y-direction (= radial direction if axisymmetric).
        if (ifaxis) then
          if (nio.eq.0) write(6,*) 'CDTM1p not implemented for ifaxis.'
          call exitt
        endif

        if(ifsplit) then
          if (nio.eq.0) write(6,*) 'CDTM1p not implemented for ifsplit.'
          call exitt
        else
!         prabal. Interpolate x to Mesh 1
          if (ndim.eq.2) then
            call mxm (ixm21,lx1,x(1,e),lx2,ta1,nyz2)
            call mxm (ta1,lx1,iytm21,ly2,wx,ly1)
          else
            call mxm(ixm21,lx1,x(1,e),lx2,ta1,nyz2)
            i1=1
            i2=1
            do iz=1,lz2
              call mxm (ta1(i1),lx1,iytm21,ly2,ta2(i2),ly1) 
              i1=i1+(lx1*ly2)
              i2=i2+(lx1*ly1)
            enddo
            call mxm (ta2,nxy1,iztm21,lz2,wx,lz1)
          endif  

!         collocate with weights
!         Jacobian goes away due to inverse jacobian of the dx/dr etc. 
          call col2(wx,w3m1,nxyz1)
        endif
C
        if (ldim.eq.2) then
          if (.not.ifdfrm(e) .and. ifalgn(e)) then

             if (      ifrsxy(e).and.isd.eq.1  .or. 
     $            .not.ifrsxy(e).and.isd.eq.2) then

!               prabal. 
                call col3 (ta1,wx,rm1(1,e),nxyz1)
                call mxm  (dxtm1,lx1,ta1,lx1,dtx(1,e),nyz1)
             else
!               prabal   
                call col3 (ta1,wx,sm1(1,e),nxyz1)
                call mxm  (ta1,lx1,dym1,ly1,dtx(1,e),ly1)

             endif
          else

!            wx*dr/dx*dv/dr
             call col3 (ta1,wx,rm1(1,e),nxyz1)
             call mxm  (dxtm1,lx1,ta1,lx1,dtx(1,e),nyz1)

!            wx*ds/dx*dv/ds
             call col3 (ta1,wx,sm1(1,e),nxyz1)
             call mxm (ta1,lx1,dym1,ly1,ta2,ly1) 

             call add2 (dtx(1,e),ta2,nxyz1)

          endif

        else
          if (ifsplit) then

            if (nio.eq.0) 
     $        write(6,*) 'CDTM1p not implemented for ifsplit.'
            call exitt
!
          else
!            wx*dr/dx*dv/dr
             call col3 (ta1,wx,rm1(1,e),nxyz1)
             call mxm  (dxtm1,lx1,ta1,lx1,dtx(1,e),nyz1)

!            wx*ds/dx*dv/ds
             call col3 (ta1,wx,sm1(1,e),nxyz1)
             i1=1
             i2=1
             do iz=1,lz1
               call mxm (ta1(i1),lx1,dym1,ly1,ta2(i2),ly1) 
               i1=i1+(lx1*ly1)
               i2=i2+(lx1*ly1)
             enddo
             call add2 (dtx(1,e),ta2,nxyz1)

!            wx*dt/dx*dv/dt
             call col3 (ta1,wx,tm1(1,e),nxyz1)
             call mxm  (ta1,nxy1,dzm1,lz1,ta2,lz1)
             call add2 (dtx(1,e),ta2,nxyz1)
          endif

        endif         
C
C     If axisymmetric, add an extra diagonal term in the radial 
C     direction (only if solving the momentum equations and ISD=2)
C     NOTE: lz1=lz2=1

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
      subroutine multdM1 (dx,x,rm1,sm1,tm1,isd,iflg)

!     Compute D*X (on Mesh 1)
!     X    : input variable, defined on M1
!     DX   : output variable, defined on M2
!     Integration done on the M1 mesh        
!     RM1 : RXM1, RYM1 or RZM1
!     SM1 : SXM1, SYM1 or SZM1
!     TM1 : TXM1, TYM1 or TZM1
!     ISD : spatial direction (x=1,y=2,z=3)
!     IFLG: OPGRAD (iflg=0) or OPDIV (iflg=1)

      
      implicit none

      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'IXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
      include 'ESOLV'

      real           dx   (lx2*ly2*lz2,lelv)
      real           x    (lx1*ly1*lz1,lelv)
      real           rm1  (lx1*ly1*lz1,lelv)
      real           sm1  (lx1*ly1*lz1,lelv)
      real           tm1  (lx1*ly1*lz1,lelv)

      real           wk1  (lx1*ly1*lz1)
      real           wk2  (lx1*ly1*lz1)

      integer isd,iflg

      real ta1,ta2,ta3
      common /ctmp1/ ta1 (lx1*ly1*lz1)
     $ ,             ta2 (lx1*ly1*lz1)
     $ ,             ta3 (lx1*ly1*lz1)

      real           duax(lx1)

      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      logical ifdfrm, iffast, ifh2, ifsolv
      include 'CTIMER'

      integer e,i1,i2,iz
      integer nxy1,nyz1,nxy2,nxyz1,nxyz2,n1,n2

C
#ifdef TIMER
      if (icalld.eq.0) tmltd=0.0
      icalld=icalld+1
      nmltd=icalld
      etime1=dnekclock()
#endif

      nxy1  = lx1*ly1
      nyz1  = ly1*lz1
      nxy2  = lx2*ly2
      nxyz1 = lx1*ly1*lz1
      nxyz2 = lx2*ly2*lz2

      n1    = lx1*ly1
      n2    = lx2*ly2

      do e=1,nelv

c        Use the appropriate derivative- and interpolation operator in 
c        the y-direction (= radial direction if axisymmetric).
         if (ifaxis) then
           if (nio.eq.0) write(6,*) 
     $       'MULTDM1 not implemented for ifaxis.'
           call exitt
          
!            ly12   = ly1*ly2
!            if (ifrzer(e)) then
!               call copy (iytm12,iatm12,ly12)
!               call copy (dytm12,datm12,ly12)
!               call copy (w3m2,w2am2,nxyz2)
!            else
!               call copy (iytm12,ictm12,ly12)
!               call copy (dytm12,dctm12,ly12)
!               call copy (w3m2,w2cm2,nxyz2)
!            endif
         endif

         if (ldim.eq.2) then
            if (.not.ifdfrm(e) .and. ifalgn(e)) then
c
               if (      ifrsxy(e).and.isd.eq.1  .or. 
     $              .not.ifrsxy(e).and.isd.eq.2) then
                  call mxm     (dxm1,lx1,x(1,e),lx1,wk1,nyz1)
!                  call mxm     (ta1,lx2,iytm12,ly1,dx(1,e),ly2)
                  call col2    (wk1,rm1(1,e),nxyz1)

!                  call mxm     (dxm12,lx2,x(1,e),lx1,ta1,nyz1)
!                  call mxm     (ta1,lx2,iytm12,ly1,dx(1,e),ly2)
!                  call col2    (dx(1,e),rm2(1,e),nxyz2)
               else
!                  call mxm     (ixm12,lx2,x(1,e),lx1,ta1,nyz1)
                  call mxm     (x(1,e),lx1,dytm1,ly1,wk1,ly1)
                  call col2    (wk1,sm1(1,e),nxyz1)

!                  call mxm     (ixm12,lx2,x(1,e),lx1,ta1,nyz1)
!                  call mxm     (ta1,lx2,dytm12,ly1,dx(1,e),ly2)
!                  call col2    (dx(1,e),sm2(1,e),nxyz2)
               endif
            else
               call mxm     (dxm1,lx1,x(1,e),lx1,wk1,nyz1)
!               call mxm     (ta1,lx2,iytm12,ly1,dx(1,e),ly2)
               call col2    (wk1,rm1(1,e),nxyz1)
!               call mxm     (ixm12,lx2,x(1,e),lx1,ta1,nyz1)
               call mxm     (x(1,e),lx1,dytm1,ly1,ta3,ly1)
               call addcol3 (wk1,ta3,sm1(1,e),nxyz1)

!               call mxm     (dxm12,lx2,x(1,e),lx1,ta1,nyz1)
!               call mxm     (ta1,lx2,iytm12,ly1,dx(1,e),ly2)
!               call col2    (dx(1,e),rm2(1,e),nxyz2)
!               call mxm     (ixm12,lx2,x(1,e),lx1,ta1,nyz1)
!               call mxm     (ta1,lx2,dytm12,ly1,ta3,ly2)
!               call addcol3 (dx(1,e),ta3,sm2(1,e),nxyz2)
            endif

         else  ! 3D

!             if (nio.eq.0) write(6,*) 
!     $         'MULTDM1 not implemented for 3D'
!             call exitt

             call mxm  (dxm1,lx1,x(1,e),lx1,ta1,nyz1)
             call col3 (wk1,ta1,rm1(1,e),nxyz1)
!
             call copy(ta3,x(1,e),nxyz1)
             i1=1
             i2=1
             do iz=1,lz1
               call mxm (ta3(i1),lx1,dytm1,ly1,ta2(i2),ly1)
               i1=i1+n1
               i2=i2+n1
             enddo
             call addcol3 (wk1,ta2,sm1(1,e),nxyz1)
!
             call copy(ta1,x(1,e),nxyz1)
             call mxm (ta1,nxy1,dztm1,lz1,ta3,lz1)
             call addcol3 (wk1,ta3,tm1(1,e),nxyz1)
         endif
C
C        Collocate with the weights on the pressure mesh


       if(ifsplit) then
!         call col2   (dx(1,e),bm1(1,1,1,e),nxyz1)
!         call invcol2(dx(1,e),jacm1(1,1,1,e),nxyz1)
       else
!        collocate with weights on Mesh 1          
         if (.not.ifaxis) call col2 (wk1,w3m1,nxyz1)
!        Using pressure test function on Mesh 1
!        integrate to get result on Mesh 2
         if (if3d) then
           call mxm(ixtm21,lx2,wk1,lx1,ta1,nyz1)
           i1=1
           i2=1
           do iz=1,lz1
             call mxm (ta1(i1),lx2,iym21,ly1,ta2(i2),ly2)
             i1=i1+(lx2*ly1)
             i2=i2+n2
           enddo
           call mxm (ta2,nxy2,izm21,lz1,dx(1,e),lz2)
         else  
           call mxm(ixtm21,lx2,wk1,lx1,ta1,lx1)
           call mxm(ta1,lx2,iym21,lx1,dx(1,e),lx2)
         endif  

!         if (.not.ifaxis) call col2 (dx(1,e),w3m2,nxyz2)
!         if (ifaxis) then
!             if (ifrzer(e)) then
!                 call col2    (dx(1,e),bm2(1,1,1,e),nxyz2)
!                 call invcol2 (dx(1,e),jacm2(1,1,1,e),nxyz2)
!             else
!                 call col2    (dx(1,e),w3m2,nxyz2)
!                 call col2    (dx(1,e),ym2(1,1,1,e),nxyz2)
!             endif
!         endif
       endif

c        If axisymmetric, add an extra diagonal term in the radial 
c        direction (ISD=2).
c        NOTE: lz1=lz2=1

!      if(ifsplit) then
!
!       if (ifaxis.and.(isd.eq.2).and.iflg.eq.1) then
!        call copy    (ta3,x(1,e),nxyz1)
!        if (ifrzer(e)) then
!           call rzero(ta3, lx1)
!           call mxm  (x(1,e),lx1,datm1,ly1,duax,1)
!           call copy (ta3,duax,lx1)
!        endif
!        call col2    (ta3,baxm1(1,1,1,e),nxyz1)
!        call add2    (dx(1,e),ta3,nxyz2)
!       endif
!
!      else
!
!       if (ifaxis.and.(isd.eq.2)) then
!            call mxm     (ixm12,lx2,x(1,e),lx1,ta1,ly1)
!            call mxm     (ta1,lx2,iytm12,ly1,ta2,ly2)
!            call col3    (ta3,bm2(1,1,1,e),ta2,nxyz2)
!            call invcol2 (ta3,ym2(1,1,1,e),nxyz2)
!            call add2    (dx(1,e),ta3,nxyz2)
!       endif
!
!      endif

      enddo
C
#ifdef TIMER
      tmltd=tmltd+(dnekclock()-etime1)
#endif
      return
      END
c-----------------------------------------------------------------------

      subroutine rho_opdivM1(outfld,inpx,inpy,inpz)
C
C     Compute OUTFLD = SUMi Di*INPi, 
C     the divergence of the vector field (INPX,INPY,INPZ)
C
C---------------------------------------------------------------------

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'SOLN'    ! \rho
      include 'TSTEP'   ! ifield

      real outfld (lx2,ly2,lz2,1)
      real inpx   (lx1,ly1,lz1,1)
      real inpy   (lx1,ly1,lz1,1)
      real inpz   (lx1,ly1,lz1,1)

      real work
      common /ctmp0/ work (lx2,ly2,lz2,lelv)
      
      integer iflg,ntot2

      iflg = 1

      ntot2 = lx2*ly2*lz2*nelv
      call rho_multdM1 (work,inpx,vtrans(1,1,1,1,ifield)
     $                  ,rxm1,sxm1,txm1,1,iflg)
      call copy  (outfld,work,ntot2)
      call rho_multdM1 (work,inpy,vtrans(1,1,1,1,ifield)
     $                  ,rym1,sym1,tym1,2,iflg)
      call add2  (outfld,work,ntot2)
      if (ldim.eq.3) then
        call rho_multdM1 (work,inpz,vtrans(1,1,1,1,ifield)
     $                    ,rzm1,szm1,tzm1,3,iflg)
         call add2  (outfld,work,ntot2)
      endif
C
      return
      end
C
c-----------------------------------------------------------------------
      subroutine rho_multdM1 (dx,x,rho,rm1,sm1,tm1,isd,iflg)

!     Compute D*X (on Mesh 1)
!     X    : input variable, defined on M1
!     DX   : output variable, defined on M2
!     Integration done on the M1 mesh        
!     RM1 : RXM1, RYM1 or RZM1
!     SM1 : SXM1, SYM1 or SZM1
!     TM1 : TXM1, TYM1 or TZM1
!     ISD : spatial direction (x=1,y=2,z=3)
!     IFLG: OPGRAD (iflg=0) or OPDIV (iflg=1)

      
      implicit none

      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'IXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
      include 'ESOLV'

      real           dx   (lx2*ly2*lz2,lelv)
      real           x    (lx1*ly1*lz1,lelv)
      real           rho  (lx1*ly1*lz1,lelv)
      real           rm1  (lx1*ly1*lz1,lelv)
      real           sm1  (lx1*ly1*lz1,lelv)
      real           tm1  (lx1*ly1*lz1,lelv)

      real           wk1  (lx1*ly1*lz1)
      real           wk2  (lx1*ly1*lz1)

      integer isd,iflg

      real ta1,ta2,ta3
      common /ctmp1/ ta1 (lx1*ly1*lz1)
     $ ,             ta2 (lx1*ly1*lz1)
     $ ,             ta3 (lx1*ly1*lz1)

      real           duax(lx1)

      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      logical ifdfrm, iffast, ifh2, ifsolv
      include 'CTIMER'

      integer e,i1,i2,iz
      integer nxy1,nyz1,nxy2,nxyz1,nxyz2,n1,n2

C
#ifdef TIMER
      if (icalld.eq.0) tmltd=0.0
      icalld=icalld+1
      nmltd=icalld
      etime1=dnekclock()
#endif

      nxy1  = lx1*ly1
      nyz1  = ly1*lz1
      nxy2  = lx2*ly2
      nxyz1 = lx1*ly1*lz1
      nxyz2 = lx2*ly2*lz2

      n1    = lx1*ly1
      n2    = lx2*ly2

      do e=1,nelv
c       Use the appropriate derivative- and interpolation operator in 
c       the y-direction (= radial direction if axisymmetric).
        if (ifaxis) then
          if (nio.eq.0) write(6,*) 
     $      'MULTDM1 not implemented for ifaxis.'
          call exitt
         
        endif

        if (ldim.eq.2) then
          if (.not.ifdfrm(e) .and. ifalgn(e)) then
            if (      ifrsxy(e).and.isd.eq.1  .or. 
     $           .not.ifrsxy(e).and.isd.eq.2) then
               call mxm     (dxm1,lx1,x(1,e),lx1,wk1,nyz1)
               call col2    (wk1,rm1(1,e),nxyz1)
            else
               call mxm     (x(1,e),lx1,dytm1,ly1,wk1,ly1)
               call col2    (wk1,sm1(1,e),nxyz1)
            endif
          else
            call mxm     (dxm1,lx1,x(1,e),lx1,wk1,nyz1)
            call col2    (wk1,rm1(1,e),nxyz1)

            call mxm     (x(1,e),lx1,dytm1,ly1,ta3,ly1)
            call addcol3 (wk1,ta3,sm1(1,e),nxyz1)
          endif

        else  ! 3D
          call mxm  (dxm1,lx1,x(1,e),lx1,ta1,nyz1)
          call col3 (wk1,ta1,rm1(1,e),nxyz1)

          call copy(ta3,x(1,e),nxyz1)
          i1=1
          i2=1
          do iz=1,lz1
            call mxm (ta3(i1),lx1,dytm1,ly1,ta2(i2),ly1)
            i1=i1+n1
            i2=i2+n1
          enddo
          call addcol3 (wk1,ta2,sm1(1,e),nxyz1)
!
          call copy(ta1,x(1,e),nxyz1)
          call mxm (ta1,nxy1,dztm1,lz1,ta3,lz1)
          call addcol3 (wk1,ta3,tm1(1,e),nxyz1)
        endif

C       Collocate with the weights on the pressure mesh

        if(ifsplit) then
!          call col2   (dx(1,e),bm1(1,1,1,e),nxyz1)
!          call invcol2(dx(1,e),jacm1(1,1,1,e),nxyz1)
        else
!         collocate with weights on Mesh 1          
          if (.not.ifaxis) call col2 (wk1,w3m1,nxyz1)
!         Collocate with density            
          call col2(wk1,rho(1,e),nxyz1)   
!         Using pressure test function on Mesh 1
!         integrate to get result on Mesh 2
          if (if3d) then
            call mxm(ixtm21,lx2,wk1,lx1,ta1,nyz1)
            i1=1
            i2=1
            do iz=1,lz1
              call mxm (ta1(i1),lx2,iym21,ly1,ta2(i2),ly2)
              i1=i1+(lx2*ly1)
              i2=i2+n2
            enddo
            call mxm (ta2,nxy2,izm21,lz1,dx(1,e),lz2)
          else  
            call mxm(ixtm21,lx2,wk1,lx1,ta1,lx1)
            call mxm(ta1,lx2,iym21,lx1,dx(1,e),lx2)
          endif  
        endif      ! ifsplit
      enddo

#ifdef TIMER
      tmltd=tmltd+(dnekclock()-etime1)
#endif
      return
      END
c-----------------------------------------------------------------------

      subroutine rho_opdiv(outfld,inpx,inpy,inpz)

C     Compute OUTFLD = SUMi Di*INPi, 
C     the divergence of the vector field (INPX,INPY,INPZ)

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'SOLN'          ! \rho

      real outfld (lx2,ly2,lz2,1)
      real inpx   (lx1,ly1,lz1,1)
      real inpy   (lx1,ly1,lz1,1)
      real inpz   (lx1,ly1,lz1,1)

      real work
      common /ctmp0/ work (lx2,ly2,lz2,lelv)

!     Should put this in some common block
      real rho2(lx2*ly2*lz2,lelv)         ! density on Mesh 2
     
      integer iflg,ntot2,e
      integer ifld

      ifld = 1
      do e=1,nelv
!       Uses common block /ctmp00/           
        call map12(rho2(1,e),vtrans(1,1,1,e,ifld),e)
      enddo      

      iflg = 1

      ntot2 = lx2*ly2*lz2*nelv
      call rho_multd (work,inpx,rho2,rxm2,sxm2,txm2,1,iflg)
      call copy  (outfld,work,ntot2)
      call rho_multd (work,inpy,rho2,rym2,sym2,tym2,2,iflg)
      call add2  (outfld,work,ntot2)
      if (ldim.eq.3) then
        call rho_multd (work,inpz,rho2,rzm2,szm2,tzm2,3,iflg)
        call add2  (outfld,work,ntot2)
      endif
C
      return
      end
C
c-----------------------------------------------------------------------

      subroutine rho_multd (dx,x,rho2,rm2,sm2,tm2,isd,iflg)

C     Compute D*X
C     X    : input variable, defined on M1
C     DX   : output variable, defined on M2 (note: D is rectangular)   
C     RM2 : RXM2, RYM2 or RZM2
C     SM2 : SXM2, SYM2 or SZM2
C     TM2 : TXM2, TYM2 or TZM2
C     ISD : spatial direction (x=1,y=2,z=3)
C     IFLG: OPGRAD (iflg=0) or OPDIV (iflg=1)

      implicit none

      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'IXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
      include 'ESOLV'

      real           dx   (lx2*ly2*lz2,lelv)
      real           x    (lx1*ly1*lz1,lelv)
      real           rho2 (lx2*ly2*lz2,lelv)    ! Density on mesh 2
      real           rm2  (lx2*ly2*lz2,lelv)
      real           sm2  (lx2*ly2*lz2,lelv)
      real           tm2  (lx2*ly2*lz2,lelv)

      real ta1,ta2,ta3
      common /ctmp1/ ta1 (lx1*ly1*lz1)
     $ ,             ta2 (lx1*ly1*lz1)
     $ ,             ta3 (lx1*ly1*lz1)

      real           duax(lx1)

      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      logical ifdfrm, iffast, ifh2, ifsolv
      include 'CTIMER'

      integer e,i1,i2,iz
      integer nxy1,nyz1,nxy2,nxyz1,nxyz2,n1,n2,ly12

      integer isd,iflg


#ifdef TIMER
      if (icalld.eq.0) tmltd=0.0
      icalld=icalld+1
      nmltd=icalld
      etime1=dnekclock()
#endif

      nyz1  = ly1*lz1
      nxy2  = lx2*ly2
      nxyz1 = lx1*ly1*lz1
      nxyz2 = lx2*ly2*lz2

      n1    = lx2*ly1
      n2    = lx2*ly2

      do e=1,nelv
c       Use the appropriate derivative- and interpolation operator in 
c       the y-direction (= radial direction if axisymmetric).
        if (ifaxis) then
           ly12   = ly1*ly2
           if (ifrzer(e)) then
              call copy (iytm12,iatm12,ly12)
              call copy (dytm12,datm12,ly12)
              call copy (w3m2,w2am2,nxyz2)
           else
              call copy (iytm12,ictm12,ly12)
              call copy (dytm12,dctm12,ly12)
              call copy (w3m2,w2cm2,nxyz2)
           endif
        endif

        if (ldim.eq.2) then
           if (.not.ifdfrm(e) .and. ifalgn(e)) then
c
              if (      ifrsxy(e).and.isd.eq.1  .or. 
     $             .not.ifrsxy(e).and.isd.eq.2) then
                 call mxm     (dxm12,lx2,x(1,e),lx1,ta1,nyz1)
                 call mxm     (ta1,lx2,iytm12,ly1,dx(1,e),ly2)
                 call col2    (dx(1,e),rm2(1,e),nxyz2)
              else
                 call mxm     (ixm12,lx2,x(1,e),lx1,ta1,nyz1)
                 call mxm     (ta1,lx2,dytm12,ly1,dx(1,e),ly2)
                 call col2    (dx(1,e),sm2(1,e),nxyz2)
              endif
           else
              call mxm     (dxm12,lx2,x(1,e),lx1,ta1,nyz1)
              call mxm     (ta1,lx2,iytm12,ly1,dx(1,e),ly2)
              call col2    (dx(1,e),rm2(1,e),nxyz2)
              call mxm     (ixm12,lx2,x(1,e),lx1,ta1,nyz1)
              call mxm     (ta1,lx2,dytm12,ly1,ta3,ly2)
              call addcol3 (dx(1,e),ta3,sm2(1,e),nxyz2)
           endif

        else  ! 3D

c          if (ifsplit) then
c
c            call mxm  (dxm12,lx2,x(1,e),lx1,dx(1,e),nyz1)
c            call col2 (dx(1,e),rm2(1,e),nxyz2)
c            i1=1
c            i2=1
c            do iz=1,lz1
c               call mxm (x(1,e),lx2,dytm12,ly1,ta1(i2),ly2)
c               i1=i1+n1
c               i2=i2+n2
c            enddo
c            call addcol3 (dx(1,e),ta1,sm2(1,e),nxyz2)
c            call mxm (x(1,e),nxy2,dztm12,lz1,ta1,lz2)
c            call addcol3 (dx(1,e),ta1,tm2(1,e),nxyz2)

c          else ! PN - PN-2

            call mxm (dxm12,lx2,x(1,e),lx1,ta1,nyz1)
            i1=1
            i2=1
            do iz=1,lz1
              call mxm (ta1(i1),lx2,iytm12,ly1,ta2(i2),ly2)
              i1=i1+n1
              i2=i2+n2
            enddo
            call mxm  (ta2,nxy2,iztm12,lz1,dx(1,e),lz2)
            call col2 (dx(1,e),rm2(1,e),nxyz2)

            call mxm  (ixm12,lx2,x(1,e),lx1,ta3,nyz1) ! reuse ta3 below
            i1=1
            i2=1
            do iz=1,lz1
              call mxm (ta3(i1),lx2,dytm12,ly1,ta2(i2),ly2)
              i1=i1+n1
              i2=i2+n2
            enddo
            call mxm     (ta2,nxy2,iztm12,lz1,ta1,lz2)
            call addcol3 (dx(1,e),ta1,sm2(1,e),nxyz2)

c           call mxm (ixm12,lx2,x(1,e),lx1,ta1,nyz1) ! reuse ta3 from above
            i1=1
            i2=1
            do iz=1,lz1
              call mxm (ta3(i1),lx2,iytm12,ly1,ta2(i2),ly2)
              i1=i1+n1
              i2=i2+n2
            enddo
            call mxm (ta2,nxy2,dztm12,lz1,ta3,lz2)
            call addcol3 (dx(1,e),ta3,tm2(1,e),nxyz2)
c          endif
        endif

C       Collocate with the weights on the pressure mesh

        if(ifsplit) then
          call col2   (dx(1,e),bm1(1,1,1,e),nxyz1)
          call invcol2(dx(1,e),jacm1(1,1,1,e),nxyz1)
        else
          if (.not.ifaxis) call col2 (dx(1,e),w3m2,nxyz2)
          if (ifaxis) then
            if (ifrzer(e)) then
              call col2    (dx(1,e),bm2(1,1,1,e),nxyz2)
              call invcol2 (dx(1,e),jacm2(1,1,1,e),nxyz2)
            else
              call col2    (dx(1,e),w3m2,nxyz2)
              call col2    (dx(1,e),ym2(1,1,1,e),nxyz2)
            endif
          endif
        endif

c       If axisymmetric, add an extra diagonal term in the radial 
c       direction (ISD=2).
c       NOTE: lz1=lz2=1

        if(ifsplit) then
          if (ifaxis.and.(isd.eq.2).and.iflg.eq.1) then
           call copy    (ta3,x(1,e),nxyz1)
           if (ifrzer(e)) then
              call rzero(ta3, lx1)
              call mxm  (x(1,e),lx1,datm1,ly1,duax,1)
              call copy (ta3,duax,lx1)
           endif
           call col2    (ta3,baxm1(1,1,1,e),nxyz1)
           call add2    (dx(1,e),ta3,nxyz2)
          endif
        else
          if (ifaxis.and.(isd.eq.2)) then
            call mxm     (ixm12,lx2,x(1,e),lx1,ta1,ly1)
            call mxm     (ta1,lx2,iytm12,ly1,ta2,ly2)
            call col3    (ta3,bm2(1,1,1,e),ta2,nxyz2)
            call invcol2 (ta3,ym2(1,1,1,e),nxyz2)
            call add2    (dx(1,e),ta3,nxyz2)
          endif
        endif

!       Collocate with density            
        call col2(dx(1,e),rho2(1,e),nxyz2)

      enddo       ! i=1,nelv

#ifdef TIMER
      tmltd=tmltd+(dnekclock()-etime1)
#endif
      return
      END
c-----------------------------------------------------------------------





