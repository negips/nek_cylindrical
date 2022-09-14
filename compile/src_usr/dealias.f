!====================================================================== 
!     Author: Prabal Negi
!     Description: Routines for Dealiasing including the density field
!
!====================================================================== 
!====================================================================== 
      subroutine advab_density

C     Eulerian scheme, add convection term to forcing function 
C     at current time step.
C     Also taking Density into account in the dealiasing step


      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      integer ntot1

      ntot1 = lx1*ly1*lz1*nelv
      call convop_rho(ta1,vx,vtrans)
!     Removing density since it is included in the dealiasing 
!     but is multiplied again later.
      call invcol2(ta1,vtrans,ntot1)

      call convop_rho(ta2,vy,vtrans)
      call invcol2(ta2,vtrans,ntot1)      ! Removing density

      call subcol3 (bfx,bm1,ta1,ntot1)
      call subcol3 (bfy,bm1,ta2,ntot1)
      if (ldim.eq.2) then
         call rzero (ta3,ntot1)
      else
         call convop_rho(ta3,vz,vtrans)
         call invcol2(ta3,vtrans,ntot1)   ! Removing density

         call subcol3 (bfz,bm1,ta3,ntot1)
      endif

      return
      end subroutine advab_density
c-----------------------------------------------------------------------
      subroutine convop_rho(conv,fi,rho)
C
C     Compute the convective term CONV for a passive scalar field FI
C     using the skew-symmetric formulation.
C     The field variable FI is defined on mesh M1 (GLL) and
C     the velocity field is assumed given.
C
C     IMPORTANT NOTE: Use the scratch-arrays carefully!!!!!
C
C     The common-block SCRNS is used in CONV1 and CONV2.
C     The common-blocks CTMP0 and CTMP1 are also used as scratch-arrays
C     since there is no direct stiffness summation or Helmholtz-solves. 


      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'DEALIAS'
      include 'MASS'
      include 'TSTEP'
!      include 'TOTAL'
      include 'CTIMER'
C
C     Use the common blocks CTMP0 and CTMP1 as work space.
C
      real cmask1,cmask2 
      common /scrch/  cmask1 (lx1,ly1,lz1,lelv)
     $ ,              cmask2 (lx1,ly1,lz1,lelv)

      real mfi,dmfi,mdmfi
      common /ctmp1/  mfi    (lx1,ly1,lz1,lelv)
     $ ,              dmfi   (lx1,ly1,lz1,lelv)
     $ ,              mdmfi  (lx1,ly1,lz1,lelv)
c
c     arrays in parameter list
c
      real    conv (lx1,ly1,lz1,1) 
      real    fi   (lx1,ly1,lz1,1)
      real    rho  (lx1,ly1,lz1,1)

      integer nxyz1,ntot1,ntotz,ntott

      if (nio.eq.0.and.loglevel.gt.2)
     $   write(6,*) 'convop_rho', ifield, ifdeal(ifield)

#ifdef TIMER
      if (icalld.eq.0) tadvc=0.0
      icalld=icalld+1
      nadvc=icalld
      etime1=dnekclock()
#endif
 
      nxyz1 = lx1*ly1*lz1
      ntot1 = lx1*ly1*lz1*nelv
      ntotz = lx1*ly1*lz1*nelfld(ifield)
      ntott = lx1*ly1*lz1*nelt
 
      call rzero  (conv,ntott)

!     Not implemented yet.
!     Somebody else can do it. 

!      if (ifdgfld(ifield)) then
!         call convect_dg (conv,fi,.false.,vxd,vyd,vzd,.true.)
!         goto 100
!      elseif (param(86).ne.0.0) then  ! skew-symmetric form
!         call convopo(conv,fi)
!         goto 100
!      endif


      if (.not. ifdeal(ifield)) then
         call conv1 (conv,fi)
         call col2(conv,rho,ntot1)   
      elseif (param(99).eq.2.or.param(99).eq.3) then
!         if (nio.eq.0) then
!           write(6,*) 'convop_rho not implemented for param(99)=',
!     $                 param(99)
!         endif
!         call exitt      
         call conv1d_rho(conv,fi,rho)
      elseif (param(99).eq.4) then
         if (ifield.eq.1) then   
           if (ifpert) then
             call convect_new_rho(conv,fi,.false.,vx,vy,vz,.false.,
     $                            rho,.false.)
           else
             call convect_new_rho(conv,fi,.false.,vxd,vyd,vzd,.true.,
     $                            rho,.false.)
           endif
         else
           if (ifpert) then
             call convect_new_rho(conv,fi,.false.,vx,vy,vz,.false.)
           else
             call convect_new_rho(conv,fi,.false.,vxd,vyd,vzd,.true.)
           endif
         endif        
         call invcol2     (conv,bm1,ntot1)  ! local mass inverse
      elseif (param(99).eq.5) then
         if (nio.eq.0) then
           write(6,*) 'convop_rho not implemented for param(99)=',
     $                 param(99)
         endif
         call exitt      

!         call convect_cons(conv,fi,.false.,vx,vy,vz,.false.)
!         call invcol2     (conv,bm1,ntot1)  ! local mass inverse
      else
         call conv1 (conv,fi)
      endif

 100  continue

#ifdef TIMER
      tadvc=tadvc+(dnekclock()-etime1)
#endif

      return
      end subroutine convop_rho
c---------------------------------------------------------------------- 
      subroutine convect_new_rho(bdu,u,ifuf,cx,cy,cz,ifcf,rho,ifrf)

C     J     - Interpolator from coarse to fine mesh
C     J^T   - Interpolator from fine to coarse mesh
C     C     - Convecting Field (on coarse mesh)
C     Bf    - Mass Matrix on fine mesh
C     \rho  - Density (on coarse mesh)
C     Compute dealiased form:  J^T (Bf * (J\rho) * (JC) .grad (Ju) ) w/ correct Jacobians

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
!      include 'TOTAL'

      real bdu(1),u(1),cx(1),cy(1),cz(1)
      logical ifuf,ifcf            ! u and/or c already on fine mesh?

      real rho(1)                  ! density 
      logical ifrf                 ! if \rho fine?

      real fx,fy,fz
      real ur,us,ut
      real tr,uf

      integer lxy,ltd
      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)
      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , tr(ltd,3),uf(ltd)


      real rhod(ltd)               ! dealiased density

      integer e,i,k
      integer iu,ic,ib,ir

      integer nxyz1,nxyzc,nxyzd,nxyzu,nxyzr


      call set_dealias_rx

      nxyz1 = lx1*ly1*lz1
      nxyzd = lxd*lyd*lzd

      nxyzu = nxyz1
      if (ifuf) nxyzu = nxyzd

      nxyzc = nxyz1
      if (ifcf) nxyzc = nxyzd

      nxyzr = nxyz1
      if (ifrf) nxyzr = nxyzd

      iu = 1    ! pointer to scalar field u
      ic = 1    ! pointer to vector field C
      ib = 1    ! pointer to scalar field Bdu
      ir = 1    ! pointer to scalar field rho

      do e=1,nelv

         if (ifcf) then

            call copy(tr(1,1),cx(ic),nxyzd)  ! already in rst form
            call copy(tr(1,2),cy(ic),nxyzd)
            if (if3d) call copy(tr(1,3),cz(ic),nxyzd)

         else  ! map coarse velocity to fine mesh (C-->F)

           call intp_rstd(fx,cx(ic),lx1,lxd,if3d,0) ! 0 --> forward
           call intp_rstd(fy,cy(ic),lx1,lxd,if3d,0) ! 0 --> forward
           if (if3d) call intp_rstd(fz,cz(ic),lx1,lxd,if3d,0) ! 0 --> forward

           if (if3d) then  ! Convert convector F to r-s-t coordinates

             do i=1,nxyzd
               tr(i,1)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)+rx(i,3,e)*fz(i)
               tr(i,2)=rx(i,4,e)*fx(i)+rx(i,5,e)*fy(i)+rx(i,6,e)*fz(i)
               tr(i,3)=rx(i,7,e)*fx(i)+rx(i,8,e)*fy(i)+rx(i,9,e)*fz(i)
             enddo

           else

             do i=1,nxyzd
               tr(i,1)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)
               tr(i,2)=rx(i,3,e)*fx(i)+rx(i,4,e)*fy(i)
             enddo

           endif

         endif

         if (ifuf) then
            call grad_rst(ur,us,ut,u(iu),lxd,if3d)
         else
            call intp_rstd(uf,u(iu),lx1,lxd,if3d,0) ! 0 --> forward
            call grad_rst(ur,us,ut,uf,lxd,if3d)
         endif

!        Interpolate density
         if (ifrf) then
           call copy(rhod,rho(ir),nxyzr)
         else 
           call intp_rstd(rhod,rho(ir),lx1,lxd,if3d,0) ! 0 --> forward
         endif   

         if (if3d) then
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
               uf(i) = (tr(i,1)*ur(i)+tr(i,2)*us(i)+tr(i,3)*ut(i))*
     $                  rhod(i) 
            enddo
         else
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
               uf(i) = (tr(i,1)*ur(i)+tr(i,2)*us(i))*rhod(i)
            enddo
         endif
         call intp_rstd(bdu(ib),uf,lx1,lxd,if3d,1) ! Project back to coarse

         ic = ic + nxyzc
         iu = iu + nxyzu
         ir = ir + nxyzr   
         ib = ib + nxyz1

      enddo

      return
      end subroutine convect_new_rho
!-----------------------------------------------------------------------
      subroutine conv1d_rho(dfi,fi,rho)

!     Compute \rho*U.D*FI 

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'GEOM'
      include 'INPUT'
      include 'DEALIAS'
      include 'TSTEP'
!      include 'TOTAL'

      real           dfi (lx1,ly1,lz1,1) 
      real           fi  (lx1,ly1,lz1,1)
      real           rho (lx1,ly1,lz1,1)

      real ta1,dfid,ta1d
      common /ctmp0/ ta1  (lx1,ly1,lz1,lelv)
     $             , dfid (lxd,lyd,lzd,lelv) 
     $             , ta1d (lxd,lyd,lzd,lelv) 

      real rhod(lxd,lyd,lzd,lelv)   ! In principle I should use a scratch array

      integer icalld
      save icalld
      data icalld /0/

      integer ntotd

      ntotd = lxd*lyd*lzd*nelv

c     interpolate ta1 and vx onto larger mesh
      call dudxyz (ta1,fi,rxm1,sxm1,txm1,jacm1,imesh,1)
      call mapw   (ta1d,lxd,ta1,lx1,1)
      call mapw   (rhod,lxd,rho,lx1,1)    ! Map density
      call col2   (ta1d,rhod,ntotd)       ! Multiply by density
      call mapw   (vxd ,lxd,vx ,lx1,1)
      call col3   (dfid,ta1d,vxd,ntotd)

c     interpolate ta1 and vy onto larger mesh
      call dudxyz  (ta1,fi,rym1,sym1,tym1,jacm1,imesh,2)
      call mapw    (ta1d,lxd,ta1,lx1,1)
      call col2    (ta1d,rhod,ntotd)      ! Multiply by density
      call mapw    (vyd ,lxd,vy ,lx1,1)
      call addcol3 (dfid,ta1d,vyd,ntotd)

      if (if3d) then

c        interpolate ta1 and vy onto larger mesh
         call dudxyz  (ta1,fi,rzm1,szm1,tzm1,jacm1,imesh,3)
         call mapw    (ta1d,lxd,ta1,lx1,1)
         call col2    (ta1d,rhod,ntotd)       ! Multiply by density
         call mapw    (vzd ,lxd,vz ,lx1,1)
         call addcol3 (dfid,ta1d,vzd,ntotd)

      endif

c     Now, *project* DFID onto mesh 1 using L2 projection

      call mapwp(dfid,lxd,dfi,lx1,-1)

      return
      end subroutine conv1d_rho
C------------------------------------------------------------------------








