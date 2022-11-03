!======================================================================
!
!     Author: Prabal Negi
!     Description: Cylindrical coordinates heat solver
!
!     Outside dependencies: 
!     generic_subs.f          : ortho_subspace()     
!
!      
!======================================================================
!---------------------------------------------------------------------- 
      subroutine cdscal_cyl (igeom)

!     Solve the convection-diffusion equation for passive scalar IPSCAL

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'MVGEOM'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'ORTHOT'

      logical          ifprint
      common  /cprint/ ifprint
      logical          ifconv

      real ta,tb
      common /scrns/ ta(lx1,ly1,lz1,lelt)
     $              ,tb(lx1,ly1,lz1,lelt)

      real h1,h2
      common /scrvh/ h1(lx1,ly1,lz1,lelt)
     $              ,h2(lx1,ly1,lz1,lelt)

      integer igeom

      integer iter,ifld1,nel,n
      integer isd,intype

      if (ifdgfld(ifield)) then
        call cdscal_dg(igeom)
        return
      endif


      ifld1 = ifield-1
      napproxt(1,ifld1) = laxtt

      nel    = nelfld(ifield)
      n   = lx1*ly1*lz1*nel

      if (igeom.eq.1) then   ! geometry at t^{n-1}
        call makeq
        call lagscal
      else                   ! geometry at t^n
         if (ifprint) then
           if (ifield.eq.2.and.nid.eq.0)
     $         write (6,*) 
     $      ' Cylindrical Temperature/Passive scalar solution'
         endif

         write(name4t,1) ifld1-1
    1    format('PS',i2)
         if (ifield.eq.2) write(name4t,'(A4)') 'TEMP'

!        New geometry

         isd = 1
         if (ifaxis.and.ifaziv.and.ifield.eq.2) isd = 2
c        if (ifaxis.and.ifmhd) isd = 2 !This is a problem if T is to be T!

         do 1000 iter=1,nmxnl ! iterate for nonlin. prob. (e.g. radiation b.c.)

         intype = 0
         if (iftran) intype = -1
         call sethlm  (h1,h2,intype)
         call bcneusc (ta,-1)
         call add2    (h2,ta,n)
         call bcdirsc (t(1,1,1,1,ifield-1))
         call axhelm  (ta,t(1,1,1,1,ifield-1),h1,h2,imesh,ISD)
         call sub3    (tb,bq(1,1,1,1,ifield-1),ta,n)
         call bcneusc (ta,1)
         call add2    (tb,ta,n)

         if (iftmsh(ifield)) then
           call hsolve  (name4t,ta,tb,h1,h2 
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxt(ifield-1),1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),bintm1)
         else
           call hsolve  (name4t,ta,tb,h1,h2 
     $                   ,tmask(1,1,1,1,ifield-1)
     $                   ,tmult(1,1,1,1,ifield-1)
     $                   ,imesh,tolht(ifield),nmxt(ifield-1),1
     $                   ,approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)
         endif 

         call add2    (t(1,1,1,1,ifield-1),ta,n)

         call cvgnlps (ifconv) ! Check convergence for nonlinear problem 
         if (ifconv) goto 2000

C        Radiation case, smooth convergence, avoid flip-flop (ER).
         call cmult (ta,0.5,n)
         call sub2  (t(1,1,1,1,ifield-1),ta,n)

 1000    continue
 2000    continue

      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine makeq_cyl

C     Generate forcing function for the solution of a passive scalar.
C     !! NOTE: Do not change the content of the array BQ until the current
      
      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MVGEOM'
      include 'INPUT'
      include 'CTIMER'
      include 'TSTEP'

      logical  if_conv_std

      real w1
      common /scruz/ w1(lx1,ly1,lz1,lelt)

      integer nxyz,ntot

      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv

      etime0 = dnekclock()      

      if (nio.eq.0.and.loglevel.gt.2)
     $   write(6,*) 'makeq', ifield

      if_conv_std = .true.
      if (ifmhd.and.ifaxis) if_conv_std = .false. ! conv. treated in induct.f

      call makeq_aux ! nekuq, etc.

      if (ifadvc(ifield) .and. if_conv_std) then
        if (ifcvfld(ifield)) then
          if (ifmvbd) then
            call sub2 (vx, wx, ntot)
            call sub2 (vy, wy, ntot)
            call sub2 (vz, wz, ntot)
          endif

          call convab

          if (ifmvbd) then
            call add2 (vx, wx, ntot)
            call add2 (vy, wy, ntot)
            call add2 (vz, wz, ntot)
          endif
        else
          if (.not.ifchar) call convab
        endif
      endif

      if (iftran) then

         if (ifcvfld(ifield)) then

           if (ifdiff(ifield)) then
              ntot = lx1*ly1*lz1*nelfld(ifield)
              call wlaplacian(w1,t(1,1,1,1,ifield-1),
     &                        vdiff(1,1,1,1,ifield),ifield)
              call add2(bq(1,1,1,1,ifield-1),w1,ntot)
           endif

         else

           if (ifmvbd.and..not.ifchar) call admesht

           call makeabq

           if (ifchar.and.ifadvc(ifield)) then
              call convch
           else
              call makebdq
           endif

         endif

      endif

      tmakq=tmakq+(dnekclock()-etime0)

      return
      end
