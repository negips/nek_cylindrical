!======================================================================
!     Poisson solver for Level Set Method Distance function.
!     Author: Prabal S. Negi
!
!====================================================================== 
!---------------------------------------------------------------------- 

      subroutine phi_corr()

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'TSTEP'

      include 'LSM'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real dtold,timeold
      real dtlagold(10)
      integer istepold
      integer nsteps_new

      integer i,n

      integer icalld
      data icalld /0/
      save icalld

!     This holds the gradients      
      real ta1(lv)
      real ta2(lv)
      real ta3(lv)
      real ta4(lv)

      common /scruz/ ta1,ta2,ta3,ta4

!     Used as work arrays 
      real tb1(lv)
      real tb2(lv)
      real tb3(lv)
      real tb4(lv)
      common /scrmg/ tb1,tb2,tb3,tb4

!     normalized phi.  0 <= p <= 1.0
      real nphi(lv)
      common /scrcg/ nphi

      integer intype

      logical ifspvd
      common /spvd/ ifspvd

      real ft

      if (istep.eq.0) return

      ifspvd = .false.

      lsm_simp = .false.

      ifadvc(3) = .false.

      istepold = istep
      dtold    = dt
      timeold  = time
      call copy(dtlagold,dtlag,10)
      ft       = filtertype

      nsteps_new = 100
      time       = 0.0

!     Generate IC
      n = lx1*ly1*lz1*nelv
      ifield = 3
      call copy(t(1,1,1,1,2),t(1,1,1,1,1),n)
      call copy(lsm_phi,t(1,1,1,1,2),n)

      intype = 1
      call lsm_gennormals(t(1,1,1,1,2))         ! normals
      call lsm_grad_nlmap(ta1,ta2,ta3,t(1,1,1,1,2))
      call normalize_phi(nphi,t(1,1,1,1,2))
!      call lsm_forc_gkreiss(ta1,ta2,ta3,ta4,    ! forcing
!     $            tb4,nphi,intype)
      call lsm_forc_gkreiss_v1(nphi,intype)
!      call lsm_forc_shukla(nphi,intype)
!      call lsm_forc_shukla(nphi,-1)
!      call outpost2(lsm_forc,vdiff(1,1,1,1,ifield),vz,pr,t,2,'cor')

      lsm_simp = .false.
      ifspvd = .false.

      param(12) = 1.0e-6
      dt        = param(12)
      filterType = 0
      do istep=1,nsteps_new
        intype = 1
        call lsm_grad_nlmap(ta1,ta2,ta3,t(1,1,1,1,2))
        call normalize_phi(nphi,t(1,1,1,1,2))
        call lsm_forc_gkreiss_v1(nphi,intype)
!        call lsm_forc_shukla(nphi,intype)
        call copy(lsm_phi,t(1,1,1,1,2),n)
        call corr_advance()
!        if (mod(istep,10).eq.0) then
!          call lsm_gennormals(t(1,1,1,1,2))
!        endif   
!        if (mod(istep,1).eq.0) then
!          call outpost2(lsm_forc,tb4,vz,pr,t,2,'cor')
!        endif  
      enddo

      call copy(tb4,t(1,1,1,1,2),n)
      call sub2(tb4,t,n)
      call outpost2(lsm_forc,tb4,vz,pr,t,2,'cor')

      call copy(t,t(1,1,1,1,1),n)

      time = timeold
      dt = dtold
      param(12) = dt
      istep = istepold
      call copy(dtlag,dtlagold,10)
      filterType = ft

      return
      end subroutine phi_corr
!---------------------------------------------------------------------- 
      subroutine corr_advance

      implicit none

      include 'SIZE'
      include 'INPUT'
!      include 'TOTAL'
      include 'CTIMER'
      include 'GEOM'          ! ifgeom

      integer igeom
      common /cgeom/ igeom


      call nekgsync

!      call setup_convect(2) ! Save conv vel

      if (iftran) call settime

      call setsolv
      call comment

!     ! PN-2/PN-2 formulation
      call setprop
      do igeom=1,ngeom

!        std. nek case
         call corr_solve(igeom)

      enddo

      return
      end subroutine corr_advance
c-----------------------------------------------------------------------

      subroutine corr_solve (igeom)
C
C     Driver for temperature or passive scalar.
C
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'DEALIAS'

      real*8 ts, dnekclock
      integer igeom

      integer intype

      ts = dnekclock()

      if (nio.eq.0 .and. igeom.eq.2) 
     &    write(*,'(13x,a)') 'Solving for \Phi Correction'

      ifield = 3
      if (idpss(ifield-1).eq.0) then      ! helmholtz
         intype        = -1
         if (.not.iftmsh(ifield)) imesh = 1
         if (     iftmsh(ifield)) imesh = 2
         call unorm
         call settolt            ! set tolerance for temp/ps
         call cdscal(igeom)      ! everything is done here.
      endif

      if (nio.eq.0 .and. igeom.eq.2)
     &   write(*,'(4x,i7,a,1p2e12.4)') 
     &   istep,'  \Phi Correction done',time,dnekclock()-ts

      return
      end subroutine corr_solve
c-----------------------------------------------------------------------

      subroutine normalize_phi(p,phi)

      implicit none

      include 'SIZE'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real p(lv)
      real phi(lv)

      integer n
      real p1,p2

      real glmin,glmax
     

      n = lx1*ly1*lz1*nelv

!     Almost the same thing as done in lsm_grad_nlmap
      p1 = glmin(phi,n)
      p2 = glmax(phi,n)

      call copy(p,phi,n)
      p1   = glmin(p,n)
      call cadd(p,-p1,n)              ! Make sure everything is positive
      p1   = glmax(p,n) 
      p2 = 1.0/p1
      call cmult(p,p2,n)              !  0<= p <= 1

      return
      end subroutine normalize_phi
!---------------------------------------------------------------------- 



