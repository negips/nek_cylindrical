!======================================================================
!     Poisson solver for Level Set Method Distance function.
!     Author: Prabal S. Negi
!
!====================================================================== 
!---------------------------------------------------------------------- 

      subroutine ps_heatflow()

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'TSTEP'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real dtold,timeold
      real dtlagold(10)
      integer istepold
      integer nsteps_new

      integer i,n
      real d,d1,dsc
      real quad_scal
      real eps1

      real coeffs(3)
      common /quadratic_coeffs/ coeffs

      integer icalld
      data icalld /0/
      save icalld

      logical ifforc
      real xin(3),yout(3)

      xin(1) = 0.4
      xin(2) = 0.45
      xin(3) = 0.5

      yout(1) = 0.4
      yout(2) = 2.0
      yout(3) = 5.0

      if (icalld.eq.0) then
        call quad_coeff(coeffs,xin,yout)
        icalld = icalld+1
        write(6,*) 'Quad Coeffs:',coeffs(1),coeffs(2),coeffs(3)
      endif  

      ifadvc(3) = .false.
      ifforc    = .true.      ! set boundary as a forcing

      istepold = istep
      dtold    = dt
      timeold  = time
      call copy(dtlagold,dtlag,10)

      nsteps_new = 1
      time       = 0.0

!     Generate IC
      n = lx1*ly1*lz1*nelv
      ifield = 3
      eps1   = 1.0e-0
      do i=1,n
        d = t(i,1,1,1,1)-0.5
        if (abs(d).gt.0.4) then
          d1  = d
          dsc = quad_scal(abs(d),coeffs(1),coeffs(2),coeffs(3))
          d   = sign(dsc,d1)
        endif  
        t(i,1,1,1,ifield-1) = exp(-(d/eps1)**2)
      enddo

      if (ifforc) then
        call rzero(t(1,1,1,1,ifield-1),n)
      endif  

      param(12) = -1.0e-3
      dt        = param(12)
      do istep=1,nsteps_new
        call ps_advance()        
      enddo

      time = timeold
      dt = dtold
      istep = istepold
      call copy(dtlag,dtlagold,10)


      return
      end subroutine ps_heatflow
!---------------------------------------------------------------------- 
      subroutine ps_advance

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

!        only for moving boundaries         
         if (ifgeom) then
            if (.not.ifrich) call gengeom (igeom)
            call geneig  (igeom)
         endif

!        std. nek case
         call ps_heat (igeom)

      enddo

      return
      end subroutine ps_advance

c-----------------------------------------------------------------------

      subroutine ps_heat (igeom)
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
     &    write(*,'(13x,a)') 'Solving for Distance Field'

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
     &   istep,'  Distance Field done',time,dnekclock()-ts

      return
      end subroutine ps_heat
c-----------------------------------------------------------------------

      real function quad_scal(x,a,b,c)

      real x
      real a,b,c

!      a = 880.0
!      b = -696.0
!      c = 138.0

      quad_scal = a*x**2 + b*x + c

      return
      end function
!---------------------------------------------------------------------- 

      subroutine quad_coeff(coeff,x,y)

      implicit none

      real x(3)
      real y(3)

      real coeff(3)

      real matA(3,3)

      real rmult(3)
      integer indr(3),indc(3),ipiv(3)

      integer i
      integer ierr

!      x(1) = 0.4
!      x(2) = 0.45
!      x(3) = 0.5
!
!      y(1) = 0.4
!      y(2) = 2.0
!      y(3) = 5.0

      do i=1,3
        matA(i,1) = x(i)**2
        matA(i,2) = x(i)
        matA(i,3) = 1.0
      enddo  

!     Invert A      
      call gaujordf(matA,3,3,indr,indc,ipiv,ierr,rmult)

      call mxm(matA,3,y,3,coeff,1) 


      return
      end subroutine quad_coeff
!---------------------------------------------------------------------- 









