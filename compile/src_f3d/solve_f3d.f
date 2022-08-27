!======================================================================
!     Routines for introducing third component in a 2d simulation
!     Author: Prabal S. Negi
!     Right now we assume the third component is homogeneous.
!     Later, a fourier dependence can be added for the linearized solve.
!
!====================================================================== 
!-----------------------------------------------------------------------

      subroutine makef_f3d

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'   ! ifield

      include 'F3D'


      integer nxyz,ntot1

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)


      nxyz = lx1*ly1*lz1
      ntot1 = nxyz*nelv

      if (.not.iff3d) then
        call rzero(bfz,ntot1)    
        return
      endif

      call mntr_log(f3d_id,log_f3d,'Generate Forcing')

!     Build user defined forcing for uz
      call makeuf_f3d

      if (filterType.eq.2) call make_hpf_f3d(iff3d)

!     Put vx,vy,vz on Dealiased grid (rst form)
!     Standard nek routine. Don't need to change anything (yet) 
!      call set_convect_new(vxd,vyd,vzd,vx,vy,vz)
      call advab_f3d

      if (ifmvbd) call admeshv_f3d

      if (iftran) call makeabf_f3d

      if ((iftran.and..not.ifchar).or.
     $    (iftran.and..not.ifnav.and.ifchar)) call makebdf_f3d


      return
      end subroutine makef_f3d

!----------------------------------------------------------------------

      subroutine makeuf_f3d

!     Compute and add: (1) user specified forcing function (FX,FY,FZ)

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'PARALLEL'
      include 'NEKUSE'

      include 'F3D'

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      integer ntot1,iel,i,j,k,ielg


      ntot1 = lx1*ly1*lz1*nelv

      time = time-dt
      call rzero(bfz,ntot1)

      do 100 iel=1,nelv
         ielg = lglel(iel)
         do 100 k=1,lz1
         do 100 j=1,ly1
         do 100 i=1,lx1
            call nekasgn (i,j,k,iel)
            call userf   (i,j,k,ielg)
            bfx(i,j,k,iel) = ffx
            bfy(i,j,k,iel) = ffy
            bfz(i,j,k,iel) = ffz
 100  continue

      call col2  (bfx,bm1,ntot1)
      call col2  (bfy,bm1,ntot1)
      call col2  (bfz,bm1,ntot1)
      time = time+dt

      return
      end subroutine makeuf_f3d

!----------------------------------------------------------------------
      subroutine advab_f3d

!     Eulerian scheme, add convection term to forcing function 
!     at current time step.

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'GEOM'

      include 'F3D'

      include 'TEST'

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      integer ntot1

      ntot1 = lx1*ly1*lz1*nelv

!      call setup_convect(2)
      ntot1 = lx1*ly1*lz1*nelv
      call convop  (ta1,vx)
      call convop  (ta2,vy)
      call convop  (ta3,vz)

      call subcol3 (bfx,bm1,ta1,ntot1)
      call subcol3 (bfy,bm1,ta2,ntot1)
      call subcol3 (bfz,bm1,ta3,ntot1)

      if (ifcyl_f3d) then
         call invers2(ta1,ym1,ntot1)            ! 1/R
             
         call convect_w_f3d(ta2,vz,vz)
         call col2(ta2,ta1,ntot1)
         call Xaddcol3(bfy,bm1,ta2,ntot1)  ! Note: This is added (on the rhs)

         call convect_w_f3d(ta3,vz,vy)
         call col2(ta3,ta1,ntot1)
         call subcol3(bfz,bm1,ta3,ntot1)   ! This is subtracted (on the rhs)
       endif  



      return
      end subroutine advab_f3d
!-----------------------------------------------------------------------
      subroutine admeshv_f3d

!     Eulerian scheme, add convection term to forcing function 
!     at current time step.

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'GEOM'
      include 'MVGEOM'

      include 'F3D'

      include 'TEST'

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      integer ntot1

      ntot1 = lx1*ly1*lz1*nelv

!     Store Original velocities (arrays from F3D)
      call copy3(ur1x,ur2x,ur3x,vx,vy,vz,ntot1)

!     This is the new convecting velocity field
      call copy3(vx,vy,vz,wx,wy,wz,ntot1)
      call chsign(vx,ntot1)
      call chsign(vy,ntot1)
      call chsign(vz,ntot1)

!      call setup_convect(2)
      call convop  (ta1,ur1x)
      call convop  (ta2,ur2x)
      call convop  (ta3,ur3x)

      call subcol3 (bfx,bm1,ta1,ntot1)
      call subcol3 (bfy,bm1,ta2,ntot1)
      call subcol3 (bfz,bm1,ta3,ntot1)

      if (ifcyl_f3d) then
        call invers2(ta1,ym1,ntot1)            ! 1/R
            
        call convect_w_f3d(ta2,vz,ur3x)        ! w_\theta*u_\theta
        call col2(ta2,ta1,ntot1)
        call Xaddcol3(bfy,bm1,ta2,ntot1)  ! Note: This is added (on the rhs)

        call convect_w_f3d(ta3,vz,vy)          ! w_\theta*u_R
        call col2(ta3,ta1,ntot1)
        call subcol3(bfz,bm1,ta3,ntot1)   ! This is subtracted (on the rhs)
      endif  

!     Restore velocities      
      call copy3(vx,vy,vz,ur1x,ur2x,ur3x,ntot1)

      return
      end subroutine admeshv_f3d
c-----------------------------------------------------------------------

      subroutine makeabf_f3d
!
!     Eulerian scheme, add convection term to forcing function 
!     at current time step.

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'
      include 'INPUT'

      include 'F3D'


      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      integer ntot1
      real ab0,ab1,ab2

      ntot1 = lx1*ly1*lz1*nelv


      ab0 = ab(1)
      ab1 = ab(2)
      ab2 = ab(3)
      call add3s2 (ta1,abx1,abx2,ab1,ab2,ntot1)
      call add3s2 (ta2,aby1,aby2,ab1,ab2,ntot1)
      call add3s2 (ta3,abz1,abz2,ab1,ab2,ntot1)

      call copy   (abx2,abx1,ntot1)
      call copy   (aby2,aby1,ntot1)
      call copy   (abz2,abz1,ntot1)

      call copy   (abx1,bfx,ntot1)
      call copy   (aby1,bfy,ntot1)
      call copy   (abz1,bfz,ntot1)

      call add2s1 (bfx,ta1,ab0,ntot1)
      call add2s1 (bfy,ta2,ab0,ntot1)
      call add2s1 (bfz,ta3,ab0,ntot1)

!     multiply by density

      if (.not.iflomach) call col2   (bfx,vtrans,ntot1)
      if (.not.iflomach) call col2   (bfy,vtrans,ntot1)
      if (.not.iflomach) call col2   (bfz,vtrans,ntot1)

      if (ldim.eq.3) then
!       Something went wrong
        write(6,*) 'Inconsistent parameter setting,ndim,iff3d',
     $              ndim,iff3d
        call exitt
      endif


      return
      end subroutine makeabf_f3d

!-----------------------------------------------------------------------
      subroutine makebdf_f3d

!     Add contributions to F from lagged BD terms.

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'

      include 'F3D'

      real ta1,ta2,ta3,tb1,tb2,tb3,h2
      common /scrns/ ta1(lx1,ly1,lz1,lelv)
     $ ,             ta2(lx1,ly1,lz1,lelv)
     $ ,             ta3(lx1,ly1,lz1,lelv)
     $ ,             tb1(lx1,ly1,lz1,lelv)
     $ ,             tb2(lx1,ly1,lz1,lelv)
     $ ,             tb3(lx1,ly1,lz1,lelv)
     $ ,             h2 (lx1,ly1,lz1,lelv)


      integer ilag,ntot1
      real const

      ntot1 = lx1*ly1*lz1*nelv
      const = 1./dt

      if(iflomach) then
        call cfill(h2,const,ntot1)
      else
        call cmult2(h2,vtrans(1,1,1,1,ifield),const,ntot1)
      endif

      call opcolv3c (tb1,tb2,tb3,vx,vy,vz,bm1,bd(2))
      call col3(tb3,vz,bm1,ntot1)
      call cmult(tb3,bd(2),ntot1)

      do 100 ilag=2,nbd
         if (ifgeom) then
            call opcolv3c(ta1,ta2,ta3,vxlag (1,1,1,1,ilag-1),
     $                                vylag (1,1,1,1,ilag-1),
     $                                vzlag (1,1,1,1,ilag-1),
     $                                bm1lag(1,1,1,1,ilag-1),bd(ilag+1))

            call col3(ta3,vzlag(1,1,1,1,ilag-1),bm1lag(1,1,1,1,ilag-1),
     $                                          ntot1)
            call cmult(ta3,bd(ilag+1),ntot1)

         else
            call opcolv3c(ta1,ta2,ta3,vxlag (1,1,1,1,ilag-1),
     $                                vylag (1,1,1,1,ilag-1),
     $                                vzlag (1,1,1,1,ilag-1),
     $                                bm1                   ,bd(ilag+1))

            call col3(ta3,vzlag(1,1,1,1,ilag-1),bm1,ntot1)
            call cmult(ta3,bd(ilag+1),ntot1)
         endif
         call opadd2  (tb1,tb2,tb3,ta1,ta2,ta3)
         call add2 (tb3,ta3,ntot1)

 100  continue
      call opadd2col(bfx,bfy,bfz,tb1,tb2,tb3,h2)
      call xaddcol3(bfz,tb3,h2,ntot1)     


      return
      end subroutine makebdf_f3d
!-----------------------------------------------------------------------

      subroutine lagvel_f3d

!     Keep old velocity field(s)

      implicit none 

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'

      integer ilag,ntot1

      ntot1 = lx1*ly1*lz1*nelv

      do 100 ilag=3-1,2,-1
         call copy (vxlag (1,1,1,1,ilag),vxlag (1,1,1,1,ilag-1),ntot1)
         call copy (vylag (1,1,1,1,ilag),vylag (1,1,1,1,ilag-1),ntot1)
         call copy (vzlag (1,1,1,1,ilag),vzlag (1,1,1,1,ilag-1),ntot1)
 100  continue

      call copy3 (vxlag,vylag,vzlag,vx,vy,vz,ntot1)

      return
      end subroutine lagvel_f3d
!----------------------------------------------------------------------
      subroutine ophx_f3d (out1,out2,out3,inp1,inp2,inp3,h1,h2)

!     OUT = (H1*A+H2*B) * INP  

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'

      include 'F3D'


      real out1 (lx1,ly1,lz1,1)
      real out2 (lx1,ly1,lz1,1)
      real out3 (lx1,ly1,lz1,1)
      real inp1 (lx1,ly1,lz1,1)
      real inp2 (lx1,ly1,lz1,1)
      real inp3 (lx1,ly1,lz1,1)
      real h1   (lx1,ly1,lz1,1)
      real h2   (lx1,ly1,lz1,1)

      integer imesh,matmod
      

      imesh = 1

      if (ifstrs) then
         matmod = 0
         call axhmsf (out1,out2,out3,inp1,inp2,inp3,h1,h2,matmod)
      else

!        the numbers are only needed for axis-symmetric formulation
!        need to come back to this later. 
         call axhelm (out1,inp1,h1,h2,imesh,1)
         call axhelm (out2,inp2,h1,h2,imesh,2)
         call axhelm (out3,inp3,h1,h2,imesh,3)
      endif

      return
      end subroutine ophx_f3d
!-----------------------------------------------------------------------
      subroutine cresvif_f3d (resv1,resv2,resv3,h1,h2)

!     Compute startresidual/right-hand-side in the velocity solver

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'

      include 'F3D'

      include 'TEST'

      real           resv1 (lx1,ly1,lz1,1)
      real           resv2 (lx1,ly1,lz1,1)
      real           resv3 (lx1,ly1,lz1,1)
      real           h1    (lx1,ly1,lz1,1)
      real           h2    (lx1,ly1,lz1,1)
      common /scruz/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)

      integer igeom
      common /cgeom/ igeom

      integer ntot1,ntot2
      integer intype

      call mntr_log(f3d_id,log_f3d,'Generate Residual')

      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv

!     for 3d solve

      intype = -1
   
      call bcneusc_f3d(w1,intype)
      call add2(h2,w1,ntot1)

      call bcdirvc_cyl(vx,vy,vz,v1mask,v2mask,v3mask)

      call extrapp (pr,prlag)
      
      call opgradt_f3d (resv1,resv2,resv3,pr)

      call rzero(resv3,ntot1)             ! homogeneous in z

      call add2_3(resv1,resv2,resv3,bfx,bfy,bfz,ntot1)

!     Ax
      call axhmsf_cyl_real(w1,w2,w3,vx,vy,vz,h1,h2)
      call sub2(resv1,w1,ntot1)
      call sub2(resv2,w2,ntot1)
      call sub2(resv3,w3,ntot1)

      return
      end subroutine cresvif_f3d
!-----------------------------------------------------------------------
      subroutine plan3_f3d (igeom)

!     Compute pressure and velocity using consistent approximation spaces.     
!     Operator splitting technique.

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'EIGEN'
      include 'SOLN'
      include 'TSTEP'

      include 'F3D'

      include 'TEST'

      real resv1,resv2,resv3
      real dv1,dv2,dv3
      real h1,h2
      common /scrns1/  resv1 (lx1,ly1,lz1,lelv)
     $ ,               resv2 (lx1,ly1,lz1,lelv)
     $ ,               resv3 (lx1,ly1,lz1,lelv)
     $ ,               dv1   (lx1,ly1,lz1,lelv)
     $ ,               dv2   (lx1,ly1,lz1,lelv)
     $ ,               dv3   (lx1,ly1,lz1,lelv)
      common /scrvh/  h1    (lx1,ly1,lz1,lelv)
     $ ,              h2    (lx1,ly1,lz1,lelv)

      integer intype
      integer igeom
      integer ntot1
      integer i,nit

      ntot1 = lx1*ly1*lz1*nelv 

      if (ifmvbd.and.ifcyl_f3d) then
        ifaxis = .false.
      endif  

      if (igeom.eq.1) then

!        old geometry
         call makef_f3d

      else

         if (igeom.eq.2) call lagvel_f3d

!        new geometry, new b.c.
         nit = 1
         do i = 1,nit
           intype = -1
           call sethlm  (h1,h2,intype)
           call cresvif_f3d (resv1,resv2,resv3,h1,h2)
           call ophinv_real(dv1,dv2,dv3,resv1,resv2,resv3,
     $                      h1,h2,tolhv,nmxv)
           if (i.eq.nit) then
             call add2_3(vx,vy,vz,dv1,dv2,dv3,ntot1)
           else
             call cmult(dv1,0.5,ntot1)
             call cmult(dv2,0.5,ntot1)
             call cmult(dv3,0.5,ntot1)
             call add2_3(vx,vy,vz,dv1,dv2,dv3,ntot1)
           endif 
         enddo  

         call incomprn_real(igeom)

!        Need to reset for reevaluation of mesh            
         if (ifmvbd.and.ifcyl_f3d) then
           ifaxis = .true.
         endif  

      endif

      return
      end subroutine plan3_f3d

!----------------------------------------------------------------------
      subroutine incomprn_real (igeom)
c
c     Project U onto the closest incompressible field
c
c     Output: updated values of U, iproj, proj; and
c             up    := pressure correction req'd to impose div U = 0
c
c
c     Dependencies: ifield ==> which "density" (vtrans) is used.
c

      implicit none

      include 'SIZE'
      include 'SOLN'          ! vxp,vyp,vzp,prp,jp
      include 'INPUT'         ! npert
      include 'TSTEP'         ! dt,ifield
      include 'CTIMER'
      include 'GEOM'          ! YM1,YM2
      include 'MASS'

      include 'F3D'

      include 'TEST'

      real h1,h2,h2inv
      common /scrvh/ h1   (lx1,ly1,lz1,lelv)
     $ ,             h2   (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv(lx1,ly1,lz1,lelv)

      real dp2
      common /scrch/ dp2(lx2,ly2,lz2,lelv)
      logical ifprjp

      real dummy
      common /scrcg/ dummy(lx1*ly1*lz1*lelt) 

      real w1,w2,w3
      real dv1,dv2,dv3,dp
      common /scrns1/ w1    (lx1,ly1,lz1,lelv)
     $ ,              w2    (lx1,ly1,lz1,lelv)
     $ ,              w3    (lx1,ly1,lz1,lelv)
     $ ,              dv1   (lx1,ly1,lz1,lelv)
     $ ,              dv2   (lx1,ly1,lz1,lelv)
     $ ,              dv3   (lx1,ly1,lz1,lelv)
     $ ,              dp    (lx2,ly2,lz2,lelv)

      integer ntot1,ntot2,intype,istart

      real bddt,bddti,const

      integer igeom

      real dnorm
      real divv,bdivv
      common /scruz/  divv (lx2,ly2,lz2,lelv)
     $ ,              bdivv(lx2,ly2,lz2,lelv)

      real glsc2

      if (igeom.eq.1) return

      if (icalld.eq.0) tpres=0.0
      icalld=icalld+1
      npres=icalld

      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv
      bddt   = bd(1)/dt
      bddti  = dt/bd(1)

      call rzero (h1,ntot1)
      call cmult2(h2,vtrans(1,1,1,1,ifield),bddt,ntot1)
      call invers2(h2inv,h2,ntot1)

!     Note: OPDIV already contains the mass matrix multiplication
      call opdiv_f3d(dp,vx,vy,vz)
      call chsign(dp,ntot2)
      call ortho (dp)

!     prabal. checking NaNs
      call copy(bdivv,dp,ntot2)
      call col3(divv,bdivv,bm2inv,ntot2)
      dnorm = sqrt(glsc2(divv,bdivv,ntot2)/volvm2)
      if (isnan(dnorm)) then
        if (nio.eq.0) write(6,*) 'NaN found'
        call copy(pr,bm2,ntot2)
        call outpost(vx,vy,vz,dp,t,'nan')

        call exitt
      endif


      ifprjp=.false.    ! project out previous pressure solutions?
      istart=param(95)  
      if (istep.ge.istart.and.istart.ne.0) ifprjp=.true.

      ! Most likely, the following can be commented out. (pff, 1/6/2010)
      if (npert.gt.1.or.ifbase)            ifprjp=.false.

      intype =  1             ! Changing integration type here.
                              ! Need to modify cdabdtp accordingly
                              ! Also need to modify uzprec

!      if (nio.eq.0.and.igeom.eq.2) write(6,5) istep,time
      call esolver(dp,h1,h2,h2inv,intype)

!   5  format(i9,1pe14.7,' Pressure Solve:')

!     Update Pressure
      call add2(pr,dp,ntot2)

!     Update Velocity      
      call opgradt_f3d(w1 ,w2 ,w3 ,dp)
      if3d = .true.
      call opbinv_f3d(dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)
      if3d = .false.
      
      call add2_3 (vx,vy,vz,dv1,dv2,dv3, ntot1)

      return
      end subroutine incomprn_real
!------------------------------------------------------------------------
      subroutine make_hpf_f3d(iff3d)

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'         ! param(110),(111)
      include 'TSTEP'         ! ifield
      include 'MASS'          ! BM1


      integer nxyz
      parameter (nxyz=lx1*ly1*lz1)
      integer n

      integer lm,lm2
      parameter (lm=40)
      parameter (lm2=lm*lm)
      real hpf_filter(lm2)

      integer hpf_kut
      real hpf_chi
      logical hpf_ifboyd

      integer nel

      integer icalld
      save icalld
      data icalld /0/

      real TA1,TA2,TA3
      COMMON /SCRUZ/ TA1 (LX1,LY1,LZ1,LELV)
     $ ,             TA2 (LX1,LY1,LZ1,LELV)
     $ ,             TA3 (LX1,LY1,LZ1,LELV)

      real hpf_op(lx1,lx1)
      save hpf_op

      logical iff3d

c---------------------------------------- 
      if(.not. iffilter(ifield)) return

      hpf_kut = int(param(101))+1
      hpf_chi = -1.0*abs(param(103))
c     Boyd transform to preserve element boundary values is 
c     linearly unstable when used as forcing.
      hpf_ifboyd = .false.      

      nel = nelv
      n = nxyz*nel

      if (hpf_chi.eq.0) return
      if(nid.eq.0 .and. loglevel.gt.2) write(6,*) 'apply hpf ',
     $                                 ifield, hpf_kut, hpf_chi

      if (icalld.eq.0) then
c       Create the filter transfer function
        call hpf_trns_fcn(hpf_filter,hpf_kut)

c       Build the matrix to apply the filter function
c       to an input field
        call build_hpf_mat(hpf_op,hpf_filter,hpf_ifboyd)

c       Only initialize once    
        icalld=icalld+1 
      endif

      if (ifield.eq.1) then
c       Apply the filter
c       to velocity fields
        if (jp.eq.0) then 
          call build_hpf_fld(ta1,vx,hpf_op,lx1,lz1)
          call build_hpf_fld(ta2,vy,hpf_op,lx1,lz1)
          if (if3d.or.iff3d) call build_hpf_fld(ta3,vz,hpf_op,lx1,lz1)

c         Multiply by filter weight (chi)
          call cmult(ta1,hpf_chi,n)    
          call cmult(ta2,hpf_chi,n)    
          if (if3d.or.iff3d) call cmult(ta3,hpf_chi,n)    

c         Multiply by Mass matrix 
c         and add to forcing term 
          call opadd2col (bfx,bfy,bfz,ta1,ta2,ta3,bm1)
          if (iff3d) call Xaddcol3(bfz,ta3,bm1,n)
            
        else
c         Apply filter on velocity perturbation fields
          call build_hpf_fld(ta1,vxp(1,jp),hpf_op,lx1,lz1)
          call build_hpf_fld(ta2,vyp(1,jp),hpf_op,lx1,lz1)
          if (if3d.or.iff3d) then
            call build_hpf_fld(ta3,vzp(1,jp),hpf_op,lx1,lz1)
          endif  

c         Multiply by filter weight (chi)
          call cmult(ta1,hpf_chi,n)    
          call cmult(ta2,hpf_chi,n)    
          if (if3d.or.iff3d) call cmult(ta3,hpf_chi,n)    

c         Multiply by Mass matrix 
c         and add to forcing term 
          call opadd2col (bfxp(1,jp),bfyp(1,jp),bfzp(1,jp),
     $                    ta1,ta2,ta3,bm1)
          if (iff3d) call Xaddcol3(bfzp(1,jp),ta3,bm1,n)

        endif        

      else

c       Apply filter to temp/passive scalar fields
        if (jp.eq.0) then
          call build_hpf_fld(ta1,t(1,1,1,1,ifield-1),
     $         hpf_op,lx1,lz1)

c         Multiply by filter weight (chi)
          call cmult(ta1,hpf_chi,n)    

c         Multiply by Mass matrix    
c         and add to source term
          call addcol3(bq(1,1,1,1,ifield-1),ta1,bm1,n)
        else
c         Apply filter on scalar perturbation field                
          call build_hpf_fld(ta1,tp(1,ifield-1,jp),hpf_op,lx1,lz1)

c         Multiply by filter weight (chi)
          call cmult(ta1,hpf_chi,n)    

c         Multiply by Mass matrix    
c         and add to source term
          call addcol3(bqp(1,ifield-1,jp),ta1,bm1,n)
        endif  

      endif

      return
      end

c----------------------------------------------------------------------

      subroutine bcneusc_f3d(s,itype)

!     Apply Neumann boundary conditions to surface of SYM boundary
!     Use IFIELD as a guide to which boundary conditions are to be applied.
!
!     If ITYPE = 1, then S is returned as the rhs contribution to the 
!                   volumetric flux.
!
!     If ITYPE =-1, then S is returned as the lhs contribution to the 
!                   diagonal of A.

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'
      include 'GEOM'
      include 'PARALLEL'
      include 'MASS'

      include 'CTIMER'
      include 'NEKUSE'

      include 'F3D'
      include 'FS_ALE'

      real s(lx1,ly1,lz1,lelt)
      real v
      common  /nekcb/ cb
      character cb*3

      integer itype
      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer iface,ie,nfaces,nxyz,nel,ntot
      integer ia,ix,iy,iz,ieg
      real alpha,beta

      real dummy
      common /scrcg/ dummy(lx1,ly1,lz1,lelt)

      real xs,xf,fsx
      integer n
      real tol

      if (icalld.eq.0) then
         tusbc=0.0
         nusbc=0
         icalld=icalld+1
      endif
      nusbc=nusbc+1
      etime1=dnekclock()

      nfaces = 2*ldim
      nxyz   = lx1*ly1*lz1
      nel    = nelfld(ifield)
      ntot   = nxyz*nel
      call rzero(s,ntot)

!     Broadcast location of the free surface
!     At the moment we only need x coord.
      call copy(dummy,xm1,ntot)
      call col2(dummy,fs_mask,ntot)
      call fgslib_gs_op(fs_gs_handle,dummy,1,1,0)     ! 1 ==> +
      call col2(dummy,fs_vmult,ntot)

      tol = 1.0e-12
      if (itype.eq.-1) then

!        Compute diagonal contributions to accomodate Robin boundary conditions
         do 1000 ie=1,nel
         do 1000 iface=1,nfaces
            ieg=lglel(ie)
            cb =cbc(iface,ie,ifield)
!            if (cb.eq.'C  ' .or. cb.eq.'C  ' .or.
!     $          cb.eq.'R  ' .or. cb.eq.'R  ') then
            if (cb.eq.'SYM') then
               ia=0
!              ia is area counter, assumes advancing fastest index first. (ix...iy...iz)
               call facind (kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,iface)
               do 100 iz=kz1,kz2
               do 100 iy=ky1,ky2
               do 100 ix=kx1,kx2
                  ia = ia + 1
                  v  = vx(ix,iy,iz,ie)
                  if (optlevel.le.2) call nekasgn (ix,iy,iz,ie)
!                  call userbc  (ix,iy,iz,iface,ieg)
                  n  = 2
                  fsx = dummy(ix,iy,iz,ie)      ! free surface position
                  xs  = fsx+slipl_f3d           ! start of blending
                  xf  = xs+blendl_f3d           ! end of blending
                  if (abs(x-fsx).lt.slipl_f3d) then
!                   Free slip within the slip length                    
                    alpha = 0.0
                    beta  = 0.0
                    hc    = 0.0    
                  else
!                   Smoothly blend from free slip to Dirichlet                    
                    alpha = ((x-xs)/(xf-xs))**n
                    if (x.le.xf) then
                      if (abs(alpha-1.0).lt.tol) then
                        alpha = 1.0-tol
                        beta  = tol
                      else
                        beta  = 1.0-alpha
                      endif
                    else
                      alpha = 1.0-tol
                      beta  = tol
                    endif  
                    hc      = alpha/beta
                  endif  
!                 We have inverse bm1 since this is just added to h2
!                 Which gets multiplied by bm1 later                  
                  s(ix,iy,iz,ie) = s(ix,iy,iz,ie) +
     $               hc*area(ia,1,iface,ie)/bm1(ix,iy,iz,ie)
  100          continue
            endif
 1000    continue
      endif
      if (itype.eq.1) then

!        add passive scalar fluxes to rhs

         do 2000 ie=1,nel
         do 2000 iface=1,nfaces
            ieg=lglel(ie)
            cb =cbc(iface,ie,ifield)
            if (cb.eq.'SYM') then
!              Add local weighted flux values to rhs, S.
!              ia is area counter, assumes advancing fastest index first. (IX...IY...IZ)
               ia=0
               call facind (kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,iface)
               do 200 iz=kz1,kz2
               do 200 iy=ky1,ky2
               do 200 ix=kx1,kx2
                  ia = ia + 1
                  v = vx(ix,iy,iz,ie)
!                 Add computed fluxes to boundary surfaces:
!                  s(ix,iy,iz,ie) = s(ix,iy,iz,ie)
!     $                           + flux*area(ia,1,iface,ie)
  200          continue
            endif
 2000    continue
      endif

      tusbc=tusbc+(dnekclock()-etime1)

      return
      end subroutine bcneusc_f3d
c-----------------------------------------------------------------------


!----------------------------------------------------------------------


c      COMMON /CTMP1/ DUMMY1(LCTMP1)
c      COMMON /CTMP0/ DUMMY0(LCTMP0)
c
c      COMMON /SCRNS/ DUMMY2(LX1,LY1,LZ1,LELT,7)
c      COMMON /SCRUZ/ DUMMY3(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCREV/ DUMMY4(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRVH/ DUMMY5(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRMG/ DUMMY6(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCRCH/ DUMMY7(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRSF/ DUMMY8(LX1,LY1,LZ1,LELT,3)
c      COMMON /SCRCG/ DUMM10(LX1,LY1,LZ1,LELT,1)

C       X             X-coordinate
C       Y             Y-coordinate
C       Z             Z-coordinate
C       UX            X-velocity
C       UY            Y-velocity
C       UZ            Z-velocity
C       TEMP          Temperature
C       PS1           Passive scalar No. 1
C       PS2           Passive scalar No. 2
C        .             .
C        .             .
C       PS9           Passive scalar No. 9
C       SI2           Strainrate invariant II
C       SI3           Strainrate invariant III
C
C     Variables to be defined by user for imposition of
C     boundary conditions :
C
C       SH1           Shear component No. 1
C       SH2           Shear component No. 2
C       TRX           X-traction
C       TRY           Y-traction
C       TRZ           Z-traction
C       SIGMA         Surface-tension coefficient
C       FLUX          Flux
C       HC            Convection heat transfer coefficient
C       HRAD          Radiation  heat transfer coefficient
C       TINF          Temperature at infinity


