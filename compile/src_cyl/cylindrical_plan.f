!======================================================================
!
!     Author: Prabal Negi
!     Description: Cylindrical coordinates solver
!
!     Routines:
!     plan_cyl                : Main driver
!     cresvif_cyl             : Create Residual for momentum solver            
!     advab_rho_cyl           : Cylindrical advection term (with variable density)
!     advab_rho               : Cylindrical advection term (standard)
!     incomprn_cyl            : Cylindrical Pressure Poisson Solver
!
!     Outside dependencies: 
!     generic_subs.f          : ortho_subspace()     
!
!      
!======================================================================
!---------------------------------------------------------------------- 
      subroutine plan_cyl (igeom)

C     Compute pressure and velocity using consistent approximation spaces.     
C     Operator splitting technique.

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'EIGEN'
      include 'SOLN'
      include 'TSTEP'
      
      real resv1,resv2,resv3
      real dv1,dv2,dv3
      common /scrns/  resv1 (lx1,ly1,lz1,lelv)
     $ ,              resv2 (lx1,ly1,lz1,lelv)
     $ ,              resv3 (lx1,ly1,lz1,lelv)
     $ ,              dv1   (lx1,ly1,lz1,lelv)
     $ ,              dv2   (lx1,ly1,lz1,lelv)
     $ ,              dv3   (lx1,ly1,lz1,lelv)

      real h1,h2
      common /scrvh/  h1    (lx1,ly1,lz1,lelv)
     $ ,              h2    (lx1,ly1,lz1,lelv)

      integer igeom
      integer intype

      if (igeom.eq.1) then

!        old geometry

         if (nio.eq.0) write(6,*) 'Cylindrical Solver'

         call makef_cyl

      else

!        new geometry, new b.c.

!        modify v1mask to account for moving interface
         call modify_mask()

         intype = -1
         call sethlm  (h1,h2,intype)
         call cresvif_cyl (resv1,resv2,resv3,h1,h2)

!        Velocity            
         call ophinv (dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxv)
         call opadd2 (vx,vy,vz,dv1,dv2,dv3)

!        Pressure
         call incomprn_cyl(vx,vy,vz,pr)

      endif

      return
      end
!---------------------------------------------------------------------- 
      subroutine cresvif_cyl (resv1,resv2,resv3,h1,h2)

!     Compute start residual/right-hand-side in the velocity solver

      implicit none

      include 'SIZE'
      include 'SOLN'

      real           resv1 (lx1,ly1,lz1,1)
      real           resv2 (lx1,ly1,lz1,1)
      real           resv3 (lx1,ly1,lz1,1)
      real           h1    (lx1,ly1,lz1,1)
      real           h2    (lx1,ly1,lz1,1)

      real w1,w2,w3
      common /scruz/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)

      integer igeom
      common /cgeom/ igeom

      integer ntot1,ntot2

      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv
      if (igeom.eq.2) call lagvel 
      call bcdirvc (vx,vy,vz,v1mask,v2mask,v3mask)
      call bcneutr

      call extrapp (pr,prlag)
      call opgradt_cyl (resv1,resv2,resv3,pr)
      call opadd2  (resv1,resv2,resv3,bfx,bfy,bfz)
      call ophx    (w1,w2,w3,vx,vy,vz,h1,h2)
      call opsub2  (resv1,resv2,resv3,w1,w2,w3)

      return
      end
!-----------------------------------------------------------------------
      subroutine makef_cyl

C     Compute and add: (1) user specified forcing function (FX,FY,FZ)
C                      (2) driving force due to natural convection
C                      (3) convection term
C     !! NOTE: Do not change the arrays BFX, BFY, BFZ until the
C              current time step is completed.

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'INPUT'
      include 'TSTEP'
      include 'CTIMER'
      include 'MVGEOM'

      etime1 = dnekclock()
                                                call makeuf
      if (filterType.eq.2)                      call make_hpf
      if (ifnav .and..not.ifchar) then
!       Fluid Convection term        
        if (ifuservp) then
          call advab_rho_cyl(1)
        else
          call advab_cyl(1)
        endif
      endif
!      if (ifmvbd.and..not.ifchar)               call admeshv
      if (ifmvbd .and..not.ifchar) then
!       Mesh convection term        
        if (ifuservp) then
          call advab_rho_cyl(2)
        else
          call advab_cyl(2)
        endif
      endif

      if (iftran)                               call makeabf
      if ((iftran.and..not.ifchar).or.
     $    (iftran.and..not.ifnav.and.ifchar))   call makebdf
      if (ifnav.and.ifchar)                     call advchar
      
      tmakf=tmakf+(dnekclock()-etime1)

      return
      end
!---------------------------------------------------------------------- 
      subroutine advab_rho_cyl(ifld)

!     Eulerian scheme, add convection term to forcing function 
!     at current time step.
!     with variable density        

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MVGEOM'
      include 'MASS'
      include 'TSTEP'

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      real rhoi         ! density inverse
      common /scrcg/ rhoi(lx1,ly1,lz1,lelt)

      integer ntot1

      integer ifld      ! ifld = 1: Fluid velocities as convecting field
                        ! ifld = 2: -(Mesh) velocities as convecting field


      ntot1 = lx1*ly1*lz1*nelv

!     bfx,bfy,bfz get multiplied by vtrans later in makeabf.
!     So I do an inverse rho multiplication here

!     rhoi = 1/\rho      
      call invers2(rhoi,vtrans(1,1,1,1,ifield),ntot1)

      if (ifld.eq.1) then
!       Convecting using Fluid velocity            
        call convect_cylindrical_rho (ta1,vtrans(1,1,1,1,ifield),vx,
     $                                vx,vy,vz)
        call convect_cylindrical_rho (ta2,vtrans(1,1,1,1,ifield),vy,
     $                                vx,vy,vz)
        if (ldim.eq.3) then
          call convect_cylindrical_rho (ta3,vtrans(1,1,1,1,ifield),vz,
     $                                  vx,vy,vz)
        endif  

        call col2(ta1,rhoi,ntot1)
        call sub2 (bfx,ta1,ntot1)

        call col2(ta2,rhoi,ntot1)      
        call sub2 (bfy,ta2,ntot1)
        if (ldim.eq.3) then
          call col2(ta3,rhoi,ntot1)
          call sub2 (bfz,ta3,ntot1)
        endif

!       Additional terms in cylindrical formulation
!       Division by R gets cancelled by the multiplication by R
!       from the Jacobian
        if (ldim.eq.3) then      
          call dealias_rho_uv(ta2,vtrans(1,1,1,1,ifield),vz,vz)
          call col2(ta2,rhoi,ntot1)      
          call add2(bfy,ta2,ntot1)
          
          call dealias_rho_uv(ta3,vtrans(1,1,1,1,ifield),vy,vz)
          call col2(ta3,rhoi,ntot1)      
          call sub2(bfz,ta3,ntot1)
        endif

      elseif (ifld.eq.2) then
!       Convecting using Mesh velocity            
        call convect_cylindrical_rho (ta1,vtrans(1,1,1,1,ifield),vx,
     $                                wx,wy,wz)
        call convect_cylindrical_rho (ta2,vtrans(1,1,1,1,ifield),vy,
     $                                wx,wy,wz)
        if (ldim.eq.3) then
          call convect_cylindrical_rho (ta3,vtrans(1,1,1,1,ifield),vz,
     $                                  wx,wy,wz)
        endif  

!       Sign is opposite that of Fluid convection term        
        call col2(ta1,rhoi,ntot1)
        call add2 (bfx,ta1,ntot1)

        call col2(ta2,rhoi,ntot1)      
        call add2 (bfy,ta2,ntot1)
        if (ldim.eq.3) then
          call col2(ta3,rhoi,ntot1)
          call add2 (bfz,ta3,ntot1)
        endif

!       Additional terms in cylindrical formulation
!       Division by R gets cancelled by the multiplication by R
!       from the Jacobian
        if (ldim.eq.3) then      
          call dealias_rho_uv(ta2,vtrans(1,1,1,1,ifield),wz,vz)
          call col2(ta2,rhoi,ntot1)      
          call sub2(bfy,ta2,ntot1)
          
          call dealias_rho_uv(ta3,vtrans(1,1,1,1,ifield),wz,vy)
          call col2(ta3,rhoi,ntot1)      
          call add2(bfz,ta3,ntot1)
        endif
      else        
        if (nio.eq.0) then
          write(6,*) 'Unknown ifld=',ifld
          write(6,*) 'Exitting in advab_rho_cyl'
        endif
        call exitt        
      endif        

      return
      end subroutine
!---------------------------------------------------------------------- 

      subroutine advab_cyl(ifld)

!     Eulerian scheme, add convection term to forcing function 
!     at current time step.

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MVGEOM'
      include 'MASS'
      include 'TSTEP'

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      integer ntot1

      integer ifld      ! ifld = 1: Fluid velocities as convecting field
                        ! ifld = 2: -(Mesh) velocities as convecting field

      ntot1 = lx1*ly1*lz1*nelv

      if (ifld.eq.1) then
!       Convecting using Fluid velocity            
        call convect_cylindrical (ta1,vx,vx,vy,vz)
        call convect_cylindrical (ta2,vy,vx,vy,vz)
        if (ldim.eq.3) then
          call convect_cylindrical (ta3,vz,vx,vy,vz)
        endif  

        call sub2 (bfx,ta1,ntot1)
        call sub2 (bfy,ta2,ntot1)
        if (ldim.eq.3) then
          call sub2 (bfz,ta3,ntot1)
        endif

!       Additional terms in cylindrical formulation
!       Division by R gets cancelled by the multiplication by R
!       from the Jacobian
        if (ldim.eq.3) then      
          call dealias_uv(ta2,vz,vz)
          call add2(bfy,ta2,ntot1)
          
          call dealias_uv(ta3,vy,vz)
          call sub2(bfz,ta3,ntot1)
        endif
      elseif (ifld.eq.2) then
!       Convecting using Mesh velocity            
        call convect_cylindrical (ta1,vx,wx,wy,wz)
        call convect_cylindrical (ta2,vy,wx,wy,wz)
        if (ldim.eq.3) then
          call convect_cylindrical (ta3,vz,wx,wy,wz)
        endif  

        call add2 (bfx,ta1,ntot1)
        call add2 (bfy,ta2,ntot1)
        if (ldim.eq.3) then
          call add2 (bfz,ta3,ntot1)
        endif

!       Additional terms in cylindrical formulation
!       Division by R gets cancelled by the multiplication by R
!       from the Jacobian
        if (ldim.eq.3) then      
          call dealias_uv(ta2,wz,vz)
          call sub2(bfy,ta2,ntot1)
          
          call dealias_uv(ta3,vz,vy)
          call add2(bfz,ta3,ntot1)
        endif

      else
        if (nio.eq.0) then
          write(6,*) 'Unknown ifld=',ifld
          write(6,*) 'Exitting in advab_cyl'
        endif
        call exitt        

      endif

      return
      end subroutine
!---------------------------------------------------------------------- 

      subroutine incomprn_cyl (ux,uy,uz,up)
c
c     Project U onto the closest incompressible field
c
c     Input:  U     := (ux,uy,uz)
c
c     Output: updated values of U, iproj, proj; and
c             up    := pressure currection req'd to impose div U = 0
c
c
c     Dependencies: ifield ==> which "density" (vtrans) is used.
c
c     Notes  1.  up is _not_ scaled by bd(1)/dt.  This should be done
c                external to incompr().
c
c            2.  up accounts _only_ for the perturbation pressure,
c                not the current pressure derived from extrapolation.

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'    ! vtrans
      include 'MASS'
      include 'TSTEP'
      include 'CTIMER'

      real ux(1),uy(1),uz(1),up(1)

      real w1,w2,w3
      real dv1,dv2,dv3,dp
      common /scrns/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
     $ ,             dv1   (lx1,ly1,lz1,lelv)
     $ ,             dv2   (lx1,ly1,lz1,lelv)
     $ ,             dv3   (lx1,ly1,lz1,lelv)
     $ ,             dp    (lx2,ly2,lz2,lelv)

      real h1,h2,h2inv
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)

      integer nset
      parameter(nset = 1 + lbelv/lelv)

      real pset
      common /orthov/ pset(lx2*ly2*lz2*lelv*mxprev,nset)
      integer nprv
      common /orthbi/ nprv(2)
      logical ifprjp

      integer ntot1,ntot2,intype,istart
      real dtbd,bdti,const,scaledt,scaledi,dtb
      integer i

      ifprjp=.false.    ! Project out previous pressure solutions?
      istart=param(95)  
      if (istep.ge.istart.and.istart.ne.0) ifprjp=.true.

      if (icalld.eq.0) tpres=0.0
      icalld = icalld+1
      npres  = icalld
      etime1 = dnekclock()

      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv

      intype = 1

      if (intype.eq.-1) then
        call sethlm(h1,h2,intype)
        call invers2 (h2inv,h2,ntot1)
      elseif (intype.eq.1) then
        call rzero   (h1,ntot1)
        call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
        call invers2 (h2inv,h2,ntot1)
      else
        if (nio.eq.0) write(6,*) 'Unknown intype', intype
        if (nio.eq.0) write(6,*) 'Exitting in incomprn_test'
        call exitt 
      endif

      call opdiv_cyl(dp,ux,uy,uz)
!      call opdiv_rho_cyl(dp,ux,uy,uz)

      if (intype.eq.1) then
        bdti = -bd(1)/dt
        call cmult   (dp,bdti,ntot2)
      else
        call cmult   (dp,-1.0,ntot2)
      endif
      call add2col2(dp,bm2,usrdiv,ntot2) ! User-defined divergence.

!      call ortho   (dp)

      i = 1 + ifield/ifldmhd
      if (ifprjp)   call setrhsp  (dp,h1,h2,h2inv,pset(1,i),nprv(i))

      if (intype.eq.1) then 
        scaledt = dt/bd(1)
        scaledi = 1./scaledt
        call cmult(dp,scaledt,ntot2)        ! scale for tol
        call esolver_cyl(dp,h1,h2,h2inv,intype)
        call cmult(dp,scaledi,ntot2)
      else
        call esolver_cyl(dp,h1,h2,h2inv,intype)
      endif 
      if (ifprjp) call gensolnp (dp,h1,h2,h2inv,pset(1,i),nprv(i))

      call add2(up,dp,ntot2)

      call opgradt_cyl(w1 ,w2 ,w3 ,dp)
      if (intype.eq.-1) then
        call ophinv (dv1,dv2,dv3,w1,w2,w3,h1,h2,tolhs,nmxv)
      elseif (intype.eq.1) then
        call opbinv (dv1,dv2,dv3,w1,w2,w3,h2inv)
      endif
      
      if (intype.eq.1) then
        dtb  = dt/bd(1)
        call opadd2cm (ux ,uy ,uz ,dv1,dv2,dv3, dtb )
      else
        call opadd2 (ux,uy,uz,dv1,dv2,dv3)
      endif

      if (ifmhd)  call chkptol      ! to avoid repetition

      tpres=tpres+(dnekclock()-etime1)

      return
      end
c-----------------------------------------------------------------------

