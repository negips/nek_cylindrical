c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)

      implicit none  
  
      include 'SIZE'
      include 'PARALLEL'
!      include 'TOTAL'
      include 'INPUT'
      include 'TSTEP'
      include 'NEKUSE'

      integer e,ix,iy,iz,ieg

      real rhog         ! density of air
      real rhow         ! density of water
      real mug          ! dynamic viscosity of air
      real muw          ! dynamic viscosity of water
      real nug          ! kinematic viscosity of air
      real nuw          ! kinematic viscosity of water
      real alpha

      real setvp

      real distn, eps

!     densities      
      rhog   = uparam(1)         ! 1.2061e-3
      rhow   = uparam(2)         ! 1.0

!     viscosities      
      mug    = uparam(3)         ! 1.5052e-6
      muw    = uparam(4)         ! 8.3e-5

      nug    = mug/rhog
      nuw    = muw/rhow

      eps    = 1.0e-1

      if (ifield.eq.1) then
        utrans = setvp(rhog,rhow,temp,eps)
        udiff  = setvp(mug,muw,temp,eps)
      endif  

      if (ifield .eq. 2) then
        e = gllel(ieg)
        utrans = 1.0   
        udiff = 1.0 ! param(8)
      endif
 
      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)

      implicit none        
  
      include 'SIZE'
      include 'NEKUSE'

      integer ix,iy,iz,ieg

      ffx = -3.00
      ffy = 0.0
      ffz = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      qvol =  0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'SOLN'
      include 'MVGEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'MASS'
      include 'PARALLEL'      ! nelgv

      include 'F3D'
      include 'FS_ALE'
      include 'TEST'
      include 'DOMAIN'

      integer lt
      parameter (lt=lx1*ly1*lz1*lelt)

      integer i,j,k,e,n,n2

      integer igeom
      character cb*3
      integer ie,iface,nfaces

      real pos(lt)
      real wght(lt)
      real x,y

      integer lxx,levb
      parameter(lxx=lx1*lx1, levb=lelv+lbelv)
      real df,sr,ss,st
      common /fastd/  df(lx1*ly1*lz1,levb)
     $             ,  sr(lxx*2,levb),ss(lxx*2,levb),st(lxx*2,levb)

!     coarse grid
      real a,acopy                                                
      common /h1crsa/ a(lcr*lcr*lelv)      ! Could use some common
     $              , acopy(lcr*lcr*lelv)  ! block for this

      real h1,h2,h2inv
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)

      integer nit

      integer lxyz
      parameter (lxyz=(lx1+2)*(ly1+2)*(lz1+2))
      integer*8 glo_num
      common /c_is1/ glo_num(lxyz*lelv)
      real vertex
      common /ivrtx/ vertex ((2**ldim)*lelt)


      n  = lx1*ly1*lz1*nelv
      n2 = lx2*ly2*lz2*nelv

      if (istep.eq.0) then

        if (param(95).gt.0) then
          param(95) = 50        ! start of projections
        endif

        call rone(vtrans(1,1,1,1,2),n)
        call rone(vdiff(1,1,1,1,2),n)

        call phi0(t)

        ifield = 1
        call vprops
        ifield = 2
        call vprops

        ifheat = .true.

        call frame_start

        ifto = .true.

        call outpost(v1mask,v2mask,v3mask,pr,tmask,'msk')

!       Reynolds number of the field
        do i=1,n
          t(i,1,1,1,2) = vtrans(i,1,1,1,1)*1.0*7.3/vdiff(i,1,1,1,1)
        enddo       
        call outpost2(vtrans,vdiff,t(1,1,1,1,2),pr,
     $                vdiff(1,1,1,1,2),1,'vis')

!       Preconditioner
        param(42)=uparam(8)       ! 0: GMRES (nonsymmetric), 1: PCG w/o weights
        param(43)=uparam(9)       ! 0: Additive multilevel (param 42=0), 1: Original 2 level
        param(44)=uparam(10)      ! 0: E based Schwartz (FEM), 1: A based Schwartz
       

12    format(A4,2x,12(E12.5,2x))


        ifield = 1

        time = 1.0
        istep = 1

        if (ifpgll) then

          call get_vert()
          call setupds(pgs_handle,lx2,ly2,lz2,nelv,nelgv,vertex,glo_num)
      
!          call exitt
        endif

!       opgradt_cyl        
        call rone(pr,n2)
        call opgradt_cyl(vx,vy,vz,pr)

!       opdiv_cyl        
        call rone(vx,n)
        call rone(vy,n)
        call rone(vz,n)
        call opdiv_cyl(pr,vx,vy,vz)

        call rone(pr,n2)
        call col2(pr,bm2,n2)
        call outpost(vx,vy,vz,pr,t,'cyl')
        call exitt

        call theta_outpost()

        call jacobi_davidson_e()

        call exitt

        call laplace_test()

        call exitt

        call gen_mapping_mvb()

      endif 

      call frame_monitor

      call chkpt_main

      ifto = .true.

      call copy(t(1,1,1,1,2),vz,n)
        
!      if (fs_iffs) call fs_mvmesh_linear()

      if (istep.eq.nsteps.or.lastep.eq.1) then
        call frame_end
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)

      implicit none

      include 'SIZE'
      include 'NEKUSE'
      include 'GEOM'

      integer ix,iy,iz,iside,ieg
      real pi

      integer jp
      common /ppointr/ jp

      real rmid

      real glmin,glmax
      integer n

      real rad1,rad2,omega1,omega2
      real a1,a2
      common /cylindrical/ rad1,rad2,omega1,omega2,a1,a2

      ux   = 0.0
      uy   = 0.0
      uz   = 0.0

      rmid = (rad1+rad2)/2
      if (y.lt.(rmid)) uz = omega1*rad1

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)

      implicit none

      include 'SIZE'
      include 'INPUT'         ! if3d
      include 'PARALLEL'
      include 'NEKUSE'
      include 'GEOM'

      include 'F3D'

      integer ix,iy,iz,ieg
      real pi

      integer jp
      common /ppointr/ jp

      real fcoeff(3)
      real xl(3)
      real mth_ran_dst

      logical ifcouette
      logical ifpoiseuille
      logical iftaylor
      logical iftestmvb

      real rad1,rad2,omega1,omega2
      real a1,a2
      common /cylindrical/ rad1,rad2,omega1,omega2,a1,a2

      real rmid,y0,z0
      real rad

      ifcouette         = .false.
      ifpoiseuille      = .false.
      iftaylor          = .true.
      iftestmvb         = .false.

      pi = 4.0*atan(1.0)


      if (ifpoiseuille) then
        ux = 1.0 - y**2
        uy = 0.
        uz = 0.0 + 0.0
      elseif (ifcouette) then
        ux = 0.0 + 1.0*y
        uy = 0.
        uz = 0.0 + 0.0
      elseif (iftaylor) then
        ux = 0.0 ! sin(pi*(y-1.0)/0.8)
        uy = 0.
        uz = a1*y + a2/y
      elseif (iftestmvb) then
        rmid = (rad1+rad2)/2.0
        ux   = -uparam(1)*exp(-((y-rmid)/0.25)**2)
        uy   = -0.1*uparam(1)*exp(-((y-rmid)/0.25)**2)
        uz   = 0.1*(a1*y + a2/y)
      endif  


      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat   ! This routine to modify element vertices

      implicit none  
  
      include 'SIZE'      ! _before_ mesh is generated, which 
      include 'INPUT'
      include 'GEOM'
      include 'TSTEP'
      include 'PARALLEL'
!      include 'TOTAL'     ! guarantees GLL mapping of mesh.

!      ifaxis = .true.   ! just for initialization
      param(42)=0       ! 0: GMRES (nonsymmetric), 1: PCG w/o weights
      param(43)=1       ! 0: Additive multilevel (param 42=0), 1: Original 2 level
      param(44)=0       ! 0: E based Schwartz (FEM), 1: A based Schwartz



      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2   ! This routine to modify mesh coordinates

      implicit none

      include 'SIZE'
      include 'INPUT'      ! cbc
      include 'PARALLEL'
      include 'GEOM'

      integer iel,ifc
      integer n

      if (if3d.and..not.ifcyclic) then
        n = lx1*ly1*lz1*nelv
!        call cmult(zm1,0.20,n)
      endif  

      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3

      implicit none        

      include 'SIZE'
      include 'SOLN'    ! tmult
      include 'INPUT'
      include 'GEOM'

      integer i,n
      real rad1,rad2,omega1,omega2
      real a1,a2
      common /cylindrical/ rad1,rad2,omega1,omega2,a1,a2

      real x,y,z

      real radius
      common /scrcg/ radius(lx1,ly1,lz1,lelt)

      real glmin,glmax

      n = lx1*ly1*lz1*nelv
      rad1 = glmin(ym1,n)
      rad2 = glmax(ym1,n)
      omega1 = 1.0/rad1
      omega2 = 0.0
      a1 = (omega2*(rad2**2) - omega1*rad1*rad1)/(rad2**2 - rad1**2)
      a2 = (omega1 - omega2)*(rad1**2)*(rad2**2)/(rad2**2 - rad1**2)

      if (nio.eq.0) write(6,*) 'Cylindrical Params:', rad1,rad2,a1,a2

      return
      end
c-----------------------------------------------------------------------
      subroutine phi0(phi)

      implicit none

      include 'SIZE'
      include 'GEOM'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real phi(lv)
     
      integer i,n
      real x,y
      real r,r0

      real wallx(2)
      real rady(2)
      common /taylor_geom/ wallx,rady

      real pi

      n=lx1*ly1*lz1*nelv
      
      pi = 4.0*atan(1.0)
!      r0 = (rady(1) + rady(2))*0.5
      r0 = 1.00
      do i=1,n
        x = xm1(i,1,1,1)
        y = ym1(i,1,1,1)
        r = sqrt(x**2 + y**2)
        phi(i) = x - r0        ! signed distance from interface
!        phi(i) = x - (r0+ 0.025*sin(2.0*pi*(y-1.0)/0.8))        ! wavy interface 
      enddo  

      return
      end subroutine phi0
!---------------------------------------------------------------------- 

      real function heavyside(phi,eps)

      real phi,eps
      real pi

      pi = 4.0*atan(1.0) 

      if (phi.lt.-eps) then
        heavyside = 0.0
      elseif (phi.gt.eps) then
        heavyside = 1.0
      else
        heavyside = 0.5*(1.0 + phi/eps + 1.0/pi*sin(pi*phi/eps))
!        heavyside = 1.0*(1.0 + phi/eps + 0.0/pi*sin(pi*phi/eps))
      endif  

      return
      end function heavyside
!---------------------------------------------------------------------- 

      real function setvp(vg,vl,phi,eps)

      real vg     ! gas property
      real vl     ! liquid property
      real phi    ! disgned distance function
      real eps    
      real heavyside    ! heavyside function

      heavy = heavyside(phi,eps)
      setvp = vg + (vl - vg)*heavy      

      return
      end function 
!---------------------------------------------------------------------- 

      subroutine heavydist(phi,eps,off)

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'SOLN'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real phi(lv)
      real phioff
      real heavyside          ! function
      
      integer i,j,k,e,n
      real hd
      real eps
      real off

      n=lx1*ly1*lz1*nelv

      do i=1,n
        phioff = phi(i)+off
        hd = heavyside(phioff,eps)
        phi(i) = hd
      enddo  

      return
      end subroutine heavydist
!---------------------------------------------------------------------- 

      subroutine theta_outpost()

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'SOLN'

      integer i,j,k,e
      integer n

      real x,y,z,theta,phi
      real pi

      n = lx1*ly1*lz1*nelv
      pi = 4.0*atan(1.0)

      do i=1,n
        x = xm1(i,1,1,1)
        y = ym1(i,1,1,1)
        z = zm1(i,1,1,1)

        phi   = atan2(y,z)*180/pi
        theta = atan2(z,y)*180/pi

        vx(i,1,1,1) = theta
        vy(i,1,1,1) = phi

      enddo

      call outpost(vx,vy,vz,pr,t,'phi')

      return
      end subroutine        
!---------------------------------------------------------------------- 





c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
