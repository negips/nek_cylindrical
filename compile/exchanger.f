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


        call setup_exchange_m2

        call exitt

12    format(A4,2x,12(E12.5,2x))


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

      integer iel,ifc


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

      subroutine extend_pr_test

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'TEST'

      integer e,i,j,k
      real p0,el
      integer n,n2
      integer gl
      integer ix,iy,iz,iz1

      n  = lx1*ly1*lz1*nelv
      n2 = lx2*ly2*lz2*nelv

      call rzero(tmp4,n2)
      call rzero(tmp8,n2)
      call rzero(tmp12,n2)

      do e=1,nelv
        gl = lglel(e) 
        do i=1,lx2*ly2*lz2
          tmp4(i,1,1,e) = gl + 0.0
        enddo
      enddo  

      call fill_interior(tmp1,tmp4)
      call dface_ext    (tmp1)
      call dssum        (tmp1,lx1,ly1,lz1)
      call dface_add1si (tmp1,-1.)

      call copy(tmp2,tmp1,n)
     

c     Fill interiors
      call fill_interior(tmp1,tmp4)
      call dface_ext    (tmp1)
      call crnr_ext     (tmp1)

      call dssum        (tmp1,lx1,ly1,lz1)
      call dface_add1si (tmp1,-1.)
      call crnr_corr    (tmp1,tmp2)

      call outpost(tmp1,tmp2,tmp3,tmp4,tmp5,'ext')


      return
      end subroutine
!---------------------------------------------------------------------- 

      subroutine test_semhat

      implicit none

      include 'SIZE'
      include 'INPUT'         ! if3d
      include 'GEOM'
      include 'MYSEMHAT'
      include 'TEST'
      include 'WZ'
      include 'DXYZ'

      integer i,j,k,e
      integer n,n2
      integer nr

!     These are populated by Swap Length      
      real l
      common /swaplengths/ l(lx1,ly1,lz1,lelv)
      real lr ,ls ,lt   ! not used by swap lengths
      real llr(lelt)    ! length of left element along "r" 
      real lls(lelt)    ! length of left element along "s"
      real llt(lelt)    ! length of left element along "t"
      real lmr(lelt)    ! length of this element along "r"
      real lms(lelt)    ! length of this element along "s"
      real lmt(lelt)    ! length of this element along "t"
      real lrr(lelt)    ! length of right element along "r"
      real lrs(lelt)    ! length of right element along "s"
      real lrt(lelt)    ! length of right element along "t"
      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr,lls,llt
     $              , lmr,lms,lmt
     $              , lrr,lrs,lrt

!     r,s,t coordinates of the 1D spectral elements.      
      real lr1(lx1,lelt),ls1(ly1,lelt),lt1(lz1,lelt)
!     Derivative mapping
      real drdx(lx1,lelt),dsdy(ly1,lelt),dtdz(lz1,lelt)
      real derivx(lx1,ly1,lelv),derivy(lx1,ly1,lelv)
      
      real sc           ! scale
      real dd           ! delta for the end points

      logical ifext


!     ah          = Laplacian
!     bh          = diagonal mass matrix
!     ch          = convection operator b*d
!     dh          = derivative matrix
!     dph         = derivative matrix (Same as dh for now)
!     jph         = interpolation matrix 
!     z           = GLL points
!
!     zglhat      = GL points
!     bgl         = diagonal mass matrix on GL
!     dgl         = derivative matrix,    mapping from velocity nodes to pressure
!     jgl         = interpolation matrix, mapping from velocity nodes to pressure
!
!     nr          = polynomial degree (velocity space)
!     wh          = Work array
      nr = lx1-1
      call mysemhat(ah,bh,ch,dh,zh,dph,dpht,jph,bgl,
     $              zglhat,dgl,jgl,nr,wh)

!      do j=1,lx1
!        k = (j-1)*lx1
!        write(6,12) 'dgll', (dph(i+k),i=1,lx1)
!      enddo  
!12    format(A4,2x,6(E12.5,2x))
     
      n  = lx1*ly1*lz1*nelv
      n2 = lx2*ly2*lz2*nelv

      call swap_lengths ! populate llr,lmr,lmt...

      ifext = .true.

!     Build 1D x,y,z      
      do e=1,nelv
        call copy(lr1(1,e),zgm1(1,1),lx1)
        sc = lmr(e)/2.0
        call cmult(lr1(1,e),sc,lx1)
        if (ifext) then
!         Add left extension        
          sc = zgm1(lx1-1,1)-zgm1(lx1,1)
          dd = sc*llr(e)/2.0
          lr1(1,e) = lr1(1,e)+dd
!         Add right extension        
          sc = zgm1(2,1)-zgm1(1,1)
          dd = sc*lrr(e)/2.0
          lr1(lx1,e) = lr1(lx1,e)+dd
        endif  

        call copy(ls1(1,e),zgm1(1,2),ly1)
        sc = lms(e)/2.0
        call cmult(ls1(1,e),sc,lx1)
        if (ifext) then
!         Add left extension        
          sc = zgm1(lx1-1,2)-zgm1(lx1,2)
          dd = sc*lls(e)/2.0
          ls1(1,e) = ls1(1,e)+dd
!         Add right extension        
          sc = zgm1(2,2)-zgm1(1,2)
          dd = sc*lrs(e)/2.0
          ls1(lx1,e) = ls1(lx1,e)+dd
        endif  

        if (if3d) then
          call copy(lt1(1,e),zgm1(1,3),lz1)
          sc = lmt(e)/2.0
          call cmult(lt1(1,e),sc,lx1)
          if (ifext) then
!           Add left extension        
            sc = zgm1(lx1-1,3)-zgm1(lx1,3)
            dd = sc*llt(e)/2.0
            lt1(1,e) = lt1(1,e)+dd
!           Add right extension        
            sc = zgm1(2,3)-zgm1(1,3)
            dd = sc*lrt(e)/2.0
            lt1(lx1,e) = lt1(lx1,e)+dd
          endif  
        endif
      enddo 


!     Geometric factors/Jacobians      
      do e=1,nelv
        call mxm(dph,lx1,lr1(1,e),lx1,drdx(1,e),1)
        call mxm(dph,lx1,ls1(1,e),lx1,dsdy(1,e),1)
        if (if3d) call mxm(dph,lx1,lt1(1,e),lx1,dtdz(1,e),1)
      enddo
      call invcol1(drdx,lx1*nelv)
      call invcol1(dsdy,lx1*nelv)
      if (if3d) call invcol1(dtdz,lx1*nelv)

!     Output Geometric factors 
      do e=1,nelv
        do k=1,lz1
        do j=1,ly1
        do i=1,lx1
          tmp1(i,j,k,e) = drdx(i,e) ! lr1(i,e)
          tmp2(i,j,k,e) = dsdy(j,e) ! ls1(j,e)
          if (if3d) tmp3(i,j,k,e) = dtdz(k,e) ! lt1(k,e)
          tmp5(i,j,1,e) = lls(e)
        enddo  
        enddo
        enddo
      enddo
      call outpost(tmp1,tmp2,tmp3,tmp4,tmp5,'tmp')

!     Derivatives without geometric factors
      call opzero(tmp1,tmp2,tmp3)
      call copy(tmp3,xm1,n)
      call copy(tmp5,ym1,n)

      do e=1,nelv
        call mxm(dph,lx1,tmp3(1,1,1,e),lx1,tmp1(1,1,1,e),lx1)
        call mxm(tmp5(1,1,1,e),lx1,dpht,lx1,tmp2(1,1,1,e),lx1)
      enddo  

      call outpost(tmp1,tmp2,tmp3,tmp4,tmp5,'tmp')

!     Derivatives with geometric factors
      do e=1,nelv
        do i=1,lx1
        do j=1,lx1
          k = (j-1)*lx1
          derivx(i,j,e) = dph(i+k)*drdx(i,e)
          derivy(i,j,e) = dpht(i+k)*dsdy(j,e)   ! transposed
        enddo
        enddo
      enddo  

      call extract_interior(tmp4,xm1)
      call extract_interior(tmp8,ym1)
      call exchange_m2(tmp1,tmp4,tmp5,tmp6,tmp7)
      call exchange_m2(tmp2,tmp8,tmp5,tmp6,tmp7)

      call outpost(tmp1,tmp2,tmp3,tmp4,tmp5,'tmp')

      do e=1,nelv
        call mxm(derivx(1,1,e),lx1,tmp1(1,1,1,e),lx1,tmp3(1,1,1,e),lx1)
        call mxm(tmp2(1,1,1,e),lx1,derivy(1,1,e),lx1,tmp5(1,1,1,e),lx1)
      enddo

      call outpost(tmp3,tmp5,tmp6,tmp4,tmp5,'tmp')


      return
      end subroutine test_semhat
!---------------------------------------------------------------------- 

      subroutine test_semhat2

      implicit none

      include 'SIZE'
      include 'INPUT'         ! if3d
      include 'GEOM'
      include 'MYSEMHAT'
      include 'TEST'
      include 'WZ'
      include 'DXYZ'

      integer i,j,k,e
      integer n,n2
      integer nr

!     These are populated by Swap Length      
      real l
      common /swaplengths/ l(lx1,ly1,lz1,lelv)
      real lr ,ls ,lt   ! not used by swap lengths
      real llr(lelt)    ! length of left element along "r" 
      real lls(lelt)    ! length of left element along "s"
      real llt(lelt)    ! length of left element along "t"
      real lmr(lelt)    ! length of this element along "r"
      real lms(lelt)    ! length of this element along "s"
      real lmt(lelt)    ! length of this element along "t"
      real lrr(lelt)    ! length of right element along "r"
      real lrs(lelt)    ! length of right element along "s"
      real lrt(lelt)    ! length of right element along "t"
      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr,lls,llt
     $              , lmr,lms,lmt
     $              , lrr,lrs,lrt

!     r,s,t coordinates of the 1D spectral elements.      
      real lr1(lx1,lelt),ls1(ly1,lelt),lt1(lz1,lelt)
!     Derivative mapping
      real drdx(lx1,lelt),dsdy(ly1,lelt),dtdz(lz1,lelt)
      real derivx(lx1,ly1,lelv),derivy(lx1,ly1,lelv)
!     Laplacians      
      real laplx(lx1,lx1,lelv),laply(ly1,ly1,lelv),laplz(lz1,lz1,lelv)
      
      real sc           ! scale
      real dd           ! delta for the end points

      logical ifext

      real tmpx(lx1,2)
      real pi
      real x,y
      character*32 str


!     ah          = Laplacian
!     bh          = diagonal mass matrix
!     ch          = convection operator b*d
!     dh          = derivative matrix
!     dph         = derivative matrix (Same as dh for now)
!     jph         = interpolation matrix 
!     z           = GLL points
!
!     zglhat      = GL points
!     bgl         = diagonal mass matrix on GL
!     dgl         = derivative matrix,    mapping from velocity nodes to pressure
!     jgl         = interpolation matrix, mapping from velocity nodes to pressure
!
!     nr          = polynomial degree (velocity space)
!     wh          = Work array
      nr = lx1-1
      call mysemhat(ah,bh,ch,dh,zh,dph,dpht,jph,bgl,
     $              zglhat,dgl,jgl,nr,wh)


!      do j=1,lx1
!        k = (j-1)*lx1
!        write(6,12) 'dgll', (dph(i+k),i=1,lx1)
!      enddo
      call blank(str,32)
      write(str,'(A8,I2,A11)') '(A4,2x,I',lx1,'(E10.4,2x))'
      write(6,*) str
12    format(A4,2x,12(E22.16,2x))
     
      n  = lx1*ly1*lz1*nelv
      n2 = lx2*ly2*lz2*nelv

!     populate llr,lmr,lmt...
!     using Mesh2 coordinate positions      
      call swap_lengths_mesh2

      ifext = .true.

!     Build 1D x,y,z based on Mesh2
      do e=1,nelv
        lr1(1,e)   = zgm1(1,1)
        lr1(lx1,e) = zgm1(lx1,1)
        call copy(lr1(2,e),zgm2(1,1),lx2)
        sc = lmr(e)/2.0
        call cmult(lr1(1,e),sc,lx1)
        if (ifext) then
!         Add left extension        
          sc = zgm2(lx2,1)-zgm1(lx1,1)
          dd = sc*llr(e)/2.0
          lr1(1,e) = lr1(1,e)+dd
!         Add right extension        
          sc = zgm2(1,1)-zgm1(1,1)
          dd = sc*lrr(e)/2.0
          lr1(lx1,e) = lr1(lx1,e)+dd
        endif  

        ls1(1,e)   = zgm1(1,2)
        ls1(ly1,e) = zgm1(ly1,2)
        call copy(ls1(2,e),zgm2(1,2),ly2)
        sc = lms(e)/2.0
        call cmult(ls1(1,e),sc,ly1)
        if (ifext) then
!         Add left extension        
          sc = zgm2(ly2,2)-zgm1(ly1,2)
          dd = sc*lls(e)/2.0
          ls1(1,e) = ls1(1,e)+dd
!         Add right extension        
          sc = zgm2(1,2)-zgm1(1,2)
          dd = sc*lrs(e)/2.0
          ls1(ly1,e) = ls1(ly1,e)+dd
        endif  

        if (if3d) then
          lt1(1,e)   = zgm1(1,3)
          lt1(lz1,e) = zgm1(lz1,3)
          call copy(lt1(2,e),zgm2(1,3),lz2)
          sc = lmt(e)/2.0
          call cmult(lt1(1,e),sc,lz1)
          if (ifext) then
!           Add left extension        
            sc = zgm2(lz2,3)-zgm1(lz1,3)
            dd = sc*llt(e)/2.0
            lt1(1,e) = lt1(1,e)+dd
!           Add right extension        
            sc = zgm2(1,3)-zgm1(1,3)
            dd = sc*lrt(e)/2.0
            lt1(lz1,e) = lt1(lz1,e)+dd
          endif  
        endif
      enddo 
!     Output 1D elements 
      do e=1,nelv
        do k=1,lz1
        do j=1,ly1
        do i=1,lx1
          tmp1(i,j,k,e) = lr1(i,e)
          tmp2(i,j,k,e) = ls1(j,e)
          if (if3d) tmp3(i,j,k,e) = lt1(k,e)
        enddo  
        enddo
        enddo
      enddo
      call outpost(tmp1,tmp2,tmp3,tmp4,tmp5,'tmp')    ! 1D elements


!     Geometric factors/Jacobians      
      do e=1,nelv
        call mxm(dph,lx1,lr1(1,e),lx1,drdx(1,e),1)
        call mxm(dph,ly1,ls1(1,e),ly1,dsdy(1,e),1)
        if (if3d) call mxm(dph,lz1,lt1(1,e),lz1,dtdz(1,e),1)
      enddo
      call invcol1(drdx,lx1*nelv)
      call invcol1(dsdy,ly1*nelv)
      if (if3d) call invcol1(dtdz,lz1*nelv)

!     Output Geometric factors 
      do e=1,nelv
        do k=1,lz1
        do j=1,ly1
        do i=1,lx1
          tmp1(i,j,k,e) = drdx(i,e) ! lr1(i,e)
          tmp2(i,j,k,e) = dsdy(j,e) ! ls1(j,e)
          if (if3d) tmp3(i,j,k,e) = dtdz(k,e) ! lt1(k,e)
          tmp5(i,j,1,e) = lls(e)
        enddo  
        enddo
        enddo
      enddo
      call outpost(tmp1,tmp2,tmp3,tmp4,tmp5,'tmp')    ! Geom Factors

!     Derivatives with geometric factors
      do e=1,nelv
        do i=1,lx1
        do j=1,lx1
          k = (j-1)*lx1
          derivx(i,j,e) = dph(i+k)*drdx(i,e)
          derivy(i,j,e) = dpht(i+k)*dsdy(j,e)   ! transposed
        enddo
        enddo
      enddo  

      pi = 4.0*atan(1.0)
      do i=1,n2
        x = xm2(i,1,1,1)
        y = ym2(i,1,1,1)
        tmp4(i,1,1,1) = sin(2*pi*(x-1.0)/3.0)
        tmp8(i,1,1,1) = sin(2*pi*(y-1.0)/3.0)
      enddo
!      call copy(tmp4,xm2,n2)
      call copy(tmp8,ym2,n2) 

      call exchange_m2(tmp1,tmp4,tmp5,tmp6,tmp7)
      call exchange_m2(tmp2,tmp8,tmp5,tmp6,tmp7)

      call outpost(tmp1,tmp2,tmp3,tmp4,tmp5,'tmp')   ! the field to derive

      do e=1,nelv
        call mxm(derivx(1,1,e),lx1,tmp1(1,1,1,e),lx1,tmp3(1,1,1,e),lx1)
        call mxm(tmp2(1,1,1,e),lx1,derivy(1,1,e),lx1,tmp5(1,1,1,e),lx1)
      enddo

      call outpost(tmp3,tmp5,tmp6,tmp4,tmp5,'tmp') ! actual derivative


      return
      end subroutine test_semhat2
!---------------------------------------------------------------------- 

      subroutine build_local_1D_laplacian(ah,dph,drdx,bm,n)

!     This is the weak Laplacian Matrix
!     ah = DT*B*D 

      implicit none

      integer n         ! lx1/ly1/lz1
      real ah(n,n)      ! Laplacian
      real dph(n,n)     ! Gradient matrix
      real drdx(n) 
      real bm(n)        ! Mass (unweighted)

      integer i,j,k
      real s


      do i=1,n
      do j=1,n
        s = 0.0
        do k=1,n
          s = s + dph(k,i)*drdx(k)*bm(k)*(1.0/drdx(k))*dph(k,j)*drdx(k)
        enddo
        ah(i,j) = s
      enddo
      enddo 


      return
      end subroutine build_local_1D_laplacian
!---------------------------------------------------------------------- 

      subroutine build_local_1D_laplacian_nek(ah,dph,drdx,bm,n)

!     This is the weak Laplacian Matrix
!     ah = DT*B*D
!     This is how nek builds it for the local_fdm_solves.
!     I think the indexing is wrong.        

      implicit none

      integer n         ! lx1/ly1/lz1
      real ah(n,n)      ! Laplacian
      real dph(n,n)     ! Gradient matrix
      real drdx(n) 
      real bm(n)        ! Mass (unweighted)

      integer i,j,k
      real s


      do i=1,n
      do j=1,n
        s = 0.0
        do k=1,n
          s = s + dph(i,k)*drdx(k)*bm(k)*(1.0/drdx(k))*dph(j,k)*drdx(k)
        enddo
        ah(i,j) = s
      enddo
      enddo 


      return
      end subroutine build_local_1D_laplacian_nek
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
