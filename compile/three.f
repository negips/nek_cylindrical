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
       
!        call rone(vtrans,n) 
!        call reset_preconditioner()

        do e=1,nelv
          do i=1,lx1*lx1
            tmp1(i,1,1,e)  = sr(i,e)
            tmp2(i,1,1,e)  = sr(i+lxx,e)
            tmp3(i,1,1,e)  = df(i,e)

            tmp5(i,1,1,e)  = ss(i,e)
            tmp6(i,1,1,e)  = ss(i+lxx,e)

            tmp9(i,1,1,e)  = st(i,e)
            tmp10(i,1,1,e) = st(i+lxx,e)
          enddo
        enddo        

        call rzero(vz,n)

        do i=1,n2
          x = xm2(i,1,1,1)  
          y = ym2(i,1,1,1)
          tmp4(i,1,1,1) = 1.0 ! bm2(i,1,1,1)
        enddo
        call mappr(t,bm2,tmp5,tmp6)
        call copy(pr,tmp4,n2)
        call outpost(tmp1,tmp2,vz,tmp4,tmp3,'prc')

12    format(A4,2x,12(E12.5,2x))


        ifield = 1

!        call rzero(tmp8,n2)
!        call local_solves_fdm(tmp8,tmp4)
!        call outpost(tmp5,tmp6,vz,tmp8,t,'prc')

!!       build local inverses based on new method        
!        call gen_fast_again3(df,sr,ss,st)
!        call local_solves_fdm3(tmp8,tmp4)
!        call outpost(tmp1,tmp2,vz,tmp8,t,'prc')

!        call rzero(tmp12,n2)
!        call crs_solve_l2(tmp12,tmp4)
!        call outpost(tmp9,tmp10,vz,tmp12,t,'prc')


!!       dd_solver2
!        call rone(h2,n)
!        call rzero(h1,n)
!        call invers2(h2inv,h2,n)

!        nit = 1
!        istep = nit
!        call rzero(tmp8,n2)
!        call rzero(tmp12,n2)
!        call dd_solver2(pr,tmp4,h1,h2,h2inv,nit,tmp8,tmp12)
!        call outpost(tmp9,tmp10,vz,pr,t,'prc')
!
!        nit = 10
!        istep = nit
!        call rzero(tmp8,n2)
!        call rzero(tmp12,n2)
!        call dd_solver2(pr,tmp4,h1,h2,h2inv,nit,tmp8,tmp12)
!        call outpost(tmp9,tmp10,vz,pr,t,'prc')
!
!        nit = 100
!        istep = nit
!        call rzero(tmp8,n2)
!        call rzero(tmp12,n2)
!        call dd_solver2(pr,tmp4,h1,h2,h2inv,nit,tmp8,tmp12)
!        call outpost(tmp9,tmp10,vz,pr,t,'prc')

        time = 1.0
        istep = 1
        call laplace_test()
      
        call cdabdtp_check()

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
      call exchange_m2(tmp1,tmp4)
      call exchange_m2(tmp2,tmp8)

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

      call exchange_m2(tmp1,tmp4)
      call exchange_m2(tmp2,tmp8)

      call outpost(tmp1,tmp2,tmp3,tmp4,tmp5,'tmp')   ! the field to derive

      do e=1,nelv
        call mxm(derivx(1,1,e),lx1,tmp1(1,1,1,e),lx1,tmp3(1,1,1,e),lx1)
        call mxm(tmp2(1,1,1,e),lx1,derivy(1,1,e),lx1,tmp5(1,1,1,e),lx1)
      enddo

      call outpost(tmp3,tmp5,tmp6,tmp4,tmp5,'tmp') ! actual derivative


      return
      end subroutine test_semhat2
!---------------------------------------------------------------------- 

      subroutine cdabdtp_check

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'DXYZ'

      integer i,j,k,e
      integer e2

      e2 = 1

      do e=1,e2
        write(6,*) 'RXM1', e
        do i=1,lx1
          write(6,13) 'rxm1', (rxm1(i,j,1,e),j=1,lx1)
        enddo
        write(6,*) ' '
      enddo    
       
      do e=1,e2
        write(6,*) 'RXM2', e
        do i=1,lx2
          write(6,13) 'rxm2', (rxm2(i,j,1,e),j=1,lx2)
        enddo
        write(6,*) ' '
      enddo    

      do e=1,1
        write(6,*) 'DXM12', e
        do i=1,lx2
          write(6,13) 'dxm12', (dxm12(i,j),j=1,lx1)
        enddo
        write(6,*) ' '
      enddo    

      do e=1,1
        write(6,*) 'DXM1', e
        do i=1,lx1
          write(6,13) 'dxm1', (dxm1(i,j),j=1,lx1)
        enddo
        write(6,*) ' '
      enddo    

13    format(A6,2x,16(E12.5,2x))



      return
      end subroutine cdabdtp_check

!----------------------------------------------------------------------       

      subroutine cdtM1p (dtx,x,rm1,sm1,tm1,isd)
C-------------------------------------------------------------
C
C     Compute DT*X (entire field)
C     Evaluated on the M1 Mesh.        
C
C-------------------------------------------------------------

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
      integer n1,n2
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

!         ly12   = ly1*ly2
!         if (ifrzer(e)) then
!            call copy (iym12,iam12,ly12)
!            call copy (dym12,dam12,ly12)
!            call copy (w3m2,w2am2,nxyz2)
!         else
!            call copy (iym12,icm12,ly12)
!            call copy (dym12,dcm12,ly12)
!            call copy (w3m2,w2cm2,nxyz2)
!         endif
       endif
C
C      Collocate with weights
C
       if(ifsplit) then
          if (nio.eq.0) write(6,*) 'CDTM1p not implemented for ifsplit.'
          call exitt

!         call col3 (wx,bm1(1,1,1,e),x(1,e),nxyz1)
!         call invcol2(wx,jacm1(1,1,1,e),nxyz1)
       else
!         if (.not.ifaxis) call col3 (wx,w3m2,x(1,e),nxyz2)

!         prabal. Interpolate x to Mesh 1
          call mxm (ixm21,lx1,x(1,e),lx2,ta1,nyz2)
          call mxm (ta1,lx1,iytm21,ly2,wx,ly1)

!         collocate with weights
!         Jacobian goes away due to inverse jacobian of the dx/dr etc. 
          call col2(wx,w3m1,nxyz1)

!         if (ifaxis) then
!            if (ifrzer(e)) then
!                call col3    (wx,x(1,e),bm2(1,1,1,e),nxyz2)
!                call invcol2 (wx,jacm2(1,1,1,e),nxyz2)
!            else
!                call col3    (wx,w3m2,x(1,e),nxyz2)
!                call col2    (wx,ym2(1,1,1,e),nxyz2)
!            endif
!         endif
       endif
C
       if (ldim.eq.2) then
         if (.not.ifdfrm(e) .and. ifalgn(e)) then

            if (      ifrsxy(e).and.isd.eq.1  .or. 
     $           .not.ifrsxy(e).and.isd.eq.2) then

!              prabal. 
               call col3 (ta1,wx,rm1(1,e),nxyz1)
               call mxm  (dxtm1,lx1,ta1,lx1,dtx(1,e),nyz1)
!               call mxm  (ta2,lx1,iym1,ly2,dtx(1,e),ly1)


!               call col3 (ta1,wx,rm2(1,e),nxyz2)
!               call mxm  (dxtm12,lx1,ta1,lx2,ta2,nyz2)
!               call mxm  (ta2,lx1,iym12,ly2,dtx(1,e),ly1)
            else
!              prabal   
               call col3 (ta1,wx,sm1(1,e),nxyz1)
!               call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)
               call mxm  (ta1,lx1,dym1,ly1,dtx(1,e),ly1)


!               call col3 (ta1,wx,sm2(1,e),nxyz2)
!               call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)
!               call mxm  (ta2,lx1,dym12,ly2,dtx(1,e),ly1)
            endif
         else

            call col3 (ta1,wx,rm1(1,e),nxyz1)
            call mxm  (dxtm1,lx1,ta1,lx1,dtx(1,e),nyz1)

            call col3 (ta1,wx,sm1(1,e),nxyz1)
            call mxm  (ta1,lx1,dym1,ly1,ta2,ly1)
          
            call add2 (dtx(1,e),ta2,nxyz1)

!            call col3 (ta1,wx,rm2(1,e),nxyz2)
!            call mxm  (dxtm12,lx1,ta1,lx2,ta2,nyz2)
!            call mxm  (ta2,lx1,iym12,ly2,dtx(1,e),ly1)
!
!            call col3 (ta1,wx,sm2(1,e),nxyz2)
!            call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)
!            call mxm  (ta2,lx1,dym12,ly2,ta1,ly1)
!
!            call add2 (dtx(1,e),ta1,nxyz1)
         endif

       else
         if (ifsplit) then

          if (nio.eq.0) write(6,*) 'CDTM1p not implemented for ifsplit.'
          call exitt

!            call col3 (ta1,wx,rm2(1,e),nxyz2)
!            call mxm  (dxtm12,lx1,ta1,lx2,dtx(1,e),nyz2)
!            call col3 (ta1,wx,sm2(1,e),nxyz2)
!            i1 = 1
!            i2 = 1
!            do iz=1,lz2
!               call mxm  (ta1(i2),lx1,dym12,ly2,ta2(i1),ly1)
!               i1 = i1 + n1
!               i2 = i2 + n2
!            enddo
!            call add2 (dtx(1,e),ta2,nxyz1)
!            call col3 (ta1,wx,tm2(1,e),nxyz2)
!            call mxm  (ta1,nxy1,dzm12,lz2,ta2,lz1)
!            call add2 (dtx(1,e),ta2,nxyz1)
!
!         else
!
!            call col3 (ta1,wx,rm2(1,e),nxyz2)
!            call mxm  (dxtm12,lx1,ta1,lx2,ta2,nyz2)
!            i1 = 1
!            i2 = 1
!            do iz=1,lz2
!               call mxm  (ta2(i2),lx1,iym12,ly2,ta1(i1),ly1)
!               i1 = i1 + n1
!               i2 = i2 + n2
!            enddo
!            call mxm  (ta1,nxy1,izm12,lz2,dtx(1,e),lz1)
!
!            call col3 (ta1,wx,sm2(1,e),nxyz2)
!            call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)
!            i1 = 1
!            i2 = 1
!            do iz=1,lz2
!               call mxm  (ta2(i2),lx1,dym12,ly2,ta1(i1),ly1)
!               i1 = i1 + n1
!               i2 = i2 + n2
!            enddo
!            call mxm  (ta1,nxy1,izm12,lz2,ta2,lz1)
!            call add2 (dtx(1,e),ta2,nxyz1)
!
!            call col3 (ta1,wx,tm2(1,e),nxyz2)
!            call mxm  (ixtm12,lx1,ta1,lx2,ta2,nyz2)
!            i1 = 1
!            i2 = 1
!            do iz=1,lz2
!               call mxm  (ta2(i2),lx1,iym12,ly2,ta1(i1),ly1)
!               i1 = i1 + n1
!               i2 = i2 + n2
!            enddo
!            call mxm  (ta1,nxy1,dzm12,lz2,ta2,lz1)
!            call add2 (dtx(1,e),ta2,nxyz1)

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
      subroutine multdM1 (dx,x,rm1,sm1,tm1,isd,iflg)
C---------------------------------------------------------------------
C
C     Compute D*X
C     X    : input variable, defined on M1
C     DX   : output variable, defined on M2
C     Integration done on the M1 mesh        
C     RM1 : RXM1, RYM1 or RZM1
C     SM1 : SXM1, SYM1 or SZM1
C     TM1 : TXM1, TYM1 or TZM1
C     ISD : spatial direction (x=1,y=2,z=3)
C     IFLG: OPGRAD (iflg=0) or OPDIV (iflg=1)
C
C---------------------------------------------------------------------
      
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

      integer e
      integer nyz1,nxy2,nxyz1,nxyz2,n1,n2

C
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

             if (nio.eq.0) write(6,*) 
     $         'MULTDM1 not implemented for ifsplit.'
             call exitt

!             call mxm (dxm12,lx2,x(1,e),lx1,ta1,nyz1)
!             i1=1
!             i2=1
!             do iz=1,lz1
!               call mxm (ta1(i1),lx2,iytm12,ly1,ta2(i2),ly2)
!               i1=i1+n1
!               i2=i2+n2
!             enddo
!             call mxm  (ta2,nxy2,iztm12,lz1,dx(1,e),lz2)
!             call col2 (dx(1,e),rm2(1,e),nxyz2)
!
!             call mxm  (ixm12,lx2,x(1,e),lx1,ta3,nyz1) ! reuse ta3 below
!             i1=1
!             i2=1
!             do iz=1,lz1
!               call mxm (ta3(i1),lx2,dytm12,ly1,ta2(i2),ly2)
!               i1=i1+n1
!               i2=i2+n2
!             enddo
!             call mxm     (ta2,nxy2,iztm12,lz1,ta1,lz2)
!             call addcol3 (dx(1,e),ta1,sm2(1,e),nxyz2)
!
!c            call mxm (ixm12,lx2,x(1,e),lx1,ta1,nyz1) ! reuse ta3 from above
!             i1=1
!             i2=1
!             do iz=1,lz1
!               call mxm (ta3(i1),lx2,iytm12,ly1,ta2(i2),ly2)
!               i1=i1+n1
!               i2=i2+n2
!             enddo
!             call mxm (ta2,nxy2,dztm12,lz1,ta3,lz2)
!             call addcol3 (dx(1,e),ta3,tm2(1,e),nxyz2)
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
         call mxm(ixtm21,lx2,wk1,lx1,ta1,lx1)
         call mxm(ta1,lx2,iym21,lx1,dx(1,e),lx2)


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
