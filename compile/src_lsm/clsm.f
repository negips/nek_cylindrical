!======================================================================
!     Conservative Level Set Method
!     Author: Prabal S. Negi
!
!====================================================================== 
!-----------------------------------------------------------------------
       subroutine lsm_gennormals(phi)

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'MASS'
      include 'SOLN'

      include 'LSM'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      integer n
      integer i,j

      real phi(lv)

!     This holds the gradients      
      real ta1(lv)
      real ta2(lv)
      real ta3(lv)
      real ta4(lv)

      common /scruz/ ta1,ta2,ta3,ta4

!     normalized phi.  0 <= p <= 1.0
      real p(lv)
      common /newscrch1/ p

      real p1,p2
      real fac,norm
      real glmax,glmin

      real viscL,viscH
      character nam*4

      logical ifregg          ! regularize gradients

      n = lx1*ly1*lz1*nelv

      ifregg = .true.
      nam    = 'NORM'

!     This holds the (non-normalized) gradients of \phi
      call lsm_grad_nlmap(ta1,ta2,ta3,phi)

!     Regularize the gradients      
      if (ifregg) then
        viscL = 1.0e-1
        viscH = 1.0e-1
        call lsm_regularize_field(ta1,viscL,viscH,nam)
        call lsm_regularize_field(ta2,viscL,viscH,nam)
        if (ndim.eq.3) call lsm_regularize_field(ta3,viscL,viscH,nam)
      endif  

!     Almost the same thing as done in lsm_grad_nlmap
      p1 = glmin(phi,n)
      p2 = glmax(phi,n)

      call copy(p,phi,n)
      p1   = glmin(p,n)
      call cadd(p,-p1,n)              ! Make sure everything is positive
      p1   = glmax(p,n) 
      p2 = 1.0/p1
      call cmult(p,p2,n)            !  0<= p <= 1

      call rzero3(lsm_nx,lsm_ny,lsm_nz,n)
      do i=1,n
        fac    = ta1(i)**2 + ta2(i)**2
        if (ndim.eq.3) fac = fac + ta3(i)**2
        norm   = sqrt(fac)
!        if (abs(norm).lt.1.0e-14) norm = 1.0e-14
        lsm_nx(i) = ta1(i)
        lsm_ny(i) = ta2(i)
        if (ndim.eq.3) lsm_nz(i) = ta3(i)
        lsm_nnorm(i) = norm
      enddo  

      return
      end subroutine lsm_gennormals

!-----------------------------------------------------------------------

      subroutine lsm_genforc(phi)

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'MASS'
      include 'SOLN'

      include 'LSM'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      integer n,n2
      integer i,j

      real phi(lv)

!     This holds the gradients      
      real ta1(lv)
      real ta2(lv)
      real ta3(lv)
      real ta4(lv)

      common /scruz/ ta1,ta2,ta3,ta4

!     Holds the normal vectors 
      real tb1(lv)
      real tb2(lv)
      real tb3(lv)
      real tb4(lv)
      common /newscrch2/ tb1,tb2,tb3,tb4

!     normalized phi.  0 <= p <= 1.0
      real p(lv)
      common /newscrch1/ p

      real p1,p2
      real fac,norm,dotp
      real glmax,glmin

      real xst,xbl
      real blend
      real fbl    ! blended field

      integer intype

      n = lx1*ly1*lz1*nelv

!     This holds the (non-normalized) gradients of \phi
      call lsm_grad_nlmap(ta1,ta2,ta3,phi)

!     Almost the same thing as done in lsm_grad_nlmap
      p1 = glmin(phi,n)
      p2 = glmax(phi,n)

      call copy(p,phi,n)
      p1   = glmin(p,n)
      call cadd(p,-p1,n)              ! Make sure everything is positive
      p1   = glmax(p,n) 
      p2 = 1.0/p1
      call cmult(p,p2,n)            !  0<= p <= 1


      call rzero3(tb1,tb2,tb3,n)

      xst  = 1.0e-8
      xbl  = 1.0e-6
      do i=1,n
        fac    = ta1(i)**2 + ta2(i)**2
        if (ndim.eq.3) fac = fac + ta3(i)**2
        norm   = sqrt(fac)
        if (abs(norm).lt.1.0e-14) norm = 1.0e-14
        
        tb1(i) = ta1(i)/norm
        tb2(i) = ta2(i)/norm
        if (ndim.eq.3) tb3(i) = ta3(i)/norm
       
      enddo  
!     tbx is now the normals

!     Store the normal vectors      
      call rzero3(lsm_nx,lsm_ny,lsm_nz,n)
      call opcopy(lsm_nx,lsm_ny,lsm_nz,tb1,tb2,tb3)

      intype = 1
      call lsm_forc_gkreiss(tb1,tb2,tb3,ta4,tb4,p,intype)
      return
      end subroutine lsm_genforc
!---------------------------------------------------------------------- 

      subroutine lsm_grad_nlmap(px,py,pz,phi)

      implicit none

      include 'SIZE'
      include 'LSM'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      integer n
      integer i,j

!     Input      
      real phi(lv)

!     Outputs      
      real px(lv)
      real py(lv)
      real pz(lv)

!     Work arrays      
      real ta1(lv)
      real ta2(lv)
      real ta3(lv)
      real ta4(lv)
      common /newscrch2/ ta1,ta2,ta3,ta4

      real p(lv)
      common /newscrch1/ p

      real a,ai            ! alpha
      real fac,norm,p1,p2

      real glmin
      real glmax

      real scal,scali

      logical ifds            ! if dssum


      ifds   = .false.

      n = lx1*ly1*lz1*nelv

      call copy(p,phi,n)
      p1   = glmin(p,n)
      call cadd(p,-p1,n)              ! Make sure everything is positive
      p1   = glmax(p,n) + 1.0e-14
      scal = 1.0/p1
      call cmult(p,scal,n)            !  0<= p <= 1
      scali = 1.0/scal

      a  = lsm_alpha
      ai = 1.0/a
      do i=1,n
        ta4(i) = (p(i)**a)/(p(i)**a + (1.0-p(i))**a)
      enddo  
      
      call rzero3(ta1,ta2,ta3,n)
      call gradm1(ta1,ta2,ta3,ta4)

!     In principle we could also regularize after the transformation      
      if (ifds) then
        call dsavg(ta1)
        call dsavg(ta2)
        if (ndim.eq.3) call dsavg(ta3)
      endif  

      do i=1,n
        fac    = ai*((p(i)*(1.0-p(i)))**(1.0-a))*
     $                ((p(i)**a + (1.0 - p(i))**a)**2)
        px(i) = ta1(i)*fac*scali
        py(i) = ta2(i)*fac*scali
        pz(i) = ta3(i)*fac*scali
      enddo  

      return
      end subroutine lsm_grad_nlmap
!---------------------------------------------------------------------- 
      subroutine lsm_forc_gkreiss_v1(p,intype)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'SOLN'
      include 'LSM'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      integer n,n2
      integer i,j

      real p(lv)        ! phi

!     Work arrays 
      real wk1(lv)
      real wk2(lv)
      real wk3(lv)
      real wk4(lv)      
      real wk5(lv)
      real wk6(lv)
      real wk7(lv)
     
      common /scrns/ wk1,wk2,wk3,wk4,wk5,wk6,wk7

      real ta1(lv)
      real ta2(lv)
      real ta3(lv)
      real ta4(lv)
      common /scruz/ ta1,ta2,ta3,ta4


      real epsh
      real s,x,y

      integer intype
      real smalldivide

      n = lx1*ly1*lz1*nelv

!     We treat Diffusion implicitly

      call rzero3(wk1,wk2,wk3,n) 
      do i=1,n
         s        = p(i)*(1.0 - p(i))       ! S = p(1.0 - p)
         wk1(i)   = s*lsm_nx(i)
         wk2(i)   = s*lsm_ny(i)
         if (ndim.eq.3) wk3(i)   = s*lsm_nz(i)
      enddo

      n2 = lx2*ly2*lz2*nelv
!     Calculate divergence (its on Mesh 2)
!     Contains BM2      
      call opdiv(wk5,wk1,wk2,wk3)
      call col2(wk5,bm2inv,n2)       ! Remove BM2
      call mappr(wk4,wk5,wk6,wk7)    ! Map back to Mesh1 (wk4) 

!     Calculate convecting field: (wk1,wk2,wk3)/|n|
      do i=1,n
        y = lsm_nnorm(i)
        x = wk1(i)
        wk1(i) = smalldivide(x,y)
        x = wk2(i)
        wk2(i) = smalldivide(x,y)
        if (ndim.eq.3) then
          x = wk3(i)
          wk3(i) = smalldivide(x,y)
        endif  
      enddo

      call convect_new
     $      (wk5,lsm_nnorm,.false.,wk1,wk2,wk3,.false.)
      call invcol2(wk5,bm1,n)  ! local mass inverse

      do i=1,n
        y = lsm_nnorm(i)
        x = wk4(i)-wk5(i)
        lsm_forc(i) = smalldivide(x,y)
      enddo


      return
      end subroutine lsm_forc_gkreiss_v1
!---------------------------------------------------------------------- 

      subroutine lsm_forc_gkreiss(px,py,pz,wk1,wk2,p,intype)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'SOLN'
      include 'LSM'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      integer n,n2
      integer i,j

      real p(lv)        ! phi

!     This holds the gradients      
      real px(lv)
      real py(lv)
      real pz(lv)
      real wk1(lv)      


!     Holds the normal vectors 
      real nx(lv)
      real ny(lv)
      real nz(lv)
      real wk2(lv)

      real epsh
      real dotp,fac

      real xst,xbl,xin,fbl
      logical ifblend

      real blend

      integer intype

      ifblend = .true.

      n = lx1*ly1*lz1*nelv

      xst  = 1.0e-10
      xbl  = 1.0e-8

      epsh = eps_corr
      do i=1,n
!       \grad\phi . n      
        dotp = px(i)*lsm_nx(i) + py(i)*lsm_ny(i) 
        if (ndim.eq.3) dotp = dotp + pz(i)*lsm_nz(i)

        if (lsm_simp) then
          if (intype.eq.1) then
!           Explicit part (forcing)
            fac = -p(i)/2.0
          else
            fac = p(i)/2.0 - p(i)*tlag(i,1,1,1,1,2)
          endif  
        else  
          fac = epsh*dotp - p(i)*(1.0-p(i))
        endif  

!       (\grad\phi . n - \phi*(1.0 - \phi))*n      
        px(i) = fac*lsm_nx(i)
        py(i) = fac*lsm_ny(i)
        if (ndim.eq.3) pz(i) = fac*lsm_nz(i)
!       Still need to divide by the norm
      enddo

      n2 = lx2*ly2*lz2*nelv
!     Calculate divergence (its on Mesh 2)
!     Contains BM2      
      call opdiv(wk1,px,py,pz)
      call col2(wk1,bm2inv,n2)          ! Remove BM2
      call mappr(lsm_forc,wk1,nx,ny)    ! wk1 is on Mesh 2

      if (ifblend) then
!       testing
        xst  = 1.0e-10
        xbl  = 1.0e-8
        do i=1,n
          fac         = p(i)*(1.0 - p(i))      ! x for blending
          fbl         = blend(0.0,lsm_forc(i),fac,xst,xbl)
          lsm_forc(i) = fbl
        enddo
      endif  


      return
      end subroutine lsm_forc_gkreiss
!---------------------------------------------------------------------- 
      subroutine lsm_forc_shukla(p,intype)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'SOLN'
      include 'LSM'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      integer n
      integer i,j

      real p(lv)        ! phi

      real wk1(lv)
      real wk2(lv)
!      real wk3(lv)
!      real wk4(lv)      
!      real wk5(lv)
!      real wk6(lv)
!      real wk7(lv)
     
!      common /scrns/ wk1,wk2 !,wk3,wk4,wk5,wk6,wk7

      real smalldivide
      real sig,logf,logn,logd

      integer intype

      n = lx1*ly1*lz1*nelv
      
      if (intype.eq.1) call rzero(lsm_forc,n)
      if (intype.eq.-1) call rzero(lsm_impf,n)
      call rzero(wk1,n)
      call rzero(wk2,n)

      do i=1,n
!       Diffusion is implicit
!        if (intype.eq.1) wk1(i) = -p(i)*(1.0-p(i))

!       This is the semi implicit version      
        if (intype.eq.1) wk1(i) = -p(i)/2.0
        if (intype.eq.-1) wk1(i) = p(i)*(0.5 - lsm_phi(i))
      enddo

      call convect_new (wk2,wk1,.false.,lsm_nx,lsm_ny,lsm_nz,.false.)
      if (intype.eq.1) call invcol2 (wk2,bm1,n)  ! local mass inverse

      do i=1,n
!        lsm_forc(i)     = 0.01*smalldivide(wk2(i),lsm_nnorm(i))
!       Semi implicit version      
        if (intype.eq.1) lsm_forc(i) = smalldivide(wk2(i),lsm_nnorm(i))
        if (intype.eq.-1) lsm_impf(i) = smalldivide(wk2(i),lsm_nnorm(i))
      enddo  

      return
      end subroutine lsm_forc_shukla
!---------------------------------------------------------------------- 

      function simson_step(x)

      real simson_step
      real y

      if (x.le.0.0) then
        simson_step = 0.0
      elseif (x.ge.1.0) then
        simson_step = 1.0
      else
        y = 1.0/(x-1.0) + 1.0/x
        simson_step = 1.0/(1.0 + exp(y))
      endif  

      return
      end function
!---------------------------------------------------------------------- 

      function blend(u1,u2,x,xst,xbl)

      real blend  

      real u1,u2        ! two variables to blend
      real x            ! at what position
      real xst          ! start of blending
      real xbl          ! width of blending
      real y

      y = (x-xst)/xbl
      blend = u1 + (u2 - u1)*simson_step(y)

      return
      end function blend
!---------------------------------------------------------------------- 

      subroutine axhelm_spd(au,u,helm1,helm2,imsh,isd)

!     In case we want the diffusion to act differently on different
!     spectral components.
!     Axhelm in the core must also be modified for this to work        

      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'MASS'

      include 'LSM'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real au(lv),u(lv)
      real helm1(lv),helm2(lv)

      integer imsh,isd

      real lowpass(lv)
      real highpass(lv)
      real temp1(lv)
      real temp2(lv)
      real temp3(lv)

      real zro(lv)

      integer n

!     work arrays 
      real ta1(lv)
      real ta2(lv)
      real ta3(lv)
      real ta4(lv)
      common /newscrch2/ ta1,ta2,ta3,ta4

!     normalized phi.  0 <= p <= 1.0
      real nphi(lv)
      common /newscrch1/ nphi

      logical ifspvd
      common /spvd/ ifspvd

      integer intype

      n = lx1*ly1*lz1*nelv
!     Set to false so the normal axhelm works      
      ifspvd = .false.
 
      if (ifield.eq.2) then
        call copy(lowpass,u,n)
        call hpf_dist(highpass,u)   ! 
        call dsavg(highpass)
        call sub2(lowpass,highpass,n)

        call rzero(zro,n)
   
!       Just get BM1*\beta*u/dt
        call axhelm(au,u,zro,helm2,imsh,isd)
        
!       High Diffusion of high wavenumber components
        call axhelm(temp3,highpass,helm1,zro,imesh,isd)
        call cmult(temp3,1000.0,n)
        call add2(au,temp3,n)

!       Regular Diffusion for the lowpass component 
        call axhelm(temp3,lowpass,helm1,zro,imesh,isd)
        call add2(au,temp3,n)
      elseif (ifield.eq.3.and.lsm_simp) then
!        intype = 1
!        call lsm_grad_nlmap(ta1,ta2,ta3,u)
!        call normalize_phi(nphi,u)
!        call lsm_forc_gkreiss(ta1,ta2,ta3,ta4,    ! forcing
!     $            temp1,nphi,intype)
!
!!       Add non-linear term to operator
!        call axhelm(au,u,helm1,helm2,imesh,isd)
!        call col2(lsm_forc,bm1,n)
!        call add2(au,lsm_forc,n)

        call axhelm_lsm(au,u,helm1,helm2,imesh,isd)

      else
!       Just regular axhelm
        call axhelm(au,u,helm1,helm2,imesh,isd)
      endif   

      ifspvd = .true.

      return
      end subroutine axhelm_spd
!---------------------------------------------------------------------- 

      subroutine lsm_regularize_field(fld,viscL,viscH,name)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'SOLN'    ! TMULT
      include 'TSTEP'   ! ifield


      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real fld(lv)
      character name*4

      real resv1(lv),resv2(lv),resv3(lv)
      real dv1(lv),dv2(lv),dv3(lv),dv4(lv)
      common /scrns/  resv1
     $ ,              resv2
     $ ,              resv3
     $ ,              dv1  
     $ ,              dv2  
     $ ,              dv3
     $ ,              dv4 

      real h1(lv),h2(lv)  
      common /scrvh/  h1      ! viscosity
     $ ,              h2      ! rho/dt

      real viscL              ! viscosity for low wavenumbers
      real viscH              ! viscosity for High wavenumbers

      integer n
      integer ifld,imsh,isd,maxit
      real tli                ! input tolerance

      logical ifspvd

!      name = 'REGU'           ! Regularization

      imsh = 1
      maxit = 6000
      isd   = 1
      tli   = 1.0e-08


      n = lx1*ly1*lz1*nelv

      call col2(fld,bm1,n)
      
      call rone(h2,n)

      call rone(resv1,n)            ! mask
      call copy(resv2,tmult,n)      ! multiplicity

      ifspvd = .false.
      if (abs(viscL-viscH).gt.1.0e-12) ifspvd = .true.

      if (ifspvd) then
!        ifspvd = .false.
        call hpf_dist(dv1,fld)   ! High pass 
        call sub2(fld,dv1,n)     ! Low pass    
        call copy(dv2,fld,n)

!       Low wavenumber solution        
        call cfill(h1,viscL,n)
        call hmholtz(name,fld,dv2,h1,h2,resv1,resv2,imsh,tli,maxit,isd)
            
!       High wavenumber solution        
        call cfill(h1,viscH,n)
        call hmholtz(name,dv2,dv1,h1,h2,resv1,resv2,imsh,tli,maxit,isd)
        call add2(fld,dv2,n)
!        ifspvd = .true.
      else
!       Same viscosity for all wavenumbers        
        call cfill(h1,viscL,n)
        call copy(dv1,fld,n)
        call hmholtz(name,fld,dv1,h1,h2,resv1,resv2,imsh,tli,maxit,isd)
      endif

      return
      end subroutine lsm_regularize_field
!---------------------------------------------------------------------- 

      subroutine axhelm_lsm(au,u,helm1,helm2,imesh,isd)
C------------------------------------------------------------------
C
C     Compute the (Helmholtz) matrix-vector product,
C     AU = helm1*[A]u + helm2*[B]u, for NEL elements.
C
C------------------------------------------------------------------

      implicit none

      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
      include 'PARALLEL'
      include 'CTIMER'

      include 'LSM'
C
      REAL WDDX,WDDYT,WDDZT
      COMMON /FASTAX/ WDDX(LX1,LX1),WDDYT(LY1,LY1),WDDZT(LZ1,LZ1)
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
C
      REAL           AU    (LX1,LY1,LZ1,1)
     $ ,             U     (LX1,LY1,LZ1,1)
     $ ,             HELM1 (LX1,LY1,LZ1,1)
     $ ,             HELM2 (LX1,LY1,LZ1,1)

      REAL DUDR,DUDS,DUDT
      REAL TMP1,TMP2,TMP3
      COMMON /CTMP1/ DUDR  (LX1,LY1,LZ1)
     $ ,             DUDS  (LX1,LY1,LZ1)
     $ ,             DUDT  (LX1,LY1,LZ1)
     $ ,             TMP1  (LX1,LY1,LZ1)
     $ ,             TMP2  (LX1,LY1,LZ1)
     $ ,             TMP3  (LX1,LY1,LZ1)

      REAL           TM1   (LX1,LY1,LZ1)
      REAL           TM2   (LX1,LY1,LZ1)
      REAL           TM3   (LX1,LY1,LZ1)
      REAL           DUAX  (LX1)
      REAL           YSM1  (LX1)
      EQUIVALENCE    (DUDR,TM1),(DUDS,TM2),(DUDT,TM3)

      integer i,j,iz,e
      integer nxy,nyz,nxz,nxyz,ntot,nel
      integer imesh,isd
      real h1
      real term1,term2

      real ndot
      real smalldivide
      real t1,t2
      integer intype

      naxhm = naxhm + 1
      etime1 = dnekclock()

      nel=nelt
      if (imesh.eq.1) nel=nelv

      NXY=lx1*ly1
      NYZ=ly1*lz1
      NXZ=lx1*lz1
      NXYZ=lx1*ly1*lz1
      NTOT=NXYZ*NEL

      IF (.NOT.IFSOLV) CALL SETFAST(HELM1,HELM2,IMESH)
      CALL RZERO (AU,NTOT)

      do 100 e=1,nel
C
        if (ifaxis) call setaxdy ( ifrzer(e) )
C
        IF (ldim.EQ.2) THEN
C
C       2-d case ...............
C
           if (iffast(e)) then
C
C          Fast 2-d mode: constant properties and undeformed element
C
           h1 = helm1(1,1,1,e)
           call mxm   (wddx,lx1,u(1,1,1,e),lx1,tm1,nyz)
           call mxm   (u(1,1,1,e),lx1,wddyt,ly1,tm2,ly1)
           call col2  (tm1,g4m1(1,1,1,e),nxyz)
           call col2  (tm2,g5m1(1,1,1,e),nxyz)
           call add3  (au(1,1,1,e),tm1,tm2,nxyz)
           call cmult (au(1,1,1,e),h1,nxyz)
C
           else
C
C
           call mxm  (dxm1,lx1,u(1,1,1,e),lx1,dudr,nyz)
           call mxm  (u(1,1,1,e),lx1,dytm1,ly1,duds,ly1)
           call col3 (tmp1,dudr,g1m1(1,1,1,e),nxyz)
           call col3 (tmp2,duds,g2m1(1,1,1,e),nxyz)
           if (ifdfrm(e)) then
              call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
              call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
           endif
!          (\grad\phi.n)n
           j = (e-1)*nxyz
           do i=1,nxyz
             ndot = tmp1(i,1,1)*lsm_nx(j+i) + tmp2(i,1,1)*lsm_ny(j+i)
             t1 = ndot*lsm_nx(j+i) 
             t2 = ndot*lsm_ny(j+i)
             tmp1(i,1,1) = smalldivide(t1,lsm_nnorm(i))
             tmp2(i,1,1) = smalldivide(t2,lsm_nnorm(i))
           enddo  
           call col2 (tmp1,helm1(1,1,1,e),nxyz)
           call col2 (tmp2,helm1(1,1,1,e),nxyz)


!          Divergence part 
           call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
           call mxm  (tmp2,lx1,dym1,ly1,tm2,ly1)
           call add2 (au(1,1,1,e),tm1,nxyz)
           call add2 (au(1,1,1,e),tm2,nxyz)

        endif
C
        else
C
C       3-d case ...............
C
           if (iffast(e)) then
C
C          Fast 3-d mode: constant properties and undeformed element
C
           h1 = helm1(1,1,1,e)
           call mxm   (wddx,lx1,u(1,1,1,e),lx1,tm1,nyz)
           do 5 iz=1,lz1
           call mxm   (u(1,1,iz,e),lx1,wddyt,ly1,tm2(1,1,iz),ly1)
 5         continue
           call mxm   (u(1,1,1,e),nxy,wddzt,lz1,tm3,lz1)
           call col2  (tm1,g4m1(1,1,1,e),nxyz)
           call col2  (tm2,g5m1(1,1,1,e),nxyz)
           call col2  (tm3,g6m1(1,1,1,e),nxyz)
           call add3  (au(1,1,1,e),tm1,tm2,nxyz)
           call add2  (au(1,1,1,e),tm3,nxyz)
           call cmult (au(1,1,1,e),h1,nxyz)
C
           else
C
C
           call mxm(dxm1,lx1,u(1,1,1,e),lx1,dudr,nyz)
           do 10 iz=1,lz1
              call mxm(u(1,1,iz,e),lx1,dytm1,ly1,duds(1,1,iz),ly1)
   10      continue
           call mxm     (u(1,1,1,e),nxy,dztm1,lz1,dudt,lz1)
           call col3    (tmp1,dudr,g1m1(1,1,1,e),nxyz)
           call col3    (tmp2,duds,g2m1(1,1,1,e),nxyz)
           call col3    (tmp3,dudt,g3m1(1,1,1,e),nxyz)
           if (ifdfrm(e)) then
              call addcol3 (tmp1,duds,g4m1(1,1,1,e),nxyz)
              call addcol3 (tmp1,dudt,g5m1(1,1,1,e),nxyz)
              call addcol3 (tmp2,dudr,g4m1(1,1,1,e),nxyz)
              call addcol3 (tmp2,dudt,g6m1(1,1,1,e),nxyz)
              call addcol3 (tmp3,dudr,g5m1(1,1,1,e),nxyz)
              call addcol3 (tmp3,duds,g6m1(1,1,1,e),nxyz)
           endif
           call col2 (tmp1,helm1(1,1,1,e),nxyz)
           call col2 (tmp2,helm1(1,1,1,e),nxyz)
           call col2 (tmp3,helm1(1,1,1,e),nxyz)
           call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
           do 20 iz=1,lz1
              call mxm(tmp2(1,1,iz),lx1,dym1,ly1,tm2(1,1,iz),ly1)
   20      continue
           call mxm  (tmp3,nxy,dzm1,lz1,tm3,lz1)
           call add2 (au(1,1,1,e),tm1,nxyz)
           call add2 (au(1,1,1,e),tm2,nxyz)
           call add2 (au(1,1,1,e),tm3,nxyz)
C
           endif
c
        endif
C
 100  continue

!       semi-implicit non-linear term
        intype = -1
        call lsm_forc_shukla(u,intype)
        call add2(au,lsm_impf,ntot)
C
      if (ifh2) call addcol4 (au,helm2,bm1,u,ntot)
C
C     If axisymmetric, add a diagonal term in the radial direction (ISD=2)
C
      if (ifaxis.and.(isd.eq.2)) then
         do 200 e=1,nel
C
            if (ifrzer(e)) then
               call mxm(u  (1,1,1,e),lx1,datm1,ly1,duax,1)
               call mxm(ym1(1,1,1,e),lx1,datm1,ly1,ysm1,1)
            endif
c
            do 190 j=1,ly1
            do 190 i=1,lx1
C               if (ym1(i,j,1,e).ne.0.) then
                  if (ifrzer(e)) then
                     term1 = 0.0
                     if(j.ne.1) 
     $             term1 = bm1(i,j,1,e)*u(i,j,1,e)/ym1(i,j,1,e)**2
                     term2 =  wxm1(i)*wam1(1)*dam1(1,j)*duax(i)
     $                       *jacm1(i,1,1,e)/ysm1(i)
                  else
                   term1 = bm1(i,j,1,e)*u(i,j,1,e)/ym1(i,j,1,e)**2
                     term2 = 0.
                  endif
                  au(i,j,1,e) = au(i,j,1,e)
     $                          + helm1(i,j,1,e)*(term1+term2)
C               endif
  190       continue
  200    continue
      endif

      taxhm=taxhm+(dnekclock()-etime1)
      return
      end subroutine axhelm_lsm
!---------------------------------------------------------------------- 

      function smalldivide(x,y)

!     Evaluate x/y for vanishing x,y

!     x/y = sign(x)/sign(y)*exp( log(|x|+\eps) - log(|y|+\eps) )

      real smalldivide

      real x,y
      real sigx,sigy,lognum,logden,logd
      real eps
      
      eps         = 1.0e-20

      sigx        = sign(1.0,x)
      sigy        = sign(1.0,y)
      lognum      = log(abs(x)+eps)
      logden      = log(abs(y)+eps)
      logd        = lognum - logden
      smalldivide = (sigx/sigy)*exp(logd)

      return
      end function smalldivide
!---------------------------------------------------------------------- 



