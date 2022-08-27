!======================================================================
!     Routines for Weak Laplacian/Weak Stress divergence 
!     in Cylindrical Coordinates.
!     Assuming Fourier in the 3rd (\theta) direction      
!     Author: Prabal S. Negi
!
!====================================================================== 
!-----------------------------------------------------------------------

      subroutine axhmsf_cyl(Au1r,Au2r,Au3r,Au1i,Au2i,Au3i,
     $                      u1r,u2r,u3r,u1i,u2i,u3i,h1,h2)

!     Fluid (MATMOD .GE. 0) :  Hij Uj = Aij*Uj + H2*B*Ui 

      implicit none

      include 'SIZE'
      include 'TSTEP'   ! nelfld
      include 'MASS'

      real u1r(1),u2r(1),u3r(1),u1i(1),u2i(1),u3i(1)
      real Au1r(1),Au2r(1),Au3r(1),Au1i(1),Au2i(1),Au3i(1)
      real h1(1),h2(1)

      integer matmod,nel,ntot1

      logical ifdfrm,iffast,ifh2,ifsolv
      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv


      matmod = 0                ! Newtonian Fluids
      nel    = nelfld(ifield)   ! ifield should be already set

      ntot1 = lx1*ly1*lz1*nel


!     Common blocks are used in succession.

!     Aij*Uj      
      call stnrate_cyl(u1r,u2r,u3r,u1i,u2i,u3i,nel,matmod)
      call stress_cyl (h1,h2,nel,matmod)
      call div_stress_cyl(Au1r,Au2r,Au3r,Au1i,Au2i,Au3i,nel)   ! aijuj

!     Add Helmholtz contributions
!     + H2*B*Ui      
      call addcol4 (Au1r,bm1,h2,u1r,ntot1)
      call addcol4 (Au2r,bm1,h2,u2r,ntot1)
      call addcol4 (Au3r,bm1,h2,u3r,ntot1)

      call addcol4 (Au1i,bm1,h2,u1i,ntot1)
      call addcol4 (Au2i,bm1,h2,u2i,ntot1)
      call addcol4 (Au3i,bm1,h2,u3i,ntot1)

      return
      end subroutine axhmsf_cyl             
!-----------------------------------------------------------------------
     
      subroutine stnrate_cyl(u1r,u2r,u3r,u1i,u2i,u3i,nel,matmod)

      implicit none

!     Compute strainrates

!     CAUTION : Stresses and strainrates share the same scratch commons

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'          ! ym1 == R
      include 'TSTEP'

      include 'F3D'

!!     Real Variables      
!      real erxt,errt,erxx,erxr,errr,ertt
!      common /scrpsn3/ erxt(lx1*ly1*lz1*lelt)      ! Er_x\theta
!     $               , errt(lx1*ly1*lz1*lelt)      ! Er_Rt
!     $               , erxx(lx1*ly1*lz1*lelt)      ! Er_xx
!     $               , erxr(lx1*ly1*lz1*lelt)      ! Er_xR
!     $               , errr(lx1*ly1*lz1*lelt)      ! Er_RR
!     $               , ertt(lx1*ly1*lz1*lelt)      ! Er_\theta\theta
!
!!     Made a new scratch array here.      
!!     Imaginary Variables
!      real eixt,eirt,eixx,eixr,eirr,eitt,wrk
!      common /scrpsn4/  eixt(lx1*ly1*lz1*lelt)      ! Ei_x\theta
!     $                , eirt(lx1*ly1*lz1*lelt)      ! Ei_Rt
!     $                , eixx(lx1*ly1*lz1*lelt)      ! Ei_xx
!     $                , eixr(lx1*ly1*lz1*lelt)      ! Ei_xR
!     $                , eirr(lx1*ly1*lz1*lelt)      ! Ei_RR
!     $                , eitt(lx1*ly1*lz1*lelt)      ! Ei_\theta\theta
!     $                , wrk(lx1*ly1*lz1*lelt)       ! work
!
!
!!     Gradients of real part of velocities      
!      real ur1x,ur1r,ur1t,ur2x,ur2r,ur2t,ur3x,ur3r,ur3t
!      common /scrpsn1/ ur1x(lx1*ly1*lz1*lelt)      ! du1dx_r
!     $               , ur1r(lx1*ly1*lz1*lelt)      ! du1dr_r
!     $               , ur2x(lx1*ly1*lz1*lelt)      ! du2dx_r
!     $               , ur2r(lx1*ly1*lz1*lelt)      ! du2dr_r
!     $               , ur3x(lx1*ly1*lz1*lelt)      ! du3dx_r
!     $               , ur3r(lx1*ly1*lz1*lelt)      ! du3dr_r
!
!
!!     Gradient of imaginary part of velocities      
!      real ui1x,ui1r,ui1t,ui2x,ui2r,ui2t,ui3x,ui3r,ui3t
!      common /scrpsn2/ ui1x(lx1*ly1*lz1*lelt)      ! du1dx_i
!     $               , ui1r(lx1*ly1*lz1*lelt)      ! du1dr_i
!     $               , ui2x(lx1*ly1*lz1*lelt)      ! du2dx_i
!     $               , ui2r(lx1*ly1*lz1*lelt)      ! du2dr_i
!     $               , ui3x(lx1*ly1*lz1*lelt)      ! du3dx_i
!     $               , ui3r(lx1*ly1*lz1*lelt)      ! du3dr_i


      real u1r,u2r,u3r
      dimension u1r(lx1,ly1,lz1,lelv)
     $        , u2r(lx1,ly1,lz1,lelv)
     $        , u3r(lx1,ly1,lz1,lelv)

      real u1i,u2i,u3i
      dimension u1i(lx1,ly1,lz1,lelv)
     $        , u2i(lx1,ly1,lz1,lelv)
     $        , u3i(lx1,ly1,lz1,lelv)

!     Work Arrays. Discarded after this subroutine
      real rinv(lx1*ly1*lz1*lelv),wk1(lx1*ly1*lz1*lelv)
      real wk2(lx1*ly1*lz1*lelv),hii(lx1*ly1*lz1*lelv)
      common /scrpsn5/ rinv,wk1,wk2,hii        ! UXYZ already uses scrsf

      integer nel,ntot1,matmod


      ntot1 = lx1*ly1*lz1*nel

!     rzero3 zeros all 3 components
!     without checking dimensionality

!     Zero real parts
      call rzero3(ur1x,ur2x,ur3x,ntot1) 
      call rzero3(ur1r,ur2r,ur3r,ntot1) 

!     Zero Imaginary parts      
      call rzero3(ui1x,ui2x,ui3x,ntot1) 
      call rzero3(ui1r,ui2r,ui3r,ntot1) 

!     Zero Real parts  
      call rzero3 (erxx,errr,ertt,ntot1)
      call rzero3 (erxr,erxt,errt,ntot1)

!     Zero imaginary parts  
      call rzero3 (eixx,eirr,eitt,ntot1)
      call rzero3 (eixr,eixt,eirt,ntot1)

!     uxyz does not zero out variables.
!     values are just added on
!     Derivatives of Real parts      
      call uxyz  (u1r,ur1x,ur1r,wk1,nel)
      call uxyz  (u2r,ur2x,ur2r,wk1,nel)
      call uxyz  (u3r,ur3x,ur3r,wk1,nel)

!     Since the division by the Jacobian is missing
      call invcol2 (ur1x,jacm1,ntot1)
      call invcol2 (ur1r,jacm1,ntot1)
      call invcol2 (ur2x,jacm1,ntot1)
      call invcol2 (ur2r,jacm1,ntot1)
      call invcol2 (ur3x,jacm1,ntot1)
      call invcol2 (ur3r,jacm1,ntot1)

!     Derivatives of Imaginary parts
      call uxyz  (u1i,ui1x,ui1r,wk1,nel)
      call uxyz  (u2i,ui2x,ui2r,wk1,nel)
      call uxyz  (u3i,ui3x,ui3r,wk1,nel)

!     Since the division by the Jacobian is missing
      call invcol2 (ui1x,jacm1,ntot1)
      call invcol2 (ui1r,jacm1,ntot1)
      call invcol2 (ui2x,jacm1,ntot1)
      call invcol2 (ui2r,jacm1,ntot1)
      call invcol2 (ui3x,jacm1,ntot1)
      call invcol2 (ui3r,jacm1,ntot1)

      
!     rinv = 1/R
      if (ifcyl_f3d) then
        call invers2(rinv,ym1,ntot1)
      else
        call rone(rinv,ntot1)  
      endif  

!!    Er_xx 
!     (Real)
      call copy(erxx,ur1x,ntot1)    ! du1/dx

!!    Ei_xx 
!     (Imaginary)
      call copy(eixx,ui1x,ntot1)    ! du1(im)/dx

!!    Er_Rx
!     (Real)      
      call copy(erxr,ur1r,ntot1)    ! du1/dr
      call add2(erxr,ur2x,ntot1)    ! + du2/dx
      call cmult(erxr,0.5,ntot1)    ! 1/2*[du1/dr + du2/dx]

!!    Ei_Rx
!     (Imaginary)      
      call copy(eixr,ui1r,ntot1)    ! du1(im)/dr
      call add2(eixr,ui2x,ntot1)    ! + du2(im)/dx
      call cmult(eixr,0.5,ntot1)    ! 1/2*[du1(im)/dr + du2(im)/dx]


!!    Er_\thetax
!     (Real)      
      call col3(erxt,u1i,rinv,ntot1)       ! u1(im)/R
      call cmult(erxt,-k_f3d,ntot1)       ! -k/R*u1(im)
      call add2(erxt,ur3x,ntot1)           ! + du3/dx
      call cmult(erxt,0.5,ntot1)           ! 1/2*[du3/dx - k/R*u1(im)]

!!    Ei_\thetax
!     (Imaginary) 
      call col3(eixt,u1r,rinv,ntot1)       ! u1/R
      call cmult(eixt,k_f3d,ntot1)         ! k*u1/R
      call add2(eixt,ui3x,ntot1)           ! + du3(im)/dx
      call cmult(eixt,0.5,ntot1)           ! 1/2*[du3(im)/dx + k/R*u1]


!!    Er_RR
!     (Real)      
      call copy(errr,ur2r,ntot1)           ! du2/dr

!!    Ei_RR
!     (Imaginary)      
      call copy(eirr,ui2r,ntot1)           ! du2(im)/dr

      
!!    Er_\thetaR
!     (Real)      
      call copy(errt,ur3r,ntot1)           ! du3/dr
      if (ifcyl_f3d) then
        call col3(wk1,u3r,rinv,ntot1)      ! u3/R
        call sub2(errt,wk1,ntot1)          ! du3/dr - u3/R
      endif        
      call col3(wk2,u2i,rinv,ntot1)        ! u2(im)/R
      call add2s2(errt,wk2,-k_f3d,ntot1)   ! du3/dr - u3/R - k/R*u2(im)
      call cmult(errt,0.5,ntot1)           ! 0.5*[du3/dr - u3/R - k/R*u2(im)]

!!    Ei_\thetaR
!     (Imaginary)      
      call copy(eirt,ui3r,ntot1)           ! du3(im)/dr
      if (ifcyl_f3d) then
        call col3(wk1,u3i,rinv,ntot1)      ! u3(im)/R
        call sub2(eirt,wk1,ntot1)          ! du3(im)/dr - u3(im)/R
      endif        
      call col3(wk2,u2r,rinv,ntot1)        ! u2/R
      call add2s2(eirt,wk2,k_f3d,ntot1)    ! du3(im)/dr - u3(im)/R + k*u2/R
      call cmult(eirt,0.5,ntot1)           ! 0.5*[du3(im)/dr - u3(im)/R + k*u2/R]


!!    Er_\theta\theta
!     (Real)      
      call col3(ertt,u3i,rinv,ntot1)       ! u3(im)/R      
      call cmult(ertt,-k_f3d,ntot1)        ! -k*u3(im)/R
      if (ifcyl_f3d) then
        call xaddcol3(ertt,u2r,rinv,ntot1) ! -k*u3(im)/R + u2/R
      endif  

!!    Ei_\theta\theta
!     (Imaginary)      
      call col3(eitt,u3r,rinv,ntot1)       ! u3/R      
      call cmult(eitt,k_f3d,ntot1)         ! k*u3/R
      if (ifcyl_f3d) then
        call xaddcol3(eitt,u2i,rinv,ntot1) ! k*u3/R + u2(im)/R
      endif


      return
      end subroutine stnrate_cyl
c-----------------------------------------------------------------------
      subroutine stress_cyl (h1,h2,nel,matmod)
C
C     MATMOD.GE.0        Fluid material models
C     MATMOD.LT.0        Solid material models
C
C     CAUTION : Stresses and strainrates share the same scratch commons

      implicit none

      include 'SIZE'
      include 'F3D'
 
      real t11,t22,t33,hii
      common /scrpsn5/ t11(lx1,ly1,lz1,lelt)
     $               , t22(lx1,ly1,lz1,lelt)
     $               , t33(lx1,ly1,lz1,lelt)
     $               , hii(lx1,ly1,lz1,lelt)


!!     Real Variables      
!      real erxt,errt,erxx,erxr,errr,ertt
!      common /scrpsn3/ erxt(lx1*ly1*lz1*lelt)      ! Er_x\theta
!     $             , errt(lx1*ly1*lz1*lelt)      ! Er_Rt
!     $             , erxx(lx1*ly1*lz1*lelt)      ! Er_xx
!     $             , erxr(lx1*ly1*lz1*lelt)      ! Er_xR
!     $             , errr(lx1*ly1*lz1*lelt)      ! Er_RR
!     $             , ertt(lx1*ly1*lz1*lelt)      ! Er_\theta\theta
!
!!     Made a new scratch array here.      
!!     Imaginary Variables
!      real eixt,eirt,eixx,eixr,eirr,eitt,wrk
!      common /scrpsn4/  eixt(lx1*ly1*lz1*lelt)      ! Ei_x\theta
!     $                , eirt(lx1*ly1*lz1*lelt)      ! Ei_Rt
!     $                , eixx(lx1*ly1*lz1*lelt)      ! Ei_xx
!     $                , eixr(lx1*ly1*lz1*lelt)      ! Ei_xR
!     $                , eirr(lx1*ly1*lz1*lelt)      ! Ei_RR
!     $                , eitt(lx1*ly1*lz1*lelt)      ! Ei_\theta\theta
!     $                , wrk(lx1*ly1*lz1*lelt)       ! work
!

      real h1(lx1,ly1,lz1,1),h2(lx1,ly1,lz1,1)
      integer nel,ntot1,matmod
      real const


      ntot1 = lx1*ly1*lz1*nel

      IF (MATMOD.EQ.0) THEN

c        newtonian fluids

         const = 2.0
         call cmult2 (hii,h1,const,ntot1)

!        Real Parts         
         call col2   (erxt,hii,ntot1)
         call col2   (errt,hii,ntot1)
         call col2   (erxx,hii,ntot1)
         call col2   (erxr,hii,ntot1)
         call col2   (errr,hii,ntot1)
         call col2   (ertt,hii,ntot1)

!        Imaginary parts            
         call col2   (eixt,hii,ntot1)
         call col2   (eirt,hii,ntot1)
         call col2   (eixx,hii,ntot1)
         call col2   (eixr,hii,ntot1)
         call col2   (eirr,hii,ntot1)
         call col2   (eitt,hii,ntot1)

      elseif (matmod.eq.-1) then

!        Elastic solids
!        Just removed it.
!        Somebody will need to implement it carefully        

         if (nid.eq.0)
     $     write(6,*) 'Elasticity not implemented for F3D'

         call exitt

      endif


      return
      end subroutine stress_cyl
!-----------------------------------------------------------------------
      subroutine div_stress_cyl (Au1r,Au2r,Au3r,Au1i,Au2i,Au3i,nel)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'GEOM'

      include 'F3D'


!!     Real Variables      
!      real erxt,errt,erxx,erxr,errr,ertt
!      common /scrpsn3/ erxt(lx1*ly1*lz1*lelt)      ! Er_x\theta
!     $               , errt(lx1*ly1*lz1*lelt)      ! Er_Rt
!     $               , erxx(lx1*ly1*lz1*lelt)      ! Er_xx
!     $               , erxr(lx1*ly1*lz1*lelt)      ! Er_xR
!     $               , errr(lx1*ly1*lz1*lelt)      ! Er_RR
!     $               , ertt(lx1*ly1*lz1*lelt)      ! Er_\theta\theta
!
!!     Made a new scratch array here.      
!!     Imaginary Variables
!      real eixt,eirt,eixx,eixr,eirr,eitt,wrk
!      common /scrpsn4/  eixt(lx1*ly1*lz1*lelt)      ! Ei_x\theta
!     $                , eirt(lx1*ly1*lz1*lelt)      ! Ei_Rt
!     $                , eixx(lx1*ly1*lz1*lelt)      ! Ei_xx
!     $                , eixr(lx1*ly1*lz1*lelt)      ! Ei_xR
!     $                , eirr(lx1*ly1*lz1*lelt)      ! Ei_RR
!     $                , eitt(lx1*ly1*lz1*lelt)      ! Ei_\theta\theta
!     $                , wrk(lx1*ly1*lz1*lelt)       ! work

      real rinv(lx1*ly1*lz1*lelv),wk1(lx1*ly1*lz1*lelv)
      real wk2(lx1*ly1*lz1*lelv),hii(lx1*ly1*lz1*lelv)
      common /scrpsn5/ rinv,wk1,wk2,hii


      real Au1r,Au2r,Au3r
      dimension Au1r(lx1,ly1,lz1,1)
     $        , Au2r(lx1,ly1,lz1,1)
     $        , Au3r(lx1,ly1,lz1,1)

      real Au1i,Au2i,Au3i
      dimension Au1i(lx1,ly1,lz1,1)
     $        , Au2i(lx1,ly1,lz1,1)
     $        , Au3i(lx1,ly1,lz1,1)
     
      integer nel,ntot1


      ntot1 = lx1*ly1*lz1*nel

!     rinv = 1/R
      if (ifcyl_f3d) then 
        call invers2(rinv,ym1,ntot1)
      else
        call rone(rinv,ntot1)
      endif  
      call col2c(rinv,bm1,-k_f3d,ntot1)   ! rinv = -k*BM1/R

!     Real Variables      
      call ttxyz(Au1r,erxx,erxr,erxt,nel)  ! [erxx*dv/dx + erxr*dv/dr]*BM1
      call col3(wk1,eixt,rinv,ntot1)       ! eixt*k*BM1/R    
      call sub2(Au1r,wk1,ntot1)

      call ttxyz(Au2r,erxr,errr,errt,nel)  ! [erxr*dv/dx + errr*dv/dr]*BM1
      call col3(wk1,eirt,rinv,ntot1)       ! eirt*k*BM1/R    
      call sub2(Au2r,wk1,ntot1)
!     Additional terms from vector gradient of test functions
      if (ifcyl_f3d) then
        call invcol2(wk1,ym1,ntot1)
        call col2c(wk1,bm1,-1.0,ntot1)       ! -BM1*v/R
        call Xaddcol3(Au2r,errt,wk1,ntot1)   ! -BM1*v/R*[errt]
      endif        

      call ttxyz(Au3r,erxt,errt,ertt,nel)  ! [erxt*dv/dx + errt*dv/dr]*BM1
      call col3(wk1,eitt,rinv,ntot1)       ! eitt*k*BM1/R    
      call sub2(Au3r,wk1,ntot1)
!     Additional terms from vector gradient of test functions
      if (ifcyl_f3d) then
        call invcol2(wk1,ym1,ntot1)
        call col2(wk1,bm1,ntot1)             ! BM1*v/R
        call Xaddcol3(Au3r,ertt,wk1,ntot1)   ! BM1*v/R*[ertt]
      endif        


!     Imaginary Variables      
      call ttxyz(Au1i,eixx,eixr,eixt,nel)  ! [eixx*dv/dx + eixr*dv/dr]*BM1
      call col3(wk1,erxt,rinv,ntot1)       ! erxt*k*BM1/R    
      call add2(Au1i,wk1,ntot1)

      call ttxyz(Au2i,eixr,eirr,eirt,nel) ! [erxr*dv/dx + errr*dv/dr]*BM1
      call col3(wk1,errt,rinv,ntot1)      ! errt*k*BM1/R    
      call add2(Au2i,wk1,ntot1)
!     Additional terms from vector gradient of test functions
      if (ifcyl_f3d) then
        call invcol2(wk1,ym1,ntot1)
        call col2c(wk1,bm1,-1.0,ntot1)       ! -BM1*v/R
        call Xaddcol3(Au2i,eirt,wk1,ntot1)   ! -BM1*v/R*[eirt]
      endif        

      call ttxyz(Au3i,eixt,eirt,eitt,nel) ! [eixt*dv/dx + eirt*dv/dr]*BM1
      call col3(wk1,ertt,rinv,ntot1)      ! ertt*k*BM1/R    
      call add2(Au3i,wk1,ntot1)
!     Additional terms from vector gradient of test functions
      if (ifcyl_f3d) then
        call invcol2(wk1,ym1,ntot1)
        call col2(wk1,bm1,ntot1)             ! BM1*v/R
        call Xaddcol3(Au3i,eitt,wk1,ntot1)   ! BM1*v/R*[eitt]
      endif        


!      if (ifaxis)    call axitzz (au2,tzz,nel)
!      if (ldim.eq.3) call ttxyz  (au3,txz,tyz,tzz,nel)

      return
      end subroutine div_stress_cyl
!-----------------------------------------------------------------------
!---------------------------------------------------------------------- 
      subroutine axhmsf_cyl_real(Au1,Au2,Au3,u1,u2,u3,h1,h2)

!     Fluid (MATMOD .GE. 0) :  Hij Uj = Aij*Uj + H2*B*Ui 

      implicit none

      include 'SIZE'
      include 'TSTEP'   ! nelfld
      include 'MASS'

      real u1(1),u2(1),u3(1)
      real Au1(1),Au2(1),Au3(1)
      real h1(1),h2(1)

      integer matmod,nel,ntot1

      logical ifdfrm,iffast,ifh2,ifsolv
      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv


      matmod = 0                ! Newtonian Fluids
      nel    = nelfld(ifield)   ! ifield should be already set

      ntot1 = lx1*ly1*lz1*nel


!     Common blocks are used in succession.

!     Aij*Uj      
      call stnrate_cyl_real(u1,u2,u3,nel,matmod)
      call stress_cyl (h1,h2,nel,matmod)
      call div_stress_cyl_real(Au1,Au2,Au3,nel)   ! aijuj

!     Add Helmholtz contributions
!     + H2*B*U 
      call addcol4 (Au1,bm1,h2,u1,ntot1)
      call addcol4 (Au2,bm1,h2,u2,ntot1)
      call addcol4 (Au3,bm1,h2,u3,ntot1)


      return
      end subroutine axhmsf_cyl_real 
!-----------------------------------------------------------------------
      subroutine stnrate_cyl_real(u1r,u2r,u3r,nel,matmod)

      implicit none

!     Compute strainrates

!     CAUTION : Stresses and strainrates share the same scratch commons

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'          ! ym1 == R
      include 'TSTEP'

      include 'F3D'

!!     Real Variables      
!      real erxt,errt,erxx,erxr,errr,ertt
!      common /scrpsn3/ erxt(lx1*ly1*lz1*lelt)      ! Er_x\theta
!     $               , errt(lx1*ly1*lz1*lelt)      ! Er_Rt
!     $               , erxx(lx1*ly1*lz1*lelt)      ! Er_xx
!     $               , erxr(lx1*ly1*lz1*lelt)      ! Er_xR
!     $               , errr(lx1*ly1*lz1*lelt)      ! Er_RR
!     $               , ertt(lx1*ly1*lz1*lelt)      ! Er_\theta\theta
!
!!     Made a new scratch array here.      
!!     Imaginary Variables
!      real eixt,eirt,eixx,eixr,eirr,eitt,wrk
!      common /scrpsn4/  eixt(lx1*ly1*lz1*lelt)      ! Ei_x\theta
!     $                , eirt(lx1*ly1*lz1*lelt)      ! Ei_Rt
!     $                , eixx(lx1*ly1*lz1*lelt)      ! Ei_xx
!     $                , eixr(lx1*ly1*lz1*lelt)      ! Ei_xR
!     $                , eirr(lx1*ly1*lz1*lelt)      ! Ei_RR
!     $                , eitt(lx1*ly1*lz1*lelt)      ! Ei_\theta\theta
!     $                , wrk(lx1*ly1*lz1*lelt)       ! work
!
!
!!     Gradients of real part of velocities      
!      real ur1x,ur1r,ur1t,ur2x,ur2r,ur2t,ur3x,ur3r,ur3t
!      common /scrpsn1/ ur1x(lx1*ly1*lz1*lelt)      ! du1dx_r
!     $               , ur1r(lx1*ly1*lz1*lelt)      ! du1dr_r
!     $               , ur2x(lx1*ly1*lz1*lelt)      ! du2dx_r
!     $               , ur2r(lx1*ly1*lz1*lelt)      ! du2dr_r
!     $               , ur3x(lx1*ly1*lz1*lelt)      ! du3dx_r
!     $               , ur3r(lx1*ly1*lz1*lelt)      ! du3dr_r
!
!
!!     Gradient of imaginary part of velocities      
!      real ui1x,ui1r,ui1t,ui2x,ui2r,ui2t,ui3x,ui3r,ui3t
!      common /scrpsn2/ ui1x(lx1*ly1*lz1*lelt)      ! du1dx_i
!     $               , ui1r(lx1*ly1*lz1*lelt)      ! du1dr_i
!     $               , ui2x(lx1*ly1*lz1*lelt)      ! du2dx_i
!     $               , ui2r(lx1*ly1*lz1*lelt)      ! du2dr_i
!     $               , ui3x(lx1*ly1*lz1*lelt)      ! du3dx_i
!     $               , ui3r(lx1*ly1*lz1*lelt)      ! du3dr_i


      real u1r,u2r,u3r
      dimension u1r(lx1,ly1,lz1,lelv)
     $        , u2r(lx1,ly1,lz1,lelv)
     $        , u3r(lx1,ly1,lz1,lelv)

      real u1i,u2i,u3i
      dimension u1i(lx1,ly1,lz1,lelv)
     $        , u2i(lx1,ly1,lz1,lelv)
     $        , u3i(lx1,ly1,lz1,lelv)

!     Work Arrays. Discarded after this subroutine
      real rinv(lx1*ly1*lz1*lelv),wk1(lx1*ly1*lz1*lelv)
      real wk2(lx1*ly1*lz1*lelv),hii(lx1*ly1*lz1*lelv)
      common /scrpsn5/ rinv,wk1,wk2,hii        ! UXYZ already uses scrsf

      integer nel,ntot1,matmod


      ntot1 = lx1*ly1*lz1*nel

!     Since we don't pass the imaginary component.
!     I create a new array and set it to zero.
!     Obviously not optimum for memory.
!     But the rest of the routine is then unchanged.      

      call rzero3(u1i,u2i,u3i,ntot1)

!     rzero3 zeros all 3 components
!     without checking dimensionality

!     Zero real parts
      call rzero3(ur1x,ur2x,ur3x,ntot1) 
      call rzero3(ur1r,ur2r,ur3r,ntot1) 

!     Zero Imaginary parts      
      call rzero3(ui1x,ui2x,ui3x,ntot1) 
      call rzero3(ui1r,ui2r,ui3r,ntot1) 

!     Zero Real parts  
      call rzero3 (erxx,errr,ertt,ntot1)
      call rzero3 (erxr,erxt,errt,ntot1)

!     Zero imaginary parts  
      call rzero3 (eixx,eirr,eitt,ntot1)
      call rzero3 (eixr,eixt,eirt,ntot1)

!     uxyz does not zero out variables.
!     values are just added on
!     Derivatives of Real parts      
      call uxyz  (u1r,ur1x,ur1r,wk1,nel)
      call uxyz  (u2r,ur2x,ur2r,wk1,nel)
      call uxyz  (u3r,ur3x,ur3r,wk1,nel)

!     Since the division by the Jacobian is missing
      call invcol2 (ur1x,jacm1,ntot1)
      call invcol2 (ur1r,jacm1,ntot1)
      call invcol2 (ur2x,jacm1,ntot1)
      call invcol2 (ur2r,jacm1,ntot1)
      call invcol2 (ur3x,jacm1,ntot1)
      call invcol2 (ur3r,jacm1,ntot1)

!     Derivatives of Imaginary parts
      call uxyz  (u1i,ui1x,ui1r,wk1,nel)
      call uxyz  (u2i,ui2x,ui2r,wk1,nel)
      call uxyz  (u3i,ui3x,ui3r,wk1,nel)

!     Since the division by the Jacobian is missing
      call invcol2 (ui1x,jacm1,ntot1)
      call invcol2 (ui1r,jacm1,ntot1)
      call invcol2 (ui2x,jacm1,ntot1)
      call invcol2 (ui2r,jacm1,ntot1)
      call invcol2 (ui3x,jacm1,ntot1)
      call invcol2 (ui3r,jacm1,ntot1)

      
!     rinv = 1/R
      if (ifcyl_f3d) then
        call invers2(rinv,ym1,ntot1)
      else
        call rone(rinv,ntot1)  
      endif  

!!    Er_xx 
!     (Real)
      call copy(erxx,ur1x,ntot1)    ! du1/dx

!!    Ei_xx 
!     (Imaginary)
      call copy(eixx,ui1x,ntot1)    ! du1(im)/dx

!!    Er_Rx
!     (Real)      
      call copy(erxr,ur1r,ntot1)    ! du1/dr
      call add2(erxr,ur2x,ntot1)    ! + du2/dx
      call cmult(erxr,0.5,ntot1)    ! 1/2*[du1/dr + du2/dx]

!!    Ei_Rx
!     (Imaginary)      
      call copy(eixr,ui1r,ntot1)    ! du1(im)/dr
      call add2(eixr,ui2x,ntot1)    ! + du2(im)/dx
      call cmult(eixr,0.5,ntot1)    ! 1/2*[du1(im)/dr + du2(im)/dx]


!!    Er_\thetax
!     (Real)      
      call col3(erxt,u1i,rinv,ntot1)       ! u1(im)/R
      call cmult(erxt,-k_f3d,ntot1)       ! -k/R*u1(im)
      call add2(erxt,ur3x,ntot1)           ! + du3/dx
      call cmult(erxt,0.5,ntot1)           ! 1/2*[du3/dx - k/R*u1(im)]

!!    Ei_\thetax
!     (Imaginary) 
      call col3(eixt,u1r,rinv,ntot1)       ! u1/R
      call cmult(eixt,k_f3d,ntot1)         ! k*u1/R
      call add2(eixt,ui3x,ntot1)           ! + du3(im)/dx
      call cmult(eixt,0.5,ntot1)           ! 1/2*[du3(im)/dx + k/R*u1]


!!    Er_RR
!     (Real)      
      call copy(errr,ur2r,ntot1)           ! du2/dr

!!    Ei_RR
!     (Imaginary)      
      call copy(eirr,ui2r,ntot1)           ! du2(im)/dr

      
!!    Er_\thetaR
!     (Real)      
      call copy(errt,ur3r,ntot1)           ! du3/dr
      if (ifcyl_f3d) then
        call col3(wk1,u3r,rinv,ntot1)      ! u3/R
        call sub2(errt,wk1,ntot1)          ! du3/dr - u3/R
      endif        
      call col3(wk2,u2i,rinv,ntot1)        ! u2(im)/R
      call add2s2(errt,wk2,-k_f3d,ntot1)   ! du3/dr - u3/R - k/R*u2(im)
      call cmult(errt,0.5,ntot1)           ! 0.5*[du3/dr - u3/R - k/R*u2(im)]

!!    Ei_\thetaR
!     (Imaginary)      
      call copy(eirt,ui3r,ntot1)           ! du3(im)/dr
      if (ifcyl_f3d) then
        call col3(wk1,u3i,rinv,ntot1)      ! u3(im)/R
        call sub2(eirt,wk1,ntot1)          ! du3(im)/dr - u3(im)/R
      endif        
      call col3(wk2,u2r,rinv,ntot1)        ! u2/R
      call add2s2(eirt,wk2,k_f3d,ntot1)    ! du3(im)/dr - u3(im)/R + k*u2/R
      call cmult(eirt,0.5,ntot1)           ! 0.5*[du3(im)/dr - u3(im)/R + k*u2/R]


!!    Er_\theta\theta
!     (Real)      
      call col3(ertt,u3i,rinv,ntot1)       ! u3(im)/R      
      call cmult(ertt,-k_f3d,ntot1)        ! -k*u3(im)/R
      if (ifcyl_f3d) then
        call xaddcol3(ertt,u2r,rinv,ntot1) ! -k*u3(im)/R + u2/R
      endif  

!!    Ei_\theta\theta
!     (Imaginary)      
      call col3(eitt,u3r,rinv,ntot1)       ! u3/R      
      call cmult(eitt,k_f3d,ntot1)         ! k*u3/R
      if (ifcyl_f3d) then
        call xaddcol3(eitt,u2i,rinv,ntot1) ! k*u3/R + u2(im)/R
      endif


      return
      end subroutine stnrate_cyl_real
!-----------------------------------------------------------------------
      subroutine div_stress_cyl_real (Au1r,Au2r,Au3r,nel)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'GEOM'

      include 'F3D'

!!     Real Variables      
!      real erxt,errt,erxx,erxr,errr,ertt
!      common /scrpsn3/ erxt(lx1*ly1*lz1*lelt)      ! Er_x\theta
!     $               , errt(lx1*ly1*lz1*lelt)      ! Er_Rt
!     $               ,  erxx(lx1*ly1*lz1*lelt)      ! Er_xx
!     $               , erxr(lx1*ly1*lz1*lelt)      ! Er_xR
!     $               , errr(lx1*ly1*lz1*lelt)      ! Er_RR
!     $               , ertt(lx1*ly1*lz1*lelt)      ! Er_\theta\theta
!
!!     Made a new scratch array here.      
!!     Imaginary Variables
!      real eixt,eirt,eixx,eixr,eirr,eitt,wrk
!      common /scrpsn4/  eixt(lx1*ly1*lz1*lelt)      ! Ei_x\theta
!     $                , eirt(lx1*ly1*lz1*lelt)      ! Ei_Rt
!     $                , eixx(lx1*ly1*lz1*lelt)      ! Ei_xx
!     $                , eixr(lx1*ly1*lz1*lelt)      ! Ei_xR
!     $                , eirr(lx1*ly1*lz1*lelt)      ! Ei_RR
!     $                , eitt(lx1*ly1*lz1*lelt)      ! Ei_\theta\theta
!     $                , wrk(lx1*ly1*lz1*lelt)       ! work

      real rinv(lx1*ly1*lz1*lelv),wk1(lx1*ly1*lz1*lelv)
      real wk2(lx1*ly1*lz1*lelv),hii(lx1*ly1*lz1*lelv)
      common /scrpsn5/ rinv,wk1,wk2,hii


      real Au1r,Au2r,Au3r
      dimension Au1r(lx1,ly1,lz1,1)
     $        , Au2r(lx1,ly1,lz1,1)
     $        , Au3r(lx1,ly1,lz1,1)
     
      integer nel,ntot1


      ntot1 = lx1*ly1*lz1*nel

!     Since we don't pass the imaginary component.
!     I force the ei** arrays to zero. Just in case precision erros
!     have piled up.

      call rzero3(eixt,eirt,eixx,ntot1)
      call rzero3(eixr,eirr,eitt,ntot1)

!     rinv = 1/R
      if (ifcyl_f3d) then 
        call invers2(rinv,ym1,ntot1)
      else
        call rone(rinv,ntot1)
      endif  
      call col2c(rinv,bm1,-k_f3d,ntot1)   ! rinv = -k*BM1/R

!     Real Variables      
      call ttxyz(Au1r,erxx,erxr,erxt,nel)  ! [erxx*dv/dx + erxr*dv/dr]*BM1
      call col3(wk1,eixt,rinv,ntot1)       ! eixt*k*BM1/R    
      call sub2(Au1r,wk1,ntot1)

      call ttxyz(Au2r,erxr,errr,errt,nel)  ! [erxr*dv/dx + errr*dv/dr]*BM1
      call col3(wk1,eirt,rinv,ntot1)       ! eirt*k*BM1/R    
      call sub2(Au2r,wk1,ntot1)
!     Additional terms from vector gradient of test functions
      if (ifcyl_f3d) then
        call invcol2(wk1,ym1,ntot1)
        call col2c(wk1,bm1,-1.0,ntot1)       ! -BM1*v/R
        call Xaddcol3(Au2r,errt,wk1,ntot1)   ! -BM1*v/R*[errt]
      endif        

      call ttxyz(Au3r,erxt,errt,ertt,nel)  ! [erxt*dv/dx + errt*dv/dr]*BM1
      call col3(wk1,eitt,rinv,ntot1)       ! eitt*k*BM1/R    
      call sub2(Au3r,wk1,ntot1)
!     Additional terms from vector gradient of test functions
      if (ifcyl_f3d) then
        call invcol2(wk1,ym1,ntot1)
        call col2(wk1,bm1,ntot1)             ! BM1*v/R
        call Xaddcol3(Au3r,ertt,wk1,ntot1)   ! BM1*v/R*[ertt]
      endif        

!     Imaginary Variables (Just deleted). 


      return
      end subroutine div_stress_cyl_real
!-----------------------------------------------------------------------







