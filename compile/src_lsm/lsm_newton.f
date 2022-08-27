!======================================================================
!     Newton solver for Level Set Method reinitialization.
!     Author: Prabal S. Negi
!
!====================================================================== 
!-----------------------------------------------------------------------
      subroutine lsm_jacobian_action(Ax,SSign,px,py,pz,phip,ifqqt)

      implicit none

      include 'SIZE'
      include 'MASS'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

!     Action of the Jacobian
      real Ax(lv)

!     Smooth Sign
      real SSign(lv)

!     gradients of perturbations      
      real px(lv)
      real py(lv)
      real pz(lv)

      real phip(lv)         ! perturbation

!     pert gradients      
      real ppx(lv)
      real ppy(lv)
      real ppz(lv)
      real ppt(lv)
!      common /scruz/ ppx,ppy,ppz,ppt

!     Work arrays      
      real tb1(lv)
      real tb2(lv)
      real tb3(lv)
      real tb4(lv)
!      common /scrmg/ tb1,tb2,tb3,tb4

      logical ifqqt           ! just local actions?

      integer i,n
     
      n = lx1*ly1*lz1*nelv 

      call gradm1(ppx,ppy,ppz,phip) 
      call rzero(ppt,n)

      do i=1,n
!       grad\phi.grad\phi     ! in principle should never be zero      
        tb1(i) = px(i)**2 + py(i)**2
        if (ndim.eq.3) tb1(i) = tb1(i)+pz(i)**2

!       grad\phi.grad\phi'
        ppt(i) = px(i)*ppx(i) + py(i)*ppy(i)
        if (ndim.eq.3) ppt(i) = ppt(i) + pz(i)*ppz(i)

        tb2(i) = ppt(i)/tb1(i)
      enddo

      call col3(tb3,SSign,tb2,n)
      call chsign(tb3,n)
      call col2(tb3,bm1,n)
      if (ifqqt) then
        call invertB(Ax,tb3)
      else
        call copy(Ax,tb3,n)
      endif  

!     Ax is now: Ax = -S0*(\grad\phi.\grad\phi')/(\grad\phi.\grad\phi)

      return
      end subroutine lsm_jacobian_action      
!---------------------------------------------------------------------- 
      
      subroutine lsm_newton_method()

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      include 'TEST'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real F0(lv)           ! this is the residual
      real SSign(lv)

      real Seps

      real phix(lv)
      real phiy(lv)
      real phiz(lv)
!      common /scrsf11/ phix,phiy,phiz

      real deltaphi(lv)

!     Work arrays      
      real tb1(lv)
      real tb2(lv)
      real tb3(lv)
      real tb4(lv)
!      common /scrmg/ tb1,tb2,tb3,tb4

      integer i,j,n

      real gdn

      real op_glsc2_wt
      real distn

      logical ifqqt
      logical ifprop          ! use propagator

      ifield = 2

      n = lx1*ly1*lz1*nelv

!     Calculate gradient of distance function
      call rzero3(phix,phiy,phiz,n) 
      call gradm1(phix,phiy,phiz,t)

!     Create Smooth sign fuction
      Seps = 1.0e-1
      do i=1,n
        gdn = sqrt(phix(i)**2 + phiy(i)**2 + phiz(i)**2)
        SSign(i) = t(i,1,1,1,1)/(sqrt(t(i,1,1,1,1)**2 
     $                           + (gdn*Seps)**2))    ! sign function
      enddo

      call copy(tmp1,t,n)     ! make a copy
      do j=1,20
!       Calculate gradient of distance function
        call rzero3(phix,phiy,phiz,n) 
        call gradm1(phix,phiy,phiz,t)

        distn = sqrt(op_glsc2_wt(phix,phiy,phiz,
     $              phix,phiy,phiz,bm1)/volvm1)
        write(6,*) '|Newton Distance|: ', j, distn
   
        if (ifprop) then
          nsteps = 100
          call copy(F0,t,n)
          call nllsm_march(F0,nsteps)
          call sub2(F0,t,n)         ! (F(x) - x)
          call chsign(F0,n)         ! -F0
        else  
          do i=1,n
            gdn = 1.0 - sqrt(phix(i)**2 + phiy(i)**2 + phiz(i)**2)
            tb1(i) = SSign(i)*gdn
          enddo
          call col2(tb1,bm1,n)
          call chsign(tb1,n)             
          call invertB(F0,tb1)          ! Residual (C0 continuous)
        endif  

        call lsm_gmres(deltaphi,F0,SSign,phix,phiy,phiz,ifprop)
        ifqqt = .true.
        call lsm_jacobian_action(tb4,SSign,phix,phiy,phiz,deltaphi,
     $                           ifqqt)
!        call outpost(phix,phiy,phiz,pr,SSign,'grd')
!        call outpost(deltaphi,F0,vz,pr,tb4,'new')
        call add2(t,deltaphi,n)
      enddo  
!     reset temperature back
      call copy(t,tmp1,n)

      return
      end subroutine lsm_newton_method
!---------------------------------------------------------------------- 

      subroutine lsm_gmres(soln,resi,SSign,px,py,pz,ifprop)

      implicit none

      include 'SIZE'
      include 'MASS'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real resi(lv)         ! Initial residual
      real SSign(lv)
      real soln(lv)         ! Solution
     
      real px(lv)
      real py(lv)
      real pz(lv)

      integer lkryl,lkryl1
      parameter (lkryl=50)
      parameter (lkryl1=lkryl+1)

      real lls_res(lkryl1)           ! residual for linear least squares
      real lls_sol(lkryl1)           ! solution for linear least squares
      real lls_wk1(lkryl1)           ! work array
      real lls_wk2(lkryl1)           ! work array

      real KrylV(lv,lkryl1)
      integer nkryl
      real scvec(lkryl),scvec1(lkryl1)
      real Hessen(lkryl1,lkryl)
      real Hes_wk(lkryl1,lkryl)
      real Hes_wk2(lkryl1,lkryl)
      real Hes_wk3(lkryl1,lkryl)

!     Work arrays      
      real tb1(lv)
      real tb2(lv)
      real tb3(lv)
      real tb4(lv)
!      common /scrmg/ tb1,tb2,tb3,tb4

      integer i,j,k,ngs,n,l
      integer ik
      
      real gl2norm      ! function. overwrites scrsf
      real gln,gln1,gln0

      real glsc3        ! function
      real sc

!     parameters for lapack      
      integer m         ! no rows in matrix
!      integer n         ! no columns in matrix       ! already declared
      integer lda       ! leading dimension of the matrix
      integer ldb       ! leading dimension of rhs vector/matrix
      integer nrhs      ! no of RHS vectors

      logical ifqqt
      logical ifprop

      integer nsteps
      
      l = lx1*ly1*lz1*lelv
      n = lx1*ly1*lz1*nelv

      call copy(tb1,resi,n)
     
      gln = gl2norm(tb1,n)

      call rzero(krylV,l*lkryl1)         ! initialize Krylov space

      call cmult(tb1,1.0/gln,n)           ! normalize starting vector
      call copy(krylV(1,1),tb1,n)         ! initialize Krylov space
      call rzero(lls_res,lkryl1)
      lls_res(1) = gln
      call copy(lls_sol,lls_res,lkryl1)
     
      nkryl = 3
      ngs   = 1         ! no of Gram-Schmidt

      call rzero(Hessen,lkryl*lkryl)
      ik = 1
      do i=1,nkryl 
!       tb2 = Jac*tb1     
        gln0 = gl2norm(tb1,n)
        if (ifprop) then
          nsteps = 100
          call copy(tb2,tb1,n)
          call rone(tb1,n)
          call linlsm_march(tb2,SSign,px,py,pz,nsteps)
          call sub2(tb2,tb1,n) 
        else
          ifqqt = .true.
          call lsm_jacobian_action(tb2,SSign,px,py,pz,tb1,ifqqt)
        endif  
        gln1 = gl2norm(tb2,n)
        do k=1,ngs
          call rzero(tb3,n)
          call rzero(scvec1,lkryl)
          do j=1,i
            sc       = glsc3(krylV(1,j),tb2,bm1,n)
            scvec(j) = sc/volvm1
            call add2s2(tb3,krylV(1,j),scvec(j),n)
            scvec1(j) = scvec1(j)+scvec(j)
          enddo
          call sub2(tb2,tb3,n)
        enddo  
        gln = gl2norm(tb2,n)

        if (gln.le.1.0e-10) then
          exit
        else
          ik = i+1
          scvec1(ik) = gln
          call copy(Hessen(1,ik-1),scvec1,ik)
          call copy(tb1,tb2,n) 
          call cmult(tb1,1.0/gln,n)
          call copy(krylV(1,ik),tb1,n)
        endif  
      enddo

      m     = ik
      n     = ik-1
      lda   = lkryl1
      ldb   = lkryl1
      nrhs  = 1


      call copy(Hes_wk,Hessen,lkryl1*lkryl)
      call copy(Hes_wk2,Hessen,lkryl1*lkryl)
    
!     Testing inverse      
      call wrp_dgeinv(Hes_wk,lda,n) 
      call mxm(Hes_wk,lda,lls_res,lda,lls_wk1,1)

      write(6,*) 'Res1:', (lls_res(i),i=1,n)
      write(6,*) 'Sol1:', (lls_wk1(i),i=1,n)
!     Temp. Get Ax      
      call mxm(Hes_wk2,lda,lls_wk1,n,lls_wk2,nrhs)
      write(6,*) 'Ax1 :', (lls_wk2(i),i=1,n)

      call mxm(krylV,l,lls_wk1,n,soln,1)



      return
      end subroutine lsm_gmres
!---------------------------------------------------------------------- 

      subroutine local_qr(q,r,ldr,m,n)

      implicit none

      integer ldr             ! array sizes
      integer n               ! no of columns in the matrix
      integer m               ! no of rows in the matrix

      real q(ldr,n)
      real r(ldr,n)           ! input matrix 
      real v(ldr)
      real w(ldr)

      real vlsc2
      real s

      integer i,j,k
      integer ngs       ! Gram-Schmidt cycles

!      n=lda
      ngs = 1

      call rzero(q,ldr*n)

      do i=1,n
        call copy(v,r(1,i),m)
        call rzero(r(1,i),m)
        do k=1,ngs
          call rzero(w,m)
          do j=1,i-1  
           s = vlsc2(q(1,j),v,m)            ! projection
           r(j,i) = r(j,i) + s
           call add2s2(w,q(1,j),s,m)        ! gather projections
          enddo
          call sub2(v,w,m)
        enddo
        s = sqrt(vlsc2(v,v,m))
        r(i,i) = s
        call cmult(v,1.0/s,m) 
        call copy(q(1,i),v,m)
      enddo        

      return
      end subroutine local_qr
!---------------------------------------------------------------------- 

      subroutine linlsm_rk4(phi,SSign,px,py,pz)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'TSTEP'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real SSign(lv)
      real dt2

!     gradients of base field
      real px(lv)
      real py(lv)
      real pz(lv)

      real phi(lv)            ! evolving field

      real ddt(lv,4)
!      common /scrmg/ ddt

      real ta4(lv,4)
!      common /scruz/ ta4

      real fil(lv) 
!      common /scrcg/ fil

      real gdn

      integer i,n
      logical iffil
      logical ifqqt


      ifield = 2

      iffil = .false.
      ifqqt = .false.

      n = lx1*ly1*lz1*nelv

      dt2 = 1.0*dt

      call copy(ta4(1,4),phi,n)

!     lsm_jacobian_action(Ax,SSign,px,py,pz,phip)
!     1st step
      call lsm_jacobian_action(ddt(1,1),SSign,px,py,pz,ta4(1,4),ifqqt)
      call col2(ddt(1,1),bm1,n)
      if (iffil) then
        call hpf_dist(fil,ta4(1,4))
        call add2col2(ddt(1,1),fil,bm1,n)           ! already has negative sign
      endif        
      call copy(fil,ddt(1,1),n)
      call invertB(ta4(1,4),fil)
      call add2s1(ta4(1,4),phi,dt2/2.0,n)

!     2nd step
      call lsm_jacobian_action(ddt(1,2),SSign,px,py,pz,ta4(1,4),ifqqt)
      call col2(ddt(1,2),bm1,n)
      if (iffil) then
        call hpf_dist(fil,ta4(1,4)) 
        call add2col2(ddt(1,2),fil,bm1,n)
      endif 
      call copy(fil,ddt(1,2),n)
      call invertB(ta4(1,4),fil)
      call add2s1(ta4(1,4),phi,dt2/2.0,n)

!     3rd step
      call lsm_jacobian_action(ddt(1,3),SSign,px,py,pz,ta4(1,4),ifqqt)
      call col2(ddt(1,3),bm1,n)
      if (iffil) then 
        call hpf_dist(fil,ta4(1,4))
        call add2col2(ddt(1,3),fil,bm1,n)
      endif  
      call copy(fil,ddt(1,3),n)
      call invertB(ta4(1,4),fil)
      call add2s1(ta4(1,4),phi,dt2,n)

!     4th step
      call lsm_jacobian_action(ddt(1,4),SSign,px,py,pz,ta4(1,4),ifqqt)
      call col2(ddt(1,4),bm1,n)
      if (iffil) then
        call hpf_dist(fil,ta4(1,4)) 
        call add2col2(ddt(1,4),fil,bm1,n)
      endif  
      do i=1,n
        fil(i) = (ddt(i,1) + 2.0*ddt(i,2) 
     $              + 2.0*ddt(i,3) + ddt(i,4))
      enddo  
      call invertB(ta4(1,4),fil)
      call add2s2(phi,ta4(1,4),dt2/6.0,n)


      return
      end subroutine linlsm_rk4
!---------------------------------------------------------------------- 

      subroutine linlsm_march(phi,SSign,px,py,pz,nsteps)

      implicit none

      include 'SIZE'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real SSign(lv)

!     gradients of base field
      real px(lv)
      real py(lv)
      real pz(lv)

      real phi(lv)      ! evolving field

      integer i,nsteps

      do i=1,nsteps
        call linlsm_rk4(phi,SSign,px,py,pz)
      enddo

      return
      end subroutine linlsm_march
!----------------------------------------------------------------------

      subroutine nllsm_march(phi,nsteps)

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real phi(lv)

      real ta1(lv)
      real ta2(lv)
      real ta3(lv)
      real ta4(lv)
!      common /scruz/ ta1,ta2,ta3,ta4

      real SSign(lv), dummy(lv)
!      common /scrch/ SSign(lv), dummy(lv)

      integer i,j
      integer n
      integer nsteps

      real Seps         ! smoothening for Sign
      real gdn


      n = lx1*ly1*lz1*nelv

!     Create Smooth sign fuction
!     This only at Step 0      
      Seps = 1.0e-1
      call gradm1(ta1,ta2,ta3,phi)
      do j=1,n
        gdn = sqrt(ta1(j)**2 + ta2(j)**2 + ta3(j)**2)
        SSign(j) = phi(j)/(sqrt(phi(j)**2 + (gdn*Seps)**2))    ! sign function
      enddo

      do i=1,nsteps
!        call lsm_rk2(phi,SSign)
        call lsm_rk4(phi,SSign)
!        call lsm_ExpEuler(phi,SSign)
      enddo  

      return
      end subroutine nllsm_march
!---------------------------------------------------------------------- 



