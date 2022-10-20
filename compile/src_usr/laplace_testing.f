!======================================================================
!     Author: Prabal Negi
!     Description: Testing Solver/Preconditioner for Laplace.
!
!======================================================================       
      subroutine pseudolaplace_arnoldi()

      implicit none
  
      include 'SIZE'
      include 'INPUT'         ! uparam/param
      include 'SOLN'          ! vtrans
      include 'MASS'          ! bm2
      include 'TSTEP'         ! ifield
      include 'GEOM'          ! xm2

      include 'ARN_ARPD'
      include 'TSTEPPERD'
 
      include 'TEST'

      real h1,h2,h2inv
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)

      integer i,ntot1,ntot2
      integer intype
      integer iflg,j

      real rnd
      real rad

      real lambda

      integer igmres
      real x,y,z

      logical ifcylindrical

      if (nio.eq.0) write(6,*) 'Pseudo Laplacian Arnoldi'

      ifield = 1
      istep  = 0

      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv

!      call rone  (vtrans,ntot1)

      intype = 1
!      call rzero   (h1,ntot1)
!      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
!      call invers2 (h2inv,h2,ntot1)
!      call reset_preconditioner() 

      do i=1,ntot2
        call random_number(rnd)
        x        = xm2(i,1,1,1)
        y        = ym2(i,1,1,1)
        z        = zm2(i,1,1,1)
        prp(i,1) = 1.0*rnd ! sin(x)*sin(y)*sin(z) ! 1.0*rnd
      enddo

      ifflow = .false.
      ifheat = .false.
!      call tst_init()   ! also calls arn_init()

      param(63) = 1     ! 8 byte output

      istep = 0

      ifcylindrical = .false.
      if (ifcylindrical) call col2(bm1,ym1,ntot1)     ! for (B^-1)

      do while (istep.lt.nsteps)

        istep = istep+1
        call settime
        call setdt
      
        if (nio.eq.0) write(6,*) 'Iteration:', istep,dt 

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

!        call rk4_advance(prp,dt)

!!       tmp = (U^T)*E*(U)*p
!        call ortho_right(prp) 
!        call cdabdtp(tmp4,prp,h1,h2,h2inv,intype)
!        call ortho_left(tmp4)
!!       p = (I + \lambda*E)*p
!        lambda = -1.0*dt
!        call add2s2(prp,tmp4,lambda,ntot2)

!!       tmp = (U^T)*E*(U)*(M^-1)*p        
!!        if (tst_vstep.eq.0.and.tst_istep.eq.0) call ortho_new(prp)
!        call eprec2_new(tmp4,prp)
!        call ortho_new(tmp4)
!        call copy(tmp8,tmp4,ntot2)
!        call cdabdtp(tmp4,tmp8,h1,h2,h2inv,intype)
!        call ortho_new(tmp4)
!!       p = (I + \lambda*E)*p
!        lambda = -1.0*dt
!        call add2s2(prp,tmp4,lambda,ntot2)


!!       tmp = E*(U^T)*(M^-1)*(U)*p
!!        call ortho_right(prp) 
!        call eprec2_new(tmp4,prp)
!!        call ortho_left(tmp4)
!!        call ortho_right(tmp4)
!        call copy(tmp8,tmp4,ntot2)
!        call cdabdtp(tmp4,tmp8,h1,h2,h2inv,intype)
!!        call ortho_left(tmp4)
!!       p = (I + \lambda*E)*p
!        lambda = -1.0*dt
!        call add2s2(prp,tmp4,lambda,ntot2)

!!       tmp = (U^T)E*(M^-1)*(U)*p
!        call ortho_new(prp) 
!        call eprec2_new(tmp4,prp)
!        call copy(tmp8,tmp4,ntot2)
!        call cdabdtp(tmp4,tmp8,h1,h2,h2inv,intype)
!        call ortho_new(tmp4)
!!       p = (I + \lambda*E)*p
!        lambda = -1.0*dt
!        call add2s2(prp,tmp4,lambda,ntot2)

!!       local_solves_fdm
!        call local_solves_fdm(tmp4,prp)
!!        call new_fdm(tmp4,prp)
!!       p = (I + \lambda*E)*p
!        lambda = -1.0*dt
!        call add2s2(prp,tmp4,lambda,ntot2)


!!       What Nek Does.
!!       tmp = E*(U)*(M^-1)*p;
!!       Remove null space only from initial vector.
!        if (tst_vstep.eq.0.and.tst_istep.eq.0) call ortho_new(prp)
!        call eprec2_new(tmp4,prp)
!        call ortho_new(tmp4)
!        call copy(tmp8,tmp4,ntot2)
!        call cdabdtp(tmp4,tmp8,h1,h2,h2inv,intype)
!!       p = (I + \lambda*E)*p
!        lambda = -1.0*dt
!        call add2s2(prp,tmp4,lambda,ntot2)



!       For inf-sup
!       tmp = E*p;
!       Remove null space only from initial vector.
!        if (tst_vstep.eq.0.and.tst_istep.eq.0) then
!          call cdabdtp(tmp4,prp,h1,h2,h2inv,intype)
!          call copy(prp,tmp4,ntot2)
!          call outpost(vx,vy,vz,prp,t,'ini')
!        endif
!        call cM1dabdtp(tmp4,prp,h1,h2,h2inv,intype)
!        call col3(tmp8,prp,bm2,ntot2)
!        call cdabdtp(tmp4,prp,h1,h2,h2inv,intype)
!        call cdabdtp_cyl(tmp4,prp,h1,h2,h2inv,intype)
!        call crhodabdtp(tmp4,prp,h1,h2,h2inv,intype)
        call cdab_full_dtp(tmp4,prp,h1,h2,h2inv,intype)
!        call invcol2(tmp4,bm2,ntot2)      ! (B^-1)E
!        call cM1dabdtp (tmp4,prp,h1,h2,h2inv,intype)
!       p = (I + \lambda*E)*p
        lambda = -uparam(7)
        call add2s2(prp,tmp4,lambda,ntot2)


!!       Preconditioned S
!        call cdabdtp(tmp4,prp,h1,h2,h2inv,intype)
!        call invcol2(tmp4,bm2,ntot2)      ! (B^-1)E
!!       p = (I + \lambda*(B^-1)*E)*p
!        lambda = -uparam(7)
!        call add2s2(prp,tmp4,lambda,ntot2)


!        if (tst_vstep.eq.0.and.tst_istep.eq.0) call tst_init()
!        call tst_solve()    

        if (lastep.eq.1) istep = nsteps

      enddo

14    format(A5,2x,16(E12.5,2x))


      return
      end
c-----------------------------------------------------------------------

      subroutine rk4_pseudolaplace(p,h1,h2,h2inv,intype,dt)

      implicit none

      include 'SIZE'
      include 'MASS'          ! bm2
      include 'PARALLEL'      ! nio

      integer lt,lt2
      parameter (lt  = lx1*ly1*lz1*lelt)
      parameter (lt2 = lx2*ly2*lz2*lelt)

      real p1,p2,p3
      real Ep,Ep1,Ep2,Ep3
      common /scrns/ p1  (lt)
     $ ,             p2  (lt)
     $ ,             p3  (lt)
     $ ,             Ep  (lt)
     $ ,             Ep1 (lt)
     $ ,             Ep2 (lt)
     $ ,             Ep3 (lt)

      real h1    (lt)
      real h2    (lt)
      real h2inv (lt)

      real p(lt2)
      integer intype

      real s0,s1,s2,s3
      real s

      integer n1,n2
      real dt
      real visc

      if (nio.eq.0) write(6,*) 'RK4', dt

      n1  = lx1*ly1*lz1*nelv
      n2  = lx2*ly2*lz2*nelv

      visc = 1.0

      s0  = 1.0
      s1  = visc*dt/2.0
      s2  = visc*dt/2.0
      s3  = visc*dt/1.0

      call cdabdtp(Ep,p,h1,h2,h2inv,intype)     ! Ep  = E*p
      call invcol2(Ep,bm2,n2)                   ! Ep  = (B^-1)*E*p
      call add3s2(p1,p,Ep,s0,s1,n2)             ! p1  = p - (dt/2)*E*p

      call cdabdtp(Ep1,p1,h1,h2,h2inv,intype)   ! Ep1 = E*p1
      call invcol2(Ep1,bm2,n2)                  ! Ep1 = (B^-1)*E*p1
      call add3s2(p2,p,Ep1,s0,s2,n2)            ! p2  = p - (dt/2)*E*p1

      call cdabdtp(Ep2,p2,h1,h2,h2inv,intype)   ! Ep2 = E*p2
      call invcol2(Ep2,bm2,n2)                  ! Ep2 = (B^-1)*E*p2
      call add3s2(p3,p,Ep2,s0,s3,n2)            ! p3  = p - (dt)*E*p2

      call cdabdtp(Ep3,p3,h1,h2,h2inv,intype)   ! Ep3 = E*p3
      call invcol2(Ep3,bm2,n2)                  ! Ep3 = (B^-1)*E*p3
   
      call add2s2(Ep,Ep1,2.0,n2)    ! Ep = E*p + 2*E*p1
      call add2s2(Ep,Ep2,2.0,n2)    ! Ep = E*p + 2*E*p1 + 2*E*p2
      call add2s2(Ep,Ep3,1.0,n2)    ! Ep = E*p + 2*E*p1 + 2*E*p2 + E*p3

      s =  visc*dt/6.0
      call add2s2(p,Ep,s,n2)
    
      return
      end subroutine rk4_pseudolaplace
!---------------------------------------------------------------------- 
      subroutine rk4_advance(v,dt)

      implicit none

      include 'SIZE'
      include 'SOLN'    ! vtrans
      include 'MASS'    ! bm2

      integer lt,lt2
      parameter (lt  = lx1*ly1*lz1*lelt)
      parameter (lt2 = lx2*ly2*lz2*lelt)

      real vi
      real Ev
      common /scrns_new/ vi (lt,4)
     $ ,                 Ev (lt,4)

      real v(1)         ! input/output

      real wk1,wk2,wk3
      common /scruz/ wk1(lt2)
     $             , wk2(lt2)
     $             , wk3(lt2)


      real s0
      real si(4),s

      integer n1,n2,n
      real dt
      real visc

      integer i

      real h1    (lt)
      real h2    (lt)
      real h2inv (lt)
      integer intype

      integer ifld

      intype = 1

      n1  = lx1*ly1*lz1*nelv
      n2  = lx2*ly2*lz2*nelv

      n   = n2

      ifld = 1
      call rzero   (h1,n1)
      call copy    (h2,vtrans(1,1,1,1,ifld),n1)
!      call rone    (h2,n1)
      call invers2 (h2inv,h2,n1)

      visc =  1.0

      s0    = 1.0

      si(1) = 2.0 
      si(2) = 2.0
      si(3) = 1.0
      si(4) = 6.0


      call rzero(Ev,lt*4)
      call rzero(vi,lt*4)

      call copy(vi(1,1),v,n)
      call copy(wk1,v,n)
!      call col2(wk1,bm2,n)         ! Mass Matrix
      
      do i=1,3

!       M*vi
        call cdabdtp(Ev(1,i),vi(1,i),h1,h2,h2inv,intype)    ! Ev_i = E*v_i
!        call cmult(Ev(1,i),visc,n)
!        call col2(Ev(1,i),bm2inv,n)

!       vi+1 = v + (dt*fac)*M*vi
        s = visc*dt/si(i)
        call add3s2(vi(1,i+1),wk1,Ev(1,i),s0,s,n)  ! v_i+1 = Bv + (dt/2)*E*v_i
!        call invcol2(vi(1,i+1),bm2,n)              ! v_i+1 = (B^-1)*v_i+1
      enddo

      call cdabdtp(Ev(1,4),vi(1,4),h1,h2,h2inv,intype)   ! Ev4 = E*p4
!      call cmult(Ev(1,4),visc,n)
!      call col2(Ev(1,4),bm2inv,n)                        ! Ev4 = (B^-1)*E*p4

      do i=1,3
        s = si(i)
        call add2s2(Ev(1,1),Ev(1,i+1),s,n)
      enddo
      
      s = visc*dt/si(4)
      call add3s2(v,wk1,Ev(1,1),s0,s,n)
!      call invcol2(v,bm2,n)               ! v_i+1 = (B^-1)*v_i+1
   
      return
      end subroutine rk4_advance
!---------------------------------------------------------------------- 
      subroutine bdfk_advance(v,istep,dt)

      implicit none

      include 'SIZE'
      include 'SOLN'    ! vtrans
      include 'MASS'    ! bm2

      integer lt,lt2
      parameter (lt  = lx1*ly1*lz1*lelt)
      parameter (lt2 = lx2*ly2*lz2*lelt)

      real vi
      real Ev
      common /SOLN_NEW/ vi (lt,4)
     $ ,                Ev (lt,4)

      real v(1)         ! input/output
      integer istep     ! for bdfk/extk

      real wk1,wk2,wk3
      common /scruz/ wk1(lt2)
     $             , wk2(lt2)
     $             , wk3(lt2)

      integer n1,n2,n
      real dt
      real visc

      integer i

      real h1    (lt)
      real h2    (lt)
      real h2inv (lt)
      integer intype

      integer ifld

      real bdfk(4)
      real extk(4)
      real const

      intype = 1

      n1  = lx1*ly1*lz1*nelv
      n2  = lx2*ly2*lz2*nelv

      n   = n2

      ifld = 1
      call rzero   (h1,n1)
      call copy    (h2,vtrans(1,1,1,1,ifld),n1)
!      call rone    (h2,n1)
      call invers2 (h2inv,h2,n1)

      visc =  1.0-6

      call copy(vi(1,1),v,n2)

      call get_bdfextk(bdfk,extk,istep)

      call cdabdtp(Ev(1,1),v,h1,h2,h2inv,intype)    ! Ev = E*v
      call cmult(Ev(1,1),visc,n2)

!     Extrapolate
      call cmult2(Ev(1,4),Ev(1,1),extk(2),n2)
      call add2s2(Ev(1,4),Ev(1,2),extk(3),n2)
      call add2s2(Ev(1,4),Ev(1,3),extk(4),n2)

!     Lag rhs
      call copy(Ev(1,3),Ev(1,2),n2)
      call copy(Ev(1,2),Ev(1,1),n2)

!     Backward difference terms
      call rzero(wk2,n2)
      do i=1,3
!        call col3(wk1,bm2,vi(1,i),n2)
!        const = -bdfk(i+1)/dt
!        call add2s2(wk2,wk1,const,n2)

        const = -bdfk(i+1)/dt
        call add2s2(wk2,vi(1,i),const,n2)
      enddo 

!     Lag soln
      call copy(vi(1,3),vi(1,2),n2)
      call copy(vi(1,2),vi(1,1),n2)

      call add2(Ev(1,4),wk2,n2)
!      call col3(v,Ev(1,4),bm2inv,n2)
      call copy(v,Ev(1,4),n2)
      const = dt/bdfk(1)
      call cmult(v,const,n2)
   
      return
      end subroutine bdfk_advance
!---------------------------------------------------------------------- 

      subroutine get_bdfextk(bdf,ext,i)

      implicit none

!     Backward difference coeff.
      real bdfk(4,3)
      real bdf(4)

!     Extrapolation Coefficients
      real extk(4,3)
      real ext(4)

      integer i      ! istep
      integer bdfo   ! bdf order
      integer exto   ! ext order

!-------------------------------------------------- 
!     Backward Difference
!     First order
      bdfk(1,1) =  1.0
      bdfk(2,1) = -1.0
      bdfk(3,1) =  0.0
      bdfk(4,1) =  0.0
      
!     Second order
      bdfk(1,2) =  3.0/2.0
      bdfk(2,2) = -4.0/2.0
      bdfk(3,2) =  1.0/2.0
      bdfk(4,2) =  0.0

!     Third order
      bdfk(1,3) =  11.0/6.0
      bdfk(2,3) = -18.0/6.0
      bdfk(3,3) =   9.0/6.0
      bdfk(4,3) = - 2.0/6.0

!     Extrapolation
!     First order Explicit 
      extk(1,1) = 0.0
      extk(2,1) = 1.0
      extk(3,1) = 0.0
      extk(4,1) = 0.0

!     Second order Explicit 
      extk(1,2) =  0.0
      extk(2,2) =  2.0
      extk(3,2) = -1.0
      extk(4,2) =  0.0

!     Third order Explicit 
      extk(1,3) =  0.0
      extk(2,3) =  3.0
      extk(3,3) = -3.0
      extk(4,3) =  1.0
!-------------------------------------------------- 

      if (i.lt.1) then
        call rzero(bdf,4)
        call rzero(ext,4)
        return
      endif

      bdfo = min(i,3)
      exto = min(i,3)
        
      call copy(bdf,bdfk(1,bdfo),4)
      call copy(ext,extk(1,exto),4)



      return
      end subroutine get_bdfextk
!---------------------------------------------------------------------- 
      subroutine laplace_test()

      implicit none
  
      include 'SIZE'
      include 'INPUT'         ! uparam/param
      include 'SOLN'          ! vtrans
      include 'MASS'          ! bm2
      include 'TSTEP'         ! ifield
      include 'GEOM'          ! xm2
      include 'DOMAIN'        ! lcr,lxc
      include 'TEST'

      real w1,w2,w3
      real dv1,dv2,dv3
      real dp
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

      integer i,ntot1,ntot2
      integer intype
      integer iflg,j

      real uc,w
      common /scrpre/ uc(lcr*lelt),w(2*lx1*ly1*lz1)

      real rand
      integer seed
      parameter (seed=86456)
      real rnd
      real rad

      integer igmres


      if (nio.eq.0) write(6,*) 'Testing Laplace'

      ifield = 1
      istep  = 2

      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv

!     Preconditioner
      param(42)=uparam(8)       ! 0: GMRES (nonsymmetric), 1: PCG w/o weights
      param(43)=uparam(9)       ! 0: Additive multilevel (param 42=0), 1: Original 2 level
      param(44)=uparam(10)      ! 0: E based Schwartz (FEM), 1: A based Schwartz
      
!      call rone(vtrans,ntot1) 
      call reset_preconditioner()

      intype = 1        ! explicit

      call rzero   (h1,ntot1)
      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call rone    (h2,ntot1)
      call invers2 (h2inv,h2,ntot1)

      do i=1,ntot2
        call random_number(rnd)
        if (ifcyclic) then
          rad = sqrt(ym2(i,1,1,1)**2 + zm2(i,1,1,1)**2)
        else
          rad = ym2(i,1,1,1)
        endif
!        dp(i,1,1,1) = (1.0e-0)*rad + 0.0*rnd
        dp(i,1,1,1) = (1.0e-0)*xm2(i,1,1,1) + 0.0*rnd
      enddo

!      call rone(tmp4,ntot2)
!      call ortho(tmp4)

      call col2(dp,bm2,ntot2) ! Mass matrix

      call crs_solve_l2(tmp4,dp)

      call outpost(tmp1,tmp2,tmp3,tmp4,tmp5,'lap')

!      call ortho(dp)
      
      call rone(tmp8,ntot2)
      call cdabdtp(tmp4,tmp8,h1,h2,h2inv,intype)

      call opgradt (w1 ,w2 ,w3 ,tmp8)
      call opbinv  (tmp1,tmp2,tmp3,w1 ,w2 ,w3 ,h2inv)
      call opdiv   (tmp4,tmp1,tmp2,tmp3)

!      call opzero(tmp1,tmp2,tmp3)
      call copy(tmp8,dp,ntot2)

      call outpost(tmp1,tmp2,tmp3,tmp8,tmp5,'lap') 

!     Solve
      igmres = 3        ! 1: weighted; 4: Left preconditioned, 3: Std
      call esolver_new (dp,h1,h2,h2inv,intype,igmres)

      call opgradt (w1 ,w2 ,w3 ,dp)
      call opbinv  (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)

      call opcopy(tmp1,tmp2,tmp3,dv1,dv2,dv3)
      call copy(tmp4,dp,ntot2)
     
      call outpost(tmp1,tmp2,tmp3,tmp4,tmp5,'lap') 

!     Solve std.
      call copy(dp,tmp8,ntot2)            ! restore dp

      igmres = 3        ! standard 
      call esolver_new (dp,h1,h2,h2inv,intype,igmres)

      call opgradt (w1 ,w2 ,w3 ,dp)
      call opbinv  (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)

      call opcopy(tmp1,tmp2,tmp3,dv1,dv2,dv3)
     
      call outpost(tmp1,tmp2,tmp3,dp,tmp5,'lap')
     
!     difference between standard and new gmres      
      call sub2(tmp4,dp,ntot2) 
      call outpost(tmp1,tmp2,tmp3,tmp4,tmp5,'lap')


14    format(A5,2x,16(E12.5,2x))


      return
      end
c-----------------------------------------------------------------------
      subroutine esolver_new (res,h1,h2,h2inv,intype,ig)
C
C     Choose E-solver
C
C--------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'ESOLV'
      INCLUDE 'INPUT'
      include 'CTIMER'
C
      real res   (lx2,ly2,lz2,lelv)
      real h1    (lx1,ly1,lz1,lelv)
      real h2    (lx1,ly1,lz1,lelv)
      real h2inv (lx1,ly1,lz1,lelv)
      common /scruz/ wk1(lx2*ly2*lz2*lelv)
     $             , wk2(lx2*ly2*lz2*lelv)
     $             , wk3(lx2*ly2*lz2*lelv)

      integer ig

      integer igmres

      if (icalld.eq.0) teslv=0.0

!      call ortho_left(res) !Ensure that residual is orthogonal to null space

      icalld=icalld+1
      neslv=icalld
      etime1=dnekclock()

      if (.not. ifsplit) then
        if (param(42).eq.1) then
          call uzawa_new(res,h1,h2,h2inv,intype,icg)
        else
          if (ig.eq.1) call uzawa_gmres_wt(res,h1,h2,h2inv,intype,icg)
!          if (ig.eq.2) call uzawa_gmres_new(res,h1,h2,h2inv,intype,icg)
          if (ig.eq.3) call uzawa_gmres_std(res,h1,h2,h2inv,intype,icg)
          if (ig.eq.4) call uzawa_gmres_lpr(res,h1,h2,h2inv,intype,icg)
          if (ig.eq.5) call uzawa_gmres_cyl(res,h1,h2,h2inv,intype,icg)
          if (ig.eq.6) then
            igmres = 3
            call uzawa_gmres_testing(res,h1,h2,h2inv,intype,icg,igmres)
          endif  
          if (ig.gt.6) then 
            write(6,*) 'Unknown GMRES. exitting in esolver_new()'
            call exitt
          endif  
        endif
      else
        write(6,*) 'error: e-solver does not exist pnpn'
        call exitt
      endif

      teslv=teslv+(dnekclock()-etime1)

      RETURN
      END
c-----------------------------------------------------------------------
      subroutine uzawa_new(rcg,h1,h2,h2inv,intype,iter)
C-----------------------------------------------------------------------
C
C     Solve the pressure equation by (nested) preconditioned 
C     conjugate gradient iteration.
C     INTYPE =  0  (steady)
C     INTYPE =  1  (explicit)
C     INTYPE = -1  (implicit)
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      COMMON  /CTOLPR/ DIVEX
      COMMON  /CPRINT/ IFPRINT
      LOGICAL          IFPRINT
      REAL             RCG  (LX2,LY2,LZ2,LELV)
      REAL             H1   (LX1,LY1,LZ1,LELV)
      REAL             H2   (LX1,LY1,LZ1,LELV)
      REAL             H2INV(LX1,LY1,LZ1,LELV)
      COMMON /SCRUZ/   WP   (LX2,LY2,LZ2,LELV)
     $ ,               XCG  (LX2,LY2,LZ2,LELV)
     $ ,               PCG  (LX2,LY2,LZ2,LELV) 
     $ ,               RPCG (LX2,LY2,LZ2,LELV)
 
      real*8 etime1,dnekclock
      integer*8 ntotg,nxyz2


      etime1 = dnekclock()
      DIVEX = 0.
      ITER  = 0
c
      CALL CHKTCG2 (TOLPS,RCG,ICONV)
      if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
     $   TOLPS = abs(param(21))
C
c      IF (ICONV.EQ.1) THEN
c         IF (NID.EQ.0) WRITE(6,9999) ITER,DIVEX,TOLPS
c         return
c      ENDIF

      nxyz2 = lx2*ly2*lz2
      ntot2 = nxyz2*nelv
      ntotg = nxyz2*nelgv

!      CALL UZPREC  (RPCG,RCG,H1,H2,INTYPE,WP)
      call eprec2_new(rpcg,rcg)
     
      RRP1 = GLSC2 (RPCG,RCG,NTOT2)
      CALL COPY    (PCG,RPCG,NTOT2)
      CALL RZERO   (XCG,NTOT2)
      if (rrp1.eq.0) return
      BETA = 0.
      div0=0.
C
      tolpss = tolps
      DO 1000 ITER=1,10000 ! NMXP
C
C        CALL CONVPR  (RCG,tolpss,ICONV,RNORM)
         call convprn (iconv,rnorm,rrp1,rcg,rpcg,tolpss)

         if (iter.eq.1)      div0   = rnorm
         if (param(21).lt.0) tolpss = abs(param(21))*div0

         ratio = rnorm/div0
         IF (IFPRINT.AND.NIO.EQ.0) 
     $   WRITE (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66    format(i5,1p4e12.5,i8,' Divergence')
c
         IF (ICONV.EQ.1.and.iter.gt.1) GOTO 9000
c        IF (ICONV.EQ.1.and.(iter.gt.1.or.istep.le.2)) GOTO 9000
c        IF (ICONV.EQ.1) GOTO 9000
c        if (ratio.le.1.e-5) goto 9000


         IF (ITER .NE. 1) THEN
            BETA = RRP1/RRP2
            CALL ADD2S1 (PCG,RPCG,BETA,NTOT2)
         ENDIF

         CALL CDABDTP  (WP,PCG,H1,H2,H2INV,INTYPE)
!         CALL CM1DABDTP  (WP,PCG,H1,H2,H2INV,INTYPE)
!         CALL CDDTP    (WP,PCG)

         call ortho_left(wp)        ! prabal         

         PAP   = GLSC2 (PCG,WP,NTOT2)

         IF (PAP.NE.0.) THEN
            ALPHA = RRP1/PAP
         ELSE
            pcgmx = glamax(pcg,ntot2)
            wp_mx = glamax(wp ,ntot2)
            ntot1 = lx1*ly1*lz1*nelv
            h1_mx = glamax(h1 ,ntot1)
            h2_mx = glamax(h2 ,ntot1)
            if (nid.eq.0) write(6,*) 'ERROR: pap=0 in uzawa.'
     $      ,iter,pcgmx,wp_mx,h1_mx,h2_mx
            call exitt
         ENDIF
         CALL ADD2S2 (XCG,PCG,ALPHA,NTOT2)
         CALL ADD2S2 (RCG,WP,-ALPHA,NTOT2)

         if (iter.eq.-1) then
            call convprn (iconv,rnrm1,rrpx,rcg,rpcg,tolpss)
            if (iconv.eq.1) then
               rnorm = rnrm1
               ratio = rnrm1/div0
               if (nio.eq.0) 
     $         write (6,66) iter,tolpss,rnrm1,div0,ratio,istep
               goto 9000
            endif
         endif

!         call ortho(rcg)

         RRP2 = RRP1
!         CALL UZPREC  (RPCG,RCG,H1,H2,INTYPE,WP)
         call eprec2_new(rpcg,rcg)
        
c        RRP1 = GLSC2 (RPCG,RCG,NTOT2)

 1000 CONTINUE
      if (nid.eq.0) WRITE (6,3001) ITER,RNORM,tolpss
c     if (istep.gt.20) CALL EMERXIT
 3001 FORMAT(I6,' **ERROR**: Failed to converge in UZAWA:',6E13.4)
 9000 CONTINUE

      divex = rnorm
      iter  = iter-1

!     prabal
!      call uzprec(rcg,xcg,h1,h2,intype,wp)
!      call copy(xcg,rcg,ntot2)

      if (iter.gt.0) call copy (rcg,xcg,ntot2)
!      call ortho(rcg)
      call ortho_right(rcg)   ! prabal

      etime1 = dnekclock()-etime1
      IF (NIO.EQ.0) WRITE(6,9999) ISTEP, '  U-Press std. ',
     &                            ITER,DIVEX,div0,tolpss,etime1
 9999 FORMAT(I11,a,I7,1p4E13.4)
19999 FORMAT(I11,'  U-Press 1.e-5: ',I7,1p4E13.4)
C
C
      return
      END
c-----------------------------------------------------------------------

      subroutine uzawa_gmres_new(res,h1,h2,h2inv,intype,iter)

c     Solve the pressure equation by right-preconditioned 
c     GMRES iteration.
c     intype =  0  (steady)
c     intype =  1  (explicit)
c     intype = -1  (implicit)

      include 'SIZE'
      include 'TOTAL'
      include 'GMRES'
      common  /ctolpr/ divex
      common  /cprint/ ifprint
      logical          ifprint
      real             res  (lx2*ly2*lz2*lelv)
      real             h1   (lx1,ly1,lz1,lelv)
      real             h2   (lx1,ly1,lz1,lelv)
      real             h2inv(lx1,ly1,lz1,lelv)

      common /scrmg/    wp (lx2,ly2,lz2,lelv)

      common /ctmp0/   wk1(lgmres),wk2(lgmres)
      common /cgmres1/ y(lgmres)

      real alpha, l, temp
      integer j,m
c
      logical iflag
      save    iflag
      data    iflag /.false./
      real    norm_fac
      save    norm_fac
c
      real*8 etime1,dnekclock
c
      if(.not.iflag) then
         iflag=.true.
         call uzawa_gmres_split0(ml_gmres,mu_gmres,bm2,bm2inv,
     $                           lx2*ly2*lz2*nelv)
         norm_fac = 1./sqrt(volvm2)
      endif
c
      etime1 = dnekclock()
      etime_p = 0.
      divex = 0.
      iter  = 0
      m = lgmres
c
      call chktcg2(tolps,res,iconv)
      if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
     $   tolps = abs(param(21))
c     if (param(21).lt.0) tolps = abs(param(21))
      if (istep.eq.0) tolps = 1.e-4
      tolpss = tolps
c
      ntot2  = lx2*ly2*lz2*nelv
c
      iconv = 0
      call rzero(x_gmres,ntot2)

      do while(iconv.eq.0.and.iter.lt.20000)

         if(iter.eq.0) then
                                                        !      -1
            call col3(r_gmres,ml_gmres,res,ntot2)       ! r = L  res
c           call copy(r_gmres,res,ntot2)
         else
            !update residual
            call copy(r_gmres,res,ntot2)                      ! r = res

            call cdabdtp(w_gmres,x_gmres,h1,h2,h2inv,intype)  ! w = A x
!            call cM1dabdtp(w_gmres,x_gmres,h1,h2,h2inv,intype)  ! w = A x
!            call cddtp(w_gmres,x_gmres)                       ! w = A x

            call ortho_left(w_gmres)

            call add2s2(r_gmres,w_gmres,-1.,ntot2)            ! r = r - w
                                                              !      -1
            call col2(r_gmres,ml_gmres,ntot2)                 ! r = L   r
         endif
                                                            !            ______
         gamma_gmres(1) = sqrt(glsc2(r_gmres,r_gmres,ntot2))! gamma  = \/ (r,r) 
                                                            !      1
         if(iter.eq.0) then
            div0 = gamma_gmres(1)*norm_fac
            if (param(21).lt.0) tolpss=abs(param(21))*div0
         endif

         !check for lucky convergence
         rnorm = 0.
         if(gamma_gmres(1) .eq. 0.) goto 9000
         temp = 1./gamma_gmres(1)
         call cmult2(v_gmres(1,1),r_gmres,temp,ntot2)! v  = r / gamma
                                                     !  1            1
         do j=1,m
            iter = iter+1
                                                           !       -1
            call col3(w_gmres,mu_gmres,v_gmres(1,j),ntot2) ! w  = U   v
                                                           !           j
            
            etime2 = dnekclock()
            if(param(43).eq.1) then
!               call uzprec(z_gmres(1,j),w_gmres,h1,h2,intype,wp)
               call eprec2_new(z_gmres(1,j),w_gmres)
            else                                        !       -1
               call hsmg_solve(z_gmres(1,j),w_gmres)    ! z  = M   w
c              call copy(z_gmres(1,j),w_gmres,ntot2)    ! z  = M   w
            endif     
            etime_p = etime_p + dnekclock()-etime2
    
!            call ortho_right(z_gmres(1,j)) 

            call cdabdtp(w_gmres,z_gmres(1,j),    ! w = A z
     $                   h1,h2,h2inv,intype)      !        j

!            call cm1dabdtp(w_gmres,z_gmres(1,j),    ! w = A z
!     $                   h1,h2,h2inv,intype)      !        j


!            call cddtp(w_gmres,z_gmres(1,j))      ! w = A z
                                                  !        j

            call ortho_left(w_gmres)
                                                  !      -1
            call col2(w_gmres,ml_gmres,ntot2)     ! w = L   w

c           !modified Gram-Schmidt
c           do i=1,j
c              h_gmres(i,j)=glsc2(w_gmres,v_gmres(1,i),ntot2) ! h    = (w,v )
c                                                             !  i,j       i
c              call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),ntot2) ! w = w - h    v
c           enddo                                                    !          i,j  i


c           2-PASS GS, 1st pass:

            do i=1,j
               h_gmres(i,j)=vlsc2(w_gmres,v_gmres(1,i),ntot2) ! h    = (w,v )
            enddo                                             !  i,j       i

            call gop(h_gmres(1,j),wk1,'+  ',j)          ! sum over P procs

            do i=1,j
               call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),ntot2) ! w = w - h    v
            enddo                                                    !          i,j  i


c           2-PASS GS, 2nd pass:
c
c           do i=1,j
c              wk1(i)=vlsc2(w,v_gmres(1,i),ntot2) ! h    = (w,v )
c           enddo                                 !  i,j       i
c                                                 !
c           call gop(wk1,wk2,'+  ',j)             ! sum over P procs
c
c           do i=1,j
c              call add2s2(w,v_gmres(1,i),-wk1(i),ntot2) ! w = w - h    v
c              h(i,j) = h(i,j) + wk1(i)                  !          i,j  i
c           enddo


            !apply Givens rotations to new column
            do i=1,j-1
               temp = h_gmres(i,j)                   
               h_gmres(i  ,j)=  c_gmres(i)*temp 
     $                        + s_gmres(i)*h_gmres(i+1,j)  
               h_gmres(i+1,j)= -s_gmres(i)*temp 
     $                        + c_gmres(i)*h_gmres(i+1,j)
            enddo
                                                              !            ______
            alpha = sqrt(glsc2(w_gmres,w_gmres,ntot2))        ! alpha =  \/ (w,w)
            rnorm = 0.
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h_gmres(j,j)*h_gmres(j,j)+alpha*alpha)
            temp = 1./l
            c_gmres(j) = h_gmres(j,j) * temp
            s_gmres(j) = alpha  * temp
            h_gmres(j,j) = l
            gamma_gmres(j+1) = -s_gmres(j) * gamma_gmres(j)
            gamma_gmres(j)   =  c_gmres(j) * gamma_gmres(j)

c            call outmat(h,m,j,' h    ',j)
            
            rnorm = abs(gamma_gmres(j+1))*norm_fac
            ratio = rnorm/div0
            if (ifprint.and.nio.eq.0) 
     $         write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66       format(i5,1p4e12.5,i8,' Divergence')

#ifndef FIXITER
            if (rnorm .lt. tolpss) goto 900  !converged
#else
            if (iter.gt.param(151)-1) goto 900
#endif
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v_gmres(1,j+1),w_gmres,temp,ntot2) ! v    = w / alpha
                                                           !  j+1            
         enddo
  900    iconv = 1
 1000    continue
         !back substitution
         !     -1
         !c = H   gamma
         do k=j,1,-1
            temp = gamma_gmres(k)
            do i=j,k+1,-1
               temp = temp - h_gmres(k,i)*c_gmres(i)
            enddo
            c_gmres(k) = temp/h_gmres(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x_gmres,z_gmres(1,i),c_gmres(i),ntot2) 
                       ! x = x + c  z
                       !          i  i
         enddo
c        if(iconv.eq.1) call dbg_write(x,lx2,ly2,lz2,nelv,'esol',3)
         call ortho_right(x_gmres)
      enddo
 9000 continue
c
      divex = rnorm
c     iter = iter - 1

      call copy(res,x_gmres,ntot2)

!      call ortho (res)  ! Orthogonalize wrt null space, if present

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep,'  U-PRES gmres  ', 
     &                            iter,divex,div0,tolpss,etime_p,etime1
c     call flush_hack
 9999 format(i11,a,I6,1p5e13.4)

      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE EPREC2_NEW(Z2,R2)
C----------------------------------------------------------------
C
C     Precondition the explicit pressure operator (E) with
C     a Neumann type (H1) Laplace operator: JT*A*J.
C     Invert A by conjugate gradient iteration or multigrid.
C
C     NOTE: SCRNS is used.
C
C----------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'SOLN'
      INCLUDE 'MASS'
      INCLUDE 'PARALLEL'
      INCLUDE 'TSTEP'
      REAL           Z2   (LX2,LY2,LZ2,LELV)
      REAL           R2   (LX2,LY2,LZ2,LELV)
      COMMON /SCRNS/ MASK (LX1,LY1,LZ1,LELV)
     $              ,R1   (LX1,LY1,LZ1,LELV)
     $              ,X1   (LX1,LY1,LZ1,LELV)
     $              ,W2   (LX2,LY2,LZ2,LELV)
     $              ,H1   (LX1,LY1,LZ1,LELV)
     $              ,H2   (LX1,LY1,LZ1,LELV)
      REAL    MASK
c
      integer icalld
      save    icalld
      data    icalld/0/
      icalld=icalld+1
c
      ntot2  = lx2*ly2*lz2*nelv
      call rzero(z2,ntot2)

c     Both local and global solver...
      call dd_solver_new (z2,r2)


c
c  Local solver only
c      call local_solves_fdm (z2,r2)
c
c
c
      return
      end
c-----------------------------------------------------------------------
      subroutine dd_solver_new(u,v)

      implicit none

      include 'SIZE'
      include 'DOMAIN'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'CTIMER'
c
      real u(1),v(1)
      real uc
      common /scrprc/ uc(lx1*ly1*lz1*lelt)

      integer ntot
      real alpha

      if (icalld.eq.0) then
         tddsl=0.0
         tcrsl=0.0
         nddsl=0
         ncrsl=0
      endif
      icalld = icalld + 1
      nddsl  = nddsl  + 1
      ncrsl  = ncrsl  + 1

      ntot  = lx2*ly2*lz2*nelv

!      if (ifpgll) then
!        call copy(u,v,ntot)
!        return
!      endif
!
!      if (lx2.lt.lx1-2) then
!        call copy(u,v,ntot)
!        return
!      endif

      call copy(u,v,ntot)

      call rzero(u,ntot)
      etime1=dnekclock()
      call local_solves_fdm    (u,v)
      tddsl=tddsl+dnekclock()-etime1

      etime1=dnekclock()
      call crs_solve_l2(uc,v)
      tcrsl=tcrsl+dnekclock()-etime1

      alpha = 1.00
c     if (param(89).ne.0.) alpha = abs(param(89))
      call add2s2(u,uc,alpha,ntot)

      return
      end
c-----------------------------------------------------------------------

      subroutine cddtp (ap,wp)

C     INTYPE= 0  Compute the matrix-vector product    D*DT*p
C     INTYPE= 1  Compute the matrix-vector product    D*DT*p
C     INTYPE=-1  Compute the matrix-vector product    D*DT*p

      implicit none

      include 'SIZE'

      REAL           AP(LX2,LY2,LZ2,LELV)
      REAL           WP(LX2,LY2,LZ2,LELV)

      REAL TA1,TA2,TA3,TB1,TB2,TB3
      COMMON /SCRNS/ TA1 (LX1,LY1,LZ1,LELV)
     $ ,             TA2 (LX1,LY1,LZ1,LELV)
     $ ,             TA3 (LX1,LY1,LZ1,LELV)
     $ ,             TB1 (LX1,LY1,LZ1,LELV)
     $ ,             TB2 (LX1,LY1,LZ1,LELV)
     $ ,             TB3 (LX1,LY1,LZ1,LELV)

      call opgradt (ta1,ta2,ta3,wp)       ! DT*p
      call opdiv   (ap,ta1,ta2,ta3)       ! D*DT*p

      return
      end
C
C-----------------------------------------------------------------------

      subroutine map_f_to_c_l2_bilin_test(uc,uf,w)

c     TRANSPOSE of L2 Iterpolation operator:                    T
c                                 (linear --> spectral GLL mesh)

      include 'SIZE'
      include 'DOMAIN'
      include 'INPUT'

      parameter (lxyz = lx2*ly2*lz2)
      real uc(nxyz_c,lelt),uf(lxyz,lelt),w(1)

      ltot22 = 2*lx2*ly2*lz2
      nx_crs = 2   ! bilinear only

      do ie=1,nelv
         call maph1_to_l2t_test(uc(1,ie),nx_crs,uf(1,ie),
     $                          lx2,if3d,w,ltot22)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------

      subroutine maph1_to_l2t_test(b,nb,a,na,if3d,w,ldw)
c
c     Input:   a
c     Output:  b
c
      real a(1),b(1),w(1)
      logical if3d
c
      parameter(lx=50)
      real za(lx),zb(lx)
c
      real iba(lx*lx),ibat(lx*lx)
      save iba,ibat
c
      integer nao,nbo
      save    nao,nbo
      data    nao,nbo  / -9, -9/
c
c
      if (na.gt.lx.or.nb.gt.lx) then
         write(6,*)'ERROR: increase lx in maph1_to_l2 to max:',na,nb
         call exitt
      endif
c
      if (na.ne.nao  .or.   nb.ne.nbo) then
         nao = na
         nbo = nb
         call zwgl (za,w,na)
         call zwgll(zb,w,nb)
!         call igllm(iba,ibat,zb,za,nb,na,nb,na)
         call iglm(iba,ibat,za,zb,na,nb,na,nb)
      endif
c
!      call specmpn(b,nb,a,na,ibat,iba,if3d,w,ldw)
      call specmpn(b,nb,a,na,iba,ibat,if3d,w,ldw)
     
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine crs_solve_l2_test(uf,vf)
c
c     Given an input vector v, this generates the H1 coarse-grid solution
c
      include 'SIZE'
      include 'DOMAIN'
      include 'ESOLV'
      include 'GEOM'
      include 'PARALLEL'
      include 'SOLN'
      include 'INPUT'
      include 'TSTEP'
      real uf(1),vf(1)
      common /scrpre/ uc(lcr*lelt),w(2*lx1*ly1*lz1)

      call map_f_to_c_l2_bilin_test(uf,vf,w)
      call fgslib_crs_solve(xxth(ifield),uc,uf)
      call map_c_to_f_l2_bilin(uf,uc,w)

      return
      end
c
c-----------------------------------------------------------------------

      subroutine new_fdm(u,v)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'DOMAIN'
      include 'PARALLEL'
c
      include 'TSTEP'
      include 'CTIMER'
c
      real u(lx2,ly2,lz2,lelv),v(lx2,ly2,lz2,lelv)
      real v1,w1,w2
      common /scrpre/ v1(lx1,ly1,lz1,lelv)
     $               ,w1(lx1,ly1,lz1),w2(lx1,ly1,lz1)

      integer lxx,levb
      parameter(lxx=lx1*lx1, levb=lelv+lbelv)

      real df,sr,ss,st
      common /fastd/  df(lx1*ly1*lz1,levb)
     $             ,  sr(lxx*2,levb),ss(lxx*2,levb),st(lxx*2,levb)
      integer e,eb,eoff

      integer ntot1,ntot2

      if (icalld.eq.0) tsolv=0.0
      icalld=icalld+1
      nsolv=icalld
c
      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv

!     Put Mesh 2 points on Mesh 1
      call exchange_m2(v1,v)

c     Now solve each subdomain problem:

      etime1=dnekclock()

      eoff  = 0
      if (ifield.gt.1) eoff  = nelv

      do e = 1,nelv
         eb = e + eoff
         call fastdm1(v1(1,1,1,e),df(1,eb)
     $                           ,sr(1,eb),ss(1,eb),st(1,eb),w1,w2)
      enddo
      tsolv=tsolv+dnekclock()-etime1
c
c     Exchange/add elemental solutions
c
!      call s_face_to_int (v1,-1.)
!      call dssum         (v1,lx1,ly1,lz1)
!      call s_face_to_int (v1, 1.)
!      if(param(42).eq.0) call do_weight_op(v1)
c
c     Map back to pressure grid (extract interior values)
c
      call extract_interior(u,v1)

      return
      end subroutine new_fdm
!---------------------------------------------------------------------- 

      subroutine incomprn_test (ux,uy,uz,up)
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
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
c
      common /scrns/ w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)
     $ ,             dv1   (lx1,ly1,lz1,lelv)
     $ ,             dv2   (lx1,ly1,lz1,lelv)
     $ ,             dv3   (lx1,ly1,lz1,lelv)
     $ ,             dp    (lx2,ly2,lz2,lelv)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)

      parameter(nset = 1 + lbelv/lelv)
      common /orthov/ pset(lx2*ly2*lz2*lelv*mxprev,nset)
      common /orthbi/ nprv(2)
      logical ifprjp

      ifprjp=.false.    ! Project out previous pressure solutions?
      istart=param(95)  
      if (istep.ge.istart.and.istart.ne.0) ifprjp=.true.

      if (icalld.eq.0) tpres=0.0
      icalld = icalld+1
      npres  = icalld
      etime1 = dnekclock()

      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv

      intype = -1

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

      call opdiv   (dp,ux,uy,uz)

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
!                    call esolver  (dp,h1,h2,h2inv,intype)  ! prabal
                    call esolver_new(dp,h1,h2,h2inv,intype,3)
                    call cmult(dp,scaledi,ntot2)
                 else
!                    call esolver  (dp,h1,h2,h2inv,intype)  ! prabal
                    call esolver_new(dp,h1,h2,h2inv,intype,1)
                 endif 
      if (ifprjp)   call gensolnp (dp,h1,h2,h2inv,pset(1,i),nprv(i))

      call add2(up,dp,ntot2)

      call opgradt  (w1 ,w2 ,w3 ,dp)
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

      if (ifmhd)  call chkptol	! to avoid repetition

      tpres=tpres+(dnekclock()-etime1)

      return
      end
c-----------------------------------------------------------------------


