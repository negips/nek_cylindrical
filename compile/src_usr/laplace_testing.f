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

      if (nio.eq.0) write(6,*) 'Pseudo Laplacian Arnoldi'

      ifield = 1
      istep  = 1

      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv

      call rone    (vtrans,ntot1)

      intype = 1        ! explicit

      call rzero   (h1,ntot1)
      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call invers2 (h2inv,h2,ntot1)

      call reset_preconditioner() 

      do i=1,ntot2
        call random_number(rnd)
        prp(i,1) = 1.0*rnd
      enddo

      call opzero(vxp,vyp,vzp)
      call cdabdtp(tmp4,prp,h1,h2,h2inv,intype)

      ifflow = .false.
      ifheat = .false.
!      call tst_init()   ! also calls arn_init()

      istep = 0
      do while (istep.lt.nsteps)

        istep = istep+1
        call settime
        call setdt
      
        if (nio.eq.0) write(6,*) 'Iteration:', istep,dt 

!        call rk4_advance(prp,dt)

!!       tmp = (U^T)*E*(U)*p
!        call ortho_right(prp) 
!        call cdabdtp(tmp4,prp,h1,h2,h2inv,intype)
!        call ortho_left(tmp4)
!!       p = (I + \lambda*E)*p
!        lambda = -1.0*dt
!        call add2s2(prp,tmp4,lambda,ntot2)

!       tmp = (U^T)*E*(U)*(M^-1)*p        
!        if (tst_vstep.eq.0.and.tst_istep.eq.0) call ortho_new(prp)
        call eprec2_new(tmp4,prp)
        call ortho_new(tmp4)
        call copy(tmp8,tmp4,ntot2)
        call cdabdtp(tmp4,tmp8,h1,h2,h2inv,intype)
        call ortho_new(tmp4)
!       p = (I + \lambda*E)*p
        lambda = -1.0*dt
        call add2s2(prp,tmp4,lambda,ntot2)


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
      call add3s2(p1,p,Ep,s0,s1,n2)            ! p1  = p - (dt/2)*E*p

      call cdabdtp(Ep1,p1,h1,h2,h2inv,intype)   ! Ep1 = E*p1
      call invcol2(Ep1,bm2,n2)                  ! Ep1 = (B^-1)*E*p1
      call add3s2(p2,p,Ep1,s0,s2,n2)           ! p2  = p - (dt/2)*E*p1

      call cdabdtp(Ep2,p2,h1,h2,h2inv,intype)   ! Ep2 = E*p2
      call invcol2(Ep2,bm2,n2)                  ! Ep2 = (B^-1)*E*p2
      call add3s2(p3,p,Ep2,s0,s3,n2)           ! p3  = p - (dt)*E*p2

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
      include 'MASS'    ! bm2

      integer lt,lt2
      parameter (lt  = lx1*ly1*lz1*lelt)
      parameter (lt2 = lx2*ly2*lz2*lelt)

      real vi
      real Ev
      common /scrns/ vi (lt,4)
     $ ,             Ev (lt,4)

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

      intype = 1

      n1  = lx1*ly1*lz1*nelv
      n2  = lx2*ly2*lz2*nelv

      n   = n2

      call rzero   (h1,n1)
      call rone    (h2,n1)
      call invers2 (h2inv,h2,n1)

      visc = -1.0

      s0    = 1.0

      si(1) = 2.0 
      si(2) = 2.0
      si(3) = 1.0
      si(4) = 6.0


      call rzero(Ev,lt*4)
      call rzero(vi,lt*4)

      call copy(vi(1,1),v,n2)
      call copy(wk1,v,n2)
!      call col2(wk1,bm2,n2)         ! Mass Matrix
      
      do i=1,3

!       M*vi
        call cdabdtp(Ev(1,i),vi(1,i),h1,h2,h2inv,intype)    ! Ev_i = E*v_i

!       vi+1 = v + (dt*fac)*M*vi
        s = visc*dt/si(i)
        call add3s2(vi(1,i+1),wk1,Ev(1,i),s0,s,n)  ! v_i+1 = Bv + (dt/2)*E*v_i
!        call invcol2(vi(1,i+1),bm2,n)              ! v_i+1 = (B^-1)*v_i+1
      enddo

      call cdabdtp(Ev(1,4),vi(1,4),h1,h2,h2inv,intype)   ! Ev4 = E*p4
!      call invcol2(Ev(1,4),bm2,n)                        ! Ev4 = (B^-1)*E*p4

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
!      call reset_preconditioner()

      intype = 1        ! explicit

      call rzero   (h1,ntot1)
!      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call rone    (h2,ntot1)
      call invers2 (h2inv,h2,ntot1)

      do i=1,ntot2
        call random_number(rnd)
        if (ifcyclic) then
          rad = sqrt(ym2(i,1,1,1)**2 + zm2(i,1,1,1)**2)
        else
          rad = ym2(i,1,1,1)
        endif
        dp(i,1,1,1) = (1.0e-0)*rad + 0.0*rnd
!        dp(i,1,1,1) = (1.0e-0)*xm2(i,1,1,1) + 0.0*rnd
      enddo

!      call rone(tmp4,ntot2)
!      call ortho(tmp4)

      call col2(dp,bm2,ntot2) ! Mass matrix

!      call crs_solve_l2(tmp4,dp)

      call outpost(tmp1,tmp2,tmp3,dp,tmp5,'lap')

!      call ortho(dp)
      
      call rone(tmp8,ntot2)
      call cdabdtp(tmp4,tmp8,h1,h2,h2inv,intype)

      call opgradt (w1 ,w2 ,w3 ,tmp8)
      call opbinv  (tmp1,tmp2,tmp3,w1 ,w2 ,w3 ,h2inv)
      call opdiv   (tmp4,tmp1,tmp2,tmp3)

!      call opzero(tmp1,tmp2,tmp3)
      call copy(tmp8,dp,ntot2)

!     Solve
      igmres = 1        ! 1: weighted; 4: Left preconditioned
      call esolver_new (dp,h1,h2,h2inv,intype,igmres)

      call opgradt (w1 ,w2 ,w3 ,dp)
      call opbinv  (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)

      call opcopy(tmp1,tmp2,tmp3,dv1,dv2,dv3)
      call copy(tmp4,dp,ntot2)
     
      call outpost(tmp1,tmp2,tmp3,tmp4,tmp5,'lap') 

!     Solve std.
      call copy(dp,tmp8,ntot2)            ! restore dp
      istep  = 3

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
          if (ig.eq.5) then
            igmres = 2
            call uzawa_gmres_testing(res,h1,h2,h2inv,intype,icg,igmres)
          endif  
          if (ig.gt.5) then 
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

      subroutine eprec2_new(z2,r2)
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
      if (.not.ifpgll) then
        call copy(u,v,ntot)
        return
      endif

      if (lx2.ne.(lx1-2)) then
        call copy(u,v,ntot)
        return
      endif

      call copy(u,v,ntot)

      call rzero(u,ntot)
      etime1=dnekclock()
      call local_solves_fdm    (u,v)
      tddsl=tddsl+dnekclock()-etime1

      etime1=dnekclock()
      call crs_solve_l2(uc,v)
      tcrsl=tcrsl+dnekclock()-etime1

      alpha = 10.00
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
      subroutine cM1dabdtp (ap,wp,h1,h2,h2inv,intype)

C     INTYPE= 0  Compute the matrix-vector product    DA(-1)DT*p
C     INTYPE= 1  Compute the matrix-vector product    D(B/DT)(-1)DT*p
C     INTYPE=-1  Compute the matrix-vector product    D(A+B/DT)(-1)DT*p

      include 'SIZE'
      include 'TOTAL'
      REAL           AP    (LX2,LY2,LZ2,1)
      REAL           WP    (LX2,LY2,LZ2,1)
      REAL           H1    (LX1,LY1,LZ1,1)
      REAL           H2    (LX1,LY1,LZ1,1)
      REAL           H2INV (LX1,LY1,LZ1,1)
C
      COMMON /SCRNS/ TA1 (LX1,LY1,LZ1,LELV)
     $ ,             TA2 (LX1,LY1,LZ1,LELV)
     $ ,             TA3 (LX1,LY1,LZ1,LELV)
     $ ,             TB1 (LX1,LY1,LZ1,LELV)
     $ ,             TB2 (LX1,LY1,LZ1,LELV)
     $ ,             TB3 (LX1,LY1,LZ1,LELV)

      call opgradtM1(ta1,ta2,ta3,wp)
      if ((intype.eq.0).or.(intype.eq.-1)) then
         tolhin=tolhs
         call ophinv (tb1,tb2,tb3,ta1,ta2,ta3,h1,h2,tolhin,nmxv)
      else
         if (ifanls) then
            dtbdi = dt/bd(1)   ! scale by dt*backwd-diff coefficient
            CALL OPBINV1(TB1,TB2,TB3,TA1,TA2,TA3,dtbdi)
         else
            CALL OPBINV (TB1,TB2,TB3,TA1,TA2,TA3,H2INV)
         endif
      endif
      call opdivM1 (ap,tb1,tb2,tb3)

      return
      end
C
C-----------------------------------------------------------------------

      subroutine opdivM1(outfld,inpx,inpy,inpz)
C---------------------------------------------------------------------
C
C     Compute OUTFLD = SUMi Di*INPi, 
C     the divergence of the vector field (INPX,INPY,INPZ)
C
C---------------------------------------------------------------------
      include 'SIZE'
      include 'GEOM'
      real outfld (lx2,ly2,lz2,1)
      real inpx   (lx1,ly1,lz1,1)
      real inpy   (lx1,ly1,lz1,1)
      real inpz   (lx1,ly1,lz1,1)
      common /ctmp0/ work (lx2,ly2,lz2,lelv)
C
      iflg = 1

      ntot2 = lx2*ly2*lz2*nelv
      call multdM1 (work,inpx,rxm1,sxm1,txm1,1,iflg)
      call copy  (outfld,work,ntot2)
      call multdM1 (work,inpy,rym1,sym1,tym1,2,iflg)
      call add2  (outfld,work,ntot2)
      if (ldim.eq.3) then
         call multdM1 (work,inpz,rzm1,szm1,tzm1,3,iflg)
         call add2  (outfld,work,ntot2)
      endif
C
      return
      end
C
c-----------------------------------------------------------------------
      subroutine opgradtM1(outx,outy,outz,inpfld)
C------------------------------------------------------------------------
C
C     Compute DTx, DTy, DTz of an input field INPFLD 
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      real outx   (lx1,ly1,lz1,1)
      real outy   (lx1,ly1,lz1,1)
      real outz   (lx1,ly1,lz1,1)
      real inpfld (lx2,ly2,lz2,1)
C
      call cdtM1p (outx,inpfld,rxm1,sxm1,txm1,1)
      call cdtM1p (outy,inpfld,rym1,sym1,tym1,2)
      if (ldim.eq.3) 
     $   call cdtM1p (outz,inpfld,rzm1,szm1,tzm1,3)
C
      return
      end
c-----------------------------------------------------------------------

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

      integer e,i1,i2,iz
      integer nxy1,nyz1,nxy2,nxyz1,nxyz2,n1,n2

C
#ifdef TIMER
      if (icalld.eq.0) tmltd=0.0
      icalld=icalld+1
      nmltd=icalld
      etime1=dnekclock()
#endif

      nxy1  = lx1*ly1
      nyz1  = ly1*lz1
      nxy2  = lx2*ly2
      nxyz1 = lx1*ly1*lz1
      nxyz2 = lx2*ly2*lz2

      n1    = lx1*ly1
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
     $         'MULTDM1 not implemented for 3D'
             call exitt

             call mxm  (dxm1,lx2,x(1,e),lx1,ta1,nyz1)
             call col3 (wk1,ta1,rm1(1,e),nxyz1)
!
             call copy(ta3,x(1,e),nxyz1)
             i1=1
             i2=1
             do iz=1,lz1
               call mxm (ta3(i1),lx1,dytm1,ly1,ta2(i2),ly1)
               i1=i1+n1
               i2=i2+n1
             enddo
             call addcol3 (wk1,ta2,sm1(1,e),nxyz1)
!
             call copy(ta1,x(1,e),nxyz1)
             call mxm (ta1,nxy1,dztm1,lz1,ta3,lz1)
             call addcol3 (wk1,ta3,tm1(1,e),nxyz1)
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
         if (if3d) then
           call mxm(ixtm21,lx2,wk1,lx1,ta1,nyz1)
           i1=1
           i2=1
           do iz=1,lz1
             call mxm (ta1(i1),lx2,iym21,ly1,ta2(i2),ly2)
             i1=i1+(lx2*ly1)
             i2=i2+n2
           enddo
           call mxm (ta2,nxy2,izm21,lz1,dx(1,e),lz2)
         else  
           call mxm(ixtm21,lx2,wk1,lx1,ta1,lx1)
           call mxm(ta1,lx2,iym21,lx1,dx(1,e),lx2)
         endif  


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

