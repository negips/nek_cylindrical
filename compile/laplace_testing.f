!======================================================================
!     Author: Prabal Negi
!     Description: Testing Solver/Preconditioner for Laplace.
!
!======================================================================       
      subroutine laplace_test()

      implicit none
  
      include 'SIZE'
      include 'INPUT'         ! uparam/param
      include 'SOLN'          ! vtrans
      include 'MASS'          ! bm2
      include 'TSTEP'         ! ifield
      include 'GEOM'          ! xm2
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

      if (nio.eq.0) write(6,*) 'Testing Laplace'

      ifield = 1

      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv

!     Preconditioner
      param(42)=uparam(8)       ! 0: GMRES (nonsymmetric), 1: PCG w/o weights
      param(43)=uparam(9)       ! 0: Additive multilevel (param 42=0), 1: Original 2 level
      param(44)=uparam(10)      ! 0: E based Schwartz (FEM), 1: A based Schwartz
      
      call rone(vtrans,ntot1) 
      call reset_preconditioner()

      intype = 1        ! explicit

      call rzero   (h1,ntot1)
      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call rone(h2,ntot1)
      call invers2 (h2inv,h2,ntot1)

      do i=1,ntot2
        dp(i,1,1,1) = (1.0e-0)*xm2(i,1,1,1)
      enddo

      call col2(dp,bm2,ntot2) ! Mass matrix

      call opzero(tmp1,tmp2,tmp3)
      call copy(tmp8,dp,ntot2)
      call outpost(tmp1,tmp2,tmp3,tmp8,tmp5,'lap') 
   
!      call col2(dp,bm2,ntot2) ! Mass matrix

      call esolver_new (dp,h1,h2,h2inv,intype)

      call opgradt  (w1 ,w2 ,w3 ,dp)
      call opbinv   (dv1,dv2,dv3,w1 ,w2 ,w3 ,h2inv)

      call opcopy(tmp1,tmp2,tmp3,dv1,dv2,dv3)
      call copy(tmp4,dp,ntot2)
     
      call outpost(tmp1,tmp2,tmp3,tmp4,tmp5,'lap') 

      call cdabdtp(tmp4,tmp8,h1,h2,h2inv,intype)

      call cdtp (tmp1,tmp8,rxm2,sxm2,txm2,1)
      iflg = 1
      call multd (tmp12,tmp1,rxm2,sxm2,txm2,1,iflg)

      call outpost(tmp1,tmp2,tmp3,tmp4,tmp5,'lap') 

      write(6,*) 'CDDTP'
      do i=1,lx2
        write(6,14) 'cddtp',(tmp12(i,j,1,1),j=1,lx2)
      enddo  

14    format(A5,2x,16(E12.5,2x))


      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE ESOLVER_NEW (RES,H1,H2,H2INV,INTYPE)
C
C     Choose E-solver
C
C--------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'ESOLV'
      INCLUDE 'INPUT'
C
      REAL RES   (LX2,LY2,LZ2,LELV)
      REAL H1    (LX1,LY1,LZ1,LELV)
      REAL H2    (LX1,LY1,LZ1,LELV)
      REAL H2INV (LX1,LY1,LZ1,LELV)
      common /scruz/ wk1(lx2*ly2*lz2*lelv)
     $             , wk2(lx2*ly2*lz2*lelv)
     $             , wk3(lx2*ly2*lz2*lelv)

      include 'CTIMER'
      real kwave2

      if (icalld.eq.0) teslv=0.0

      call ortho(res) !Ensure that residual is orthogonal to null space

      icalld=icalld+1
      neslv=icalld
      etime1=dnekclock()

      if (.not. ifsplit) then
         if (param(42).eq.1) then
            CALL UZAWA_NEW(RES,H1,H2,H2INV,INTYPE,ICG)
         else
            call uzawa_gmres_new(res,h1,h2,h2inv,intype,icg)
         endif
      else
         WRITE(6,*) 'ERROR: E-solver does not exist PnPn'
         CALL EXITT
      ENDIF

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
      DO 1000 ITER=1,500 ! NMXP
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
!         CALL CDDTP    (WP,PCG)
                                        

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

         call ortho(rcg)

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
      call ortho(rcg)

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

      do while(iconv.eq.0.and.iter.lt.1000)

         if(iter.eq.0) then
                                                        !      -1
            call col3(r_gmres,ml_gmres,res,ntot2)       ! r = L  res
c           call copy(r_gmres,res,ntot2)
         else
            !update residual
            call copy(r_gmres,res,ntot2)                      ! r = res
            call cdabdtp(w_gmres,x_gmres,h1,h2,h2inv,intype)  ! w = A x
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
     
            call cdabdtp(w_gmres,z_gmres(1,j),    ! w = A z
     $                   h1,h2,h2inv,intype)      !        j

!            call cddtp(w_gmres,z_gmres(1,j))      ! w = A z
                                                  !        j


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
      enddo
 9000 continue
c
      divex = rnorm
c     iter = iter - 1

      call copy(res,x_gmres,ntot2)

      call ortho (res)  ! Orthogonalize wrt null space, if present

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

c  Both local and global solver...
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
c
      include 'SIZE'
      include 'DOMAIN'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'CTIMER'
c
      real u(1),v(1)
      common /scrprc/ uc(lx1*ly1*lz1*lelt)
c
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
      call copy(u,v,ntot)

!      call rzero(u,ntot)
      etime1=dnekclock()
!      call local_solves_fdm    (u,v)
      tddsl=tddsl+dnekclock()-etime1

      etime1=dnekclock()
      call crs_solve_l2 (uc,v)
      tcrsl=tcrsl+dnekclock()-etime1

      alpha = 0.
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







