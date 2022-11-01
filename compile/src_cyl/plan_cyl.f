!======================================================================
!
!     Author: Prabal Negi
!     Description: Cylindrical coordinates solver
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

         intype = -1
         call sethlm  (h1,h2,intype)
         call cresvif_cyl (resv1,resv2,resv3,h1,h2)

         call ophinv (dv1,dv2,dv3,resv1,resv2,resv3,h1,h2,tolhv,nmxv)
         call opadd2 (vx,vy,vz,dv1,dv2,dv3)

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
        if (ifuservp) then
          call advab_rho_cyl()
        else
          call advab_cyl()
        endif
      endif
      if (ifmvbd.and..not.ifchar)               call admeshv
      if (iftran)                               call makeabf
      if ((iftran.and..not.ifchar).or.
     $    (iftran.and..not.ifnav.and.ifchar))   call makebdf
      if (ifnav.and.ifchar)                     call advchar
      
      tmakf=tmakf+(dnekclock()-etime1)

      return
      end
!---------------------------------------------------------------------- 
      subroutine advab_rho_cyl()

!     Eulerian scheme, add convection term to forcing function 
!     at current time step.
!     with variable density        

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      real rhoi         ! density inverse
      common /scrcg/ rhoi(lx1,ly1,lz1,lelt)

      integer ntot1

      ntot1 = lx1*ly1*lz1*nelv

!     bfx,bfy,bfz get multiplied by vtrans later in makeabf.
!     So I do an inverse rho multiplication here

!     rhoi = 1/\rho      
      call invers2(rhoi,vtrans(1,1,1,1,ifield),ntot1)

      call convect_cylindrical_rho (ta1,vtrans(1,1,1,1,ifield),vx,
     $                              vx,vy,vz)
      call convect_cylindrical_rho (ta2,vtrans(1,1,1,1,ifield),vy,
     $                              vx,vy,vz)
      if (ldim.eq.3) then
        call convect_cylindrical_rho (ta3,vtrans(1,1,1,1,ifield),vz,
     $                                vx,vy,vz)
      endif  

      call col2(ta1,rhoi,ntot1)
      call sub2 (bfx,ta1,ntot1)

      call col2(ta2,rhoi,ntot1)      
      call sub2 (bfy,ta2,ntot1)
      if (ldim.eq.3) then
        call col2(ta3,rhoi,ntot1)
        call sub2 (bfz,ta3,ntot1)
      endif

!     Additional terms in cylindrical formulation
!     Division by R gets cancelled by the multiplication by R
!     from the Jacobian
      if (ldim.eq.3) then      
        call dealias_rho_uv(ta2,vtrans(1,1,1,1,ifield),vz,vz)
        call col2(ta2,rhoi,ntot1)      
        call add2(bfy,ta2,ntot1)
        
        call dealias_rho_uv(ta3,vtrans(1,1,1,1,ifield),vy,vz)
        call col2(ta3,rhoi,ntot1)      
        call sub2(bfz,ta3,ntot1)
      endif  

      return
      end subroutine
!---------------------------------------------------------------------- 

      subroutine advab_cyl()

!     Eulerian scheme, add convection term to forcing function 
!     at current time step.

      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'TSTEP'

      real ta1,ta2,ta3
      common /scruz/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)

      integer ntot1

      ntot1 = lx1*ly1*lz1*nelv

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

!     Additional terms in cylindrical formulation
!     Division by R gets cancelled by the multiplication by R
!     from the Jacobian
      if (ldim.eq.3) then      
        call dealias_uv(ta2,vz,vz)
        call add2(bfy,ta2,ntot1)
        
        call dealias_uv(ta3,vy,vz)
        call sub2(bfz,ta3,ntot1)
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

      integer gmtype

      gmtype = 5        ! std cylindrical

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
!         call esolver  (dp,h1,h2,h2inv,intype)  ! prabal
         call esolver_new(dp,h1,h2,h2inv,intype,gmtype)
         call cmult(dp,scaledi,ntot2)
      else
!         call esolver  (dp,h1,h2,h2inv,intype)  ! prabal
         call esolver_new(dp,h1,h2,h2inv,intype,gmtype)
      endif 
      if (ifprjp)   call gensolnp (dp,h1,h2,h2inv,pset(1,i),nprv(i))

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

      if (ifmhd)  call chkptol	! to avoid repetition

      tpres=tpres+(dnekclock()-etime1)

      return
      end
c-----------------------------------------------------------------------
      subroutine esolver_cyl (res,h1,h2,h2inv,intype)

!     Choose E-solver

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
          call uzawa_cyl(res,h1,h2,h2inv,intype,icg)
        else
          call uzawa_gmres_cyl(res,h1,h2,h2inv,intype,icg)
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

      subroutine uzawa_gmres_cyl(res,h1,h2,h2inv,intype,iter)

c     Solve the cylindrical pressure equation by right-preconditioned 
c     GMRES iteration.
c     intype =  0  (steady)
c     intype =  1  (explicit)
c     intype = -1  (implicit)

      implicit none

      include 'SIZE'
      include 'INPUT'   ! param
      include 'MASS'    ! bm2
      include 'TSTEP'
      include 'GMRES'

      integer lt1,lt2
      parameter(lt1 = lx1*ly1*lz1*lelv)
      parameter(lt2 = lx2*ly2*lz2*lelv)

      real divex
      common  /ctolpr/ divex
      
      logical          ifprint
      common  /cprint/ ifprint

      real             res  (lt2)
      real             h1   (lt1)
      real             h2   (lt1)
      real             h2inv(lt1)
      
      real wp
      common /scrmg/   wp (lt2)

      real wk1,wk2
      common /ctmp0/   wk1(lgmres),wk2(lgmres)

      real y
      common /cgmres1/ y(lgmres)

      real alpha, l, temp
      integer j,m

      logical iflag
      save    iflag
      data    iflag /.false./
      real    norm_fac
      save    norm_fac

      real*8 etime1,dnekclock

      integer ntot2
      real glsc2,glsc3,vlsc2,vlsc3
      integer iconv,intype
      real tolpss,div0
      integer i,k,iter
      real etime2,etime_p,ratio,rnorm

      integer maxiter

      logical ifwgt           ! If Weighted orthogonalization
      integer ngs             ! No of Gram-Schmid

      ifwgt = .false.
      ngs   = 1

!     I've removed ortho from earlier calls and call it here at
!     the beginning      
      call ortho(res)
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
!      if (istep.eq.0) tolps = 1.e-4
      tolpss = tolps

      ntot2  = lx2*ly2*lz2*nelv

      iconv = 0
      call rzero(x_gmres,ntot2)

      do while(iconv.eq.0.and.iter.lt.1000)

         if(iter.eq.0) then
                                                        !      -1
            call col3(r_gmres,ml_gmres,res,ntot2)       ! r = L  res
         else
            !update residual
            call copy(r_gmres,res,ntot2)                      ! r = res
            call cdabdtp_cyl(w_gmres,x_gmres,h1,h2,h2inv,intype)  ! w = A x
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
               call uzprec(z_gmres(1,j),w_gmres,h1,h2,intype,wp)
            else                                        !       -1
               call hsmg_solve(z_gmres(1,j),w_gmres)    ! z  = M   w
            endif     
            etime_p = etime_p + dnekclock()-etime2
     
            call cdabdtp_cyl(w_gmres,z_gmres(1,j),    ! w = A z
     $                   h1,h2,h2inv,intype)      !        j
     
                                                  !      -1
            call col2(w_gmres,ml_gmres,ntot2)     ! w = L   w

!           Gram-Schmidt:
            call ortho_subspace(w_gmres,ntot2,h_gmres(1,j),v_gmres,
     $            lt2,j,bm2,ifwgt,ngs,wk1,wk2)

!           Apply Givens rotations to new column
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
         enddo          ! j=1,m

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

      enddo           ! while(iconv.eq.0.and.iter.lt.1000)
 9000 continue
c
      divex = rnorm

      call copy(res,x_gmres,ntot2)

      call ortho (res)  ! Orthogonalize wrt null space, if present

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep,'  U-PRES gmres (CYL) ', 
     &                            iter,divex,div0,tolpss,etime_p,etime1
c     call flush_hack
 9999 format(i11,a,I6,1p5e13.4)

      return
      end

c-----------------------------------------------------------------------
      subroutine uzawa_cyl (rcg,h1,h2,h2inv,intype,iter)

C     Solve the pressure equation by (nested) preconditioned 
C     conjugate gradient iteration.
C     INTYPE =  0  (steady)
C     INTYPE =  1  (explicit)
C     INTYPE = -1  (implicit)
      
      implicit none

      include 'SIZE'
      include 'MASS'
      include 'INPUT'
      include 'PARALLEL'
      include 'TSTEP' 

      integer lt1,lt2
      parameter(lt1 = lx1*ly1*lz1*lelv)
      parameter(lt2 = lx2*ly2*lz2*lelv)

      real divex
      common  /ctolpr/ divex

      logical          ifprint
      common  /cprint/ ifprint
      real             rcg  (lt2)
      real             h1   (lt1)
      real             h2   (lt1)
      real             h2inv(lt1)

      real wp,xcg,pcg,rpcg
      common /scruz/   wp   (lt2)
     $ ,               xcg  (lt2)
     $ ,               pcg  (lt2) 
     $ ,               rpcg (lt2)
 
      real*8 etime1,dnekclock
      integer*8 ntotg,nxyz2

      integer intype,iter,iconv,ntot1,ntot2
      real alpha,beta,rrp1,rrp2,pap,div0,ratio
      real rnrm1,rrpx,rnorm,tolpss
      real h1_mx,h2_mx,wp_mx

      real glamax,glsc2
      real pcgmx


      etime1 = dnekclock()
      DIVEX = 0.
      ITER  = 0

      call chktcg2 (tolps,rcg,iconv)
      if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
     $   tolps = abs(param(21))

      nxyz2 = lx2*ly2*lz2
      ntot2 = nxyz2*nelv
      ntotg = nxyz2*nelgv

      call uzprec  (rpcg,rcg,h1,h2,intype,wp)
      rrp1 = glsc2 (rpcg,rcg,ntot2)
      call copy    (pcg,rpcg,ntot2)
      call rzero   (xcg,ntot2)
      if (rrp1.eq.0) return
      beta = 0.
      div0=0.

      tolpss = tolps
      do 1000 iter=1,nmxp

         call convprn (iconv,rnorm,rrp1,rcg,rpcg,tolpss)

         if (iter.eq.1)      div0   = rnorm
         if (param(21).lt.0) tolpss = abs(param(21))*div0

         ratio = rnorm/div0
         if (ifprint.and.nio.eq.0) 
     $   write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66    format(i5,1p4e12.5,i8,' Divergence')
c
         if (iconv.eq.1.and.iter.gt.1) goto 9000

         if (iter .ne. 1) then
            beta = rrp1/rrp2
            call add2s1 (pcg,rpcg,beta,ntot2)
         endif

         call cdabdtp_cyl(wp,pcg,h1,h2,h2inv,intype)
         pap   = glsc2 (pcg,wp,ntot2)

         if (pap.ne.0.) then
            alpha = rrp1/pap
         else
            pcgmx = glamax(pcg,ntot2)
            wp_mx = glamax(wp ,ntot2)
            ntot1 = lx1*ly1*lz1*nelv
            h1_mx = glamax(h1 ,ntot1)
            h2_mx = glamax(h2 ,ntot1)
            if (nid.eq.0) write(6,*) 'ERROR: pap=0 in uzawa.'
     $      ,iter,pcgmx,wp_mx,h1_mx,h2_mx
            call exitt
         endif
         call add2s2 (xcg,pcg,alpha,ntot2)
         call add2s2 (rcg,wp,-alpha,ntot2)

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

         rrp2 = rrp1
         call uzprec  (rpcg,rcg,h1,h2,intype,wp)

 1000 continue
      if (nid.eq.0) write (6,3001) iter,rnorm,tolpss
 3001 format(I6,' **ERROR**: Failed to converge in UZAWA:',6E13.4)
 9000 continue

      divex = rnorm
      iter  = iter-1

      if (iter.gt.0) call copy (rcg,xcg,ntot2)
      call ortho(rcg)

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep, '  U-Press std (CYL). ',
     &                            iter,divex,div0,tolpss,etime1
 9999 format(I11,a,I7,1p4E13.4)
19999 format(I11,'  U-Press 1.e-5: ',I7,1p4E13.4)


      return
      end subroutine uzawa_cyl
c-----------------------------------------------------------------------

