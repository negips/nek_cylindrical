!======================================================================
!     Author: Prabal Negi      
!     Description: Routines for 3D cylindrical solve implementation
!     Routines:   cdabdtp_cyl()           : Pressure Pseudolaplacian
!                 opdiv_cyl()             : Cylindrical Divergence
!                 multd_cyl()             : D*u = q*(D*u)
!                 opgradt_cyl()           : Pressure gradient term
!                 cdtp_cyl()              : (D^T)*p = p*(D*v)
!                 convect_cylindircal_rho : Cylindrical convective term with density
!                 convect_cylindircal     : Cylindrical convective term
!                 dealias_rho_uv          : Dealiased \rho*u*v      
!                 dealias_uv              : Dealiased u*v
!                 axhmsf_cyl              : Ax (cylindrical)
!                 stnrate_cyl             : 1/2*(grad(u) + grad(u)^T)      
!                 stress_cyl              : Sij = 2*mu*Eij
!                 div_stress_cyl          : Grad(v)\dot Sij
!                 ttxyz_cyl               : Grad(v)\dot      
!
!====================================================================== 

      subroutine rhob_inv_dense(o1,o2,o3,i1,i2,i3,h1,h2,tolh,nmxhi)

C     Outfld = (H2*B)-1 * inpfld  (implicit)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'VPROJ'
      include 'TSTEP'
 
      real o1 (lx1,ly1,lz1,1)
      real o2 (lx1,ly1,lz1,1)      
      real o3 (lx1,ly1,lz1,1)
      real i1 (lx1,ly1,lz1,1) 
      real i2 (lx1,ly1,lz1,1) 
      real i3 (lx1,ly1,lz1,1) 
      real h1 (lx1,ly1,lz1,1)
      real h2 (lx1,ly1,lz1,1)

      real tolh
      integer nmxhi
 
      imesh = 1

      if (ifield.eq.1) then
        call cggo_dense(o1,i1,h1,h2,v1mask,vmult,imesh,tolh,nmxhi,1,
     $            h2,'B-1x')

        call cggo_dense(o2,i2,h1,h2,v1mask,vmult,imesh,tolh,nmxhi,1,
     $            h2,'B-1y')
      
        if (if3d)
     $     call cggo_dense(o3,i3,h1,h2,v1mask,vmult,imesh,tolh,nmxhi,1,
     $            h2,'B-1z')


      endif

      return
      end
c--------------------------------------------------------------------
      subroutine cggo_dense(x,f,h1,h2,mask,mult,imsh,tin,maxit,isd,
     $                      rho,name)

C     Solve the Helmholtz equation, H*U = RHS,
C     using preconditioned conjugate gradient iteration.
C     Preconditioner: diag(H).

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'   ! ifield
!      include 'TOTAL'
      include 'MASS'
      include 'FDMH1'
 
      COMMON  /CPRINT/ IFPRINT, IFHZPC
      LOGICAL          IFPRINT, IFHZPC

      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      logical ifdfrm, iffast, ifh2, ifsolv

      logical ifmcor,ifprint_hmh
 
      real x(1),f(1),h1(1),h2(1),mask(1),mult(1)
      integer lg
      parameter (lg=lx1*ly1*lz1*lelt)

      real d,scalar,r,w,p,z
      COMMON /SCRCG/ d (lg) , scalar(2)
      common /SCRMG/ r (lg) , w (lg) , p (lg) , z (lg)

      real rho(lg)      ! on Mesh 1

      integer maxcg
      parameter (maxcg=900)

      real diagt,upper
      common /tdarray/ diagt(maxcg),upper(maxcg)

      integer niterhm
      common /iterhm/ niterhm
      character*4 name

      integer maxit,imsh,isd
      real tin,vol

      real alpha,beta,sigma,smean,rtz2,rtz1,skmin,rmean
      integer iter,niter,krylov,iter_max

      integer n,nel,nxyz
      real fmax,h2max,rbn0,rbn2,sigma0,tol

!     functions      
      real glmax,glamax,glmin,vlsc3,vlsc32,glsc2,glsc3,glsum


c **  zero out stuff for Lanczos eigenvalue estimator
      call rzero(diagt,maxcg)
      call rzero(upper,maxcg)
      sigma = 0.00
C
C     Initialization
C
      NXYZ   = lx1*ly1*lz1
      NEL    = NELV
      VOL    = VOLVM1
      IF (IMSH.EQ.2) NEL=NELT
      IF (IMSH.EQ.2) VOL=VOLTM1
      n      = NEL*NXYZ

      tol=abs(tin)

c     overrule input tolerance
      if (restol(ifield).ne.0) tol=restol(ifield)

      if (tin.lt.0) tol=abs(tin)
      niter = 10 !min(maxit,maxcg)

      if (.not.ifsolv) then
!        call setfast(h1,h2,imsh)
        ifsolv = .true.
      endif

C     Set up diag preconditioner.
      call copy(D,bm1,n)
      call col2(D,rho,n)
      call dssum(D,lx1,ly1,lz1)
      call invcol1(D,n)

      call copy (r,f,n)
      call rzero(x,n)
      call rzero(p,n)

      fmax = glamax(f,n)
      if (fmax.eq.0.0) return

c     Check for non-trivial null-space

      ifmcor = .false.
      h2max = glmax(h2  ,n)
      skmin = glmin(mask,n)

      krylov = 0
      rtz1=1.0
      niterhm = 0

      do iter=1,niter

!        Apply Diagonal preconditioner         
         call col3(z,r,d,n)

         rtz2=rtz1
         scalar(1)=vlsc3 (z,r,mult,n)
         scalar(2)=vlsc32(r,mult,binvm1,n)
         call gop(scalar,w,'+  ',2)
         rtz1=scalar(1)
         rbn2=sqrt(scalar(2)/vol)
         if (iter.eq.1) rbn0 = rbn2
         if (param(22).lt.0) tol=abs(param(22))*rbn0
         if (tin.lt.0)       tol=abs(tin)*rbn0

         ifprint_hmh = .false.
         if (nio.eq.0.and.ifprint.and.param(74).ne.0) ifprint_hmh=.true.
         if (nio.eq.0.and.istep.eq.1)                 ifprint_hmh=.true.

         if (ifprint_hmh)
     &      write(6,3002) istep,' Lanczos ' // name,
     &                    iter,rbn2,h1(1),tol,h2(1),ifmcor


         iter_max = 10000
         if (iter.gt.iter_max) then
           NITER = ITER-1
           if (nio.eq.0)
     &        write(6,3000) istep,' Lanczos ' // name,
     &                      iter_max,rbn2,rbn0,tol
           goto 9999
         endif
c
         beta = rtz1/rtz2
         if (iter.eq.1) beta=0.0
         call add2s1 (p,z,beta,n)

!         call axhelm (w,p,h1,h2,imsh,isd)
         call rho_b(w,p,rho)

         call dssum  (w,lx1,ly1,lz1)
         call col2   (w,mask,n)
c
         sigma0 = sigma
         sigma  = glsc3(w,p,mult,n)
         alpha=rtz1/sigma
         call add2s2(x,p ,alpha,n)
         call add2s2(r,w ,-alpha,n)

c        Generate tridiagonal matrix for Lanczos scheme
         if (iter.eq.1) then
            krylov = krylov+1
            diagt(iter) = sigma/rtz1
         elseif (iter.le.maxcg) then
            krylov = krylov+1
            diagt(iter)    = (beta**2 * sigma0 + sigma ) / rtz1
            upper(iter-1)  = -beta * sigma0 / sqrt(rtz2 * rtz1)
         endif
 1000 enddo
      niter = iter-1
c
      if (nio.eq.0) write (6,3001) istep, '  Error Hmholtz ' // name,
     &                             niter,rbn2,rbn0,tol


 3000 format(i11,a,1x,I7,1p4E13.4)
 3001 format(i11,a,1x,I7,1p4E13.4)
 3002 format(i11,a,1x,I7,1p4E13.4,l4)
 9999 continue
      niterhm = niter
      ifsolv = .false.

      return
      end
!---------------------------------------------------------------------- 

      subroutine rho_b(Au,u,rho)

      implicit none

      include 'SIZE'
      include 'GEOM'    ! JACM1
      
      integer lxb
      parameter (lxb=2*lx2)   ! Arbitrarily set for now

!     Interpolation matrices for dense mass matrices
      real jgl(lxb,lx1)
      real jglt(lx1,lxb)
      real wght(lxb)
      common /densemass/ jgl,jglt,wght

      integer e
      real u(lx1*ly1*lz1,lelv)
      real Au(lx1*ly1*lz1,lelv)
      real rho(lx1*ly1*lz1,lelv)

      real wk1lxb(lxb*lxb*lxb)
      real wk2lxb(lxb*lxb*lxb)
      real ulxb(lxb*lxb*lxb)

      integer lxb3

      integer i,j,k,ii
      integer i1,i2,iz

      integer icalld
      save icalld
      data icalld /0/

      logical ifcylindrical

      ifcylindrical = .false.

      if (icalld.eq.0) then
        call setup_interp_lxb(jgl,jglt,wght,lxb)
        icalld = icalld + 1
      endif  

      lxb3 = lxb*lxb*lxb

      do e=1,nelv

!       u_lx1 -> u_lxb          
        call intp_lxb(ulxb,u(1,e),jgl,jglt,lxb) 

!       \rho*u
        call intp_lxb(wk1lxb,rho(1,e),jgl,jglt,lxb)
        call col2(ulxb,wk1lxb,lxb3)

!       J*\rho*u
        call intp_lxb(wk2lxb,jacm1(1,1,1,e),jgl,jglt,lxb)
        call col2(ulxb,wk2lxb,lxb3)

        if (ifcylindrical) then
!         R*J*\rho*u
          call intp_lxb(wk2lxb,ym1(1,1,1,e),jgl,jglt,lxb)
          call col2(ulxb,wk2lxb,lxb3)
        endif  

!       W*J*\rho*u
        do k=1,lxb
        do j=1,lxb
        do i=1,lxb
          ii = i + (j-1)*lxb + (k-1)*lxb*lxb
          ulxb(ii) = ulxb(ii)*wght(i)*wght(j)*wght(k)
        enddo
        enddo
        enddo

!       (I^T)*W*J*\rho*u
        call mxm (jglt,lx1,ulxb,lxb,wk1lxb,lxb*lxb)
        i1=1
        i2=1
        do iz=1,lxb
          call mxm (wk1lxb(i1),lx1,jgl,lxb,wk2lxb(i2),ly1)
          i1=i1+lx1*lxb
          i2=i2+lx1*ly1
        enddo
        call mxm  (wk2lxb,lx1*ly1,jgl,lxb,Au(1,e),lz1)

      enddo

      return
      end subroutine rho_b
!---------------------------------------------------------------------- 
      subroutine setup_interp_lxb(jgl,jglt,wght,lxb)

      implicit none

c     Set up interpolator from lx1 to lxb mesh 

      include 'SIZE'
      include 'WZ'

      integer lxb

      real jgl(lxb,lx1)
      real jglt(lx1,lxb)

      real wght(lxb)
      real z(lxb)

!     Gauss-Legendre Mesh 
      call zwgl (z,wght,lxb)

!     Interpolator from M1 mesh to lxb Mesh      
      call igllm (jgl,jglt,zgm1(1,1),z,lx1,lxb,lx1,lxb)

      return
      end
c-----------------------------------------------------------------------

      subroutine intp_lxb(fldf,fld,jgl,jglt,lxb)

      implicit none

      include 'SIZE'
      include 'GEOM'

      integer lxb

      real fld  (lx1*ly1*lz1)
      real fldf (lxb*lxb*lxb)

      real jgl(lxb,lx1)
      real jglt(lx1,lxb)

      real wk1(lxb*lxb*lxb)
      real wk2(lxb*lxb*lxb)

      integer iz,i1,i2


      if (ndim.eq.2) then
!       I12*u        
        call mxm (jgl,lxb,fld,lx1,wk1,ly1)
        call mxm (wk1,lxb,jglt,ly1,fldf,lxb)
      else        
!       I12*u        
        call mxm (jgl,lxb,fld,lx1,wk1,ly1*lz1)
        i1=1
        i2=1
        do iz=1,lz1
          call mxm (wk1(i1),lxb,jglt,ly1,wk2(i2),lxb)
          i1=i1+lxb*ly1
          i2=i2+lxb*lxb
        enddo
        call mxm  (wk2,lxb*lxb,jglt,lz1,fldf,lxb)
      endif        


      return
      end subroutine intp_lxb
!---------------------------------------------------------------------- 

      subroutine cyl_gmres(res,h2,mask,maxiter)

!     using right-preconditioned GMRES iteration.
!     to solve \rho*B*u = f

!     We solve:   \rho*B*(M^-1)(Mu) = f,
!              => \rho*B*(M^-1)v    = f,
!              => u = (M^-1)v

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'GMRES_M1'

      integer lt1,lt2
      parameter (lt1 = lx1*ly1*lz1*lelv)
      parameter (lt2 = lx2*ly2*lz2*lelv)

      logical          ifprint
      common  /cprint/ ifprint

      real             res  (lt1)
      real             h2   (lt1)
      real             mask (lt1)

      real theta        ! Approximated eigenvalue

      real wp
      common /scrmg_m1/    wp (lt1)

      real wk1,wk2
      common /ctmp0_m1/   wk1(lgmres),wk2(lgmres)

      real y
      common /cgmres1_m1/ y(lgmres)

      real alpha, l, temp
      real t1,t2
      integer j,m

      logical iflag
      save    iflag
      data    iflag /.false./
      real    norm_fac
      save    norm_fac

      real*8 etime1,dnekclock

      integer ntot
      real glsc2,glsc3,vlsc2,vlsc3
      integer iconv,intype
      real tolpss,div0
      integer i,k,iter
      real etime2,etime_p,ratio,rnorm

      integer maxiter

      logical ifwgt           ! If weighted Gmres       
      logical ifprec

      integer ngs             ! No of Gram-Schmidt
      character*132 str

      real D(lt1)

      ntot        = lx1*ly1*lz1*nelv

      ifprint     = .true.

      ifwgt       = .false.           ! Weighted Orthogonalization
      ifprec      = .true.           ! Use preconditioner
      ngs         = 2

      if (ifwgt) then
        norm_fac = 1./sqrt(volvm1)
      else
        call rone(wp,ntot)
        alpha = sqrt(glsc2(wp,wp,ntot))
        norm_fac = 1.0/alpha
      endif

      norm_fac = 1.0 

C     Set up diag preconditioner.
      if (ifprec) then
        call copy(D,bm1,ntot)
        call col2(D,h2,ntot)
        call dssum(D,lx1,ly1,lz1)
        call invcol1(D,ntot)
      else
        call rone(D,ntot)  
      endif  

      etime1 = dnekclock()
      etime_p = 0.
      iter  = 0
      m = min(maxiter,lgmres)

      tolpss = 1.0e-12

      iconv = 0
      call rzero  (x_gmres,ntot)
      call dssum  (res,lx1,ly1,lz1)
      call col2   (res,mask,ntot)

      do while(iconv.eq.0.and.iter.lt.maxiter)            ! prabal

         if(iter.eq.0) then
            call copy(r_gmres,res,ntot)                ! r = res
         else
!           update residual
            call copy(r_gmres,res,ntot)                ! r = res

!           w = A*x
            call rho_b  (w_gmres,x_gmres,h2)
            call dssum  (w_gmres,lx1,ly1,lz1)
            call col2   (w_gmres,mask,ntot)

!           r = r - w
            call add2s2(r_gmres,w_gmres,-1.,ntot)

         endif

         if (ifwgt) then
!          Weighted inner product                                !            ______
           gamma_gmres(1) = sqrt(glsc3(r_gmres,r_gmres,bm1,ntot))! gamma  = \/(Br,r)
         else    
!          Un-weighted inner product                         !            ______
           gamma_gmres(1) = sqrt(glsc2(r_gmres,r_gmres,ntot))! gamma  = \/(r,r)
         endif   
                                                           
         if(iter.eq.0) then
           div0 = gamma_gmres(1)*norm_fac
         endif

!        check for lucky convergence
         rnorm = 0.
         if (gamma_gmres(1) .eq. 0.) goto 9000
         temp = 1./gamma_gmres(1)
         call cmult2(v_gmres(1,1),r_gmres,temp,ntot)  ! v  = r / gamma

         do j=1,m
            iter = iter+1

            call copy(w_gmres,v_gmres(1,j),ntot)      ! w  = v_j
 
            etime2 = dnekclock()
!           z = (M^-1)w      
            if (ifprec) then
              call col3(z_gmres(1,j),w_gmres,D,ntot)
            else
              call copy(z_gmres(1,j),w_gmres,ntot)
            endif

            etime_p = etime_p + dnekclock()-etime2

!           r = (M^-1)w      
            call copy(r_gmres,z_gmres(1,j),ntot)

!           w = A*(M^-1)w
            call rho_b(w_gmres,r_gmres,h2) 
            call dssum  (w_gmres,lx1,ly1,lz1)
            call col2   (w_gmres,mask,ntot)

!           Gram-Schmidt:
            call ortho_subspace(w_gmres,ntot,h_gmres(1,j),v_gmres,
     $            lt1,j,bm1,ifwgt,ngs,wk1,wk2)
!            call ortho_subspace_fullM1(w_gmres,ntot,h_gmres(1,j),
!     $            v_gmres,lt1,j,bm1,ifwgt,ngs,wk1,wk2,.true.)

!           Apply Givens rotations to new column
            do i=1,j-1
              t1             = h_gmres(i,j)                   
              t2             = h_gmres(i+1,j)                   
              h_gmres(i  ,j) =  c_gmres(i)*t1 + s_gmres(i)*t2 
              h_gmres(i+1,j) = -s_gmres(i)*t1 + c_gmres(i)*t2
            enddo
!                      ______
!           alpha =  \/(Bw,w) 
            if (ifwgt) then
              alpha = sqrt(glsc3(w_gmres,w_gmres,bm1,ntot))    
            else
              alpha = sqrt(glsc2(w_gmres,w_gmres,ntot))        
            endif
            if(alpha.eq.0.) goto 900  !converged

!           Calculate new Givens rotation to convert H to Triangular
            l = sqrt(h_gmres(j,j)*h_gmres(j,j)+alpha*alpha)
            temp = 1./l
            c_gmres(j) = h_gmres(j,j) * temp
            s_gmres(j) = alpha  * temp
            h_gmres(j,j) = l

!           Apply New Givens rotation to the rhs
            t1 =  c_gmres(j) * gamma_gmres(j)
            t2 = -s_gmres(j) * gamma_gmres(j)
            gamma_gmres(j)   = t1 
            gamma_gmres(j+1) = t2 

            rnorm = abs(gamma_gmres(j+1))*norm_fac
            ratio = rnorm/div0
            if (ifprint.and.nio.eq.0) 
     $         write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66       format(i5,1p4e12.5,i8,' (rho*B)^-1 ')

            if (rnorm .lt. tolpss) goto 900  !converged

            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v_gmres(1,j+1),w_gmres,temp,ntot) ! v = w / alpha
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

!        sum up Arnoldi vectors
!        x_gmres = x_gmres + (M^-1)*(V*c)
!     => x_gmres = x_gmres + (M^-1 * V)*c = X + Z*c
         do i=1,j
            call add2s2(x_gmres,z_gmres(1,i),c_gmres(i),ntot) 
         enddo

      enddo
 9000 continue

      call copy(res,x_gmres,ntot)


      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep,' (rho*B)^-1 gmres  ', 
     &                            iter,rnorm,div0,tolpss,etime_p,etime1

 9999 format(i11,a,I6,1p5e13.4)

      return
      end

c-----------------------------------------------------------------------

      subroutine ortho_subspace_fullM1(r,nt,h,V,ldv,k,wgt,ifwgt,ngs
     $                               ,wk1,wk2,iffull)

      implicit none

      integer nt              ! Length of the vector r
      integer ldv             ! Leading dimension of V
      integer k               ! No of Columns in V
      real r(nt)              ! Vector to orthogonalize
      real V(ldv,k)           ! Orthogonalizing Space
      real wgt(nt)            ! Weights
      logical ifwgt           ! If Weighted orthogonalization
      integer ngs             ! No. of Gram-Schmidt
      real h(k)               ! Projections on V
      logical iffull          ! Full Mass matrix for orthogonalization

!     Work Arrays      
      real wk1(k)
      real wk2(k)

      integer igs,i

      real vlsc2,vlsc3        ! Functions
      real vlsc2_fm

!     Zero projections      
      call rzero(h,k)

      do igs = 1,ngs
!       Gram-Schmidt:
        do i=1,k
          if (ifwgt) then
            if (iffull) then
              wk1(i) = vlsc2_fm(r,V(1,i))
            else  
              wk1(i)=vlsc3(r,V(1,i),wgt,nt)     ! wk1 = (Bw,V )
            endif  
          else
            wk1(i)=vlsc2(r,V(1,i),nt)           ! wk1 = (w,V )
          endif
        enddo                                             
        call gop(wk1,wk2,'+  ',k)               ! sum over all procs

        do i=1,k
          call add2s2(r,V(1,i),-wk1(i),nt)      ! r = r - V*wk1
          h(i) = h(i) + wk1(i)                  ! h = h + wk1 
        enddo
      enddo       ! igs 

      return
      end subroutine ortho_subspace_fullM1
!----------------------------------------------------------------------

      function glsc2_full_M1(u,v)

      implicit none

      include 'SIZE'

      real u(1),v(1)
      real sc,tmp
      real glsc2_full_M1

      call local_inner_prod_full_M1(sc,u,v)
      call gop(sc,tmp,'+  ',1)

      glsc2_full_M1 = sc

      return 
      end function 

!---------------------------------------------------------------------- 
      subroutine local_inner_prod_full_M1(s,u,v)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'GEOM'

      real s
      real u(lx1*ly1*lz1,lelt)
      real v(lx1*ly1*lz1,lelt)

      integer lxb
      parameter (lxb=2*lx1)   ! Arbitrarily set for now

!     Interpolation matrices for dense mass matrix inner product
      real jgl(lxb,lx1)
      real jglt(lx1,lxb)
      real wght(lxb)
      common /denseinp/ jgl,jglt,wght     ! dense inner product

      real wk1(lxb*lxb*lxb)
      real wk2(lxb*lxb*lxb)

      integer e,i,j,k,ii
      integer lxb3

      real vlsum

      integer icalld
      save icalld
      data icalld /0/

      if (icalld.eq.0) then
        call setup_interp_lxb(jgl,jglt,wght,lxb)
        icalld = icalld + 1
      endif  

      lxb3 = lxb*lxb*lxb

      s = 0.0
      do e=1,nelv

!       u_lx1 -> u_lxb          
        call intp_lxb(wk1,u(1,e),jgl,jglt,lxb) 

!       v_lx1 -> v_lxb
        call intp_lxb(wk2,v(1,e),jgl,jglt,lxb)
        call col2(wk1,wk2,lxb3)

!       J_lx1 -> J_lxb
        call intp_lxb(wk2,jacm1(1,1,1,e),jgl,jglt,lxb)
        call col2(wk1,wk2,lxb3)

!       W*J*\rho*u
        do k=1,lxb
        do j=1,lxb
        do i=1,lxb
          ii = i + (j-1)*lxb + (k-1)*lxb*lxb
          wk1(ii) = wk1(ii)*wght(i)*wght(j)*wght(k)
        enddo
        enddo
        enddo

        s = s + vlsum(wk1,lxb3) 

      enddo

      return
      end subroutine local_inner_prod_full_M1
!---------------------------------------------------------------------- 

      subroutine op_gmres(r1,r2,r3,h2,maxiter)

!     using right-preconditioned GMRES iteration.
!     to solve \rho*B*u = f

!     We solve:   \rho*B*(M^-1)(Mu) = f,
!              => \rho*B*(M^-1)v    = f,
!              => u = (M^-1)v

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'OP_GMRES_M1'

      integer lt1,lt2
      parameter (lt1 = lx1*ly1*lz1*lelv)
      parameter (lt2 = lx2*ly2*lz2*lelv)

      logical          ifprint
      common  /cprint/ ifprint

      real             r1    (lt1)
      real             r2    (lt1)
      real             r3    (lt1)
      real             h2    (lt1)

      real theta        ! Approximated eigenvalue

      real wghts
      common /scrmg_m1/    wghts (lt1*3)

      real wk1,wk2
      common /ctmp0_m1/   wk1(lgmres),wk2(lgmres)

      real y
      common /cgmres1_m1/ y(lgmres)

      real alpha, l, temp
      integer j,m

      logical iflag
      save    iflag
      data    iflag /.false./
      real    norm_fac
      save    norm_fac

      real*8 etime1,dnekclock

      integer nt1,nt3,lt3
      real glsc2,glsc3,vlsc2,vlsc3,op_glsc2_wt
      integer iconv,intype
      real tolpss,div0
      integer i,k,iter
      real etime2,etime_p,ratio,rnorm

      integer maxiter

      logical ifwgt           ! If weighted Gmres       
      logical ifprec

      integer ngs             ! No of Gram-Schmidt
      character*132 str

      real D1(lt1),D2(lt1),D3(lt1)

      integer i1,i2,i3

      nt1         = lx1*ly1*lz1*nelv
      nt3         = nt1*3
      lt3         = lt1*3

      i1          = 1
      i2          = 1 + nt1
      i3          = 1 + 2*nt1

      ifprint     = .true.

      ifwgt       = .false.           ! Weighted Orthogonalization
      ifprec      = .true.           ! Use preconditioner
      ngs         = 1

!     I use this for weights
      if (ifwgt) then
        call opcopy(wghts(i1),wghts(i2),wghts(i3),bm1,bm1,bm1)
      else
        call rone(wghts,nt3)
      endif  
      alpha = sqrt(glsc2(wghts,wghts,nt3))
      norm_fac = 1.0/alpha

C     Set up diag preconditioner.
      if (ifprec) then
        call col3(D1,bm1,h2,nt1)
        call col3(D2,bm1,h2,nt1)
        call col3(D3,bm1,h2,nt1)
        call opdssum(D1,D2,D3)
        call invcol1(D1,nt1)
        call invcol1(D2,nt1)
        call invcol1(D3,nt1)
      else
        call rone(D1,nt1)  
        call rone(D2,nt1)  
        call rone(D3,nt1)  
      endif  

      etime1 = dnekclock()
      etime_p = 0.
      iter  = 0
      m = min(maxiter,lgmres)

      tolpss = 1.0e-12

      iconv = 0
      call rzero   (x_gmres,lt3)
      call opdssum (r1,r2,r3)
      call opmask  (r1,r2,r3)

      do while(iconv.eq.0.and.iter.lt.maxiter)

         if(iter.eq.0) then
            call opcopy(r_gmres(i1),r_gmres(i2),r_gmres(i3),
     $                  r1,r2,r3)

         else
            !update residual
            call opcopy(r_gmres(i1),r_gmres(i2),r_gmres(i3),
     $                  r1,r2,r3)

!           w = A*x
            call rho_b(w_gmres(i1),x_gmres(i1),h2)
            call rho_b(w_gmres(i2),x_gmres(i2),h2)
            call rho_b(w_gmres(i3),x_gmres(i3),h2)

            call opdssum(w_gmres(i1),w_gmres(i2),w_gmres(i3))
            call opmask(w_gmres(i1),w_gmres(i2),w_gmres(i3))

!           r = r - w
            call add2s2(r_gmres,w_gmres,-1.,nt3)

         endif

         if (ifwgt) then
!          Weighted inner product                                
!                     ______
!          gamma  = \/(Br,r) 
           gamma_gmres(1) = sqrt(glsc3(r_gmres,r_gmres,wghts,nt3))
         else    
!          Un-weighted inner product
!                     ______
!          gamma  = \/(r,r)                          
           gamma_gmres(1) = sqrt(glsc2(r_gmres,r_gmres,nt3))
         endif   
                                                           
         if(iter.eq.0) then
           div0 = gamma_gmres(1)*norm_fac
           if (param(21).lt.0) tolpss=abs(param(21))*div0
         endif

!        check for lucky convergence
         rnorm = 0.
         if (gamma_gmres(1) .eq. 0.) goto 9000
         temp = 1./gamma_gmres(1)
         call cmult2(v_gmres(1,1),r_gmres,temp,nt3)! v  = r / gamma

         do j=1,m
            iter = iter+1

            call copy(w_gmres,v_gmres(1,j),nt3) ! w  = v_j
 
            etime2 = dnekclock()
!           z = (M^-1)w      
            if (ifprec) then
              call col3(z_gmres(i1,j),w_gmres(i1),D1,nt1)
              call col3(z_gmres(i2,j),w_gmres(i2),D2,nt1)
              call col3(z_gmres(i3,j),w_gmres(i3),D3,nt1)
            else
              call copy(z_gmres(1,j),w_gmres,nt3)
            endif

            etime_p = etime_p + dnekclock()-etime2

!           r = (M^-1)w 
            call copy(r_gmres,z_gmres(1,j),nt3)

!           w = A*(M^-1)w
            call rho_b(w_gmres(i1),r_gmres(i1),h2)
            call rho_b(w_gmres(i2),r_gmres(i2),h2)
            call rho_b(w_gmres(i3),r_gmres(i3),h2)

            call opdssum  (w_gmres(i1),w_gmres(i2),w_gmres(i3))
            call opmask   (w_gmres(i1),w_gmres(i2),w_gmres(i3))

!           Gram-Schmidt:
            call ortho_subspace(w_gmres,nt3,h_gmres(1,j),v_gmres,
     $            lt3,j,wghts,ifwgt,ngs,wk1,wk2)

!           Apply Givens rotations to new column
            do i=1,j-1
              temp = h_gmres(i,j)                   
              h_gmres(i  ,j)=  c_gmres(i)*temp 
     $                       + s_gmres(i)*h_gmres(i+1,j)  
              h_gmres(i+1,j)= -s_gmres(i)*temp 
     $                       + c_gmres(i)*h_gmres(i+1,j)
            enddo

            if (ifwgt) then
!                        ______
!             alpha =  \/(Bw,w) 
              alpha = sqrt(glsc3(w_gmres,w_gmres,wghts,nt3))
            else
!                        ______
!             alpha =  \/(w,w) 
              alpha = sqrt(glsc2(w_gmres,w_gmres,nt3))        
            endif
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
   66       format(i5,1p4e12.5,i8,' (rho*B)^-1 ')

            if (rnorm .lt. tolpss) goto 900  !converged

            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v_gmres(1,j+1),w_gmres,temp,nt3) ! v = w / alpha
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
!        sum up Arnoldi vectors
!        x_gmres = (M^-1)*(V*c)
!     => x_gmres = (M^-1 * V)*c = Z*c
         do i=1,j
!          x = x + Z*c
           call add2s2(x_gmres,z_gmres(1,i),c_gmres(i),nt3) 
         enddo

      enddo
 9000 continue

      call copy(r1,x_gmres(i1),nt1)
      call copy(r2,x_gmres(i2),nt1)
      call copy(r3,x_gmres(i3),nt1)

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep,' (rho*B)^-1 gmres  ', 
     &                            iter,rnorm,div0,tolpss,etime_p,etime1

 9999 format(i11,a,I6,1p5e13.4)

      return
      end

c-----------------------------------------------------------------------





















