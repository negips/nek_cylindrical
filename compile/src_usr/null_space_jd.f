!======================================================================
!
!     Author: Prabal Negi
!     Description: Find the Null Space of the E operator
!                  Using Jacobi Davidson
!
!======================================================================
!---------------------------------------------------------------------- 

      subroutine jacobi_davidson_e()

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'SOLN'          ! vtrans
      include 'TSTEP'         ! ifield,DT

      include 'TEST'

      integer jdlv
      parameter (jdlv = 30)

      integer lt2
      parameter(lt2=lx2*ly2*lz2*lelv)

      real jd_V(lt2,jdlv)           ! Space V
      real jd_W(lt2,jdlv)           ! Space W = A*V
      real jd_H(jdlv,jdlv)

      real jdv(lt2)
      real jdw(lt2)
      real jdr(lt2)
      real jdu0(lt2)                ! Last eigen direction

!     For Eigen Decomposition
      real ev_H(jdlv,jdlv)
      real ev_vr(jdlv,jdlv)
      real ev_vl(jdlv,jdlv)
      real ev_wr(jdlv)
      real ev_wi(jdlv)

      real e0r,e0i

      real diff


      integer jdnv
      integer nt1,nt2

      real hii
      real theta

      real rnd
      real vnorm,vnormi
      real glsc2,glsc3

      real h1,h2,h2inv
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)

      real h(jdlv),ht(jdlv),wk1(jdlv),wk2(jdlv)

      integer intype

      logical ifwgt           ! If Weighted orthogonalization
      integer ngs             ! No of Gram-Schmidts

      integer i,j,k,ijd,e0
      integer maxiter         ! Max. GMRES iterations

      integer jd_restarts     ! Max Jacobi Davidson restarts
      logical ifconv

      real eigtol
      real mu                 ! temporary (I - mu*E)

      nt1 = lx1*ly1*lz1*nelv
      nt2 = lx2*ly2*lz2*nelv

      ifwgt       = .true.
      ngs         = 2
      maxiter     = 1000
      jd_restarts = 200
      eigtol      = 1.0e-12

      mu          = -dt

      ifield = 1
      intype = 1

      call rzero   (h1,nt1)
!      call copy    (h2,vtrans(1,1,1,1,ifield),nt1)
      call rone    (h2,nt1)
      call invers2 (h2inv,h2,nt1)

      call rzero(jd_V,lt2*jdlv)
      call rzero(jd_W,lt2*jdlv)
      call rzero(jd_H,jdlv*jdlv)

      jdnv = jdlv - 1

!     Generate random initial vector      
      do i=1,nt2
        call random_number(rnd)
        jdv(i) = rnd
      enddo

!     Normalize initial vector
      if (ifwgt) then 
        vnorm = sqrt(glsc3(jdv,jdv,bm2,nt2))
      else
        vnorm = sqrt(glsc2(jdv,jdv,nt2))
      endif

      vnormi = 1.0/vnorm
      call cmult(jdv,vnormi,nt2)
      call copy(jd_V(1,1),jdv,nt2)

!     w = E*v      
      call jd_Eop(jdw,jdv,h1,h2,h2inv,intype) 
      call copy(jd_W(1,1),jdw,nt2)

!     h = V'*B2*w
      if (ifwgt) then      
        hii = glsc3(jd_V(1,1),jdw,bm2,nt2)
      else
        hii = glsc2(jd_V(1,1),jdw,nt2)
      endif
      jd_H(1,1) = hii

!     Orthogonalizing direction
!     Also stores the best approximation of our wanted eigenvector      
      call copy(jdu0,jdv,nt2) 

      theta = 1.0 ! 1.0e-12

!     r = w - \theta*v      
      call add3s2(jdr,jdw,jdv,1.0,theta,nt2)

      ijd = 0
      ifconv = .false.

      if (abs(hii-theta).lt.eigtol) then
        ifconv = .true.
      endif  

      istep = 0
      do while(ijd.lt.jd_restarts.and..not.ifconv)

        ijd = ijd+1

        do j=1,jdnv
!         r = -r      
          call chsign(jdr,nt2)
          
          istep = istep+1
!         Correction.
!         (I - u0*u0')*(A - \theta*I)*(I - u0*u0')t = -r 
          call jd_gmres(jdr,jdu0,theta,h1,h2,h2inv,intype,maxiter)

!         Orthogonalize t with respect to V
!         h = V'*t        
!         t = t - V*h
          call ortho_subspace(jdr,nt2,h,jd_V,lt2,j,
     $                        bm2,ifwgt,ngs,wk1,wk2)

!         Normalize t
          if (ifwgt) then 
            vnorm = sqrt(glsc3(jdr,jdr,bm2,nt2))
          else
            vnorm = sqrt(glsc2(jdr,jdr,nt2))
          endif
          vnormi = 1.0/vnorm

          call copy(jdv,jdr,nt2)
          call cmult(jdv,vnormi,nt2)
          call copy(jd_V(1,j+1),jdv,nt2)    ! Extend subspace V

!         w = A*v
          call jd_Eop(jdw,jdv,h1,h2,h2inv,intype) 
          call copy(jd_W(1,j+1),jdw,nt2)    ! Extend subspace W

!         Get Projections of w on to V to get the column H(:,j+1)
!         h = V'*w
          call ortho_subspace(jdw,nt2,jd_H(1,j+1),jd_V,lt2,j+1,
     $                        bm2,ifwgt,ngs,wk1,wk2)

!         Get Projections of v on to W to get the rows H(j+1,1:j)
!         ht = v'*W
          call ortho_subspace(jdv,nt2,ht,jd_W,lt2,j,
     $                        bm2,ifwgt,ngs,wk1,wk2)

          do i=1,j
            jd_H(j+1,i) = ht(i)
          enddo

!         Calculate Eigen values/Eigenvectors      
          call copy(ev_H,jd_H,jdlv*jdlv)       ! Gets overwritten in dgeev
          call wrp_dgeev(ev_H,jdlv,j+1,ev_wr,ev_wi,ev_vl,jdlv,
     $                   ev_vr,jdlv)

          diff = 1e+5
          e0 = 1
          do i=1,j+1
            wk1(i) = sqrt(ev_wr(i)**2 + ev_wi(i)**2)
            if (abs(wk1(i) - theta).lt.diff) then
              diff = abs(wk1(i)-theta)
              e0 = i
            endif  
          enddo
          e0r = (ev_wr(e0) - 1.0)/mu
          e0i = (ev_wi(e0) - 0.0)/mu 
          if (nio.eq.0) write(6,16) 'EIG0',istep,
     $                               ev_wr(e0), ev_wi(e0),e0r,e0i
          if (nio.eq.0) write(6,*) ' '

!         Taking the eigenvector corresponding to the
!         nearest eigenvalue as u0. 
          call mxm(jd_V,lt2,ev_vr(1,e0),j+1,jdu0,1)

!         Normalize new eigenvector
          if (ifwgt) then 
            vnorm = sqrt(glsc3(jdu0,jdu0,bm2,nt2))
          else
            vnorm = sqrt(glsc2(jdu0,jdu0,nt2))
          endif
          vnormi = 1.0/vnorm
          call cmult(jdu0,vnormi,nt2)

          call copy(jdv,jdu0,nt2)
          call jd_Eop(jdw,jdv,h1,h2,h2inv,intype) 

!         r = w - \theta*v      
          call add3s2(jdr,jdw,jdv,1.0,theta,nt2)

          diff = abs(ev_wr(e0)-theta)
          if (diff.lt.eigtol) then
            ifconv = .true.
            if (nio.eq.0) write(6,*) 'Jacobi Davison Converged, res=',
     $                                diff             
          endif  

!         Transform Eigenvalues.
          do i=1,j+1
            ev_wr(i) = (ev_wr(i)-1.0)/mu
            ev_wi(i) = (ev_wi(i)-0.0)/mu
            if (nio.eq.0) write(6,16) 'JDEIG',i,ev_wr(i), ev_wi(i)
          enddo

          if (ifconv) exit 
        enddo       ! j=1,jdnv

!       Outpost current null space estimate
        time  = ev_wr(i)
        istep = ijd
        call outpost(vx,vy,vz,jdu0,t,'est')

        if (.not.ifconv) then
!         Reinitialize matrices
          call copy(jd_V(1,1),jdv,nt2)
          call copy(jd_W(1,1),jdw,nt2)

          call rzero(jd_H,jdlv*jdlv)

!         h = V'*B2*w
          if (ifwgt) then      
            hii = glsc3(jd_V(1,1),jdw,bm2,nt2)
          else
            hii = glsc2(jd_V(1,1),jdw,nt2)
          endif
          jd_H(1,1) = hii
        endif       ! if (.not.ifconv)

      enddo         ! while (ijd.lt.jd_restarts.and..not.ifconv)


!     Output eigenvectors
      k = jdnv
      if (ifconv) then
        k = j 
        istep = i
        call mxm(jd_V,lt2,ev_vr(1,e0),k+1,pr,1)
        call outpost(vx,vy,vz,pr,t,'tst')
      else
        k = jdnv
        do i=1,k+1
          istep = i
          time  = ev_wr(i)
          call mxm(jd_V,lt2,ev_vr(1,i),k+1,pr,1)
          call outpost(vx,vy,vz,pr,t,'tst')
        enddo
      endif 

      do i=1,k+1
!!        We have already scaled these      
!         ev_wr(i) = (ev_wr(i)-1.0)/mu
!         ev_wi(i) = (ev_wi(i)-0.0)/mu
         if (nio.eq.0) write(6,16) 'JDEIG',i,ev_wr(i), ev_wi(i)
      enddo  

16    format(A5,2x,I5,2x,16(E22.14,2x))
17    format(A5,2x,16(E22.14,2x))



      return
      end subroutine jacobi_davidson_e
!---------------------------------------------------------------------- 

      subroutine jd_gmres(res,u0,theta,h1,h2,h2inv,intype,maxiter)

c     Solve the Jacobi-Davidson correction equation
c     for the E operator
c     using right-preconditioned GMRES iteration.
c     intype =  0  (steady)
c     intype =  1  (explicit)
c     intype = -1  (implicit)

!     We solve:   (U^T)*E*U*(M^-1)(Mp) = f,
!              => (U^T)*E*U*(M^-1)q    = f,
!              => p = (M^-1)q.
!     Where,      U = (I - u0*(u0^T)*B),
!     is the subspace orthogonal to u0.
!
!     ortho_left(r)  = (U^T)*r
!     ortho_right(r) = U*r

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'GMRES'

      integer lt1,lt2
      parameter(lt1 = lx1*ly1*lz1*lelv)
      parameter(lt2 = lx2*ly2*lz2*lelv)

      logical          ifprint
      common  /cprint/ ifprint

      real             res  (lt2)
      real             h1   (lt1)
      real             h2   (lt1)
      real             h2inv(lt1)

!     Orthogonalizing vector      
      real             u0   (lt2)

      real theta        ! Approximated eigenvalue

      real wp
      common /scrmg/    wp (lt2)

      real wk1,wk2
      common /ctmp0/   wk1(lgmres),wk2(lgmres)

      real y
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

      integer ntot2
      real glsc2,glsc3,vlsc2,vlsc3
      integer iconv,intype
      real tolpss,div0
      integer i,k,iter
      real etime2,etime_p,ratio,rnorm

      integer maxiter

      logical ifwgt           ! If weighted Gmres       
      logical ifwgtortho      ! If weighted ortho
      logical ifprec

      integer ngs             ! No of Gram-Schmidt
      character*132 str


      ntot2  = lx2*ly2*lz2*nelv

      ifprint = .true.

      ifwgt       = .false.            ! Weighted Ortho for GMRES
      ifwgtortho  = .false.
      ifprec      = .false.           ! Use preconditioner
      ngs         = 1

!     Orthogonalize w.r.t constant vector      
      call jd_ortho(res,u0,wp,0,ifwgt)

      if (ifwgt) then
        norm_fac = 1./sqrt(volvm2)
      else
        call rone(wp,ntot2)
        alpha = sqrt(glsc2(wp,wp,ntot2))
        norm_fac = 1.0/alpha
      endif

      etime1 = dnekclock()
      etime_p = 0.
      iter  = 0
      m = min(maxiter,lgmres)

      if (ifwgt) then
!       Weighted inner product               !            ______
        div0 = sqrt(glsc3(res,res,bm2,ntot2))! gamma  = \/(Br,r)
      else    
!       Un-weighted inner product        !            ______
        div0 = sqrt(glsc2(res,res,ntot2))! gamma  = \/(r,r)
      endif   
      div0 = div0*norm_fac

!      call chktcg2(tolps,res,iconv)
!      tolpss = tolps
      tolpss = 1.0e-8*div0

      iconv = 0
      call rzero(x_gmres,ntot2)

      do while(iconv.eq.0.and.iter.lt.maxiter)            ! prabal

         if(iter.eq.0) then
            call copy(r_gmres,res,ntot2)                ! r = res
         else
            !update residual
            call copy(r_gmres,res,ntot2)                ! r = res

!           Using z_gmres(1,1) as temporary variable
            call copy(z_gmres,x_gmres,ntot2)            ! z = x

!           U*z
            call jd_ortho(z_gmres,u0,wp,1,ifwgtortho)

!           w = A*U*z            
            call jd_Ecorr(w_gmres,z_gmres,h1,h2,h2inv,intype,theta)

!           w = (U^T)*A*U*x
            call jd_ortho(w_gmres,u0,wp,0,ifwgtortho)

            call add2s2(r_gmres,w_gmres,-1.,ntot2)    ! r = r - w

         endif

         if (ifwgt) then
!          Weighted inner product                                 !            ______
           gamma_gmres(1) = sqrt(glsc3(r_gmres,r_gmres,bm2,ntot2))! gamma  = \/(Br,r)
         else    
!          Un-weighted inner product                          !            ______
           gamma_gmres(1) = sqrt(glsc2(r_gmres,r_gmres,ntot2))! gamma  = \/(r,r)
         endif   
                                                           
         if(iter.eq.0) then
           div0 = gamma_gmres(1)*norm_fac
           if (param(21).lt.0) tolpss=abs(param(21))*div0
         endif

!        check for lucky convergence
         rnorm = 0.
         if (gamma_gmres(1) .eq. 0.) goto 9000
         temp = 1./gamma_gmres(1)
         call cmult2(v_gmres(1,1),r_gmres,temp,ntot2)! v  = r / gamma

         do j=1,m
            iter = iter+1

            call copy(w_gmres,v_gmres(1,j),ntot2) ! w  = v_j
 
            etime2 = dnekclock()
            if (ifprec) then
              if(param(43).eq.1) then
!                 call uzprec(z_gmres(1,j),w_gmres,h1,h2,intype,wp)
                 call eprec2_new(z_gmres(1,j),w_gmres)
              else
                call blank(str,132)
                str = "Multigrid not implented for 
     $ Jacobi Davidson. param(43)= $"
                call exitti(str,int(param(43)))
!                 call hsmg_solve(z_gmres(1,j),w_gmres)    ! z  = (M^-1)w
              endif
            else
              call copy(z_gmres(1,j),w_gmres,ntot2)
            endif

            etime_p = etime_p + dnekclock()-etime2

!           r = (M^-1)w      
            call copy(r_gmres,z_gmres(1,j),ntot2)

!           r = U*(M^-1)*w
            call jd_ortho(r_gmres,u0,wp,1,ifwgtortho)

!           w = A*U*(M^-1)w    
            call jd_Ecorr(w_gmres,r_gmres,h1,h2,h2inv,intype,theta)

!           w = (U^T)*A*U*(M^-1)*w
            call jd_ortho(w_gmres,u0,wp,0,ifwgtortho)

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

!                      ______
!           alpha =  \/(Bw,w) 
            if (ifwgt) then
              alpha = sqrt(glsc3(w_gmres,w_gmres,bm2,ntot2))    
            else
              alpha = sqrt(glsc2(w_gmres,w_gmres,ntot2))        
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
   66       format(i5,1p4e12.5,i8,' JD Correction')

#ifndef FIXITER
            if (rnorm .lt. tolpss) goto 900  !converged
#else
            if (iter.gt.param(151)-1) goto 900
#endif
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v_gmres(1,j+1),w_gmres,temp,ntot2) ! v = w / alpha
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
            call add2s2(x_gmres,z_gmres(1,i),c_gmres(i),ntot2) 
                       ! x = x + Z*c
         enddo

      enddo
 9000 continue

      call copy(res,x_gmres,ntot2)


      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep,' JD-Corr. gmres  ', 
     &                            iter,rnorm,div0,tolpss,etime_p,etime1

 9999 format(i11,a,I6,1p5e13.4)

      return
      end

c-----------------------------------------------------------------------

      subroutine jd_ortho(res,u0,wk,iflg,ifwgt)

!     iflg=0: Left  Ortho: res = (I - (B^T)*(u0*u0^T))*res        
!     iflg=1: Right Ortho: res = (I - (u0*u0^T)*B)*res        

!     Orthogonalize with respect vector u0

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'          ! bm2

      real res (lx2,ly2,lz2,lelv)
      real u0  (lx2,ly2,lz2,lelv)   ! Orthogonalizing vector
      real wk  (lx2,ly2,lz2,lelv)

      integer iflg

      integer n2
      real sc
      real glsc2,glsc3

      real vol2          ! volume squared
      logical ifwgt      ! if weighted

      n2    = nx2*ny2*nz2*nelv

!      ifwgt = .true.

      if (ifwgt) then
        vol2 = glsc3(u0,u0,bm2,n2)
      else  
        vol2 = glsc2(u0,u0,n2)
      endif

      if (iflg.eq.0) then
        if (ifwgt) then
          sc = glsc2(u0,res,n2)/vol2      ! sc  = u0'*res/(vol^2)
          call col3(wk,u0,bm2,n2)         ! wk  = BM2*u0
          call add2s2(res,wk,-sc,n2)      ! res = res - sc*BM2*u0
        else  
          sc = glsc2(u0,res,n2)/vol2      ! sc  = u0'*res/(vol^2) 
          call add2s2(res,u0,-sc,n2)      ! res = res - sc*u0
        endif  
      elseif (iflg.eq.1) then
        if (ifwgt) then
          sc = glsc3(u0,res,bm2,n2)/vol2  ! sc  = u0'*BM2*res/(vol^2)
          call add2s2(res,u0,-sc,n2)      ! res = res - sc*u0 
        else  
          sc = glsc2(u0,res,n2)/vol2      ! sc  = u0'*res/(vol^2)
          call add2s2(res,u0,-sc,n2)      ! res = res - sc*u0
        endif  
      else
         call exitti('jd_ortho: unknown iflg: $',iflg)
      endif

      return
      end
c------------------------------------------------------------------------


      subroutine jd_Eop(w,v,h1,h2,h2inv,intype) 

      implicit none

      include 'SIZE'
      include 'TSTEP'   ! dt
      include 'MASS'    ! BM2inv

      real w(1)
      real v(1)
      real h1(1)
      real h2(1)
      real h2inv(1)
      integer intype

      integer nt2
      real mu

      nt2 = nx2*ny2*nz2*nelv

!     w = Ev      
!      call cdabdtp(w,v,h1,h2,h2inv,intype)
      call fm_cdabdtp(w,v,h1,h2,h2inv,intype)

!      call col2(w,bm2inv,nt2)

!     w = (I + \mu*E)v
      mu = -dt
      call add2s1(w,v,mu,nt2)


      return
      end subroutine jd_Eop
!---------------------------------------------------------------------- 

      subroutine jd_Ecorr(w,v,h1,h2,h2inv,intype,theta) 

      implicit none

      include 'SIZE'
      include 'MASS'          ! BM2

      real w(1)
      real v(1)
      real h1(1)
      real h2(1)
      real h2inv(1)
      integer intype

      integer i,nt2
      real theta

      nt2 = nx2*ny2*nz2*nelv

      call jd_Eop(w,v,h1,h2,h2inv,intype)

      call add2s2(w,v,-theta,nt2)


      return
      end subroutine jd_Ecorr
!---------------------------------------------------------------------- 

