!---------------------------------------------------------------------- 
      subroutine uzawa_gmres_lpr(res,h1,h2,h2inv,intype,iter)

c     Solve the pressure equation by right-preconditioned
c     And also with left preconditioned
c     GMRES iteration.
c     intype =  0  (steady)
c     intype =  1  (explicit)
c     intype = -1  (implicit)

!     We solve:   L*(U^T)*E*U*(M^-1)(Mp) = L*f,
!              => L*(U^T)*E*U*(M^-1)q    = L*f,
!              => p = (M^-1)q.
!     Where,      U = (I - p0*(p0^T)*B),
!     is the Subspace without the constant pressure field.
!
!     ortho_left(r)  = (U^T)*r
!     ortho_right(r) = U*r

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'GMRES'

      real divex
      common  /ctolpr/ divex

      logical          ifprint
      common  /cprint/ ifprint

      real             res  (lx2*ly2*lz2*lelv)
      real             h1   (lx1,ly1,lz1,lelv)
      real             h2   (lx1,ly1,lz1,lelv)
      real             h2inv(lx1,ly1,lz1,lelv)

      real wp
      common /scrmg/    wp (lx2,ly2,lz2,lelv)

      real wk1,wk2
      common /ctmp0/   wk1(lgmres),wk2(lgmres)

      real y
      common /cgmres1/ y(lgmres)

      real lpr(lx2,ly2,lz2,lelv)

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
      integer e,i,k,iter
      real etime2,etime_p,ratio,rnorm
      integer nxyz2
      real localvol
      real vlsum

      logical ifweighted      ! if weighted inner products for Gram-Schmidt
      logical iftwopass       ! if 2-pass Gram-Schmidt
      logical iforthowt       ! if weighted null space orthogonalization
      logical iflpr           ! if left preconditioner

      nxyz2  = lx2*ly2*lx2
      ntot2  = lx2*ly2*lz2*nelv

      ifweighted = .false.
      iftwopass  = .false.
      iforthowt  = .true.
      iflpr      = .true.

!     Orthogonalize w.r.t constant vector      
      if (iforthowt) then
        call ortho_left(res)
      else
        call ortho_new(res)
      endif     

!     Create Left preconditioner
      if (iflpr) then
        do e=1,nelv
          localvol = vlsum(bm2(1,1,1,e),nxyz2)
          call cfill(lpr(1,1,1,e),localvol,nxyz2)
        enddo
      else
        call rone(lpr,ntot2)
      endif  

      if (iflpr) call col2(res,lpr,ntot2) 

      if(.not.iflag) then
         iflag=.true.
         call uzawa_gmres_split0(ml_gmres,mu_gmres,bm2,bm2inv,
     $                           lx2*ly2*lz2*nelv)

!         call rone(ml_gmres,ntot2)
!         call rone(mu_gmres,ntot2) 
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
      iconv = 0
      call rzero(x_gmres,ntot2)


      do while(iconv.eq.0.and.iter.lt.1000)            ! prabal

         if(iter.eq.0) then
                                                        !      -1
            call col3(r_gmres,ml_gmres,res,ntot2)       ! r = L  res
         else
            !update residual
            call copy(r_gmres,res,ntot2)                      ! r = res

!           U*x
            if (iforthowt) then
              call ortho_right(x_gmres)
            else
              call ortho_new(x_gmres)
            endif     

            call cdabdtp(w_gmres,x_gmres,h1,h2,h2inv,intype)  ! w = A x

!           (U^T)*A*U*x
            if (iforthowt) then
              call ortho_left(w_gmres)
            else
              call ortho_new(w_gmres)
            endif     

!           Left preconditioner: L*(U^T)*A*U*x
            if (iflpr) call col2(w_gmres,lpr,ntot2) 
           
            call add2s2(r_gmres,w_gmres,-1.,ntot2)            ! r = r - w

                                                              !      -1
            call col2(r_gmres,ml_gmres,ntot2)                 ! r = L   r
         endif

         if (ifweighted) then
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
                                                           !            j

 
            etime2 = dnekclock()
            if(param(43).eq.1) then
!               call uzprec(z_gmres(1,j),w_gmres,h1,h2,intype,wp)
               call eprec2_new(z_gmres(1,j),w_gmres)
!               call ortho(z_gmres(1,j))   ! done in uzprec
            else
               call hsmg_solve(z_gmres(1,j),w_gmres)    ! z  = (M^-1)w
            endif 
            etime_p = etime_p + dnekclock()-etime2

!           z = U*(M^-1)*w             ! Remove constant vector
            if (iforthowt) then
              call ortho_right(z_gmres(1,j))
            else
              call ortho_new(z_gmres(1,j))
            endif     
     
            call cdabdtp(w_gmres,z_gmres(1,j),    ! w = A*U*(M^-1)w
     $                   h1,h2,h2inv,intype)     

!           w = (U^T)*A*U*(M^-1)*w
            if (iforthowt) then
              call ortho_right(w_gmres)
            else
              call ortho_new(w_gmres)
            endif    

!           Left preconditioner: L*(U^T)*A*U*x
            if (iflpr) call col2(w_gmres,lpr,ntot2) 

                                                  !      -1
            call col2(w_gmres,ml_gmres,ntot2)     ! w = L   w


c           Gram-Schmidt, 1st pass:
            do i=1,j
               if (ifweighted) then
                 h_gmres(i,j)=vlsc3(w_gmres,v_gmres(1,i),bm2,ntot2) ! h = (Bw,V )
               else
                 h_gmres(i,j)=vlsc2(w_gmres,v_gmres(1,i),ntot2)     ! h = (w,V )
               endif
            enddo                                             
            call gop(h_gmres(1,j),wk1,'+  ',j)          ! sum over P procs

            do i=1,j
               call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),ntot2) ! w = w - Vh
            enddo


!           gram-Schmidt, 2nd pass:
            if (iftwopass) then
              do i=1,j
                if (ifweighted) then
                  wk1(i)=vlsc3(w_gmres,v_gmres(1,i),bm2,ntot2) ! g = (Bw,V )
                else
                  wk1(i)=vlsc2(w_gmres,v_gmres(1,i),ntot2)     ! g = (w,V )
                endif
              enddo                                             
              call gop(wk1,wk2,'+  ',j)          ! sum over P procs

              do i=1,j
                call add2s2(w_gmres,v_gmres(1,i),-wk1(i),ntot2) ! w = w - Vg
                h_gmres(i,j) = h_gmres(i,j) + wk1(i)            ! h = h + g 
              enddo
            endif

            !apply Givens rotations to new column
            do i=1,j-1
               temp = h_gmres(i,j)                   
               h_gmres(i  ,j)=  c_gmres(i)*temp 
     $                        + s_gmres(i)*h_gmres(i+1,j)  
               h_gmres(i+1,j)= -s_gmres(i)*temp 
     $                        + c_gmres(i)*h_gmres(i+1,j)
            enddo

!                      ______
!           alpha =  \/(Bw,w) 
            if (ifweighted) then
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
   66       format(i5,1p4e12.5,i8,' Divergence')

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
c        if(iconv.eq.1) call dbg_write(x,lx2,ly2,lz2,nelv,'esol',3)
!         call ortho_right(x_gmres)              ! prabal
      enddo
 9000 continue
c
      divex = rnorm
c     iter = iter - 1
      call copy(res,x_gmres,ntot2)

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep,'  U-PRES gmres  ', 
     &                            iter,divex,div0,tolpss,etime_p,etime1
c     call flush_hack
 9999 format(i11,a,I6,1p5e13.4)

      return
      end

c-----------------------------------------------------------------------
      subroutine uzawa_gmres_wt(res,h1,h2,h2inv,intype,iter)

c     Solve the pressure equation by right-preconditioned 
c     GMRES iteration.
c     intype =  0  (steady)
c     intype =  1  (explicit)
c     intype = -1  (implicit)

!     We solve:   (U^T)*E*U*(M^-1)(Mp) = f,
!              => (U^T)*E*U*(M^-1)q    = f,
!              => p = (M^-1)q.
!     Where,      U = (I - p0*(p0^T)*B),
!     is the Subspace without the constant pressure field.
!
!     ortho_left(r)  = (U^T)*r
!     ortho_right(r) = U*r

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'GMRES'

      real divex
      common  /ctolpr/ divex

      logical          ifprint
      common  /cprint/ ifprint

      real             res  (lx2*ly2*lz2*lelv)
      real             h1   (lx1,ly1,lz1,lelv)
      real             h2   (lx1,ly1,lz1,lelv)
      real             h2inv(lx1,ly1,lz1,lelv)

      real wp
      common /scrmg/    wp (lx2,ly2,lz2,lelv)

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

      logical ifweighted
      logical iftwopass
      logical iforthowt

      ntot2  = lx2*ly2*lz2*nelv

      ifweighted = .true.
      iftwopass  = .false.
      iforthowt  = .true.

!     Orthogonalize w.r.t constant vector      
      if (iforthowt) then
        call ortho_left(res)
      else
        call ortho_new(res)
      endif     

      if(.not.iflag) then
         iflag=.true.
         call uzawa_gmres_split0(ml_gmres,mu_gmres,bm2,bm2inv,
     $                           lx2*ly2*lz2*nelv)

!         call rone(ml_gmres,ntot2)
!         call rone(mu_gmres,ntot2) 
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
      iconv = 0
      call rzero(x_gmres,ntot2)

      do while(iconv.eq.0.and.iter.lt.2000)            ! prabal

         if(iter.eq.0) then
            call col3(r_gmres,ml_gmres,res,ntot2)       ! r = (L^-1)*res
         else
            !update residual
            call copy(r_gmres,res,ntot2)                ! r = res

!           U*x
            if (iforthowt) then
              call ortho_right(x_gmres)
            else
              call ortho_new(x_gmres)
            endif     

            call cdabdtp(w_gmres,x_gmres,h1,h2,h2inv,intype)  ! w = A x

!           (U^T)*A*U*x
            if (iforthowt) then
              call ortho_left(w_gmres)
            else
              call ortho_new(w_gmres)
            endif     

            call add2s2(r_gmres,w_gmres,-1.,ntot2)    ! r = r - w

            call col2(r_gmres,ml_gmres,ntot2)         ! r = (L^-1)r
         endif

         if (ifweighted) then
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

         !check for lucky convergence
         rnorm = 0.
         if(gamma_gmres(1) .eq. 0.) goto 9000
         temp = 1./gamma_gmres(1)
         call cmult2(v_gmres(1,1),r_gmres,temp,ntot2)! v  = r / gamma

         do j=1,m
            iter = iter+1

            call col3(w_gmres,mu_gmres,v_gmres(1,j),ntot2) ! w  = (U^-1)v_j

 
            etime2 = dnekclock()
            if(param(43).eq.1) then
!               call uzprec(z_gmres(1,j),w_gmres,h1,h2,intype,wp)
               call eprec2_new(z_gmres(1,j),w_gmres)
!               call ortho(z_gmres(1,j))   ! done in uzprec
            else
               call hsmg_solve(z_gmres(1,j),w_gmres)    ! z  = (M^-1)w
            endif 
            etime_p = etime_p + dnekclock()-etime2

!           r = (M^-1)w      
            call copy(r_gmres,z_gmres(1,j),ntot2)

!           r = U*(M^-1)*w             ! Remove constant vector
            if (iforthowt) then
              call ortho_right(r_gmres)
            else
              call ortho_new(r_gmres)
            endif     

!           w = A*U*(M^-1)w    
            call cdabdtp(w_gmres,r_gmres,h1,h2,h2inv,intype)

!           w = (U^T)*A*U*(M^-1)*w
            if (iforthowt) then
              call ortho_right(w_gmres)
            else
              call ortho_new(w_gmres)
            endif    

            call col2(w_gmres,ml_gmres,ntot2)     ! w = (L^-1)w

c           Gram-Schmidt, 1st pass:
            do i=1,j
               if (ifweighted) then
                 h_gmres(i,j)=vlsc3(w_gmres,v_gmres(1,i),bm2,ntot2) ! h = (Bw,V )
               else
                 h_gmres(i,j)=vlsc2(w_gmres,v_gmres(1,i),ntot2)     ! h = (w,V )
               endif
            enddo                                             
            call gop(h_gmres(1,j),wk1,'+  ',j)          ! sum over P procs

            do i=1,j
               call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),ntot2) ! w = w - Vh
            enddo


!           Gram-Schmidt, 2nd pass:
            if (iftwopass) then
              do i=1,j
                if (ifweighted) then
                  wk1(i)=vlsc3(w_gmres,v_gmres(1,i),bm2,ntot2) ! g = (Bw,V )
                else
                  wk1(i)=vlsc2(w_gmres,v_gmres(1,i),ntot2)     ! g = (w,V )
                endif
              enddo                                             
              call gop(wk1,wk2,'+  ',j)          ! sum over P procs

              do i=1,j
                call add2s2(w_gmres,v_gmres(1,i),-wk1(i),ntot2) ! w = w - Vg
                h_gmres(i,j) = h_gmres(i,j) + wk1(i)            ! h = h + g 
              enddo
            endif

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
            if (ifweighted) then
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
   66       format(i5,1p4e12.5,i8,' Divergence')

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
!         call ortho_right(x_gmres)              ! prabal
      enddo
 9000 continue
c
      divex = rnorm
c     iter = iter - 1
      call copy(res,x_gmres,ntot2)

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep,'  U-PRES gmres  ', 
     &                            iter,divex,div0,tolpss,etime_p,etime1
c     call flush_hack
 9999 format(i11,a,I6,1p5e13.4)

      return
      end

c-----------------------------------------------------------------------

      subroutine uzawa_gmres_std(res,h1,h2,h2inv,intype,iter)

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

!     prabal. I've removed ortho from earlier calls and call it here at
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
               call ortho(z_gmres(1,j))
            else                                        !       -1
               call hsmg_solve(z_gmres(1,j),w_gmres)    ! z  = M   w
c              call copy(z_gmres(1,j),w_gmres,ntot2)    ! z  = M   w
            endif     
            etime_p = etime_p + dnekclock()-etime2
     
            call cdabdtp(w_gmres,z_gmres(1,j),    ! w = A z
     $                   h1,h2,h2inv,intype)      !        j
     
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
c
c     DIAGNOSTICS
c      call copy   (w,x,ntot2)
       call ortho  (w_gmres) ! Orthogonalize wrt null space, if present
c      call copy(r,res,ntot2) !r = res
c      call cdabdtp(r,w,h1,h2,h2inv,intype)  ! r = A w
c      do i=1,ntot2
c         r(i) = res(i) - r(i)               ! r = res - r
c      enddo
c      call uzawa_gmres_temp(r,bm2inv,ntot2)
c                                               !            ______
c      gamma(1) = sqrt(glsc2(r,r,ntot2)/volvm2) ! gamma  = \/ (r,r) 
c                                               !      1
c      print *, 'GMRES end resid:',gamma(1)
c     END DIAGNOSTICS
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
      subroutine uzawa_gmres_testing(res,h1,h2,h2inv,intype,iter,igmres)

c     Solve the pressure equation by right-preconditioned 
c     GMRES iteration.
c     intype =  0  (steady)
c     intype =  1  (explicit)
c     intype = -1  (implicit)

!     We solve:   (U^T)*E*U*(M^-1)*(M*p) = f,
!              => (U^T)*E*U*(M^-1)*q     = f,
!              => p = (M^-1)*q.
!     Where,      U = (I - p0*(p0^T)*B),
!     is the Subspace without the constant pressure field.
!
!     ortho_left(r)  = (U^T)*r
!     ortho_right(r) = U*r

!     IGMRES:    1) (U^T)*E*(M^-1) *M*U*p               = (U^T)*f,
!                2) (U^T)*E*U*(M^-1) *(M*p)             = (U^T)*f,
!                3) (U^T)*E*(M^-1)*U *(M*p)             = (U^T)*f,

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'GMRES'

      real divex
      common  /ctolpr/ divex

      logical          ifprint
      common  /cprint/ ifprint

      real             res  (lx2*ly2*lz2*lelv)
      real             h1   (lx1,ly1,lz1,lelv)
      real             h2   (lx1,ly1,lz1,lelv)
      real             h2inv(lx1,ly1,lz1,lelv)

      real wp
      common /scrmg/    wp (lx2,ly2,lz2,lelv)

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

      logical ifweighted
      logical iftwopass
      logical iforthowt

      integer igmres

      ntot2  = lx2*ly2*lz2*nelv

      ifweighted = .false.
      iftwopass  = .false.
      iforthowt  = .false.

      if (ifweighted) iforthowt = .true.

!     Orthogonalize w.r.t constant vector
 
      if (igmres.gt.0.and.igmres.le.3) then
        if (iforthowt) then
          call ortho_left(res)
        else
          call ortho_new(res)
        endif
      endif 

      if(.not.iflag) then
         iflag=.true.
         call uzawa_gmres_split0(ml_gmres,mu_gmres,bm2,bm2inv,
     $                           lx2*ly2*lz2*nelv)

!         call rone(ml_gmres,ntot2)
!         call rone(mu_gmres,ntot2) 
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
      iconv = 0
      call rzero(x_gmres,ntot2)

      do while(iconv.eq.0.and.iter.lt.2000)            ! prabal

         if(iter.eq.0) then
            call col3(r_gmres,ml_gmres,res,ntot2)       ! r = (L^-1)*res
         else
            !update residual
            call copy(r_gmres,res,ntot2)                ! r = res

!           igmres==1: x
!           igmres==2: Ux
            if (igmres.eq.2.or.igmres.eq.3) then
              if (iforthowt) then
                call ortho_right(x_gmres)
              else
                call ortho_new(x_gmres)
              endif
            endif 

            call cdabdtp(w_gmres,x_gmres,h1,h2,h2inv,intype)  ! w = A x

!           (U^T)*A*U*x
            if (iforthowt) then
              call ortho_left(w_gmres)
            else
              call ortho_new(w_gmres)
            endif     

            call add2s2(r_gmres,w_gmres,-1.,ntot2)    ! r = r - w

            call col2(r_gmres,ml_gmres,ntot2)         ! r = (L^-1)r
         endif

         if (ifweighted) then
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

         !check for lucky convergence
         rnorm = 0.
         if(gamma_gmres(1) .eq. 0.) goto 9000
         temp = 1./gamma_gmres(1)
         call cmult2(v_gmres(1,1),r_gmres,temp,ntot2)! v  = r / gamma

         do j=1,m
            iter = iter+1

            call col3(w_gmres,mu_gmres,v_gmres(1,j),ntot2) ! w  = (U^-1)v_j

!           w = U*w      
            if (igmres.eq.3) then
              if (iforthowt) then
                call ortho_right(w_gmres)
              else
                call ortho_new(w_gmres)
              endif
            endif 

            etime2 = dnekclock()
            if(param(43).eq.1) then
!               call uzprec(z_gmres(1,j),w_gmres,h1,h2,intype,wp)
               call eprec2_new(z_gmres(1,j),w_gmres)
            else
               call hsmg_solve(z_gmres(1,j),w_gmres)    ! z  = (M^-1)w
            endif 
            etime_p = etime_p + dnekclock()-etime2

!           r = (M^-1)w      
            call copy(r_gmres,z_gmres(1,j),ntot2)
            if (igmres.eq.2) then
!             r = U*(M^-1)*w 
              if (iforthowt) then
                call ortho_right(r_gmres)
              else
                call ortho_new(r_gmres)
              endif
            endif 

!           w = A*U*(M^-1)w    
            call cdabdtp(w_gmres,r_gmres,h1,h2,h2inv,intype)

!           w = (U^T)*A*U*(M^-1)*w
            if (igmres.eq.1.or.igmres.eq.2.or.igmres.eq.3) then
              if (iforthowt) then
                call ortho_right(w_gmres)
              else
                call ortho_new(w_gmres)
              endif
            endif 

            call col2(w_gmres,ml_gmres,ntot2)     ! w = (L^-1)w

c           Gram-Schmidt, 1st pass:
            do i=1,j
               if (ifweighted) then
                 h_gmres(i,j)=vlsc3(w_gmres,v_gmres(1,i),bm2,ntot2) ! h = (Bw,V )
               else
                 h_gmres(i,j)=vlsc2(w_gmres,v_gmres(1,i),ntot2)     ! h = (w,V )
               endif
            enddo                                             
            call gop(h_gmres(1,j),wk1,'+  ',j)          ! sum over P procs

            do i=1,j
               call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),ntot2) ! w = w - Vh
            enddo


!           Gram-Schmidt, 2nd pass:
            if (iftwopass) then
              do i=1,j
                if (ifweighted) then
                  wk1(i)=vlsc3(w_gmres,v_gmres(1,i),bm2,ntot2) ! g = (Bw,V )
                else
                  wk1(i)=vlsc2(w_gmres,v_gmres(1,i),ntot2)     ! g = (w,V )
                endif
              enddo                                             
              call gop(wk1,wk2,'+  ',j)          ! sum over P procs

              do i=1,j
                call add2s2(w_gmres,v_gmres(1,i),-wk1(i),ntot2) ! w = w - Vg
                h_gmres(i,j) = h_gmres(i,j) + wk1(i)            ! h = h + g 
              enddo
            endif

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
            if (ifweighted) then
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
   66       format(i5,1p4e12.5,i8,' Divergence')

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
!         call ortho_right(x_gmres)              ! prabal

!        x = U*x
         if (igmres.eq.2) then   
           if (iforthowt) then
             call ortho_right(x_gmres)
           else
             call ortho_new(x_gmres)
           endif
         endif   
      enddo
 9000 continue
c
      divex = rnorm
c     iter = iter - 1
      call copy(res,x_gmres,ntot2)

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep,'  U-PRES gmres  ', 
     &                            iter,divex,div0,tolpss,etime_p,etime1
c     call flush_hack
 9999 format(i11,a,I6,1p5e13.4)

      return
      end subroutine uzawa_gmres_testing

c-----------------------------------------------------------------------




