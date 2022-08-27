c-----------------------------------------------------------------------

      subroutine lsm_march_phi2

      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'MASS'
      include 'SOLN'
      include 'DEALIAS'

      include 'LSM'

      include 'TEST'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real ta1(lv)
      real ta2(lv)
      real ta3(lv)
      real ta4(lv)
      common /scruz/ ta1,ta2,ta3,ta4

      integer i,j
      integer n

      real eps
      real eps_ss

      logical ifdistort
      character nam*3

!      ifdistort = .false.
!      call phi0_ellipse(t,ifdistort)

      n = lx1*ly1*lz1*nelv

      ifield = 3

      if (nio.eq.0) write(6,*) 'Distance Correction on ifield: ', ifield

      lsm_noncon  = .false.         ! no convection.
      lsm_dealias = .true.          ! use dealiasing for nonlinear term
      lsm_ifcf    = .false.
      lsm_fil     = .true.

      if (lsm_noncon) lsm_dealias = .false.
      if (.not.lsm_dealias) lsm_ifcf = .false.

      call eval_normals(lsm_nx,lsm_ny,lsm_nz,lsm_nnorm,t)

      eps_ss = 1.0e-2
      call SmoothSign(lsm_ssign,t,eps_ss,n)      ! Smooth Sign! SmoothSign is the forcing term
      nam = 'SIGN'
!      call lsm_regularize_field(lsm_ssign,1.0e-4,1.0e-4,nam)
      call copy(lsm_forc,lsm_ssign,n)

      if (lsm_ifcf) then
        call set_convect_sc(vxd,vyd,vzd,lsm_nx,lsm_ny,lsm_nz,lsm_ssign)
      endif  
!      set_convect_sc(cr,cs,ct,ux,uy,uz,sc)

      if (.not.lsm_ifcf) then
        call col2(lsm_nx,lsm_ssign,n)  
        call col2(lsm_ny,lsm_ssign,n)  
        if (ndim.eq.3) call col2(lsm_nz,lsm_ssign,n)
      endif

      call opcopy(vx,vy,vz,lsm_nx,lsm_ny,lsm_nz)

      ifdistort = .true.
      call phi0_ellipse(t,ifdistort)

      call copy(t(1,1,1,1,2),lsm_ssign,n)
      call outpost2(vx,vy,vz,pr,t,2,'  ')

      do i=1,nsteps

         istep = 100+i
!        call lsm_rk2(t,SSign)
!        call lsm_rk4(t,SSign)
        call lsm_ssprk3_v2(t)
!        call lsm_ExpEuler(t,SSign)

        if (mod(i,iostep).eq.0) then
!          call eval_normals(lsm_nx,lsm_ny,lsm_nz,lsm_nnorm,t)
!          call SmoothSign(lsm_ssign,t,eps_ss,n)      ! Smooth Sign! SmoothSign is the forcing term
!          nam = 'SIGN'           
!          call lsm_regularize_field(lsm_ssign,1.0e-2,1.0e-2,nam)
!          call copy(lsm_forc,lsm_ssign,n)

!          if (lsm_ifcf)  
!     $      call set_convect_sc(vxd,vyd,vzd,vx,vy,vz,lsm_ssign)

!          if (.not.lsm_ifcf) then
!            call col2(lsm_nx,lsm_ssign,n)  
!            call col2(lsm_ny,lsm_ssign,n)  
!            if (ndim.eq.3) call col2(lsm_nz,lsm_ssign,n)
!          endif
!          call opcopy(vx,vy,vz,lsm_nx,lsm_ny,lsm_nz)

!          call copy(t(1,1,1,1,2),lsm_ssign,n)
        endif
        if (mod(i,iostep).eq.0) then
          call outpost2(vx,vy,vz,pr,t,2,'  ')
        endif  
      enddo  

      return
      end subroutine lsm_march_phi2
!---------------------------------------------------------------------- 
      subroutine lsm_ssprk3_v2(d0)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'TSTEP'
      include 'SOLN'
      include 'INPUT'
      include 'DEALIAS'

      include 'LSM'

      include 'TEST'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real d0(lv)

      real ta1(lv)
      real ta2(lv)
      real ta3(lv)
      real ta4(lv)
      common /scruz/ ta1,ta2,ta3,ta4

      real fil(lv) 
      common /scrcg/ fil

      real tmpsol(lv,3)
      common /scrsf/ tmpsol

      real dt2,dti

      real h1(lv)
      real h2(lv)

      real gdn

      integer i,n
      logical iffil,ifhmh,ifcf,ifss,nocon

      real distn
      real op_glsc2_wt,glmax,glamax
      real gmax,gmaxdt

      real s1,s2

      integer imsh,maxit,isd
      real tli

      iffil = lsm_fil
      ifhmh = .false.
      ifcf  = lsm_ifcf
      ifss  = .false.    ! multiply by smooth sign
      if (ifcf) ifss = .false.
      nocon = lsm_noncon

      n = lx1*ly1*lz1*nelv

      dt2 = 1.0*dt
      dti = 1.0/dt2

      ifield = 3

!      istep = 20

      if (ifhmh) then
        call cfill(h2,dti,n)
        call cfill(h1,cpfld(ifield,1),n)

        imsh = 1
        maxit = 2000
        isd  = 1
        tli  = 1.0e-10

!        call copy(t(1,1,1,1,2),h1,n)
      endif  

!     1st step
!--------------------------------------------------       
      if (lsm_dealias) then
        if (ifcf) then
          call convect_new(ta1,d0,.false.,vxd,vyd,vzd,.true.)
        else
          call convect_new(ta1,d0,.false.,vx,vy,vz,.false.)
        endif
        call invcol2(ta1,bm1,n)
      else
        call gradm1(ta2,ta3,ta4,d0)
        call dsavg(ta2)
        call dsavg(ta3)
        if (ndim.eq.3) call dsavg(ta4)
        if (nocon) then
          do i=1,n
            s1 = ta2(i)**2 + ta3(3)**2
            if (ndim.eq.3) s1 = s1 + ta4(i)**2
            ta2(i) = sqrt(s1)
          enddo
          call dealias_uv(ta1,ta2,lsm_ssign)
          call invcol2(ta1,bm1,n) 
        else  
          if (ndim.eq.3) then
            call vdot3(ta1,vx,vy,vz,ta2,ta3,ta4,n)
          else
            call vdot2(ta1,vx,vy,ta2,ta3,n)
          endif
        endif   
      endif
      if (ifss) then
         call col2(ta1,lsm_ssign,n)
      endif   
      call add2s1(ta1,lsm_forc,-1.0,n)
      if (iffil) then
        call hpf_dist(fil,d0)
        call add2(ta1,fil,n)           ! already has negative sign
      endif
      call col2(ta1,bm1,n)
      if (ifhmh) then
!       Euler Solution
        call col3(ta3,d0,bm1,n)
        call add2s2(ta1,ta3,dti,n)
        call hmholtz('CORR',ta2,ta1,h1,h2,tmask,tmult,
     $                imsh,tli,maxit,isd)

      else  
        call invertB(ta2,ta1)
        call add2s1(ta2,d0,dt2,n)
      endif  
      call copy(tmpsol(1,1),ta2,n)             ! 1st soln: f1 = Euler(f0)

!     2nd step
!--------------------------------------------------       
      if (lsm_dealias) then
        if (ifcf) then
          call convect_new(ta1,tmpsol(1,1),.false.,vxd,vyd,vzd,.true.)
        else
          call convect_new(ta1,tmpsol(1,1),.false.,vx,vy,vz,.false.)
        endif
        call invcol2(ta1,bm1,n)
      else
        call gradm1(ta2,ta3,ta4,tmpsol(1,1))
        call dsavg(ta2)
        call dsavg(ta3)
        if (ndim.eq.3) call dsavg(ta4)

        if (nocon) then
          do i=1,n
            s1 = ta2(i)**2 + ta3(3)**2
            if (ndim.eq.3) s1 = s1 + ta4(i)**2
            ta2(i) = sqrt(s1)
          enddo
          call dealias_uv(ta1,ta2,lsm_ssign)
          call invcol2(ta1,bm1,n) 
        else  
          if (ndim.eq.3) then
            call vdot3(ta1,vx,vy,vz,ta2,ta3,ta4,n)
          else
            call vdot2(ta1,vx,vy,ta2,ta3,n)
          endif
        endif  
      endif  
      if (ifss) then
        call col2(ta1,lsm_ssign,n)
      endif   
    
      call add2s1(ta1,lsm_forc,-1.0,n)
      if (iffil) then
        call hpf_dist(fil,tmpsol(1,1))
        call add2(ta1,fil,n)           ! already has negative sign
      endif
      call col2(ta1,bm1,n)
      if (ifhmh) then
!       Euler Solution            
        call col3(ta3,tmpsol(1,1),bm1,n)
        call add2s2(ta1,ta3,dti,n)
        call hmholtz('CORR',ta2,ta1,h1,h2,tmask,tmult,
     $                imsh,tli,maxit,isd)

      else  
        call invertB(ta2,ta1)
        call add2s1(ta2,tmpsol(1,1),dt2,n)   ! Euler solution
      endif  

      s2    = 1.0/4.0
      call cmult(ta2,s2,n)
      call copy(tmpsol(1,2),ta2,n)
      s2    = 3.0/4.0
      call add2s2(tmpsol(1,2),d0,s2,n)    ! 2nd soln: f2 = 3/4*f0 + 1/4*Euler(f1)

!     3rd step
!--------------------------------------------------       
      if (lsm_dealias) then 
        if (ifcf) then
          call convect_new(ta1,tmpsol(1,2),.false.,vxd,vyd,vzd,.true.)
        else
          call convect_new(ta1,tmpsol(1,2),.false.,vx,vy,vz,.false.)
        endif
        call invcol2(ta1,bm1,n)
      else
        call gradm1(ta2,ta3,ta4,tmpsol(1,2))
        call dsavg(ta2)
        call dsavg(ta3)
        if (ndim.eq.3) call dsavg(ta4)
        if (nocon) then
          do i=1,n
            s1 = ta2(i)**2 + ta3(3)**2
            if (ndim.eq.3) s1 = s1 + ta4(i)**2
            ta2(i) = sqrt(s1)
          enddo
          call dealias_uv(ta1,ta2,lsm_ssign)
          call invcol2(ta1,bm1,n) 
        else  
          if (ndim.eq.3) then
            call vdot3(ta1,vx,vy,vz,ta2,ta3,ta4,n)
          else
            call vdot2(ta1,vx,vy,ta2,ta3,n)
          endif
        endif  
      endif  
      if (ifss) then
        call col2(ta1,lsm_ssign,n)
      endif   
    
      call add2s1(ta1,lsm_forc,-1.0,n)
      if (iffil) then
        call hpf_dist(fil,tmpsol(1,2))
        call add2(ta1,fil,n)           ! already has negative sign
      endif
      call col2(ta1,bm1,n)
      if (ifhmh) then
!       Euler Solution            
        call col3(ta3,tmpsol(1,2),bm1,n)
        call add2s2(ta1,ta3,dti,n)
        call hmholtz('CORR',ta2,ta1,h1,h2,tmask,tmult,
     $                imsh,tli,maxit,isd)

      else  
        call invertB(ta2,ta1)
        call add2s1(ta2,tmpsol(1,2),dt2,n)   ! Euler solution
      endif  

      s2    = 2.0/3.0
      call cmult(ta2,s2,n)
      call copy(tmpsol(1,3),ta2,n)
      s2    = 1.0/3.0
      call add2s2(tmpsol(1,3),d0,s2,n)    ! Final soln: fn+1 = 1/3*f0 + 2/3*Euler(f2)

      call copy(d0,tmpsol(1,3),n)


      return
      end subroutine lsm_ssprk3_v2
!---------------------------------------------------------------------- 
      subroutine set_convect_sc(cr,cs,ct,ux,uy,uz,sc)
C
C     Put vxd,vyd,vzd into rst form on fine mesh
C
C     For rst form, see eq. (4.8.5) in Deville, Fischer, Mund (2002).
C
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
!      include 'TOTAL'

      integer lxy,ltd
      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)

      real cr(ltd,1),cs(ltd,1),ct(ltd,1)
      real ux(lxy,1),uy(lxy,1),uz(lxy,1)
      real sc(lxy,1)

      real fx,fy,fz
      real ur,us,ut
      real tr,uf
      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , tr(ltd,3),uf(ltd)
      real fsc(ltd)     ! scalar on fine mesh

      integer i,e
      integer ic,nxyz1,nxyzd

      call set_dealias_rx

      nxyz1 = lx1*ly1*lz1
      nxyzd = lxd*lyd*lzd

      ic = 1    ! pointer to vector field C

      do e=1,nelv 

c        Map coarse velocity to fine mesh (C-->F)

         call intp_rstd(fx,ux(1,e),lx1,lxd,if3d,0) ! 0 --> forward
         call intp_rstd(fy,uy(1,e),lx1,lxd,if3d,0) ! 0 --> forward
         if (if3d) call intp_rstd(fz,uz(1,e),lx1,lxd,if3d,0) ! 0 --> forward

!        Interpolate the scalar variable   
         call intp_rstd(fsc,sc(1,e),lx1,lxd,if3d,0) ! 0 --> forward
!        Scale the convecting velocities with the 
!        scalar on the fine mesh            
         call col2(fx,fsc,nxyzd)
         call col2(fy,fsc,nxyzd)
         if (if3d) call col2(fz,fsc,nxyzd)

c        Convert convector F to r-s-t coordinates

         if (if3d) then

           do i=1,nxyzd
              cr(i,e)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)+rx(i,3,e)*fz(i)
              cs(i,e)=rx(i,4,e)*fx(i)+rx(i,5,e)*fy(i)+rx(i,6,e)*fz(i)
              ct(i,e)=rx(i,7,e)*fx(i)+rx(i,8,e)*fy(i)+rx(i,9,e)*fz(i)
           enddo

         else

           do i=1,nxyzd
              cr(i,e)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)
              cs(i,e)=rx(i,3,e)*fx(i)+rx(i,4,e)*fy(i)
           enddo

         endif
      enddo

      return
      end subroutine set_convect_sc
c-----------------------------------------------------------------------

      subroutine eval_normals(px,py,pz,pn,phi)

      implicit none

      include 'SIZE'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real phi(lv)

      real px(lv)
      real py(lv)
      real pz(lv)
      real pn(lv)

      logical ifregg
      real viscL,viscH

      real gn

      integer i,n

      character nam*4

      n = lx1*ly1*lz1*nelv

      ifregg = .true.
      nam    = 'NORM'

      call rzero3(px,py,pz,n)
      call gradm1(px,py,pz,phi)

!     Gradient Norm    
      do i=1,n
        gn = px(i)**2 + py(i)**2
        if (ndim.eq.3) gn = gn + pz(i)**2
        pn(i) = sqrt(gn)    
      enddo 

!     Normalize gradients (= Normals)        
      call invcol2(px,pn,n)  
      call invcol2(py,pn,n)  
      if (ndim.eq.3) call invcol2(pz,pn,n)

      if (ifregg) then
        viscL = 1.0e-0
        viscH = 1.0e-0
        call lsm_regularize_field(px,viscL,viscH,nam)
        call lsm_regularize_field(py,viscL,viscH,nam)
        if (ndim.eq.3) call lsm_regularize_field(pz,viscL,viscH,nam)
      endif
!     Gradient Norm    
      do i=1,n
        gn = px(i)**2 + py(i)**2
        if (ndim.eq.3) gn = gn + pz(i)**2
        pn(i) = sqrt(gn)    
      enddo 

!     Normalize gradients (= Normals)        
      call invcol2(px,pn,n)  
      call invcol2(py,pn,n)  
      if (ndim.eq.3) call invcol2(pz,pn,n)

!!     Trying the new dealiasing version    
!      call col2(lsm_nx,lsm_ssign,n)  
!      call col2(lsm_ny,lsm_ssign,n)  
!      if (ndim.eq.3) call col2(lsm_nz,lsm_ssign,n)
 
!      call opcopy(vx,vy,vz,lsm_nx,lsm_ny,lsm_nz)


      return
      end subroutine eval_normals
!----------------------------------------------------------------------       
      subroutine eval_normals2(px,py,pz,pn,phi)

      implicit none

      include 'SIZE'

      integer lv
      parameter (lv=lx1*ly1*lz1*lelv)

      real phi(lv)

      real px(lv)
      real py(lv)
      real pz(lv)
      real pn(lv)

      logical ifregg
      real viscL,viscH

      real gn

      integer i,n

      character nam*4

      n = lx1*ly1*lz1*nelv

      ifregg = .true.
      nam    = 'NORM'

      call rzero3(px,py,pz,n)
      call gradm1(px,py,pz,phi)
     
      if (ifregg) then
        viscL = 1.0e-4
        viscH = 1.0e-4
        call lsm_regularize_field(px,viscL,viscH,nam)
        call lsm_regularize_field(py,viscL,viscH,nam)
        if (ndim.eq.3) call lsm_regularize_field(pz,viscL,viscH,nam)
      endif
!     Gradient Norm    
      do i=1,n
        gn = px(i)**2 + py(i)**2
        if (ndim.eq.3) gn = gn + pz(i)**2
        pn(i) = sqrt(gn)    
      enddo 

!     Normalize gradients (= Normals)        
      call invcol2(px,pn,n)  
      call invcol2(py,pn,n)  
      if (ndim.eq.3) call invcol2(pz,pn,n)

!     Gradient Norm    
      do i=1,n
        gn = px(i)**2 + py(i)**2
        if (ndim.eq.3) gn = gn + pz(i)**2
        px(i) = px(i)/sqrt(gn)
        py(i) = py(i)/sqrt(gn)
        if (ndim.eq.3) pz(i) = pz(i)/sqrt(gn)
      enddo 

!!     Trying the new dealiasing version    
!      call col2(lsm_nx,lsm_ssign,n)  
!      call col2(lsm_ny,lsm_ssign,n)  
!      if (ndim.eq.3) call col2(lsm_nz,lsm_ssign,n)
 
!      call opcopy(vx,vy,vz,lsm_nx,lsm_ny,lsm_nz)


      return
      end subroutine eval_normals2
!----------------------------------------------------------------------       
      subroutine SmoothSign(ss,phi,eps,n)

      implicit none

      real ss(1),phi(1)
      real eps

      integer i,n


      do i=1,n
        if (phi(i).lt.-eps) then
          ss(i) = -1.0
        elseif(phi(i).gt.eps) then
          ss(i) = 1.0
        else
          ss(i) = phi(i)/(sqrt(phi(i)**2 + eps**2))  
        endif
      enddo


      return
      end subroutine SmoothSign
!---------------------------------------------------------------------- 

      subroutine dealias_uv(cku,u,v)

!     Compute dealiased form:  J^T Bf *Jv .Ju w/ correct Jacobians

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'

      real cku(1),u(1),cx(1),cy(1),v(1)
      logical ifuf,ifcf            ! u and/or c already on fine mesh?

      integer lxy,ltd
      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)

      real fx,fy,fz
      real ur,us,ut
      real tr,uf,wd2,jacm1d
      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , jacm1d(ltd),tr(ltd)
     $             , wd2(ltd),uf(ltd)

      integer e,k
      integer iu,iv,ijc,ick

      integer nxyz1,nxyzv,nxyzd,nxyzu,nxyzj

      real zd,wd
      common /dealias1/ zd(lxd),wd(lxd)

      integer i,j,l

      nxyz1 = lx1*ly1*lz1
      nxyzd = lxd*lyd*lzd

      nxyzu = nxyz1
!      if (ifuf) nxyzu = nxyzd

      nxyzv = nxyz1
!      if (ifcf) nxyzc = nxyzd

      nxyzj = nxyz1


      iu  = 1    ! pointer to scalar field u
      iv  = 1    ! pointer to vector field v
      ijc = 1    ! pointer to scalar JACM1
      ick = 1    ! pointer to scalar cku 

      call zwgl (zd,wd,lxd)  ! zwgl -- NOT zwgll!

      if (if3d) then
        do k=1,lzd
        do j=1,lyd
        do i=1,lxd
           l = (k-1)*lyd*lxd + (j-1)*lxd + i 
           wd2(l) = wd(i)*wd(j)*wd(k)
        enddo
        enddo
        enddo
      else
        do j=1,lyd
        do i=1,lxd
           l = (j-1)*lxd + i 
           wd2(l) = wd(i)*wd(j)
        enddo
        enddo

      endif


      do e=1,nelv

!       Interpolate Convecting Field   
        call intp_rstd(fz,v(iv),lx1,lxd,if3d,0) ! 0 --> forward

!       Interpolate Convected Field   
        call intp_rstd(uf,u(iu),lx1,lxd,if3d,0) ! 0 --> forward

!       Interpolate Jacobian (Probably only needs to be done once) 
        call intp_rstd(jacm1d,jacm1(ijc,1,1,1),lx1,lxd,if3d,0) ! 0 --> forward

        do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
           tr(i) = wd2(i)*jacm1d(i)*uf(i)*fz(i)
        enddo
        call intp_rstd(cku(ick),tr,lx1,lxd,if3d,1) ! Project back to coarse

        iv  = iv  + nxyzv
        iu  = iu  + nxyzu
        ijc = ijc + nxyzj
        ick = ick + nxyz1

      enddo

      return
      end subroutine dealias_uv
!-----------------------------------------------------------------------
