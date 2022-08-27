!----------------------------------------------------------------------
!     Author: Prabal Negi
!     Description: Routines for Generating global basis function.
!                  for surface representation
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      subroutine fs_smooth_meshmv(wx,wy,wz)

      implicit none

      include 'SIZE'
      include 'FS_ALE'
!      include 'SOLN'

      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)

      integer icalld
      save icalld
      data icalld /0/

      integer i

!     Generate the global basis anyway      
      if (icalld.eq.0) then
        call fs_global_basis
        icalld = icalld+1
      endif

!     Zero out tangential component of mesh velocity
      call fs_mvmeshn2(wx,wy,wz)


      if (fs_ifgsm) then

        call fs_gllo_xyz

        call fs_gllo_flds(wx,wy,wz)
        call fs_intp_setup
        call fs_get_localpts          ! SEM -> Global

!       Filtering here requires BOYD transformation to be active.
!       Otherwise the boundary points move.
!       Alternately, one could correct boundary points.      
!!       Filter normal velocities
!        if (fs_iffil) then
!          do i=1,ndim
!            call fs_glfilter(fld_fs(1,1,i))
!          enddo  
!        endif
        call fs_get_globalpts         ! Global -> SEM
        call fs_restore_int(wx,wy,wz,'Norm')
      endif       ! fs_ifgsm 

!     Correction for tangential movement to minimize
!     mesh deformation
      if (fs_iftc) then
        call fs_tang_corr(wx,wy,wz)
      endif

!     make sure fluid does not come off the wall      
      call fs_fixcorners(wx,wy,wz) 

!     Free the handles
      if (fs_ifgsm) then 
        call fgslib_findpts_free(intgh_fs)
        call fgslib_findpts_free(intlh_fs)

        call mntr_log(fs_id,fs_log,'Interface Smoothening Done')
      endif  

      return
      end subroutine fs_smooth_meshmv
!----------------------------------------------------------------------
      subroutine fs_smoothmv(wx,wy,wz)

!     Mesh motion using smooth normals/Tangentials        

      implicit none

      include 'SIZE'
      include 'FS_ALE'

      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)

      integer icalld
      save icalld
      data icalld /0/

      if (icalld.eq.0) then
        call fs_global_basis
        icalld = icalld+1
      endif

      call fs_smcoor_mv(wx,wy,wz)

      call mntr_log(fs_id,fs_log,'Interface Smoothening Done')

      return
      end subroutine fs_smoothmv
!----------------------------------------------------------------------
      subroutine fs_smcoor_mv(wx,wy,wz)

!     Here I approximate the coordinates using global polynomials and
!     calculate the normals and tangential directions, which are
!     globally smooth along the interface.

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'FS_ALE'

!     Temporary arrays
      real wk1,wk2,wk3,wk4,wk5,wk6,wk7
      common /scrns3/ wk1(lx1,ly1,lz1,lelt),    ! Tangengial directions     
     $               wk2(lx1,ly1,lz1,lelt),     ! Tangengial directions
     $               wk3(lx1,ly1,lz1,lelt),     ! Tangengial directions
     $               wk4(lx1,ly1,lz1,lelt),     ! Displacement x 
     $               wk5(lx1,ly1,lz1,lelt),     ! Displacement y
     $               wk6(lx1,ly1,lz1,lelt),     ! Displacement z
     $               wk7(lx1,ly1,lz1,lelt)

      integer ix,iy,ntot
      integer i,e,j1,j2,j3

      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)
      
      real wallvx(lx1*lelt),wallvy(lx1*lelt)

      real xs,ys,ss
      real xr,yr,rr

      real s,p,r,mu,dy,tol

      integer tangcorr

!     Save corner point velocities      
      do i = 1,fs_nsymo
        e  = fs_ie(i)
        j1 = fs_ix(i)
        j2 = fs_iy(i)
        j3 = fs_iz(i)
        wallvx(i) = wx(j1,j2,j3,e)
        wallvy(i) = 0.0
      enddo  

      call opcopy(wk1,wk2,wk3,xm1,ym1,zm1)

      call fs_gllo_xyz
      call fs_gllo_flds(wk1,wk2,wk3)
      call fs_intp_setup
      call fs_get_localpts          ! SEM -> Global

!     We could filter coordinates
!     Or we could filter normals/tangentials      
      if (fs_iffil) then
        do i=1,ndim
          call fs_glfilter(fld_fs(1,1,i))
        enddo  
      endif

      call mxm (dx_fs,lxfs,fld_fs(1,1,1),lxfs,xr_fs,lyfs)   ! dx/dr
      call mxm (dx_fs,lxfs,fld_fs(1,1,2),lxfs,yr_fs,lyfs)   ! dy/dr

      call mxm (fld_fs(1,1,1),lxfs,dyt_fs,lyfs,xs_fs,lyfs)  ! dx/ds
      call mxm (fld_fs(1,1,2),lxfs,dyt_fs,lyfs,ys_fs,lyfs)  ! dy/ds

      do 100 iy=1,lyfs
      do 100 ix=1,lxfs
         xs  = xs_fs(ix,iy)
         ys  = ys_fs(ix,iy)
         xr  = xr_fs(ix,iy)
         yr  = yr_fs(ix,iy)

         ss  = sqrt( xs**2 + ys**2 )
         rr  = sqrt( xr**2 + yr**2 )

         t1x_fs(ix,iy)  =  xs / ss
         t1y_fs(ix,iy)  =  ys / ss
         unx_fs(ix,iy)  =  t1y_fs(ix,iy)
         uny_fs(ix,iy)  = -t1x_fs(ix,iy)

         if (ndim.eq.3) then
           t2x_fs(ix,iy)  =  xr / rr
           t2y_fs(ix,iy)  =  yr / rr
!          Calculate normals for 3D
         endif

  100 continue


      call copy(fld_fs(1,1,1),unx_fs,lxfs*lyfs)
      call copy(fld_fs(1,1,2),uny_fs,lxfs*lyfs)
!!     Remove high wavenumbers      
!      if (fs_iffil) then
!        do i=1,ndim
!          call fs_glfilter(fld_fs(1,1,i))
!        enddo  
!      endif

      call fs_get_globalpts         ! Global -> SEM
      call fs_restore_int(wk1,wk2,wk3,'XYZN')
      ntot = lx1*ly1*lz1*nelv
      call unitvec(wk1,wk2,wk3,ntot)
      if (ndim.eq.2) then
        call vdot2(wk4,wk1,wk2,wx,wy,ntot)      ! dot product with normal dir
        call opcopy(wx,wy,wz,wk1,wk2,wk3)       ! copy normals 
        call col2(wx,wk4,ntot)                  ! normal velocities 
        call col2(wy,wk4,ntot)                  ! normal velocities
      else
        call vdot3(wk4,wk1,wk2,wk3,wx,wy,wz,ntot) 

      endif

!     Tangential correction
      tangcorr = 1            ! which tangential correction
      if (fs_iftc) then
          call copy(fld_fs(1,1,1),t1x_fs,lxfs*lyfs)
          call copy(fld_fs(1,1,2),t1y_fs,lxfs*lyfs)
!!         Remove high wavenumbers      
!          if (fs_iffil) then
!            do i=1,ndim
!              call fs_glfilter(fld_fs(1,1,i))
!            enddo  
!          endif

          call fs_get_globalpts         ! Global -> SEM
          call fs_restore_int(wk1,wk2,wk3,'XYZT')
          call unitvec(wk1,wk2,wk3,ntot)

        if (tangcorr.eq.1) then
          tol = 1.0e-14
          do i=1,ntot
            p  = wy(i,1,1,1)    ! current velocity along y
            s  = wk2(i,1,1,1)   ! tangential velocity along y
            if (abs(s).gt.tol) then
              r  = p/s
            else
              r  = 0.0
            endif  
            wx(i,1,1,1) = wx(i,1,1,1) - r*wk1(i,1,1,1)
            wy(i,1,1,1) = wy(i,1,1,1) - r*wk2(i,1,1,1)
            if (ndim.eq.3) wz(i,1,1,1) = wz(i,1,1,1) 
     $                                  - r*wk3(i,1,1,1)
          enddo
        elseif (tangcorr.eq.2) then
!         Local point displacements      
          call opcopy(wk4,wk5,wk6,xm1,ym1,zm1)
          call opsub2(wk4,wk5,wk6,xm0_fs,ym0_fs,zm0_fs) ! dx,dy,dz
     
          tol = 1.0e-14
          mu  = 0.1
          do i=1,ntot
            dy = wk5(i,1,1,1)   ! displacement along y
            p  = wk2(i,1,1,1)   ! tangential velocity along y
            if (abs(p).gt.tol) then
              s  = -p/abs(p)
            else
              s  = 0.0
            endif 
            wk7(i,1,1,1) = s*mu*dy
          enddo
        endif

      endif       ! fs_iftc 

!     Restore corner point velocities      
      do i=1,fs_nsymo
        e  = fs_ie(i)
        j1 = fs_ix(i)
        j2 = fs_iy(i)
        j3 = fs_iz(i)

        wx(j1,j2,j3,e) = wallvx(i)
        wy(j1,j2,j3,e) = wallvy(i)
      enddo   

!     Free the handles      
      call fgslib_findpts_free(intgh_fs)
      call fgslib_findpts_free(intlh_fs)


      return
      end subroutine fs_smcoor_mv

!---------------------------------------------------------------------- 
      subroutine fs_mvmeshn(ux,uy,uz)
!     Only in 2D for now
!     Project velocities on to Normal directions

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'MASS'

      include 'FS_ALE'

      real ux(lx1,ly1,lz1,lelv)
      real uy(lx1,ly1,lz1,lelv)
      real uz(lx1,ly1,lz1,lelv)

      integer i,n,nface
      integer js1,js2,jf1,jf2,jskip1,jskip2,ifc,e
      integer j1,j2,j3,nxyz
      integer ifld

      real rnor,rtn1

      character cb*3

      real dummy1,dummy2,dummy3
      common /scrsf/ dummy1(lx1,ly1,lz1,lelt),
     $               dummy2(lx1,ly1,lz1,lelt),
     $               dummy3(lx1,ly1,lz1,lelt)

      real wallvx(lx1*lelt),wallvy(lx1*lelt)
      real tol
      integer icalld
      save icalld
      data icalld /0/

      ifld  = 1
      nxyz  = lx1*ly1*lz1
      nface = 2*ndim

      do i = 1,fs_nsymo
        e  = fs_ie(i)
        j1 = fs_ix(i)
        j2 = fs_iy(i)
        j3 = fs_iz(i)
        wallvx(i) = ux(j1,j2,j3,e)
        wallvy(i) = 0.0
      enddo  
        

      do 200 e=1,nelv
      do 200 ifc=1,nface
        cb  = fs_cbc(ifc,e)
        if (cb.eq.'INT') then
          call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,ifc)
          i = 0
          do 220 j2=js2,jf2,jskip2
          do 220 j1=js1,jf1,jskip1
             i = i + 1
!            normal component         
             rnor = ( ux(j1,j2,1,e)*unx(i,1,ifc,e) +
     $                uy(j1,j2,1,e)*uny(i,1,ifc,e) )
!            remove tangential component
             ux(j1,j2,1,e) = rnor*unx(i,1,ifc,e)
             uy(j1,j2,1,e) = rnor*uny(i,1,ifc,e)


  220      continue
        endif                 
  200 continue

      call dsavg(ux)
      call dsavg(uy)

      do i=1,fs_nsymo
        e  = fs_ie(i)
        j1 = fs_ix(i)
        j2 = fs_iy(i)
        j3 = fs_iz(i)

        ux(j1,j2,j3,e) = wallvx(i)
        uy(j1,j2,j3,e) = wallvy(i)
      enddo   


      return
      end subroutine fs_mvmeshn        
!---------------------------------------------------------------------- 

      subroutine fs_mvmeshn2(ux,uy,uz)
!     Only in 2D for now
!     Project velocities on to Normal directions

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'MASS'

      include 'FS_ALE'

      real ux(lx1,ly1,lz1,lelv)
      real uy(lx1,ly1,lz1,lelv)
      real uz(lx1,ly1,lz1,lelv)

      integer i,n,nface
      integer js1,js2,jf1,jf2,jskip1,jskip2,ifc,e
      integer j1,j2,j3,nxyz
      integer ifld

      integer iop
      real norx,nory,norz

      real rnor,rtn1

      character cb*3

!     Temporary arrays
      real wk1,wk2,wk3,wk4,wk5,wk6,wk7
      common /scrns/ wk1(lx1,ly1,lz1,lelt),    ! Normal directions     
     $               wk2(lx1,ly1,lz1,lelt),     ! Normal directions
     $               wk3(lx1,ly1,lz1,lelt),     ! Normal directions
     $               wk4(lx1,ly1,lz1,lelt),     ! Tangential directions 
     $               wk5(lx1,ly1,lz1,lelt),     ! Tangential directions
     $               wk6(lx1,ly1,lz1,lelt),     ! Tangential directions
     $               wk7(lx1,ly1,lz1,lelt)

      real wallvx(lx1*lelt),wallvy(lx1*lelt)
      real tol
      integer icalld
      save icalld
      data icalld /0/

      ifld  = 1
      nxyz  = lx1*ly1*lz1
      n     = nxyz*nelv
      nface = 2*ndim

      do i = 1,fs_nsymo
        e  = fs_ie(i)
        j1 = fs_ix(i)
        j2 = fs_iy(i)
        j3 = fs_iz(i)
        wallvx(i) = ux(j1,j2,j3,e)
        wallvy(i) = 0.0
      enddo  

      call rzero3(wk1,wk2,wk3,n)

      nface = 2*ndim
      iop   = 1 
      do i=1,fs_nel
        e   = fs_elno(i)
        ifc = fs_iface(i)
        
!       Get Normal directions      
        call facexv(unx(1,1,ifc,e),uny(1,1,ifc,e),unz(1,1,ifc,e),
     $              wk1(1,1,1,e),wk2(1,1,1,e),wk3(1,1,1,e),ifc,iop)
      enddo
!      call dsavg(wk1)
!      call dsavg(wk2)
!      if (ndim.eq.3) call dsavg(wk3)
!     Renormalize after dssum
!      call unitvec(wk1,wk2,wk3,n)
!      call fs_int_project(wk1,wk2,wk3)

      do 200 i=1,fs_nel
        e   = fs_elno(i)
        ifc = fs_iface(i)
        call facind2 (js1,jf1,jskip1,js2,jf2,jskip2,ifc)
        do 220 j2=js2,jf2,jskip2
        do 220 j1=js1,jf1,jskip1
!         normal component
          norx  =  wk1(j1,j2,1,e)
          nory  =  wk2(j1,j2,1,e)
          norz  =  0.0
          if (ndim.eq.3) norz  =  wk3(j1,j2,1,e)

          rnor = ux(j1,j2,1,e)*norx + uy(j1,j2,1,e)*nory
          if (ndim.eq.3) rnor = rnor + uz(j1,j2,1,e)*norz 

!         remove tangential component
          ux(j1,j2,1,e) = rnor*norx
          uy(j1,j2,1,e) = rnor*nory
          if (ndim.eq.3) uz(j1,j2,1,e) = rnor*norz

  220   continue                 
  200 continue

      call dsavg(ux)
      call dsavg(uy)
      if (ndim.eq.3) call dsavg(uz)


      do i=1,fs_nsymo
        e  = fs_ie(i)
        j1 = fs_ix(i)
        j2 = fs_iy(i)
        j3 = fs_iz(i)

        ux(j1,j2,j3,e) = wallvx(i)
        uy(j1,j2,j3,e) = wallvy(i)
      enddo   


      return
      end subroutine fs_mvmeshn2 
!---------------------------------------------------------------------- 
      subroutine fs_tang_corr(wx,wy,wz)

      implicit none

      include 'SIZE'
      include 'FS_ALE'
      include 'GEOM'
      include 'MASS'

      include 'TEST'
      include 'SOLN'

      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)
     
      integer i,n
      integer e,j1,j2,j3
      real p,s,r        ! 
      real mu           ! exponential decay
      real dy, tol

      integer ifc,nface,iop
      integer js1,js2,jf1,jf2,jskip1,jskip2

      character cb*3

!     Temporary arrays
      real wk1,wk2,wk3,wk4,wk5,wk6,wk7
      common /scrns3/ wk1(lx1,ly1,lz1,lelt),    ! Tangengial directions     
     $               wk2(lx1,ly1,lz1,lelt),     ! Tangengial directions
     $               wk3(lx1,ly1,lz1,lelt),     ! Tangengial directions
     $               wk4(lx1,ly1,lz1,lelt),     ! Displacement x 
     $               wk5(lx1,ly1,lz1,lelt),     ! Displacement y
     $               wk6(lx1,ly1,lz1,lelt),     ! Displacement z
     $               wk7(lx1,ly1,lz1,lelt)


      n = lx1*ly1*lz1*nelv

      call rzero3(wk1,wk2,wk3,n)
      nface = 2*ndim
      iop   = 1 
      do i=1,fs_nel
        e   = fs_elno(i)
        ifc = fs_iface(i)
        call facexv(t1x(1,1,ifc,e),t1y(1,1,ifc,e),t1z(1,1,ifc,e),
     $              wk1(1,1,1,e),wk2(1,1,1,e),wk3(1,1,1,e),ifc,iop)

      enddo
!      call dsavg(wk1)
!      call dsavg(wk2)
!      if (ndim.eq.3) call dsavg(wk3)
!!     Renormalize after dssum
!      call unitvec(wk1,wk2,wk3,n)
!      call fs_int_project(wk1,wk2,wk3)

      if (fs_ifgsm) then
        call fs_gllo_flds(wk1,wk2,wk3)
        call fs_get_localpts          ! SEM -> Global
!       Remove high wavenumbers      
        if (fs_iffil) then
          do i=1,ndim
            call fs_glfilter(fld_fs(1,1,i))
          enddo  
        endif
        call fs_get_globalpts         ! Global -> SEM
        call fs_restore_int(wk1,wk2,wk3,'Tang')    ! Tangential velocities
        call unitvec(wk1,wk2,wk3,n)
      endif       ! fs_ifgsm  

!      call opcopy(wk4,wk5,wk6,xm1,ym1,zm1)
!      call opsub2(wk4,wk5,wk6,xm0_fs,ym0_fs,zm0_fs) ! dx,dy,dz
     
      tol = 1.0e-14
      do i=1,n
        p  = wy(i,1,1,1)    ! current velocity along y
        s  = wk2(i,1,1,1)   ! tangential velocity along y
        if (abs(s).gt.tol) then
          r  = p/s
        else
          r  = 0.0
        endif  
        wk7(i,1,1,1) = r    
      enddo
      call col2(wk1,wk7,n)
      call col2(wk2,wk7,n)
      if (ndim.eq.3) call col2(wk3,wk7,n)
      
!      call localfilterfld(wk1)
!      call localfilterfld(wk2)
!      if (ndim.eq.3) call localfilterfld(wk3)
      call dsavg(wk1)
      call dsavg(wk2)
      if (ndim.eq.3) call dsavg(wk3)

!     No correction for SYM/O points
      do i=1,fs_nsymo
        e  = fs_ie(i)
        j1 = fs_ix(i)
        j2 = fs_iy(i)
        j3 = fs_iz(i)
        wk1(j1,j2,j3,e) = 0.
        wk2(j1,j2,j3,e) = 0.
        wk3(j1,j2,j3,e) = 0.
      enddo  

!      call fs_int_project(wk1,wk2,wk3)
      call opsub2(wx,wy,wz,wk1,wk2,wk3)

      return
      end subroutine fs_tang_corr
!----------------------------------------------------------------------
      subroutine fs_tang_corr4(wx,wy,wz)

      implicit none

      include 'SIZE'
      include 'FS_ALE'
      include 'GEOM'
      include 'MASS'

      include 'TEST'
      include 'SOLN'

      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)
     
      integer i,n
      integer e,j1,j2,j3
      real p,s,r        ! 
      real mu           ! exponential decay
      real dy, tol

      integer ifc,nface,iop
      integer js1,js2,jf1,jf2,jskip1,jskip2

      character cb*3

!     Temporary arrays
      real wk1,wk2,wk3,wk4,wk5,wk6,wk7
      common /scrns3/ wk1(lx1,ly1,lz1,lelt),    ! Tangengial directions     
     $               wk2(lx1,ly1,lz1,lelt),     ! Tangengial directions
     $               wk3(lx1,ly1,lz1,lelt),     ! Tangengial directions
     $               wk4(lx1,ly1,lz1,lelt),     ! Displacement x 
     $               wk5(lx1,ly1,lz1,lelt),     ! Displacement y
     $               wk6(lx1,ly1,lz1,lelt),     ! Displacement z
     $               wk7(lx1,ly1,lz1,lelt)



      n = lx1*ly1*lz1*nelv

      call rzero3(wk1,wk2,wk3,n)
      nface = 2*ndim
      iop   = 1 
      do i=1,fs_nel
        e   = fs_elno(i)
        ifc = fs_iface(i)
        call facexv(t1x(1,1,ifc,e),t1y(1,1,ifc,e),t1z(1,1,ifc,e),
     $              wk1(1,1,1,e),wk2(1,1,1,e),wk3(1,1,1,e),ifc,iop)

      enddo
      call dsavg(wk1)
      call dsavg(wk2)
      if (ndim.eq.3) call dsavg(wk3)
!     Renormalize after dssum
      call unitvec(wk1,wk2,wk3,n)
      call fs_int_project(wk1,wk2,wk3)

!      call opcopy(wk4,wk5,wk6,xm1,ym1,zm1)
!      call opsub2(wk4,wk5,wk6,xm0_fs,ym0_fs,zm0_fs) ! dx,dy,dz
     
      tol = 1.0e-14
      do i=1,n
        p  = wy(i,1,1,1)    ! current velocity along y
        s  = wk2(i,1,1,1)   ! tangential velocity along y
        if (abs(s).gt.tol) then
          r  = p/s
        else
          r  = 0.0
        endif  
        wk7(i,1,1,1) = r    
      enddo
      call fs_int_project(wk7,wk7,wk7) 
!      call localfilterfld(wk7)
!      call dsavg(wk7) 

      call col2(wk1,wk7,n)
      call col2(wk2,wk7,n)
      if (ndim.eq.3) call col2(wk3,wk7,n)

!!     No correction for SYM/O points
!      do i=1,fs_nsymo
!        e  = fs_ie(i)
!        j1 = fs_ix(i)
!        j2 = fs_iy(i)
!        j3 = fs_iz(i)
!        wk1(j1,j2,j3,e) = 0.
!        wk2(j1,j2,j3,e) = 0.
!        wk3(j1,j2,j3,e) = 0.
!      enddo  

      if (fs_ifgsm) then
        call fs_gllo_flds(wk1,wk2,wk3)
        call fs_get_localpts          ! SEM -> Global
!       Remove high wavenumbers      
        if (fs_iffil) then
          do i=1,ndim
            call fs_glfilter(fld_fs(1,1,i))
          enddo  
        endif
        call fs_get_globalpts         ! Global -> SEM
        call fs_restore_int(wk1,wk2,wk3,'Tang')    ! Tangential velocities
        call unitvec(wk1,wk2,wk3,n)
      endif       ! fs_ifgsm  

      call dsavg(wk1)
      call dsavg(wk2)
      if (ndim.eq.3) call dsavg(wk3)

!!     No correction for SYM/O points
!      do i=1,fs_nsymo
!        e  = fs_ie(i)
!        j1 = fs_ix(i)
!        j2 = fs_iy(i)
!        j3 = fs_iz(i)
!        wk1(j1,j2,j3,e) = 0.
!        wk2(j1,j2,j3,e) = 0.
!        wk3(j1,j2,j3,e) = 0.
!      enddo  

      call fs_int_project(wk1,wk2,wk3)

      call opsub2(wx,wy,wz,wk1,wk2,wk3)

      return
      end subroutine fs_tang_corr4
!----------------------------------------------------------------------
      subroutine fs_tang_corr2(wx,wy,wz)

      implicit none

      include 'SIZE'
      include 'FS_ALE'
      include 'GEOM'
      include 'MASS'

      include 'TEST'
      include 'SOLN'

      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)
     
      integer i,n
      integer e,j1,j2,j3
      real p,s,r        ! 
      real mu           ! exponential decay
      real dy, tol

      integer ifc,nface,iop
      integer js1,js2,jf1,jf2,jskip1,jskip2

      character cb*3

!     Temporary arrays
      real wk1,wk2,wk3,wk4,wk5,wk6,wk7
      common /scrns3/ wk1(lx1,ly1,lz1,lelt),    ! Tangengial directions     
     $               wk2(lx1,ly1,lz1,lelt),     ! Tangengial directions
     $               wk3(lx1,ly1,lz1,lelt),     ! Tangengial directions
     $               wk4(lx1,ly1,lz1,lelt),     ! Displacement x 
     $               wk5(lx1,ly1,lz1,lelt),     ! Displacement y
     $               wk6(lx1,ly1,lz1,lelt),     ! Displacement z
     $               wk7(lx1,ly1,lz1,lelt)



      n = lx1*ly1*lz1*nelv

      call rzero3(wk1,wk2,wk3,n)
      nface = 2*ndim
      iop   = 1 
      do i=1,fs_nel
        e   = fs_elno(i)
        ifc = fs_iface(i)
        call facexv(t1x(1,1,ifc,e),t1y(1,1,ifc,e),t1z(1,1,ifc,e),
     $              wk1(1,1,1,e),wk2(1,1,1,e),wk3(1,1,1,e),ifc,iop)

      enddo
!      call dsavg(wk1)
!      call dsavg(wk2)
!      if (ndim.eq.3) call dsavg(wk3)
!!     Renormalize after dssum
!      call unitvec(wk1,wk2,wk3,n)
!      call fs_int_project(wk1,wk2,wk3)

!     If global smoothening      
      if (fs_ifgsm) then
        call fs_gllo_flds(wk1,wk2,wk3)
        call fs_get_localpts          ! SEM -> Global
!       Remove high wavenumbers      
        if (fs_iffil) then
          do i=1,ndim
            call fs_glfilter(fld_fs(1,1,i))
          enddo  
        endif
        call fs_get_globalpts         ! Global -> SEM
        call fs_restore_int(wk1,wk2,wk3,'Tang')    ! Tangential velocities
        call unitvec(wk1,wk2,wk3,n)
      endif       ! fs_ifgsm
!     wk1,wk2,wk3 should now be smooth.      

!!     No corrections for SYM/O points
!      do i=1,fs_nsymo
!        e  = fs_ie(i)
!        j1 = fs_ix(i)
!        j2 = fs_iy(i)
!        j3 = fs_iz(i)
!        wk1(j1,j2,j3,e) = 0.
!        wk2(j1,j2,j3,e) = 0.
!        wk3(j1,j2,j3,e) = 0.
!      enddo  


!     Local point displacements      
      call opcopy(wk4,wk5,wk6,xm1,ym1,zm1)
      call opsub2(wk4,wk5,wk6,xm0_fs,ym0_fs,zm0_fs) ! dx,dy,dz
     
      tol = 1.0e-14
      mu  = 0.1
      do i=1,n
        dy = wk5(i,1,1,1)   ! displacement along y
        p  = wk2(i,1,1,1)   ! tangential velocity along y
        if (abs(p).gt.tol) then
          s  = -p/abs(p)
        else
          s  = 0.0
        endif 
        wk7(i,1,1,1) = s*mu*dy
      enddo

      call fs_int_project(wk7,wk7,wk7)
!      call localfilterfld(wk7)
      call dsavg(wk7) 

      call addcol3(wx,wk1,wk7,n)
      call addcol3(wy,wk2,wk7,n)
      if (ndim.eq.3) call addcol3(wz,wk3,wk7,n)

      return
      end subroutine fs_tang_corr2
!----------------------------------------------------------------------

      subroutine fs_tang_corr3(wx,wy,wz)

      implicit none

      include 'SIZE'
      include 'FS_ALE'
      include 'GEOM'
      include 'MASS'

      include 'TEST'
      include 'SOLN'

      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)
     
      integer i,n
      integer e,j1,j2,j3
      real p,s,r        ! 
      real mu           ! exponential decay
      real dy, tol

      integer ifc,nface,iop
      integer js1,js2,jf1,jf2,jskip1,jskip2

      character cb*3

!     Temporary arrays
      real wk1,wk2,wk3,wk4,wk5,wk6,wk7
      common /scrns3/ wk1(lx1,ly1,lz1,lelt),    ! Tangengial directions     
     $               wk2(lx1,ly1,lz1,lelt),     ! Tangengial directions
     $               wk3(lx1,ly1,lz1,lelt),     ! Tangengial directions
     $               wk4(lx1,ly1,lz1,lelt),     ! Displacement x 
     $               wk5(lx1,ly1,lz1,lelt),     ! Displacement y
     $               wk6(lx1,ly1,lz1,lelt),     ! Displacement z
     $               wk7(lx1,ly1,lz1,lelt)



      n = lx1*ly1*lz1*nelv

      call rzero3(wk1,wk2,wk3,n)
      nface = 2*ndim
      iop   = 1 
      do i=1,fs_nel
        e   = fs_elno(i)
        ifc = fs_iface(i)
        call facexv(t1x(1,1,ifc,e),t1y(1,1,ifc,e),t1z(1,1,ifc,e),
     $              wk1(1,1,1,e),wk2(1,1,1,e),wk3(1,1,1,e),ifc,iop)

      enddo
      call dsavg(wk1)
      call dsavg(wk2)
      if (ndim.eq.3) call dsavg(wk3)
!     Renormalize after dssum
      call unitvec(wk1,wk2,wk3,n)
      call fs_int_project(wk1,wk2,wk3)

!     Local point displacements      
      call opcopy(wk4,wk5,wk6,xm1,ym1,zm1)
      call opsub2(wk4,wk5,wk6,xm0_fs,ym0_fs,zm0_fs) ! dx,dy,dz
     
      tol = 1.0e-14
      mu  = 0.1
      do i=1,n
        dy = wk5(i,1,1,1)   ! displacement along y
        p  = wk2(i,1,1,1)   ! tangential velocity along y
        if (abs(p).gt.tol) then
          s  = -p/abs(p)
        else
          s  = 0.0
        endif 
        wk7(i,1,1,1) = s*mu*dy
      enddo
      call fs_int_project(wk7,wk7,wk7)

      call col2(wk1,wk7,n)
      call col2(wk2,wk7,n)
      if (ndim.eq.3) call col2(wk3,wk7,n)

!!     No (normal) corrections for SYM/O points
!      do i=1,fs_nsymo
!        e  = fs_ie(i)
!        j1 = fs_ix(i)
!        j2 = fs_iy(i)
!        j3 = fs_iz(i)
!        wk1(j1,j2,j3,e) = 0.
!        wk2(j1,j2,j3,e) = 0.
!        wk3(j1,j2,j3,e) = 0.
!      enddo  

!     If global smoothening      
      if (fs_ifgsm) then
        call fs_gllo_flds(wk1,wk2,wk3)
        call fs_get_localpts          ! SEM -> Global
!       Remove high wavenumbers      
        if (fs_iffil) then
          do i=1,ndim
            call fs_glfilter(fld_fs(1,1,i))
          enddo  
        endif
        call fs_get_globalpts         ! Global -> SEM
        call fs_restore_int(wk1,wk2,wk3,'Tang')    ! Tangential velocities
        call unitvec(wk1,wk2,wk3,n)
      endif       ! fs_ifgsm  

!!     No corrections for SYM/O points
!      do i=1,fs_nsymo
!        e  = fs_ie(i)
!        j1 = fs_ix(i)
!        j2 = fs_iy(i)
!        j3 = fs_iz(i)
!        wk1(j1,j2,j3,e) = 0.
!        wk2(j1,j2,j3,e) = 0.
!        wk3(j1,j2,j3,e) = 0.
!      enddo  

      call opadd2(wx,wy,wz,wk1,wk2,wk3)

      return
      end subroutine fs_tang_corr3
!----------------------------------------------------------------------
      subroutine fs_global_basis

      implicit none

      include 'SIZE'
      include 'FS_ALE'

      integer ix,iy,iz

      real xs2,ys2,xs4,ys4

      call rzero(zx_fs,lxfs)
      call rzero(zy_fs,lxfs)
      call rzero(wx_fs,lxfs)
      call rzero(wy_fs,lxfs)

      call zwgll(zx_fs,wx_fs,lxfs)
      call zwgll(zy_fs,wy_fs,lyfs)
   
      if (ndim.eq.3) then
        do ix=1,lxfs
        do iy=1,lyfs
          w2_fs(ix,iy) = wx_fs(ix)*wy_fs(iy)
        enddo
        enddo
      else
        call copy(w2_fs,wx_fs,lxfs)
      endif

!     Derivative matrices      
      call dgll (dx_fs,dxt_fs,zx_fs,lxfs,lxfs)
      call dgll (dy_fs,dyt_fs,zy_fs,lyfs,lyfs)

      return
      end subroutine fs_global_basis
!----------------------------------------------------------------------
      subroutine fs_gllo_xyz

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'          ! bm1, temporary
      include 'PARALLEL'
      include 'WZ'

      include 'FS_ALE'

      include 'GFLDR'

!      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
!      real wz(lx1,ly1,lz1,lelv)

      integer e,ifc,n,ne,nfaces

      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz,ia,ii,jj

      character cb*3
      real xm1_min,xm1_max,ym1_min,ym1_max,zm1_min,zm1_max
      real glmin,glmax

      real s(3)   ! surface normals
      real vn     ! normal velocity      

      xm1_min = glmin(xm1,lx1*ly1*lz1*nelv)
      xm1_max = glmax(xm1,lx1*ly1*lz1*nelv)
      ym1_min = glmin(ym1,lx1*ly1*lz1*nelv)
      ym1_max = glmax(ym1,lx1*ly1*lz1*nelv)
      zm1_min = glmin(zm1,lx1*ly1*lz1*nelv)
      zm1_max = glmax(zm1,lx1*ly1*lz1*nelv)

      call rzero(xg_fs,lxfs*lyfs*2)
      call rzero(yg_fs,lxfs*lyfs*2)
      call rzero(zg_fs,lxfs*lyfs*2)

      do iy=1,lyfs
        do ix=1,lxfs
          xg_fs(ix,iy,1)=(zx_fs(ix)+1.0)*(zm1_max-zm1_min)/2.0 + zm1_min
          yg_fs(ix,iy,1)=(zy_fs(iy)+1.0)*(ym1_max-ym1_min)/2.0 + ym1_min
          if (ndim.eq.2) xg_fs(ix,iy,1) = zx_fs(ix)
        enddo
      enddo

!     Get the surface x,y,z
      nfaces = 2*ndim
      ne     = 0              ! total number of interface elements
      ii     = 0
      do e=1,nelv
        jj     = 0
        do ifc=1,nfaces
          cb  = fs_cbc(ifc,e)
          if (cb.eq.'INT') then
            ne = ne+1
            call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,ifc)
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
!              call getSnormal(s,ix,iy,iz,ifc,e)   ! surface normals
              ii = ii + 1
              xi_fs(ii) = xm1(ix,iy,iz,e)
              yi_fs(ii) = ym1(ix,iy,iz,e)
              if (ndim.eq.3) zi_fs(ii) = zm1(ix,iy,iz,e)
              if (ndim.eq.2) xi_fs(ii) = zgm1(ix,1)

              if (ndim.eq.3) then
                jj = jj + 1
                xm1_fs(jj,1,ne) = xm1(ix,iy,iz,e)
                ym1_fs(jj,1,ne) = ym1(ix,iy,iz,e)
                xm1_fs(jj,1,ne) = zm1(ix,iy,iz,e)
              else
!               ndim.eq.2 needs very special treatment
                if (kx1.eq.kx2) then
                  do ia = 1,lx1
                    xm1_fs(ia,iy,ne)     = zgm1(ia,1)
                    ym1_fs(ia,iy,ne)     = ym1(ix,iy,iz,e)
                  enddo
                elseif (ky1.eq.ky2) then
                  do ia = 1,ly1
                    xm1_fs(ix,ia,ne)     = zgm1(ia,1)
                    ym1_fs(ix,ia,ne)     = ym1(ix,iy,iz,e)
                  enddo
                endif     ! kx1.eq.kx2
              endif       ! ndim.eq.2
            enddo         ! ix
            enddo         ! iy
            enddo         ! iz
          endif           ! cb.eq.INT
        enddo             ! ifc
      enddo               ! e

      
!!     debugging
!      call nekgsync()
!      if (nid.eq.0) then
!        do e=1,fs_nel
!          do iy=1,ly1
!            write(6,*) e, (ym1_fs(ix,iy,e),ix=1,lx1)
!          enddo
!        enddo
!      endif
!      call nekgsync()
!      call sleep(1)

      return
      end subroutine fs_gllo_xyz
!----------------------------------------------------------------------
       subroutine fs_gllo_flds(wx,wy,wz)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'          ! bm1, temporary
      include 'PARALLEL'
      include 'WZ'

      include 'FS_ALE'

      include 'GFLDR'

      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)

      integer e,ifc,n,ne,nfaces

      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz,ia,ii

      character cb*3

      real s(3)   ! surface normals
      real vn     ! normal velocity      

!     Get the surface x,y,z
      nfaces = 2*ndim
      ne     = 0              ! total number of interface elements
      do e=1,nelv
      do ifc=1,nfaces
        cb  = fs_cbc(ifc,e)
        if (cb.eq.'INT') then
          ne = ne+1
          ii     = 0
          call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,ifc)
          do iz=kz1,kz2
          do iy=ky1,ky2
          do ix=kx1,kx2
            ii = ii+1
!           call getSnormal(s,ix,iy,iz,ifc,e)   ! surface normals
            if (ndim.eq.3) then
              lfld_fs(ii,1,ne,1)= wx(ix,iy,iz,e)
              lfld_fs(ii,1,ne,2)= wy(ix,iy,iz,e)
              lfld_fs(ii,1,ne,3)= wz(ix,iy,iz,e)
            else  
!             ndim.eq.2 needs very special treatment
              if (kx1.eq.kx2) then
                do ia = 1,lx1
                  lfld_fs(ia,iy,ne,1) = wx(ix,iy,iz,e)
                  lfld_fs(ia,iy,ne,2) = wy(ix,iy,iz,e)
                enddo
              elseif (ky1.eq.ky2) then
                do ia = 1,ly1
                  lfld_fs(ix,ia,ne,1) = wx(ix,iy,iz,e)
                  lfld_fs(ix,ia,ne,2) = wy(ix,iy,iz,e)
                enddo
              endif     ! kx1.eq.kx2
            endif       ! ndim.eq.2
          enddo         ! ix
          enddo         ! iy
          enddo         ! iz
        endif           ! cb.eq.INT
      enddo             ! ifc
      enddo             ! e

!!     debugging
!      call nekgsync()
!      if (nid.eq.0) then
!        do e=1,fs_nel
!          do iy=1,ly1
!            write(6,*) e, (lfld_fs(ix,iy,e,1),ix=1,lx1)
!          enddo
!        enddo
!      endif
!      call nekgsync()
!      call sleep(1)


      return
      end subroutine fs_gllo_flds
!----------------------------------------------------------------------
     
      subroutine fs_intp_setup

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'GFLDR'
      include 'PARALLEL'

      include 'FS_ALE'

      integer i

!     testing interpolation
      integer nxf,nyf,nzf
      integer nhash,nmax
      integer ldim2

      integer nintp           ! no of interpolation points
      
      integer nidd,npp,nekcomm,nekgroup,nekreal
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal


!     initialize interpolation tool using global mesh
      nxf   = 2*lx1
      nyf   = 2*ly1
      nzf   = 2*1
      nhash = fs_nel*lx1*ly1
      if (nhash.eq.0) then
        nhash = lx1*ly1
      endif  
      nmax  = 128
!     We do all Global calculations on nid 0
      if (nid.eq.0) then
        nels  = 1
      else
        nels  = 0
      endif  

      ldim2 = 2


!     Interpolation handle for Global mesh.      
      call fgslib_findpts_setup(intgh_fs,nekcomm,np,ldim2,
     &                          xg_fs,yg_fs,zg_fs,lxfs,lyfs,lzfs,
     &                          nels,nxf,nyf,nzf,bb_t,
     &                          nhash,nhash,nmax,tol)

      call mntr_log(fs_id,fs_log,'Global Interpolation Setup Done')

!     initialize interpolation tool using local sem mesh
      nxf   = 2*lx1
      nyf   = 2*ly1
      nzf   = 2*1
      nhash = 1*lx1*ly1
      nmax  = 128
!     We do all Global calculations on nid 0
      nels  = fs_nel
     
!     Interpolation handle for SEM surface mesh. 
      call fgslib_findpts_setup(intlh_fs,nekcomm,np,ldim2,
     &                          xm1_fs,ym1_fs,zm1_fs,lx1,ly1,1,
     &                          nels,nxf,nyf,nzf,bb_t,
     &                          nhash,nhash,nmax,tol)

      call mntr_log(fs_id,fs_log,'Local Interpolation Setup Done')

      return
      end subroutine fs_intp_setup
!---------------------------------------------------------------------- 
      subroutine fs_get_globalpts

!     Find the local sem surface points within the global mesh.

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'GFLDR'
      include 'PARALLEL'

      include 'FS_ALE'

      integer i

!     testing interpolation
      integer nxf,nyf,nzf
      integer nhash,nmax
      integer ldim2
      real xin,yin,zin,fldout

      integer nintp           ! no of interpolation points
      integer*8 nfail
      integer*8 nfail_sum
      integer*8 i8glsum

      real toldist

      integer ix,iy

      nfail = 0
      toldist = 1e-12
!
      if (ndim.eq.2) then      
        nintp  = fs_nel*lx1*ly1
      else
        nintp  = fs_nel*lx1*ly1
      endif
      ldim2  = 2

      call fgslib_findpts(intgh_fs,
     &                    grcode,1,
     &                    gproc,1,
     &                    gelid,1,
     &                    grst,ldim2,
     &                    gdist,1,
     &                    xm1_fs,1,
     &                    ym1_fs,1,
     &                    zm1_fs,1,nintp)


      do i=1,nintp
         if(grcode(i).eq.1 .and. sqrt(gdist(i)).gt.toldist)
     &     nfail = nfail + 1
         if(grcode(i).eq.2) nfail = nfail + 1
      enddo

      nfail_sum = i8glsum(nfail,1)
      if(nfail_sum.gt.0) then
        if(nio.eq.0) write(6,*)
     &    ' WARNING: Unable to find all mesh points in source fld ',
     &    nfail_sum
        if(nio.eq.0) write(6,*)
     &    ' WARNING: in fs_get_globalpts()'
      endif

      call rzero(fs_gfldout,lx1*ly1*lelv*3)
!     evaluate input field at sem points
      do i=1,ndim
        call fgslib_findpts_eval(intgh_fs,
     &                           fs_gfldout(1,1,1,i),1,
     &                           grcode,1,
     &                           gproc,1,
     &                           gelid,1,
     &                           grst,ldim2,nintp,
     &                           fld_fs(1,1,i))
      enddo  

      call mntr_log(fs_id,fs_log,'Global -> SEM: Done')

      return
      end subroutine fs_get_globalpts
!----------------------------------------------------------------------
      subroutine fs_get_localpts

!     Find the Global mesh points within the local sem surface mesh.

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'GFLDR'
      include 'PARALLEL'

      include 'FS_ALE'

      integer i,j

!     testing interpolation
      integer nxf,nyf,nzf
      integer nhash,nmax
      integer ldim2
      real xin,yin,zin,fldout

      integer nintp           ! no of interpolation points
      integer*8 nfail
      integer*8 nfail_sum
      integer*8 i8glsum

      real toldist

      nfail = 0
      toldist = 5e-14

      if (nid.eq.0) then
        nintp  = lxfs*lyfs
      else
        nintp  = 0  
      endif
      ldim2  = 2

      call fgslib_findpts(intlh_fs,
     &                    grcode,1,
     &                    gproc,1,
     &                    gelid,1,
     &                    grst,ldim2,
     &                    gdist,1,
     &                    xg_fs,1,
     &                    yg_fs,1,
     &                    zg_fs,1,nintp)


      do i=1,nintp
         if(grcode(i).eq.1 .and. sqrt(gdist(i)).gt.toldist)
     &     nfail = nfail + 1
         if(grcode(i).eq.2) nfail = nfail + 1
      enddo

      nfail_sum = i8glsum(nfail,1)
      if(nfail_sum.gt.0) then
        if(nio.eq.0) write(6,*)
     &    ' WARNING: Unable to find all mesh points in source fld ',
     &    nfail_sum
        if(nio.eq.0) write(6,*)
     &    ' WARNING: in fs_get_localpts()'
      endif

!     Evaluate fields at global mesh points      
      do i=1,ndim
!       evaluate inut field at given points
        call fgslib_findpts_eval(intlh_fs,
     &                           fs_lfldout(1,1,i),1,
     &                           grcode,1,
     &                           gproc,1,
     &                           gelid,1,
     &                           grst,ldim2,nintp,
     &                           lfld_fs(1,1,1,i))

!        write(6,*) 'Fldout', (fs_lfldout(j,1,i),j=1,nintp)

!       This is now the globally smooth field from which we interpolate
!       the local sem points 
        if (nintp.gt.0) then
          call copy(fld_fs(1,1,i),fs_lfldout(1,1,i),nintp)
        else
          call rzero(fld_fs(1,1,i),lxfs*lyfs)
        endif  
      enddo  

      call mntr_log(fs_id,fs_log,'SEM -> Global: Done')

      return
      end subroutine fs_get_localpts
!----------------------------------------------------------------------

      subroutine fs_restore_int(wx,wy,wz,dirc)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'          ! bm1, temporary
      include 'PARALLEL'
      include 'TSTEP'

      include 'FS_ALE'

      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)

      integer e,ifc,n,ne,nfaces

      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz,ia,ii

      character cb*3
      real xm1_min,xm1_max,ym1_min,ym1_max
      real glmin,glmax,glsum

      real erx,ery,erz,err         ! errors due to smoothening
      real ar,arsum                ! area

      character str*120
      integer loglev

      character dirc*4

!     Get the surface x,y,z
      nfaces = 2*ndim

      ar  = 0.0
      ne  = 0
      erx = 0.0
      ery = 0.0
      erz = 0.0
      do e=1,nelv
      do ifc=1,nfaces
        cb  = fs_cbc(ifc,e)
        if (cb.eq.'INT') then
          ne = ne + 1
          call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,ifc)
          ii = 0
          do iz=kz1,kz2
          do iy=ky1,ky2
          do ix=kx1,kx2
            ii = ii+1
            erx = erx + 
     $      area(ii,1,ifc,e)*(fs_gfldout(ix,iy,ne,1)-wx(ix,iy,iz,e))**2
            ery = ery + 
     $      area(ii,1,ifc,e)*(fs_gfldout(ix,iy,ne,2)-wy(ix,iy,iz,e))**2
            if (ndim.eq.3)
     $        erz = erz + 
     $      area(ii,1,ifc,e)*(fs_gfldout(ix,iy,ne,3)-wz(ix,iy,iz,e))**2
          
            ar = ar +  area(ii,1,ifc,e)

!            wx(ix,iy,iz,e) = fs_gfldout(ii,1,ne,1)
!            wy(ix,iy,iz,e) = fs_gfldout(ii,1,ne,2)
!            if (ndim.eq.3) wz(ix,iy,iz,e) = fs_gfldout(ii,1,ne,3)

            if (ndim.eq.3) then
              wx(ix,iy,iz,e) = fs_gfldout(ii,1,ne,1) 
              wy(ix,iy,iz,e) = fs_gfldout(ii,1,ne,2) 
              wz(ix,iy,iz,e) = fs_gfldout(ii,1,ne,3) 
            else  
!             ndim.eq.2 needs very special treatment
              wx(ix,iy,iz,e) = lfld_fs(ix,iy,ne,1)
              wy(ix,iy,iz,e) = lfld_fs(ix,iy,ne,2)
            endif       ! ndim.eq.2
          enddo
          enddo
          enddo
        endif
      enddo
      enddo  

      err = (erx + ery + erz)
      erx = glsum(err,1)      ! now contains total error
      arsum = glsum(ar,1)
      err = sqrt(erx/arsum)

      write(str,'(A4,1x,A24,1x,E11.4E2)') dirc,
     $      'Smooth projection error:',err

!     We (almost) always want to see this error      
      loglev = fs_log
      call mntr_log(fs_id,loglev,str)

      return
      end subroutine fs_restore_int
!----------------------------------------------------------------------

      subroutine fs_glfilter(fld)

      implicit none

      include 'SIZE'
!      include 'SOLN'
!      include 'INPUT'         ! param(110),(111)
      include 'TSTEP'         ! ifield
      include 'MASS'          ! BM1
      include 'FS_ALE'

      real fld(1)       ! Field to filter

      integer lm2
      parameter (lm2=lxym*lxym)
      real filterfcn(lm2)
      real wght
      integer kut

      integer icalld
      save icalld
      data icalld /0/

      kut = int(lxfs/4)+0

      if (icalld.eq.0) then
!       Create the filter transfer function
        wght = 1.00
        call fs_trns_fcn(filterfcn,kut,wght)

!       Build the matrix to apply the filter function
!       to an input field
        call fs_filtermat(glfiltop_fs,filterfcn)

!       Only initialize once    
        icalld=icalld+1 
      endif

!     Apply the filter
      call fs_filterfld(fld,glfiltop_fs,lxfs)

      return
      end subroutine fs_glfilter
!---------------------------------------------------------------------- 
      subroutine fs_trns_fcn(diag,kut,wght)

      implicit none

      include 'SIZE'
      include 'FS_ALE'

      real diag(lxym*lxym)
      integer nx,k0,kut,kk,k

      real amp,wght

      nx = lxym
      call ident   (diag,nx)

      k0 = nx-kut
      do k=k0+1,nx
        kk = k+nx*(k-1)
        amp = wght*((k-k0)*(k-k0)+0.)/(kut*kut+0.)     ! Normalized amplitude. quadratic growth
        diag(kk) = 1.-amp
      enddo

c     Output normalized transfer function
!      k0 = lxym+1
!      if (nio.eq.0) then
!        write(6,6) 'HPF :',((1.-diag(k)), k=1,lx1*lx1,k0)
!   6    format(a8,16f9.6,6(/,8x,16f9.6))
!      endif

      return
      end subroutine fs_trns_fcn

!---------------------------------------------------------------------- 

      subroutine fs_filtermat(op_mat,f_filter)

c     Builds the operator for high pass filtering
c     Transformation matrix from nodal to modal space.
c     Applies f_filter to the the legendre coefficients
c     Transforms back to nodal space
c     Operation: V * f_filter * V^(-1)
c     Where V is the transformation matrix from modal to nodal space

      implicit none

      include 'SIZE'
      include 'FS_ALE'

      integer lm2
      parameter (lm2=lxym*lxym)

      real f_filter(lm2)
      real op_mat(lxym,lxym)

      real ref_xmap(lm2)
      real wk_xmap(lm2)

      real wk1(lm2),wk2(lm2)
      real indr(lxym),ipiv(lxym),indc(lxym)

      real rmult(lxym)
      integer ierr

      call fs_speccoeff_init(ref_xmap)
      
      call copy(wk_xmap,ref_xmap,lm2)
      call copy(wk1,wk_xmap,lm2)

      call gaujordf  (wk1,lxym,lxym,indr,indc,ipiv,ierr,rmult)  ! xmap inverse

      call mxm  (f_filter,lxym,wk1,lxym,wk2,lxym)        !          -1
      call mxm  (wk_xmap,lxym,wk2,lxym,op_mat,lxym)      !     V D V

      return
      end

!---------------------------------------------------------------------- 
     
      subroutine fs_speccoeff_init(ref_xmap)
c     Initialise spectral coefficients
c     For legendre transform

      implicit none

      include 'SIZE'
      include 'WZ'
      include 'FS_ALE'

      integer lm2
      parameter (lm2=lxym*lxym)

c     local variables
      integer i, j, k, n, nx, kj
c     Legendre polynomial
      real plegx(lxym)
      real z
      real ref_xmap(lm2)
      real pht(lm2)

c     Change of basis
      logical ifboyd

!     We don't want to change boundary values      
      ifboyd = .false.

      nx = lxym
      kj = 0
      n  = nx-1
      do j=1,nx
        z = zx_fs(j)
        call legendre_poly(plegx,z,n)
        kj = kj+1
        pht(kj) = plegx(1)
        kj = kj+1
        pht(kj) = plegx(2)

        if (ifboyd) then        ! change basis to preserve element boundary values
          do k=3,nx
             kj = kj+1
             pht(kj) = plegx(k)-plegx(k-2)
          enddo
        else                    ! legendre basis    
          do k=3,nx
             kj = kj+1
             pht(kj) = plegx(k)
          enddo         
        endif
      enddo

      call transpose (ref_xmap,nx,pht,nx)

      return
      end subroutine fs_speccoeff_init
!---------------------------------------------------------------------- 
      subroutine fs_filterfld(v,f,nx)

c     Appies the operator f to field u
c     using tensor operations
c     v = f*u

      implicit none

      include 'SIZE'
      include 'FS_ALE'

      integer nxyz 
      parameter (nxyz=lxfs*lyfs*lzfs)

      real w1(nxyz),w2(nxyz)              ! work arrays
      real v(nxyz)                        ! output fld
      real u(nxyz)                        ! input fld

      integer nx

      real f(nx,nx),ft(nx,nx)             ! operator f and its transpose

      integer e,i,j,k
      integer nel

!      call copy(v,u,nxyz)

      call transpose(ft,nx,f,nx)

!     Filter
      call copy(w1,v,nxyz)
      call mxm(f ,nx,w1,nx,w2,nx)
      call mxm(w2,nx,ft,nx,w1,nx)         ! w1 is low pass filtered 

!      call sub3(w2,u,w1,nxyz)
!      call copy(v,w2,nxyz)
      call copy(v,w1,nxyz) 

      return
      end subroutine fs_filterfld

!---------------------------------------------------------------------- 
      subroutine fs_fixcorners(wx,wy,wz)

      implicit none

      include 'SIZE'
      include 'FS_ALE'

      integer i,e,j1,j2,j3
      real wx(lx1,ly1,lz1,lelv),wy(lx1,ly1,lz1,lelv)
      real wz(lx1,ly1,lz1,lelv)

!     Save corner point velocities      
      do i = 1,fs_nsymo
        e  = fs_ie(i)
        j1 = fs_ix(i)
        j2 = fs_iy(i)
        j3 = fs_iz(i)
!        wallvx(i) = wx(j1,j2,j3,e)
!        wallvy(i) = 0.0
        wy(j1,j2,j3,e) = 0.0
      enddo  


      return
      end subroutine fs_fixcorners        
!---------------------------------------------------------------------- 




