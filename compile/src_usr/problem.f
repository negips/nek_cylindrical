!======================================================================
!     Author: Prabal Negi      
!     Description: Problem specific routines
!     Routines:   modify_mask()              : Modify v1mask
!
!======================================================================

      
      subroutine modify_mask()

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'GEOM'
      include 'NEKUSE'
      include 'CYLINDRICAL'
      include 'FS_ALE'

      integer ie,iface,nfaces,ifld
      integer nxyz,nx,ny,nz
      real radmean

      integer kx1,kx2,ky1,ky2,kz1,kz2
      integer ix,iy,iz

      real slipl        ! slip length
      real blendl       ! blend length
      real alpha,beta
      real intdis

      real xs,xf,fsx
      integer n
      real mu

      common  /nekcb/ cb
      character cb*3

      nfaces = 2*ndim
      ifld   = 1
      nxyz   = lx1*ly1*lz1
      nx     = lx1
      ny     = ly1
      nz     = lz1

      slipl  = 0.01
      blendl = 0.025
      n      = 4
      mu     = 0.025

!     Get interface position
      call fs_get_intpos(fs_intpos)

      do ie=1,nelv
        do iface=1,nfaces
          cb  = cbc(iface,ie,ifld)
          call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface)

!          call surface_int(sint,sarea,ym1,ie,iface)
!          sint = sint/sarea
!          if (cb.eq.'O  ') then
          if (cb.eq.'v  ') then
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
              if (optlevel.le.2) call nekasgn (ix,iy,iz,ie)
              fsx      = fs_intpos(ix,iy,iz,ie)  ! free surface position
              intdis   = abs(x-fsx)
              xs       = slipl               ! start of blending
              xf       = slipl+blendl        ! end of blending
!             Smoothly blend from free slip to Dirichlet                    
              if (intdis.le.xs) then
                v1mask(ix,iy,iz,ie) = 1.0        ! Free slip in this region
              elseif (intdis.le.xf) then
!                alpha = ((intdis-xs)/(xf-xs))**n
                intdis = intdis - xs
                alpha  = exp(-(intdis/mu)**2)
                v1mask(ix,iy,iz,ie) = alpha      ! Blended mask in this region
              else
                v1mask(ix,iy,iz,ie) = 0.0        ! Dirichlet
              endif                   
            enddo      ! ix
            enddo      ! iy
            enddo      ! iz
          endif        ! cb.eq.v            
        enddo          ! iface
      enddo            ! ie

      return
      end subroutine
!---------------------------------------------------------------------- 
