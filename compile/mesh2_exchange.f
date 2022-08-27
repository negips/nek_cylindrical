!======================================================================
!     Author: Prabal Negi      
!     Setting up communicator to exchange Mesh2
!     data across neighboring elements.
!
!====================================================================== 
      subroutine setup_exchange_m2()

!     Exchange interior values of neighboring elements
!     On the Pressure mesh        

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      include 'TEST'          ! for debugging

      integer lt,lt2
      parameter (lt=lx1*ly1*lz1*lelv)
      parameter (lt2=lx2*ly2*lz2*lelv)

      integer n,n2,nxyz,nxyz2

      real scr1
      real scpr               ! global point numbering on Mesh2 
      common /scrch/ scr1(lt),scpr(lx2,ly2,lz2,lelt)

      integer i,j,k,e
      integer gl,gl0

      n  = lx1*ly1*lz1*nelv
      n2 = lx2*ly2*lz2*nelv

      nxyz  = lx1*ly1*lz1
      nxyz2 = lx2*ly2*lz2
     
!     Create a global numbering.
!     This is easy since there are no coincident points on Mesh2
!     All points are unique.      

      do e=1,nelv
        gl   = lglel(e)             ! global element number
        gl0  = (gl-1)*nxyz2
        do i=1,nxyz2
          scpr(i,1,1,e) = gl0 + i + 0.0
        enddo
      enddo
      call copy(tmp4,scpr,n2)

      call rzero(scr1,n)
      call fill_interior(scr1,scpr)       ! Put M2 on M1 interior
      call dface_ext   (scr1)             ! Extend to face
      call edge_ext    (scr1)             ! Extend to Edges
      call crnr_ext    (scr1)             ! Extend to Corners
      call dssum       (scr1,lx1,ly1,lz1) ! Add up all contributions
      call dface_add1si(scr1,-1.)         ! Here we have exchanged face values
      call edge_corr   (scr1)
      call crnr_corr   (scr1)

      call copy(tmp1,scr1,n)

      call outpost(tmp1,tmp2,tmp3,tmp4,tmp3,'exc')

      return
      end subroutine 
!---------------------------------------------------------------------- 

      subroutine exchange_m2(x1,x2)

!     Exchange interior values of neighboring elements
!     On the Pressure mesh        

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      integer n,n2

      real x1(lx1,ly1,lz1,lelv)     ! output
      real x2(lx2,ly2,lz2,lelv)     ! input

      n  = lx1*ly1*lz1*nelv
      n2 = lx2*ly2*lz2*nelv

      call rzero(x1,n)
      call fill_interior(x1,x2)       ! Put M2 on M1 interior
      call dface_ext   (x1)             ! Extend to face
      call edge_ext    (x1)             ! Extend to Edges
      call crnr_ext    (x1)             ! Extend to Corners
      call dssum       (x1,lx1,ly1,lz1) ! Add up all contributions
      call dface_add1si(x1,-1.)         ! Here we have exchanged face values
      call edge_corr   (x1)             ! Correct the edge values
      call crnr_corr   (x1)             ! Correct the corner values

      return
      end subroutine 
!---------------------------------------------------------------------- 

      subroutine edge_ext(x)

      implicit none        

      include 'SIZE'
      include 'INPUT'

      real x(lx1,ly1,lz1,1)
      integer ie,i,j,k
      integer ix,iy,iz

      do ie=1,nelv

         if (if3d) then
           do k=1,nz1,(nz1-1)
           do j=1,ny1,(ny1-1)
           do i=2,nx1-1
             ix = i+1
             if (ix.gt.lx1) ix=lx1-1
             iy = j+1
             if (iy.gt.ly1) iy=ly1-1
             iz = k+1
             if (iz.gt.lz1) iz=lz1-1
             x(i,j,k,ie)      = x(i,iy,k,ie)
             x(i,j,k,ie)      = x(i,j,iz,ie)
           enddo
           enddo
           enddo 

           do k=1,nz1,(nz1-1)
           do j=2,ny1-1
           do i=1,nx1,(nx1-1)
             ix = i+1
             if (ix.gt.lx1) ix=lx1-1
             iy = j+1
             if (iy.gt.ly1) iy=ly1-1
             iz = k+1
             if (iz.gt.lz1) iz=lz1-1
             x(i,j,k,ie)      = x(ix,j,k,ie)
             x(i,j,k,ie)      = x(i,j,iz,ie)
           enddo
           enddo
           enddo 

           do k=2,nz1-1
           do j=1,ny1,(ny1-1)
           do i=1,nx1,(nx1-1)
             ix = i+1
             if (ix.gt.lx1) ix=lx1-1
             iy = j+1
             if (iy.gt.ly1) iy=ly1-1
             iz = k+1
             if (iz.gt.lz1) iz=lz1-1
             x(i,j,k,ie)      = x(ix,j,k,ie)
             x(i,j,k,ie)      = x(i,iy,k,ie)
           enddo
           enddo
           enddo 

!           do i=2,lx1-1
!             x(i,1,1,ie)      = x(i,2,1,ie)
!             x(i,ly1,1,ie)    = x(i,ly1-1,1,ie)
!             x(i,1,lz1,ie)    = x(i,2,lz1,ie)
!             x(i,ly1,lz1,ie)  = x(i,ly1-1,lz1,ie)
!           enddo  
!
!           
!           do j=2,ly1-1
!             x(1,j,1,ie)      = x(1,j,2,ie)
!             x(1,j,lz1,ie)    = x(1,j,lz1-1,ie)
!             x(lx1,j,1,ie)    = x(lx1,j,2,ie)
!             x(lx1,j,lz1,ie)  = x(lx1,j,lz1-1,ie)
!           enddo  
!
!           do k=2,lz1-1
!             x(1,1,k,ie)      = x(2,1,k,ie)
!             x(lx1,1,k,ie)    = x(lx1-1,1,k,ie)
!             x(1,ly1,k,ie)    = x(2,ly1,k,ie)
!             x(lx1,ly1,k,ie)  = x(lx1-1,ly1,k,ie)
!           enddo  
         
         else
!          No Edge extensions in 2D.
!          Only Faces And Corners           

        endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine edge_corr(x)

      implicit none

      include 'SIZE'
      include 'INPUT'

      integer ie

      real x(lx1,ly1,lz1,1)
      real s
      real f1,f2,f3           ! face correction
      real lc                 ! local element correction
      integer i,j,k
      integer ix,iy,iz

      do ie=1,nelv

         if (if3d) then

           do k=1,nz1,(nz1-1)
           do j=1,ny1,(ny1-1)
           do i=2,nx1-1
             ix = i+1
             if (ix.gt.lx1) ix=lx1-1
             iy = j+1
             if (iy.gt.ly1) iy=ly1-1
             iz = k+1
             if (iz.gt.lz1) iz=lz1-1
             s           = x(i,j,k,ie)          ! all contributions
             f1          = x(i,iy,k,ie)         ! face 1
             f2          = x(i,j,iz,ie)         ! face 2
             lc          = x(i,iy,iz,ie)        ! internal (local)
             x(i,j,k,ie) = s - (f1+f2) - lc
           enddo
           enddo
           enddo 

           do k=1,nz1,(nz1-1)
           do j=2,ny1-1
           do i=1,nx1,(nx1-1)
             ix = i+1
             if (ix.gt.lx1) ix=lx1-1
             iy = j+1
             if (iy.gt.ly1) iy=ly1-1
             iz = k+1
             if (iz.gt.lz1) iz=lz1-1
             s           = x(i,j,k,ie)          ! all contributions
             f1          = x(ix,j,k,ie)         ! face 1
             f2          = x(i,j,iz,ie)         ! face 2
             lc          = x(ix,j,iz,ie)        ! internal (local)
             x(i,j,k,ie) = s - (f1+f2) - lc
           enddo
           enddo
           enddo 

           do k=2,nz1-1
           do j=1,ny1,(ny1-1)
           do i=1,nx1,(nx1-1)
             ix = i+1
             if (ix.gt.lx1) ix=lx1-1
             iy = j+1
             if (iy.gt.ly1) iy=ly1-1
             iz = k+1
             if (iz.gt.lz1) iz=lz1-1
             s           = x(i,j,k,ie)          ! all contributions
             f1          = x(ix,j,k,ie)         ! face 1
             f2          = x(i,iy,k,ie)         ! face 2
             lc          = x(ix,iy,k,ie)        ! internal (local)
             x(i,j,k,ie) = s - (f1+f2) - lc
           enddo
           enddo
           enddo 

         else
!           No Edge Correction for 2D            

        endif
      enddo

      return
      end
c-----------------------------------------------------------------------

      subroutine crnr_ext(x)

      implicit none

      include 'SIZE'
      include 'INPUT'

      real x(lx1,ly1,lz1,lelt)
      integer ie,i,j,k
      integer ix,iy,iz
      real s
c
      do ie=1,nelv
c
         if (if3d) then

           do k=1,nz1,(nz1-1)
           do j=1,ny1,(ny1-1)
           do i=1,nx1,(nx1-1)
             ix = i+1
             if (ix.gt.lx1) ix=lx1-1
             iy = j+1
             if (iy.gt.ly1) iy=ly1-1
             iz = k+1
             if (iz.gt.lz1) iz=lz1-1
             s           = x(ix,iy,iz,ie)
             x(i,j,k,ie) = s
           enddo
           enddo
           enddo

         else

           do j=1,ny1,(ny1-1)
           do i=1,nx1,(nx1-1)
             ix = i+1
             if (ix.gt.lx1) ix=lx1-1
             iy = j+1
             if (iy.gt.ly1) iy=ly1-1
             s           = x(ix,iy,1,ie)
             x(i,j,1,ie) = s
           enddo
           enddo

        endif
      enddo
      return
      end
c-----------------------------------------------------------------------

      subroutine crnr_corr(x)

      implicit none

      include 'SIZE'
      include 'INPUT'

      integer ie

      real x(lx1,ly1,lz1,1)
      real s
      real f1,f2,f3           ! face correction
      real e1,e2,e3           ! edge correction
      real lc                 ! local element correction
      integer i,j,k
      integer ix,iy,iz

      do ie=1,nelv
c
         if (if3d) then

           do k=1,nz1,(nz1-1)
           do j=1,ny1,(ny1-1)
           do i=1,nx1,(nx1-1)
             ix = i+1
             if (ix.gt.lx1) ix=lx1-1
             iy = j+1
             if (iy.gt.ly1) iy=ly1-1
             iz = k+1
             if (iz.gt.lz1) iz=lz1-1
!            all contributions           
             s                   = x(i,j,k,ie)
!            face correction            
             f1                  = x(ix,iy,k,ie)
             f2                  = x(ix,j,iz,ie)
             f3                  = x(i,iy,iz,ie)
!            edge correction 
             e1                  = x(ix,j,k,ie)
             e2                  = x(i,iy,k,ie)
             e3                  = x(i,j,iz,ie)
!            local correction 
             lc                  = x(ix,iy,iz,ie)
!            Corrected value 
             x(i,j,k,ie)         = s - (f1+f2+f3) - (e1+e2+e3) - lc
           enddo
           enddo
           enddo


         else

           do j=1,ny1,(ny1-1)
           do i=1,nx1,(nx1-1)
             ix = i+1
             if (ix.gt.lx1) ix=lx1-1
             iy = j+1
             if (iy.gt.ly1) iy=ly1-1
!            all contributions           
             s                   = x(i,j,1,ie)
!            face correction            
             f1                  = x(ix,j,1,ie)
             f2                  = x(i,iy,1,ie)
!            No edge corrections in 2D. 
!            local correction 
             lc                  = x(ix,iy,1,ie)
!            Corrected value 
             x(i,j,1,ie)         = s - (f1+f2) - lc
           enddo
           enddo

        endif
      enddo

      return
      end
c-----------------------------------------------------------------------

      subroutine fill_interior(x,y)

!     Put Mesh2 node values on to the interior
!     Mesh1 nodes.

      implicit none

      include 'SIZE'

      real x(lx1,ly1,lz1,lelv)
      real y(lx2,ly2,lz2,lelv)

      integer iz1,e,ix,iy,iz
      integer n

      n = lx1*ly1*lz1*nelv

      call rzero(x,n)
!     Fill interiors
      iz1 = 0
      if (ndim.eq.3) iz1=1
      do e=1,nelv
         do iz=1,lz2
         do iy=1,ly2
         do ix=1,lx2
           x(ix+1,iy+1,iz+iz1,e) = y(ix,iy,iz,e)
         enddo
         enddo
         enddo
      enddo

      return
      end subroutine 
!---------------------------------------------------------------------- 

       subroutine extract_interior(x,y)

      implicit none

      include 'SIZE'

      real y(lx1,ly1,lz1,lelv)
      real x(lx2,ly2,lz2,lelv)

      integer iz1,e,ix,iy,iz
      integer n
c
c     Map back to pressure grid (extract interior values)
c
      do e=1,nelv
         do iz=1,lz2
         do iy=1,ly2
         do ix=1,lx2
            x(ix,iy,iz,e) = y(ix+1,iy+1,iz+iz1,e)
         enddo
         enddo
         enddo
      enddo


      return
      end subroutine 
!---------------------------------------------------------------------- 

