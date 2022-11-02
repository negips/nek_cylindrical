!====================================================================== 
!     Author: Prabal Negi
!     Description: FDM preconditioner for cylindrical solve.
!      
!      
!====================================================================== 
      subroutine set_up_fast_1D_sem_cyl(s,lam,n,lbc,rbc,ll,lm,lr,
     $                                    rho,rad,ie,isd)

      implicit none

      include 'SIZE'
      include 'SEMHAT'
c
      common /fast1dsem/ g(lr2),w(lr2)
c
      real g,w
      real s(1),lam(1),ll,lm,lr
      integer lbc,rbc
      
      integer bb0,bb1,eb0,eb1,n,n1
      logical l,r

      integer ie
      real rho(lx1)     ! density
      real rad(lx1)     ! radius
      real rmean        ! mean radius
      real mri2         ! mean radius inverse squared
      real vlsc2

      real dummy(lx1)
      integer isd
     
      n=lx1-1
!     BCs on E are from normal vel component
      if(lbc.eq.2 .or. lbc.eq.3) then !wall,sym - dirichlet velocity
         eb0=1
      else !outflow,element - neumann velocity
         eb0=0
      endif
      if(rbc.eq.2 .or. rbc.eq.3) then !wall,sym - dirichlet velocity
         eb1=n-1
      else !outflow,element - neumann velocity
         eb1=n
      endif

!     BCs on B are from tangent vel component
      if(lbc.eq.2) then !wall - dirichlet velocity
         bb0=1
      else !outflow,element,sym - neumann velocity
         bb0=0
      endif
      if(rbc.eq.2) then !wall - dirichlet velocity
         bb1=n-1
      else !outflow,element,sym - neumann velocity
         bb1=n
      endif

      l = (lbc.eq.0)
      r = (rbc.eq.0)

c     calculate E tilde operator
      call set_up_fast_1D_sem_op_cyl(s,eb0,eb1,l,r,ll,lm,lr,
     $                               bh,jgl,dgl,bgl,rho,rad,isd)

c     calculate B tilde operator
      call set_up_fast_1D_sem_mass_cyl(g,bb0,bb1,l,r,ll,lm,lr,
     $                               bh,jgl,bgl,rad,isd)

      n=n+1
      call generalev(s,g,lam,n,w)

!     In cylindrical form, opz is divided by <R>
!     In order to make an approximate tensor product form      
      if (isd.eq.3) then
        rmean = vlsc2(rad,bh,lx1)/2.0
        mri2   = (1.0/rmean)**2
        call cmult(lam,mri2,lx1)
      endif  

      if(.not.l) call row_zero(s,n,n,1)
      if(.not.r) call row_zero(s,n,n,n)
      call transpose(s(n*n+1),n,s,n) ! compute the transpose of s

      return
      end
c-----------------------------------------------------------------------

      subroutine set_up_fast_1D_sem_op_cyl(s,e0,e1,l,r,
     $           ll,lm,lr,bh,jgl,dgl,bgl,rho,radm,isd)

!                  -1 T
!     S = D (rho*B)  D
!
!
!     gives the inexact restriction of this matrix to
!     an element plus one node on either side

      implicit none

      include 'SIZE'

      real s(lx1,lx1)               ! Pseudo Laplacian

      real bh(lx1)                  ! Reference mass matrix
      real jgl(lx2,lx1)             ! Interpolation operator
      real dgl(lx2,lx1)             ! Differential operator
      real bgl(lx2)                 ! Diagonal Mass Matrix on M2
      integer e0,e1                 ! The range for Bhat indices for
                                    ! s (enforces b.c.)
      logical l                     ! If connected to left element
      logical r                     ! If connected to right element
      real ll                       ! Length of left element
      real lm                       ! Length of middle element
      real lr                       ! Length of right element

      real rho(lx1)                 ! density
      real radm(lx1)                ! Radius (middle element) on Mesh 1

      integer isd                   ! Direction

      real bl(lx1)                  ! Mass matrix (inverse) of left element
      real bm(lx1)                  ! Mass matrix (inverse) of middle element
      real br(lx1)                  ! Mass matrix (inverse) of right element

!     Radii Mesh 2      
      real radm2(lx2)               ! Radius (middle element) on Mesh 2
      real radl2(lx2)               ! Radius (left element)   on Mesh 2
      real radr2(lx2)               ! Radius (right element)  on Mesh 2

!     Geometric factors      
      real gl                       ! Geometric factor left 
      real gm                       ! Geometric factor middle
      real gr                       ! Geometric factor right
      real gll                      ! Geometric factor left*left
      real glm                      ! Geometric factor left*middle
      real gmm                      ! Geometric factor middle*middle
      real gmr                      ! Geometric factor middle*right
      real grr                      ! Geometric factor right*right

      common /fast1dwork/ bl,bm,br,radl2,radm2,radr2,
     $                    gl,gm,gr,gll,glm,gmm,gmr,grr 


      integer n
      integer i0,i1,i,j,k

      real sm
      real zr
      real rmean
      real vlsc2

      real d,dt                     ! pointwise values of D,DT


      n=lx1

c     compute the scale factors for J      
      gl=0.0*0.5*ll
      gm=0.0*0.5*lm
      gr=0.0*0.5*lr

      gll = gl*gl
      glm = gl*gm
      gmm = gm*gm
      gmr = gm*gr
      grr = gr*gr

!     Since I have shifted array indicies by 1      
      i0 = e0+1
      i1 = e1+1

!     compute the summed inverse mass matrices for
!     the middle, left, and right elements
      call rzero(bm(1),lx1)
      do i=2,lx1-1
        bm(i) = 2.0/(lm*rho(i)*bh(i))
      enddo
      if (i0.eq.1) then
        if (l) then
          bm(1) = rho(1)*0.5*(ll*bh(lx1) + lm*bh(1))
        else
          bm(1) = 0.5*lm*bh(1)
        endif  
        bm(1)   = 1.0/bm(1)
      endif

      if (i1.eq.lx1) then
        if (r) then
          bm(lx1)= rho(lx1)*0.5*(lr*bh(1) + lm*bh(lx1))
        else
          bm(lx1)= rho(lx1)*0.5*lm*bh(lx1)
        endif
        bm(lx1)  = 1.0/bm(lx1)
      endif

!     note that in computing bl for the left element,
!     bl(1) is missing the contribution from its left neighbor
      if (l) then
        do i=1,lx1-1
          bl(i)=2.0/(ll*rho(i)*bh(i))
        enddo
        bl(lx1)=bm(1)
      endif
!     note that in computing br for the right element,
!     br(n) is missing the contribution from its right neighbor
      if (r) then
        br(1)=bm(lx1)
        do i=2,lx1
          br(i)=2.0/(lr*rho(i)*bh(i))
        enddo
      endif

!     Mass Matrix scaled by R for isd=2
!     I use 's' as a work array here      
      if (isd.eq.2) then
        rmean = vlsc2(radm,bh,lx1)/2.0 
        do i=1,lx1
          bm(i) = bm(i)/radm(i)
          zr    = (radm(i)-rmean)*2.0/lm               ! reference coordinate

!         Scale left element entries            
          s(i,1)  = radm(1) - ((1.0 - zr)*0.5)*ll        ! left element radius(i)
          bl(i)   = bl(i)/s(i,1)

!         Scale right element entries  
          s(i,2)  = radm(lx1) + ((1.0 + zr)*0.5)*lr      ! right element radius(i)
          br(i)   = br(i)/s(i,2)
        enddo
!       Interpolate middle Radius to Mesh 2        
        call mxm(jgl,lx2,radm,lx1,radm2,1)

!       Interpolate Left Radius to Mesh 2        
        call mxm(jgl,lx2,s(1,1),lx1,radl2,1)

!       Interpolate Right Radius to Mesh 2        
        call mxm(jgl,lx2,s(1,2),lx1,radr2,1)

!       jgl is multiplied by mass. We need to remove that factor
        call invcol2(radl2,bgl,lx2)
        call invcol2(radm2,bgl,lx2)
        call invcol2(radr2,bgl,lx2)

      endif  

!     Initialize operator      
      call rzero(s,lx1*lx1)
!     Here we build the interior of the matrix      
      do j=1,lx2
      do i=1,lx2
        sm = 0.0
        do k=i0,i1
          if (isd.eq.1) then
            dt = dgl(j,k)
            d  = dgl(i,k)
          elseif (isd.eq.2) then 
            dt = (radm2(j)*dgl(j,k) + jgl(j,k)*gm)
            d  = (radm2(i)*dgl(i,k) + jgl(i,k)*gm)
          else
            dt = dgl(j,k)
            d  = dgl(i,k)
          endif
!         D*(B^-1)*(D^T)
          sm = sm + d*bm(k)*dt
        enddo
        s(i+1,j+1) = sm
      enddo
      enddo
      
!     Left element contributions      
      if (l) then
        do i=1,lx2
          if (isd.eq.1) then  
            dt = dgl(lx2,lx1)
            d  = dgl(i,1)
          elseif (isd.eq.2) then
            dt = (radl2(lx2)*dgl(lx2,lx1) + jgl(lx2,lx1)*gl)
            d  = (radm2(i)*dgl(i,1)       + jgl(i,1)*gm)
          else
            dt = dgl(lx2,lx1)
            d  = dgl(i,1)
          endif
          s(i+1,1) = d*bm(1)*dt
          s(1,i+1) = s(i+1,1)
        enddo
!       the following is inexact
!       the neighbors bc's are ignored, and the contribution
!       from the neighbor's neighbor is left out
!       that is, bl(0) could be off as noted above
!       or maybe i should go from 1 to n
        do i=1,lx1
          if (isd.eq.1) then
            dt = dgl(lx2,i)
            d  = dgl(lx2,i)
          elseif (isd.eq.2) then
            dt = (radl2(lx2)*dgl(lx2,i) + jgl(lx2,i)*gl)
            d  = (radl2(lx2)*dgl(lx2,i) + jgl(lx2,i)*gl)
          else
            dt = dgl(lx2,i)
            d  = dgl(lx2,i)
          endif
!         D*(B^-1)*(D^T)
          s(1,1) = s(1,1) + d*bl(i)*dt
        enddo
      else
        s(1,1)=1.
      endif

!     Right element contributions      
      if (r) then
        do i=1,lx2
          if (isd.eq.1) then
            dt = dgl(1,1)
            d  = dgl(i,lx1)
          elseif (isd.eq.2) then
            dt = (radr2(1)*dgl(1,1) + jgl(1,1)*gr)
            d  = (radm2(i)*dgl(i,lx1) + jgl(i,lx1)*gm)
          else  
            dt = dgl(1,1)
            d  = dgl(i,lx1)
          endif
          s(i+1,lx1) = d*bm(lx1)*dt
          s(lx1,i+1) = s(i+1,lx1)
        enddo
!       the following is inexact
!       the neighbors bc's are ignored, and the contribution
!       from the neighbor's neighbor is left out
!       that is, br(lx1) could be off as noted above
!       or maybe i should go from 0 to n-1
        do i=1,lx1
          if (isd.eq.1) then
            dt = dgl(1,i)
            d  = dgl(1,i)
          elseif (isd.eq.2) then
!           rad(1)_right = rad(lx1)_middle            
            dt = (radr2(1)*dgl(1,i) + jgl(1,i)*gr)
            d  = (radr2(1)*dgl(1,i) + jgl(1,i)*gr)
          else
            dt = dgl(1,i)
            d  = dgl(1,i)
          endif
!         D*(B^-1)*(D^T)
          s(lx1,lx1) = s(lx1,lx1) + d*br(i)*dt
        enddo
      else
        s(lx1,lx1)=1.
      endif
     
      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_fast_1D_sem_mass_cyl(s,b0,b1,l,r,
     $           ll,lm,lr,bh,jgl,bgl,radm,isd)

c              -1 T
c     S = J (B)  J

c
c     gives the inexact restriction of this matrix to
c     an element plus one node on either side

      implicit none

      include 'SIZE'

      real s(lx1,lx1)               ! Mass Matrix

      real bh(lx1)                  ! Reference mass matrix
      real jgl(lx2,lx1)             ! Interpolation operator
      real bgl(lx2)                 ! Diagonal Mass Matrix on M2
      integer b0,b1                 ! The range for Bhat indices for
                                    ! g (enforces b.c.)
      logical l                     ! If connected to left element
      logical r                     ! If connected to right element
      real ll                       ! Length of left element
      real lm                       ! Length of middle element
      real lr                       ! Length of right element

      real radm(lx1)                ! Radius (middle element) on Mesh 1

      integer isd                   ! Direction

      real bl(lx1)                  ! Mass matrix (inverse) of left element
      real bm(lx1)                  ! Mass matrix (inverse) of middle element
      real br(lx1)                  ! Mass matrix (inverse) of right element

!     Radii Mesh 2      
      real radm2(lx2)               ! Radius (middle element) on Mesh 2
      real radl2(lx2)               ! Radius (left element)   on Mesh 2
      real radr2(lx2)               ! Radius (right element)  on Mesh 2

!     Geometric factors      
      real gl                       ! Geometric factor left 
      real gm                       ! Geometric factor middle
      real gr                       ! Geometric factor right
      real gll                      ! Geometric factor left*left
      real glm                      ! Geometric factor left*middle
      real gmm                      ! Geometric factor middle*middle
      real gmr                      ! Geometric factor middle*right
      real grr                      ! Geometric factor right*right

      common /fast1dwork/ bl,bm,br,radl2,radm2,radr2,
     $                    gl,gm,gr,gll,glm,gmm,gmr,grr 

      integer n
      integer i0,i1,i,j,k

      real sm    
      real zr
      real rmean
      real vlsc2

      real d,dt                     ! pointwise values of J,JT

      n=lx1

c     compute the scale factors for J      
      gl=0.5*ll
      gm=0.5*lm
      gr=0.5*lr

      gll = gl*gl
      glm = gl*gm
      gmm = gm*gm
      gmr = gm*gr
      grr = gr*gr

!     Since I have shifted array indicies by 1      
      i0 = b0+1
      i1 = b1+1

!     compute the summed inverse mass matrices for
!     the middle, left, and right elements
      do i=2,lx1-1
        bm(i) = 2.0/(lm*bh(i))
      enddo
      if (i0.eq.1) then
        if (l) then
          bm(1) = 0.5*(ll*bh(lx1) + lm*bh(1))
        else
          bm(1) = 0.5*lm*bh(1)
        endif  
        bm(1)   = 1.0/bm(1)
      endif

      if (i1.eq.lx1) then
        if (r) then
          bm(lx1)= 0.5*(lr*bh(1) + lm*bh(lx1))
        else
          bm(lx1)= 0.5*lm*bh(lx1)
        endif
        bm(lx1)  = 1.0/bm(lx1)
      endif

!     note that in computing bl for the left element,
!     bl(1) is missing the contribution from its left neighbor
      if (l) then
        do i=1,lx1-1
          bl(i)=2.0/(ll*bh(i))
        enddo
        bl(lx1)=bm(1)
      endif
!     note that in computing br for the right element,
!     br(n) is missing the contribution from its right neighbor
      if (r) then
        br(1)=bm(lx1)
        do i=2,lx1
          br(i)=2.0/(lr*bh(i))
        enddo
      endif

!     Mass Matrix scaled by R for isd=2
!     I use 's' as a work array here      
      if (isd.eq.2) then
        rmean = vlsc2(radm,bh,lx1)/2.0 
        do i=1,lx1
          bm(i) = bm(i)/radm(i)
          zr    = (radm(i)-rmean)*2.0/lm               ! reference coordinate

!         Scale left element entries            
          s(i,1)  = radm(1) - ((1.0 - zr)*0.5)*ll        ! left element radius(i)
          bl(i)   = bl(i)/s(i,1)

!         Scale right element entries  
          s(i,2)  = radm(lx1) + ((1.0 + zr)*0.5)*lr      ! right element radius(i)
          br(i)   = br(i)/s(i,2)
        enddo
!       Interpolate middle Radius to Mesh 2        
        call mxm(jgl,lx2,radm,lx1,radm2,1)

!       Interpolate Left Radius to Mesh 2        
        call mxm(jgl,lx2,s(1,1),lx1,radl2,1)

!       Interpolate Right Radius to Mesh 2        
        call mxm(jgl,lx2,s(1,2),lx1,radr2,1)
      endif  


!     Initialize operator      
      call rzero(s,lx1*lx1)
!     Here we build the interior of the matrix      
      do j=1,lx2
      do i=1,lx2
        sm = 0.0
        do k=i0,i1
          if (isd.eq.1) then
            dt = jgl(j,k)*gm
            d  = jgl(i,k)*gm
          elseif (isd.eq.2) then 
            dt = (radm2(j)*jgl(j,k)*gm)
            d  = (radm2(i)*jgl(i,k)*gm)
          else
            dt = jgl(j,k)*gm
            d  = jgl(i,k)*gm
          endif
!         J*(B^-1)*(J^T)
          sm = sm + d*bm(k)*dt
        enddo
        s(i+1,j+1) = sm
      enddo
      enddo

!     Left element contributions      
      if (l) then
        do i=1,lx2
          if (isd.eq.1) then  
            dt = jgl(lx2,lx1)*gl
            d  = jgl(i,1)*gm
          elseif (isd.eq.2) then
            dt = radl2(lx2)*jgl(lx2,lx1)*gl
            d  = radm2(i)*jgl(i,1)*gm
          else
            dt = jgl(lx2,lx1)*gl
            d  = jgl(i,1)*gm
          endif
          s(i+1,1) = d*bm(1)*dt
          s(1,i+1) = s(i+1,1)
        enddo
!       the following is inexact
!       the neighbors bc's are ignored, and the contribution
!       from the neighbor's neighbor is left out
!       that is, bl(0) could be off as noted above
!       or maybe i should go from 1 to n
        do i=1,lx1
          if (isd.eq.1) then
            dt = jgl(lx2,i)*gl
            d  = jgl(lx2,i)*gl
          elseif (isd.eq.2) then
            dt = radl2(lx2)*jgl(lx2,i)*gl
            d  = radl2(lx2)*jgl(lx2,i)*gl
          else
            dt = jgl(lx2,i)*gl
            d  = jgl(lx2,i)*gl
          endif
!         J*(B^-1)*(J^T)
          s(1,1) = s(1,1) + d*bl(i)*dt
        enddo
      else
        s(1,1)=1.
      endif

!     Right element contributions      
      if (r) then
        do i=1,lx2
          if (isd.eq.1) then
            dt = jgl(1,1)*gr
            d  = jgl(i,lx1)*gm
          elseif (isd.eq.2) then
            dt = radr2(1)*jgl(1,1)*gr
            d  = radm2(i)*jgl(i,lx1)*gm
          else  
            dt = jgl(1,1)*gr
            d  = jgl(i,lx1)*gm
          endif
          s(i+1,lx1) = d*bm(lx1)*dt
          s(lx1,i+1) = s(i+1,lx1)
        enddo
!       the following is inexact
!       the neighbors bc's are ignored, and the contribution
!       from the neighbor's neighbor is left out
!       that is, br(lx1) could be off as noted above
!       or maybe i should go from 0 to n-1
        do i=1,lx1
          if (isd.eq.1) then
            dt = jgl(1,i)*gr
            d  = jgl(1,i)*gr
          elseif (isd.eq.2) then
            dt = radr2(1)*jgl(1,i)*gr
            d  = radr2(1)*jgl(1,i)*gr
          else
            dt = jgl(1,i)*gr
            d  = jgl(1,i)*gr
          endif
!         J*(B^-1)*(J^T)
          s(lx1,lx1) = s(lx1,lx1) + d*br(i)*dt
        enddo
      else
        s(lx1,lx1)=1.
      endif
     
      return
      end
c-----------------------------------------------------------------------






























