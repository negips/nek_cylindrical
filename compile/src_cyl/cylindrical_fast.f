!====================================================================== 
!     Author: Prabal Negi
!     Description: Resetting the FDM preconditioner for cylindrical
!                  solve.
!      
!      
!====================================================================== 
      subroutine set_up_fast_1D_sem_cyl(s,lam,n,lbc,rbc,ll,lm,lr,
     $                                    rho,ie,isd)

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
      real rho(lx1),dummy(lx1)
      real rhobh(lx1)     ! density*mass

      integer isd

      n=lx1
      call col3(rhobh,bh,rho,n)
     
      n=lx1-1
      !bcs on E are from normal vel component
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
      !bcs on B are from tangent vel component
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
c
      l = (lbc.eq.0)
      r = (rbc.eq.0)
c
c     calculate E tilde operator
      call set_up_fast_1D_sem_op_cyl(s,eb0,eb1,l,r,ll,lm,lr,bh,dgl,
     $                                 rho,0)

c     call outmat(s,n+1,n+1,'  Et  ',ie)
c     calculate B tilde operator
      call set_up_fast_1D_sem_cyl(g,bb0,bb1,l,r,ll,lm,lr,bh,jgl,1)
     
      n=n+1
      call generalev(s,g,lam,n,w)
      if(.not.l) call row_zero(s,n,n,1)
      if(.not.r) call row_zero(s,n,n,n)
      call transpose(s(n*n+1),n,s,n) ! compute the transpose of s

      return
      end
c-----------------------------------------------------------------------

      subroutine set_up_fast_1D_sem_op_cyl(s,e0,e1,l,r,
     $           ll,lm,lr,bh,jgl,dgl,jgl,rho,rad)

c                  -1 T
c     S = D (rho*B)  D
c
c                  -1 T
c     G = J (rho*B)  J

c
c     gives the inexact restriction of this matrix to
c     an element plus one node on either side

!     rho - density

      implicit none

      include 'SIZE'

      real s(lx1,lx1)               ! Pseudo Laplacian

      real bh(lx1)                  ! Reference mass matrix
      real jgl(lx2,lx1)             ! Interpolation operator
      real dgl(lx2,lx1)             ! Differential operator
      real ll                       ! Length of left element
      real lm                       ! Length of current element
      real lr                       ! Length of right element
      integer e0,e1                 ! The range for Bhat indices for
                                    ! s (enforces b.c.)
      integer b0,b1                 ! The range for Bhat indices for
                                    ! g (enforces b.c.)
      logical l                     ! If connected to left element
      logical r                     ! If connected to right element

      real bl(lx1)                  ! Mass matrix (inverse) of left element
      real bm(lx1)                  ! Mass matrix (inverse) of current element
      real br(lx1)                  ! Mass matrix (inverse) of right element
!     Geometric factors      
      real gl                       ! Geometric factor left 
      real gm                       ! Geometric factor middle/current
      real gr                       ! Geometric factor right
      real gll                      ! Geometric factor left*left
      real glm                      ! Geometric factor left*middle
      real gmm                      ! Geometric factor middle*middle
      real gmr                      ! Geometric factor middle*right
      real grr                      ! Geometric factor right*right
      real rho(lx1)                 ! density
      real rhobh(lx1)               ! density*mass
      real rad(lx1)                 ! Radius

      integer n
      integer i0,i1

      real sm    

      real d,dt                     ! pointwise values of D,DT
      real j,jt                     ! pointwise values of J,JT

      n=lx1

!     (middle element density)*(reference mass)
      call col3(rhobh,bh,rho,lx1)

c     compute the scale factors for J      
      gl=0.5*ll
      gm=0.5*lm
      gr=0.5*lr

      gll = gl*gl
      glm = gl*gm
      gmm = gm*gm
      gmr = gm*gr
      grr = gr*gr

!     compute the summed inverse mass matrices for
!     the middle, left, and right elements
      do i=2,lx1-1
        bm(i) = 2.0/(lm*rhobh(i))
      enddo
      if (e0.eq.0) then
        if (l) then
          bm(1) = 0.5*(ll*rhobh(lx1)+lm*rhobh(1))
        else
          bm(1) = 0.5*lm*rhobh(1)
        endif  
        bm(1)   = 1.0/bm(1)
      endif

      if (e1.eq.(lx1-1)) then
        bm(lx1) = 0.5*lm*rhobh(n)
        if (r) then
          bm(lx1)= 0.5*(lr*rhobh(1)+lm*rhobh(lx1))
        else
          bm(lx1)= 0.5*lm*rhobh(lx1)
        endif
        bm(lx1)  =1.0/bm(lx1)
      endif

!     note that in computing bl for the left element,
!     bl(1) is missing the contribution from its left neighbor
      if (l) then
        do i=1,lx1-1
          bl(i)=2.0/(ll*rhobh(i))
        enddo
        bl(lx1)=bm(1)
      endif
!     note that in computing br for the right element,
!     br(n) is missing the contribution from its right neighbor
      if (r) then
        br(1)=bm(lx1)
        do i=2,lx1
          br(i)=2.0/(lr*rhobh(i))
        enddo
      endif

!     Initialize operator      
      call rzero(s,lx1*lx1)
!     Since I have shifted array indicies by 1      
      i0 = e0+1
      i1 = e1+1
!     Here we build the interior of the matrix      
      do j=1,lx2
      do i=1,lx2
        sm = 0.0
        do k=i0,i1
          if (isd.eq.1) then
            dt = dgl(j,k)
            d  = dgl(i,k)
          elseif (isd.eq.2) then 
            dt = (rad(j)*dgl(j,k) + jgl(j,k)*gm)
            d  = (rad(i)*dgl(i,k) + jgl(i,k)*gm)
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
      
      if (l) then
        do i=2,lx1-1
          if (isd.eq.1) then  
            dt = dgl(lx2,lx1)
            d  = dgl(i,1)
          elseif (isd.eq.2) then
            dt = (rad(lx2)*dgl(lx2,lx1) + jgl(lx2,lx1)*gl)
            d  = (rad(i)*dgl(i,1) + jgl(i,1)*gm)
          else
            dt = dgl(lx2,lx1)
            d  = dgl(i,1)
          endif
          s(i,1) = d*bm(1)*dt
          s(1,i) = s(i,1)
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
            dt = (rad(lx2)*dgl(lx2,i) + jgl(lx2,i)*gl)
            d  = (rad(lx2)*dgl(lx2,i) + jgl(lx2,i)*gl)
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
      
      if (r) then
        do i=2,lx1-1
          if (isd.eq.1) then
            dt = dgl(1,1)
            d  = dgl(i,lx1)
          elseif (isd.eq.2) then
            dt = (rad(1)*dgl(1,1) + jgl(1,1)*gr)
            d  = (rad(i)*dgl(i,lx1) + jgl(i,lx1)*gm)
          else  
            dt = dgl(1,1)
            d  = dgl(i,lx1)
          endif
          s(i,lx1) = d*bm(lx1)*dt
          s(lx1,i) = s(i,lx1)
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
            dt = (rad(1)*dgl(1,i) + jgl(1,i)*gr)
            d  = (rad(1)*dgl(1,i) + jgl(1,i)*gr)
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
     $           ll,lm,lr,bh,jgl,dgl,jgl,rho,rad)

c              -1 T
c     S = J (B)  J

c
c     gives the inexact restriction of this matrix to
c     an element plus one node on either side

!     rho - density

      implicit none

      include 'SIZE'

      real s(lx1,lx1)               ! Mass Matrix

      real bh(lx1)                  ! Reference mass matrix
      real jgl(lx2,lx1)             ! Interpolation operator
      real dgl(lx2,lx1)             ! Differential operator
      real ll                       ! Length of left element
      real lm                       ! Length of current element
      real lr                       ! Length of right element

      integer b0,b1                 ! The range for Bhat indices for
                                    ! g (enforces b.c.)
      logical l                     ! If connected to left element
      logical r                     ! If connected to right element

      real bl(lx1)                  ! Mass matrix (inverse) of left element
      real bm(lx1)                  ! Mass matrix (inverse) of current element
      real br(lx1)                  ! Mass matrix (inverse) of right element
!     Geometric factors      
      real gl                       ! Geometric factor left 
      real gm                       ! Geometric factor middle/current
      real gr                       ! Geometric factor right
      real gll                      ! Geometric factor left*left
      real glm                      ! Geometric factor left*middle
      real gmm                      ! Geometric factor middle*middle
      real gmr                      ! Geometric factor middle*right
      real grr                      ! Geometric factor right*right
      real rho(lx1)                 ! density
      real rhobh(lx1)               ! density*mass
      real rad(lx1)                 ! Radius

      integer n
      integer i0,i1

      real sm    

      real d,dt                     ! pointwise values of D,DT
      real j,jt                     ! pointwise values of J,JT

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
          bm(1) = 0.5*(ll*bh(lx1)+lm*bh(1))
        else
          bm(1) = 0.5*lm*bh(1)
        endif  
        bm(1)   = 1.0/bm(1)
      endif

      if (i1.eq.lx1) then
        bm(lx1) = 0.5*lm*bh(n)
        if (r) then
          bm(lx1)= 0.5*(lr*bh(1)+lm*bh(lx1))
        else
          bm(lx1)= 0.5*lm*bh(lx1)
        endif
        bm(lx1)  =1.0/bm(lx1)
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
            dt = (rad(j)*dgl(j,k) + jgl(j,k)*gm)
            d  = (rad(i)*dgl(i,k) + jgl(i,k)*gm)
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
      
      if (l) then
        do i=2,lx1-1
          if (isd.eq.1) then  
            dt = dgl(lx2,lx1)
            d  = dgl(i,1)
          elseif (isd.eq.2) then
            dt = (rad(lx2)*dgl(lx2,lx1) + jgl(lx2,lx1)*gl)
            d  = (rad(i)*dgl(i,1) + jgl(i,1)*gm)
          else
            dt = dgl(lx2,lx1)
            d  = dgl(i,1)
          endif
          s(i,1) = d*bm(1)*dt
          s(1,i) = s(i,1)
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
            dt = (rad(lx2)*dgl(lx2,i) + jgl(lx2,i)*gl)
            d  = (rad(lx2)*dgl(lx2,i) + jgl(lx2,i)*gl)
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
      
      if (r) then
        do i=2,lx1-1
          if (isd.eq.1) then
            dt = dgl(1,1)
            d  = dgl(i,lx1)
          elseif (isd.eq.2) then
            dt = (rad(1)*dgl(1,1) + jgl(1,1)*gr)
            d  = (rad(i)*dgl(i,lx1) + jgl(i,lx1)*gm)
          else  
            dt = dgl(1,1)
            d  = dgl(i,lx1)
          endif
          s(i,lx1) = d*bm(lx1)*dt
          s(lx1,i) = s(i,lx1)
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
            dt = (rad(1)*dgl(1,i) + jgl(1,i)*gr)
            d  = (rad(1)*dgl(1,i) + jgl(1,i)*gr)
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

