!======================================================================
!     Author: Prabal Negi      
!     Description: Algorithms from the Book
!                  Implementing Spectral Methods for 
!                  Partial Differential Equations: 
!                  Algorithms for Scientists and Engineers (2009)       
!                  - David A. Kopriva        
!     Routines: BarycentricWeights              : Algorithm 30
!               LagrangeInterpolation           : Algorithm 31      
!               PolynomialInterpolationMatrix   : Algorithm 32
!               LagrangeInterpolatingPolynomial : Algorithm 34     
!               LagrangeInterpolantDerivative   : Algorithm 36
!               LagrangeDerivativeMatrix        : Custom        
!               PolynomialDerivativeMatrix      : Algorithm 37      
!               HigherOrdPolyDerivMatrix        : Algorithm 38
!
!
!     Function: AlmostEqual                     : Algorithm 131
!
!====================================================================== 
!---------------------------------------------------------------------- 
      function AlmostEqual(a,b)

!     Algorithm 131
      implicit none

      real a,b
      logical AlmostEqual
      real eps

      eps = epsilon(a)

      if ((a.eq.0.0).or.(b.eq.0.0)) then
        if (abs(a-b).le.2*eps) then
          AlmostEqual = .true.
        else
          AlmostEqual = .false.
        endif
      else
        if ((abs(a-b).le.eps*abs(a)).and.
     $      (abs(a-b).le.eps*abs(b))) then
          AlmostEqual = .true.
        else
          AlmostEqual = .false.
        endif
      endif  

      return
      end function AlmostEqual
!---------------------------------------------------------------------- 

      subroutine BarycentricWeights(w,x,n)

!     Algorithm 30: BarycentricWeights: Weights for Lagrange Interpolation
      implicit none

      integer n   ! No of Points
      real x(n)   ! Locations of Lagrange Interpolating points
      real w(n)   ! Baryweights

      integer j,k

      do j=1,n
        w(j) = 1.0
      enddo  

      do j=2,n
      do k=1,j-1
        w(k) = w(k)*(x(k) - x(j))
        w(j) = w(j)*(x(j) - x(k))
      enddo
      enddo  

      do j=1,n
        w(j) = 1.0/w(j)
      enddo  

      return
      end subroutine BarycentricWeights
!---------------------------------------------------------------------- 
      subroutine LagrangeInterpolation(fx,f,xi,x,w,n)

!     Algorithm 31: LagrangeInterpolation        
      implicit none

      real xi     ! input point
      real fx     ! output value f(xi)

      integer n   ! no of points in the interpolating polynomial
      real x(n)   ! Lagrange interpolating points
      real f(n)   ! Function values at x
      real w(n)   ! Baryweights

      logical AlmostEqual     ! function
      real t,num,den

      integer j

      num = 0.0
      den = 0.0

      do j=1,n
        if (AlmostEqual(xi,x(j))) then
          fx = f(j)
          return
        endif
        t   = w(j)/(xi - x(j))
        num = num + t*f(j)
        den = den + t
      enddo

      fx = num/den

      return
      end subroutine LagrangeInterpolation
!----------------------------------------------------------------------       
      subroutine LagrangeInterpolatingPolynomial
     $           (lx,xi,x,w,n)

!     Algorithm 34: LagrangeInterpolatingPolynomials

      implicit none

!     Input:
      integer n   ! no of points in the interpolating polynomial

      real xi     ! input point
      real x(n)   ! Lagrange interpolating points
      real w(n)   ! Baryweights

!     Output:      
      real lx(n)  ! Interpolating polynomial (coefficients)

      logical AlmostEqual     ! function
      logical xMatchesNode
      real t,s

      integer j

      xMatchesNode = .false.
      do j=1,n
        lx(j) = 0.0
        if (AlmostEqual(xi,x(j))) then
          lx(j) = 1.0
          xMatchesNode = .true.
        endif
      enddo

      if (xMatchesNode) return

      s = 0.0
      do j=1,n
        t     = w(j)/(xi - x(j))
        lx(j) = t
        s = s + t
      enddo

      do j=1,n
        lx(j) = lx(j)/s
      enddo  

      return
      end subroutine LagrangeInterpolatingPolynomial
!----------------------------------------------------------------------       

      subroutine PolynomialInterpolationMatrix
     $           (T,y,m,x,w,n)

!     Algorithm 32: PolynomialInterpolationMatrix

      implicit none

      integer m   ! no of points to interpolate to
      real y(m)   ! input locations

      integer n
      real x(n)   ! Lagrange interpolating points
      real w(n)   ! Baryweights

      real T(m,n) ! Interpolation matrix

      logical AlmostEqual     ! function
      logical rowhasmatch

      integer j,k
      real s,r

!     Initialize
      do j=1,n
      do k=1,m
        T(k,j) = 0.0
      enddo
      enddo

      do k=1,m
        rowhasmatch = .false.
        do j=1,n
          if (AlmostEqual(y(k),x(j))) then
            rowhasmatch = .true.
            T(k,j) = 1.0
          endif  
        enddo

        if (.not.rowhasmatch) then
          s = 0.0
          do j=1,n
            r      = w(j)/(y(k) - x(j))
            s      = s + r
            T(k,j) = r
          enddo
          do j=1,n
            T(k,j) = T(k,j)/s
          enddo
        endif  
      enddo

      return
      end subroutine 
!----------------------------------------------------------------------  

      subroutine LagrangeInterpolantDerivative
     $           (df,f,xi,x,w,n)         

!     Algorithm 36: LagrangeInterpolatDerivative:
!     Direct Computation of the Polynomial Derivative 
!     in Barycentric form        
      implicit none

      real xi     ! input point
      real df     ! output value df/dx(xi)

      integer n   ! no of points in the interpolating polynomial
      real x(n)   ! Lagrange interpolating points
      real f(n)   ! Function values at x
      real w(n)   ! Baryweights

      logical AlmostEqual     ! function
      logical atNode
      real t,num,den

      real p      ! f(xi) using LagrangeInterpolation 


      integer i,j

      atNode = .false.
      num = 0.0

      do j=1,n
        if (AlmostEqual(xi,x(j))) then
          atNode = .true.
          p   = f(j)
          den = -w(j)
          i   = j  
        endif
      enddo

      if (atNode) then
        do j=1,n
          if (j.ne.i) then
            num = num + w(j)*(p-f(j))/(xi - x(j))
          endif
        enddo
      else
        den = 0.0  
        call LagrangeInterpolation(p,f,xi,x,w,n)
        do j=1,n
          t = w(j)/(xi - x(j))
          num = num + t*(p-f(j))/(xi - x(j))
          den = den + t
        enddo
      endif

      df = num/den  

      return
      end subroutine
!----------------------------------------------------------------------       

      subroutine LagrangeDerivativeMatrix
     $           (D,y,m,x,w,n)         

!     Derivative matrix at points y(m) via the 
!     Lagrange interpolating polynomial at x(n)

!     Custom algorithm built from Eqn. (3.45) and (3.46)
!     in Kopriva (2006) Implementing Spectral Methods for PDEs
      implicit none

!     Input      
      integer m   ! no of derivative points
      real y(m)   ! Points of derivative evaluation

      integer n   ! no of points in the interpolating polynomial
      real x(n)   ! Lagrange interpolating points
      real w(n)   ! Baryweights

!     Output      
      real D(m,n)

      real lx(n)  ! Lagrange Interpolating Polynomial

      logical AlmostEqual     ! function
      logical atNode
      real t,s,s2

      real p      ! f(xi) using LagrangeInterpolation 

      integer i,j,k


!     Initialize
      do j=1,n
      do k=1,m
        D(k,j) = 0.0
      enddo
      enddo

      do k=1,m
        atNode = .false.
        do j=1,n
          if (AlmostEqual(y(k),x(j))) then
            atNode = .true.
            s2  = -w(j)
            i   = j  
          endif
        enddo

        if (atNode) then
          s = 0.0
          do j=1,n
            if (j.ne.i) then
              t = w(j)/(y(k) - x(j))
              s = s + t 
              D(k,j) = -t
            endif
          enddo
          D(k,i) = s
          do j=1,n
            D(k,j) = D(k,j)/s2
          enddo
        else
          call LagrangeInterpolatingPolynomial
     $           (lx,y(k),x,w,n)
          s  = 0.0
          s2 = 0.0
          do j=1,n
            t      = w(j)/(y(k) - x(j))
            D(k,j) = -t/(y(k) - x(j))
            s      = s + t
            s2     = s2 + t/(y(k) - x(j))
          enddo
            
          do j=1,n
            D(k,j) = D(k,j) + s2*lx(j)
            D(k,j) = D(k,j)/s
          enddo  
          
        endif

      enddo  

      return
      end subroutine
!----------------------------------------------------------------------       

      subroutine PolynomialDerivativeMatrix
     $           (D,x,w,n)

!     Algorithm 37: PolynomialDerivativeMatrix
!     Derivative matrix for the Nodal points.        

      implicit none

!     Input      
      integer n   ! No of Lagrange Interpolating points
      real x(n)   ! Lagrange interpolating points
      real w(n)   ! Baryweights

!     Output      
      real D(n,n) ! Derivative matrix

      integer i,j

!     Initialize
      do j=1,n
      do i=1,n
        D(i,j) = 0.0
      enddo
      enddo

      do i=1,n
        D(i,i) = 0.0
        do j=1,n
          if (j.ne.i) then
            D(i,j) = w(j)/(w(i)*(x(i) - x(j)))
            D(i,i) = D(i,i) - D(i,j)
          endif  
        enddo
      enddo

      return
      end subroutine 
!----------------------------------------------------------------------  
     
      subroutine HigherOrdPolyDerivMatrix
     $           (D,x,w,n,m)

!     Algorithm 38: PolynomialDerivativeMatrix
!     Derivative matrix for the Nodal points.        

      implicit none

!     Input      
      integer m   ! Order of derivatives required
      integer n   ! No of Lagrange Interpolating points
      real x(n)   ! Lagrange interpolating points
      real w(n)   ! Baryweights

!     Output      
      real D(n,n,m) ! Derivative matrix

      real t

      integer i,j,k

!     Initialize
      do k=1,m
      do j=1,n
      do i=1,n
        D(i,j,k) = 0.0
      enddo
      enddo
      enddo

      call PolynomialDerivativeMatrix(D(1,1,k),x,w,n)

      if (m.eq.1) return

      do k=2,m 
        do i=1,n
          D(i,i,k) = 0.0
          do j=1,n
            if (j.ne.i) then
              t        = k/(x(i) - x(j))
              D(i,j,k) = t*(w(j)*D(i,i,k-1)/w(i) - D(i,j,k-1))
              D(i,i,k) = D(i,i,k) - D(i,j,k)
            endif  
          enddo   ! i
        enddo     ! j
      enddo       ! k

      return
      end subroutine 
!----------------------------------------------------------------------  










