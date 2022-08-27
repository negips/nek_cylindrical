!> @file dneupd.doc.f
!! @ingroup arn_arp
!! @brief Description of dneupd subroutine in arpack 
!! @author Prabal Negi
!! @date June 23, 2021
!
!=======================================================================
!\Name: dneupd
!
!\Description: 
!
!  This subroutine returns the converged approximations to eigenvalues
!  of A*z = lambda*B*z and (optionally):
!
!      (1) The corresponding approximate eigenvectors;
!
!      (2) An orthonormal basis for the associated approximate
!          invariant subspace;
!
!      (3) Both.
!
!  There is negligible additional cost to obtain eigenvectors.  An orthonormal
!  basis is always computed.  There is an additional storage cost of n*nev
!  if both are requested (in this case a separate array Z must be supplied).
!
!  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
!  are derived from approximate eigenvalues and eigenvectors of
!  of the linear operator OP prescribed by the MODE selection in the
!  call to DNAUPD.  DNAUPD must be called before this routine is called.
!  These approximate eigenvalues and vectors are commonly called Ritz
!  values and Ritz vectors respectively.  They are referred to as such
!  in the comments that follow.  The computed orthonormal basis for the
!  invariant subspace corresponding to these Ritz values is referred to as a
!  Schur basis.
!
!  See documentation in the header of the subroutine DNAUPD for 
!  definition of OP as well as other terms and the relation of computed
!  Ritz values and Ritz vectors of OP with respect to the given problem
!  A*z = lambda*B*z.  For a brief description, see definitions of 
!  IPARAM(7), MODE and WHICH in the documentation of DNAUPD.
!
!\Usage:
!  call dneupd 
!     ( RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI, WORKEV, BMAT, 
!       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, 
!       LWORKL, INFO )
!
!\Arguments:
!  RVEC    LOGICAL  (INPUT) 
!          Specifies whether a basis for the invariant subspace corresponding 
!          to the converged Ritz value approximations for the eigenproblem 
!          A*z = lambda*B*z is computed.
!
!             RVEC = .FALSE.     Compute Ritz values only.
!
!             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.
!                                See Remarks below. 
! 
!  HOWMNY  Character*1  (INPUT) 
!          Specifies the form of the basis for the invariant subspace 
!          corresponding to the converged Ritz values that is to be computed.
!
!          = 'A': Compute NEV Ritz vectors; 
!          = 'P': Compute NEV Schur vectors;
!          = 'S': compute some of the Ritz vectors, specified
!                 by the logical array SELECT.
!
!  SELECT  Logical array of dimension NCV.  (INPUT)
!          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
!          computed. To select the Ritz vector corresponding to a
!          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE.. 
!          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace.
!
!  DR      Double precision array of dimension NEV+1.  (OUTPUT)
!          If IPARAM(7) = 1,2 or 3 and SIGMAI=0.0  then on exit: DR contains 
!          the real part of the Ritz  approximations to the eigenvalues of 
!          A*z = lambda*B*z. 
!          If IPARAM(7) = 3, 4 and SIGMAI is not equal to zero, then on exit:
!          DR contains the real part of the Ritz values of OP computed by 
!          DNAUPD. A further computation must be performed by the user
!          to transform the Ritz values computed for OP by DNAUPD to those
!          of the original system A*z = lambda*B*z. See remark 3 below.
!
!  DI      Double precision array of dimension NEV+1.  (OUTPUT)
!          On exit, DI contains the imaginary part of the Ritz value 
!          approximations to the eigenvalues of A*z = lambda*B*z associated
!          with DR.
!
!          NOTE: When Ritz values are complex, they will come in complex 
!                conjugate pairs.  If eigenvectors are requested, the 
!                corresponding Ritz vectors will also come in conjugate 
!                pairs and the real and imaginary parts of these are 
!                represented in two consecutive columns of the array Z 
!                (see below).
!
!  Z       Double precision N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. (OUTPUT)
!          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of 
!          Z represent approximate eigenvectors (Ritz vectors) corresponding 
!          to the NCONV=IPARAM(5) Ritz values for eigensystem 
!          A*z = lambda*B*z. 
! 
!          The complex Ritz vector associated with the Ritz value 
!          with positive imaginary part is stored in two consecutive 
!          columns.  The first column holds the real part of the Ritz 
!          vector and the second column holds the imaginary part.  The 
!          Ritz vector associated with the Ritz value with negative 
!          imaginary part is simply the complex conjugate of the Ritz vector 
!          associated with the positive imaginary part.
!
!          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
!
!          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
!          the array Z may be set equal to first NEV+1 columns of the Arnoldi
!          basis array V computed by DNAUPD.  In this case the Arnoldi basis
!          will be destroyed and overwritten with the eigenvector basis.
!
!  LDZ     Integer.  (INPUT)
!          The leading dimension of the array Z.  If Ritz vectors are
!          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.
!
!  SIGMAR  Double precision  (INPUT)
!          If IPARAM(7) = 3 or 4, represents the real part of the shift. 
!          Not referenced if IPARAM(7) = 1 or 2.
!
!  SIGMAI  Double precision  (INPUT)
!          If IPARAM(7) = 3 or 4, represents the imaginary part of the shift. 
!          Not referenced if IPARAM(7) = 1 or 2. See remark 3 below.
!
!  WORKEV  Double precision work array of dimension 3*NCV.  (WORKSPACE)
!
!  **** The remaining arguments MUST be the same as for the   ****
!  **** call to DNAUPD that was just completed.               ****
!
!  NOTE: The remaining arguments
!
!           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
!           WORKD, WORKL, LWORKL, INFO
!
!         must be passed directly to DNEUPD following the last call
!         to DNAUPD.  These arguments MUST NOT BE MODIFIED between
!         the the last call to DNAUPD and the call to DNEUPD.
!
!  Three of these parameters (V, WORKL, INFO) are also output parameters:
!
!  V       Double precision N by NCV array.  (INPUT/OUTPUT)
!
!          Upon INPUT: the NCV columns of V contain the Arnoldi basis
!                      vectors for OP as constructed by DNAUPD .
!
!          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
!                       contain approximate Schur vectors that span the
!                       desired invariant subspace.  See Remark 2 below.
!
!          NOTE: If the array Z has been set equal to first NEV+1 columns
!          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
!          Arnoldi basis held by V has been overwritten by the desired
!          Ritz vectors.  If a separate array Z has been passed then
!          the first NCONV=IPARAM(5) columns of V will contain approximate
!          Schur vectors that span the desired invariant subspace.
!
!  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          WORKL(1:ncv*ncv+3*ncv) contains information obtained in
!          dnaupd.  They are not changed by dneupd.
!          WORKL(ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) holds the
!          real and imaginary part of the untransformed Ritz values,
!          the upper quasi-triangular matrix for H, and the
!          associated matrix representation of the invariant subspace for H.
!
!          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
!          of the above information computed by dneupd.
!          -------------------------------------------------------------
!          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
!                     original system.
!          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
!                     the original system.
!          IPNTR(11): pointer to the NCV corresponding error bounds.
!          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
!                     Schur matrix for H.
!          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
!                     of the upper Hessenberg matrix H. Only referenced by
!                     dneupd if RVEC = .TRUE. See Remark 2 below.
!          -------------------------------------------------------------
!
!  INFO    Integer.  (OUTPUT)
!          Error flag on output.
!
!          =  0: Normal exit.
!
!          =  1: The Schur form computed by LAPACK routine dlahqr
!                could not be reordered by LAPACK routine dtrsen.
!                Re-enter subroutine dneupd with IPARAM(5)=NCV and 
!                increase the size of the arrays DR and DI to have 
!                dimension at least dimension NCV and allocate at least NCV 
!                columns for Z. NOTE: Not necessary if Z and V share 
!                the same space. Please notify the authors if this error
!                occurs.
!
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV-NEV >= 2 and less than or equal to N.
!          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work WORKL array is not sufficient.
!          = -8: Error return from calculation of a real Schur form.
!                Informational error from LAPACK routine dlahqr.
!          = -9: Error return from calculation of eigenvectors.
!                Informational error from LAPACK routine dtrevc.
!          = -10: IPARAM(7) must be 1,2,3,4.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
!          = -12: HOWMNY = 'S' not yet implemented
!          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
!          = -14: DNAUPD did not find any eigenvalues to sufficient
!                 accuracy.
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
!     Real Matrices", Linear Algebra and its Applications, vol 88/89,
!     pp 575-595, (1987).
!
!\Routines called:
!     ivout   ARPACK utility routine that prints integers.
!     dmout   ARPACK utility routine that prints matrices
!     dvout   ARPACK utility routine that prints vectors.
!     dgeqr2  LAPACK routine that computes the QR factorization of 
!             a matrix.
!     dlacpy  LAPACK matrix copy routine.
!     dlahqr  LAPACK routine to compute the real Schur form of an
!             upper Hessenberg matrix.
!     dlamch  LAPACK routine that determines machine constants.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     dlaset  LAPACK matrix initialization routine.
!     dorm2r  LAPACK routine that applies an orthogonal matrix in 
!             factored form.
!     dtrevc  LAPACK routine to compute the eigenvectors of a matrix
!             in upper quasi-triangular form.
!     dtrsen  LAPACK routine that re-orders the Schur form.
!     dtrmm   Level 3 BLAS matrix times an upper triangular matrix.
!     dger    Level 2 BLAS rank one update to a matrix.
!     dcopy   Level 1 BLAS that copies one vector to another .
!     ddot    Level 1 BLAS that computes the scalar product of two vectors.
!     dnrm2   Level 1 BLAS that computes the norm of a vector.
!     dscal   Level 1 BLAS that scales a vector.
!
!\Remarks
!
!  1. Currently only HOWMNY = 'A' and 'P' are implemented.
!
!     Let X' denote the transpose of X.
!
!  2. Schur vectors are an orthogonal representation for the basis of
!     Ritz vectors. Thus, their numerical properties are often superior.
!     If RVEC = .TRUE. then the relationship
!             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
!     V(:,1:IPARAM(5))' * V(:,1:IPARAM(5)) = I are approximately satisfied.
!     Here T is the leading submatrix of order IPARAM(5) of the real 
!     upper quasi-triangular matrix stored workl(ipntr(12)). That is,
!     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; 
!     each 2-by-2 diagonal block has its diagonal elements equal and its
!     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
!     diagonal block is a complex conjugate pair of Ritz values. The real
!     Ritz values are stored on the diagonal of T.
!
!  3. If IPARAM(7) = 3 or 4 and SIGMAI is not equal zero, then the user must
!     form the IPARAM(5) Rayleigh quotients in order to transform the Ritz
!     values computed by DNAUPD for OP to those of A*z = lambda*B*z. 
!     Set RVEC = .true. and HOWMNY = 'A', and
!     compute 
!           Z(:,I)' * A * Z(:,I) if DI(I) = 0.
!     If DI(I) is not equal to zero and DI(I+1) = - D(I), 
!     then the desired real and imaginary parts of the Ritz value are
!           Z(:,I)' * A * Z(:,I) +  Z(:,I+1)' * A * Z(:,I+1),
!           Z(:,I)' * A * Z(:,I+1) -  Z(:,I+1)' * A * Z(:,I), respectively.
!     Another possibility is to set RVEC = .true. and HOWMNY = 'P' and
!     compute V(:,1:IPARAM(5))' * A * V(:,1:IPARAM(5)) and then an upper
!     quasi-triangular matrix of order IPARAM(5) is computed. See remark
!     2 above.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University 
!     Chao Yang                    Houston, Texas
!     Dept. of Computational &
!     Applied Mathematics          
!     Rice University           
!     Houston, Texas            

