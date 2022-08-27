!> @file zneupd.doc.f
!! @ingroup arn_arp
!! @brief Description of zneupd subroutine in arpack 
!! @author Prabal Negi
!! @date June 23, 2021
!
!=======================================================================
!\Name: zneupd 
! 
!\Description: 
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
!  call to ZNAUPD.  ZNAUPD must be called before this routine is called.
!  These approximate eigenvalues and vectors are commonly called Ritz
!  values and Ritz vectors respectively.  They are referred to as such 
!  in the comments that follow.   The computed orthonormal basis for the 
!  invariant subspace corresponding to these Ritz values is referred to as a 
!  Schur basis. 
! 
!  The definition of OP as well as other terms and the relation of computed
!  Ritz values and vectors of OP with respect to the given problem
!  A*z = lambda*B*z may be found in the header of ZNAUPD.  For a brief 
!  description, see definitions of IPARAM(7), MODE and WHICH in the
!  documentation of ZNAUPD.
!
!\Usage:
!  call zneupd 
!     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, WORKEV, BMAT, 
!       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, 
!       WORKL, LWORKL, RWORK, INFO )
!
!\Arguments:
!  RVEC    LOGICAL  (INPUT)
!          Specifies whether a basis for the invariant subspace corresponding
!          to the converged Ritz value approximations for the eigenproblem 
!          A*z = lambda*B*z is computed.
!
!             RVEC = .FALSE.     Compute Ritz values only.
!
!             RVEC = .TRUE.      Compute Ritz vectors or Schur vectors.
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
!          computed. To select the  Ritz vector corresponding to a
!          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
!          If HOWMNY = 'A' or 'P', SELECT need not be initialized 
!          but it is used as internal workspace.
!
!  D       Complex*16 array of dimension NEV+1.  (OUTPUT)
!          On exit, D contains the  Ritz  approximations 
!          to the eigenvalues lambda for A*z = lambda*B*z.
!
!  Z       Complex*16 N by NEV array    (OUTPUT)
!          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of 
!          Z represents approximate eigenvectors (Ritz vectors) corresponding 
!          to the NCONV=IPARAM(5) Ritz values for eigensystem
!          A*z = lambda*B*z.
!
!          If RVEC = .FALSE. or HOWMNY = 'P', then Z is NOT REFERENCED.
!
!          NOTE: If if RVEC = .TRUE. and a Schur basis is not required, 
!          the array Z may be set equal to first NEV+1 columns of the Arnoldi 
!          basis array V computed by ZNAUPD.  In this case the Arnoldi basis 
!          will be destroyed and overwritten with the eigenvector basis.
!
!  LDZ     Integer.  (INPUT)
!          The leading dimension of the array Z.  If Ritz vectors are
!          desired, then  LDZ .ge.  max( 1, N ) is required.  
!          In any case,  LDZ .ge. 1 is required.
!
!  SIGMA   Complex*16  (INPUT)
!          If IPARAM(7) = 3 then SIGMA represents the shift. 
!          Not referenced if IPARAM(7) = 1 or 2.
!
!  WORKEV  Complex*16 work array of dimension 2*NCV.  (WORKSPACE)
!
!  **** The remaining arguments MUST be the same as for the   ****
!  **** call to ZNAUPD that was just completed.               ****
!
!  NOTE: The remaining arguments 
!
!           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, 
!           WORKD, WORKL, LWORKL, RWORK, INFO 
!
!         must be passed directly to ZNEUPD following the last call 
!         to ZNAUPD.  These arguments MUST NOT BE MODIFIED between
!         the the last call to ZNAUPD and the call to ZNEUPD.
!
!  Three of these parameters (V, WORKL and INFO) are also output parameters:
!
!  V       Complex*16 N by NCV array.  (INPUT/OUTPUT)
!
!          Upon INPUT: the NCV columns of V contain the Arnoldi basis
!                      vectors for OP as constructed by ZNAUPD .
!
!          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
!                       contain approximate Schur vectors that span the
!                       desired invariant subspace.
!
!          NOTE: If the array Z has been set equal to first NEV+1 columns
!          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
!          Arnoldi basis held by V has been overwritten by the desired
!          Ritz vectors.  If a separate array Z has been passed then
!          the first NCONV=IPARAM(5) columns of V will contain approximate
!          Schur vectors that span the desired invariant subspace.
!
!  WORKL   Complex*16 work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          WORKL(1:ncv*ncv+2*ncv) contains information obtained in
!          znaupd.  They are not changed by zneupd.
!          WORKL(ncv*ncv+2*ncv+1:3*ncv*ncv+4*ncv) holds the
!          untransformed Ritz values, the untransformed error estimates of 
!          the Ritz values, the upper triangular matrix for H, and the
!          associated matrix representation of the invariant subspace for H.
!
!          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
!          of the above information computed by zneupd.
!          -------------------------------------------------------------
!          IPNTR(9):  pointer to the NCV RITZ values of the
!                     original system.
!          IPNTR(10): Not used
!          IPNTR(11): pointer to the NCV corresponding error estimates.
!          IPNTR(12): pointer to the NCV by NCV upper triangular
!                     Schur matrix for H.
!          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
!                     of the upper Hessenberg matrix H. Only referenced by
!                     zneupd if RVEC = .TRUE. See Remark 2 below.
!          -------------------------------------------------------------
!
!  INFO    Integer.  (OUTPUT)
!          Error flag on output.
!          =  0: Normal exit.
!
!          =  1: The Schur form computed by LAPACK routine csheqr
!                could not be reordered by LAPACK routine ztrsen.
!                Re-enter subroutine zneupd with IPARAM(5)=NCV and
!                increase the size of the array D to have
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
!          = -8: Error return from LAPACK eigenvalue calculation.
!                This should never happened.
!          = -9: Error return from calculation of eigenvectors.
!                Informational error from LAPACK routine ztrevc.
!          = -10: IPARAM(7) must be 1,2,3
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
!          = -12: HOWMNY = 'S' not yet implemented
!          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
!          = -14: ZNAUPD did not find any eigenvalues to sufficient
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
!  3. B. Nour-Omid, B. N. Parlett, T. Ericsson and P. S. Jensen,
!     "How to Implement the Spectral Transformation", Math Comp.,
!     Vol. 48, No. 178, April, 1987 pp. 664-673. 
!
!\Routines called:
!     ivout   ARPACK utility routine that prints integers.
!     zmout   ARPACK utility routine that prints matrices
!     zvout   ARPACK utility routine that prints vectors.
!     zgeqr2  LAPACK routine that computes the QR factorization of 
!             a matrix.
!     zlacpy  LAPACK matrix copy routine.
!     zlahqr  LAPACK routine that computes the Schur form of a
!             upper Hessenberg matrix.
!     zlaset  LAPACK matrix initialization routine.
!     ztrevc  LAPACK routine to compute the eigenvectors of a matrix
!             in upper triangular form.
!     ztrsen  LAPACK routine that re-orders the Schur form.
!     zunm2r  LAPACK routine that applies an orthogonal matrix in 
!             factored form.
!     dlamch  LAPACK routine that determines machine constants.
!     ztrmm   Level 3 BLAS matrix times an upper triangular matrix.
!     zgeru   Level 2 BLAS rank one update to a matrix.
!     zcopy   Level 1 BLAS that copies one vector to another .
!     zscal   Level 1 BLAS that scales a vector.
!     zdscal  Level 1 BLAS that scales a complex vector by a real number.
!     dznrm2  Level 1 BLAS that computes the norm of a complex vector.
!
!\Remarks
!
!  1. Currently only HOWMNY = 'A' and 'P' are implemented. 
!
!  2. Schur vectors are an orthogonal representation for the basis of
!     Ritz vectors. Thus, their numerical properties are often superior.
!     If RVEC = .true. then the relationship
!             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
!     V(:,1:IPARAM(5))' * V(:,1:IPARAM(5)) = I are approximately satisfied.
!     Here T is the leading submatrix of order IPARAM(5) of the 
!     upper triangular matrix stored workl(ipntr(12)). 
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Chao Yang                    Houston, Texas 
!     Dept. of Computational & 
!     Applied Mathematics 
!     Rice University 
!     Houston, Texas
!

