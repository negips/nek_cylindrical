!> @file dnaupd.doc.f
!! @ingroup arn_arp
!! @brief Description of dnaupd subroutine in arpack 
!! @warning There is no restart option for serial ARPACK version. It is
!!   supported by parallel PARPACK only.
!! @author Prabal Negi
!! @date June 23, 2021
!
!=======================================================================
!\Name: dnaupd
!
!\Description: 
!  Reverse communication interface for the Implicitly Restarted Arnoldi
!  iteration. This subroutine computes approximations to a few eigenpairs 
!  of a linear operator "OP" with respect to a semi-inner product defined by 
!  a symmetric positive semi-definite real matrix B. B may be the identity 
!  matrix. NOTE: If the linear operator "OP" is real and symmetric 
!  with respect to the real positive semi-definite symmetric matrix B, 
!  i.e. B*OP = (OP')*B, then subroutine ssaupd should be used instead.
!
!  The computed approximate eigenvalues are called Ritz values and
!  the corresponding approximate eigenvectors are called Ritz vectors.
!
!  dnaupd is usually called iteratively to solve one of the 
!  following problems:
!
!  Mode 1:  A*x = lambda*x.
!           ===> OP = A  and  B = I.
!
!  Mode 2:  A*x = lambda*M*x, M symmetric positive definite
!           ===> OP = inv[M]*A  and  B = M.
!           ===> (If M can be factored see remark 3 below)
!
!  Mode 3:  A*x = lambda*M*x, M symmetric semi-definite
!           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M. 
!           ===> shift-and-invert mode (in real arithmetic)
!           If OP*x = amu*x, then 
!           amu = 1/2 * [ 1/(lambda-sigma) + 1/(lambda-conjg(sigma)) ].
!           Note: If sigma is real, i.e. imaginary part of sigma is zero;
!                 Real_Part{ inv[A - sigma*M]*M } == inv[A - sigma*M]*M 
!                 amu == 1/(lambda-sigma). 
!  
!  Mode 4:  A*x = lambda*M*x, M symmetric semi-definite
!           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M. 
!           ===> shift-and-invert mode (in real arithmetic)
!           If OP*x = amu*x, then 
!           amu = 1/2i * [ 1/(lambda-sigma) - 1/(lambda-conjg(sigma)) ].
!
!  Both mode 3 and 4 give the same enhancement to eigenvalues close to
!  the (complex) shift sigma.  However, as lambda goes to infinity,
!  the operator OP in mode 4 dampens the eigenvalues more strongly than
!  does OP defined in mode 3.
!
!  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
!        should be accomplished either by a direct method
!        using a sparse matrix factorization and solving
!
!           [A - sigma*M]*w = v  or M*w = v,
!
!        or through an iterative method for solving these
!        systems.  If an iterative method is used, the
!        convergence test must be more stringent than
!        the accuracy requirements for the eigenvalue
!        approximations.
!
!\Usage:
!  call dnaupd
!     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
!       IPNTR, WORKD, WORKL, LWORKL, INFO )
!
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.  IDO must be zero on the first 
!          call to dnaupd.  IDO will be set internally to
!          indicate the type of operation to be performed.  Control is
!          then given back to the calling routine which has the
!          responsibility to carry out the requested operation and call
!          dnaupd with the result.  The operand is given in
!          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    This is for the initialization phase to force the
!                    starting vector into the range of OP.
!          IDO =  1: compute  Y = OP * Z  and Z = B * X where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y,
!                    IPNTR(3) is the pointer into WORKD for Z.
!          IDO =  2: compute  Y = M * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!          IDO =  3: compute and return the shifts in the first 
!                    NP locations of WORKL.
!          IDO =  4: compute Z = OP * X
!          IDO = 99: done
!
!          After the initialization phase, when the routine is used in 
!          the "shift-and-invert" mode, the vector B * X is already 
!          available and does not need to be recomputed in forming OP*X.
!             
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B that defines the
!          semi-inner product for the operator OP.
!          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
!          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
!
!  N       Integer.  (INPUT)
!          Dimension of the eigenproblem.
!
!  WHICH   Character*2.  (INPUT)
!          'LM' -> want the NEV eigenvalues of largest magnitude.
!          'SM' -> want the NEV eigenvalues of smallest magnitude.
!          'LR' -> want the NEV eigenvalues of largest real part.
!          'SR' -> want the NEV eigenvalues of smallest real part.
!          'LI' -> want the NEV eigenvalues of largest imaginary part.
!          'SI' -> want the NEV eigenvalues of smallest imaginary part.
!
!  NEV     Integer.  (INPUT)
!          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
!
!  TOL     Double precision scalar.  (INPUT)
!          Stopping criterion: the relative accuracy of the Ritz value 
!          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
!          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
!          DEFAULT = DLAMCH('EPS')  (machine precision as computed
!                    by the LAPACK auxiliary subroutine DLAMCH).
!
!  RESID   Double precision array of length N.  (INPUT/OUTPUT)
!          On INPUT: 
!          If INFO .EQ. 0, a random initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          On OUTPUT:
!          RESID contains the final residual vector.
!
!  NCV     Integer.  (INPUT)
!          Number of columns of the matrix V. NCV must satisfy the two
!          inequalities 2 <= NCV-NEV and NCV <= N.
!          This will indicate how many Arnoldi vectors are generated 
!          at each iteration.  After the startup phase in which NEV 
!          Arnoldi vectors are generated, the algorithm generates 
!          approximately NCV-NEV Arnoldi vectors at each subsequent update 
!          iteration. Most of the cost in generating each Arnoldi vector is 
!          in the matrix-vector operation OP*x. 
!          NOTE: 2 <= NCV-NEV in order that complex conjugate pairs of Ritz 
!          values are kept together. (See remark 4 below)
!
!  V       Double precision array N by NCV.  (OUTPUT)
!          Contains the final set of Arnoldi basis vectors. 
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling program.
!
!  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
!          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
!          The shifts selected at each iteration are used to restart
!          the Arnoldi iteration in an implicit fashion.
!          -------------------------------------------------------------
!          ISHIFT = 0: the shifts are provided by the user via
!                      reverse communication.  The real and imaginary
!                      parts of the NCV eigenvalues of the Hessenberg
!                      matrix H are returned in the part of the WORKL 
!                      array corresponding to RITZR and RITZI. See remark 
!                      5 below.
!          ISHIFT = 1: exact shifts with respect to the current
!                      Hessenberg matrix H.  This is equivalent to 
!                      restarting the iteration with a starting vector
!                      that is a linear combination of approximate Schur
!                      vectors associated with the "wanted" Ritz values.
!          -------------------------------------------------------------
!
!          IPARAM(2) = No longer referenced.
!
!          IPARAM(3) = MXITER
!          On INPUT:  maximum number of Arnoldi update iterations allowed. 
!          On OUTPUT: actual number of Arnoldi update iterations taken. 
!
!          IPARAM(4) = NB: blocksize to be used in the recurrence.
!          The code currently works only for NB = 1.
!
!          IPARAM(5) = NCONV: number of "converged" Ritz values.
!          This represents the number of Ritz values that satisfy
!          the convergence criterion.
!
!          IPARAM(6) = IUPD
!          No longer referenced. Implicit restarting is ALWAYS used.  
!
!          IPARAM(7) = MODE
!          On INPUT determines what type of eigenproblem is being solved.
!          Must be 1,2,3,4; See under \Description of dnaupd for the 
!          four modes available.
!
!          IPARAM(8) = NP
!          When ido = 3 and the user provides shifts through reverse
!          communication (IPARAM(1)=0), dnaupd returns NP, the number
!          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
!          5 below.
!
!          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
!          OUTPUT: NUMOP  = total number of OP*x operations,
!                  NUMOPB = total number of B*x operations if BMAT='G',
!                  NUMREO = total number of steps of re-orthogonalization.        
!
!  IPNTR   Integer array of length 14.  (OUTPUT)
!          Pointer to mark the starting locations in the WORKD and WORKL
!          arrays for matrices/vectors used by the Arnoldi iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X in WORKD.
!          IPNTR(2): pointer to the current result vector Y in WORKD.
!          IPNTR(3): pointer to the vector B * X in WORKD when used in 
!                    the shift-and-invert mode.
!          IPNTR(4): pointer to the next available location in WORKL
!                    that is untouched by the program.
!          IPNTR(5): pointer to the NCV by NCV upper Hessenberg matrix
!                    H in WORKL.
!          IPNTR(6): pointer to the real part of the ritz value array 
!                    RITZR in WORKL.
!          IPNTR(7): pointer to the imaginary part of the ritz value array
!                    RITZI in WORKL.
!          IPNTR(8): pointer to the Ritz estimates in array WORKL associated
!                    with the Ritz values located in RITZR and RITZI in WORKL.
!
!          Note: IPNTR(9:13) is only referenced by dneupd. See Remark 2 below.
!
!          IPNTR(9): pointer to the real part of the NCV RITZ values of the 
!                    original system.
!          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of 
!                     the original system.
!          IPNTR(11): pointer to the NCV corresponding error bounds.
!          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
!                     Schur matrix for H.
!          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
!                     of the upper Hessenberg matrix H. Only referenced by
!                     dneupd if RVEC = .TRUE. See Remark 2 below.
!          Note: IPNTR(9:13) is only referenced by dneupd. See Remark 2 below.
!          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
!          -------------------------------------------------------------
!          
!  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The user should not use WORKD 
!          as temporary workspace during the iteration. Upon termination
!          WORKD(1:N) contains B*RESID(1:N). If an invariant subspace
!          associated with the converged Ritz values is desired, see remark
!          2 below, subroutine dneupd uses this output.
!          See Data Distribution Note below.  
!
!  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  See Data Distribution Note below.
!
!  LWORKL  Integer.  (INPUT)
!          LWORKL must be at least 3*NCV**2 + 6*NCV.
!
!  INFO    Integer.  (INPUT/OUTPUT)
!          If INFO .EQ. 0, a randomly initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          Error flag on output.
!          =  0: Normal exit.
!          =  1: Maximum number of iterations taken.
!                All possible eigenvalues of OP has been found. IPARAM(5)  
!                returns the number of wanted converged Ritz values.
!          =  2: No longer an informational error. Deprecated starting
!                with release 2 of ARPACK.
!          =  3: No shifts could be applied during a cycle of the 
!                Implicitly restarted Arnoldi iteration. One possibility 
!                is to increase the size of NCV relative to NEV. 
!                See remark 4 below.
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV-NEV >= 2 and less than or equal to N.
!          = -4: The maximum number of Arnoldi update iteration 
!                must be greater than zero.
!          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work array is not sufficient.
!          = -8: Error return from LAPACK eigenvalue calculation;
!          = -9: Starting vector is zero.
!          = -10: IPARAM(7) must be 1,2,3,4.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
!          = -12: IPARAM(1) must be equal to 0 or 1.
!          = -9999: Could not build an Arnoldi factorization.
!                   IPARAM(5) returns the size of the current Arnoldi
!                   factorization.
!
!\Remarks
!  1. The computed Ritz values are approximate eigenvalues of OP. The
!     selection of WHICH should be made with this in mind when
!     Mode = 3 and 4.  After convergence, approximate eigenvalues of the
!     original problem may be obtained with the ARPACK subroutine dneupd.
!
!  2. If a basis for the invariant subspace corresponding to the converged Ritz 
!     values is needed, the user must call dneupd immediately following 
!     completion of dnaupd. This is new starting with release 2 of ARPACK.
!
!  3. If M can be factored into a Cholesky factorization M = LL'
!     then Mode = 2 should not be selected.  Instead one should use
!     Mode = 1 with  OP = inv(L)*A*inv(L').  Appropriate triangular 
!     linear systems should be solved with L and L' rather
!     than computing inverses.  After convergence, an approximate
!     eigenvector z of the original problem is recovered by solving
!     L'z = x  where x is a Ritz vector of OP.
!
!  4. At present there is no a-priori analysis to guide the selection of NCV 
!     relative to NEV.  The only formal requirement is that NCV > NEV + 2.
!     However, it is recommended that NCV .ge. 2*NEV+1.  If many problems of
!     the same type are to be solved, one should experiment with increasing
!     NCV while keeping NEV fixed for a given test problem.  This will 
!     usually decrease the required number of OP*x operations but it
!     also increases the work and storage required to maintain the orthogonal
!     basis vectors.  The optimal "cross-over" with respect to CPU time
!     is problem dependent and must be determined empirically. 
!     See Chapter 8 of Reference 2 for further information.
!
!  5. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the 
!     NP = IPARAM(8) real and imaginary parts of the shifts in locations 
!         real part                  imaginary part
!         -----------------------    --------------
!     1   WORKL(IPNTR(14))           WORKL(IPNTR(14)+NP)
!     2   WORKL(IPNTR(14)+1)         WORKL(IPNTR(14)+NP+1)
!                        .                          .
!                        .                          .
!                        .                          .
!     NP  WORKL(IPNTR(14)+NP-1)      WORKL(IPNTR(14)+2*NP-1).
!
!     Only complex conjugate pairs of shifts may be applied and the pairs 
!     must be placed in consecutive locations. The real part of the 
!     eigenvalues of the current upper Hessenberg matrix are located in 
!     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1) and the imaginary part 
!     in WORKL(IPNTR(7)) through WORKL(IPNTR(7)+NCV-1). They are ordered
!     according to the order defined by WHICH. The complex conjugate
!     pairs are kept together and the associated Ritz estimates are located in
!     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
!
!-----------------------------------------------------------------------




