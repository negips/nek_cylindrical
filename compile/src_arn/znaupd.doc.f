!> @file znaupd.doc.f
!! @ingroup arn_arp
!! @brief Description of znaupd subroutine in arpack 
!! @author Prabal Negi
!! @date June 17, 2021
!
!=======================================================================
!\Name: znaupd
!
!\Description: 
!  Reverse communication interface for the Implicitly Restarted Arnoldi
!  iteration. This is intended to be used to find a few eigenpairs of a 
!  complex linear operator OP with respect to a semi-inner product defined 
!  by a hermitian positive semi-definite real matrix B. B may be the identity 
!  matrix.  NOTE: if both OP and B are real, then dsaupd or dnaupd should
!  be used.
!
!
!  The computed approximate eigenvalues are called Ritz values and
!  the corresponding approximate eigenvectors are called Ritz vectors.
!
!  znaupd is usually called iteratively to solve one of the 
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
!           ===> OP =  inv[A - sigma*M]*M   and  B = M. 
!           ===> shift-and-invert mode 
!           If OP*x = amu*x, then lambda = sigma + 1/amu.
!  
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
!  call znaupd
!     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
!       IPNTR, WORKD, WORKL, LWORKL, RWORK, INFO )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.  IDO must be zero on the first 
!          call to znaupd.  IDO will be set internally to
!          indicate the type of operation to be performed.  Control is
!          then given back to the calling routine which has the
!          responsibility to carry out the requested operation and call
!          znaupd with the result.  The operand is given in
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
!          -------------------------------------------------------------
!          After the initialization phase, when the routine is used in 
!          the "shift-and-invert" mode, the vector M * X is already 
!          available and does not need to be recomputed in forming OP*X.
!             
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B that defines the
!          semi-inner product for the operator OP.
!          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
!          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x
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
!  TOL     Double precision  scalar.  (INPUT)
!          Stopping criteria: the relative accuracy of the Ritz value 
!          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
!          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
!          DEFAULT = dlamch('EPS')  (machine precision as computed
!                    by the LAPACK auxiliary subroutine dlamch).
!
!  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
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
!  V       Complex*16 array N by NCV.  (OUTPUT)
!          Contains the final set of Arnoldi basis vectors. 
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling program.
!
!  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
!          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
!          The shifts selected at each iteration are used to filter out
!          the components of the unwanted eigenvector.
!          -------------------------------------------------------------
!          ISHIFT = 0: the shifts are to be provided by the user via
!                      reverse communication.  The NCV eigenvalues of 
!                      the Hessenberg matrix H are returned in the part
!                      of WORKL array corresponding to RITZ.
!          ISHIFT = 1: exact shifts with respect to the current
!                      Hessenberg matrix H.  This is equivalent to 
!                      restarting the iteration from the beginning 
!                      after updating the starting vector with a linear
!                      combination of Ritz vectors associated with the 
!                      "wanted" eigenvalues.
!          ISHIFT = 2: other choice of internal shift to be defined.
!          -------------------------------------------------------------
!
!          IPARAM(2) = No longer referenced 
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
!          Must be 1,2,3,4; See under \Description of znaupd for the 
!          four modes available.
!
!          IPARAM(8) = NP
!          When ido = 3 and the user provides shifts through reverse
!          communication (IPARAM(1)=0), _naupd returns NP, the number
!          of shifts the user is to provide. 0 < NP < NCV-NEV.
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
!          IPNTR(5): pointer to the NCV by NCV upper Hessenberg
!                    matrix H in WORKL.
!          IPNTR(6): pointer to the  ritz value array  RITZ
!          IPNTR(7): pointer to the (projected) ritz vector array Q
!          IPNTR(8): pointer to the error BOUNDS array in WORKL.
!          Note: IPNTR(9:13) is only referenced by zneupd. See Remark 2 below.
!          IPNTR(9): pointer to the NCV RITZ values of the 
!                    original system.
!          IPNTR(10): Not Used
!          IPNTR(11): pointer to the NCV corresponding error bounds.
!          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
!          -------------------------------------------------------------
!          
!  WORKD   Complex*16 work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The user should not use WORKD 
!          as temporary workspace during the iteration !!!!!!!!!!
!          See Data Distribution Note below.  
!
!  WORKL   Complex*16 work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  See Data Distribution Note below.
!
!  LWORKL  Integer.  (INPUT)
!          LWORKL must be at least 3*NCV**2 + 5*NCV.
!
!  RWORK   Double precision  work array of length NCV (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.
!
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
!          = -10: IPARAM(7) must be 1,2,3.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
!          = -12: IPARAM(1) must be equal to 0 or 1.
!          = -9999: Could not build an Arnoldi factorization.
!                   User input error highly likely.  Please
!                   check actual array dimensions and layout.
!                   IPARAM(5) returns the size of the current Arnoldi
!                   factorization.
!
!\Remarks
!  1. The computed Ritz values are approximate eigenvalues of OP. The
!     selection of WHICH should be made with this in mind when using
!     Mode = 3.  When operating in Mode = 3 setting WHICH = 'LM' will
!     compute the NEV eigenvalues of the original problem that are
!     closest to the shift SIGMA . After convergence, approximate eigenvalues 
!     of the original problem may be obtained with the ARPACK subroutine zneupd.
!
!  2. If a basis for the invariant subspace corresponding to the converged Ritz 
!     values is needed, the user must call zneupd immediately following 
!     completion of znaupd. This is new starting with release 2 of ARPACK.
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
!     NP = IPARAM(8) complex shifts in locations
!     WORKL(IPNTR(14)), WORKL(IPNTR(14)+1), ... , WORKL(IPNTR(14)+NP).
!     Eigenvalues of the current upper Hessenberg matrix are located in
!     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1). They are ordered
!     according to the order defined by WHICH.  The associated Ritz estimates
!     are located in WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... ,
!     WORKL(IPNTR(8)+NCV-1).
!
!-----------------------------------------------------------------------



