!====================================================================== 
!     Author: Prabal Negi
!     Description: Routines for Continuous pressure 
!                  approximation space
!
!======================================================================       


      SUBROUTINE DGLLGLL (D,DT,ZM1,ZM2,IM12,NZM1,NZM2,ND1,ND2)
C-----------------------------------------------------------------------
C
C     Compute the (one-dimensional) derivative matrix D and its
C     transpose DT associated with taking the derivative of a variable
C     expanded on a Gauss-Lobatto Legendre mesh (M1), and evaluate its
C     derivative on a Guass Legendre mesh (M2).
C     Need the one-dimensional interpolation operator IM12
C     (see subroutine IGLLGL).
C     Note: D and DT are rectangular matrices.
C
C-----------------------------------------------------------------------
      REAL D(ND2,ND1), DT(ND1,ND2), ZM1(ND1), ZM2(ND2), IM12(ND2,ND1)
      IF (NZM1.EQ.1) THEN
        D (1,1) = 0.
        DT(1,1) = 0.
        RETURN
      ENDIF
      EPS = 1.E-6
      NM1 = NZM1-1

      FN = NM1
      d0 = FN*(FN+1.)/4.

      DO 10 IP = 1, NZM2
         DO 10 JQ = 1, NZM1
            ZP = ZM2(IP)
            ZQ = ZM1(JQ)
            IF ((ABS(ZP) .LT. EPS).AND.(ABS(ZQ) .LT. EPS)) THEN
                D(IP,JQ) = 0.
            ELSEIF ((IP.EQ.1).AND.(JQ.EQ.1)) THEN
                D(IP,JQ) = -d0
            ELSEIF ((IP.EQ.NZM2).AND.(JQ.EQ.NZM1)) THEN
                D(IP,JQ) =  d0
            ELSE
                D(IP,JQ) = (PNLEG(ZP,NM1)/PNLEG(ZQ,NM1)
     $                     -IM12(IP,JQ))/(ZP-ZQ)
            ENDIF
            DT(JQ,IP) = D(IP,JQ)
 10   CONTINUE
      RETURN
      END SUBROUTINE DGLLGLL
!---------------------------------------------------------------------- 
