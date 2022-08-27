!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!     Description: Complex arithmetic for Arnoldi
!     Author:      Prabal Negi
!
!----------------------------------------------------------------------

      subroutine czero(a,n)

      implicit none

      integer i,n
      complex*16 a(1)

      do i=1,n
        a(i) = cmplx(0.,0.)
      enddo

      return
      end subroutine czero
!---------------------------------------------------------------------- 
      subroutine ccopy(a,b,n)

      implicit none

      integer i,n
      complex*16 a(1),b(1)

      do i=1,n
        a(i) = b(i)
      enddo

      return
      end subroutine ccopy
!---------------------------------------------------------------------- 

      subroutine copytocomplex(c,a,b,n)

      implicit none

      integer i,n
      complex*16 c(1)
      real a(1),b(1)

      do i=1,n
        c(i) = cmplx(a(i),b(i))
      enddo

      return
      end subroutine copytocomplex
!----------------------------------------------------------------------
      subroutine copytoreal(a,b,c,n)

      implicit none

      integer i,n
      complex*16 c(1)
      real a(1),b(1)

      do i=1,n
        a(i) = real(c(i))
        b(i) = aimag(c(i))
      enddo

      return
      end subroutine copytoreal
!----------------------------------------------------------------------

      subroutine col2_cr(a,b,n)

      implicit none

      integer i,n
      complex*16 a(1)
      real b(1)

      do i=1,n
        a(i) = a(i)*b(i)
      enddo

      return
      end subroutine col2_cr
!---------------------------------------------------------------------- 
      subroutine cmult_cr(a,c,n)

      implicit none

      integer i,n
      complex*16 a(1)
      real c

      do i=1,n
        a(i) = a(i)*c
      enddo

      return
      end subroutine cmult_cr
!---------------------------------------------------------------------- 



