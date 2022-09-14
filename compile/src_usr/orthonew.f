!====================================================================== 
!
!     Author: Prabal negi
!     Description: Slightly different way of removing mean pressure.
!
!
!======================================================================       

      subroutine ortho_new(respr)

C     Orthogonalize the residual in the pressure solver with respect 
C     to (1,1,...,1)T  (only if all Dirichlet b.c.).

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'TSTEP'
      include 'MASS'          ! bm2

      real respr (lx2,ly2,lz2,lelv)
      integer*8 ntotg,nxyz2

      integer ntot
      real rlam
      real glsum

      real glsc3,glsc2

      nxyz2 = lx2*ly2*lz2
      ntot  = nxyz2*nelv
      ntotg = nxyz2*nelgv

      if (ifield.eq.1) then
         if (ifvcor) then
            rlam  = glsum (respr,ntot)/ntotg
            call cadd (respr,-rlam,ntot)

!            rlam  = glsum (respr,ntot)/volvm2
!            rlam  = glsum (respr,ntot)/ntotg
!            call cadd (respr,-rlam,ntot)
!            call add2s2(respr,bm2,-rlam,ntot)

!            rlam = glsc2(respr,bm2,ntot)/volvm2
!            call cadd(respr,-rlam,ntot)

         endif
       elseif (ifield.eq.ifldmhd) then
         if (ifbcor) then
            rlam = glsum (respr,ntot)/ntotg
            call cadd (respr,-rlam,ntot)
         endif
       else
         call exitti('ortho: unaccounted ifield = $',ifield)
      endif

      return
      end
c------------------------------------------------------------------------

      subroutine ortho_left(respr)

!     portho = (I - (B^T)*(p*p^T))*respr


C     Orthogonalize the residual in the pressure solver with respect 
C     to (1,1,...,1)T  (only if all Dirichlet b.c.).

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'TSTEP'
      include 'MASS'          ! bm2

      real respr (lx2,ly2,lz2,lelv)
      integer*8 ntotg,nxyz2

      integer ntot
      real rlam
      real glsum

      real glsc3,glsc2

      nxyz2 = lx2*ly2*lz2
      ntot  = nxyz2*nelv
      ntotg = nxyz2*nelgv

      if (ifield.eq.1) then
         if (ifvcor) then
!            rlam  = glsum (respr,ntot)/ntotg
!            call cadd (respr,-rlam,ntot)

            rlam  = glsum (respr,ntot)/volvm2
!            rlam  = glsum (respr,ntot)/ntotg
!            call cadd (respr,-rlam,ntot)
            call add2s2(respr,bm2,-rlam,ntot)

!            rlam = glsc2(respr,bm2,ntot)/volvm2
!            call cadd(respr,-rlam,ntot)

         endif
       elseif (ifield.eq.ifldmhd) then
         if (ifbcor) then
            rlam = glsum (respr,ntot)/ntotg
            call cadd (respr,-rlam,ntot)
         endif
       else
         call exitti('ortho: unaccounted ifield = $',ifield)
      endif

      return
      end
c------------------------------------------------------------------------
      subroutine ortho_right(respr)


!     portho = (I - (p*p^T)B)*respr

C     Orthogonalize the residual in the pressure solver with respect 
C     to (1,1,...,1)T  (only if all Dirichlet b.c.).

      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'TSTEP'
      include 'MASS'          ! bm2

      real respr (lx2,ly2,lz2,lelv)
      integer*8 ntotg,nxyz2

      integer ntot
      real rlam
      real glsum

      real glsc3,glsc2

      nxyz2 = lx2*ly2*lz2
      ntot  = nxyz2*nelv
      ntotg = nxyz2*nelgv

      if (ifield.eq.1) then
         if (ifvcor) then
!            rlam  = glsum (respr,ntot)/ntotg
!            call cadd (respr,-rlam,ntot)

!            rlam  = glsum (respr,ntot)/volvm2
!            rlam  = glsum (respr,ntot)/ntotg
!            call cadd (respr,-rlam,ntot)
!            call add2s2(respr,bm2,-rlam,ntot)

            rlam = glsc2(respr,bm2,ntot)/volvm2
            call cadd(respr,-rlam,ntot)

         endif
       elseif (ifield.eq.ifldmhd) then
         if (ifbcor) then
            rlam = glsum (respr,ntot)/ntotg
            call cadd (respr,-rlam,ntot)
         endif
       else
         call exitti('ortho: unaccounted ifield = $',ifield)
      endif

      return
      end
c------------------------------------------------------------------------

