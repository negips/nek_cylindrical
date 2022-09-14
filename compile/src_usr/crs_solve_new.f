!====================================================================== 
!     Author: Prabal Negi
!     Description: Coarse solver using a modified M2 -> M1 mapping.
!
!
!======================================================================       
!---------------------------------------------------------------------- 
      subroutine crs_solve_l2_test(uf,vf)
c
c     Given an input vector v, this generates the H1 coarse-grid solution
c
      include 'SIZE'
      include 'DOMAIN'
      include 'ESOLV'
      include 'GEOM'
      include 'PARALLEL'
      include 'SOLN'
      include 'INPUT'
      include 'TSTEP'
      real uf(1),vf(1)
      common /scrpre/ uc(lcr*lelt),w(2*lx1*ly1*lz1)

      call map_f_to_c_l2_bilin_test(uf,vf,w)
      call fgslib_crs_solve(xxth(ifield),uc,uf)
      call map_c_to_f_l2_bilin(uf,uc,w)

      return
      end
c
c-----------------------------------------------------------------------

      subroutine map_f_to_c_l2_bilin_test(uc,uf,w)

c     TRANSPOSE of L2 Iterpolation operator:                    T
c                                 (linear --> spectral GLL mesh)

      include 'SIZE'
      include 'DOMAIN'
      include 'INPUT'

      parameter (lxyz = lx2*ly2*lz2)
      real uc(nxyz_c,lelt),uf(lxyz,lelt),w(1)

      ltot22 = 2*lx2*ly2*lz2
      nx_crs = 2   ! bilinear only

      do ie=1,nelv
         call maph1_to_l2t_test(uc(1,ie),nx_crs,uf(1,ie),
     $                          lx2,if3d,w,ltot22)
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------

      subroutine maph1_to_l2t_test(b,nb,a,na,if3d,w,ldw)
c
c     Input:   a
c     Output:  b
c
      real a(1),b(1),w(1)
      logical if3d
c
      parameter(lx=50)
      real za(lx),zb(lx)
c
      real iba(lx*lx),ibat(lx*lx)
      save iba,ibat
c
      integer nao,nbo
      save    nao,nbo
      data    nao,nbo  / -9, -9/
c
c
      if (na.gt.lx.or.nb.gt.lx) then
         write(6,*)'ERROR: increase lx in maph1_to_l2 to max:',na,nb
         call exitt
      endif
c
      if (na.ne.nao  .or.   nb.ne.nbo) then
         nao = na
         nbo = nb
         call zwgl (za,w,na)
         call zwgll(zb,w,nb)
!         call igllm(iba,ibat,zb,za,nb,na,nb,na)
         call iglm(iba,ibat,za,zb,na,nb,na,nb)
      endif
c
!      call specmpn(b,nb,a,na,ibat,iba,if3d,w,ldw)
      call specmpn(b,nb,a,na,iba,ibat,if3d,w,ldw)
     
c
      return
      end
c
c-----------------------------------------------------------------------

