!====================================================================== 
!     Author: Prabal Negi
!     Description: Resetting the h1_crs preconditioner in Nek.
!
!====================================================================== 
      subroutine reset_preconditioner()

      implicit none

      include 'SIZE'
      include 'INPUT'   ! uparam
      include 'TSTEP'   ! ifield
      include 'SOLN'    ! vtrans

      integer n
      real h1(lx1,ly1,lz1,lelv)
      real h2(lx1,ly1,lz1,lelv)
      integer ist

      logical iffdmr    ! if fdm reset
      logical ifcrsr    ! if coarse grid reset

      real*8 t0
      real*8 dnekclock

      t0 = dnekclock()

      n = lx1*ly1*lz1*nelv

      iffdmr = .false.
      ifcrsr = .false.

      if (uparam(5).eq.1.0) iffdmr = .true.
      if (uparam(6).eq.1.0) ifcrsr = .true.

!     The multigrid has not been modified yet
      param(43) = 1
      ifmgrid   = .false.

      ifield = 1
      if (iffdmr) call set_overlap_again

      if (ifcrsr) then
        call invers2(h1,vtrans(1,1,1,1,ifield),n)
        call rzero(h2,n)
        call set_up_h1_crs_again(h1,h2)     ! reversed
      endif  
      t0 = dnekclock() - t0

      if (iffdmr.or.ifcrsr) then
        if (nio.eq.0) then
           write(6,*) 'Preconditioner Resets ',t0, ' sec'
        endif
      endif  


      return
      end subroutine

!----------------------------------------------------------------------       

      subroutine set_up_h1_crs_again(h1_usr,h2_usr)

      include 'SIZE'
      include 'GEOM'
      include 'DOMAIN'
      include 'INPUT'
      include 'PARALLEL'
      include 'TSTEP'
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      common /ivrtx/ vertex ((2**ldim)*lelt)
      integer*8 vertex

      integer gs_handle
      integer null_space,e

      character*3 cb
      common /scrxxt/ cmlt(lcr,lelv),mask(lcr,lelv)
      common /scrxxti/ ia(lcr,lcr,lelv), ja(lcr,lcr,lelv)
      real mask
      integer ia,ja
      real z

      common /scrch/ iwork(2,lx1*ly1*lz1*lelv)
      common /scrns/ w(7*lx1*ly1*lz1*lelv)
!      common /vptsol/ a(27*lx1*ly1*lz1*lelv)   ! Common block used for
                                                ! solutions
      real a,acopy                                                
      common /h1crsa/ a(lcr*lcr*lelv)      ! Could use some common
     $              , acopy(lcr*lcr*lelv)  ! block for this

      integer w
      real wr(1)
      equivalence (wr,w)

      common /scrvhx/ h1(lx1*ly1*lz1*lelv),h2(lx1*ly1*lz1*lelv)
      common /scrmgx/ w1(lx1*ly1*lz1*lelv),w2(lx1*ly1*lz1*lelv)

      integer*8 ngv
      character*132 amgfile_c
      character*1   fname1(132)
      equivalence  (fname1,amgfile_c)
      integer nnamg

      integer i,j

      real h1_usr(lx1*ly1*lz1*lelv)
      real h2_usr(lx1*ly1*lz1*lelv)

      t0 = dnekclock()

c     nxc is order of coarse grid space + 1, nxc=2, linear, 3=quad,etc.
c     nxc=param(82)
c     if (nxc.gt.lxc) then
c        nxc=lxc
c        write(6,*) 'WARNING :: coarse grid space too large',nxc,lxc 
c     endif
c     if (nxc.lt.2) nxc=2

      nxc     = 2
      nx_crs  = nxc

!      if(nio.eq.0) write(6,*) 'setup h1 coarse grid (again), nx_crs=', 
!     $                         nx_crs

      ncr     = nxc**ldim
      nxyz_c  = ncr
c
c     Set SEM_to_GLOB
c
      call get_vertex
      call set_vert(se_to_gcrs,ngv,nxc,nelv,vertex,.true.)

c     Set mask
      z=0
      ntot=nelv*nxyz_c
      nzc=1
      if (if3d) nzc=nxc
      call rone(mask,ntot)
      call rone(cmlt,ntot)
      nfaces=2*ldim
c     ifield=1			!c? avo: set in set_overlap through 'TSTEP'?

      if (ifield.eq.1) then
         do ie=1,nelv
         do iface=1,nfaces
            cb=cbc(iface,ie,ifield)
            if (cb.eq.'o  '  .or.  cb.eq.'on '  .or. 
     $          cb.eq.'O  '  .or.  cb.eq.'ON '  .or.  cb.eq.'MM '  .or.
     $          cb.eq.'mm '  .or.  cb.eq.'ms '  .or.  cb.eq.'MS ')
     $           call facev(mask,ie,iface,z,nxc,nxc,nzc) ! 'S* ' & 's* ' ?avo?
         enddo
         enddo
      elseif (ifield.eq.ifldmhd) then   ! no ifmhd ?avo?
         do ie=1,nelv
         do iface=1,nfaces
            cb=cbc(iface,ie,ifield)
            if (cb.eq.'ndd'  .or.  cb.eq.'dnd'  .or.  cb.eq.'ddn')
     $          call facev(mask,ie,iface,z,nxc,nxc,nzc)
         enddo
         enddo
      endif

c     Set global index of dirichlet nodes to zero; xxt will ignore them


      call fgslib_gs_setup(gs_handle,se_to_gcrs,ntot,nekcomm,mp)
      call fgslib_gs_op   (gs_handle,mask,1,2,0)  !  "*"
      call fgslib_gs_op   (gs_handle,cmlt,1,1,0)  !  "+"
      call fgslib_gs_free (gs_handle)             ! Free up memory
      call set_jl_crs_mask(ntot,mask,se_to_gcrs)

      call invcol1(cmlt,ntot)

c     Setup local SEM-based Neumann operators (for now, just full...)

c      if (param(51).eq.1) then     ! old coarse grid
c         nxyz1=lx1*ly1*lz1
c         lda = 27*nxyz1*lelt
c         ldw =  7*nxyz1*lelt
c         call get_local_crs(a,lda,nxc,h1,h2,w,ldw)
c      else
c        NOTE: a(),h1,...,w2() must all be large enough
         n = lx1*ly1*lz1*nelv
!         call rone  (h1,n)
!         call rzero (h2,n)
         call copy  (h1,h1_usr,n)
         call copy  (h2,h2_usr,n)
         call get_local_crs_galerkin(a,ncr,nxc,h1,h2,w1,w2)
!         call get_wt_local_crs_galerkin(a,ncr,nxc,h1,h2,w1,w2)

         call copy(acopy,a,ncr*ncr*nelv)  ! for debugging
c      endif

      call set_mat_ij(ia,ja,ncr,nelv)
      null_space=0
      if (ifield.eq.1) then
         if (ifvcor)  null_space=1
      elseif (ifield.eq.ifldmhd) then
         if (ifbcor)  null_space=1
      endif

      nz=ncr*ncr*nelv
      isolver = param(40)

      call blank(fname1,132)
      lamgn = ltrunc(amgfile,len(amgfile))
      call chcopy(fname1,amgfile,lamgn)
      call chcopy(fname1(lamgn+1),char(0),1)

      ierr = 0
      call fgslib_crs_free(xxth(ifield))  ! free up the memory
      call fgslib_crs_setup(xxth(ifield),isolver,nekcomm,mp,ntot,
     $     se_to_gcrs,nz,ia,ja,a, null_space, crs_param, 
     $     amgfile_c,ierr)
      ierr = iglmax(ierr,1)
      if (ifneknek) ierr = iglmax_ms(ierr,1)
      if (ierr.eq.1) then
         call exitt
      endif

      t0 = dnekclock()-t0
      if (nio.eq.0) then
         write(6,*) 'h1 coarse grid (reset) ',t0, ' sec'
      endif

      return
      end
c
c-----------------------------------------------------------------------

      subroutine get_wt_local_crs_galerkin(a,ncl,nxc,h1,h2,w1,w2)

c     This routine generates Nelv submatrices of order ncl using
c     Weighted Galerkin projection

      include 'SIZE'
      include 'MASS'          ! BM1 (weight)

      real    a(ncl,ncl,1),h1(1),h2(1)
      real    w1(lx1*ly1*lz1,nelv),w2(lx1*ly1*lz1,nelv)

      parameter (lcrd=lx1**ldim)
      common /ctmp1z/ b(lcrd,8)

      integer e

      do j=1,ncl
         call gen_crs_basis(b(1,j),j) ! bi- or tri-linear interpolant
      enddo

      isd  = 1
      imsh = 1

      nxyz = lx1*ly1*lz1
      do j = 1,ncl
         do e = 1,nelv
            call copy(w1(1,e),b(1,j),nxyz)
         enddo

         call axhelm (w2,w1,h1,h2,imsh,isd)        ! A^e * bj
         call col2(w2,bm1,nxyz*nelv)               ! Multiply by weight

         do e = 1,nelv
         do i = 1,ncl
            a(i,j,e) = vlsc2(b(1,i),w2(1,e),nxyz)  ! bi^T * B * A^e * bj
         enddo
         enddo

      enddo

      return
      end
c-----------------------------------------------------------------------

      subroutine crs_solve_l2_wt(uf,vf)
c
c     Given an input vector v, this generates the H1 coarse-grid solution
c     Using Weighted projection
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
      common /scrpre/ uc(lcr*lelt),w1(2*lx1*ly1*lz1),w2(2*lx1*ly1*lz1)

      call map_f_to_c_l2_bilin_wt(uf,vf,w1,w2)
      call fgslib_crs_solve(xxth(ifield),uc,uf)
      call map_c_to_f_l2_bilin(uf,uc,w1)

      return
      end
c
c-----------------------------------------------------------------------
      subroutine map_f_to_c_l2_bilin_wt(uc,uf,w1,w2)

c     TRANSPOSE of L2 Iterpolation operator:                    T
c                                 (linear --> spectral GLL mesh)

      implicit none

      include 'SIZE'
      include 'DOMAIN'
      include 'MASS'          ! BM1
      include 'INPUT'

      integer lxyz
      parameter (lxyz = lx2*ly2*lz2)
      real uc(nxyz_c,lelt),uf(lxyz,lelt)
      real w1(1),w2(1)

      integer ltot22
      integer ie


      ltot22 = 2*lx2*ly2*lz2
      nx_crs = 2   ! bilinear only

      do ie=1,nelv
!       Using BM2 as the weight. Integration is done on the M2 mesh.      
        call maph1_to_l2t_wt(uc(1,ie),nx_crs,uf(1,ie),
     $                         lx2,if3d,w1,w2,ltot22,bm2(1,1,1,ie))
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine maph1_to_l2t_wt(b,nb,a,na,if3d,w1,w2,ldw,wt)

!     Instead of interpolation we do a weighted projection

c     Input:   a
c     Output:  b
c
      real a(1),b(1),wt(1)
      real w1(1),w2(1)        ! work arrays
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
         nao = na       ! lx2
         nbo = nb       ! nx_crs
         call zwgl (za,w,na)
         call zwgll(zb,w,nb)
         call igllm(iba,ibat,zb,za,nb,na,nb,na)
      endif

!     na    - lx2
!     nb    - nx_crs      
!     iba   - interpolation matrix of Linear lagrange interpolants on
!             GLL grid to GL points
!     ibat  - transpose of iba      

      call specmpn_wt(b,nb,a,na,ibat,iba,if3d,w1,w2,ldw,wt)
c
      return
      end
c
c-----------------------------------------------------------------------

      subroutine specmpn_wt(b,nb,a,na,ba,ab,if3d,w1,w2,ldw,wt)
C
C     -  Spectral (weighted) projection from A to B via tensor products
C     -  scratch arrays: w(na*na*nb + nb*nb*na)
C
C     5/3/00  -- this routine replaces specmp in navier1.f, which
c                has a potential memory problem
C

      implicit none

      logical if3d

      integer nb,na,ldw

      real b(nb,nb,nb),a(na,na,na),wt(na,na,na)
      real w1(ldw),w2(ldw)       ! work

      real ba(1)        ! - interpolation matrix of Linear lagrange interpolants on
                        !   GLL grid to GL points (lx2,nx_crs)
      real ab(1)        ! - transpose of ab (nx_crs,lx2)

      integer ltest,nab,nbb
      integer k,l,iz
c

!     na    - lx2
!     nb    - nx_crs      
!
!     ab    - interpolation matrix of Linear lagrange interpolants on
!             GLL grid to GL points (lx2,nx_crs)
!     ba    - transpose of ab (nx_crs,lx2)
!
!     wt    - BM2      

      ltest = na*nb
      if (if3d) ltest = na*na*nb + nb*na*na
      if (ldw.lt.ltest) then
         write(6,*) 'ERROR specmp:',ldw,ltest,if3d
         call exitt
      endif

!     Collocate with weight      
      if (if3d) then
        call col3(w2,a,wt,na*na*na)
      else
        call col3(w2,a,wt,na*na*1)
      endif
c
      if (if3d) then
         nab = na*nb
         nbb = nb*nb
         call mxm(ba,nb,w2,na,w1,na*na)
         k=1
         l=na*na*nb + 1
         do iz=1,na
            call mxm(w1(k),nb,ab,na,w1(l),nb)
            k=k+nab
            l=l+nbb
         enddo
         l=na*na*nb + 1
         call mxm(w1(l),nbb,ab,na,b,nb)
      else
         call mxm(ba,nb,w2,na,w1,na)
         call mxm(w1,nb,ab,na,b,nb)
      endif

      return
      end

!-----------------------------------------------------------------------










