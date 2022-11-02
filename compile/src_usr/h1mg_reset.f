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

      real*8 t0
      real*8 dnekclock

      t0 = dnekclock()

      n = lx1*ly1*lz1*nelv

!     The multigrid has not been modified yet
      param(43) = 1
      ifmgrid   = .false.

      ifield = 1
      if (uparam(5).eq.1.0) then
        call set_overlap_again
      endif 

      if (uparam(6).eq.1.0) then
        call invers2(h1,vtrans(1,1,1,1,ifield),n)
!        call rone(h1,n)    
        call rzero(h2,n)
        call set_up_h1_crs_again(h1,h2)
      endif  
      t0 = dnekclock() - t0
      if (nio.eq.0) then
         write(6,*) 'Preconditioner Resets ',t0, ' sec'
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
      real a(lcr*lcr*lelv)                      ! Could use some common
                                                ! block for this

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

      ncr     = nxc**ldim
      nxyz_c  = ncr

c     Set SEM_to_GLOB

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
         call copy  (h1,h1_usr,n)
         call copy  (h2,h2_usr,n)
         call get_local_crs_galerkin(a,ncr,nxc,h1,h2,w1,w2)
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













