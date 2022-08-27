!     Major modifications for local_solve_fdm

!======================================================================
      subroutine local_solves_fdm3(u,v)
c
c     Given an input vector v, this returns the additive Schwarz solution
c
c                       -1
c                  u = M   v
c
c                      T     -1  
c     where M = sum ( R_i A_i    R_i )
c                i
c
c     The R_i's are simply index_set restriction operators.
c
c     The local solves are performed with the fast diagonalization method.
c
      include 'SIZE'
      include 'INPUT'
      include 'DOMAIN'
      include 'PARALLEL'
      include 'SOLN'          ! vmult
c
      include 'TSTEP'
      include 'CTIMER'
      include 'TEST'          ! temporary
c
      real u(lx2,ly2,lz2,lelv),v(lx2,ly2,lz2,lelv)
      common /scrpre/ v1(lx1,ly1,lz1,lelv)
     $               ,w1(lx1,ly1,lz1),w2(lx1,ly1,lz1)
      common /scrover/ ar(lelv)

      parameter(lxx=lx1*lx1, levb=lelv+lbelv)
      common /fastd/  df(lx1*ly1*lz1,levb)
     $             ,  sr(lxx*2,levb),ss(lxx*2,levb),st(lxx*2,levb)
      integer e,eb,eoff

      real wk1(lx1,ly1,lz1,lelv)
      real wk2(lx1,ly1,lz1,lelv)
      real wk3(lx1,ly1,lz1,lelv)

      integer nit


      nit = 1

c
      if (icalld.eq.0) tsolv=0.0
      icalld=icalld+1
      nsolv=icalld
c
      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv

      call copy(u,v,ntot2)

      do i=1,nit
!       Exchange data across elements
        call exchange_m2(v1,u,wk1,wk2,wk3)
        call copy(tmp1,v1,ntot1)
c
c       Now solve each subdomain problem:
c
        etime1=dnekclock()

        eoff  = 0
        if (ifield.gt.1) eoff  = nelv

        do e = 1,nelv
           eb = e + eoff
           call fastdm1(v1(1,1,1,e),df(1,eb)
     $                             ,sr(1,eb),ss(1,eb),st(1,eb),w1,w2)
        enddo
        tsolv=tsolv+dnekclock()-etime1
        call copy(tmp2,v1,ntot1)
c
c       Exchange/add elemental solutions
c

c       Copy back to pressure grid
c       lx1 -> lx2
        call extract_interior(u,v1)
      enddo  

c
      return
      end
c-----------------------------------------------------------------------


      subroutine local_solves_fdm2(u,v)
c
c     Given an input vector v, this returns the additive Schwarz solution
c
c                       -1
c                  u = M   v
c
c                      T     -1  
c     where M = sum ( R_i A_i    R_i )
c                i
c
c     The R_i's are simply index_set restriction operators.
c
c     The local solves are performed with the fast diagonalization method.
c
      include 'SIZE'
      include 'INPUT'
      include 'DOMAIN'
      include 'PARALLEL'
      include 'SOLN'          ! vmult
c
      include 'TSTEP'
      include 'CTIMER'
c
      real u(lx2,ly2,lz2,lelv),v(lx2,ly2,lz2,lelv)
      common /scrpre/ v1(lx1,ly1,lz1,lelv)
     $               ,w1(lx1,ly1,lz1),w2(lx1,ly1,lz1)
      common /scrover/ ar(lelv)

      parameter(lxx=lx1*lx1, levb=lelv+lbelv)
      common /fastd/  df(lx1*ly1*lz1,levb)
     $             ,  sr(lxx*2,levb),ss(lxx*2,levb),st(lxx*2,levb)
      integer e,eb,eoff

      real wk1(lx1,ly1,lz1,lelv)

c
      if (icalld.eq.0) tsolv=0.0
      icalld=icalld+1
      nsolv=icalld
c
      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv
c
!     Fill interiors
      call fill_interior(v1,v)
      call dface_ext    (v1)
      call dssum        (v1,lx1,ly1,lz1)
      call dface_add1si (v1,-1.)

      call copy(wk1,v1,ntot1)

c     Fill interiors
      call fill_interior(v1,v)
      call dface_ext    (v1)
      call crnr_ext     (v1)

      call dssum        (v1,lx1,ly1,lz1)
      call dface_add1si (v1,-1.)
      call crnr_corr    (v1)
c
c     Now solve each subdomain problem:
c
      etime1=dnekclock()

      eoff  = 0
      if (ifield.gt.1) eoff  = nelv

      do e = 1,nelv
         eb = e + eoff
         call fastdm1(v1(1,1,1,e),df(1,eb)
     $                           ,sr(1,eb),ss(1,eb),st(1,eb),w1,w2)
      enddo
      tsolv=tsolv+dnekclock()-etime1
c
c     Exchange/add elemental solutions
c
!      call s_face_to_int (v1,-1.)
!      call dssum         (v1,lx1,ly1,lz1)
!      call s_face_to_int (v1, 1.)
!      if(param(42).eq.0) call do_weight_op(v1)

c     Copy back to pressure grid
c     lx1 -> lx2
      call extract_interior(u,v1)

c
      return
      end
c-----------------------------------------------------------------------

      subroutine dd_solver2(u,v,h1,h2,h2inv,nit,wk1,wk2)

      implicit none

      include 'SIZE'
      include 'DOMAIN'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'CTIMER'

      real h1(1),h2(1),h2inv(1)
      real u(1),v(1),wk1(1),wk2(1)
      real uc
      common /scrprc/ uc(lx1*ly1*lz1*lelt)

      real alpha
      integer nit,i
      integer ntot
      integer intype
      real sigma
      real gl2norm2
      real nor
c
      if (icalld.eq.0) then
         tddsl=0.0
         tcrsl=0.0
         nddsl=0
         ncrsl=0
      endif
      icalld = icalld + 1
      nddsl  = nddsl  + 1
      ncrsl  = ncrsl  + 1

      ntot  = lx2*ly2*lz2*nelv
      call rzero(u,ntot)
      call copy(wk1,v,ntot)
      intype = 1         ! explicit
!      nit    = 50       ! Iterations
      etime1=dnekclock()
      sigma  = 1.00
      do i=1,nit
!        call local_solves_fdm(wk2,wk1)          ! wk2 = inv(M)*wk1
        call local_solves_fdm2(wk2,wk1)
        call add2s2(u,wk2,sigma,ntot)           ! u   = u + \sigma*wk2
        call cdabdtp(wk1,u,h1,h2,h2inv,intype)  ! wk1 = A u
        call sub2(wk1,v,ntot)                   ! wk1 = (Au - g)
        call chsign(wk1,ntot)                   ! wk1 = (g - Au)
      enddo  

      nor = gl2norm2(wk1,ntot)
      if (nio.eq.0) write(6,*) 'Residual Norm',nit,nor


      tddsl=tddsl+dnekclock()-etime1

      etime1=dnekclock()
!      call crs_solve_l2 (uc,v)
      tcrsl=tcrsl+dnekclock()-etime1

      alpha = 10.
c     if (param(89).ne.0.) alpha = abs(param(89))
!      call add2s2(u,uc,alpha,ntot)

      return
      end
c-----------------------------------------------------------------------
     
      subroutine mysemhat(a,b,c,d,z,dgll,dgllt,jgll,bgl,zgl,dgl,jgl,n,w)

      implicit none

c     Generate matrices for single element, 1D operators:
c
c     a        = Laplacian
c     b        = diagonal mass matrix
c     c        = convection operator b*d
c     d        = derivative matrix
c     dgll     = derivative matrix
c     dgllt    = derivative matrix transpose
c     jgll     = interpolation matrix 
c     z        = GLL points
c
c     zgl      = GL points
c     bgl      = diagonal mass matrix on GL
c     dgl      = derivative matrix,    mapping from velocity nodes to pressure
c     jgl      = interpolation matrix, mapping from velocity nodes to pressure
c
c     n        = polynomial degree (velocity space)
c     w        = work array of size 2*n+2


      integer i,j,k,n,np,nm,n2

      real a(0:n,0:n),b(0:n),c(0:n,0:n),d(0:n,0:n),z(0:n)
      real dgll(0:n,0:n),dgllt(0:n,0:n),jgll(0:n,0:n)
c
      real bgl(1:n-1),zgl(1:n-1)
      real dgl(1:n-1,0:n),jgl(1:n-1,0:n)
c
      real w(0:2*n+1)
c
      np = n+1    ! No of Points on M1 (lx1)
      nm = n-1    ! No of Points on M2 (lx2)
      n2 = n-2    ! Polynomial order of M2 (lx2-1)
c
      call zwgll (z,b,np)
c
      do i=0,n
         call fd_weights_full(z(i),z(0),n,1,w)
         do j=0,n
            d(i,j) = w(j+np)                   !  Derivative matrix
         enddo
      enddo

      if (n.eq.1) return                       !  No interpolation for n=1

!     Weights and Interpolation on GLL points      
      do i=0,n
!        Using n instead of n2
!        as polynomial order      
         call fd_weights_full(z(i),z(0),n,1,w(0))
         do j=0,n
            jgll(i,j) = w(j   )                  !  Interpolation matrix
            dgll(i,j) = w(j+np)                  !  Derivative    matrix
         enddo
      enddo

!     prabal
      call transpose(dgllt,np,dgll,np)
c
      call rzero(a,np*np)
      do j=0,n
      do i=0,n
         do k=0,n
            a(i,j) = a(i,j) + d(k,i)*b(k)*d(k,j)
         enddo
         c(i,j) = b(i)*d(i,j)
      enddo
      enddo
c
      call zwgl (zgl,bgl,nm)
c
      do i=1,n-1
         call fd_weights_full(zgl(i),z,n,1,w)
         do j=0,n
            jgl(i,j) = w(j   )                  !  Interpolation matrix
            dgl(i,j) = w(j+np)                  !  Derivative    matrix
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------

      subroutine swap_lengths_mesh2

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'WZ'

      real l
      common /swaplengths/ l(lx1,ly1,lz1,lelv)

      real lr, ls, lt   ! not used by swap lengths
      real llr(lelt)    ! length of left element along "r" 
      real lls(lelt)    ! length of left element along "s"
      real llt(lelt)    ! length of left element along "t"
      real lmr(lelt)    ! length of this element along "r"
      real lms(lelt)    ! length of this element along "s"
      real lmt(lelt)    ! length of this element along "t"
      real lrr(lelt)    ! length of right element along "r"
      real lrs(lelt)    ! length of right element along "s"
      real lrt(lelt)    ! length of right element along "t"
      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr,lls,llt
     $              , lmr,lms,lmt
     $              , lrr,lrs,lrt

      integer e,i,j,k,f,nfaces
      integer n,n2,nz0,nzn,nx

!     Coordinates
!     M2 -> M1 (Interior Points just copied) 
!     Face points from M1 mesh      
      real x21(lx1,ly1,lz1,lelt)
      real y21(lx1,ly1,lz1,lelt)
      real z21(lx1,ly1,lz1,lelt)

      common /scrsf/ x21,y21,z21

!     Does not matter
!     We get the size of the whole element      
!!     Fill interior values      
!      call fill_interior(x21,xm2)
!      call fill_interior(y21,ym2)
!      call rzero(z21,lx1*ly1*lz1*nelv)
!      if (if3d) call fill_interior(z21,zm2)
!
!!     Copy all face values      
!      nfaces = 2*ndim
!      do e=1,nelv
!      do f=1,nfaces
!        call facec(x21,xm1,e,f,nx1,ny1,nz1,nelv)
!        call facec(y21,ym1,e,f,nx1,ny1,nz1,nelv)
!        if (if3d) call facec(z21,zm1,e,f,nx1,ny1,nz1,nelv)
!      enddo
!      enddo

      call opcopy(x21,y21,z21,xm1,ym1,zm1)
      if (.not.if3d) call rzero(z21,lx1*ly1*lz1*nelv)

      n2 = lx1-1
      nz0 = 1
      nzn = 1
      nx  = lx1-2
      if (if3d) then
         nz0 = 0
         nzn = n2
      endif
      call plane_space(lmr,lms,lmt,0,n2,wxm1,x21,y21,z21,nx,n2,nz0,nzn)

      n=n2+1
      if (if3d) then
         do e=1,nelv
         do j=2,n2
         do k=2,n2
            l(1,k,j,e) = lmr(e)
            l(n,k,j,e) = lmr(e)
            l(k,1,j,e) = lms(e)
            l(k,n,j,e) = lms(e)
            l(k,j,1,e) = lmt(e)
            l(k,j,n,e) = lmt(e)
         enddo
         enddo
         enddo
         call dssum(l,n,n,n)
         do e=1,nelv
            llr(e) = l(1,2,2,e)-lmr(e)
            lrr(e) = l(n,2,2,e)-lmr(e)
            lls(e) = l(2,1,2,e)-lms(e)
            lrs(e) = l(2,n,2,e)-lms(e)
            llt(e) = l(2,2,1,e)-lmt(e)
            lrt(e) = l(2,2,n,e)-lmt(e)
         enddo
      else
         do e=1,nelv
         do j=2,n2
            l(1,j,1,e) = lmr(e)
            l(n,j,1,e) = lmr(e)
            l(j,1,1,e) = lms(e)
            l(j,n,1,e) = lms(e)
         enddo
         enddo
         call dssum(l,n,n,1)
         do e=1,nelv
            llr(e) = l(1,2,1,e)-lmr(e)
            lrr(e) = l(n,2,1,e)-lmr(e)
            lls(e) = l(2,1,1,e)-lms(e)
            lrs(e) = l(2,n,1,e)-lms(e)
         enddo
      endif
      return
      end
c----------------------------------------------------------------------

      subroutine gen_fast_again3(df,sr,ss,st)
c
c     Generate fast diagonalization matrices for each element
c
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'WZ'

      integer lxx
      parameter(lxx=lx1*lx1)
      real df(lx1*ly1*lz1,1),sr(lxx*2,1),ss(lxx*2,1),st(lxx*2,1)

      real lr ,ls ,lt   ! not used by swap lengths
      real llr(lelt)    ! length of left element along "r" 
      real lls(lelt)    ! length of left element along "s"
      real llt(lelt)    ! length of left element along "t"
      real lmr(lelt)    ! length of this element along "r"
      real lms(lelt)    ! length of this element along "s"
      real lmt(lelt)    ! length of this element along "t"
      real lrr(lelt)    ! length of right element along "r"
      real lrs(lelt)    ! length of right element along "s"
      real lrt(lelt)    ! length of right element along "t"
      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr,lls,llt
     $              , lmr,lms,lmt
     $              , lrr,lrs,lrt


!     Boundary conditions
      integer lbr,rbr,lbs,rbs,lbt,rbt,e

!     1D density along r,s,t      
      real fldr(lx1,lelv),flds(lx1,lelv),fldt(lx1,lelv)
      real fld(lx1,ly1,lz1,lelv)
      common /ctmpf_rho/ fldr,flds,fldt,fld

!     r,s,t coordinates of the 1D spectral elements.      
      real lr1(lx1,lelt),ls1(ly1,lelt),lt1(lz1,lelt)
!     Derivative mapping
      real drdx(lx1,lelt),dsdy(ly1,lelt),dtdz(lz1,lelt)
      common /geom1D/lr1,ls1,lt1,drdx,dsdy,dtdz

      integer ifld
      integer i,j,k,l,n,nr,ns,nt
      integer ierr,ierrmx,iglmax

      real vlmax
      real diag,eps
      

      ierr = 0

      if (param(44).eq.1) then
        if (nio.eq.0) then
          write(6,*) 'Density in FEM local solves not implemented'
          write(6,*) 'Exitting in gen_fast_again()'
        endif  

        call exitt
      endif

      ifld = 1
      n = lx1*ly1*lz1*nelv
      call copy(fld,vtrans(1,1,1,1,ifld),n)
      call get_1D_fld(fldr,flds,fldt,fld)
!     fldr -> Density along "r" for each element     
!     flds -> Density along "s" for each element     
!     fldt -> Density along "t" for each element     

!     lmr,drdx,...
!     lr1,ls1,lt1
!     drdx,dsdy,dtdz      
      call local_1D_geom()

      do e=1,nelv

         if (param(44).eq.1) then
!           call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,2,ierr)
         else
           call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,3,ierr)
         endif

         if (param(44).eq.1) then
!           call set_up_fast_1D_fem( sr(1,e),lr,nr ,lbr,rbr
!     $                      ,llr(e),lmr(e),lrr(e),zgm2(1,1),lx2,e)
         else
           call set_up_1D_gev( sr(1,e),lr,nr ,lbr,rbr
     $                 ,llr(e),lmr(e),lrr(e),fldr(1,e),drdx(1,e),e)
         endif
         if (ifaxis) then
           if (nio.eq.0) write(6,*) 'Not implemented for ifaxis yet.'
           if (nio.eq.0) write(6,*) 'Exitting in gen_fast_again3'
           call exitt

         else
            if (param(44).eq.1) then
!               call set_up_fast_1D_fem( ss(1,e),ls,ns ,lbs,rbs
!     $                      ,lls(e),lms(e),lrs(e),zgm2(1,2),ly2,e)
            else
               call set_up_1D_gev( ss(1,e),ls,ns ,lbs,rbs
     $                 ,lls(e),lms(e),lrs(e),flds(1,e),dsdy(1,e),e)
            endif
         endif
         if (if3d) then
            if (param(44).eq.1) then
!               call set_up_fast_1D_fem( st(1,e),lt,nt ,lbt,rbt
!     $                      ,llt(e),lmt(e),lrt(e),zgm2(1,3),lz2,e)
            else
               call set_up_1D_gev( st(1,e),lt,nt ,lbt,rbt
     $                 ,llt(e),lmt(e),lrt(e),fldt(1,e),dtdz(1,e),e)
            endif
         endif
c
c        Set up diagonal inverse
c
         if (if3d) then
            eps = 1.e-5 * (vlmax(lr(2),nr-2)
     $                  +  vlmax(ls(2),ns-2) + vlmax(lt(2),nt-2))
            l   = 1
            do k=1,nt
            do j=1,ns
            do i=1,nr
               diag = lr(i) + ls(j) + lt(k)
               if (diag.gt.eps) then
                  df(l,e) = 1.0/diag
               else
                  df(l,e) = 0.0
               endif
               l = l+1
            enddo
            enddo
            enddo
         else
            eps = 1.e-5*(vlmax(lr(2),nr-2) + vlmax(ls(2),ns-2))
            l   = 1
            do j=1,ns
            do i=1,nr
               diag = lr(i) + ls(j)
               if (diag.gt.eps) then
                  df(l,e) = 1.0/diag
               else
                  df(l,e) = 0.0
               endif
               l = l+1
            enddo
            enddo
         endif
c
c        Next element ....
c
      enddo

      ierrmx = iglmax(ierr,1)
      if (ierrmx.gt.0) then
         if (ierr.gt.0) write(6,*) nid,ierr,' BC FAIL'
         call exitti('E INVALID BC FOUND in genfast$',ierrmx)
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_1D_gev(s,lam,n,lbc,rbc,ll,lm,lr,
     $                                    rho,drdx,ie)

!     Calculate the Generalized Eigenvalues

      implicit none

      include 'SIZE'
      include 'MYSEMHAT'
      
      real g,w
      common /fast1dsem/ g(lr2),w(lr2)

      real s(1),lam(1),ll,lm,lr
      integer lbc,rbc
      
      integer bb0,bb1,eb0,eb1,n,n1
      logical l,r

      integer ie
      real rho(lx1)
      real drdx(lx1)
      integer i,j,j0      ! prabal

      l = (lbc.eq.0)
      r = (rbc.eq.0)
      n=lx1

!     calculate E tilde, B tilde operator
      call set_up_1D_op(s,g,bh,dph,drdx,jph,rho)


      call row_zero(s,n,n,1)
      call rzero(s(1),n)    ! Column zero
      s(1) = 1.0
      call row_zero(s,n,n,n)
      call rzero(s((n-1)*n + 1),n)    ! Column zero
      s(n*n) = 1.0

      call row_zero(g,n,n,1)
      call rzero(g(1),n)    ! Column zero
      g(1) = 1.0
      call row_zero(g,n,n,n)
      call rzero(g((n-1)*n + 1),n)    ! Column zero
      g(n*n) = 1.0

      do i=1,lx1
        write(6,12) 'lapM', (s(j),j=i,lx1*lx1,lx1)
      enddo

      do i=1,lx1
        write(6,12) 'MasM', (g(j),j=i,lx1*lx1,lx1)
      enddo

12    format(A4,2x,16(E12.5,2x))

!      if (l) then
!        call row_zero(s,n,n,1)
!        s(1) = 1.0
!        call row_zero(g,n,n,1)
!        g(1) = 1.0
!      endif  
!
!      if (r) then
!        call row_zero(s,n,n,n)
!        s(n*n) = 1.0
!        call row_zero(g,n,n,n)
!        g(n*n) = 1.0
!      endif  
       
      n=lx1
      call generalev(s,g,lam,n,w)
      if(.not.l) call row_zero(s,n,n,1)
      if(.not.r) call row_zero(s,n,n,n)
      call transpose(s(n*n+1),n,s,n) ! compute the transpose of s

      write(6,12) 'EigM', (lam(i),i=1,lx1) 


      return
      end
c-----------------------------------------------------------------------

      subroutine set_up_1D_op(s,g,bh,dgl,drdx,jgl,rho)

!                 -1  T
!     S = D B*(rho)  D
!
!              T
!     S = J B J
!
!     rho - density

      implicit none

      include 'SIZE'

      real s(lx1,lx1)
      real g(lx1,lx1)
      real bh(lx1)
      real dgl(lx1,lx1)
      real jgl(lx1,lx1)
      real drdx(lx1)

      real rho(lx1)       ! density
      real rhoi(lx1)      ! density inverse

      integer n
      integer i,j

      n=lx1

      call invers2(rhoi,rho,n)
      
      call build_1D_laplacian(s,dgl,drdx,bh,rhoi,n)
   
      call local_1D_mass(g,jgl,drdx,bh,n)


      return
      end
c-----------------------------------------------------------------------
      subroutine build_1D_laplacian(ah,dph,drdx,bm,rhoi,n)

!     This is the weak Laplacian Matrix
!     ah = DT*B*D 

      implicit none

      integer n         ! lx1/ly1/lz1
      real ah(n,n)      ! Laplacian
      real dph(n,n)     ! Gradient matrix
      real drdx(n) 
      real bm(n)        ! Mass (unweighted)
      real rhoi(n)

      integer i,j,k
      real s,s1

      logical ifbdry    ! Add boundary terms?


      do i=1,n
      do j=1,n
        s = 0.0
        do k=1,n
          s1 = dph(k,i)*drdx(k)*rhoi(k)*bm(k)*(1.0/drdx(k))*
     $         dph(k,j)*drdx(k)
          s  = s + s1
!          s = s + dph(k,i)*drdx(k)*bm(k)*(1.0/drdx(k))*dph(k,j)*drdx(k)
        enddo
        ah(i,j) = s
      enddo
      enddo 

      ifbdry = .false.
      do j=1,n
!       Left Boundary
        k  = 1
        s  = 1.0*rhoi(k)*bm(k)*(1.0/drdx(k))*
     $       dph(k,j)*drdx(k)
        ah(1,j) = ah(1,j) - s
!       Right Boundary        
        k  = n
        s  = 1.0*rhoi(k)*bm(k)*(1.0/drdx(k))*
     $       dph(k,j)*drdx(k)
        ah(n,j) = ah(n,j) - s
      enddo  

      return
      end subroutine build_1D_laplacian
!---------------------------------------------------------------------- 

      subroutine local_1D_mass(ah,jgl,drdx,bm,n)

!     This is the Local 1D Mass Matrix
!     Built along with the interpolation operator        
!     ah = J*B*JT

      implicit none

      integer n         ! lx1/ly1/lz1
      real ah(n,n)      ! Laplacian
      real jgl(n,n)     ! Interpolation matrix
      real drdx(n) 
      real bm(n)        ! Mass (unweighted)
      real mu(n)

      integer i,j,k
      real s

      call rzero(ah,n*n)
      do i=1,n
      do j=1,n
        s = 0.0
        do k=1,n
          s = s + jgl(k,i)*bm(k)*(1.0/drdx(k))*jgl(k,j)
!          s = s + jgl(k,i)*jgl(k,j)
        enddo
        ah(i,j) = s
      enddo
      enddo 


      return
      end
!---------------------------------------------------------------------- 

      subroutine local_1D_geom()

      implicit none

      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'MYSEMHAT'

      integer e,nr
      real sc,dd,l0

!     These are populated by swap_Lengths_mesh2
      real l
      common /swaplengths/ l(lx1,ly1,lz1,lelv)
      real lr ,ls ,lt   ! not used by swap lengths
      real llr(lelt)    ! length of left element along "r" 
      real lls(lelt)    ! length of left element along "s"
      real llt(lelt)    ! length of left element along "t"
      real lmr(lelt)    ! length of this element along "r"
      real lms(lelt)    ! length of this element along "s"
      real lmt(lelt)    ! length of this element along "t"
      real lrr(lelt)    ! length of right element along "r"
      real lrs(lelt)    ! length of right element along "s"
      real lrt(lelt)    ! length of right element along "t"
      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr,lls,llt
     $              , lmr,lms,lmt
     $              , lrr,lrs,lrt

!     r,s,t coordinates of the 1D spectral elements.      
      real lr1(lx1,lelt),ls1(ly1,lelt),lt1(lz1,lelt)
!     Derivative mapping
      real drdx(lx1,lelt),dsdy(ly1,lelt),dtdz(lz1,lelt)
      common /geom1D/lr1,ls1,lt1,drdx,dsdy,dtdz

!     ah          = Laplacian
!     bh          = diagonal mass matrix
!     ch          = convection operator b*d
!     dh          = derivative matrix
!     dph         = derivative matrix (Same as dh for now)
!     jph         = interpolation matrix 
!     z           = GLL points
!
!     zglhat      = GL points
!     bgl         = diagonal mass matrix on GL
!     dgl         = derivative matrix,    mapping from velocity nodes to pressure
!     jgl         = interpolation matrix, mapping from velocity nodes to pressure
!
!     nr          = polynomial degree (velocity space)
!     wh          = Work array

      nr = lx1-1
      call mysemhat(ah,bh,ch,dh,zh,dph,dpht,jph,bgl,
     $              zglhat,dgl,jgl,nr,wh)

!     populate llr,lmr,lmt...
!     using Mesh2 coordinate positions      
!      call swap_lengths_mesh2
      call swap_lengths
     
!     Build 1D x,y,z based on Mesh2
      do e=1,nelv
        lr1(1,e)   = zgm1(1,1)
        lr1(lx1,e) = zgm1(lx1,1)
        call copy(lr1(2,e),zgm2(1,1),lx2)
        sc = lmr(e)/2.0
        call cmult(lr1(1,e),sc,lx1)
!       Add left extension        
        sc = zgm2(lx2,1)-zgm1(lx1,1)
        dd = sc*llr(e)/2.0
        lr1(1,e) = lr1(1,e)+dd
!       Add right extension        
        sc = zgm2(1,1)-zgm1(1,1)
        dd = sc*lrr(e)/2.0
        lr1(lx1,e) = lr1(lx1,e)+dd

        ls1(1,e)   = zgm1(1,2)
        ls1(ly1,e) = zgm1(ly1,2)
        call copy(ls1(2,e),zgm2(1,2),ly2)
        sc = lms(e)/2.0
        call cmult(ls1(1,e),sc,ly1)
!       Add left extension        
        sc = zgm2(ly2,2)-zgm1(ly1,2)
        dd = sc*lls(e)/2.0
        ls1(1,e) = ls1(1,e)+dd
!       Add right extension        
        sc = zgm2(1,2)-zgm1(1,2)
        dd = sc*lrs(e)/2.0
        ls1(ly1,e) = ls1(ly1,e)+dd

        if (ndim.eq.3) then
          lt1(1,e)   = zgm1(1,3)
          lt1(lz1,e) = zgm1(lz1,3)
          call copy(lt1(2,e),zgm2(1,3),lz2)
          sc = lmt(e)/2.0
          call cmult(lt1(1,e),sc,lz1)
!         Add left extension        
          sc = zgm2(lz2,3)-zgm1(lz1,3)
          dd = sc*llt(e)/2.0
          lt1(1,e) = lt1(1,e)+dd
!         Add right extension        
          sc = zgm2(1,3)-zgm1(1,3)
          dd = sc*lrt(e)/2.0
          lt1(lz1,e) = lt1(lz1,e)+dd
        endif
      enddo 

!     Geometric factors/Jacobians      
      do e=1,nelv
        call mxm(dph,lx1,lr1(1,e),lx1,drdx(1,e),1)                ! dx/dr
        call mxm(dph,ly1,ls1(1,e),ly1,dsdy(1,e),1)                ! dy/ds
        if (ndim.eq.3) call mxm(dph,lz1,lt1(1,e),lz1,dtdz(1,e),1) ! dz/dt
      enddo
      call invcol1(drdx,lx1*nelv)                     ! dr/dx
      call invcol1(dsdy,ly1*nelv)                     ! ds/dy
      if (ndim.eq.3) call invcol1(dtdz,lz1*nelv)      ! dt/dz

      return
      end subroutine 
!---------------------------------------------------------------------- 


