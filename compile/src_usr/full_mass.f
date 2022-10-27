!======================================================================
!     Author: Prabal Negi      
!     Description: Evaluations using consistent integration
!     Routines:   fm_setup            : lx1 -> lxfm setup
!                 fm_intp             : Interpolate lx1 -> lxfm
!                 fm_tensor_op        : Apply op as tensor product      
!                 fm_tensor3_op       : opx,opy,opz as tensor operators
!                 fm_tensorx_op       : opx as tensor operator
!                 fm_tensory_op       : opy as tensor operator
!                 fm_tensorz_op       : opz as tensor operator
!                 fm_col_weights      : collocate with weights      
!                 fm_cdtp             : DT*x in the lxfm mesh.
!                 fm_multd            : D*u in the lxfm mesh.
!
!      
!     Functions:  vlsc2_fm            : Local inner product 
!                 glsc2_fm            : Global inner product      
!
!     Dependency: kopriva.f      
!====================================================================== 

      subroutine fm_setup()

      implicit none

!     Set up weights/interpolator from lx1 to lxfm mesh
!     For now I am using a fine Gauss-Legendre mesh since,
!     DGLLGL routine is available to build derivative operators.
!     Possibly could be made more flexible in the future.      

      include 'SIZE'
      include 'WZ'

      include 'FULLMASS'

      integer icalld
      save icalld
      data icalld /0/

      integer i,j,n,m

!     If already initialized      
      if (icalld.gt.0) return

!     Gauss-Legendre Mesh 
      call zwgl (fm_z,fm_wght,lxfm)

!!     Interpolator from M1 mesh to lxfm Mesh      
!      call igllm  (fm_jgl,fm_jglt,zgm1(1,1),fm_z,lx1,lxfm,lx1,lxfm)
!
!!     Derivative of M1 mesh variables on lxfm Mesh      
!      call dgllgl (fm_dgl,fm_dglt,zgm1(1,1),fm_z,fm_jgl,
!     $                                       lx1,lxfm,lx1,lxfm)
!!     Interpolator from M2 mesh to lxfm Mesh      
!      call iglm   (fm_jgl2,fm_jglt2,zgm2(1,1),fm_z,lx2,lxfm,lx2,lxfm)
!
!      Can't build Mesh 2 -> lxfm derivatives from speclib


!     Mesh 1 -> lxfm
      call BaryCentricWeights(fm_bw1,zgm1,lx1)

!     Interpolation Mesh 1 -> lxfm
      call PolynomialInterpolationMatrix
     $            (fm_jgl,fm_z,lxfm,zgm1,fm_bw1,lx1)
      call transpose(fm_jglt,lx1,fm_jgl,lxfm)

!     Derivative Mesh 1 -> lxfm 
      call LagrangeDerivativeMatrix(fm_dgl,fm_z,lxfm,zgm1,fm_bw1,lx1)
      call transpose(fm_dglt,lx1,fm_dgl,lxfm)


!     Aparently we don't seem to have a DGLGL sort of routine
!     So we are doing it through the barycentric algorithms
!     in kopriva.f

!     Mesh 2 -> lxfm      
      call BaryCentricWeights(fm_bw2,zgm2,lx2)

!     Interpolation Mesh 2 -> lxfm
      call PolynomialInterpolationMatrix
     $            (fm_jgl2,fm_z,lxfm,zgm2,fm_bw2,lx2)
      call transpose(fm_jglt2,lx2,fm_jgl2,lxfm)

!     Derivative Mesh 2 -> lxfm 
      call LagrangeDerivativeMatrix(fm_dgl2,fm_z,lxfm,zgm2,fm_bw2,lx2)
      call transpose(fm_dglt2,lx2,fm_dgl2,lxfm)

      icalld = icalld+1

      return
      end
c-----------------------------------------------------------------------

      subroutine fm_intp(fldf,fld,imsh)

      implicit none

      include 'SIZE'
      include 'INPUT'         ! if3d
      include 'GEOM'

      include 'FULLMASS'

      real fld  (1)
      real fldf (lxfm**ldim)

      integer imsh

      integer iz,i1,i2
      integer ldw

      real fm_op_wk(lxwk**ldim,2)
      common /fm_op_work/ fm_op_wk

      ldw = (lxwk**ldim)*2

      if (imsh.eq.1) then
!       Velocity to lxfm         
        call specmpn(fldf,lxfm,fld,lx1,fm_jgl,fm_jglt,
     $               if3d,fm_op_wk,ldw)
      elseif (imsh.eq.2) then
!       Pressure to lxfm        
        call specmpn(fldf,lxfm,fld,lx2,fm_jgl2,fm_jglt2,
     $               if3d,fm_op_wk,ldw)
      else
        if (nio.eq.0) write(6,*) 'fm_intp: Unknown imsh', imsh
        call exitt  
      endif

      return
      end subroutine fm_intp
!---------------------------------------------------------------------- 
      subroutine fm_tensor_op(fldf,fld,nx,ny,nz,op1d,op1dt,idir)

      implicit none

      include 'SIZE'
      include 'INPUT'         ! if3d
      include 'GEOM'

      include 'FULLMASS'

      integer nx,ny,nz

      real fld  (nx*ny*nz)
      real fldf (lxfm**ldim)

!     One dimensional operator      
      real op1d(lxfm*nx)
!     Transpose of the one dimensional operator      
      real op1dt(nx*lxfm)

      integer idir
      integer ldw

      real fm_op_wk(lxwk**ldim,2)
      common /fm_op_work/ fm_op_wk

      ldw = (lxfm**ldim)*2

      if (idir.eq.0) then
        call specmpn(fldf,lxfm,fld,nx,op1d,op1dt,if3d,fm_op_wk,ldw)
      elseif (idir.eq.1) then
        call specmpn(fld,nx,fldf,lxfm,op1dt,op1d,if3d,fm_op_wk,ldw)
      else
        if (nio.eq.0) write(6,*) 'fm_tensor_op: Invalid idir:', idir
        call exitt 
      endif

      return
      end subroutine fm_tensor_op
!---------------------------------------------------------------------- 
      subroutine fm_tensor3_op(fldf,fld,nx,ny,nz,opx,opyt,opzt,mx,my,mz)

      implicit none

      include 'SIZE'
      include 'INPUT'         ! if3d
      include 'GEOM'

      include 'FULLMASS'

      integer nx,ny,nz
      integer mx  ! No. of rows of opx
      integer my  ! No. of columns of opyt
      integer mz  ! No. of columns of opzt

      real fld  (nx*ny*nz)
      real fldf (mx*my*mz)

!     x operator is applied direction      
      real opx(mx*nx)

!     y,z operators apply as transposed operators 
      real opyt(ny*my)
      real opzt(nz*mz)

      real fm_op_wk(lxwk**ldim,2)
      common /fm_op_work/ fm_op_wk

      if ((mx*my*mz.gt.(lxwk**ldim)) .and. 
     $    (nx*ny*nz.gt.(lxwk**ldim))) then
        if (nio.eq.0) write(6,*) 
     $    'Inadequate workspace size in fm_tensor3_op'
        call exitt
      endif  

      call fm_tensorx_op(fm_op_wk(1,1),fld,nx,ny,nz,opx,mx)
      if (if3d) then
        call fm_tensory_op(fm_op_wk(1,2),fm_op_wk(1,1),mx,ny,nz,opyt,my)
        call fm_tensorz_op(fldf,fm_op_wk(1,2),mx,my,nz,opzt,mz)
      else
        call fm_tensory_op(fldf,fm_op_wk(1,1),mx,ny,nz,opyt,my)
      endif  


      return
      end subroutine fm_tensor3_op
!---------------------------------------------------------------------- 

      subroutine fm_tensorx_op(fldo,fld,nx,ny,nz,opx,mx)

      implicit none

      include 'SIZE'
      include 'INPUT'         ! if3d

      integer nx,ny,nz
      integer mx              ! no of rows

      real fld  (1)
      real fldo (1)

      real opx(mx*nx)

      call mxm  (opx,mx,fld,nx,fldo,ny*nz)

      return
      end subroutine fm_tensorx_op
!---------------------------------------------------------------------- 

      subroutine fm_tensory_op(fldo,fld,nx,ny,nz,opy,my)

      implicit none

      include 'SIZE'
      include 'INPUT'         ! if3d

      integer nx,ny,nz
      integer my              ! no of columns

      real fld  (1)
      real fldo (1)

      real opy(ny*my)
      integer i,j,k

      if (if3d) then
        i = 1
        j = 1
        do k = 1,nz
          call mxm  (fld(i),nx,opy,ny,fldo(j),my)
          i = i + nx*ny
          j = j + nx*my
        enddo  
      else
        call mxm  (fld,nx,opy,ny,fldo,my)
      endif  


      return
      end subroutine fm_tensory_op
!---------------------------------------------------------------------- 

      subroutine fm_tensorz_op(fldo,fld,nx,ny,nz,opz,mz)

      implicit none

      include 'SIZE'
      include 'INPUT'         ! if3d

      integer nx,ny,nz
      integer mz              ! no of columns

      real fld  (1)
      real fldo (1)

      real opz(nz*mz)

      if (if3d) then
        call mxm  (fld,nx*ny,opz,nz,fldo,mz)
      endif  


      return
      end subroutine fm_tensorz_op
!---------------------------------------------------------------------- 

      function vlsc2_fm(u,v)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'GEOM'

      include 'FULLMASS'

      real vlsc2_fm
      real s
      real u(lx1*ly1*lz1,lelt)
      real v(lx1*ly1*lz1,lelt)

      integer e,i,j,k,ii
      integer lxfm3
      integer imsh

      real vlsum

      call fm_setup()

      lxfm3 = lxfm**ndim

      s = 0.0
      imsh = 1
      do e=1,nelv

!       u_lx1 -> u_lxfm          
        call fm_intp(fm_wk3,u(1,e),imsh) 

!       v_lx1 -> v_lxfm
        call fm_intp(fm_wk4,v(1,e),imsh)
        call col2(fm_wk3,fm_wk4,lxfm3)

!       J_lx1 -> J_lxfm
        call fm_intp(fm_wk4,jacm1(1,1,1,e),imsh)
        call col2(fm_wk3,fm_wk4,lxfm3)

!       W*J*\rho*u
        call fm_col_weights(fm_wk3)

        s = s + vlsum(fm_wk3,lxfm3) 

      enddo
      
      vlsc2_fm = s


      return
      end function vlsc2_fm
!---------------------------------------------------------------------- 

      function glsc2_fm(u,v)

      implicit none

      include 'SIZE'

      real u(1),v(1)
      real sc,tmp
      real glsc2_fm
      real vlsc2_fm

      sc = vlsc2_fm(u,v)
      call gop(sc,tmp,'+  ',1)

      glsc2_fm = sc

      return 
      end function glsc2_fm 

!---------------------------------------------------------------------- 

      subroutine fm_col_weights(x)

      implicit none

      include 'SIZE'
      include 'FULLMASS'

      real x(1)

      integer i,j,k,ii

!     W*x
      ii = 0
      if (ndim.eq.2) then
        do j=1,lyfm
        do i=1,lxfm
          ii = ii + 1
          x(ii) = x(ii)*fm_wght(i)*fm_wght(j)
        enddo
        enddo
      else          
        do k=1,lzfm
        do j=1,lyfm
        do i=1,lxfm
          ii = ii + 1
          x(ii) = x(ii)*fm_wght(i)*fm_wght(j)*fm_wght(k)
        enddo
        enddo
        enddo
      endif  

      return
      end subroutine

!----------------------------------------------------------------------         

      subroutine fm_cdtp (dtx,x,rm1,sm1,tm1,isd)

!     Compute DT*X (entire field)
!     Evaluated on the lxfm mesh 

      implicit none

      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'IXYZ'
      include 'GEOM'
      include 'MASS'
      include 'INPUT'
      include 'ESOLV'
      include 'CTIMER'

      include 'FULLMASS'

      real dtx  (lx1*ly1*lz1,lelv)
      real x    (lx2*ly2*lz2,lelv)

      real rm1  (lx1*ly1*lz1,lelv)
      real sm1  (lx1*ly1*lz1,lelv)
      real tm1  (lx1*ly1*lz1,lelv)

      real wx,ta1,ta2,ta3
      common /ctmp1/ wx  (lx1*ly1*lz1)
     $ ,             ta1 (lx1*ly1*lz1)
     $ ,             ta2 (lx1*ly1*lz1)
     $ ,             ta3 (lx1*ly1,lz1)

      integer isd
      integer e
      integer nxyz1
      integer imsh1,imsh2

#ifdef TIMER
      if (icalld.eq.0) tcdtp=0.0
      icalld=icalld+1
      ncdtp=icalld
      etime1=dnekclock()
#endif

      if (ifaxis) then
        if (nio.eq.0) 
     $    write(6,*) 'FM_CDTp not implemented for ifaxis.'
        call exitt
      endif

      if(ifsplit) then
        if (nio.eq.0) 
     $    write(6,*) 'FM_CDTp not implemented for ifsplit.'
        call exitt
      endif


      nxyz1 = lx1*ly1*lz1
      imsh1 = 1
      imsh2 = 2

      do e=1,nelv

!       Interpolate x to lxfm Mesh
!       wk1 = x            
        call fm_intp(fm_wk1,x(1,e),imsh2)         
!       Collocate with weights
!       Jacobian goes away due to inverse jacobian of the dx/dr etc.
!       wk1 = W*x
        call fm_col_weights(fm_wk1)

!       wk2 = dr/dx_i
        call fm_intp(fm_wk2,rm1(1,e),imsh1)
!       wk3 = W*x*dr/dx
        call col3(fm_wk3,fm_wk1,fm_wk2,lxfm*lyfm*lzfm)
!       dtx = (dv/dr)*(W*x*dr/dx_i)
        call fm_tensor3_op(dtx(1,e),fm_wk3,lxfm,lyfm,lzfm,
     $       fm_dglt,fm_jgl,fm_jgl,lx1,lx1,lx1)

!       wk2 = ds/dx_i
        call fm_intp(fm_wk2,sm1(1,e),imsh1)
!       wk3 = W*x*ds/dx_i
        call col3(fm_wk3,fm_wk1,fm_wk2,lxfm*lyfm*lzfm)
!       ta1 = (dv/ds)*(W*x*ds/dx_i)
        call fm_tensor3_op(ta1,fm_wk3,lxfm,lyfm,lzfm,
     $       fm_jglt,fm_dgl,fm_jgl,lx1,lx1,lx1)
        
!       dtx = dtx + ta1
        call add2(dtx(1,e),ta1,nxyz1)

        if (ndim.eq.3) then
!         wk2 = dt/dx_i
          call fm_intp(fm_wk2,tm1(1,e),imsh1)
!         wk3 = W*x*dt/dx_i
          call col3(fm_wk3,fm_wk1,fm_wk2,lxfm*lyfm*lzfm)
!         ta1 = (dv/dt)*(W*x*dt/dx_i)
          call fm_tensor3_op(ta1,fm_wk3,lxfm,lyfm,lzfm,
     $         fm_jglt,fm_jgl,fm_dgl,lx1,lx1,lx1)

!         dtx = dtx + ta1
          call add2(dtx(1,e),ta1,nxyz1)
        endif  

!       If axisymmetric, add an extra diagonal term in the radial 
!       direction (only if solving the momentum equations and ISD=2)
!       NOTE: lz1=lz2=1

      enddo
!
#ifdef TIMER
      tcdtp=tcdtp+(dnekclock()-etime1)
#endif
      return
      end subroutine fm_cdtp
!---------------------------------------------------------------------- 

      subroutine fm_multd (du,u,rm1,sm1,tm1,isd)

!     Compute D*U (on Mesh 1)
!     U    : input variable, defined on M1
!     DU   : output variable, defined on M2
!     Integration done on the lxfm mesh        
!     RM1 : RXM1, RYM1 or RZM1
!     SM1 : SXM1, SYM1 or SZM1
!     TM1 : TXM1, TYM1 or TZM1
!     ISD : spatial direction (x=1,y=2,z=3)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'ESOLV'
      include 'CTIMER'

      include 'FULLMASS'

      real           du   (lx2*ly2*lz2,lelv)
      real           u    (lx1*ly1*lz1,lelv)
      real           rm1  (lx1*ly1*lz1,lelv)
      real           sm1  (lx1*ly1*lz1,lelv)
      real           tm1  (lx1*ly1*lz1,lelv)

      integer isd

      integer e
      integer imsh1


#ifdef TIMER
      if (icalld.eq.0) tmltd=0.0
      icalld=icalld+1
      nmltd=icalld
      etime1=dnekclock()
#endif

      if (ifaxis) then
        if (nio.eq.0) write(6,*) 
     $    'FM_MULTD not implemented for ifaxis.'
        call exitt
      endif

      imsh1 = 1
      do e=1,nelv

!       wk1 = du/dr
        call fm_tensor3_op(fm_wk1,u(1,e),lx1,ly1,lz1,
     $       fm_dgl,fm_jglt,fm_jglt,lxfm,lyfm,lzfm)

!       wk2 = (dr/dx_i)
        call fm_intp(fm_wk2,rm1(1,e),imsh1)
!       wk1 = (du/dr)*(dr/dx_i)
        call col2(fm_wk1,fm_wk2,lxfm*lyfm*lzfm)

!       wk3 = du/ds
        call fm_tensor3_op(fm_wk3,u(1,e),lx1,ly1,lz1,
     $       fm_jgl,fm_dglt,fm_jglt,lxfm,lyfm,lzfm)
!       wk2 = (ds/dx_i)
        call fm_intp(fm_wk2,sm1(1,e),imsh1)
!       wk3 = (du/ds)*(ds/dx_i)
        call col2(fm_wk3,fm_wk2,lxfm*lyfm*lzfm)

!       wk1 = wk1 + (ds/dx_i)*(du/dt)
        call add2(fm_wk1,fm_wk3,lxfm*lyfm*lzfm)

        if (ldim.eq.3) then
!         wk3 = du/dt
          call fm_tensor3_op(fm_wk3,u(1,e),lx1,ly1,lz1,
     $         fm_jgl,fm_jglt,fm_dglt,lxfm,lyfm,lzfm)
!         wk2 = (dt/dx_i)
          call fm_intp(fm_wk2,tm1(1,e),imsh1)
!         wk3 = (du/dt)*(dt/dx_i)
          call col2(fm_wk3,fm_wk2,lxfm*lyfm*lzfm)

!         wk1 = wk1 + (dt/dx_i)*(du/dt)
          call add2(fm_wk1,fm_wk3,lxfm*lyfm*lzfm)
        endif  

!       Collocate with the weights on the pressure mesh
!       wk1 = W*(du/dx_i)
        call fm_col_weights(fm_wk1)

!       du = q*W*(du/dx_i)
        call fm_tensor3_op(du(1,e),fm_wk1,lxfm,lyfm,lzfm,
     $         fm_jglt2,fm_jgl2,fm_jgl2,lx2,ly2,lz2)


      enddo

#ifdef TIMER
      tmltd=tmltd+(dnekclock()-etime1)
#endif
      return
      end subroutine fm_multd
!-----------------------------------------------------------------------

      subroutine fm_opdiv(outfld,inpx,inpy,inpz)

!     Compute OUTFLD = SUMi Di*INPi, 
!     the divergence of the vector field (INPX,INPY,INPZ)
!     Integrated on the lxfm mesh        

      implicit none

      include 'SIZE'
      include 'GEOM'

      real outfld (lx2,ly2,lz2,1)
      real inpx   (lx1,ly1,lz1,1)
      real inpy   (lx1,ly1,lz1,1)
      real inpz   (lx1,ly1,lz1,1)

      real work
      common /ctmp0/ work (lx2,ly2,lz2,lelv)
      
      integer iflg,ntot2

      iflg = 1

      ntot2 = lx2*ly2*lz2*nelv
      call fm_multd (work,inpx,rxm1,sxm1,txm1,1)
      call copy  (outfld,work,ntot2)
      call fm_multd (work,inpy,rym1,sym1,tym1,2)
      call add2  (outfld,work,ntot2)
      if (ldim.eq.3) then
         call fm_multd (work,inpz,rzm1,szm1,tzm1,3)
         call add2  (outfld,work,ntot2)
      endif

      return
      end

!-----------------------------------------------------------------------

      subroutine fm_opgradt(outx,outy,outz,inpfld)

!     Compute DTx, DTy, DTz of an input field INPFLD
!     Evaluated on lxfm 

      implicit none

      include 'SIZE'
      include 'GEOM'

      real outx   (lx1,ly1,lz1,1)
      real outy   (lx1,ly1,lz1,1)
      real outz   (lx1,ly1,lz1,1)
      real inpfld (lx2,ly2,lz2,1)

      call fm_cdtp (outx,inpfld,rxm1,sxm1,txm1,1)
      call fm_cdtp (outy,inpfld,rym1,sym1,tym1,2)
      if (ldim.eq.3) 
     $   call fm_cdtp (outz,inpfld,rzm1,szm1,tzm1,3)

      return
      end
!-----------------------------------------------------------------------
      subroutine fm_cdabdtp (ap,wp,h1,h2,h2inv,intype)

C     INTYPE= 0  Compute the matrix-vector product    DA(-1)DT*p
C     INTYPE= 1  Compute the matrix-vector product    D(B/DT)(-1)DT*p
C     INTYPE=-1  Compute the matrix-vector product    D(A+B/DT)(-1)DT*p

      implicit none

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'         ! ifanls

      real           ap    (lx2,ly2,lz2,1)
      real           wp    (lx2,ly2,lz2,1)
      real           h1    (lx1,ly1,lz1,1)
      real           h2    (lx1,ly1,lz1,1)
      real           h2inv (lx1,ly1,lz1,1)

      real ta1,ta2,ta3,tb1,tb2,tb3
      common /scrns/ ta1 (lx1,ly1,lz1,lelv)
     $ ,             ta2 (lx1,ly1,lz1,lelv)
     $ ,             ta3 (lx1,ly1,lz1,lelv)
     $ ,             tb1 (lx1,ly1,lz1,lelv)
     $ ,             tb2 (lx1,ly1,lz1,lelv)
     $ ,             tb3 (lx1,ly1,lz1,lelv)

      integer intype
      real tolhin,dtbdi

      call fm_opgradt(ta1,ta2,ta3,wp)
      if ((intype.eq.0).or.(intype.eq.-1)) then
        tolhin=tolhs
        call ophinv (tb1,tb2,tb3,ta1,ta2,ta3,h1,h2,tolhin,nmxv)
      else
        if (ifanls) then
          dtbdi = dt/bd(1)   ! scale by dt*backwd-diff coefficient
          call opbinv1(tb1,tb2,tb3,ta1,ta2,ta3,dtbdi)
        else
          call opbinv (tb1,tb2,tb3,ta1,ta2,ta3,h2inv)
        endif
      endif
      call fm_opdiv (ap,tb1,tb2,tb3)

      return
      end
C
C-----------------------------------------------------------------------




















