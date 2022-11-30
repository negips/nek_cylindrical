!======================================================================
!
!     Author: Prabal Negi
!     Description: Cylindrical coordinates solver
!     
!     Routines:
!     esolver_cyl       : Choose E solver
!     uzawa_gmres_cyl   : GMRES solver
!     uzawa_cyl         : CG solver
!     cggosf_cyl        : CG solver (velocity)
!     setprec_cyl       : Diagonal Preconditioner (Velocity)       
!      
!
!     Outside dependencies: 
!     generic_subs.f          : ortho_subspace()     
!
!      
!======================================================================
!-----------------------------------------------------------------------
      subroutine esolver_cyl (res,h1,h2,h2inv,intype)

!     Choose E-solver

      INCLUDE 'SIZE'
      INCLUDE 'ESOLV'
      INCLUDE 'INPUT'
      include 'CTIMER'
C
      real res   (lx2,ly2,lz2,lelv)
      real h1    (lx1,ly1,lz1,lelv)
      real h2    (lx1,ly1,lz1,lelv)
      real h2inv (lx1,ly1,lz1,lelv)
      common /scruz/ wk1(lx2*ly2*lz2*lelv)
     $             , wk2(lx2*ly2*lz2*lelv)
     $             , wk3(lx2*ly2*lz2*lelv)

      integer ig

      integer igmres

      if (icalld.eq.0) teslv=0.0

!     Moved to inside gmres      
!      call ortho(res) !Ensure that residual is orthogonal to null space

      icalld=icalld+1
      neslv=icalld
      etime1=dnekclock()

      if (.not. ifsplit) then
        if (param(42).eq.1) then
          call uzawa_cyl(res,h1,h2,h2inv,intype,icg)
        else
          call uzawa_gmres_cyl(res,h1,h2,h2inv,intype,icg)
          if (ig.gt.6) then 
            write(6,*) 'Unknown GMRES. exitting in esolver_new()'
            call exitt
          endif  
        endif
      else
        write(6,*) 'error: e-solver does not exist pnpn'
        call exitt
      endif

      teslv=teslv+(dnekclock()-etime1)

      RETURN
      END
c-----------------------------------------------------------------------

      subroutine uzawa_gmres_cyl(res,h1,h2,h2inv,intype,iter)

c     Solve the cylindrical pressure equation by right-preconditioned 
c     GMRES iteration.
c     intype =  0  (steady)
c     intype =  1  (explicit)
c     intype = -1  (implicit)

      implicit none

      include 'SIZE'
      include 'INPUT'   ! param
      include 'MASS'    ! bm2
      include 'TSTEP'
      include 'GMRES'

      integer lt1,lt2
      parameter(lt1 = lx1*ly1*lz1*lelv)
      parameter(lt2 = lx2*ly2*lz2*lelv)

      real divex
      common  /ctolpr/ divex
      
      logical          ifprint
      common  /cprint/ ifprint

      real             res  (lt2)
      real             h1   (lt1)
      real             h2   (lt1)
      real             h2inv(lt1)
      
      real wp
      common /scrmg/   wp (lt2)

      real wk1,wk2
      common /ctmp0/   wk1(lgmres),wk2(lgmres)

      real y
      common /cgmres1/ y(lgmres)

      real alpha, l, temp
      integer j,m

      logical iflag
      save    iflag
      data    iflag /.false./
      real    norm_fac
      save    norm_fac

      real*8 etime1,dnekclock

      integer ntot2
      real glsc2,glsc3,vlsc2,vlsc3
      integer iconv,intype
      real tolpss,div0
      integer i,k,iter
      real etime2,etime_p,ratio,rnorm

      integer maxiter

      logical ifwgt           ! If Weighted orthogonalization
      integer ngs             ! No of Gram-Schmid

      ifwgt = .false.
      ngs   = 1

      ntot2  = lx2*ly2*lz2*nelv

!     I've removed ortho from earlier calls and call it here at
!     the beginning      
      call ortho(res)
c
      if(.not.iflag) then
         iflag=.true.
         call uzawa_gmres_split0(ml_gmres,mu_gmres,bm2,bm2inv,
     $                           lx2*ly2*lz2*nelv)
         norm_fac = 1./sqrt(volvm2)
      endif
c
      etime1 = dnekclock()
      etime_p = 0.
      divex = 0.
      iter  = 0
      m = lgmres
c
      call chktcg2(tolps,res,iconv)
      if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
     $   tolps = abs(param(21))
!      if (istep.eq.0) tolps = 1.e-4
      tolpss = tolps


      iconv = 0
      call rzero(x_gmres,ntot2)

      do while(iconv.eq.0.and.iter.lt.1000)

         if(iter.eq.0) then
                                                        !      -1
            call col3(r_gmres,ml_gmres,res,ntot2)       ! r = L  res
         else
            !update residual
            call copy(r_gmres,res,ntot2)                      ! r = res
            call cdabdtp_cyl(w_gmres,x_gmres,h1,h2,h2inv,intype)  ! w = A x
            call add2s2(r_gmres,w_gmres,-1.,ntot2)            ! r = r - w
                                                              !      -1
            call col2(r_gmres,ml_gmres,ntot2)                 ! r = L   r
         endif
                                                            !            ______
         gamma_gmres(1) = sqrt(glsc2(r_gmres,r_gmres,ntot2))! gamma  = \/ (r,r) 
                                                            !      1
         if(iter.eq.0) then
            div0 = gamma_gmres(1)*norm_fac
            if (param(21).lt.0) tolpss=abs(param(21))*div0
         endif

         !check for lucky convergence
         rnorm = 0.
         if(gamma_gmres(1) .eq. 0.) goto 9000
         temp = 1./gamma_gmres(1)
         call cmult2(v_gmres(1,1),r_gmres,temp,ntot2)! v  = r / gamma
                                                     !  1            1
         do j=1,m
            iter = iter+1
                                                           !       -1
            call col3(w_gmres,mu_gmres,v_gmres(1,j),ntot2) ! w  = U   v
                                                           !           j
            
            etime2 = dnekclock()
            if(param(43).eq.1) then
               call uzprec(z_gmres(1,j),w_gmres,h1,h2,intype,wp)
            else                                        !       -1
               call hsmg_solve(z_gmres(1,j),w_gmres)    ! z  = M   w
            endif     
            etime_p = etime_p + dnekclock()-etime2
     
            call cdabdtp_cyl(w_gmres,z_gmres(1,j),    ! w = A z
     $                   h1,h2,h2inv,intype)      !        j
     
                                                  !      -1
            call col2(w_gmres,ml_gmres,ntot2)     ! w = L   w

!           Gram-Schmidt:
            call ortho_subspace(w_gmres,ntot2,h_gmres(1,j),v_gmres,
     $            lt2,j,bm2,ifwgt,ngs,wk1,wk2)

!           Apply Givens rotations to new column
            do i=1,j-1
              temp = h_gmres(i,j)                   
              h_gmres(i  ,j) =  c_gmres(i)*temp 
     $                         +s_gmres(i)*h_gmres(i+1,j)  
              h_gmres(i+1,j) = -s_gmres(i)*temp 
     $                         +c_gmres(i)*h_gmres(i+1,j)
            enddo
                                                              !            ______
            alpha = sqrt(glsc2(w_gmres,w_gmres,ntot2))        ! alpha =  \/ (w,w)
            rnorm = 0.
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h_gmres(j,j)*h_gmres(j,j)+alpha*alpha)
            temp = 1./l
            c_gmres(j) = h_gmres(j,j) * temp
            s_gmres(j) = alpha  * temp
            h_gmres(j,j) = l
            gamma_gmres(j+1) = -s_gmres(j) * gamma_gmres(j)
            gamma_gmres(j)   =  c_gmres(j) * gamma_gmres(j)

            
            rnorm = abs(gamma_gmres(j+1))*norm_fac
            ratio = rnorm/div0
            if (ifprint.and.nio.eq.0) 
     $         write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66       format(i5,1p4e12.5,i8,' Divergence')

#ifndef FIXITER
            if (rnorm .lt. tolpss) goto 900  !converged
#else
            if (iter.gt.param(151)-1) goto 900
#endif
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v_gmres(1,j+1),w_gmres,temp,ntot2) ! v    = w / alpha
                                                           !  j+1            
         enddo          ! j=1,m

  900    iconv = 1
 1000    continue
         !back substitution
         !     -1
         !c = H   gamma
         do k=j,1,-1
            temp = gamma_gmres(k)
            do i=j,k+1,-1
               temp = temp - h_gmres(k,i)*c_gmres(i)
            enddo
            c_gmres(k) = temp/h_gmres(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x_gmres,z_gmres(1,i),c_gmres(i),ntot2) 
                       ! x = x + c  z
                       !          i  i
         enddo

      enddo           ! while(iconv.eq.0.and.iter.lt.1000)
 9000 continue
c
      divex = rnorm

      call copy(res,x_gmres,ntot2)

      call ortho (res)  ! Orthogonalize wrt null space, if present

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep,'  U-PRES gmres (CYL) ', 
     &                            iter,divex,div0,tolpss,etime_p,etime1
c     call flush_hack
 9999 format(i11,a,I6,1p5e13.4)

      return
      end

c-----------------------------------------------------------------------
      subroutine uzawa_cyl (rcg,h1,h2,h2inv,intype,iter)

C     Solve the pressure equation by (nested) preconditioned 
C     conjugate gradient iteration.
C     INTYPE =  0  (steady)
C     INTYPE =  1  (explicit)
C     INTYPE = -1  (implicit)
      
      implicit none

      include 'SIZE'
      include 'MASS'
      include 'INPUT'
      include 'PARALLEL'
      include 'TSTEP' 

      integer lt1,lt2
      parameter(lt1 = lx1*ly1*lz1*lelv)
      parameter(lt2 = lx2*ly2*lz2*lelv)

      real divex
      common  /ctolpr/ divex

      logical          ifprint
      common  /cprint/ ifprint
      real             rcg  (lt2)
      real             h1   (lt1)
      real             h2   (lt1)
      real             h2inv(lt1)

      real wp,xcg,pcg,rpcg
      common /scruz/   wp   (lt2)
     $ ,               xcg  (lt2)
     $ ,               pcg  (lt2) 
     $ ,               rpcg (lt2)
 
      real*8 etime1,dnekclock
      integer*8 ntotg,nxyz2

      integer intype,iter,iconv,ntot1,ntot2
      real alpha,beta,rrp1,rrp2,pap,div0,ratio
      real rnrm1,rrpx,rnorm,tolpss
      real h1_mx,h2_mx,wp_mx

      real glamax,glsc2
      real pcgmx


      etime1 = dnekclock()
      divex = 0.
      iter  = 0

      call chktcg2 (tolps,rcg,iconv)
      if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
     $   tolps = abs(param(21))

      nxyz2 = lx2*ly2*lz2
      ntot2 = nxyz2*nelv
      ntotg = nxyz2*nelgv

      call uzprec  (rpcg,rcg,h1,h2,intype,wp)
      rrp1 = glsc2 (rpcg,rcg,ntot2)
      call copy    (pcg,rpcg,ntot2)
      call rzero   (xcg,ntot2)
      if (rrp1.eq.0) return
      beta = 0.
      div0=0.

      tolpss = tolps
      do 1000 iter=1,nmxp

         call convprn (iconv,rnorm,rrp1,rcg,rpcg,tolpss)

         if (iter.eq.1)      div0   = rnorm
         if (param(21).lt.0) tolpss = abs(param(21))*div0

         ratio = rnorm/div0
         if (ifprint.and.nio.eq.0) 
     $   write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66    format(i5,1p4e12.5,i8,' Divergence')
c
         if (iconv.eq.1.and.iter.gt.1) goto 9000

         if (iter .ne. 1) then
            beta = rrp1/rrp2
            call add2s1 (pcg,rpcg,beta,ntot2)
         endif

         call cdabdtp_cyl(wp,pcg,h1,h2,h2inv,intype)
         pap   = glsc2 (pcg,wp,ntot2)

         if (pap.ne.0.) then
            alpha = rrp1/pap
         else
            pcgmx = glamax(pcg,ntot2)
            wp_mx = glamax(wp ,ntot2)
            ntot1 = lx1*ly1*lz1*nelv
            h1_mx = glamax(h1 ,ntot1)
            h2_mx = glamax(h2 ,ntot1)
            if (nid.eq.0) write(6,*) 'ERROR: pap=0 in uzawa.'
     $      ,iter,pcgmx,wp_mx,h1_mx,h2_mx
            call exitt
         endif
         call add2s2 (xcg,pcg,alpha,ntot2)
         call add2s2 (rcg,wp,-alpha,ntot2)

         if (iter.eq.-1) then
            call convprn (iconv,rnrm1,rrpx,rcg,rpcg,tolpss)
            if (iconv.eq.1) then
               rnorm = rnrm1
               ratio = rnrm1/div0
               if (nio.eq.0) 
     $         write (6,66) iter,tolpss,rnrm1,div0,ratio,istep
               goto 9000
            endif
         endif

         call ortho(rcg)

         rrp2 = rrp1
         call uzprec  (rpcg,rcg,h1,h2,intype,wp)

 1000 continue
      if (nid.eq.0) write (6,3001) iter,rnorm,tolpss
 3001 format(I6,' **ERROR**: Failed to converge in UZAWA:',6E13.4)
 9000 continue

      divex = rnorm
      iter  = iter-1

      if (iter.gt.0) call copy (rcg,xcg,ntot2)
      call ortho(rcg)

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep, '  U-Press std (CYL). ',
     &                            iter,divex,div0,tolpss,etime1
 9999 format(I11,a,I7,1p4E13.4)
19999 format(I11,'  U-Press 1.e-5: ',I7,1p4E13.4)


      return
      end subroutine uzawa_cyl
c-----------------------------------------------------------------------
      subroutine cggosf_cyl (u1,u2,u3,r1,r2,r3,h1,h2,rmult,binv,
     $                   vol,tin,maxit,matmod)

C     Conjugate gradient iteration for solution of coupled 
C     Helmholtz equations 

      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
      include 'FDMH1'

      common /screv/  dpc(lx1*ly1*lz1*lelt)
     $     ,          p1 (lx1*ly1*lz1*lelt)
      common /scrch/  p2 (lx1*ly1*lz1*lelt)
     $     ,          p3 (lx1*ly1*lz1*lelt)
      common /scrsl/  qq1(lx1*ly1*lz1*lelt)
     $     ,          qq2(lx1*ly1*lz1*lelt)
     $     ,          qq3(lx1*ly1*lz1*lelt)
      common /scrmg/  pp1(lx1*ly1*lz1*lelt)
     $     ,          pp2(lx1*ly1*lz1*lelt)
     $     ,          pp3(lx1*ly1*lz1*lelt)
     $     ,          wa (lx1*ly1*lz1*lelt)
      real ap1(1),ap2(1),ap3(1)
      equivalence (ap1,pp1),(ap2,pp2),(ap3,pp3)

      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      common /cprint/ ifprint
      logical ifdfrm, iffast, ifh2, ifsolv, ifprint

      real u1(1),u2(1),u3(1),
     $     r1(1),r2(1),r3(1),h1(1),h2(1),rmult(1),binv(1)


      logical iffdm,ifcrsl

      iffdm  = .true.
      iffdm  = .false.
c     ifcrsl = .true.
      ifcrsl = .false.

      nel   = nelfld(ifield)
      nxyz  = lx1*ly1*lz1
      n     = nxyz*nel

      if (istep.le.1.and.iffdm) call set_fdm_prec_h1A

      tol  = tin

c     overrule input tolerance
      if (restol(ifield).ne.0) tol=restol(ifield)

      if (ifcrsl) call set_up_h1_crs_strs(h1,h2,ifield,matmod)

      if ( .not.ifsolv ) then           !     Set logical flags
         call setfast (h1,h2,imesh)
         ifsolv = .true.
      endif

      call opdot (wa,r1,r2,r3,r1,r2,r3,n)
      rbnorm = glsc3(wa,binv,rmult,n)
      rbnorm = sqrt ( rbnorm / vol )
      if (rbnorm .lt. tol**2) then
         iter = 0
         r0 = rbnorm
c        if ( .not.ifprint )  goto 9999
         if (matmod.ge.0.and.nio.eq.0) write (6,3000) 
     $                                 istep,iter,rbnorm,r0,tol
         if (matmod.lt.0.and.nio.eq.0) write (6,3010) 
     $                                 istep,iter,rbnorm,r0,tol
         goto 9999
      endif

C     Evaluate diagional pre-conidtioner for fluid solve
      call setprec_cyl(qq1,h1,h2,imesh,1)
      call setprec_cyl(qq2,h1,h2,imesh,2)
      if (ldim.eq.3) then
         call setprec_cyl(qq3,h1,h2,imesh,3)
      endif

      if (iffdm) then
!         call set_fdm_prec_h1b(dpc,h1,h2,nel)
!         call fdm_h1a (pp1,r1,dpc,nel,ktype(1,1,1),wa)
!         call fdm_h1a (pp2,r2,dpc,nel,ktype(1,1,2),wa)
!         call fdm_h1a (pp3,r3,dpc,nel,ktype(1,1,3),wa)
!         call rmask   (pp1,pp2,pp3,nel)
!         call opdssum (pp1,pp2,pp3)
      else
        call col3 (pp1,qq1,r1,n)
        call col3 (pp2,qq2,r2,n)
        if (if3d) call col3 (pp3,qq3,r3,n)
      endif
      if (ifcrsl) then
!        call crs_strs(p1,p2,p3,r1,r2,r3)
!        call rmask   (p1,p2,p3,nel)
      else
        call opzero(p1,p2,p3)
      endif
      call opadd2       (p1,p2,p3,pp1,pp2,pp3)
      rpp1 = op_glsc2_wt(p1,p2,p3,r1,r2,r3,rmult)

      maxit=200
      do 1000 iter=1,maxit
!        call axhmsf  (ap1,ap2,ap3,p1,p2,p3,h1,h2,matmod)
        call axhmsf_cyl (ap1,ap2,ap3,p1,p2,p3,h1,h2,matmod)      ! prabal
        call rmask   (ap1,ap2,ap3,nel)
        call opdssum (ap1,ap2,ap3)
        pap   = op_glsc2_wt(p1,p2,p3,ap1,ap2,ap3,rmult)
        alpha = rpp1 / pap

        call opadds (u1,u2,u3,p1 ,p2 ,p3 , alpha,n,2)
        call opadds (r1,r2,r3,ap1,ap2,ap3,-alpha,n,2)

        call opdot  (wa,r1,r2,r3,r1,r2,r3,n)
        rbnorm = glsc3(wa,binv,rmult,n)
        rbnorm = sqrt (rbnorm/vol)

        if (iter.eq.1) r0 = rbnorm

        if (rbnorm.lt.tol) then
           ifin = iter
           if (nio.eq.0) then
              if (matmod.ge.0) write(6,3000) istep,ifin,rbnorm,r0,tol
              if (matmod.lt.0) write(6,3010) istep,ifin,rbnorm,r0,tol
           endif
           goto 9999
        endif

        if (iffdm) then
!           call fdm_h1a (pp1,r1,dpc,nel,ktype(1,1,1),wa)
!           call fdm_h1a (pp2,r2,dpc,nel,ktype(1,1,2),wa)
!           call fdm_h1a (pp3,r3,dpc,nel,ktype(1,1,3),wa)
!           call rmask   (pp1,pp2,pp3,nel)
!           call opdssum (pp1,pp2,pp3)
        else
           call col3 (pp1,qq1,r1,n)
           call col3 (pp2,qq2,r2,n)
           if (if3d) call col3 (pp3,qq3,r3,n)
        endif

        if (ifcrsl) then
!          call crs_strs(qq1,qq2,qq3,r1,r2,r3)
!          call rmask   (qq1,qq2,qq3,nel)
!          call opadd2  (pp1,pp2,pp3,qq1,qq2,qq3)
        endif

        call opdot (wa,r1,r2,r3,pp1,pp2,pp3,n)

        rpp2 = rpp1
        rpp1 = glsc2(wa,rmult,n)
        beta = rpp1/rpp2
        call opadds (p1,p2,p3,pp1,pp2,pp3,beta,n,1)

 1000 continue
      if (matmod.ge.0.and.nio.eq.0) write (6,3001) 
     $                              istep,iter,rbnorm,r0,tol
      if (matmod.lt.0.and.nio.eq.0) write (6,3011) 
     $                              istep,iter,rbnorm,r0,tol

 9999 continue
      ifsolv = .false.


 3000 format(i11,'  Helmh3 fluid (CYL) ',I6,1p3E13.4)
 3010 format(i11,'  Helmh3 mesh   ',I6,1p3E13.4)
 3001 format(i11,'  Helmh3 fluid unconverged! ',I6,1p3E13.4)
 3011 format(i11,'  Helmh3 mesh unconverged! ',I6,1p3E13.4)

      return
      end
c-----------------------------------------------------------------------
      subroutine setprec_cyl (dpcm1,helm1,helm2,imsh,isd)

!     Generate diagonal preconditioner for the Helmholtz operator.

      implicit none

      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'MASS'
      include 'CYLINDRICAL'

      real            dpcm1 (lx1,ly1,lz1,1)
      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      logical ifdfrm, iffast, ifh2, ifsolv
      real            helm1(lx1,ly1,lz1,1), helm2(lx1,ly1,lz1,1)
      real ysm1(ly1)

      integer ie,iq,iz,iy,ix,ntot
      integer i,j,e,iel,ijk
      integer nel,imsh,isd
      real term1,term2

      real const
      real r

      real rinv(lx1,ly1,lz1,lelv)
      real rinv2(lx1,ly1,lz1,lelv)

      nel=nelt
      if (imsh.eq.1) nel=nelv

      ntot = nel*lx1*ly1*lz1

      call rzero(dpcm1,ntot)
      do 1000 ie=1,nel

        if (ifaxis) call setaxdy ( ifrzer(ie) )

        const = 1.0
        if (isd.eq.1) const=2.0

        do 320 iq=1,lx1
        do 320 iz=1,lz1
        do 320 iy=1,ly1
        do 320 ix=1,lx1
           r                  = cyl_radius(ix,iy,iz,ie)
           dpcm1(ix,iy,iz,ie) = dpcm1(ix,iy,iz,ie) + 
     $                    r*const*g1m1(iq,iy,iz,ie) * dxtm1(ix,iq)**2
  320   continue

        const = 1.0
        if (isd.eq.2) const=2.0
        do 340 iq=1,ly1
        do 340 iz=1,lz1
        do 340 iy=1,ly1
        do 340 ix=1,lx1
           r                  = cyl_radius(ix,iy,iz,ie)
           dpcm1(ix,iy,iz,ie) = dpcm1(ix,iy,iz,ie) + 
     $                     r*const*g2m1(ix,iq,iz,ie) * dytm1(iy,iq)**2
  340   continue
        if (ldim.eq.3) then
           const = 1.0
           if (isd.eq.2) const=2.0
         
           do 360 iq=1,lz1
           do 360 iz=1,lz1
           do 360 iy=1,ly1
           do 360 ix=1,lx1
              r                  = cyl_radius(ix,iy,iz,ie)
              dpcm1(ix,iy,iz,ie) = dpcm1(ix,iy,iz,ie) + const* 
     $     g3m1(ix,iy,iq,ie)*(dztm1(iz,iq)**2)/r
  360      continue
c
c          add cross terms if element is deformed.
c
!          I am not caring about this since we don't have deformed
!          elements in the cylindrical formulation (yet)
           if (ifdfrm(ie)) then
              do 600 iy=1,ly1,ly1-1
              do 600 iz=1,lz1,max(1,lz1-1)
              dpcm1(1,iy,iz,ie) = dpcm1(1,iy,iz,ie)
     $            + g4m1(1,iy,iz,ie) * dxtm1(1,1)*dytm1(iy,iy)
     $            + g5m1(1,iy,iz,ie) * dxtm1(1,1)*dztm1(iz,iz)
              dpcm1(lx1,iy,iz,ie) = dpcm1(lx1,iy,iz,ie)
     $            + g4m1(lx1,iy,iz,ie) * dxtm1(lx1,lx1)*dytm1(iy,iy)
     $            + g5m1(lx1,iy,iz,ie) * dxtm1(lx1,lx1)*dztm1(iz,iz)
  600         continue
              do 700 ix=1,lx1,lx1-1
              do 700 iz=1,lz1,max(1,lz1-1)
                 dpcm1(ix,1,iz,ie) = dpcm1(ix,1,iz,ie)
     $            + g4m1(ix,1,iz,ie) * dytm1(1,1)*dxtm1(ix,ix)
     $            + g6m1(ix,1,iz,ie) * dytm1(1,1)*dztm1(iz,iz)
                 dpcm1(ix,ly1,iz,ie) = dpcm1(ix,ly1,iz,ie)
     $            + g4m1(ix,ly1,iz,ie) * dytm1(ly1,ly1)*dxtm1(ix,ix)
     $            + g6m1(ix,ly1,iz,ie) * dytm1(ly1,ly1)*dztm1(iz,iz)
  700         continue
              do 800 ix=1,lx1,lx1-1
              do 800 iy=1,ly1,ly1-1
                 dpcm1(ix,iy,1,ie) = dpcm1(ix,iy,1,ie)
     $                + g5m1(ix,iy,1,ie) * dztm1(1,1)*dxtm1(ix,ix)
     $                + g6m1(ix,iy,1,ie) * dztm1(1,1)*dytm1(iy,iy)
                 dpcm1(ix,iy,lz1,ie) = dpcm1(ix,iy,lz1,ie)
     $                + g5m1(ix,iy,lz1,ie) * dztm1(lz1,lz1)*dxtm1(ix,ix)
     $                + g6m1(ix,iy,lz1,ie) * dztm1(lz1,lz1)*dytm1(iy,iy)
  800         continue
           endif

        else  ! 2d

           iz=1
           if (ifdfrm(ie)) then
              do 602 iy=1,ly1,ly1-1
                 dpcm1(1,iy,iz,ie) = dpcm1(1,iy,iz,ie)
     $                + g4m1(1,iy,iz,ie) * dxtm1(1,1)*dytm1(iy,iy)
                 dpcm1(lx1,iy,iz,ie) = dpcm1(lx1,iy,iz,ie)
     $                + g4m1(lx1,iy,iz,ie) * dxtm1(lx1,lx1)*dytm1(iy,iy)
  602         continue
              do 702 ix=1,lx1,lx1-1
                 dpcm1(ix,1,iz,ie) = dpcm1(ix,1,iz,ie)
     $                + g4m1(ix,1,iz,ie) * dytm1(1,1)*dxtm1(ix,ix)
                 dpcm1(ix,ly1,iz,ie) = dpcm1(ix,ly1,iz,ie)
     $                + g4m1(ix,ly1,iz,ie) * dytm1(ly1,ly1)*dxtm1(ix,ix)
  702         continue
           endif

        endif
 1000 continue

!     Add additional terms for cylindrical formulation

      call invers2(rinv,cyl_radius,ntot)        ! rinv = 1/R
      call copy(rinv2,rinv,ntot)
      call invcol2(rinv2,rinv,ntot)             ! rinv2 = 1/R^2

      if (isd.eq.1) then
!        call cmult(rinv2,const,ntot)            ! (k^2)/R^2
!        call Xaddcol3(dpcm1,rinv2,bm1,ntot)     ! D = D + BM1*(k^2)/R^2

      elseif (isd.eq.2) then
        const = 2.0
        call cmult(rinv2,const,ntot)            ! (2/R^2)
        call Xaddcol3(dpcm1,rinv2,bm1,ntot)     ! D = D + BM1*(2/R^2)
      elseif (isd.eq.3) then
!       Again, ignoring cross derivatives for now        
        do 380 ie=1,nelv
        do 380 iq=1,ly1
        do 380 iz=1,lz1
        do 380 iy=1,ly1
        do 380 ix=1,lx1
           dpcm1(ix,iy,iz,ie) = dpcm1(ix,iy,iz,ie)
     $     - sym1(ix,iq,iz,ie)*w3m1(ix,iy,iz)*dym1(iy,iq)
     $     - sym1(ix,iq,iz,ie)*w3m1(ix,iy,iz)*dytm1(iy,iq)
     $     + bm1(ix,iy,iz,ie)*rinv2(ix,iy,iz,ie)             
  380   continue
      else

        if (nid.eq.0) then
          write(6,*) 'Unrecognized ISD in setprec_cyl', ISD
        endif
        call exitt
      endif  

      call col2    (dpcm1,helm1,ntot)
      call addcol3 (dpcm1,helm2,bm1,ntot)

      call dssum (dpcm1,lx1,ly1,lz1)
      call invcol1 (dpcm1,ntot)

      return
      end
!---------------------------------------------------------------------- 




