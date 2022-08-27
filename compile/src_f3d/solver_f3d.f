!======================================================================
!     Solver Routines to take care of complex variables 
!     in Cylindrical Coordinates.
!     Author: Prabal S. Negi
!
!====================================================================== 
!-----------------------------------------------------------------------

      subroutine ophinv_cyl(or1,or2,or3,oi1,oi2,oi3,
     $                      ir1,ir2,ir3,ii1,ii2,ii3,
     $                      h1,h2,tolh,nmxhi)

!     Ok = (H1*A+H2*B)-1 * Ik  (implicit)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'VPROJ'
      include 'TSTEP'

!     Real outputs
      real or1(lx1,ly1,lz1,1),or2(lx1,ly1,lz1,1),or3(lx1,ly1,lz1,1)
!     Imaginary outputs      
      real oi1(lx1,ly1,lz1,1),oi2(lx1,ly1,lz1,1),oi3(lx1,ly1,lz1,1)
     
!     Real Inputs      
      real ir1(lx1,ly1,lz1,1),ir2(lx1,ly1,lz1,1),ir3(lx1,ly1,lz1,1)
!     Imaginary Outputs      
      real ii1(lx1,ly1,lz1,1),ii2(lx1,ly1,lz1,1),ii3(lx1,ly1,lz1,1)
     
      real h1 (lx1,ly1,lz1,1) , h2 (lx1,ly1,lz1,1)

      integer mtmp,i,matmod,nmxhi
      real tolh
 
c      ifproj = .false.
c      if (param(94).gt.0)    ifproj = .true.
c      if (ifprojfld(ifield)) ifproj = .true.
c 
c      if (.not.ifproj) then
c         if (ifield.eq.1) call ophinv
c     $      (o1,o2,o3,i1,i2,i3,h1,h2,tolh,nmxhi)
c         if (ifield.eq.ifldmhd) call ophinv
c     $      (o1,o2,o3,i1,i2,i3,h1,h2,tolh,nmxhi)
c         return
c      endif
 
      mtmp = param(93)
      do i=1,2*ldim
         ivproj(1,i) = min(mxprev,mtmp) - 1
      enddo
 
      imesh = 1
 
      if (ifstrs) then
         matmod = 0
!        'NOMG' = No multigrid
         call hmhzsf_cyl('NOMG',or1,or2,or3,oi1,oi2,oi3,
     $                  ir1,ir2,ir3,ii1,ii2,ii3,h1,h2,
     $                  v1mask,v2mask,v3mask,vmult,
     $                  tolh,nmxhi,matmod)
      else
         if (ifield.eq.1) then
           if (nid.eq.0) then 
             write(6,*) 'IFSTRS = ', ifstrs
             write(6,*) 'Cylindrical Solver needs coupled solve'
             write(6,*) 'Set IFSTRS = TRUE'
           endif  

           call exitt

         elseif (ifield.eq.ifldmhd) then  ! B-field

           if (nid.eq.0) 
     $       write(6,*) 'Cylindrical Solver not implemented for MHD'

           call exitt
         endif
      endif
C
      return
      end subroutine ophinv_cyl 
c--------------------------------------------------------------------

      subroutine hmhzsf_cyl(name,ur1,ur2,ur3,ui1,ui2,ui3,
     $                   rr1,rr2,rr3,ri1,ri2,ri3,h1,h2,
     $                   rmask1,rmask2,rmask3,rmult,
     $                   tol,maxit,matmod)

!     Solve coupled Helmholtz equations (stress formulation)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'   ! For outpost diagnostic call
      include 'TSTEP'
      include 'ORTHOSTRS'
      include 'CTIMER'

      real ur1(1),ur2(1),ur3(1),ui1(1),ui2(1),ui3(1)
      real rr1(1),rr2(1),rr3(1),ri1(1),ri2(1),ri3(1)
      real rmask1(1),rmask2(1),rmask3(1),rmult(1)
      real h1(1),h2(1)
      character name*4

      integer iproj,nel,n,maxit,matmod
      real tol
      real vol

#ifdef TIMER
      nhmhz = nhmhz + 1
      etime1 = dnekclock()
#endif

      nel = nelfld(ifield)
      vol = volfld(ifield)
      n   = lx1*ly1*lz1*nel

      napproxstrs(1) = 0
      iproj = 0
      if (ifprojfld(ifield)) iproj = param(94)
      if (iproj.gt.0.and.istep.ge.iproj) napproxstrs(1)=param(93)
      napproxstrs(1)=min(napproxstrs(1),mxprev)

!     Just apply the mask.
!     More complicated for Symmetry/On boundary conditions
!     Which is not implemented yet.      
!      call rmask   (rr1,rr2,rr3,nel)
!      call rmask   (ri1,ri2,ri3,nel)
      call col2_3(rr1,rr2,rr3,v1mask,v2mask,v3mask,n)
      call col2_3(ri1,ri2,ri3,v1mask,v2mask,v3mask,n)
    
      call dssum3 (rr1,rr2,rr3)
      call dssum3 (ri1,ri2,ri3)

      call rzero3  (ur1,ur2,ur3,n)
      call rzero3  (ui1,ui2,ui3,n)

      if (imesh.eq.1) then
!         call chktcgs (r1,r2,r3,rmask1,rmask2,rmask3,rmult,binvm1
!     $                ,vol,tol,nel)

!         call strs_project_a(r1,r2,r3,h1,h2,rmult,ifield,ierr,matmod)

!         call cggosf  (u1,u2,u3,r1,r2,r3,h1,h2,rmult,binvm1
!     $                ,vol,tol,maxit,matmod)

         call cggosf_cyl(ur1,ur2,ur3,ui1,ui2,ui3,
     $                   rr1,rr2,rr3,ri1,ri2,ri3,h1,h2,
     $                   rmult,binvm1,vol,tol,maxit,matmod)

!         call strs_project_b(u1,u2,u3,h1,h2,rmult,ifield,ierr)

      else

!         call chktcgs (r1,r2,r3,rmask1,rmask2,rmask3,rmult,bintm1
!     $                ,vol,tol,nel)
!         call cggosf  (u1,u2,u3,r1,r2,r3,h1,h2,rmult,bintm1
!     $                ,vol,tol,maxit,matmod)

      endif

#ifdef TIMER
      thmhz=thmhz+(dnekclock()-etime1)
#endif

      return
      end subroutine hmhzsf_cyl
!----------------------------------------------------------------------       

      subroutine cggosf_cyl(u1r,u2r,u3r,u1i,u2i,u3i,
     $                      r1r,r2r,r3r,r1i,r2i,r3i,h1,h2,
     $                      rmult,binv,vol,tin,maxit,matmod)

!     Conjugate gradient iteration for solution of coupled 
!     Helmholtz equations, with complex entires
!     Actual entiries treated with real variables 

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'          ! V?MASK
!      include 'TOTAL'
      include 'MASS'
      include 'DOMAIN'
      include 'FDMH1'

      include 'F3D'           ! k_f3d

      include 'TEST'

!     r?r,r?i           - Residuals
!     p?r,p?i           - Search directions
!     Ap?r, Ap?i        - A*p
!     qq?               - Diagonal Preconditioner      

      real dpc,p1r,p2r,p3r
      common /screv11/  dpc(lx1*ly1*lz1*lelt)
     $     ,          p1r (lx1*ly1*lz1*lelt)
      common /scrch11/  p2r (lx1*ly1*lz1*lelt)
     $     ,          p3r (lx1*ly1*lz1*lelt)

!     This is not how its done in the original routines
!     But I'm going to use qq1,qq22,qq3 as arrays for diagonal
!     preconditioning.
!     Originally it has been used for coarse grid solves.
!     But we don't do any coarse grid solve.      
      real qq1,qq2,qq3
      common /scrsl11/  qq1(lx1*ly1*lz1*lelt)
     $     ,          qq2(lx1*ly1*lz1*lelt)
     $     ,          qq3(lx1*ly1*lz1*lelt)

      real Ap1r,Ap2r,Ap3r,wa
      common /scrmg11/  Ap1r(lx1*ly1*lz1*lelt)
     $     ,          Ap2r(lx1*ly1*lz1*lelt)
     $     ,          Ap3r(lx1*ly1*lz1*lelt)
     $     ,          wa (lx1*ly1*lz1*lelt)

      real p1i,p2i,p3i,Ap1i,Ap2i,Ap3i
      common /scruz11/ p1i(lx1*ly1*lz1*lelt),
     $               p2i(lx1*ly1*lz1*lelt),
     $               p3i(lx1*ly1*lz1*lelt),
     $               Ap1i(lx1*ly1*lz1*lelt)
      common /scrvh11/ Ap2i(lx1*ly1*lz1*lelt),
     $               Ap3i(lx1*ly1*lz1*lelt)



!      real Ap1r(1),Ap2r(1),Ap3r(1)
!      real Ap1i(1),Ap2i(1),Ap3i(1)
!     Why do we do this equivalence?
!      equivalence (Ap1r,pp1r),(Ap2r,pp2r),(Ap3r,pp3r)
!      equivalence (Ap1i,pp1i),(Ap2i,pp2i),(Ap3i,pp3i)

      logical ifdfrm, iffast, ifh2, ifsolv, ifprint
      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      common /cprint/ ifprint

      real u1r(1),u2r(1),u3r(1),u1i(1),u2i(1),u3i(1),
     $     r1r(1),r2r(1),r3r(1),r1i(1),r2i(1),r3i(1),
     $     h1(1),h2(1),rmult(1),binv(1)

      integer maxit,matmod,nel,nxyz,n,iter,ifin
      logical iffdm,ifcrsl

      real tin,tol,vol
      real alphar,alphai,betar,betai,papr,papi
      real rpp1r,rpp1i,rpp2r,rpp2i

      real tmpval

      real r0,rbnorm

      real glsc2,glsc3,op_glsc2_wt
      real opnorm2_wt_comp

!     This common block is used in weak laplacian calculations
!     In this routine this is only used as a work array
!     Since the entries are not required after the weak laplacian
!     has been evaluated      
      real wk1r,wk2r,wk3r,wk1i,wk2i,wk3i,wk4
      common /scrns/   wk1r(lx1*ly1*lz1*lelt)
     $               , wk2r(lx1*ly1*lz1*lelt)
     $               , wk3r(lx1*ly1*lz1*lelt)
     $               , wk1i(lx1*ly1*lz1*lelt)
     $               , wk2i(lx1*ly1*lz1*lelt)
     $               , wk3i(lx1*ly1*lz1*lelt)
     $               , wk4(lx1*ly1*lz1*lelt) 

!      integer lkryl
!      parameter (lkryl = 5)
!      real krylr(lx1*ly1*lz1*lelv,3,lkryl)
!      real kryli(lx1*ly1*lz1*lelv,3,lkryl)
!      real orthor(lkryl,lkryl),orthoi(lkryl,lkryl)
!      complex*16 orthoc(lkryl,lkryl)

      logical ifjacobi        ! Apply Jacobi preconditioner?

      integer n1,n2,i,j

!     No Fast Diagonalization Method            
      iffdm  = .false.
!     Jacobi preconditioner 
      ifjacobi = .true.
!     No Coarse grid
      ifcrsl = .false.

      nel   = nelfld(ifield)
      nxyz  = lx1*ly1*lz1
      n     = nxyz*nel

      if (istep.le.1.and.iffdm) call set_fdm_prec_h1A

      tol  = tin

!     overrule input tolerance
      if (restol(ifield).ne.0) tol=restol(ifield)
      if (ifcrsl) call set_up_h1_crs_strs(h1,h2,ifield,matmod)

      if ( .not.ifsolv ) then           !     Set logical flags
         ifaxis = .true.
         call setfast (h1,h2,imesh)
         ifsolv = .true.
         ifaxis = .false.
      endif

!     call opdot (wa,r1r,r2r,r3r,r1r,r2r,r3r,n)
!     rbnorm = glsc3(wa,binv,rmult,n)
!     rbnorm = sqrt ( rbnorm / vol )
      call col3(wa,binv,rmult,n)
      rbnorm  = opnorm2_wt_comp(r1r,r2r,r3r,r1i,r2i,r3i,wa,n)
      rbnorm = sqrt(rbnorm/vol) 

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

!     prabal. Not sure why preconditioners were added together for
!     different components.        
!     Evaluate diagional pre-conidtioner for fluid solve
      if (ifjacobi) then
        call setprec_cyl (qq1,h1,h2,k_f3d,imesh,1)
        call setprec_cyl (qq2,h1,h2,k_f3d,imesh,2)
        call setprec_cyl (qq3,h1,h2,k_f3d,imesh,3)
      else
        call rone(qq1,n)
        call rone(qq2,n)
        call rone(qq3,n)
      endif  

      if (iffdm) then
!         call set_fdm_prec_h1b(dpc,h1,h2,nel)
!         call fdm_h1a (pp1,r1,dpc,nel,ktype(1,1,1),wa)
!         call fdm_h1a (pp2,r2,dpc,nel,ktype(1,1,2),wa)
!         call fdm_h1a (pp3,r3,dpc,nel,ktype(1,1,3),wa)
!         call rmask   (pp1,pp2,pp3,nel)
!         call opdssum (pp1,pp2,pp3)
      elseif (ifjacobi) then
!        Apply Diagonal preconditioners                     
!        Real
         call col3 (Ap1r,qq1,r1r,n)
         call col3 (Ap2r,qq2,r2r,n)
         call col3 (Ap3r,qq3,r3r,n)

!        Imaginary         
         call col3 (Ap1i,qq1,r1i,n)
         call col3 (Ap2i,qq2,r2i,n)
         call col3 (Ap3i,qq3,r3i,n)
      else
         call copy3(Ap1r,Ap2r,Ap3r,r1r,r2r,r3r,n) 
         call copy3(Ap1i,Ap2i,Ap3i,r1i,r2i,r3i,n) 
      endif

      if (ifcrsl) then
!         call crs_strs(p1,p2,p3,r1,r2,r3)
!         call rmask   (p1,p2,p3,nel)
      else
         call rzero3(p1r,p2r,p3r,n)
         call rzero3(p1i,p2i,p3i,n)
      endif

      call add2_3(p1r,p2r,p3r,Ap1r,Ap2r,Ap3r,n)
      call add2_3(p1i,p2i,p3i,Ap1i,Ap2i,Ap3i,n)

!      rpp1 = r*rmult*(D*r)
      call opglsc2_wt_comp(rpp1r,rpp1i,r1r,r2r,r3r,r1i,r2i,r3i,
     $                     Ap1r,Ap2r,Ap3r,Ap1i,Ap2i,Ap3i,rmult,n)

!     Imaginary part must be zero.
!     Otherwise the Diagonal is not purely real, which means 
!     the problem is not symmetric and PCG will not work.
!     Should probably print out and verify
!     Setting it zero so we don't accumulate round off errors
      rpp1i = 0.      

      maxit=500
      do 1000 iter=1,maxit
         call axhmsf_cyl(Ap1r,Ap2r,Ap3r,Ap1i,Ap2i,Ap3i,
     $                      p1r,p2r,p3r,p1i,p2i,p3i,h1,h2)


!!     prabal
!        call col3(Ap1r,p1r,bm1,n)
!        call col3(Ap2r,p2r,bm1,n)
!        call col3(Ap3r,p3r,bm1,n)
!        call col3(Ap1i,p1i,bm1,n)
!        call col3(Ap2i,p2i,bm1,n)
!        call col3(Ap3i,p3i,bm1,n)

!         call rmask(ap1r,ap2r,ap3r,nel)
!         call rmask(ap1i,ap2i,ap3i,nel)
         call col2_3(Ap1r,Ap2r,Ap3r,v1mask,v2mask,v3mask,n)
         call col2_3(Ap1i,Ap2i,Ap3i,v1mask,v2mask,v3mask,n)

         call dssum3(Ap1r,Ap2r,Ap3r)
         call dssum3(Ap1i,Ap2i,Ap3i)

!        pAp = pn-1*Apn-1            
         call opglsc2_wt_comp(pApr,pApi,p1r,p2r,p3r,p1i,p2i,p3i,
     $                        Ap1r,Ap2r,Ap3r,Ap1i,Ap2i,Ap3i,rmult,n)


!        \alpha = (rn-1*rn-1)/(pn-1*Apn-1) 
         tmpval = (pApr*pApr + pApi*pApi)
         alphar = (rpp1r*pApr/tmpval + rpp1i*pApi/tmpval)
         alphai = (rpp1i*pApr/tmpval - rpp1r*pApi/tmpval)

!        Update Soln: u = u + \alpha*p
!        ur = ur + \alphar*pr - \alphai*pi
         call opadds_3(u1r,u2r,u3r,p1r,p2r,p3r,alphar,n,2)
         call opadds_3(u1r,u2r,u3r,p1r,p2r,p3r,-alphai,n,2)

!        ui = ui + \alphar*pi + \alphai*pr
         call opadds_3(u1i,u2i,u3i,p1i,p2i,p3i,alphar,n,2)
         call opadds_3(u1i,u2i,u3i,p1r,p2r,p3r,alphai,n,2)


!        Update Residual: r = r - \alpha*p 
!        rr = rr - \alphar*Apr + \alphai*Api
         call opadds_3(r1r,r2r,r3r,Ap1r,Ap2r,Ap3r,-alphar,n,2)
         call opadds_3(r1r,r2r,r3r,Ap1i,Ap2i,Ap3i,alphai,n,2)

!        ri = ri - \alphar*Api - \alphai*Apr
         call opadds_3(r1i,r2i,r3i,Ap1i,Ap2i,Ap3i,-alphar,n,2)
         call opadds_3(r1i,r2i,r3i,Ap1r,Ap2r,Ap3r,-alphai,n,2)

!        rbnorm = ||r||                
         call col3(wa,binv,rmult,n)
         rbnorm = opnorm2_wt_comp(r1r,r2r,r3r,r1i,r2i,r3i,wa,n)
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
         elseif (ifjacobi) then
!          Apply Diagonal preconditioners                     
!          Real           
           call col3 (Ap1r,qq1,r1r,n)
           call col3 (Ap2r,qq2,r2r,n)
           call col3 (Ap3r,qq3,r3r,n)

!          Imaginary         
           call col3 (Ap1i,qq1,r1i,n)
           call col3 (Ap2i,qq2,r2i,n)
           call col3 (Ap3i,qq3,r3i,n)
         else
           call copy3(Ap1r,Ap2r,Ap3r,r1r,r2r,r3r,n) 
           call copy3(Ap1i,Ap2i,Ap3i,r1i,r2i,r3i,n) 
         endif

         if (ifcrsl) then
!           call crs_strs(qq1,qq2,qq3,r1,r2,r3)
!           call rmask   (qq1,qq2,qq3,nel)
!           call opadd2  (pp1,pp2,pp3,qq1,qq2,qq3)
         endif

!        Old Residual
         rpp2r = rpp1r
         rpp2i = rpp1i

!        New Residual         
!        rpp1 = r*rmult*(D*r)
!        Here Ap == D*r            
         call opglsc2_wt_comp(rpp1r,rpp1i,r1r,r2r,r3r,r1i,r2i,r3i,
     $                        Ap1r,Ap2r,Ap3r,Ap1i,Ap2i,Ap3i,rmult,n)

!        Again, imaginary part must be zero.
!        Otherwise the Diagonal is not purely real, which means 
!        the problem is not symmetric and PCG will not work.
!        Setting it zero so we don't accumulate round off errors
         rpp1i = 0. 

         tmpval = (rpp2r*rpp2r + rpp2i*rpp2i)

         betar = rpp1r*rpp2r/tmpval + rpp1i*rpp2i/tmpval
         betai = rpp1i*rpp2r/tmpval - rpp1r*rpp2i/tmpval

!        Again, betai must also be exactly zero.
!        Setting it to zero
         betai = 0.         

         call copy3(wk1r,wk2r,wk3r,p1r,p2r,p3r,n)
         call copy3(wk1i,wk2i,wk3i,p1i,p2i,p3i,n)

!        pr = D*rr + \betar*pr - \betai*pi
         call copy3(p1r,p2r,p3r,Ap1r,Ap2r,Ap3r,n)
         call opadds_3(p1r,p2r,p3r,wk1r,wk2r,wk3r,betar,n,2)   
!         call opadds_3(p1r,p2r,p3r,wk1i,wk2i,wk3i,-betai,n,2)   

!        pi = D*ri + \betar*pi + \betai*pr
         call copy3(p1i,p2i,p3i,Ap1i,Ap2i,Ap3i,n)
         call opadds_3(p1i,p2i,p3i,wk1i,wk2i,wk3i,betar,n,2)   



 1000 continue
      if (matmod.ge.0.and.nio.eq.0) write (6,3001) 
     $                              istep,iter,rbnorm,r0,tol
      if (matmod.lt.0.and.nio.eq.0) write (6,3011) 
     $                              istep,iter,rbnorm,r0,tol

 9999 continue
      ifsolv = .false.


 3000 format(i11,'  Helmh3 fluid  ',I6,1p3E13.4)
 3010 format(i11,'  Helmh3 mesh   ',I6,1p3E13.4)
 3001 format(i11,'  Helmh3 fluid unconverged! ',I6,1p3E13.4)
 3011 format(i11,'  Helmh3 mesh unconverged! ',I6,1p3E13.4)


      return
      end subroutine cggosf_cyl
!-----------------------------------------------------------------------

      real function opnorm2_wt_comp(u1r,u2r,u3r,u1i,u2i,u3i,wt,n)

!                _       
!     Calculate: u *Wt* u
!     Assuming the Weight itself is real        

      implicit none

!      real opnorm2_wt_comp

      integer i,n

      real u1r(n),u2r(n),u3r(n),u1i(n),u2i(n),u3i(n)
      real wt(n)

      real a1,a2,wk
      real glsum

      a1 = 0.
      a2 = 0.
      do i=1,n
        a1 = a1 + wt(i)*(u1r(i)*u1r(i) + u2r(i)*u2r(i) + u3r(i)*u3r(i))
        a2 = a2 + wt(i)*(u1i(i)*u1i(i) + u2i(i)*u2i(i) + u3i(i)*u3i(i))
      enddo
     
      wk = a1 + a2

      opnorm2_wt_comp = glsum(wk,1)

!      opnorm2_wt_comp = wk

      return
      end function opnorm2_wt_comp

!-----------------------------------------------------------------------

      subroutine opglsc2_wt_comp(scr,sci,u1r,u2r,u3r,u1i,u2i,u3i,
     $                              v1r,v2r,v3r,v1i,v2i,v3i,wt,n)

!                _       
!     Calculate: u *Wt* v
!     Assuming the Weight itself is real        

      implicit none

      integer i,n

      real u1r(n),u2r(n),u3r(n),u1i(n),u2i(n),u3i(n)
      real v1r(n),v2r(n),v3r(n),v1i(n),v2i(n),v3i(n)
      real wt(n)

      real scr,sci,a1,a2
      real glsum

      a1 = 0.
      a2 = 0.
      do i=1,n
        a1 = a1 + wt(i)*(u1r(i)*v1r(i) + u2r(i)*v2r(i) + u3r(i)*v3r(i))
        a1 = a1 + wt(i)*(u1i(i)*v1i(i) + u2i(i)*v2i(i) + u3i(i)*v3i(i))

        a2 = a2 + wt(i)*(u1r(i)*v1i(i) + u2r(i)*v2i(i) + u3r(i)*v3i(i))
        a2 = a2 - wt(i)*(u1i(i)*v1r(i) + u2i(i)*v2r(i) + u3i(i)*v3r(i))
      enddo

!     a1 = ur*Wt*vr + ui*Wt*vi
!     a2 = ur*Wt*vi - ui*Wt*vr
!     sc = a1 + i*a2      

      scr = a1
      sci = a2

!     Sum over all processors 
!      call glsum(scr,1)
!      call glsum(sci,1)
      scr = glsum(a1,1)
      sci = glsum(a2,1)


      return
      end subroutine opglsc2_wt_comp

!-----------------------------------------------------------------------

      subroutine copy3(a1,a2,a3,b1,b2,b3,n)

      implicit none  

      integer i,n

      real a1(1),a2(1),a3(1)
      real b1(1),b2(1),b3(1)

      do i=1,n
         a1(i)=b1(i)
         a2(i)=b2(i)
         a3(i)=b3(i)
      enddo


      return
      end subroutine copy3
!-----------------------------------------------------------------------

      subroutine dssum3(a1,a2,a3)

      implicit none

      include 'SIZE'

      real a1(1),a2(1),a3(1)

      call dssum(a1,nx1,ny1,nz1)
      call dssum(a2,nx1,ny1,nz1)      
      call dssum(a3,nx1,ny1,nz1)

      return  
      end subroutine dssum3        
!-----------------------------------------------------------------------

      subroutine opadds_3 (a1,a2,a3,b1,b2,b3,const,n,isc)

      implicit none

      integer n,isc
      real const

      real a1(1),a2(1),a3(1),b1(1),b2(1),b3(1)

      if (isc.eq.1) then
         call add2s1 (a1,b1,const,n)
         call add2s1 (a2,b2,const,n)
         call add2s1 (a3,b3,const,n)
      elseif (isc.eq.2) then
         call add2s2 (a1,b1,const,n)
         call add2s2 (a2,b2,const,n)
         call add2s2 (a3,b3,const,n)
      endif

      return
      end subroutine opadds_3
!-----------------------------------------------------------------------
      subroutine add2_3 (a1,a2,a3,b1,b2,b3,n)

      implicit none

      integer n

      real a1(1),a2(1),a3(1),b1(1),b2(1),b3(1)

      call add2(a1,b1,n)
      call add2(a2,b2,n)
      call add2(a3,b3,n)

      return
      end subroutine add2_3
!-----------------------------------------------------------------------
      subroutine col2_3 (a1,a2,a3,b1,b2,b3,n)

      implicit none

      integer n

      real a1(1),a2(1),a3(1),b1(1),b2(1),b3(1)

      call col2(a1,b1,n)
      call col2(a2,b2,n)
      call col2(a3,b3,n)

      return
      end subroutine col2_3
!-----------------------------------------------------------------------
      subroutine setprec_cyl (dpcm1,helm1,helm2,k_f3d,imsh,isd)

!     Generate diagonal preconditioner for the Helmholtz operator.

      implicit none

      include 'SIZE'
      include 'WZ'
      include 'DXYZ'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'MASS'

      real            dpcm1 (lx1,ly1,lz1,1)
      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      logical ifdfrm, iffast, ifh2, ifsolv
      real            helm1(lx1,ly1,lz1,1), helm2(lx1,ly1,lz1,1)
      real ysm1(ly1)

      integer ie,iq,iz,iy,ix,ntot
      integer i,j,e,iel,ijk
      integer nel,imsh,isd
      real term1,term2

      real k_f3d                ! passed as argument
      real const

      real rinv(lx1,ly1,lz1,lelv)
      real rinv2(lx1,ly1,lz1,lelv)

      nel=nelt
      if (imsh.eq.1) nel=nelv

      ntot = nel*lx1*ly1*lz1

      call rzero(dpcm1,ntot)
      do 1000 ie=1,nel

        if (ifaxis) call setaxdy ( ifrzer(ie) )

        do 320 iq=1,lx1
        do 320 iz=1,lz1
        do 320 iy=1,ly1
        do 320 ix=1,lx1
           dpcm1(ix,iy,iz,ie) = dpcm1(ix,iy,iz,ie) + 
     $                          g1m1(iq,iy,iz,ie) * dxtm1(ix,iq)**2
  320   continue
        do 340 iq=1,ly1
        do 340 iz=1,lz1
        do 340 iy=1,ly1
        do 340 ix=1,lx1
           dpcm1(ix,iy,iz,ie) = dpcm1(ix,iy,iz,ie) + 
     $                          g2m1(ix,iq,iz,ie) * dytm1(iy,iq)**2
  340   continue
        if (ldim.eq.3) then
           do 360 iq=1,lz1
           do 360 iz=1,lz1
           do 360 iy=1,ly1
           do 360 ix=1,lx1
              dpcm1(ix,iy,iz,ie) = dpcm1(ix,iy,iz,ie) + 
     $                             g3m1(ix,iy,iq,ie) * dztm1(iz,iq)**2
  360      continue
c
c          add cross terms if element is deformed.
c
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

!     Add additional terms for ifcyl and Fourier

      call invers2(rinv,ym1,ntot)               ! rinv = 1/R
      call copy(rinv2,rinv,ntot)
      call invcol2(rinv2,rinv,ntot)             ! rinv2 = 1/R^2
      const = k_f3d*k_f3d                       ! k^2

      if (isd.eq.1) then
        call cmult(rinv2,const,ntot)            ! (k^2)/R^2
        call Xaddcol3(dpcm1,rinv2,bm1,ntot)     ! D = D + BM1*(k^2)/R^2

      elseif (isd.eq.2) then
        call cmult(rinv2,const,ntot)            ! (k^2)/R^2
        call Xaddcol3(dpcm1,rinv2,bm1,ntot)     ! D = D + BM1*(k^2)/R^2

      elseif (isd.eq.3) then
        const = 2.0*const      
        call cmult(rinv2,const,ntot)            ! (2*k^2)/R^2
        call Xaddcol3(dpcm1,rinv2,bm1,ntot)     ! D = D + BM1*(2*k^2)/R^2

!       Second term
        call col2(rinv,bm1,ntot)                ! BM1/R

!       Using rinv2 as a temporary array        
!       dvi(ri)/dr*dr/dy + dvj(sj)/ds*ds/dy
        do e = 1,nel
        do j = 1,ly1
        do i = 1,lx1
          ijk = i + (j-1)*lx1
          rinv2(i,j,1,e) = dxm1(i,i)*jacmi(ijk,e)*rym1(i,j,1,e) +
     $                     dym1(j,j)*jacmi(ijk,e)*sym1(i,j,1,e)

        enddo
        enddo
        enddo
        call col2(rinv2,rinv,ntot)              ! BM1/R*dv/dy
        call sub2(dpcm1,rinv2,ntot)             ! D = D + BM1/R*dv/dy
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

      subroutine uzawa_f3d(rcgr,rcgi,h1,h2,h2inv,intype,iter)

!     Solve the pressure equation by (nested) preconditioned 
!     conjugate gradient iteration.
!     INTYPE =  0  (steady)
!     INTYPE =  1  (explicit)
!     INTYPE = -1  (implicit)

      implicit none

      include 'SIZE'
      include 'MASS'
      include 'INPUT'
      include 'TSTEP'

!      include 'TOTAL'

      real divex
      common  /ctolpr/ divex

      common  /cprint/ ifprint
      logical          ifprint
      real             rcgr (lx2,ly2,lz2,lelv)
      real             rcgi (lx2,ly2,lz2,lelv)
      real             h1   (lx1,ly1,lz1,lelv)
      real             h2   (lx1,ly1,lz1,lelv)
      real             h2inv(lx1,ly1,lz1,lelv)
 
      real rpcgr,rpcgi,pcgr,pcgi,xcgr,xcgi 
      common /scrpsn1/ rpcgr(lx2*ly2*lz2*lelt)
     $               , rpcgi(lx2*ly2*lz2*lelt)
     $               , pcgr(lx2*ly2*lz2*lelt)
     $               , pcgi(lx2*ly2*lz2*lelt)
     $               , xcgr(lx2*ly2*lz2*lelt)
     $               , xcgi(lx2*ly2*lz2*lelt)

      real Apr,Api,wpr,wpi,wk1,wk2
      common /scrpsn2/   Apr(lx2,ly2,lz2,lelv)
     $ ,                 Api(lx2,ly2,lz2,lelv)
     $ ,                 wpr(lx2,ly2,lz2,lelv) 
     $ ,                 wpi(lx2,ly2,lz2,lelv)
     $ ,                 wk1(lx2,ly2,lz2,lelv)
     $ ,                 wk2(lx2,ly2,lz2,lelv)


      real*8 etime1,dnekclock
      integer*8 ntotg,nxyz2

      real div0,h1_mx,h2_mx,pap,pcgmx,ratio,rnorm
      real rpx,rpy,tolpss,rnrm1,rrp1,rrp2,wp_mx
      integer iter,iconv,nelgv,ntot1,ntot2
      integer intype

      real glsc2,glamax

      real pApr,pApi
      real alphar,alphai
      real betar,betai
      real rrp1r,rrp1i
      real rrp2r,rrp2i
      real const

      logical ifprec

      etime1 = dnekclock()
      divex = 0.
      iter  = 0

      ifprec = .true.

!      call chktcg2 (tolps,rcg,iconv)
      if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
     $   tolps = abs(param(21))

      nxyz2 = lx2*ly2*lz2
      ntot2 = nxyz2*nelv
      ntotg = nxyz2*nelgv

!     x0 = 0.      
      call rzero(xcgr,ntot2)
      call rzero(xcgi,ntot2)

!     r0 = b - Ax0 = b

!     z0 = (M^-1)r0      ==> rpcg = (M^-1)r0
      if (ifprec) then
        call uzprec(rpcgr,rcgr,h1,h2,intype,wpr)
        call uzprec(rpcgi,rcgi,h1,h2,intype,wpi)
      else  
        call copy(rpcgr,rcgr,ntot2)
        call copy(rpcgi,rcgi,ntot2)
      endif  

!     p0 = z0 
      call copy(pcgr,rpcgr,ntot2)
      call copy(pcgi,rpcgi,ntot2)

!     Check Convergence. Also calculate rrp1 = (r^T)*z
      call convprn_f3d(iconv,rnorm,rrp1r,rrp1i,
     $                 rcgr,rcgi,rpcgr,rpcgi,tolpss)
      div0  = rnorm
      betar = 0.
      betai = 0.
      if (param(21).lt.0) tolpss = abs(param(21))*div0

      if (iconv.eq.1) goto 9000

      rrp1 = sqrt(rrp1r*rrp1r + rrp1i*rrp1i)
      if (rrp1.eq.0) return

      tolpss = tolps
      do 1000 iter=1,500 !nmxp

!       rrp1 = (r^T)*z
        call glsc2_comp(rrp1r,rrp1i,rcgr,rcgi,rpcgr,rpcgi,ntot2)

!       Ap = A*p            
        call cdabdtp_f3d(Apr,pcgr,h1,h2,h2inv,intype)
        call cdabdtp_f3d(Api,pcgi,h1,h2,h2inv,intype)

!       pAp = (p^T)*Ap         
        call glsc2_comp(pApr,pApi,pcgr,pcgi,Apr,Api,ntot2)

        pAp = pApr*pApr + pApi*pApi

!       \alpha = ((r^T)*z)/((p^T)*Ap)
        alphar = (rrp1r*pApr + rrp1i*pApi)/pAp
        alphai = (rrp1i*pApr - rrp1r*pApi)/pAp
        
        if (pAp.eq.0) then
           pcgmx = glamax(pcgr,ntot2)
           wp_mx = glamax(Apr,ntot2)
           ntot1 = lx1*ly1*lz1*nelv
           h1_mx = glamax(h1 ,ntot1)
           h2_mx = glamax(h2 ,ntot1)
           if (nid.eq.0) write(6,*) 'error: pap=0 in uzawa.'
     $     ,iter,pcgmx,wp_mx,h1_mx,h2_mx
           call exitt
        endif

!        write(6,*) 'alpha', alphar,alphai
        alphai = 0.

!       x = x + \alpha*p
!       xr = xr + \alpha_r*pr
        call add2s2(xcgr,pcgr,alphar,ntot2)
!       xr = xr - \alpha_i*pi
        call add2s2(xcgr,pcgi,-alphai,ntot2)

!       xi = xi + \alpha_r*pi
        call add2s2(xcgi,pcgi,alphar,ntot2)
!       xi = xi + \alpha_i*pr
        call add2s2(xcgi,pcgr,alphai,ntot2)


!       r = r - \alpha*Ap
!       rr = rr - \alpha_r*Apr
        const = -alphar
        call add2s2(rcgr,Apr,const,ntot2)
!       rr = rr + \alpha_i*Api
        const = alphai
        call add2s2(rcgr,Api,const,ntot2)

!       ri = ri - \alpha_r*Api
        const = -alphar
        call add2s2(rcgi,Api,const,ntot2)
!       ri = ri - \alpha_i*Apr
        const = -alphai
        call add2s2(rcgi,Apr,const,ntot2)

        rrp2r = rrp1r
        rrp2i = rrp1i
        rrp2 = rrp2r*rrp2r + rrp2i*rrp2i

!       Check Convergence        
        call convprn_f3d(iconv,rnorm,rpx,rpy,
     $                   rcgr,rcgi,rpcgr,rpcgi,tolpss)

        ratio = rnorm/div0
        if (ifprint.and.nio.eq.0) 
     $  write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66   format(i5,1p4e12.5,i8,' Divergence')
        if (iconv.eq.1) goto 9000

        call ortho(rcgr)
        call ortho(rcgi)

!        rrp2r = rrp1r
!        rrp2i = rrp1i

!       z = (M^-1)r
        if (ifprec) then 
          call uzprec(rpcgr,rcgr,h1,h2,intype,wpr)
          call uzprec(rpcgi,rcgi,h1,h2,intype,wpi)
        else
          call copy(rpcgr,rcgr,ntot2)
          call copy(rpcgi,rcgi,ntot2)
        endif  

!       rrp1 = (r^T)*z
        call glsc2_comp(rrp1r,rrp1i,rcgr,rcgi,rpcgr,rpcgi,ntot2)

!       \beta = rrp1/rrp2 
        betar = (rrp1r*rrp2r + rrp1i*rrp2i)/rrp2
        betai = (rrp1i*rrp2r - rrp1r*rrp2i)/rrp2
        
!        write(6,*) 'Beta', betar,betai   
        betai = 0.

        call copy(wk1,pcgr,ntot2)
        call copy(wk2,pcgi,ntot2)

!       p = z + \beta*p
!       pr = zr + \beta_r*pr - \beta_i*pi
        call copy(pcgr,rpcgr,ntot2)         ! pcgr = zr
        const = betar
        call add2s2(pcgr,wk1,const,ntot2)
        const = -betai
        call add2s2(pcgr,wk2,const,ntot2)

!       pi = zi + \beta_r*pi + \beta_i*pr
        call copy(pcgi,rpcgi,ntot2)         ! pcgi = zi
        call add2s2(pcgi,wk2,betar,ntot2)
        call add2s2(pcgi,wk1,betai,ntot2)

 1000 continue
      if (nid.eq.0) write (6,3001) iter,rnorm,tolpss
c     if (istep.gt.20) call emerxit
 3001 format(I6,' **ERROR**: Failed to converge in UZAWA:',6E13.4)
 9000 continue

      divex = rnorm
      iter  = iter-1

      if (iter.gt.0) then
        call copy (rcgr,xcgr,ntot2)
        call copy (rcgi,xcgi,ntot2)
      endif  

      call ortho(rcgr)
      call ortho(rcgi)

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,9999) istep, '  U-Press std. ',
     &                            iter,divex,div0,tolpss,etime1
 9999 format(I11,a,I7,1p4E13.4)
19999 format(I11,'  U-Press 1.e-5: ',I7,1p4E13.4)

      return
      end subroutine uzawa_f3d
!-----------------------------------------------------------------------
      subroutine convprn_f3d(iconv,rbnorm,rrptr,rrpti,
     $                       resr,resi,zr,zi,tol)

!                                               T
!     Convergence test for the pressure step;  r z

      implicit none

      include 'SIZE'
      include 'MASS'

      real resr(1),resi(1)
      real zr(1),zi(1)
      real wrk1(4),wrk2(4)
      real rrptr,rrpti,rbnorm
      real tol

      integer iconv
      integer ntot2
      real vlsc21
      real vlsc2

      ntot2   = lx2*ly2*lz2*nelv

      wrk1(1) = vlsc21(resr,bm2inv,ntot2)  !  resr*bm2inv*resr
      wrk1(2) = vlsc21(resi,bm2inv,ntot2)  !  resi*bm2inv*resi

      call gop(wrk1,wrk2,'+  ',2)
      rbnorm  = sqrt((wrk1(1)+wrk1(2))/volvm2) ! sqrt((rr² + ri²)/vol)
      
      call glsc2_comp(rrptr,rrpti,resr,resi,zr,zi,ntot2)

!      wrk1(1) = vlsc2(resr,zr,ntot2)     !  resr*zr
!      wrk1(2) = vlsc2(resi,zi,ntot2)     !  resi*zi
!      wrk1(3) = vlsc2(resr,zi,ntot2)     !  resr*zi
!      wrk1(4) = vlsc2(resi,zr,ntot2)     !  resi*zr
!
!      call gop(wrk1,wrk2,'+  ',4)
!      rrptr = wrk1(1) + wrk1(2)
!      rrpti = wrk1(3) - wrk1(4)

      iconv  = 0
      if (rbnorm.lt.tol) iconv=1
      return
      end subroutine convprn_f3d
!---------------------------------------------------------------------- 

      subroutine glsc2_comp(scr,sci,ur,ui,vr,vi,n)
!                _       
!     Calculate: u * v
!     Assuming the Weight itself is real        

      implicit none

      integer i,n

      real ur(n),ui(n)
      real vr(n),vi(n)

      real scr,sci,wk(2),wk2(2)
      real glsum

      wk(1) = 0.
      wk(2) = 0.
      do i=1,n
        wk(1) = wk(1) + ur(i)*vr(i) + ui(i)*vi(i)
        wk(2) = wk(2) + ur(i)*vi(i) - ui(i)*vr(i)
      enddo

!     Sum over all processors
      call gop(wk,wk2,'+  ',2)
      scr = wk(1)
      sci = wk(2)

      return
      end subroutine glsc2_comp

!-----------------------------------------------------------------------
      subroutine esolver_f3d(resr,resi,h1,h2,h2inv,intype)

!     Choose E-solver

      implicit none        

      include 'SIZE'
      include 'ESOLV'
      include 'INPUT'
      include 'CTIMER'

      real resr(lx2,ly2,lz2,lelv)
      real resi(lx2,ly2,lz2,lelv)
     
      real h1(lx1,ly1,lz1,lelv)
      real h2(lx1,ly1,lz1,lelv)
      real h2inv(lx1,ly1,lz1,lelv)

      integer intype
      integer icg

      if (icalld.eq.0) teslv=0.0

      call ortho(resr) ! Ensure that residual is orthogonal to null space
      call ortho(resi) ! Ensure that residual is orthogonal to null space

      icalld=icalld+1
      neslv=icalld
      etime1=dnekclock()

      if (.not. ifsplit) then
        if (param(42).eq.1) then
          call uzawa_f3d(resr,resi,h1,h2,h2inv,intype,icg)
        else
!          call uzawa_gmres(res,h1,h2,h2inv,intype,icg)
          write(6,*) 'ERROR: E-solver (F3D) not implemented for GMRES'
          call exitt
        endif
      else
         write(6,*) 'ERROR: E-solver does not exist PnPn'
         call exitt
      endif

      teslv=teslv+(dnekclock()-etime1)

      return
      end subroutine esolver_f3d
!-----------------------------------------------------------------------
      subroutine ophinv_real(o1,o2,o3,i1,i2,i3,h1,h2,tolh,nmxhi)
C
C     Ok = (H1*A+H2*B)-1 * Ik  (implicit)
C
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'VPROJ'
      include 'TSTEP'

      include 'F3D'
c      logical ifproj
 
      real o1 (lx1,ly1,lz1,1) , o2 (lx1,ly1,lz1,1) , o3 (lx1,ly1,lz1,1)
      real i1 (lx1,ly1,lz1,1) , i2 (lx1,ly1,lz1,1) , i3 (lx1,ly1,lz1,1)
      real h1 (lx1,ly1,lz1,1) , h2 (lx1,ly1,lz1,1)
      
      integer mtmp,i,matmod,nmxhi
      real tolh

      mtmp = param(93)
      do i=1,2*ldim
         ivproj(1,i) = min(mxprev,mtmp) - 1
      enddo
 
      imesh = 1
 
      if (ifstrs) then
         matmod = 0
         call hmhzsf_real ('NOMG',o1,o2,o3,i1,i2,i3,h1,h2,
     $                  v1mask,v2mask,v3mask,vmult,
     $                  tolh,nmxhi,matmod)
      else

        if (ifield.eq.1) then
          if (nid.eq.0) then 
            write(6,*) 'IFSTRS = ', ifstrs
            write(6,*) 'Cylindrical Solver needs coupled solve'
            write(6,*) 'Set IFSTRS = TRUE'
          endif  

          call exitt
        endif  

      endif
C
      return
      end subroutine ophinv_real
!----------------------------------------------------------------------
      subroutine hmhzsf_real (name,u1,u2,u3,r1,r2,r3,h1,h2,
     $                   rmask1,rmask2,rmask3,rmult,
     $                   tol,maxit,matmod)

c     Solve coupled Helmholtz equations (stress formulation)

      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'   ! For outpost diagnostic call
      include 'TSTEP'
      include 'ORTHOSTRS'
      include 'CTIMER'

      real u1(1),u2(1),u3(1),r1(1),r2(1),r3(1),h1(1),h2(1)
      real rmask1(1),rmask2(1),rmask3(1),rmult(1)
      character name*4

      integer iproj,nel,n,maxit,matmod
      real tol
      real vol

#ifdef TIMER
      nhmhz = nhmhz + 1
      etime1 = dnekclock()
#endif

      nel = nelfld(ifield)
      vol = volfld(ifield)
      n   = lx1*ly1*lz1*nel

      napproxstrs(1) = 0
      iproj = 0
      if (ifprojfld(ifield)) iproj = param(94)
      if (iproj.gt.0.and.istep.ge.iproj) napproxstrs(1)=param(93)
      napproxstrs(1)=min(napproxstrs(1),mxprev)

!      call rmask   (r1,r2,r3,nel)
!      call opdssum (r1,r2,r3)
      call col2_3  (r1,r2,r3,v1mask,v2mask,v3mask,n)
      call dssum3  (r1,r2,r3)

      call rzero3  (u1,u2,u3,n)

      if (imesh.eq.1) then
!         call chktcgs (r1,r2,r3,rmask1,rmask2,rmask3,rmult,binvm1
!     $                ,vol,tol,nel)

!         call strs_project_a(r1,r2,r3,h1,h2,rmult,ifield,ierr,matmod)

         call cggosf_real(u1,u2,u3,r1,r2,r3,h1,h2,rmult,binvm1
     $                ,vol,tol,maxit,matmod)

!         call strs_project_b(u1,u2,u3,h1,h2,rmult,ifield,ierr)

      else

      endif

#ifdef TIMER
      thmhz=thmhz+(dnekclock()-etime1)
#endif

      return
      end subroutine hmhzsf_real

!--------------------------------------------------------------------

      subroutine cggosf_real (u1,u2,u3,r1,r2,r3,h1,h2,rmult,binv,
     $                   vol,tin,maxit,matmod)

C     Conjugate gradient iteration for solution of coupled 
C     Helmholtz equations 

      implicit none

      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
      include 'FDMH1'
      include 'F3D'

      include 'TEST'    ! just testing

      real dpc,p1,p2,p3
      common /screv11/  dpc(lx1*ly1*lz1*lelt)
     $     ,            p1 (lx1*ly1*lz1*lelt)
      common /scrch11/  p2 (lx1*ly1*lz1*lelt)
     $     ,            p3 (lx1*ly1*lz1*lelt)
      real qq1,qq2,qq3
      common /scrsl11/  qq1(lx1*ly1*lz1*lelt)
     $     ,            qq2(lx1*ly1*lz1*lelt)
     $     ,            qq3(lx1*ly1*lz1*lelt)
      real pp1,pp2,pp3,wa
      common /scrmg11/  pp1(lx1*ly1*lz1*lelt)
     $     ,            pp2(lx1*ly1*lz1*lelt)
     $     ,            pp3(lx1*ly1*lz1*lelt)
     $     ,            wa (lx1*ly1*lz1*lelt)
      real ap1(1),ap2(1),ap3(1)
      equivalence (ap1,pp1),(ap2,pp2),(ap3,pp3)

      common /fastmd/ ifdfrm(lelt), iffast(lelt), ifh2, ifsolv
      common /cprint/ ifprint
      logical ifdfrm, iffast, ifh2, ifsolv, ifprint

      real u1(1),u2(1),u3(1),
     $     r1(1),r2(1),r3(1),h1(1),h2(1),rmult(1),binv(1)

      logical iffdm,ifcrsl,ifjacobi

      integer maxit,matmod,nel,nxyz,n,iter,ifin

      real tin,tol,vol
      real alpha,beta,pap
      real rpp1,rpp2

      real tmpval

      real r0,rbnorm

      real glsc2,glsc3,op_glsc2_wt
      real opnorm2_wt_comp

!     No Fast Diagonalization Method            
      iffdm  = .false.
!     Jacobi preconditioner 
      ifjacobi = .true.
!     No Coarse grid
      ifcrsl = .false.

      if (istep.eq.1) ifprint = .true.

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

!      call opdot (wa,r1,r2,r3,r1,r2,r3,n)
      if (iff3d) then
        call vdot3 (wa,r1,r2,r3,r1,r2,r3,n)
      else
        call opdot (wa,r1,r2,r3,r1,r2,r3,n)
      endif  

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
      if (ifjacobi) then
        if (ifcyl_f3d) then
          call setprec_cyl (qq1,h1,h2,k_f3d,imesh,1)
          call setprec_cyl (qq2,h1,h2,k_f3d,imesh,2)
          call setprec_cyl (qq3,h1,h2,k_f3d,imesh,3)
        else
          call setprec (qq1,h1,h2,imesh,1)
          call setprec (qq2,h1,h2,imesh,1)
          call setprec (qq3,h1,h2,imesh,1)
        endif
      else
        call rone(qq1,n)
        call rone(qq2,n)
        call rone(qq3,n)
      endif
  
!      call setprec (dpc,h1,h2,imesh,1)
!      call setprec (wa ,h1,h2,imesh,2)
!      call add2    (dpc,wa,n)
!      if (iff3d) then
!         call setprec (wa,h1,h2,imesh,3)
!         call add2    (dpc,wa,n)
!      endif

      if (iffdm) then
!         call set_fdm_prec_h1b(dpc,h1,h2,nel)
!         call fdm_h1a (pp1,r1,dpc,nel,ktype(1,1,1),wa)
!         call fdm_h1a (pp2,r2,dpc,nel,ktype(1,1,2),wa)
!         call fdm_h1a (pp3,r3,dpc,nel,ktype(1,1,3),wa)
!         call rmask   (pp1,pp2,pp3,nel)
!         call opdssum (pp1,pp2,pp3)
      elseif (ifjacobi) then
         call col3 (pp1,qq1,r1,n)
         call col3 (pp2,qq2,r2,n)
         call col3 (pp3,qq3,r3,n)
      else
         call copy3(pp1,pp2,pp3,r1,r2,r3,n)   
      endif

      if (ifcrsl) then
!         call crs_strs(p1,p2,p3,r1,r2,r3)
!         call rmask   (p1,p2,p3,nel)
      else
         call rzero3(p1,p2,p3,n)
      endif
!      call opadd2       (p1,p2,p3,pp1,pp2,pp3)
      call add2_3(p1,p2,p3,pp1,pp2,pp3,n)
      if3d = .true.
      rpp1 = op_glsc2_wt(p1,p2,p3,r1,r2,r3,rmult)
      if3d = .false.

      maxit=500
      do 1000 iter=1,maxit
         call axhmsf_cyl_real(ap1,ap2,ap3,p1,p2,p3,h1,h2)
!         call axhmsf(ap1,ap2,ap3,p1,p2,p3,h1,h2,matmod)

!         call rzero3(tmp1,tmp2,tmp3,n)
!         call rzero3(tmp5,tmp6,tmp7,n)   
!         call axhmsf_cyl(ap1,ap2,ap3,tmp1,tmp2,tmp3,
!     $                   p1,p2,p3,tmp5,tmp6,tmp7,h1,h2)

!         call rmask   (ap1,ap2,ap3,nel)
!         call opdssum (ap1,ap2,ap3)
         call col2_3(ap1,ap2,ap3,v1mask,v2mask,v3mask,n)
         call dssum3 (ap1,ap2,ap3)

         if3d  = .true.
         pap   = op_glsc2_wt(p1,p2,p3,ap1,ap2,ap3,rmult)
         if3d  = .false.
         alpha = rpp1 / pap

!         call opadds (u1,u2,u3,p1 ,p2 ,p3 , alpha,n,2)
!         call opadds (r1,r2,r3,ap1,ap2,ap3,-alpha,n,2)
         call opadds_3 (u1,u2,u3,p1,p2,p3,alpha,n,2)
         call opadds_3 (r1,r2,r3,ap1,ap2,ap3,-alpha,n,2)

!         call opdot  (wa,r1,r2,r3,r1,r2,r3,n)
         if (iff3d) then
           call vdot3 (wa,r1,r2,r3,r1,r2,r3,n)
         else
           call opdot (wa,r1,r2,r3,r1,r2,r3,n)
         endif  

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
         elseif (ifprint) then
            if (matmod.ge.0) write(6,3000) istep,iter,rbnorm,r0,tol
         endif

         if (iffdm) then
!            call fdm_h1a (pp1,r1,dpc,nel,ktype(1,1,1),wa)
!            call fdm_h1a (pp2,r2,dpc,nel,ktype(1,1,2),wa)
!            call fdm_h1a (pp3,r3,dpc,nel,ktype(1,1,3),wa)
!            call rmask   (pp1,pp2,pp3,nel)
!            call opdssum (pp1,pp2,pp3)
         elseif (ifjacobi) then
           call col3 (pp1,qq1,r1,n)
           call col3 (pp2,qq2,r2,n)
           call col3 (pp3,qq3,r3,n)
         else  
           call copy3 (pp1,pp2,pp3,r1,r2,r3,n)
         
!            call col3 (pp1,dpc,r1,n)
!            call col3 (pp2,dpc,r2,n)
!            if (iff3d) call col3 (pp3,dpc,r3,n)
         endif

         if (ifcrsl) then
!           call crs_strs(qq1,qq2,qq3,r1,r2,r3)
!           call rmask   (qq1,qq2,qq3,nel)
!           call opadd2  (pp1,pp2,pp3,qq1,qq2,qq3)
         endif

!         call opdot (wa,r1,r2,r3,pp1,pp2,pp3,n)
         if (iff3d) then
           call vdot3 (wa,r1,r2,r3,pp1,pp2,pp3,n)
         else
           call opdot (wa,r1,r2,r3,pp1,pp2,pp3,n)
         endif  

         rpp2 = rpp1
         rpp1 = glsc2(wa,rmult,n)
         beta = rpp1/rpp2
!         call opadds (p1,p2,p3,pp1,pp2,pp3,beta,n,1)
         call opadds_3 (p1,p2,p3,pp1,pp2,pp3,beta,n,1)

 1000 continue
      if (matmod.ge.0.and.nio.eq.0) write (6,3001) 
     $                              istep,iter,rbnorm,r0,tol
      if (matmod.lt.0.and.nio.eq.0) write (6,3011) 
     $                              istep,iter,rbnorm,r0,tol

 9999 continue
      ifsolv = .false.


 3000 format(i11,'  Helmh3 fluid  ',I6,1p3E13.4)
 3002 format(i11,a,1x,I7,1p4E13.4,l4)
 3010 format(i11,'  Helmh3 mesh   ',I6,1p3E13.4)
 3001 format(i11,'  Helmh3 fluid unconverged! ',I6,1p3E13.4)
 3011 format(i11,'  Helmh3 mesh unconverged! ',I6,1p3E13.4)

      return
      end subroutine cggosf_real
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


