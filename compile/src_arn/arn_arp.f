!> @file arn_arp.f
!! @ingroup arn_arp
!! @brief Set of subroutines to solve eigenvalue problem with Arnoldi
!!   algorithm using PARPACK/ARPACK
!! @warning There is no restart option for serial ARPACK version. It is
!!   supported by parallel PARPACK only.
!! @author Adam Peplinski
!! @date Mar 7, 2016
!
!
! To define ARPACK mode: direct or inverse:
! Notice that simulation with temperature or passive scalars has to be
! performed in inverse mode due to speciffic inner product.
!
! Direct eigenvalue problem A*x = lambda*x
!#define ARPACK_DIRECT
! Generalized (inverse) eigenvalue problem A*x = lambda*B*x
#undef ARPACK_DIRECT
!=======================================================================
!> @brief Register Arnoldi ARPACK module
!! @ingroup arn_arp
!! @note This interface is called by @ref tst_register
      subroutine stepper_register()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'TSTEPPERD'
      include 'ARN_ARPD'

      ! local variables
      integer lpmid, il
      real ltim
      character*2 str

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

      ! check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid,arna_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(arna_name)//'] already registered')
         return
      endif

      ! find parent module
      call mntr_mod_is_name_reg(lpmid,tst_name)
      if (lpmid.le.0) then
         lpmid = 1
         call mntr_abort(lpmid,
     $        'parent module ['//trim(tst_name)//'] not registered')
      endif

      ! register module
      call mntr_mod_reg(arna_id,lpmid,arna_name,
     $      'Arnoldi ARPACK spectra calculation')

      ! register timers
      ! initialisation
      call mntr_tmr_reg(arna_tmr_ini_id,tst_tmr_ini_id,arna_id,
     $     'ARNA_INI','Arnoldi ARPACK initialisation time',.true.)
      ! submodule operation
      call mntr_tmr_reg(arna_tmr_evl_id,tst_tmr_evl_id,arna_id,
     $     'ARNA_EVL','Arnoldi ARPACK evolution time',.true.)

      ! register and set active section
      call rprm_sec_reg(arna_sec_id,arna_id,'_'//adjustl(arna_name),
     $     'Runtime paramere section for Arnoldi ARPACK module')
      call rprm_sec_set_act(.true.,arna_sec_id)

      ! register parameters
      call rprm_rp_reg(arna_nkrl_id,arna_sec_id,'NKRL',
     $     'Krylov space size',rpar_int,50,0.0,.false.,' ')

      call rprm_rp_reg(arna_negv_id,arna_sec_id,'NEGV',
     $     'Number of eigenvalues',rpar_int,10,0.0,.false.,' ')

      call rprm_rp_reg(arna_ifpr_id,arna_sec_id,'IFPR',
     $     'IF Pressure?',rpar_log,0,0.0,.false.,' ')

      call rprm_rp_reg(arna_ifcomplex_id,arna_sec_id,'IFCOMPLEX',
     $     'IF Complex?',rpar_log,0,0.0,.false.,' ')

      ! set initialisation flag
      arna_ifinit=.false.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(arna_tmr_ini_id,1,ltim)

      return
      end subroutine
!=======================================================================
!> @brief Initilise Arnoldi ARPACK module
!! @ingroup arn_arp
!! @note This interface is called by @ref tst_init
      subroutine stepper_init()
      implicit none

      include 'SIZE'            ! NIO, NDIM, NPERT
      include 'TSTEP'           ! NSTEPS
      include 'INPUT'           ! IF3D, IFHEAT
      include 'SOLN'            ! V?MASK, TMASK, V[XYZ]P, TP
      include 'FRAMELP'
      include 'CHKPOINTD'
      include 'TSTEPPERD'
      include 'ARN_ARPD'

      include 'F3D'

      ! ARPACK include file
      INCLUDE 'debug.h'

      ! local variables
      integer itmp, il
      real rtmp, ltim
      logical ltmp
      character*20 ctmp

      ! functions
      real dnekclock

      integer i
!-----------------------------------------------------------------------
      ! check if the module was already initialised
      if (arna_ifinit) then
         call mntr_warn(arna_id,
     $        'module ['//trim(arna_name)//'] already initiaised.')
         return
      endif

      ! timing
      ltim = dnekclock()

      ! get runtime parameters
!     if pressure
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,arna_ifpr_id,rpar_log)
      arna_ifpr = ltmp

!     if complex
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,arna_ifcomplex_id,rpar_log)
      arna_ifcomplex = ltmp

!     NKRYL
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,arna_nkrl_id,rpar_int)
      arna_nkrl = itmp

!     NEV      
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,arna_negv_id,rpar_int)
      arna_negv = itmp

!     get restart options
!     if restart      
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,chpt_ifrst_id,rpar_log)
      arna_ifrst = ltmp

!     restart file number      
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,chpt_fnum_id,rpar_int)
      arna_fnum = itmp

      if (.not.tst_iftst) return

      ! check simulation parameters
#ifdef ARPACK_DIRECT
      ! standard eigenvalue problem A*x = lambda*x
      if(IFHEAT) call mntr_abort(arna_id,
     $   'IFHEAT requires #undef ARPACK_DIRECT')
#endif

      if (arna_nkrl.gt.arna_lkrl) call mntr_abort(arna_id,
     $   'arna_nkrl bigger than arna_lkrl')

      if (arna_negv.ge.(arna_nkrl/2)) call mntr_abort(arna_id,
     $   'arna_negv > arna_nkrl/2')

      ! make sure NSTEPS is bigger than the possible number of iteration in arnoldi
      ! multiplication by 2 for OIC
      if (tst_iftst) then
        NSTEPS = max(NSTEPS,tst_step*arna_nkrl*tst_cmax*2+10)
      endif  

      ! related to restart
      nparp = 0
      ncarp = 0
      rnmarp= 0.0

      ! initialise ARPACK parameters
#ifdef ARPACK_DIRECT
      ! direct eigenvalue problem A*x = lambda*x
      bmatarp='I'
#else
      ! generalised eigenvalue problem A*x = lambda*B*x
      bmatarp='G'
#endif
      ! eigenvalues of largest magnitude
      whicharp='LM'

      call izero(iparp,11)
      call izero(ipntarp,14)
      ! exact shifts with respect to the current Hessenberg matrix
      iparp(1)=1
      ! maximum number of Arnoldi update iterations allowed
      iparp(3)=tst_cmax
#ifdef ARPACK_DIRECT
      ! A*x = lambda*x
      iparp(7)=1
#else
      ! A*x = lambda*M*x, M symmetric positive definite; bmatarp='G'
      iparp(7)=2
#endif
      ! used size of workla
      nwlarp = (3*arna_nkrl+6)*arna_nkrl

      ! user supplied initial conditions
      infarp=1

      ! get eigenvectors
      rvarp=.true.
      ! compute Ritz vectors
      howarp='A'
      ! select should be specified for howarp='S'

      ! no shift
      if (arna_ifcomplex) then
        sigarp(1) = cmplx(0.,0.)
      else  
        sigarp(1) = 0.0
        sigarp(2) = 0.0
      endif  

      ! vector lengths
      ! single vector length in Krylov space
      ! velocity
      if (iff3d) then
        arna_ns = tst_nv*3 !          ! velocity
      else
        arna_ns = tst_nv*ndim
      endif
      if (arna_ifpr) arna_ns = arna_ns + tst_np    ! pressure
      ! temperature
      if(IFHEAT) then
         arna_ns = arna_ns + tst_nt ! temperature
      endif
      if (arna_ns.gt.arna_ls) then
        call mntr_abort(arna_id,
     $   'arna_ns too big; arna_ns > arna_ls')
      endif  

      ! initialise arrays
      if (arna_ifcomplex) then
        call czero(workda,wddima)
        call czero(workla,wldima)
        call czero(workea,wedima)
        call czero(vbasea,arna_ls*arna_lkrl)
        call czero(resida,arna_ls)

!       This is real        
        call rzero(driarp,arna_lkrl*4)
      else
        call rzero(workda,wddima)
        call rzero(workla,wldima)
        call rzero(workea,wedima)
        call rzero(vbasea,arna_ls*arna_lkrl)
        call rzero(resida,arna_ls)
        call rzero(driarp,arna_lkrl*4)
      endif  

!     info level from ARPACK
      ndigit = -3
      logfil = 6

!     No Outputs from Arpack on all processes 
!     Symmetric codes
      msgets = 0
      msaitr = 0
      msapps = 0
      msaupd = 0
      msaup2 = 0
      mseupd = 0

!     Non-Symmetric codes      
      mngets = 0
      mnaitr = 0
      mnapps = 0
      mnaupd = 0
      mnaup2 = 0
      mneupd = 0

!     Complex subroutine codes      
      mcgets = 0
      mcaitr = 0
      mcapps = 0
      mcaupd = 0
      mcaup2 = 0
      mceupd = 0

!     Only process 0 outputs 
      if (nio.eq.0) then
!       Symmetric codes
        msgets = 0
        msaitr = 2
        msapps = 0
        msaupd = 2
        msaup2 = 3
        mseupd = 3

!       Non-Symmetric codes      
        mngets = 0
        mnaitr = 2
        mnapps = 0
        mnaupd = 2
        mnaup2 = 3
        mneupd = 3

!       Complex subroutine codes      
        mcgets = 0
        mcaitr = 2
        mcapps = 0
        mcaupd = 2
        mcaup2 = 3
        mceupd = 3
      endif  

      ! restart
      if (arna_ifrst) then
         ! read checkpoint
         call arn_rst_read
      else
         ! if no restart fill RESIDA with initial conditions
         ! V?MASK removes points at the wall and inflow
#ifdef ARPACK_DIRECT
         ! A*x = lambda*x
         ! velocity
         if (arna_ifcomplex) then
           i = 1
           call copytocomplex(resida(i),vxp(1,1),vxp(1,2),tst_nv)
           call col2_cr(resida(i),v1mask,tst_nv) 
           i = i + tst_nv

           call copytocomplex(resida(i),vyp(1,1),vyp(1,2),tst_nv)
           call col2_cr(resida(i),v2mask,tst_nv) 
           i = i + tst_nv
           
           if (iff3d) then
             call copytocomplex(resida(i),vzp(1,1),vzp(1,2),tst_nv)
             call col2_cr(resida(i),v3mask,tst_nv)
             i = i + tst_nv
           endif  
            
           if (arna_ifpr) then
             call copytocomplex(resida(i),prp(1,1),prp(1,2),tst_np)
             i = i + tst_np
           endif  
         else           ! real arithmetic
           i = 1
           call col3(resida(i),VXP,V1MASK,tst_nv)
           i = i + tst_nv
           call col3(resida(i),VYP,V2MASK,tst_nv)
           i = i + tst_nv
           if (IF3D.or.iff3d) then
             call col3(resida(i),VZP,V3MASK,tst_nv)
             i = i + tst_nv
           endif
           if (arna_ifpr) then 
             call copy(resida(i),prp,tst_np)
             i = i + tst_np
           endif  
         endif          ! arna_ifcomplex             

         ! no temperature here
#else
         ! A*x = lambda*M*x
         if (arna_ifcomplex) then   
           i = 1
           call copytocomplex(resida(i),vxp(1,1),vxp(1,2),tst_nv)
           call col2_cr(resida(i),v1mask,tst_nv) 
           i = i + tst_nv

           call copytocomplex(resida(i),vyp(1,1),vyp(1,2),tst_nv)
           call col2_cr(resida(i),v2mask,tst_nv) 
           i = i + tst_nv
           
           if (iff3d) then
             call copytocomplex(resida(i),vzp(1,1),vzp(1,2),tst_nv)
             call col2_cr(resida(i),v3mask,tst_nv)
             i = i + tst_nv
           endif  
            
           if (arna_ifpr) then
             call copytocomplex(resida(i),prp(1,1),prp(1,2),tst_np)
             i = i + tst_np
           endif

           if (ifheat) then
             call copytocomplex(resida(i),tp(1,1,1),tp(1,1,2),tst_nt)
             call col2_cr(resida(i),tmask,tst_nt)
             i = i + tst_nt
           endif 

         else           ! real arithmetic
           ! velocity
           i = 1
           call col3(resida(i),VXP,v1mask,tst_nv)
           i = i + tst_nv
           call col3(resida(i),VYP,v2mask,tst_nv)
           i = i + tst_nv
           if (IF3D.or.iff3d) then
             call col3(resida(i),VZP,v3mask,tst_nv)
             i = i + tst_nv
           endif
           if (arna_ifpr) then
             call copy(resida(i),PRP,tst_np)
             i = i + tst_np
           endif  
           ! temperature
           if (IFHEAT) then
             call col3(resida(i),TP,tmask,tst_nt)
             i = i + tst_nt
           endif  
         endif          ! arna_ifcomplex

#endif

         ! initialise rest of variables
         ! first call
         idoarp=0
      endif

      ! ARPACK interface
      call arn_naupd

      ! we should start stepper here
      if (idoarp.ne.-1.and.idoarp.ne.1) then
         write(ctmp,*) idoarp
         call mntr_abort(arna_id,
     $   'stepper_init; error with arn_naupd, ido = '//trim(ctmp))
      endif

      ! print info
      call mntr_log(arna_id,lp_prd,'ARPACK initialised')
      call mntr_log(arna_id,lp_prd,'Parameters:')
      call mntr_log(arna_id,lp_prd,'BMAT = '//trim(bmatarp))
      call mntr_log(arna_id,lp_prd,'WHICH = '//trim(whicharp))
      call mntr_logr(arna_id,lp_prd,'TOL = ',tst_tol)
      call mntr_logi(arna_id,lp_prd,'NEV = ',arna_negv)
      call mntr_logi(arna_id,lp_prd,'NCV = ',arna_nkrl)
      call mntr_logi(arna_id,lp_prd,'IPARAM(1) = ',iparp(1))
      call mntr_logi(arna_id,lp_prd,'IPARAM(3) = ',iparp(3))
      call mntr_logi(arna_id,lp_prd,'IPARAM(7) = ',iparp(7))
      call mntr_logl(arna_id,lp_prd,'RVEC = ',rvarp)
      call mntr_log(arna_id,lp_prd,'HOWMNY = '//trim(howarp))

      ! everything is initialised
      arna_ifinit=.true.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(arna_tmr_ini_id,1,ltim)

      return
      end subroutine
!=======================================================================
!> @brief Check if module was initialised
!! @ingroup arn_arp
!! @return stepper_is_initialised
      logical function stepper_is_initialised()
      implicit none

      include 'SIZE'
      include 'ARN_ARPD'
!-----------------------------------------------------------------------
      stepper_is_initialised = arna_ifinit

      return
      end function
!=======================================================================
!> @brief Create Krylov space, get Ritz values and restart
!!  stepper phase.
!! @ingroup arn_arp
!! @note This interface is called by @ref tst_solve
      subroutine stepper_vsolve
      implicit none

      include 'SIZE'            ! NIO
      include 'INPUT'           ! IFHEAT, IF3D
      include 'SOLN'            ! V[XYZ]P, TP, V?MASK
      include 'MASS'            ! BM1
      include 'FRAMELP'
      include 'TSTEPPERD'
      include 'ARN_ARPD'

      include 'F3D'

      ! local variables
      real  ltim        ! timing
      character(20) str

      integer i

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! timing
      ltim=dnekclock()

      ! fill work array with velocity
      ! V?MASK removes points at the boundary
#ifdef ARPACK_DIRECT
      ! no temperature here
      ! A*x = lambda*x
      ! velocity
      if (arna_ifcomplex) then
        i = ipntarp(2)
        call copytocomplex(workda(i),vxp(1,1),vxp(1,2),tst_nv)
        call col2_cr(workda(i),v1mask,tst_nv)
        i = i + tst_nv
        call copytocomplex(workda(i),vyp(1,1),vyp(1,2),tst_nv)
        call col2_cr(workda(i),v2mask,tst_nv)
        i = i + tst_nv
        if (iff3d) then
          call copytocomplex(workda(i),vzp(1,1),vzp(1,2),tst_nv)
          call col2_cr(workda(i),v3mask,tst_nv)
          i = i + tst_nv
        endif
        if (arna_ifpr) then
          call copytocomplex(workda(i),prp(1,1),prp(1,2),tst_np)
          i = i + tst_np
        endif  
      else 
        i = ipntarp(2) 
        call col3(workda(i),VXP,V1MASK,tst_nv)
        i = i + tst_nv
        call col3(workda(i),VYP,V2MASK,tst_nv)
        i = i + tst_nv
        if (IF3D.or.iff3d) then 
          call col3(workda(i),VZP,V3MASK,tst_nv)
          i = i + tst_nv
        endif  
      endif  
#else
      ! velocity
      ! A*x = lambda*M*x
      if (arna_ifcomplex) then
        i = ipntarp(2)
        call copytocomplex(workda(i),vxp(1,1),vxp(1,2),tst_nv)
        call col2_cr(workda(i),v1mask,tst_nv)
        i = i + tst_nv
        call copytocomplex(workda(i),vyp(1,1),vyp(1,2),tst_nv)
        call col2_cr(workda(i),v2mask,tst_nv)
        i = i + tst_nv
        if (iff3d) then
          call copytocomplex(workda(i),vzp(1,1),vzp(1,2),tst_nv)
          call col2_cr(workda(i),v3mask,tst_nv)
          i = i + tst_nv
        endif
        if (arna_ifpr) then
          call copytocomplex(workda(i),prp(1,1),prp(1,2),tst_np)
          i = i + tst_np
        endif
        if (ifheat) then
          call copytocomplex(workda(i),tp(1,1,1),tp(1,1,2),tst_nt)
          call col2_cr(workda(i),tmask,tst_nt)
          i = i + tst_nt
        endif
      else
        i = ipntarp(2) 
        call col3(workda(i),VXP,V1MASK,tst_nv)
        i = i + tst_nv
        call col3(workda(i),VYP,V2MASK,tst_nv)
        i = i + tst_nv
        if (IF3D.or.iff3d) then
          call col3(workda(i),VZP,V3MASK,tst_nv)
          i = i + tst_nv
        endif
!       pressure
        if (arna_ifpr) then
          call copy(workda(i),prp,tst_np)        
          i = i + tst_np
        endif  
!       temperature
        if (IFHEAT) then
          call col3(workda(i),TP,TMASK,tst_nt)
          i = i + tst_nt
        endif  
        ! this may be not necessary, but ARPACK manual is not clear about it
        !call col3(workda(ipntarp(1)),VXP,BM1,tst_nv)
        !call col3(workda(ipntarp(1)+tst_nv),VYP,BM1,tst_nv)
        !if (IF3D) call col3(workda(ipntarp(1)+2*tst_nv),VZP,BM1,tst_nv)
      endif       ! arna_ifcomplex 
#endif

!     ARPACK interface
      call arn_naupd

      if (idoarp.eq.-2) then
         ! checkpoint
         call arn_rst_save
      elseif (idoarp.eq.99) then
         ! finalise
         call arn_esolve
      elseif (idoarp.eq.-1.or.idoarp.eq.1) then
         ! stepper restart, nothing to do
      else
         write(str,*) idoarp
         call mntr_abort(arna_id,
     $    'stepper_vsolve; error with arn_naupd, ido = '//trim(str))
      endif

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(arna_tmr_evl_id,1,ltim)

      return
      end
!=======================================================================
!> @brief ARPACK postprocessing
!! @ingroup arn_arp
      subroutine arn_esolve
      implicit none

      include 'SIZE'            ! NIO, NID, LDIMT1
      include 'TSTEP'           ! ISTEP, DT, LASTEP
      include 'SOLN'            ! VX, VY, VZ, VMULT, V?MASK
      include 'INPUT'           ! IFXYO,IFPO,IFVO,IFTO,IFPSO,IF3D,IFHEAT
      include 'FRAMELP'
      include 'TSTEPPERD'
      include 'ARN_ARPD'

      include 'F3D'             ! iff3d
      ! local variables
      integer il, iunit, ierror
      real dumm
      logical lifxyo, lifpo, lifvo, lifto, lifpso(LDIMT1)
      character(20) str

      ! global comunication in nekton
      integer NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL
      common /nekmpi/ NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL

      integer i
      integer nps

      real eigr, eigi
      real ritzr, ritzi
     
!-----------------------------------------------------------------------
      if (idoarp.eq.99) then

         call mntr_logi(arna_id,lp_prd,
     $       'Postprocessing converged eigenvectors NV= ',iparp(5))

#ifdef MPI
         if (arna_ifcomplex) then            
           call pzneupd(NEKCOMM,rvarp,howarp,selarp,driarp,
     $       vbasea,arna_ls,sigarp,workea,bmatarp,arna_ns,
     $       whicharp,arna_negv,tst_tol,resida,arna_nkrl,vbasea,
     $       arna_ls,iparp,ipntarp,workda,workla,nwlarp,workra,ierrarp)
         else
           call pdneupd(NEKCOMM,rvarp,howarp,selarp,driarp,driarp(1,2),
     $       vbasea,arna_ls,sigarp(1),sigarp(2),workea,bmatarp,arna_ns,
     $       whicharp,arna_negv,tst_tol,resida,arna_nkrl,vbasea,
     $       arna_ls,iparp,ipntarp,workda,workla,nwlarp,ierrarp)
         endif  
#else
         if (arna_ifcomplex) then
!          From Documentation:              
!          zneupd( RVEC, HOWMNY, SELECT, D,
!                  Z, LDZ, SIGMA, WORKEV, BMAT, N,
!                  WHICH, NEV, TOL, RESID, NCV, V, 
!                  LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, RWORK, INFO )

           call zneupd(rvarp,howarp,selarp,driarp(1,1),
     $       vbasea,arna_ls,sigarp,workea,bmatarp,arna_ns,
     $       whicharp,arna_negv,tst_tol,resida,arna_nkrl,vbasea,
     $       arna_ls,iparp,ipntarp,workda,workla,nwlarp,workra,ierrarp)
         else
!          From Documentation: 
!          dneupd( RVEC, HOWMNY, SELECT, DR, DI, 
!                  Z, LDZ, SIGMAR, SIGMAI, WORKEV, BMAT, N, 
!                  WHICH, NEV, TOL, RESID, NCV, V, 
!                  LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO )

           call dneupd(rvarp,howarp,selarp,driarp,driarp(1,2),
     $       vbasea,arna_ls,sigarp(1),sigarp(2),workea,bmatarp,arna_ns,
     $       whicharp,arna_negv,tst_tol,resida,arna_nkrl,vbasea,
     $       arna_ls,iparp,ipntarp,workda,workla,nwlarp,ierrarp)
         endif  
#endif

         if (ierrarp.eq.0) then
            call mntr_log(arna_id,lp_prd,
     $       'Writing eigenvalues and eigenvectors')
            ierror=0
            ! open file
            if (NID.eq.0) then
               ! find free unit
               call io_file_freeid(iunit, ierror)
               if (ierror.eq.0) then
                  open (unit=iunit,file='eigenvalues.txt',
     $                 action='write', iostat=ierror)
                  write(unit=iunit,FMT=410,iostat=ierror)
 410              FORMAT(5x,'I',13x,'re(RITZ)',13x,'im(RITZ)',13x,
     $                 'ln|RITZ|',12x,'arg(RITZ)')
               endif
            endif
            ! error check
            call  mntr_check_abort(arna_id,ierror,
     $       'Error opening eigenvalue file.')

            ! integration time
            dumm = DT*tst_step
            dumm = 1.0/dumm

            ! copy and set output parameters
            lifxyo = IFXYO
            IFXYO = .TRUE.
            lifpo= IFPO
            IFPO = .false.
            if (arna_ifpr) IFPO = .true.
            lifvo= IFVO
            IFVO = .true.
            lifto= IFTO
            if (IFHEAT) then
               IFTO = .TRUE.
            else
               IFTO = .FALSE.
            endif
            if (iff3d) IFTO = .true.
            do il=1,LDIMT1
               lifpso(il)= IFPSO(il)
               IFPSO(il) = .false.
            enddo

!           We have to take into account storage of imaginary and real
!           parts of eigenvectors in arpack.
!           The complex Ritz vector associated with the Ritz value
!           with positive imaginary part is stored in two consecutive
!           columns.  The first column holds the real part of the Ritz
!           vector and the second column holds the imaginary part.  The
!           Ritz vector associated with the Ritz value with negative
!           imaginary part is simply the complex conjugate of the Ritz
!           vector associated with the positive imaginary part.

            ierror=0
            do il=1,IPARP(5)
               ! possible place to test error
               ! get growth rate; get eigenvalues of Operator
               if (arna_ifcomplex) then
                 ritzr = real(driarp(il,1))
                 ritzi = aimag(driarp(il,1))
                 eigr  = log(abs(driarp(il,1)))*dumm
                 eigi  = atan2(ritzi,ritzr)*dumm
                 driarp(il,2) = cmplx(eigr,eigi)
               else
                 ritzr = driarp(il,1)
                 ritzi = driarp(il,2)
                 eigr = log(sqrt(ritzr**2 +
     $              ritzi**2))*dumm
                 eigi = atan2(ritzi,ritzr)*dumm
                 driarp(il,3) = eigr
                 driarp(il,4) = eigi
               endif  

!               if (NID.eq.0)  write(unit=iunit,fmt=*,iostat=ierror)
!     $         il,driarp(il,1),driarp(il,2),driarp(il,3),driarp(il,4)
               if (NID.eq.0)  write(unit=iunit,fmt=411,iostat=ierror)
     $         il,ritzr,ritzi,eigr,eigi

 411           FORMAT(3x,I3,4(4x,E17.9E3))

               istep = il
               time  = eigi 
               !copy eigenvectors to perturbation variables
               if (arna_ifcomplex) then
                 i = 1
                 call copytoreal(vxp(1,1),vxp(1,2),vbasea(i,il),tst_nv)
                 i = i + tst_nv
                 call copytoreal(vyp(1,1),vyp(1,2),vbasea(i,il),tst_nv)
                 i = i + tst_nv
                 if (iff3d) then
                   call copytoreal(vzp(1,1),vzp(1,2),
     $                             vbasea(i,il),tst_nv)
                   i = i + tst_nv
                   call copy(tp(1,1,1),vzp(1,1),tst_nv)
                   call copy(tp(1,1,2),vzp(1,2),tst_nv)
                 endif  
                 if (arna_ifpr) then
                   call copytoreal(prp(1,1),prp(1,2),
     $                             vbasea(i,il),tst_np)
                   i = i + tst_np
                 endif  

!                real part                 
                 call outpost(vxp(1,1),vyp(1,1),vzp(1,1),prp(1,1),
     $                        tp(1,1,1),'egv')
!                Imaginary part                  
                 call outpost(vxp(1,2),vyp(1,2),vzp(1,2),prp(1,2),
     $                        tp(1,1,2),'egv')
               else
                 i = 1
                 call copy(VXP,vbasea(i,il),tst_nv)
                 i = i + tst_nv
                 call copy(VYP,vbasea(i,il),tst_nv)
                 i = i + tst_nv
                 if (IF3D.or.iff3d) then
                   call copy(VZP,vbasea(i,il),tst_nv)
                   i = i + tst_nv
                 endif
                 if (arna_ifpr) then
                   call copy(PRP,vbasea(i,il),tst_np)
                   i = i + tst_np
                 endif  
                  
                 nps = 0 
                 if (IFHEAT) then
                   nps = nps + 1
                   call copy(TP(1,nps,1),vbasea(i,il),tst_nt)
                   i = i + tst_nt
                 endif
                 if (iff3d.and..not.if3d) then
                   nps = nps + 1
                   call copy(TP(1,nps,1),vzp,tst_nv)
                   ifpso(nps-1) = .true.
                 endif 
                 call outpost2(VXP,VYP,VZP,PRP,TP,nps,'egv')
               endif          ! arna_ifcomplex

            enddo
            ! error check
            call  mntr_check_abort(arna_id,ierror,
     $       'Error writing eigenvalue file.')

            ! put output variables back
            IFXYO = lifxyo
            IFPO = lifpo
            IFVO = lifvo
            IFTO = lifto
            do il=1,LDIMT1
               IFPSO(il) = lifpso(il)
            enddo

            ! close eigenvalue file
            if (NID.eq.0)  close(unit=iunit)

         else                   ! ierrarp
            write(str,*) ierrarp
            call  mntr_abort(arna_id,
     $       'arn_esolve; error with _neupd, info = '//trim(str))
         endif                  ! ierrarp

         ! finish run
         LASTEP=1
      endif

      return
      end
!=======================================================================
!> @brief Interface to pdnaupd
!! @ingroup arn_arp
      subroutine arn_naupd
      implicit none

      include 'SIZE'            ! NIO, NDIM, N[XYZ]1
      include 'INPUT'           ! IF3D, IFHEAT
      include 'SOLN'            ! V?MASK, TMASK, V[XYZ]P, TP
      include 'MASS'            ! BM1
      include 'FRAMELP'
      include 'TSTEPPERD'
      include 'ARN_ARPD'

      include 'F3D'

      ! local variables
      character(20) str

      ! global comunication in nekton
      integer NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL
      common /nekmpi/ NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL

      integer i,ix,iy
!-----------------------------------------------------------------------
#ifdef MPI
      if (arna_ifcomplex) then
        call pznaupd(NEKCOMM,idoarp,bmatarp,arna_ns,whicharp,arna_negv,
     $    tst_tol,resida,arna_nkrl,vbasea,arna_ls,iparp,ipntarp,workda,
     $    workla,nwlarp,workra,infarp,nparp,rnmarp,ncarp)
      else
        call pdnaupd(NEKCOMM,idoarp,bmatarp,arna_ns,whicharp,arna_negv,
     $    tst_tol,resida,arna_nkrl,vbasea,arna_ls,iparp,ipntarp,workda,
     $    workla,nwlarp,infarp,nparp,rnmarp,ncarp)
      endif  
#else
      if (arna_complex) then
!     From Documentation:        
!     znaupd(IDO, BMAT, N, WHICH, NEV, 
!            TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, 
!            WORKL, LWORKL, RWORK, INFO )
       
        call znaupd(Idoarp,bmatarp,arna_ns,whicharp,arna_negv,
     $    tst_tol,resida,arna_nkrl,vbasea,arna_ls,iparp,ipntarp,workda,
     $    workla,nwlarp,workra,infarp)
      else
!     From Documentation:
!     dnaupd( IDO, BMAT, N, WHICH, NEV, 
!             TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, 
!             WORKL, LWORKL, INFO )
       
        call dnaupd(Idoarp,bmatarp,arna_ns,whicharp,arna_negv,
     $    tst_tol,resida,arna_nkrl,vbasea,arna_ls,iparp,ipntarp,workda,
     $    workla,nwlarp,infarp)
      endif  
#endif

      ! error check
      if (infarp.lt.0) then
         write(str,*) infarp
         call  mntr_abort(arna_id,
     $       'arn_naupd; error with _naupd, info = '//trim(str))
      endif

      if (idoarp.eq.2) then
!        From Documentation:         
!        IDO =  2: compute  Y = M * X  where
!                  IPNTR(1) is the pointer into WORKD for X,
!                  IPNTR(2) is the pointer into WORKD for Y.

!        multiply by weights and masks
         do
           if (arna_ifcomplex) then
             ix = ipntarp(1)
             iy = ipntarp(2) 
!            velocities 
             call ccopy(workda(iy),workda(ix),arna_ns)
             i = iy
             call col2_cr(workda(i),bm1,tst_nv)
             call col2_cr(workda(i),v1mask,tst_nv)
             i = i+tst_nv
             call col2_cr(workda(i),bm1,tst_nv)
             call col2_cr(workda(i),v2mask,tst_nv)
             i = i+tst_nv
             if (iff3d) then
               call col2_cr(workda(i),bm1,tst_nv)
               call col2_cr(workda(i),v3mask,tst_nv)
               i = i+tst_nv
             endif
!            pressure                  
             if (arna_ifpr) then
               call czero(workda(i),tst_np)
               i = i+tst_np
             endif
!            temperature
             if (ifheat) then
               call col2_cr(workda(i),bm1,tst_nt)
               call col2_cr(workda(i),tmask,tst_nt)
!              rescale coefficient of temperature
               call cmult_cr(workda(i),1.0,tst_nt)
               i = i+tst_nt
             endif
           else
             ix = ipntarp(1)
             iy = ipntarp(2) 
             i  = iy
!            velocity
             call col3(workda(i),BM1,V1MASK,tst_nv)
             i = i + tst_nv
             call col3(workda(i),BM1,V2MASK,tst_nv)
             i = i + tst_nv
             if (IF3D.or.iff3d) then
               call col3(workda(i),BM1,V3MASK,tst_nv)
               i = i + tst_nv
             endif
!            Pressure             
             if (arna_ifpr) then
               call rzero(workda(i),tst_np)
               i = i + tst_np
             endif     

!            Temperature
             if(IFHEAT) then
                call col3(workda(i),BM1,TMASK,tst_nt)
                i = i + tst_nt
!               Temperature coefficients
                call cht_weight_fun (workda(ipntarp(2)),
     $               workda(ipntarp(2)+tst_nv),
     $               workda(ipntarp(2)+2*tst_nv),
     $               workda(ipntarp(2)+NDIM*tst_nv),1.0)
             endif

             call col2(workda(iy),workda(ix),arna_ns)
           endif        ! arna_ifcomplex

#ifdef MPI
           if (arna_ifcomplex) then            
             call pznaupd(NEKCOMM,idoarp,bmatarp,arna_ns,whicharp,
     $         arna_negv,tst_tol,resida,arna_nkrl,vbasea,arna_ls,
     $         iparp,ipntarp,workda,workla,nwlarp,workra,infarp,
     $         nparp,rnmarp,ncarp)
           else
             call pdnaupd(NEKCOMM,idoarp,bmatarp,arna_ns,whicharp,
     $         arna_negv,tst_tol,resida,arna_nkrl,vbasea,arna_ls,
     $         iparp,ipntarp,workda,workla,nwlarp,infarp,nparp,rnmarp,
     $         ncarp)
           endif  
#else
           if (arna_ifcomplex) then
             call znaupd(idoarp,bmatarp,arna_ns,whicharp,arna_negv,
     $         tst_tol,resida,arna_nkrl,vbasea,arna_ls,iparp,ipntarp,
     $         workda,workla,nwlarp,workra,infarp)
           else
             call dnaupd(idoarp,bmatarp,arna_ns,whicharp,arna_negv,
     $         tst_tol,resida,arna_nkrl,vbasea,arna_ls,iparp,ipntarp,
     $         workda,workla,nwlarp,infarp)
           endif  
#endif

           ! error check
           if (infarp.lt.0) then
               write(str,*) infarp
               call  mntr_abort(arna_id,
     $  'arn_naupd; inner prod. error with _naupd, info = '//trim(str))
           endif
           if (idoarp.ne.2) exit
         enddo
      endif                     ! idoarp.eq.2

      ! restart stepper
      if (idoarp.eq.-1.or.idoarp.eq.1) then
         call mntr_log(arna_id,lp_prd,'Restarting stepper')

!        move renormed data back to nekton
         if (arna_ifcomplex) then
           i = ipntarp(1)
           call copytoreal(vxp(1,1),vxp(1,2),workda(i),tst_nv)
           i = i+tst_nv
           call copytoreal(vyp(1,1),vyp(1,2),workda(i),tst_nv)
           i = i+tst_nv
           if (iff3d) then 
             call copytoreal(vzp(1,1),vzp(1,2),workda(i),tst_nv)
             i = i+tst_nv
           endif
           if (arna_ifpr) then 
             call copytoreal(prp(1,1),prp(1,2),workda(i),tst_np)
             i = i+tst_np
           endif
           if (ifheat) then
             call copytoreal(tp(1,1,1),tp(1,1,2),workda(i),tst_nt)
             i = i+tst_nt
           endif  
         else  
           ! velocity
           i = ipntarp(1)
           call copy(VXP,workda(i),tst_nv)
           i = i + tst_nv
           call copy(VYP,workda(i),tst_nv)
           i = i + tst_nv
           if (IF3D.or.iff3d) then 
             call copy(VZP,workda(i),tst_nv)
             i = i + tst_nv
           endif
!          Pressure 
           if (arna_ifpr) then
             call copy(prp,workda(i),tst_np)
             i = i + tst_np
           endif  
!          Temperature
           if (IFHEAT) then
             call copy(TP,workda(i),tst_nt)
             i = i + tst_nt
           endif  

!          make sure the velocity and temperature fields are continuous at
!          element faces and edges
         endif          ! arna_ifcomplex
         call tst_dssum
      endif                     ! idoarp.eq.-1.or.idoarp.eq.1

      return
      end
!=======================================================================




