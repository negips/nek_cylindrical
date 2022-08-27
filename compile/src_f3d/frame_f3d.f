!======================================================================
!     Routines for introducing third (Fourier) component in a 2d simulation
!     Author: Prabal S. Negi
!     Subroutines to make the module compatible with Framework
!      
!====================================================================== 
!-----------------------------------------------------------------------
!> @brief Register F3D module
!! @note This routine should be called in frame_usr_register

      subroutine frame_register_f3d()     
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'F3D'

      integer lpmid
      real ltim

      real dnekclock

!     timing
      ltim = dnekclock()

!     Check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid,f3d_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(f3d_name)//'] already registered')
         return
      endif

!     Find parent module
      call mntr_mod_is_name_reg(lpmid,'FRAME')
      if (lpmid.le.0) then
         lpmid = 1
         call mntr_abort(lpmid,
     $        'Parent module ['//'FRAME'//'] not registered')
      endif

!     Register module
      call mntr_mod_reg(f3d_id,lpmid,f3d_name,
     $          'Fourier 3D Solve.')

!     Register timers
      call mntr_tmr_is_name_reg(lpmid,'FRM_TOT')
!     Total time
      call mntr_tmr_reg(f3d_tmr_tot_id,lpmid,f3d_id,
     $     'F3D_TOT','F3D module total time',.false.)
!     Initialisation time
      call mntr_tmr_reg(f3d_tmr_ini_id,lpmid,f3d_id,
     $     'F3D_INI','F3D initialisation time',.true.)
       
!     Register and set active section
      call rprm_sec_reg(f3d_sec_id,f3d_id,
     $     '_'//adjustl(f3d_name),
     $     'Runtime parameter section for F3D module')
      call rprm_sec_set_act(.true.,f3d_sec_id)

!     Register parameters
!     IFF3D
      call rprm_rp_reg(f3d_iff3d_id,f3d_sec_id,'IFF3D',
     $     'Enable F3D? ',
     $     rpar_log,0,0.0,.false.,'  ')
!     IFCYL_F3D
      call rprm_rp_reg(f3d_ifcyl_id,f3d_sec_id,'IFCYL_F3D',
     $     'Cylindrical Coordinates? ',
     $     rpar_log,0,0.0,.false.,'  ')
     
!     K_F3D
      call rprm_rp_reg(f3d_k_id,f3d_sec_id,'K_F3D',
     $     'F3D wavenumber ',
     $     rpar_real,0,0.0,.false.,'  ')

!     F3D_SLIPL
      call rprm_rp_reg(f3d_slipl_id,f3d_sec_id,'SLIPL_F3D',
     $     'Slip Length ',
     $     rpar_real,0,0.0,.false.,'  ')

!     F3D_BLENDL
      call rprm_rp_reg(f3d_blendl_id,f3d_sec_id,'BLENDL_F3D',
     $     'Blending Length ',
     $     rpar_real,0,1.0,.false.,'  ')

!     LOG_F3D
      call rprm_rp_reg(f3d_log_id,f3d_sec_id,'LOG_F3D',
     $     'Log Level ',
     $     rpar_int,1,0.0,.false.,'  ')

      ! set initialisation flag
!      otd_ifinit=.false.

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(f3d_tmr_tot_id,1,ltim)

      return
      end subroutine frame_register_f3d

!---------------------------------------------------------------------- 

      subroutine frame_get_param_f3d()

      implicit none

      include 'SIZE'
      include 'INPUT'               ! param(59), initc
      include 'TSTEP'               ! time
      include 'FRAMELP'
      include 'F3D'

      ! local variables
      integer       itmp
      real          rtmp, ltim
      logical       ltmp
      character*20  ctmp
      character*2   str1, str2
      character*200 lstring

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

!      ! check if the module was already initialised
!      if (otd_ifinit) then
!         call mntr_warn(otd_id,
!     $        'module ['//trim(otd_name)//'] already initiaised.')
!         return
!      endif

      ! get runtime parameters
!     iff3d
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,f3d_iff3d_id,rpar_log)
      iff3d = ltmp
!     ifcyl_f3d
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,f3d_ifcyl_id,rpar_log)
      ifcyl_f3d = ltmp
!     k_f3d
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,f3d_k_id,rpar_real)
      k_f3d = rtmp
!     slipl_f3d
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,f3d_slipl_id,rpar_real)
      slipl_f3d = rtmp
!     blendl_f3d
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,f3d_blendl_id,rpar_real)
      blendl_f3d = rtmp
!     log_f3d
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,f3d_log_id,rpar_int)
      log_f3d = itmp

      return
      end subroutine frame_get_param_f3d        
!---------------------------------------------------------------------- 


!-----------------------------------------------------------------------











