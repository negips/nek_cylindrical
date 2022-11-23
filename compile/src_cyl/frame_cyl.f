!======================================================================
!     Description: Routines for Cylindrical Solver
!     Author: Prabal S. Negi
!     Subroutines to make the module compatible with Framework
!      
!====================================================================== 
!-----------------------------------------------------------------------
!> @brief Register FS_ALE module
!! @note This routine should be called in frame_usr_register

      subroutine frame_register_cyl()     
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'CYLINDRICAL'

      integer lpmid
      real ltim

      real dnekclock

!     timing
      ltim = dnekclock()

!     Check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid,cyl_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(cyl_name)//'] already registered')
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
      call mntr_mod_reg(cyl_id,lpmid,cyl_name,
     $          'Cylindrical Formulation')

!     Register timers
      call mntr_tmr_is_name_reg(lpmid,'FRM_TOT')
!     Total time
      call mntr_tmr_reg(cyl_tmr_tot_id,lpmid,cyl_id,
     $     'CYL_TOT','CYLINDRICAL module total time',.false.)
!     Initialisation time
      call mntr_tmr_reg(cyl_tmr_ini_id,lpmid,cyl_id,
     $     'CYL_INI','CYLINDRICAL initialisation time',.true.)
       
!     Register and set active section
      call rprm_sec_reg(cyl_sec_id,cyl_id,
     $     '_'//adjustl(cyl_name),
     $     'Runtime parameter section for CYLINDRICAL module')
      call rprm_sec_set_act(.true.,cyl_sec_id)

!     Register parameters
!     CYL_IFCYL
      call rprm_rp_reg(cyl_ifcyl_id,cyl_sec_id,'CYL_IFCYL',
     $     'Enable CYLINDRICAL? ',
     $     rpar_log,0,0.0,.false.,'  ')

!     CYL_LOG
      call rprm_rp_reg(cyl_log_id,cyl_sec_id,'CYL_LOG',
     $     'Log Level ',
     $     rpar_int,0,0.0,.false.,'  ')

!     CYL_OMEGA1
      call rprm_rp_reg(cyl_omega1_id,cyl_sec_id,'CYL_OMEGA1',
     $     'Inner Omega ',
     $     rpar_real,0,0.0,.false.,'  ')

!     CYL_OMEGA2
      call rprm_rp_reg(cyl_omega2_id,cyl_sec_id,'CYL_OMEGA2',
     $     'Outer Omega ',
     $     rpar_real,0,0.0,.false.,'  ')

!     timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(cyl_tmr_tot_id,1,ltim)

      return
      end subroutine frame_register_cyl

!---------------------------------------------------------------------- 

      subroutine frame_get_param_cyl()

      implicit none

      include 'SIZE'
      include 'INPUT'               ! param(59), initc
      include 'TSTEP'               ! time
      include 'FRAMELP'
      include 'CYLINDRICAL'

      ! local variables
      integer       itmp
      real          rtmp, ltim
      logical       ltmp
      character*20  ctmp
      character*2   str1, str2
      character*200 lstring

      ! functions
      real dnekclock

!     Timing
      ltim = dnekclock()

!     Get runtime parameters
!     cyl_ifcyl
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,cyl_ifcyl_id,rpar_log)
      cyl_ifcyl = ltmp

!     cyl_log
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,cyl_log_id,rpar_int)
      cyl_log = itmp

!     cyl_omega1
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,cyl_omega1_id,rpar_real)
      cyl_omega(1) = rtmp

!     cyl_omega2
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,cyl_omega2_id,rpar_real)
      cyl_omega(2) = rtmp

!     Timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(cyl_tmr_tot_id,1,ltim)

      return
      end subroutine frame_get_param_cyl 
!---------------------------------------------------------------------- 






