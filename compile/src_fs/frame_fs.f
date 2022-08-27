!======================================================================
!     Description: Routines for free surface simulation
!     Author: Prabal S. Negi
!     Subroutines to make the module compatible with Framework
!      
!====================================================================== 
!-----------------------------------------------------------------------
!> @brief Register FS_ALE module
!! @note This routine should be called in frame_usr_register

      subroutine frame_register_fs()     
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'FRAMELP'
      include 'FS_ALE'

      integer lpmid
      real ltim

      real dnekclock

!     timing
      ltim = dnekclock()

!     Check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid,fs_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(fs_name)//'] already registered')
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
      call mntr_mod_reg(fs_id,lpmid,fs_name,
     $          'Free Surface ALE.')

!     Register timers
      call mntr_tmr_is_name_reg(lpmid,'FRM_TOT')
!     Total time
      call mntr_tmr_reg(fs_tmr_tot_id,lpmid,fs_id,
     $     'FS_TOT','FS module total time',.false.)
!     Initialisation time
      call mntr_tmr_reg(fs_tmr_ini_id,lpmid,fs_id,
     $     'FS_INI','FS initialisation time',.true.)
       
!     Register and set active section
      call rprm_sec_reg(fs_sec_id,fs_id,
     $     '_'//adjustl(fs_name),
     $     'Runtime parameter section for FS module')
      call rprm_sec_set_act(.true.,fs_sec_id)

!     Register parameters
!     IFFS
      call rprm_rp_reg(fs_iffs_id,fs_sec_id,'FS_IFFS',
     $     'Enable FS? ',
     $     rpar_log,0,0.0,.false.,'  ')

!     Register parameters
!     FS_IFGSM
      call rprm_rp_reg(fs_ifgsm_id,fs_sec_id,'FS_IFGSM',
     $     'Global Smoothening of FS? ',
     $     rpar_log,0,0.0,.false.,'  ')

!     FS_IFTC
      call rprm_rp_reg(fs_iftc_id,fs_sec_id,'FS_IFTC',
     $     'Tangential Correction of FS? ',
     $     rpar_log,0,0.0,.false.,'  ')

!     FS_IFGH
      call rprm_rp_reg(fs_ifgh_id,fs_sec_id,'FS_IFGH',
     $     'Gordan Hall Correction? ',
     $     rpar_log,0,0.0,.false.,'  ')

!     FS_IFFIL
      call rprm_rp_reg(fs_iffil_id,fs_sec_id,'FS_IFFIL',
     $     'Filter Global Interpolation? ',
     $     rpar_log,0,0.0,.false.,'  ')

!     FS_IFGRID
      call rprm_rp_reg(fs_ifgrid_id,fs_sec_id,'FS_IFGRID',
     $     'Global Approx. of Grid? ',
     $     rpar_log,0,0.0,.false.,'  ')
     
!     FS_OFST
      call rprm_rp_reg(fs_ofst_id,fs_sec_id,'FS_OFST',
     $     'Damping Offset ',
     $     rpar_real,0,0.0,.false.,'  ')

!     FS_LOG
      call rprm_rp_reg(fs_log_id,fs_sec_id,'FS_LOG',
     $     'Log Level ',
     $     rpar_int,0,0.0,.false.,'  ')

      ! timing
      ltim = dnekclock() - ltim
      call mntr_tmr_add(fs_tmr_tot_id,1,ltim)

      return
      end subroutine frame_register_fs

!---------------------------------------------------------------------- 

      subroutine frame_get_param_fs()

      implicit none

      include 'SIZE'
      include 'INPUT'               ! param(59), initc
      include 'TSTEP'               ! time
      include 'FRAMELP'
      include 'FS_ALE'

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

      ! get runtime parameters
!     iffs
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,fs_iffs_id,rpar_log)
      fs_iffs = ltmp
!     fs_ofst
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,fs_ofst_id,rpar_real)
      fs_ofst = rtmp
!     ifgsm
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,fs_ifgsm_id,rpar_log)
      fs_ifgsm = ltmp
!     iftc
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,fs_iftc_id,rpar_log)
      fs_iftc = ltmp
!     ifgh
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,fs_ifgh_id,rpar_log)
      fs_ifgh = ltmp
!     iffil
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,fs_iffil_id,rpar_log)
      fs_iffil = ltmp
!     ifgrid
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,fs_ifgrid_id,rpar_log)
      fs_ifgrid = ltmp

!     log_fs
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,fs_log_id,rpar_int)
      fs_log = itmp

!!     fs_slipl
!      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,fs_slipl_id,rpar_real)
!      fs_slipl = rtmp
!!     fs_blendl
!      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,fs_blendl_id,rpar_real)
!      fs_blendl = rtmp

      return
      end subroutine frame_get_param_fs 
!---------------------------------------------------------------------- 

