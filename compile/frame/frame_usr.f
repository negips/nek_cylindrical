!> @brief Register user specified modules
!====================================================================== 
      subroutine frame_usr_register
      implicit none

      include 'SIZE'
      include 'FRAMELP'

!     register modules
      call io_register
      call chkpt_register
!      call frame_register_f3d
      call frame_register_fs
      call map2D_register
      call stat_register
!      call frame_register_f3d
!      call tst_register

      return
      end subroutine
!======================================================================
!> @brief Initialise user specified modules
      subroutine frame_usr_init
      implicit none

      include 'SIZE'
      include 'FRAMELP'

!     initialise modules
      call chkpt_init
!      call map2D_init
!      call stat_init
!      call frame_get_param_f3d
      call frame_get_param_fs
!      call tst_init

      return
      end subroutine
!======================================================================
!> @brief Finalise user specified modules
      subroutine frame_usr_end
      implicit none

      include 'SIZE'
      include 'FRAMELP'

      
      return
      end subroutine

!-----------------------------------------------------------------------

