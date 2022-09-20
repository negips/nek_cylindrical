
SRC_INC     := ./src_inc
SRC_USR     := ./src_usr
SRC_FRAME   := ./frame
SRC_ARN     := ./src_arn
SRC_MATF    := ./src_matf

SRC_F3D     := ./src_f3d
SRC_FS      := ./src_fs
SRC_LSM     := ./src_lsm
SRC_STT     := ./stat
SRC_MAP     := ./map2D
SRC_CHP     := ./checkpoint

$(OBJDIR)/frame_usr.o           :$(SRC_FRAME)/frame_usr.f;              $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/frame.o               :$(SRC_FRAME)/frame.f;                  $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/mntrlog_block.o       :$(SRC_FRAME)/mntrlog_block.f;          $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/mntrlog.o             :$(SRC_FRAME)/mntrlog.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/mntrtmr_block.o       :$(SRC_FRAME)/mntrtmr_block.f;          $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/mntrtmr.o             :$(SRC_FRAME)/mntrtmr.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/rprm_block.o          :$(SRC_FRAME)/rprm_block.f;             $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/rprm.o                :$(SRC_FRAME)/rprm.f;                   $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/io_tools.o            :$(SRC_FRAME)/io_tools.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/io_tools_block.o      :$(SRC_FRAME)/io_tools_block.f;         $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/math_tools.o          :$(SRC_FRAME)/math_tools.f;             $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/solve_f3d.o           :$(SRC_F3D)/solve_f3d.f ;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/solvep_f3d.o          :$(SRC_F3D)/solvep_f3d.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/weak_lapl_cyl.o       :$(SRC_F3D)/weak_lapl_cyl.f;            $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/solver_f3d.o          :$(SRC_F3D)/solver_f3d.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/frame_f3d.o           :$(SRC_F3D)/frame_f3d.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/tstepper.o            :$(SRC_ARN)/tstepper.f;                 $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/conht_tools.o         :$(SRC_ARN)/conht_tools.f;              $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/arn_arp.o             :$(SRC_ARN)/arn_arp.f;                  $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/arn_arp_io.o          :$(SRC_ARN)/arn_arp_io.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/arn_usr.o             :$(SRC_ARN)/arn_usr.f;                  $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/arn_cmath.o           :$(SRC_ARN)/arn_cmath.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/chkpoint.o            :$(SRC_CHP)/chkpoint.f;                 $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/chkpt_mstp.o          :$(SRC_CHP)/chkpt_mstp.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/fs_ale.o              :$(SRC_FS)/fs_ale.f;                    $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/fs_sort.o             :$(SRC_FS)/fs_sort.f;                   $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/fs_global_surf.o      :$(SRC_FS)/fs_global_surf.f;            $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/frame_fs.o            :$(SRC_FS)/frame_fs.f;                  $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/localfilter.o         :$(SRC_USR)/localfilter.f;              $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/wrapper_lapack.o      :$(SRC_MATF)/wrapper_lapack.f;          $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/matf_krylov.o         :$(SRC_MATF)/matf_krylov.f;             $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/matrix_functions.o    :$(SRC_MATF)/matrix_functions.f;        $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/lsm.o                 :$(SRC_LSM)/lsm.f;                      $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/lsm_newton.o          :$(SRC_LSM)/lsm_newton.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/lsm_heat.o            :$(SRC_LSM)/lsm_heat.f;                 $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/clsm.o                :$(SRC_LSM)/clsm.f;                     $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/clsm_corr.o           :$(SRC_LSM)/clsm_corr.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/lsm_zhoge.o           :$(SRC_LSM)/lsm_zhoge.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/dealias.o             :$(SRC_USR)/dealias.f;                  $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/stat.o                :$(SRC_STT)/stat.f;                     $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/stat_IO.o             :$(SRC_STT)/stat_IO.f;                  $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/stat_usr.o            :$(SRC_STT)/stat_usr.f;                 $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/map2D.o               :$(SRC_MAP)/map2D.f;                    $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/map2D_usr.o           :$(SRC_MAP)/map2D_usr.f;                $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/fast_reset.o          :$(SRC_USR)/fast_reset.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/h1mg_reset.o          :$(SRC_USR)/h1mg_reset.f;               $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/laplace_testing.o     :$(SRC_USR)/laplace_testing.f;          $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/crs_solve_new.o       :$(SRC_USR)/crs_solve_new.f;            $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/orthonew.o            :$(SRC_USR)/orthonew.f;                 $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/rotate_cyclic.o       :$(SRC_USR)/rotate_cyclic.f;            $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/gmres_solver.o        :$(SRC_USR)/gmres_solver.f;             $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/mesh2_exchange.o      :$(SRC_USR)/mesh2_exchange.f;           $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/subs_ifpgll.o         :$(SRC_USR)/subs_ifpgll.f;              $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@
$(OBJDIR)/null_space_jd.o       :$(SRC_USR)/null_space_jd.f;            $(F77) -c $(FL2) -I$(SRC_INC)/ $< -o $@



