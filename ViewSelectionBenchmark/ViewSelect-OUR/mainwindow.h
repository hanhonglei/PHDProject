#pragma  once

#include "ui_mainwindow.h"
#include "glwidget.h"
#include "myrender.h"
#include "log.h"

class MainWindow : public QMainWindow, private Ui::MainWindow
{
	Q_OBJECT

public:
	MainWindow(QMainWindow *parent = 0);

private:
	MyRender *myrender;
	GLWidget *glw;

public slots:
	void on_btn_load_model_2d_clicked();
	void on_btn_compute_potential_2d_clicked();
	void on_btn_compute_skeleton_2d_clicked();
	void on_btn_load_model_3d_clicked();
	void on_btn_compute_potential_3d_clicked();
	void on_btn_compute_skeleton_3d_clicked();
	void on_btn_save_potential_3d_clicked();
	void on_btn_load_potential_3d_clicked();
	void on_rb_model_2d_toggled(bool checked);
	void on_rb_primary_skeleton_2d_toggled(bool checked);
	void on_rb_final_skeleton_2d_toggled(bool checked);
	void on_rb_whole_skeleton_2d_toggled(bool checked);
	void on_rb_model_3d_toggled(bool checked);
	void on_rb_model_voxels_3d_toggled(bool checked);
	void on_rb_segmentation_3d_toggled(bool checked);
	void on_rb_primary_skeleton_3d_toggled(bool checked);
	void on_rb_final_skeleton_3d_toggled(bool checked);
	void on_rb_whole_skeleton_3d_toggled(bool checked);
	void on_rb_view_dependent_curvature_toggled(bool checked);
	void on_cb_select_parts_by_mouse_stateChanged(int state);
	void on_btn_build_segmentation_clicked();
	void on_btn_compute_entropy_clicked();
	void on_btn_compute_entropy_revised_clicked();
	void on_btn_viewpoint_selection_clicked();
	void on_btn_move_clicked();
	void on_btn_testing_a_round_clicked();
	void on_btn_show_A_clicked();
	void on_btn_show_B1_clicked();
	void on_btn_show_B2_clicked();
	void on_btn_show_B3_clicked();
	void on_rb_viewing_sphere_toggled(bool checked);
	void on_btn_compute_mesh_saliency_3d_clicked();
	void on_rb_mesh_saliency_3d_toggled(bool checked);
	void on_sb_primary_skeleton_2d_valueChanged(int i);
	void on_sb_final_skeleton_2d_valueChanged(int i);
	void on_sb_primary_skeleton_3d_valueChanged(int i);
	void on_sb_final_skeleton_3d_valueChanged(int i);
	void on_btn_load_segmentation_3d_clicked();
	void on_btn_save_result_clicked();
	void on_btn_save_viewpoint_clicked();
};