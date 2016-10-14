#pragma  once

//#include "QStreamRedirector.h"
#include "ui_mainwindow.h"
#include "myrender.h"
#include "glwidget.h"

class MainWindow : public QMainWindow, private Ui::MainWindow
{
    Q_OBJECT

public:
    MainWindow(QMainWindow *parent = 0);

private:
    MyRender *myrender;
    GLWidget *glw;
   // QStreamRedirector *stream_redirector;

public slots:
    void on_action_load_model_triggered();
    void on_action_exit_triggered();

    void on_rbtn_original_image_toggled(bool checked);
    void on_rbtn_depth_buffer_image_toggled(bool checked);
    void on_rbtn_radial_curvature_image_toggled(bool checked);
    void on_rbtn_model_space_curvature_image_toggled(bool checked);
    void on_rbtn_view_dependent_curvature_image_toggled(bool checked);
    void on_rbtn_mesh_saliency_toggled(bool checked);

    void on_cmb_adaptive_box_size_currentIndexChanged(int index);
    void on_cmb_number_of_histogram_intervals_currentIndexChanged(int index);

    void on_btn_current_shannon_entropy_clicked();
    void on_btn_current_shannon_entropy_II_clicked();
    void on_btn_current_entropy_modified_II_clicked();

    void on_btn_sampling_using_shannon_entropy_clicked();
    void on_btn_sampling_using_shannon_entropy_II_clicked();
    void on_btn_sampling_using_revised_entropy_clicked();
    void on_btn_sampling_using_revised_entropy_II_clicked();

    void on_btn_current_mesh_saliency_clicked();
    void on_btn_sampling_using_mesh_saliency_clicked();

    void on_btn_testing_a_round_clicked();
    void on_btn_save_result_clicked();
};