#include <iostream>
#include <QtGui>
#include <QtCore/QString>

#include "mainwindow.h"
#include "cmath"

MainWindow::MainWindow(QMainWindow *parent /* = 0 */):QMainWindow(parent)
{
    setupUi(this);
    //setWindowFlags(Qt::Window | Qt::WindowMinimizeButtonHint | Qt::MSWindowsFixedSizeDialogHint);

    // setup OpenGL widget and Render
    myrender = new MyRender();
    
    glw = new GLWidget(myrender, gl_container);
    glw->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

    QGridLayout *gl_container_layout = new QGridLayout();
    gl_container_layout->setMargin(0);
    gl_container_layout->setSpacing(0);
    gl_container_layout->addWidget(glw,0,0);
    gl_container->setLayout(gl_container_layout);

    // setup log system
    //stream_redirector = new QStreamRedirector(edit_log);
}

void MainWindow::on_action_load_model_triggered()
{
    QString filename = QFileDialog::getOpenFileName(this);
    if (!filename.isEmpty())
        myrender->load_model(filename.toStdWString());
    glw->updateGL();
}

void MainWindow::on_action_exit_triggered()
{
    qApp->exit(0);
}

void MainWindow::on_rbtn_original_image_toggled(bool checked)
{
    if (checked)
    {
        myrender->set_showing_image_type(MyRender::original_image);
        glw->updateGL();
    }
}

void MainWindow::on_rbtn_depth_buffer_image_toggled(bool checked)
{
    if (checked)
    {
        myrender->set_showing_image_type(MyRender::depth_buffer_image);
        glw->updateGL();
    }
}

void MainWindow::on_rbtn_radial_curvature_image_toggled(bool checked)
{
    if (checked)
    {
        myrender->set_showing_image_type(MyRender::radial_curvature_image);
        glw->updateGL();
    }
}

void MainWindow::on_rbtn_model_space_curvature_image_toggled(bool checked)
{
    if (checked)
    {
        myrender->set_showing_image_type(MyRender::model_space_curvature_image);
        glw->updateGL();
    }
}

void MainWindow::on_rbtn_view_dependent_curvature_image_toggled(bool checked)
{
    if (checked)
    {
        myrender->set_showing_image_type(MyRender::view_dependent_curvature_image);
        glw->updateGL();
    }
}

void MainWindow::on_rbtn_mesh_saliency_toggled(bool checked)
{
    if (checked)
    {
        myrender->set_showing_image_type(MyRender::mesh_saliency);
        glw->updateGL();
    }
}

void MainWindow::on_cmb_adaptive_box_size_currentIndexChanged(int index)
{
    int size = (int)std::pow(2.0,index);
    std::cout
        <<"Current adaptive box size of entropy computing is changed to :\t"
        <<size
        <<"\n";
    myrender->set_adaptive_box_size(size);
}

void MainWindow::on_cmb_number_of_histogram_intervals_currentIndexChanged(int index)
{
    int size = (int)std::pow(2.0,index+1);
    std::cout
        <<"Current number of histogram intervals is changed to :\t"
        <<size
        <<"\n";
    myrender->set_number_of_histogram_intervals(size);
}

void MainWindow::on_btn_current_shannon_entropy_clicked()
{
    myrender->compute_shannon_entropy();
}

void MainWindow::on_btn_current_shannon_entropy_II_clicked()
{
    myrender->compute_shannon_entropy_II();
}

void MainWindow::on_btn_current_entropy_modified_II_clicked()
{
    myrender->compute_revised_entropy_II();
}

void MainWindow::on_btn_sampling_using_shannon_entropy_clicked()
{
    myrender->sample_using_shannon_entropy();
    glw->updateGL();
}

void MainWindow::on_btn_sampling_using_shannon_entropy_II_clicked()
{
    myrender->sample_using_shannon_entropy_II();
    glw->updateGL();
}

void MainWindow::on_btn_sampling_using_revised_entropy_clicked()
{
    myrender->sample_using_revised_entropy();
    glw->updateGL();
}

void MainWindow::on_btn_sampling_using_revised_entropy_II_clicked()
{
    myrender->sample_using_revised_entropy_II();
    glw->updateGL();
}

void MainWindow::on_btn_current_mesh_saliency_clicked()
{
    myrender->compute_mesh_saliency();
}

void MainWindow::on_btn_sampling_using_mesh_saliency_clicked()
{
    myrender->sample_using_mesh_saliency();
    glw->updateGL();
}

void MainWindow::on_btn_testing_a_round_clicked()
{
	//LOG("testing");
    myrender->testing_a_round();
    glw->updateGL();
}

void MainWindow::on_btn_save_result_clicked()
{
    myrender->save_result();
    glw->update();
}