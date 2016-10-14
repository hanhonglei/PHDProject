#include <QtGui>
#include <QString>

#include "mainwindow.h"

MainWindow::MainWindow(QMainWindow *parent /* = 0 */):QMainWindow(parent)
{
	setupUi(this);
	//setWindowFlags(Qt::Window | Qt::WindowMinimizeButtonHint | Qt::MSWindowsFixedSizeDialogHint);

	// setup log system

	_LOG = edit_log ;
	// setup OpenGL widget and initiate

	myrender = new MyRender();
	glw = new GLWidget(myrender, widget_gl_container);
	glw->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);

	QGridLayout *layout_gl_container = new QGridLayout();
	layout_gl_container->setMargin(0);
	layout_gl_container->setSpacing(0);
	layout_gl_container->addWidget(glw,0,0);
	widget_gl_container->setLayout(layout_gl_container);
}

void MainWindow::on_btn_load_model_2d_clicked()
{
	static QString filename;
	filename = QFileDialog::getOpenFileName(this);
	if (!filename.isEmpty())
	{
		myrender->load_model_2d(filename.toStdWString());
	}
	glw->updateGL();
}

void MainWindow::on_btn_compute_potential_2d_clicked()
{
	LOG("Computing 2D potential field, it maybe need a few minutes...");
	//std::cout<<std::endl; 
	// for immediate showing the msg
	myrender->compute_potential_2d();
	LOG("Computing 2D potential field, Done!");
	glw->updateGL();
}

void MainWindow::on_btn_compute_skeleton_2d_clicked()
{
	myrender->set_WOS_2D(sb_weight_of_skeleton_2d->value());

	myrender->extract_primary_skeleton_2d();
	myrender->build_skeleton_tree_2d();
	myrender->extract_final_skeleton_2d();
	myrender->segmentation_by_skeleton_2d();

	sb_primary_skeleton_2d->setValue(0);
	sb_primary_skeleton_2d->setRange(0, myrender->primary_skls2d.size()-1);
	myrender->set_showing_primary_skeleton_number_2d(0);

	sb_final_skeleton_2d->setValue(0);
	sb_final_skeleton_2d->setRange(0, myrender->final_skls2d.size()-1);
	myrender->set_showing_final_skeleton_number_2d(0);

	LOG("Extract 2D skeleton, Done!");
	glw->updateGL();
}

void MainWindow::on_btn_load_model_3d_clicked()
{
	static QString filename;
	filename = QFileDialog::getOpenFileName(
		this,"Open Model File", filename, ("Model Files(*.ply *.obj *.off)"));
	if (!filename.isEmpty())
	{
		freopen("out.txt","w",stdout);
		std::wstring ss=L"";
		for (int i=0;i<filename.size();i++) ss+=filename[i].toAscii();
		//cout<<ss;
		//printf("%s",filename.toStdString().c_str());
		//LOG("file : " + QString(filename.toStdString()));
		myrender->load_model_3d(ss);
	}
	glw->updateGL();
}

void MainWindow::on_btn_load_segmentation_3d_clicked()
{
	static QString filename;
	filename = QFileDialog::getOpenFileName(
		this,"Load Segmentation", filename, ("Segmentation Files(*.seg)"));
	if (!filename.isEmpty())
	{
		myrender->load_segmentation(filename.toStdWString());
	}
	glw->updateGL();
}

void MainWindow::on_btn_compute_potential_3d_clicked()
{
	LOG("Computing 3D potential field, it maybe need half a hour or more time...");
   // std::cout<<std::endl;
	myrender->compute_potential_3d();
	LOG("Computing 3D potential field, Done!");
	glw->updateGL();
}

void MainWindow::on_btn_compute_skeleton_3d_clicked()
{
	myrender->set_WOS_3D(sb_weight_of_skeleton_3d->value());

	myrender->extract_primary_skeleton_3d();
	myrender->build_skeleton_tree_3d();
	myrender->extract_final_skeleton_3d();
	myrender->segmentation_by_skeleton_3d();
	//myrender->map_to_mesh_segmentation();

	sb_primary_skeleton_3d->setValue(0);
	sb_primary_skeleton_3d->setRange(0, myrender->primary_skls3d.size()-1);
	myrender->set_showing_primary_skeleton_number_3d(0);

	sb_final_skeleton_3d->setValue(0);
	sb_final_skeleton_3d->setRange(0, myrender->final_skls3d.size()-1);
	myrender->set_showing_final_skeleton_number_3d(0);
	
	LOG("Extract 3D skeleton, Done!");
	glw->updateGL();
}

void MainWindow::on_btn_load_potential_3d_clicked()
{
	static QString filename;
	filename = QFileDialog::getOpenFileName(
		this,"Load Potential Field", filename, ("Potential field (*.pf)"));

	if (!filename.isEmpty())
	{
		myrender->load_potential_3d(filename.toStdWString());
	}
	glw->updateGL();
}

void MainWindow::on_btn_save_potential_3d_clicked()
{
	static QString filename;
	filename = QFileDialog::getSaveFileName(
		this,"Save Potential Field", filename, ("Potential field (*.pf)"));
	if (!filename.isEmpty())
	{
		myrender->save_potential_3d(filename.toStdWString());
	}
	glw->updateGL();
}

void MainWindow::on_rb_primary_skeleton_2d_toggled(bool checked)
{
	if (checked)
	{
		myrender->showing_type = MyRender::SHOWING_PRIMARY_SKELETON;
		myrender->showing_primary_skeleton_number_2d = sb_primary_skeleton_2d->value();
		glw->updateGL();
	}
}

void MainWindow::on_rb_final_skeleton_2d_toggled(bool checked)
{
	if (checked)
	{
		myrender->showing_type = MyRender::SHOWING_FINAL_SKELETON;
		myrender->showing_final_skeleton_number_2d = sb_final_skeleton_2d->value();
		glw->updateGL();
	}
}

void MainWindow::on_rb_whole_skeleton_2d_toggled(bool checked)
{
	if (checked)
	{
		myrender->showing_type = MyRender::SHOWING_WHOLE_SKELETON;
		glw->updateGL();
	}
}

void MainWindow::on_rb_model_2d_toggled(bool checked)
{
	if (checked)
	{
		myrender->showing_type = MyRender::SHOWING_MODEL;
		glw->updateGL();
	}
}

void MainWindow::on_rb_primary_skeleton_3d_toggled(bool checked)
{
	if (checked)
	{
		myrender->showing_type = MyRender::SHOWING_PRIMARY_SKELETON;
		myrender->showing_primary_skeleton_number_3d = sb_primary_skeleton_3d->value();
		glw->updateGL();
	}
}

void MainWindow::on_rb_final_skeleton_3d_toggled(bool checked)
{
	if (checked)
	{
		myrender->showing_type = MyRender::SHOWING_FINAL_SKELETON;
		myrender->showing_final_skeleton_number_3d = sb_final_skeleton_3d->value();
		glw->updateGL();
	}
}

void MainWindow::on_rb_whole_skeleton_3d_toggled(bool checked)
{
	if (checked)
	{
		myrender->showing_type = MyRender::SHOWING_WHOLE_SKELETON;
		glw->updateGL();
	}
}

void MainWindow::on_rb_view_dependent_curvature_toggled(bool checked)
{
	if(checked)
	{
		myrender->showing_type = MyRender::SHOWING_VIEW_DEPENDENT_CURVATURE;
		glw->updateGL();
	}
}

void MainWindow::on_rb_viewing_sphere_toggled(bool checked)
{

	if(checked)
	{
		myrender->showing_type = MyRender::SHOWING_VIEW_SPHERE_MAP;
		glw->updateGL();
	}
}

void MainWindow::on_btn_build_segmentation_clicked()
{
	std::vector<int> parts;
	std::string strbuf = le_parts->text().toStdString();
	if (strbuf.size() == 0) // restore to whole model
	{
		myrender->build_segmentation(parts);
		glw->updateGL();
		return;
	}
	std::stringstream ssbuf(strbuf);
	int i=0;
	while (!ssbuf.eof())
	{
		ssbuf>>i;
		parts.push_back(i++);
	}
	myrender->build_segmentation(parts);
	glw->updateGL();
}

void MainWindow::on_cb_select_parts_by_mouse_stateChanged(int state)
{
	if (state == Qt::Checked)
	{
		myrender->change_mouse_mode(1);
		glw->selecting_parts_by_mouse = true;
	}
	else
	{
		myrender->change_mouse_mode(0);
		glw->selecting_parts_by_mouse = false;
	}
	glw->updateGL();
}

void MainWindow::on_btn_compute_entropy_clicked()
{
	//myrender->compute_shannon_entropy_II();
	myrender->compute_entropy();
}

void MainWindow::on_btn_compute_entropy_revised_clicked()
{
	myrender->compute_revised_entropy();
}

void MainWindow::on_btn_move_clicked()
{
	while(1)
	{
		glw->mouseMove();
	}
}

void MainWindow::on_btn_viewpoint_selection_clicked()
{
#if USE_PROJECTING_COVERED_VIEWPLANE
	myrender->viewpoint_selection2();
#else
	myrender->sample_using_revised_entropy_II();
#endif
	glw->updateGL();
}

void MainWindow::on_rb_model_3d_toggled(bool checked)
{
	if (checked)
	{
		myrender->showing_type = MyRender::SHOWING_MODEL;
		glw->updateGL();
	}
}

void MainWindow::on_rb_model_voxels_3d_toggled(bool checked)
{
	if (checked)
	{
		myrender->showing_type = MyRender::SHOWING_MODEL_VOXELS;
		glw->updateGL();
	}
}

void MainWindow::on_rb_segmentation_3d_toggled(bool checked)
{
	if (checked)
	{
		myrender->showing_type = MyRender::SHOWING_SEGMENTATION;
		glw->updateGL();
	}
}

void MainWindow::on_rb_mesh_saliency_3d_toggled(bool checked)
{
	if (checked)
	{
		myrender->showing_type = MyRender::SHOWING_MESH_SALIENCY;
		glw->updateGL();
	}
}

void MainWindow::on_btn_compute_mesh_saliency_3d_clicked()
{
	myrender->compute_mesh_saliency_3d();
	glw->updateGL();
}

void  MainWindow::on_sb_primary_skeleton_2d_valueChanged(int i)
{
	myrender->set_showing_primary_skeleton_number_2d(i);
	glw->updateGL();
}

void  MainWindow::on_sb_final_skeleton_2d_valueChanged(int i)
{
	myrender->set_showing_final_skeleton_number_2d(i);
	glw->updateGL();
}

void  MainWindow::on_sb_primary_skeleton_3d_valueChanged(int i)
{
	myrender->set_showing_primary_skeleton_number_3d(i);
	glw->updateGL();
}

void  MainWindow::on_sb_final_skeleton_3d_valueChanged(int i)
{
	myrender->set_showing_final_skeleton_number_3d(i);
	glw->updateGL();
}

void MainWindow::on_btn_save_result_clicked()
{
	static QString filename;
	filename = QFileDialog::getSaveFileName(
		this,"Save The Viewpoint Selection Result", filename, ("Viewpoint selection result (*.ppm)"));
	if (filename.isEmpty()) return;

	// Find first non-used filename
	FILE *f;
	string fn = "";
	for (int i=0;i<filename.size();i++) fn+=filename[i].toAscii();
	f = fopen(fn.c_str(), "wb");
	// Read pixels
	int width = glw->size().width(), height = glw->size().height(); 
	char *buf = new char[width*height*3];
//	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	//ust change the packing to ensure no overruns!
	glPixelStorei(GL_PACK_ALIGNMENT, 1);

	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, buf);

	// Flip top-to-bottom
	for (int i = 0; i < height/2; i++) {
		char *row1 = buf + 3 * width * i;
		char *row2 = buf + 3 * width * (height - 1 - i);
		for (int j = 0; j < 3 * width; j++)
			swap(row1[j], row2[j]);
	}

	// Write out file
	fprintf(f, "P6\n%d %d\n255\n", width, height);
	fwrite(buf, width*height*3, 1, f);
	fclose(f);
	delete [] buf;
	// restore default alignment
	glPixelStorei(GL_PACK_ALIGNMENT, 4);

	LOG("The result has been saved.\n");
}

void MainWindow::on_btn_save_viewpoint_clicked()
{
	myrender->save_viewpoint();
}

void MainWindow::on_btn_testing_a_round_clicked()
{
	myrender->testing_a_round();
}

void MainWindow::on_btn_show_A_clicked()
{
	myrender->show_A();
	glw->updateGL();
}

void MainWindow::on_btn_show_B1_clicked()
{
	myrender->show_B1();
	glw->updateGL();
}

void MainWindow::on_btn_show_B2_clicked()
{
	myrender->show_B2();
	glw->updateGL();
}

void MainWindow::on_btn_show_B3_clicked()
{
	myrender->show_B3();
	glw->updateGL();
}