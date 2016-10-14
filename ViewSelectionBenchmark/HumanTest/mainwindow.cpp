#include <iostream>
#include <QtGui>
#include <QString>
#include <QAction>
#include <QFile>
#include <QTime>

#include "mainwindow.h"

QString LOC2UCS(const char *locstr)
{  
    QTextCodec *gbk_codec = QTextCodec::codecForLocale();
    return gbk_codec->toUnicode(locstr);
}
#define L2U LOC2UCS

MainWindow::MainWindow(QMainWindow *parent /* = 0 */):QMainWindow(parent)
{
    setupUi(this);
    //icon.addPixmap(QPixmap(L2U(":/HumanTest.ico")), QIcon::Normal, QIcon::Off);
    //setWindowIcon(icon);
    setWindowFlags(Qt::Window | Qt::WindowMinimizeButtonHint | Qt::MSWindowsFixedSizeDialogHint);

    // setup OpenGL widget and Render
    myrender = new MyRender();
    
    glw = new GLWidget(myrender, gl_container);
    //glw->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    QSizePolicy sp;
    sp.setHeightForWidth(true);
    glw->setSizePolicy(sp);

    QGridLayout *gl_container_layout = new QGridLayout();
    gl_container_layout->setMargin(0);
    gl_container_layout->setSpacing(0);
    gl_container_layout->addWidget(glw,0,0);
    gl_container->setLayout(gl_container_layout);

    // setup log system
    //stream_redirector = new QStreamRedirector(edit_log);

    // get model list
    QFile ml(":/modellist"); 
    ml.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream mltx(&ml);
    while(!mltx.atEnd())
    {
        QString x = mltx.readLine();
        models.push_back(x.toStdWString());
    }
    ml.close();

    // prepare result views
    // if view.dat doesn't exists or in wrong format, init all to [0,0,0,1], or load previous data
    result_quats.resize(models.size()*4,0.0);
    if (QFile::exists("./view.dat"))
    {
        QFile viewdat("./view.dat");
        viewdat.open(QIODevice::ReadOnly | QIODevice::Text);
        QTextStream viewdattx(&viewdat);
        // parse and load view data
        for (unsigned int i=0; i<models.size(); i++)
        {
            viewdattx
                >>result_quats[i*4+0]
                >>result_quats[i*4+1]
                >>result_quats[i*4+2]
                >>result_quats[i*4+3];
        }
        viewdat.close();
    }
    else
    {
        for(unsigned int i=0; i<models.size(); i++)
        {
            result_quats[i*4+3] = 1; // set default quat to [0,0,0,1]
        }
    }

    // set start position
    current_no = 0;
    myrender->load_model(models[current_no]);
    trackball_quat[0] = result_quats[current_no*4 + 0];
    trackball_quat[1] = result_quats[current_no*4 + 1];
    trackball_quat[2] = result_quats[current_no*4 + 2];
    trackball_quat[3] = result_quats[current_no*4 + 3];
    myrender->set_trackball_quat(trackball_quat);
    lbl_no->setText(QString::number(current_no+1) + QString("/") + QString::number(models.size()));
}

void MainWindow::load_model()
{
    filename = QFileDialog::getOpenFileName(
        this,"Open Model File", filename, ("Model Files(*.ply *.obj *.off)"));
    if (!filename.isEmpty())
        myrender->load_model(filename.toStdWString());
    glw->updateGL();
}

void MainWindow::save_result()
{
}

void MainWindow::closeEvent(QCloseEvent *event)
{
    // get current trackball quat
    myrender->get_trackball_quat(trackball_quat);
    result_quats[current_no*4 + 0] = trackball_quat[0];
    result_quats[current_no*4 + 1] = trackball_quat[1];
    result_quats[current_no*4 + 2] = trackball_quat[2];
    result_quats[current_no*4 + 3] = trackball_quat[3];

    // save all result
    QFile viewdat("./view.dat");
    viewdat.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream viewdattx(&viewdat);
    for (unsigned int i=0; i<models.size(); i++)
    {
        viewdattx
            <<result_quats[i*4+0]<<' '
            <<result_quats[i*4+1]<<' '
            <<result_quats[i*4+2]<<' '
            <<result_quats[i*4+3]<<'\n';
    }
    viewdat.close();
    event->accept();
}

void MainWindow::on_btn_next_model_clicked()
{
    btn_next_model->setDisabled(true);
    btn_previous_model->setDisabled(true);

    if (current_no == models.size()-1)
    {
        QMessageBox msgBox;
        msgBox.setIcon(QMessageBox::Question);
        msgBox.setWindowTitle(L2U("完成？"));
        msgBox.setWindowIcon(icon);
        msgBox.setText(L2U("这是最后一个模型。"));
        msgBox.setInformativeText(L2U("如果您对所有的模型都选择好了观察视点，请按\"Yes\"结束。\n或按\"No\"返回修改。\n\n全部修改完成后，退出程序，请您将程序目录下的\"view.dat\"文件反馈给我，非常感谢！"));
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::No);
        int ret = msgBox.exec();
        if (ret == QMessageBox::Yes)
        {
            this->close();
        }
    }
    else
    {
        // save current track ball quat
        myrender->get_trackball_quat(trackball_quat);
        result_quats[current_no*4 + 0] = trackball_quat[0];
        result_quats[current_no*4 + 1] = trackball_quat[1];
        result_quats[current_no*4 + 2] = trackball_quat[2];
        result_quats[current_no*4 + 3] = trackball_quat[3];

        current_no++;
        myrender->load_model(models[current_no]);

        // load previous track ball quat
        trackball_quat[0] = result_quats[current_no*4 + 0];
        trackball_quat[1] = result_quats[current_no*4 + 1];
        trackball_quat[2] = result_quats[current_no*4 + 2];
        trackball_quat[3] = result_quats[current_no*4 + 3];
        myrender->set_trackball_quat(trackball_quat);
        glw->updateGL();
        lbl_no->setText(QString::number(current_no+1) + QString("/") + QString::number(models.size()));
    }

    btn_next_model->setDisabled(false);
    btn_previous_model->setDisabled(false);
}

void MainWindow::on_btn_previous_model_clicked()
{
    btn_next_model->setDisabled(true);
    btn_previous_model->setDisabled(true);

    if (current_no == 0)
    {
        // do nothing
        QMessageBox msgBox;
        msgBox.setIcon(QMessageBox::Information);
        msgBox.setWindowTitle(L2U("第一个模型"));
        msgBox.setWindowIcon(icon);
        msgBox.setText(L2U("这是第一个模型。"));
        msgBox.setStandardButtons(QMessageBox::Yes);
        msgBox.exec();
    }
    else
    {
        // save current track ball quat
        myrender->get_trackball_quat(trackball_quat);
        result_quats[current_no*4 + 0] = trackball_quat[0];
        result_quats[current_no*4 + 1] = trackball_quat[1];
        result_quats[current_no*4 + 2] = trackball_quat[2];
        result_quats[current_no*4 + 3] = trackball_quat[3];

        current_no--;
        myrender->load_model(models[current_no]);

        // load previous track ball quat
        trackball_quat[0] = result_quats[current_no*4 + 0];
        trackball_quat[1] = result_quats[current_no*4 + 1];
        trackball_quat[2] = result_quats[current_no*4 + 2];
        trackball_quat[3] = result_quats[current_no*4 + 3];
        myrender->set_trackball_quat(trackball_quat);
        glw->updateGL();

        // set init
        lbl_no->setText(QString::number(current_no+1) + QString("/") + QString::number(models.size()));
    }

    btn_next_model->setDisabled(false);
    btn_previous_model->setDisabled(false);
}