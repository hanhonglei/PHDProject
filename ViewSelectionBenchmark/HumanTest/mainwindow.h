#pragma  once

#include <QButtonGroup>
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
    //QStreamRedirector *stream_redirector;

    QIcon icon;
    QString filename;

    int current_no;
    vector<wstring> models;
    vector<float> result_quats; // 4*sizeof(models)
    QTime timer;

    float trackball_quat[4];

protected:
    void closeEvent(QCloseEvent *event);

public slots:
    void load_model();
    void save_result();

    void on_btn_next_model_clicked();
    void on_btn_previous_model_clicked();
};