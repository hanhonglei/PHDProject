/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Tue Feb 14 12:58:23 2012
**      by: Qt User Interface Compiler version 4.8.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QComboBox>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QTextEdit>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *action_load_model;
    QAction *action_exit;
    QAction *action_testing_a_round_original_entropy;
    QAction *action_testing_a_round_modified_entropy;
    QAction *action_testing_a_round_mesh_saliency;
    QWidget *central_widget;
    QHBoxLayout *hboxLayout;
    QHBoxLayout *certral_layout;
    QWidget *gl_container;
    QGridLayout *gridLayout;
    QGroupBox *groupBox_1;
    QGridLayout *gridLayout_3;
    QGridLayout *gridLayout_2;
    QRadioButton *rbtn_original_image;
    QRadioButton *rbtn_depth_buffer_image;
    QRadioButton *rbtn_view_dependent_curvature_image;
    QRadioButton *rbtn_model_space_curvature_image;
    QRadioButton *rbtn_radial_curvature_image;
    QRadioButton *rbtn_mesh_saliency;
    QGroupBox *groupBox_2;
    QGridLayout *gridLayout_5;
    QComboBox *cmb_adaptive_box_size;
    QGroupBox *groupBox_3;
    QGridLayout *gridLayout_7;
    QGridLayout *gridLayout_6;
    QPushButton *btn_current_shannon_entropy;
    QPushButton *btn_current_mesh_saliency;
    QPushButton *btn_sampling_using_mesh_saliency;
    QPushButton *btn_current_shannon_entropy_II;
    QPushButton *btn_sampling_using_revised_entropy_II;
    QPushButton *btn_sampling_using_shannon_entropy;
    QPushButton *btn_sampling_using_shannon_entropy_II;
    QPushButton *btn_sampling_using_revised_entropy;
    QPushButton *btn_current_entropy_modified_II;
    QPushButton *btn_testing_a_round;
    QPushButton *btn_save_result;
    QGroupBox *groupBox_4;
    QGridLayout *gridLayout_8;
    QTextEdit *edit_log;
    QGroupBox *groupBox;
    QVBoxLayout *verticalLayout;
    QComboBox *cmb_number_of_histogram_intervals;
    QMenuBar *menuBar;
    QMenu *menu_File;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(438, 752);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/images/Giraffe.png"), QSize(), QIcon::Normal, QIcon::Off);
        MainWindow->setWindowIcon(icon);
        action_load_model = new QAction(MainWindow);
        action_load_model->setObjectName(QString::fromUtf8("action_load_model"));
        action_exit = new QAction(MainWindow);
        action_exit->setObjectName(QString::fromUtf8("action_exit"));
        action_testing_a_round_original_entropy = new QAction(MainWindow);
        action_testing_a_round_original_entropy->setObjectName(QString::fromUtf8("action_testing_a_round_original_entropy"));
        action_testing_a_round_modified_entropy = new QAction(MainWindow);
        action_testing_a_round_modified_entropy->setObjectName(QString::fromUtf8("action_testing_a_round_modified_entropy"));
        action_testing_a_round_mesh_saliency = new QAction(MainWindow);
        action_testing_a_round_mesh_saliency->setObjectName(QString::fromUtf8("action_testing_a_round_mesh_saliency"));
        central_widget = new QWidget(MainWindow);
        central_widget->setObjectName(QString::fromUtf8("central_widget"));
        hboxLayout = new QHBoxLayout(central_widget);
        hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
        certral_layout = new QHBoxLayout();
        certral_layout->setObjectName(QString::fromUtf8("certral_layout"));
        gl_container = new QWidget(central_widget);
        gl_container->setObjectName(QString::fromUtf8("gl_container"));

        certral_layout->addWidget(gl_container);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        groupBox_1 = new QGroupBox(central_widget);
        groupBox_1->setObjectName(QString::fromUtf8("groupBox_1"));
        gridLayout_3 = new QGridLayout(groupBox_1);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        rbtn_original_image = new QRadioButton(groupBox_1);
        rbtn_original_image->setObjectName(QString::fromUtf8("rbtn_original_image"));
        rbtn_original_image->setChecked(true);

        gridLayout_2->addWidget(rbtn_original_image, 0, 0, 1, 1);

        rbtn_depth_buffer_image = new QRadioButton(groupBox_1);
        rbtn_depth_buffer_image->setObjectName(QString::fromUtf8("rbtn_depth_buffer_image"));

        gridLayout_2->addWidget(rbtn_depth_buffer_image, 0, 1, 1, 1);

        rbtn_view_dependent_curvature_image = new QRadioButton(groupBox_1);
        rbtn_view_dependent_curvature_image->setObjectName(QString::fromUtf8("rbtn_view_dependent_curvature_image"));

        gridLayout_2->addWidget(rbtn_view_dependent_curvature_image, 2, 0, 1, 1);

        rbtn_model_space_curvature_image = new QRadioButton(groupBox_1);
        rbtn_model_space_curvature_image->setObjectName(QString::fromUtf8("rbtn_model_space_curvature_image"));

        gridLayout_2->addWidget(rbtn_model_space_curvature_image, 1, 0, 1, 1);

        rbtn_radial_curvature_image = new QRadioButton(groupBox_1);
        rbtn_radial_curvature_image->setObjectName(QString::fromUtf8("rbtn_radial_curvature_image"));

        gridLayout_2->addWidget(rbtn_radial_curvature_image, 1, 1, 1, 1);

        rbtn_mesh_saliency = new QRadioButton(groupBox_1);
        rbtn_mesh_saliency->setObjectName(QString::fromUtf8("rbtn_mesh_saliency"));

        gridLayout_2->addWidget(rbtn_mesh_saliency, 2, 1, 1, 1);


        gridLayout_3->addLayout(gridLayout_2, 0, 0, 1, 1);


        gridLayout->addWidget(groupBox_1, 0, 0, 1, 1);

        groupBox_2 = new QGroupBox(central_widget);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        gridLayout_5 = new QGridLayout(groupBox_2);
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        cmb_adaptive_box_size = new QComboBox(groupBox_2);
        cmb_adaptive_box_size->setObjectName(QString::fromUtf8("cmb_adaptive_box_size"));

        gridLayout_5->addWidget(cmb_adaptive_box_size, 0, 0, 1, 1);


        gridLayout->addWidget(groupBox_2, 1, 0, 1, 1);

        groupBox_3 = new QGroupBox(central_widget);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        gridLayout_7 = new QGridLayout(groupBox_3);
        gridLayout_7->setObjectName(QString::fromUtf8("gridLayout_7"));
        gridLayout_6 = new QGridLayout();
        gridLayout_6->setObjectName(QString::fromUtf8("gridLayout_6"));
        btn_current_shannon_entropy = new QPushButton(groupBox_3);
        btn_current_shannon_entropy->setObjectName(QString::fromUtf8("btn_current_shannon_entropy"));

        gridLayout_6->addWidget(btn_current_shannon_entropy, 0, 0, 1, 1);

        btn_current_mesh_saliency = new QPushButton(groupBox_3);
        btn_current_mesh_saliency->setObjectName(QString::fromUtf8("btn_current_mesh_saliency"));

        gridLayout_6->addWidget(btn_current_mesh_saliency, 4, 0, 1, 1);

        btn_sampling_using_mesh_saliency = new QPushButton(groupBox_3);
        btn_sampling_using_mesh_saliency->setObjectName(QString::fromUtf8("btn_sampling_using_mesh_saliency"));

        gridLayout_6->addWidget(btn_sampling_using_mesh_saliency, 4, 1, 1, 1);

        btn_current_shannon_entropy_II = new QPushButton(groupBox_3);
        btn_current_shannon_entropy_II->setObjectName(QString::fromUtf8("btn_current_shannon_entropy_II"));

        gridLayout_6->addWidget(btn_current_shannon_entropy_II, 2, 0, 1, 1);

        btn_sampling_using_revised_entropy_II = new QPushButton(groupBox_3);
        btn_sampling_using_revised_entropy_II->setObjectName(QString::fromUtf8("btn_sampling_using_revised_entropy_II"));

        gridLayout_6->addWidget(btn_sampling_using_revised_entropy_II, 3, 1, 1, 1);

        btn_sampling_using_shannon_entropy = new QPushButton(groupBox_3);
        btn_sampling_using_shannon_entropy->setObjectName(QString::fromUtf8("btn_sampling_using_shannon_entropy"));

        gridLayout_6->addWidget(btn_sampling_using_shannon_entropy, 0, 1, 1, 1);

        btn_sampling_using_shannon_entropy_II = new QPushButton(groupBox_3);
        btn_sampling_using_shannon_entropy_II->setObjectName(QString::fromUtf8("btn_sampling_using_shannon_entropy_II"));

        gridLayout_6->addWidget(btn_sampling_using_shannon_entropy_II, 2, 1, 1, 1);

        btn_sampling_using_revised_entropy = new QPushButton(groupBox_3);
        btn_sampling_using_revised_entropy->setObjectName(QString::fromUtf8("btn_sampling_using_revised_entropy"));

        gridLayout_6->addWidget(btn_sampling_using_revised_entropy, 1, 1, 1, 1);

        btn_current_entropy_modified_II = new QPushButton(groupBox_3);
        btn_current_entropy_modified_II->setObjectName(QString::fromUtf8("btn_current_entropy_modified_II"));

        gridLayout_6->addWidget(btn_current_entropy_modified_II, 3, 0, 1, 1);

        btn_testing_a_round = new QPushButton(groupBox_3);
        btn_testing_a_round->setObjectName(QString::fromUtf8("btn_testing_a_round"));

        gridLayout_6->addWidget(btn_testing_a_round, 5, 0, 1, 1);

        btn_save_result = new QPushButton(groupBox_3);
        btn_save_result->setObjectName(QString::fromUtf8("btn_save_result"));

        gridLayout_6->addWidget(btn_save_result, 5, 1, 1, 1);


        gridLayout_7->addLayout(gridLayout_6, 0, 0, 1, 1);


        gridLayout->addWidget(groupBox_3, 3, 0, 1, 1);

        groupBox_4 = new QGroupBox(central_widget);
        groupBox_4->setObjectName(QString::fromUtf8("groupBox_4"));
        gridLayout_8 = new QGridLayout(groupBox_4);
        gridLayout_8->setObjectName(QString::fromUtf8("gridLayout_8"));
        edit_log = new QTextEdit(groupBox_4);
        edit_log->setObjectName(QString::fromUtf8("edit_log"));
        edit_log->viewport()->setProperty("cursor", QVariant(QCursor(Qt::ArrowCursor)));
        edit_log->setAcceptDrops(false);
        edit_log->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
        edit_log->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        edit_log->setReadOnly(true);

        gridLayout_8->addWidget(edit_log, 0, 0, 1, 1);


        gridLayout->addWidget(groupBox_4, 4, 0, 1, 1);

        groupBox = new QGroupBox(central_widget);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        verticalLayout = new QVBoxLayout(groupBox);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        cmb_number_of_histogram_intervals = new QComboBox(groupBox);
        cmb_number_of_histogram_intervals->setObjectName(QString::fromUtf8("cmb_number_of_histogram_intervals"));

        verticalLayout->addWidget(cmb_number_of_histogram_intervals);


        gridLayout->addWidget(groupBox, 2, 0, 1, 1);


        certral_layout->addLayout(gridLayout);


        hboxLayout->addLayout(certral_layout);

        MainWindow->setCentralWidget(central_widget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 438, 19));
        menu_File = new QMenu(menuBar);
        menu_File->setObjectName(QString::fromUtf8("menu_File"));
        MainWindow->setMenuBar(menuBar);
        QWidget::setTabOrder(rbtn_original_image, rbtn_depth_buffer_image);
        QWidget::setTabOrder(rbtn_depth_buffer_image, rbtn_view_dependent_curvature_image);
        QWidget::setTabOrder(rbtn_view_dependent_curvature_image, btn_current_shannon_entropy);
        QWidget::setTabOrder(btn_current_shannon_entropy, edit_log);

        menuBar->addAction(menu_File->menuAction());
        menu_File->addAction(action_load_model);
        menu_File->addSeparator();
        menu_File->addAction(action_exit);

        retranslateUi(MainWindow);

        cmb_adaptive_box_size->setCurrentIndex(0);
        cmb_number_of_histogram_intervals->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "View Selection", 0, QApplication::UnicodeUTF8));
        action_load_model->setText(QApplication::translate("MainWindow", "&Load Model", 0, QApplication::UnicodeUTF8));
        action_load_model->setShortcut(QApplication::translate("MainWindow", "Ctrl+O", 0, QApplication::UnicodeUTF8));
        action_exit->setText(QApplication::translate("MainWindow", "&Exit", 0, QApplication::UnicodeUTF8));
        action_exit->setShortcut(QApplication::translate("MainWindow", "Esc", 0, QApplication::UnicodeUTF8));
        action_testing_a_round_original_entropy->setText(QApplication::translate("MainWindow", "Turn a round and output original entropy", 0, QApplication::UnicodeUTF8));
        action_testing_a_round_modified_entropy->setText(QApplication::translate("MainWindow", "Turn a round and output modified entropy", 0, QApplication::UnicodeUTF8));
        action_testing_a_round_mesh_saliency->setText(QApplication::translate("MainWindow", "Turn a round and output mesh saliency", 0, QApplication::UnicodeUTF8));
        groupBox_1->setTitle(QApplication::translate("MainWindow", "Image Type", 0, QApplication::UnicodeUTF8));
        rbtn_original_image->setText(QApplication::translate("MainWindow", "Original", 0, QApplication::UnicodeUTF8));
        rbtn_depth_buffer_image->setText(QApplication::translate("MainWindow", "Depth buffer", 0, QApplication::UnicodeUTF8));
        rbtn_view_dependent_curvature_image->setText(QApplication::translate("MainWindow", "View-Dependent Curvature", 0, QApplication::UnicodeUTF8));
        rbtn_model_space_curvature_image->setText(QApplication::translate("MainWindow", "Model Space Curvature", 0, QApplication::UnicodeUTF8));
        rbtn_radial_curvature_image->setText(QApplication::translate("MainWindow", "Radial Curvature", 0, QApplication::UnicodeUTF8));
        rbtn_mesh_saliency->setText(QApplication::translate("MainWindow", "Mesh Saliency", 0, QApplication::UnicodeUTF8));
        groupBox_2->setTitle(QApplication::translate("MainWindow", "Adaptive Averaging Curvature Box Size", 0, QApplication::UnicodeUTF8));
        cmb_adaptive_box_size->clear();
        cmb_adaptive_box_size->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "1", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "2", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "4", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "8", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "16", 0, QApplication::UnicodeUTF8)
        );
        groupBox_3->setTitle(QApplication::translate("MainWindow", "Operation", 0, QApplication::UnicodeUTF8));
        btn_current_shannon_entropy->setText(QApplication::translate("MainWindow", "Current Shannon Entropy", 0, QApplication::UnicodeUTF8));
        btn_current_mesh_saliency->setText(QApplication::translate("MainWindow", "Current Mesh Saliency", 0, QApplication::UnicodeUTF8));
        btn_sampling_using_mesh_saliency->setText(QApplication::translate("MainWindow", "Sampling Using Mesh Saliency", 0, QApplication::UnicodeUTF8));
        btn_current_shannon_entropy_II->setText(QApplication::translate("MainWindow", "Current Shannon Entropy II", 0, QApplication::UnicodeUTF8));
        btn_sampling_using_revised_entropy_II->setText(QApplication::translate("MainWindow", "Sampling Using Entropy Modified II", 0, QApplication::UnicodeUTF8));
        btn_sampling_using_shannon_entropy->setText(QApplication::translate("MainWindow", "Sampling Using Shannon Entropy", 0, QApplication::UnicodeUTF8));
        btn_sampling_using_shannon_entropy_II->setText(QApplication::translate("MainWindow", "Sampling Using Shannon Entropy II", 0, QApplication::UnicodeUTF8));
        btn_sampling_using_revised_entropy->setText(QApplication::translate("MainWindow", "Sampling Using Entropy Modified", 0, QApplication::UnicodeUTF8));
        btn_current_entropy_modified_II->setText(QApplication::translate("MainWindow", "Current Entropy Modified II", 0, QApplication::UnicodeUTF8));
        btn_testing_a_round->setText(QApplication::translate("MainWindow", "Testing A Round", 0, QApplication::UnicodeUTF8));
        btn_save_result->setText(QApplication::translate("MainWindow", "Append Current Result to File", 0, QApplication::UnicodeUTF8));
        groupBox_4->setTitle(QApplication::translate("MainWindow", "Log", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("MainWindow", "Number of Histogram Intervals", 0, QApplication::UnicodeUTF8));
        cmb_number_of_histogram_intervals->clear();
        cmb_number_of_histogram_intervals->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "2", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "4", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "8", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "16", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "32", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "64", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "128", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "256", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "512", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "1024", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "2048", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "4096", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "8192", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "16384", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "32768", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "65536", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "131072", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "262144", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "524288", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "1048576", 0, QApplication::UnicodeUTF8)
        );
        menu_File->setTitle(QApplication::translate("MainWindow", "&File", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
