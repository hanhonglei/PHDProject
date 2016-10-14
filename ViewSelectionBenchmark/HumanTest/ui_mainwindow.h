/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Tue Feb 14 12:55:25 2012
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
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QMainWindow>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *action_load_model;
    QAction *action_exit;
    QAction *action_save_result;
    QAction *actionA;
    QWidget *central_widget;
    QHBoxLayout *hboxLayout;
    QHBoxLayout *certral_layout;
    QWidget *gl_container;
    QGridLayout *gridLayout;
    QGroupBox *groupBox_1;
    QHBoxLayout *horizontalLayout_2;
    QGridLayout *gridLayout_2;
    QPushButton *btn_previous_model;
    QPushButton *btn_next_model;
    QSpacerItem *verticalSpacer_2;
    QLabel *lbl_no;
    QLabel *label;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(486, 461);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/HumanTest.ico"), QSize(), QIcon::Normal, QIcon::Off);
        MainWindow->setWindowIcon(icon);
        action_load_model = new QAction(MainWindow);
        action_load_model->setObjectName(QString::fromUtf8("action_load_model"));
        action_exit = new QAction(MainWindow);
        action_exit->setObjectName(QString::fromUtf8("action_exit"));
        action_save_result = new QAction(MainWindow);
        action_save_result->setObjectName(QString::fromUtf8("action_save_result"));
        actionA = new QAction(MainWindow);
        actionA->setObjectName(QString::fromUtf8("actionA"));
        central_widget = new QWidget(MainWindow);
        central_widget->setObjectName(QString::fromUtf8("central_widget"));
        hboxLayout = new QHBoxLayout(central_widget);
        hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
        certral_layout = new QHBoxLayout();
        certral_layout->setObjectName(QString::fromUtf8("certral_layout"));
        gl_container = new QWidget(central_widget);
        gl_container->setObjectName(QString::fromUtf8("gl_container"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(gl_container->sizePolicy().hasHeightForWidth());
        gl_container->setSizePolicy(sizePolicy);

        certral_layout->addWidget(gl_container);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        groupBox_1 = new QGroupBox(central_widget);
        groupBox_1->setObjectName(QString::fromUtf8("groupBox_1"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(groupBox_1->sizePolicy().hasHeightForWidth());
        groupBox_1->setSizePolicy(sizePolicy1);
        groupBox_1->setMinimumSize(QSize(169, 0));
        horizontalLayout_2 = new QHBoxLayout(groupBox_1);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        btn_previous_model = new QPushButton(groupBox_1);
        btn_previous_model->setObjectName(QString::fromUtf8("btn_previous_model"));
        sizePolicy.setHeightForWidth(btn_previous_model->sizePolicy().hasHeightForWidth());
        btn_previous_model->setSizePolicy(sizePolicy);

        gridLayout_2->addWidget(btn_previous_model, 1, 0, 1, 1);

        btn_next_model = new QPushButton(groupBox_1);
        btn_next_model->setObjectName(QString::fromUtf8("btn_next_model"));
        sizePolicy.setHeightForWidth(btn_next_model->sizePolicy().hasHeightForWidth());
        btn_next_model->setSizePolicy(sizePolicy);

        gridLayout_2->addWidget(btn_next_model, 2, 0, 1, 1);

        verticalSpacer_2 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_2->addItem(verticalSpacer_2, 3, 0, 1, 1);

        lbl_no = new QLabel(groupBox_1);
        lbl_no->setObjectName(QString::fromUtf8("lbl_no"));
        QSizePolicy sizePolicy2(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(lbl_no->sizePolicy().hasHeightForWidth());
        lbl_no->setSizePolicy(sizePolicy2);
        QFont font;
        font.setFamily(QString::fromUtf8("\351\273\221\344\275\223"));
        font.setPointSize(32);
        lbl_no->setFont(font);
        lbl_no->setLayoutDirection(Qt::LeftToRight);
        lbl_no->setAlignment(Qt::AlignCenter);

        gridLayout_2->addWidget(lbl_no, 0, 0, 1, 1);

        label = new QLabel(groupBox_1);
        label->setObjectName(QString::fromUtf8("label"));
        label->setWordWrap(true);

        gridLayout_2->addWidget(label, 4, 0, 2, 1);


        horizontalLayout_2->addLayout(gridLayout_2);


        gridLayout->addWidget(groupBox_1, 0, 0, 1, 1);


        certral_layout->addLayout(gridLayout);


        hboxLayout->addLayout(certral_layout);

        MainWindow->setCentralWidget(central_widget);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Viewselection Benchmark Human Data", 0, QApplication::UnicodeUTF8));
        action_load_model->setText(QApplication::translate("MainWindow", "\350\243\205\350\275\275\346\250\241\345\236\213(&L)", 0, QApplication::UnicodeUTF8));
        action_load_model->setShortcut(QApplication::translate("MainWindow", "Ctrl+O", 0, QApplication::UnicodeUTF8));
        action_exit->setText(QApplication::translate("MainWindow", "\351\200\200\345\207\272(&E)", 0, QApplication::UnicodeUTF8));
        action_exit->setShortcut(QApplication::translate("MainWindow", "Esc", 0, QApplication::UnicodeUTF8));
        action_save_result->setText(QApplication::translate("MainWindow", "\344\277\235\345\255\230\347\273\223\346\236\234", 0, QApplication::UnicodeUTF8));
        action_save_result->setShortcut(QApplication::translate("MainWindow", "Ctrl+S", 0, QApplication::UnicodeUTF8));
        actionA->setText(QApplication::translate("MainWindow", "a", 0, QApplication::UnicodeUTF8));
        groupBox_1->setTitle(QApplication::translate("MainWindow", "Operations", 0, QApplication::UnicodeUTF8));
        btn_previous_model->setText(QApplication::translate("MainWindow", "\344\270\212\344\270\200\344\270\252\346\250\241\345\236\213(&P)", 0, QApplication::UnicodeUTF8));
        btn_next_model->setText(QApplication::translate("MainWindow", "\344\270\213\344\270\200\344\270\252\346\250\241\345\236\213(&N)", 0, QApplication::UnicodeUTF8));
        lbl_no->setText(QApplication::translate("MainWindow", "1", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("MainWindow", "1, \350\257\267\345\257\271\346\257\217\344\270\252\346\250\241\345\236\213\351\200\211\346\213\251\344\270\200\344\270\252\346\202\250\350\247\211\345\276\227\346\234\200\344\271\240\346\203\257\343\200\201\346\234\200\345\245\275\347\232\204\350\247\202\345\257\237\350\247\222\345\272\246\343\200\202\n"
"2, \347\273\223\346\236\234\350\207\252\345\212\250\344\277\235\345\255\230\357\274\214\345\256\214\346\210\220\345\220\216\350\257\267\347\233\264\346\216\245\351\200\200\345\207\272\343\200\202\n"
"3, \351\200\200\345\207\272\344\271\213\345\220\216\357\274\214\346\202\250\350\277\230\345\217\257\344\273\245\345\206\215\346\211\223\345\274\200\347\250\213\345\272\217\344\277\256\346\224\271\345\216\237\346\235\245\347\232\204\347\273\223\346\236\234\343\200\202\n"
"4, \346\234\200\345\220\216\357\274\214\350\257\267\345\260\206view.dat\345\217\215\351\246\210\347\273\231\346\210\221\357\274\214\345\215\201\345\210\206\350\260\242\350\260\242\357\274\201", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
