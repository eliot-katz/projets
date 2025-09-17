/********************************************************************************
** Form generated from reading UI file 'MainWindow.ui'
**
** Created by: Qt User Interface Compiler version 5.15.8
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QPlainTextEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralwidget;
    QGridLayout *gridLayout_12;
    QTabWidget *tabWidget;
    QWidget *MarketInstru;
    QGridLayout *gridLayout_13;
    QGridLayout *gridLayout_11;
    QGroupBox *groupBox_2;
    QVBoxLayout *verticalLayout_2;
    QGridLayout *gridLayout_2;
    QDoubleSpinBox *sbK;
    QLabel *label_7;
    QLabel *label_5;
    QLabel *label_6;
    QDoubleSpinBox *sbT;
    QHBoxLayout *horizontalLayout_2;
    QComboBox *cbType;
    QGridLayout *gridLayout_10;
    QPushButton *btnLoadDefaults;
    QPushButton *btnSaveProject;
    QPushButton *btnLoadProject;
    QPushButton *btnReset;
    QGroupBox *groupBox;
    QVBoxLayout *verticalLayout;
    QGridLayout *gridLayout;
    QDoubleSpinBox *sbR;
    QLabel *label_4;
    QLabel *label_3;
    QDoubleSpinBox *sbQ;
    QDoubleSpinBox *sbS0;
    QLabel *label_2;
    QLabel *label;
    QDoubleSpinBox *sbSigma;
    QWidget *Pricing;
    QGridLayout *gridLayout_14;
    QHBoxLayout *horizontalLayout_3;
    QPushButton *btnRunMC;
    QCheckBox *cbLog;
    QGroupBox *groupBox_3;
    QVBoxLayout *verticalLayout_3;
    QGridLayout *gridLayout_3;
    QDoubleSpinBox *sbNTarget;
    QLabel *label_8;
    QDoubleSpinBox *sbBatch;
    QLabel *label_9;
    QDoubleSpinBox *dsbTol;
    QLabel *label_11;
    QLabel *label_10;
    QLabel *label_12;
    QSpinBox *sbNSteps;
    QSpinBox *sbSeed;
    QGroupBox *groupBox_7;
    QVBoxLayout *verticalLayout_8;
    QGridLayout *gridLayout_7;
    QLabel *lblNEff;
    QLabel *lblCIHigh;
    QLabel *lblPrice;
    QLabel *lblCILow;
    QLabel *lblElapsed;
    QLabel *lblSE;
    QLabel *label_27;
    QLabel *label_28;
    QLabel *label_29;
    QLabel *label_30;
    QLabel *label_31;
    QLabel *label_32;
    QWidget *chartContainer;
    QGroupBox *grpLegendPricing;
    QGridLayout *gridLayout_9;
    QVBoxLayout *legendLayoutPricing;
    QWidget *Greeks;
    QGridLayout *gridLayout_15;
    QGroupBox *groupBox_5;
    QVBoxLayout *verticalLayout_5;
    QVBoxLayout *verticalLayout_6;
    QRadioButton *rbBRV;
    QRadioButton *rbPW;
    QRadioButton *rbLRM;
    QCheckBox *cbCRN;
    QCheckBox *cbAnti;
    QGroupBox *groupBox_6;
    QVBoxLayout *verticalLayout_7;
    QGridLayout *gridLayout_4;
    QDoubleSpinBox *dsbBumpAbsSigma;
    QLabel *label_13;
    QLabel *label_14;
    QDoubleSpinBox *dsbBumpAbsR;
    QDoubleSpinBox *dsbBumpRelS0;
    QLabel *label_15;
    QLabel *label_16;
    QDoubleSpinBox *dsbBumpAbsT;
    QPushButton *btnRunGreeks;
    QGroupBox *groupBox_8;
    QVBoxLayout *verticalLayout_9;
    QGridLayout *gridLayout_6;
    QLabel *label_25;
    QLabel *label_26;
    QLabel *label_22;
    QLabel *label_23;
    QLabel *label_24;
    QLabel *lblDelta;
    QLabel *lblVega;
    QLabel *lblGamma;
    QLabel *lblRho;
    QLabel *lblTheta;
    QWidget *greeksChartContainer;
    QWidget *Stress;
    QVBoxLayout *verticalLayout_4;
    QGroupBox *gbStressSummary;
    QHBoxLayout *horizontalLayout_5;
    QGridLayout *gridLayout_16;
    QLabel *label_33;
    QLabel *lblPLError;
    QLabel *lblPL;
    QLabel *lblPriceStressed;
    QLabel *label_21;
    QLabel *label_34;
    QLabel *lblPLTaylor;
    QLabel *label_36;
    QLabel *label_37;
    QLabel *lblElapsedStressed;
    QHBoxLayout *horizontalLayout;
    QWidget *widget;
    QHBoxLayout *horizontalLayout_4;
    QGroupBox *gbStressControls;
    QHBoxLayout *horizontalLayout_6;
    QGridLayout *gridLayout_8;
    QLabel *label_18;
    QDoubleSpinBox *sbSigmaStress;
    QSlider *slRStress;
    QLabel *label_19;
    QLabel *label_17;
    QSlider *slSigmaStress;
    QLabel *label_20;
    QDoubleSpinBox *sbS0Stress;
    QSlider *slS0Stress;
    QDoubleSpinBox *sbRStress;
    QDoubleSpinBox *sbTStress;
    QSlider *slTStress;
    QCheckBox *cbStressAutoRun;
    QPushButton *btnRunStress;
    QGroupBox *gbStressChart;
    QHBoxLayout *horizontalLayout_7;
    QVBoxLayout *verticalLayout_10;
    QWidget *stressChartContainer;
    QWidget *gbStressLegend;
    QHBoxLayout *horizontalLayout_8;
    QHBoxLayout *legendLayoutStress;
    QWidget *Validation;
    QHBoxLayout *horizontalLayout_9;
    QGridLayout *gridLayout_5;
    QGroupBox *gbTestsControls;
    QHBoxLayout *horizontalLayout_10;
    QVBoxLayout *verticalLayout_11;
    QPushButton *btnRunAllTests;
    QSpacerItem *verticalSpacer;
    QVBoxLayout *verticalLayout_13;
    QTableWidget *twTests;
    QSpacerItem *verticalSpacer_2;
    QPlainTextEdit *pteLog;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(915, 593);
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(centralwidget->sizePolicy().hasHeightForWidth());
        centralwidget->setSizePolicy(sizePolicy);
        gridLayout_12 = new QGridLayout(centralwidget);
        gridLayout_12->setObjectName(QString::fromUtf8("gridLayout_12"));
        tabWidget = new QTabWidget(centralwidget);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        MarketInstru = new QWidget();
        MarketInstru->setObjectName(QString::fromUtf8("MarketInstru"));
        sizePolicy.setHeightForWidth(MarketInstru->sizePolicy().hasHeightForWidth());
        MarketInstru->setSizePolicy(sizePolicy);
        gridLayout_13 = new QGridLayout(MarketInstru);
        gridLayout_13->setObjectName(QString::fromUtf8("gridLayout_13"));
        gridLayout_11 = new QGridLayout();
        gridLayout_11->setObjectName(QString::fromUtf8("gridLayout_11"));
        gridLayout_11->setHorizontalSpacing(150);
        gridLayout_11->setVerticalSpacing(100);
        gridLayout_11->setContentsMargins(-1, -1, 50, -1);
        groupBox_2 = new QGroupBox(MarketInstru);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        verticalLayout_2 = new QVBoxLayout(groupBox_2);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        sbK = new QDoubleSpinBox(groupBox_2);
        sbK->setObjectName(QString::fromUtf8("sbK"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(sbK->sizePolicy().hasHeightForWidth());
        sbK->setSizePolicy(sizePolicy1);
        sbK->setDecimals(6);
        sbK->setMinimum(0.010000000000000);
        sbK->setMaximum(1000000000.000000000000000);
        sbK->setValue(100.000000000000000);

        gridLayout_2->addWidget(sbK, 0, 1, 1, 1);

        label_7 = new QLabel(groupBox_2);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_7, 0, 0, 1, 1);

        label_5 = new QLabel(groupBox_2);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        QSizePolicy sizePolicy2(QSizePolicy::Fixed, QSizePolicy::Preferred);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(label_5->sizePolicy().hasHeightForWidth());
        label_5->setSizePolicy(sizePolicy2);
        label_5->setLayoutDirection(Qt::RightToLeft);
        label_5->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_5, 2, 0, 1, 1);

        label_6 = new QLabel(groupBox_2);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_6, 1, 0, 1, 1);

        sbT = new QDoubleSpinBox(groupBox_2);
        sbT->setObjectName(QString::fromUtf8("sbT"));
        sbT->setMaximum(100.000000000000000);
        sbT->setSingleStep(0.010000000000000);
        sbT->setValue(1.000000000000000);

        gridLayout_2->addWidget(sbT, 1, 1, 1, 1);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        cbType = new QComboBox(groupBox_2);
        cbType->addItem(QString());
        cbType->addItem(QString());
        cbType->setObjectName(QString::fromUtf8("cbType"));
        cbType->setEditable(false);

        horizontalLayout_2->addWidget(cbType);


        gridLayout_2->addLayout(horizontalLayout_2, 2, 1, 1, 1);


        verticalLayout_2->addLayout(gridLayout_2);


        gridLayout_11->addWidget(groupBox_2, 2, 1, 1, 1);

        gridLayout_10 = new QGridLayout();
        gridLayout_10->setObjectName(QString::fromUtf8("gridLayout_10"));
        gridLayout_10->setSizeConstraint(QLayout::SetMinimumSize);
        btnLoadDefaults = new QPushButton(MarketInstru);
        btnLoadDefaults->setObjectName(QString::fromUtf8("btnLoadDefaults"));
        sizePolicy1.setHeightForWidth(btnLoadDefaults->sizePolicy().hasHeightForWidth());
        btnLoadDefaults->setSizePolicy(sizePolicy1);

        gridLayout_10->addWidget(btnLoadDefaults, 0, 1, 1, 1);

        btnSaveProject = new QPushButton(MarketInstru);
        btnSaveProject->setObjectName(QString::fromUtf8("btnSaveProject"));
        sizePolicy1.setHeightForWidth(btnSaveProject->sizePolicy().hasHeightForWidth());
        btnSaveProject->setSizePolicy(sizePolicy1);

        gridLayout_10->addWidget(btnSaveProject, 1, 0, 1, 1);

        btnLoadProject = new QPushButton(MarketInstru);
        btnLoadProject->setObjectName(QString::fromUtf8("btnLoadProject"));
        sizePolicy1.setHeightForWidth(btnLoadProject->sizePolicy().hasHeightForWidth());
        btnLoadProject->setSizePolicy(sizePolicy1);

        gridLayout_10->addWidget(btnLoadProject, 1, 1, 1, 1);

        btnReset = new QPushButton(MarketInstru);
        btnReset->setObjectName(QString::fromUtf8("btnReset"));
        sizePolicy1.setHeightForWidth(btnReset->sizePolicy().hasHeightForWidth());
        btnReset->setSizePolicy(sizePolicy1);

        gridLayout_10->addWidget(btnReset, 0, 0, 1, 1);


        gridLayout_11->addLayout(gridLayout_10, 1, 0, 1, 1);

        groupBox = new QGroupBox(MarketInstru);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        verticalLayout = new QVBoxLayout(groupBox);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        sbR = new QDoubleSpinBox(groupBox);
        sbR->setObjectName(QString::fromUtf8("sbR"));
        sbR->setMinimum(-1.000000000000000);
        sbR->setMaximum(1.000000000000000);
        sbR->setSingleStep(0.000100000000000);
        sbR->setValue(0.020000000000000);

        gridLayout->addWidget(sbR, 1, 1, 1, 1);

        label_4 = new QLabel(groupBox);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout->addWidget(label_4, 3, 0, 1, 1);

        label_3 = new QLabel(groupBox);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout->addWidget(label_3, 2, 0, 1, 1);

        sbQ = new QDoubleSpinBox(groupBox);
        sbQ->setObjectName(QString::fromUtf8("sbQ"));
        sbQ->setMinimum(-1.000000000000000);
        sbQ->setMaximum(1.000000000000000);
        sbQ->setSingleStep(0.000100000000000);

        gridLayout->addWidget(sbQ, 2, 1, 1, 1);

        sbS0 = new QDoubleSpinBox(groupBox);
        sbS0->setObjectName(QString::fromUtf8("sbS0"));
        sizePolicy1.setHeightForWidth(sbS0->sizePolicy().hasHeightForWidth());
        sbS0->setSizePolicy(sizePolicy1);
        sbS0->setDecimals(6);
        sbS0->setMinimum(0.010000000000000);
        sbS0->setMaximum(1000000000.000000000000000);
        sbS0->setValue(100.000000000000000);

        gridLayout->addWidget(sbS0, 0, 1, 1, 1);

        label_2 = new QLabel(groupBox);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout->addWidget(label_2, 1, 0, 1, 1);

        label = new QLabel(groupBox);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout->addWidget(label, 0, 0, 1, 1);

        sbSigma = new QDoubleSpinBox(groupBox);
        sbSigma->setObjectName(QString::fromUtf8("sbSigma"));
        sbSigma->setMaximum(5.000000000000000);
        sbSigma->setSingleStep(0.000100000000000);
        sbSigma->setValue(0.020000000000000);

        gridLayout->addWidget(sbSigma, 3, 1, 1, 1);


        verticalLayout->addLayout(gridLayout);


        gridLayout_11->addWidget(groupBox, 1, 1, 1, 1);

        gridLayout_11->setColumnStretch(0, 1);
        gridLayout_11->setColumnStretch(1, 3);

        gridLayout_13->addLayout(gridLayout_11, 0, 0, 1, 1);

        tabWidget->addTab(MarketInstru, QString());
        Pricing = new QWidget();
        Pricing->setObjectName(QString::fromUtf8("Pricing"));
        gridLayout_14 = new QGridLayout(Pricing);
        gridLayout_14->setObjectName(QString::fromUtf8("gridLayout_14"));
        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        btnRunMC = new QPushButton(Pricing);
        btnRunMC->setObjectName(QString::fromUtf8("btnRunMC"));

        horizontalLayout_3->addWidget(btnRunMC);

        cbLog = new QCheckBox(Pricing);
        cbLog->setObjectName(QString::fromUtf8("cbLog"));

        horizontalLayout_3->addWidget(cbLog);


        gridLayout_14->addLayout(horizontalLayout_3, 0, 0, 1, 1);

        groupBox_3 = new QGroupBox(Pricing);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        verticalLayout_3 = new QVBoxLayout(groupBox_3);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        gridLayout_3 = new QGridLayout();
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        sbNTarget = new QDoubleSpinBox(groupBox_3);
        sbNTarget->setObjectName(QString::fromUtf8("sbNTarget"));
        sbNTarget->setDecimals(0);
        sbNTarget->setMinimum(1.000000000000000);
        sbNTarget->setMaximum(100000000.000000000000000);
        sbNTarget->setValue(400000.000000000000000);

        gridLayout_3->addWidget(sbNTarget, 0, 1, 1, 1);

        label_8 = new QLabel(groupBox_3);
        label_8->setObjectName(QString::fromUtf8("label_8"));

        gridLayout_3->addWidget(label_8, 3, 0, 1, 1);

        sbBatch = new QDoubleSpinBox(groupBox_3);
        sbBatch->setObjectName(QString::fromUtf8("sbBatch"));
        sbBatch->setDecimals(0);
        sbBatch->setMinimum(1.000000000000000);
        sbBatch->setMaximum(100000000.000000000000000);
        sbBatch->setSingleStep(1.000000000000000);
        sbBatch->setValue(100000.000000000000000);

        gridLayout_3->addWidget(sbBatch, 1, 1, 1, 1);

        label_9 = new QLabel(groupBox_3);
        label_9->setObjectName(QString::fromUtf8("label_9"));

        gridLayout_3->addWidget(label_9, 2, 0, 1, 1);

        dsbTol = new QDoubleSpinBox(groupBox_3);
        dsbTol->setObjectName(QString::fromUtf8("dsbTol"));
        dsbTol->setDecimals(6);
        dsbTol->setMinimum(-1.000000000000000);
        dsbTol->setMaximum(1.000000000000000);
        dsbTol->setSingleStep(0.000001000000000);
        dsbTol->setValue(-1.000000000000000);

        gridLayout_3->addWidget(dsbTol, 2, 1, 1, 1);

        label_11 = new QLabel(groupBox_3);
        label_11->setObjectName(QString::fromUtf8("label_11"));

        gridLayout_3->addWidget(label_11, 0, 0, 1, 1);

        label_10 = new QLabel(groupBox_3);
        label_10->setObjectName(QString::fromUtf8("label_10"));

        gridLayout_3->addWidget(label_10, 1, 0, 1, 1);

        label_12 = new QLabel(groupBox_3);
        label_12->setObjectName(QString::fromUtf8("label_12"));

        gridLayout_3->addWidget(label_12, 4, 0, 1, 1);

        sbNSteps = new QSpinBox(groupBox_3);
        sbNSteps->setObjectName(QString::fromUtf8("sbNSteps"));
        sbNSteps->setMinimum(1);
        sbNSteps->setMaximum(10000);

        gridLayout_3->addWidget(sbNSteps, 3, 1, 1, 1);

        sbSeed = new QSpinBox(groupBox_3);
        sbSeed->setObjectName(QString::fromUtf8("sbSeed"));
        sbSeed->setMaximum(999999999);
        sbSeed->setValue(42);

        gridLayout_3->addWidget(sbSeed, 4, 1, 1, 1);


        verticalLayout_3->addLayout(gridLayout_3);


        gridLayout_14->addWidget(groupBox_3, 0, 1, 2, 1);

        groupBox_7 = new QGroupBox(Pricing);
        groupBox_7->setObjectName(QString::fromUtf8("groupBox_7"));
        verticalLayout_8 = new QVBoxLayout(groupBox_7);
        verticalLayout_8->setObjectName(QString::fromUtf8("verticalLayout_8"));
        gridLayout_7 = new QGridLayout();
        gridLayout_7->setObjectName(QString::fromUtf8("gridLayout_7"));
        lblNEff = new QLabel(groupBox_7);
        lblNEff->setObjectName(QString::fromUtf8("lblNEff"));

        gridLayout_7->addWidget(lblNEff, 4, 1, 1, 1);

        lblCIHigh = new QLabel(groupBox_7);
        lblCIHigh->setObjectName(QString::fromUtf8("lblCIHigh"));

        gridLayout_7->addWidget(lblCIHigh, 3, 1, 1, 1);

        lblPrice = new QLabel(groupBox_7);
        lblPrice->setObjectName(QString::fromUtf8("lblPrice"));

        gridLayout_7->addWidget(lblPrice, 0, 1, 1, 1);

        lblCILow = new QLabel(groupBox_7);
        lblCILow->setObjectName(QString::fromUtf8("lblCILow"));

        gridLayout_7->addWidget(lblCILow, 2, 1, 1, 1);

        lblElapsed = new QLabel(groupBox_7);
        lblElapsed->setObjectName(QString::fromUtf8("lblElapsed"));

        gridLayout_7->addWidget(lblElapsed, 5, 1, 1, 1);

        lblSE = new QLabel(groupBox_7);
        lblSE->setObjectName(QString::fromUtf8("lblSE"));

        gridLayout_7->addWidget(lblSE, 1, 1, 1, 1);

        label_27 = new QLabel(groupBox_7);
        label_27->setObjectName(QString::fromUtf8("label_27"));

        gridLayout_7->addWidget(label_27, 0, 0, 1, 1);

        label_28 = new QLabel(groupBox_7);
        label_28->setObjectName(QString::fromUtf8("label_28"));

        gridLayout_7->addWidget(label_28, 1, 0, 1, 1);

        label_29 = new QLabel(groupBox_7);
        label_29->setObjectName(QString::fromUtf8("label_29"));

        gridLayout_7->addWidget(label_29, 2, 0, 1, 1);

        label_30 = new QLabel(groupBox_7);
        label_30->setObjectName(QString::fromUtf8("label_30"));

        gridLayout_7->addWidget(label_30, 3, 0, 1, 1);

        label_31 = new QLabel(groupBox_7);
        label_31->setObjectName(QString::fromUtf8("label_31"));

        gridLayout_7->addWidget(label_31, 4, 0, 1, 1);

        label_32 = new QLabel(groupBox_7);
        label_32->setObjectName(QString::fromUtf8("label_32"));

        gridLayout_7->addWidget(label_32, 5, 0, 1, 1);


        verticalLayout_8->addLayout(gridLayout_7);


        gridLayout_14->addWidget(groupBox_7, 1, 0, 2, 1);

        chartContainer = new QWidget(Pricing);
        chartContainer->setObjectName(QString::fromUtf8("chartContainer"));

        gridLayout_14->addWidget(chartContainer, 2, 1, 2, 1);

        grpLegendPricing = new QGroupBox(Pricing);
        grpLegendPricing->setObjectName(QString::fromUtf8("grpLegendPricing"));
        gridLayout_9 = new QGridLayout(grpLegendPricing);
        gridLayout_9->setObjectName(QString::fromUtf8("gridLayout_9"));
        legendLayoutPricing = new QVBoxLayout();
        legendLayoutPricing->setObjectName(QString::fromUtf8("legendLayoutPricing"));

        gridLayout_9->addLayout(legendLayoutPricing, 0, 0, 1, 1);


        gridLayout_14->addWidget(grpLegendPricing, 3, 0, 1, 1);

        gridLayout_14->setRowStretch(2, 2);
        gridLayout_14->setRowStretch(3, 1);
        gridLayout_14->setColumnStretch(1, 2);
        tabWidget->addTab(Pricing, QString());
        Greeks = new QWidget();
        Greeks->setObjectName(QString::fromUtf8("Greeks"));
        gridLayout_15 = new QGridLayout(Greeks);
        gridLayout_15->setObjectName(QString::fromUtf8("gridLayout_15"));
        groupBox_5 = new QGroupBox(Greeks);
        groupBox_5->setObjectName(QString::fromUtf8("groupBox_5"));
        verticalLayout_5 = new QVBoxLayout(groupBox_5);
        verticalLayout_5->setObjectName(QString::fromUtf8("verticalLayout_5"));
        verticalLayout_6 = new QVBoxLayout();
        verticalLayout_6->setObjectName(QString::fromUtf8("verticalLayout_6"));
        rbBRV = new QRadioButton(groupBox_5);
        rbBRV->setObjectName(QString::fromUtf8("rbBRV"));
        rbBRV->setChecked(true);

        verticalLayout_6->addWidget(rbBRV);

        rbPW = new QRadioButton(groupBox_5);
        rbPW->setObjectName(QString::fromUtf8("rbPW"));

        verticalLayout_6->addWidget(rbPW);

        rbLRM = new QRadioButton(groupBox_5);
        rbLRM->setObjectName(QString::fromUtf8("rbLRM"));

        verticalLayout_6->addWidget(rbLRM);

        cbCRN = new QCheckBox(groupBox_5);
        cbCRN->setObjectName(QString::fromUtf8("cbCRN"));
        cbCRN->setChecked(true);

        verticalLayout_6->addWidget(cbCRN);

        cbAnti = new QCheckBox(groupBox_5);
        cbAnti->setObjectName(QString::fromUtf8("cbAnti"));

        verticalLayout_6->addWidget(cbAnti);


        verticalLayout_5->addLayout(verticalLayout_6);


        gridLayout_15->addWidget(groupBox_5, 1, 0, 1, 1);

        groupBox_6 = new QGroupBox(Greeks);
        groupBox_6->setObjectName(QString::fromUtf8("groupBox_6"));
        verticalLayout_7 = new QVBoxLayout(groupBox_6);
        verticalLayout_7->setObjectName(QString::fromUtf8("verticalLayout_7"));
        gridLayout_4 = new QGridLayout();
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        dsbBumpAbsSigma = new QDoubleSpinBox(groupBox_6);
        dsbBumpAbsSigma->setObjectName(QString::fromUtf8("dsbBumpAbsSigma"));
        dsbBumpAbsSigma->setDecimals(4);
        dsbBumpAbsSigma->setMinimum(0.000100000000000);
        dsbBumpAbsSigma->setMaximum(1.000000000000000);
        dsbBumpAbsSigma->setSingleStep(0.000100000000000);
        dsbBumpAbsSigma->setValue(0.010000000000000);

        gridLayout_4->addWidget(dsbBumpAbsSigma, 1, 1, 1, 1);

        label_13 = new QLabel(groupBox_6);
        label_13->setObjectName(QString::fromUtf8("label_13"));

        gridLayout_4->addWidget(label_13, 3, 0, 1, 1);

        label_14 = new QLabel(groupBox_6);
        label_14->setObjectName(QString::fromUtf8("label_14"));

        gridLayout_4->addWidget(label_14, 2, 0, 1, 1);

        dsbBumpAbsR = new QDoubleSpinBox(groupBox_6);
        dsbBumpAbsR->setObjectName(QString::fromUtf8("dsbBumpAbsR"));
        dsbBumpAbsR->setDecimals(4);
        dsbBumpAbsR->setMinimum(0.000000000000000);
        dsbBumpAbsR->setMaximum(0.100000000000000);
        dsbBumpAbsR->setSingleStep(0.000100000000000);
        dsbBumpAbsR->setValue(0.000100000000000);

        gridLayout_4->addWidget(dsbBumpAbsR, 2, 1, 1, 1);

        dsbBumpRelS0 = new QDoubleSpinBox(groupBox_6);
        dsbBumpRelS0->setObjectName(QString::fromUtf8("dsbBumpRelS0"));
        dsbBumpRelS0->setDecimals(4);
        dsbBumpRelS0->setMinimum(0.000100000000000);
        dsbBumpRelS0->setMaximum(0.050000000000000);
        dsbBumpRelS0->setSingleStep(0.010000000000000);
        dsbBumpRelS0->setValue(0.000100000000000);

        gridLayout_4->addWidget(dsbBumpRelS0, 0, 1, 1, 1);

        label_15 = new QLabel(groupBox_6);
        label_15->setObjectName(QString::fromUtf8("label_15"));

        gridLayout_4->addWidget(label_15, 1, 0, 1, 1);

        label_16 = new QLabel(groupBox_6);
        label_16->setObjectName(QString::fromUtf8("label_16"));

        gridLayout_4->addWidget(label_16, 0, 0, 1, 1);

        dsbBumpAbsT = new QDoubleSpinBox(groupBox_6);
        dsbBumpAbsT->setObjectName(QString::fromUtf8("dsbBumpAbsT"));
        dsbBumpAbsT->setDecimals(6);
        dsbBumpAbsT->setMaximum(1.000000000000000);
        dsbBumpAbsT->setSingleStep(0.000100000000000);
        dsbBumpAbsT->setValue(0.002739000000000);

        gridLayout_4->addWidget(dsbBumpAbsT, 3, 1, 1, 1);


        verticalLayout_7->addLayout(gridLayout_4);


        gridLayout_15->addWidget(groupBox_6, 2, 0, 1, 1);

        btnRunGreeks = new QPushButton(Greeks);
        btnRunGreeks->setObjectName(QString::fromUtf8("btnRunGreeks"));

        gridLayout_15->addWidget(btnRunGreeks, 0, 0, 1, 1);

        groupBox_8 = new QGroupBox(Greeks);
        groupBox_8->setObjectName(QString::fromUtf8("groupBox_8"));
        verticalLayout_9 = new QVBoxLayout(groupBox_8);
        verticalLayout_9->setObjectName(QString::fromUtf8("verticalLayout_9"));
        gridLayout_6 = new QGridLayout();
        gridLayout_6->setObjectName(QString::fromUtf8("gridLayout_6"));
        label_25 = new QLabel(groupBox_8);
        label_25->setObjectName(QString::fromUtf8("label_25"));

        gridLayout_6->addWidget(label_25, 3, 0, 1, 1);

        label_26 = new QLabel(groupBox_8);
        label_26->setObjectName(QString::fromUtf8("label_26"));

        gridLayout_6->addWidget(label_26, 4, 0, 1, 1);

        label_22 = new QLabel(groupBox_8);
        label_22->setObjectName(QString::fromUtf8("label_22"));

        gridLayout_6->addWidget(label_22, 0, 0, 1, 1);

        label_23 = new QLabel(groupBox_8);
        label_23->setObjectName(QString::fromUtf8("label_23"));

        gridLayout_6->addWidget(label_23, 1, 0, 1, 1);

        label_24 = new QLabel(groupBox_8);
        label_24->setObjectName(QString::fromUtf8("label_24"));

        gridLayout_6->addWidget(label_24, 2, 0, 1, 1);

        lblDelta = new QLabel(groupBox_8);
        lblDelta->setObjectName(QString::fromUtf8("lblDelta"));

        gridLayout_6->addWidget(lblDelta, 0, 1, 1, 1);

        lblVega = new QLabel(groupBox_8);
        lblVega->setObjectName(QString::fromUtf8("lblVega"));

        gridLayout_6->addWidget(lblVega, 1, 1, 1, 1);

        lblGamma = new QLabel(groupBox_8);
        lblGamma->setObjectName(QString::fromUtf8("lblGamma"));

        gridLayout_6->addWidget(lblGamma, 2, 1, 1, 1);

        lblRho = new QLabel(groupBox_8);
        lblRho->setObjectName(QString::fromUtf8("lblRho"));

        gridLayout_6->addWidget(lblRho, 3, 1, 1, 1);

        lblTheta = new QLabel(groupBox_8);
        lblTheta->setObjectName(QString::fromUtf8("lblTheta"));

        gridLayout_6->addWidget(lblTheta, 4, 1, 1, 1);


        verticalLayout_9->addLayout(gridLayout_6);

        greeksChartContainer = new QWidget(groupBox_8);
        greeksChartContainer->setObjectName(QString::fromUtf8("greeksChartContainer"));

        verticalLayout_9->addWidget(greeksChartContainer);


        gridLayout_15->addWidget(groupBox_8, 0, 1, 3, 1);

        gridLayout_15->setColumnStretch(1, 1);
        tabWidget->addTab(Greeks, QString());
        Stress = new QWidget();
        Stress->setObjectName(QString::fromUtf8("Stress"));
        verticalLayout_4 = new QVBoxLayout(Stress);
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        gbStressSummary = new QGroupBox(Stress);
        gbStressSummary->setObjectName(QString::fromUtf8("gbStressSummary"));
        horizontalLayout_5 = new QHBoxLayout(gbStressSummary);
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        gridLayout_16 = new QGridLayout();
        gridLayout_16->setObjectName(QString::fromUtf8("gridLayout_16"));
        label_33 = new QLabel(gbStressSummary);
        label_33->setObjectName(QString::fromUtf8("label_33"));

        gridLayout_16->addWidget(label_33, 1, 0, 1, 1);

        lblPLError = new QLabel(gbStressSummary);
        lblPLError->setObjectName(QString::fromUtf8("lblPLError"));

        gridLayout_16->addWidget(lblPLError, 2, 1, 1, 1);

        lblPL = new QLabel(gbStressSummary);
        lblPL->setObjectName(QString::fromUtf8("lblPL"));

        gridLayout_16->addWidget(lblPL, 0, 1, 1, 1);

        lblPriceStressed = new QLabel(gbStressSummary);
        lblPriceStressed->setObjectName(QString::fromUtf8("lblPriceStressed"));

        gridLayout_16->addWidget(lblPriceStressed, 3, 1, 1, 1);

        label_21 = new QLabel(gbStressSummary);
        label_21->setObjectName(QString::fromUtf8("label_21"));

        gridLayout_16->addWidget(label_21, 0, 0, 1, 1);

        label_34 = new QLabel(gbStressSummary);
        label_34->setObjectName(QString::fromUtf8("label_34"));

        gridLayout_16->addWidget(label_34, 2, 0, 1, 1);

        lblPLTaylor = new QLabel(gbStressSummary);
        lblPLTaylor->setObjectName(QString::fromUtf8("lblPLTaylor"));

        gridLayout_16->addWidget(lblPLTaylor, 1, 1, 1, 1);

        label_36 = new QLabel(gbStressSummary);
        label_36->setObjectName(QString::fromUtf8("label_36"));

        gridLayout_16->addWidget(label_36, 3, 0, 1, 1);

        label_37 = new QLabel(gbStressSummary);
        label_37->setObjectName(QString::fromUtf8("label_37"));

        gridLayout_16->addWidget(label_37, 4, 0, 1, 1);

        lblElapsedStressed = new QLabel(gbStressSummary);
        lblElapsedStressed->setObjectName(QString::fromUtf8("lblElapsedStressed"));

        gridLayout_16->addWidget(lblElapsedStressed, 4, 1, 1, 1);


        horizontalLayout_5->addLayout(gridLayout_16);


        verticalLayout_4->addWidget(gbStressSummary);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        widget = new QWidget(Stress);
        widget->setObjectName(QString::fromUtf8("widget"));
        horizontalLayout_4 = new QHBoxLayout(widget);
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        gbStressControls = new QGroupBox(widget);
        gbStressControls->setObjectName(QString::fromUtf8("gbStressControls"));
        QSizePolicy sizePolicy3(QSizePolicy::Preferred, QSizePolicy::Expanding);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(gbStressControls->sizePolicy().hasHeightForWidth());
        gbStressControls->setSizePolicy(sizePolicy3);
        gbStressControls->setMinimumSize(QSize(260, 0));
        horizontalLayout_6 = new QHBoxLayout(gbStressControls);
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        gridLayout_8 = new QGridLayout();
        gridLayout_8->setObjectName(QString::fromUtf8("gridLayout_8"));
        label_18 = new QLabel(gbStressControls);
        label_18->setObjectName(QString::fromUtf8("label_18"));

        gridLayout_8->addWidget(label_18, 1, 0, 1, 1);

        sbSigmaStress = new QDoubleSpinBox(gbStressControls);
        sbSigmaStress->setObjectName(QString::fromUtf8("sbSigmaStress"));
        sbSigmaStress->setDecimals(1);
        sbSigmaStress->setMinimum(-10.000000000000000);
        sbSigmaStress->setMaximum(10.000000000000000);
        sbSigmaStress->setSingleStep(0.500000000000000);

        gridLayout_8->addWidget(sbSigmaStress, 1, 2, 1, 1);

        slRStress = new QSlider(gbStressControls);
        slRStress->setObjectName(QString::fromUtf8("slRStress"));
        slRStress->setMinimum(-100);
        slRStress->setMaximum(100);
        slRStress->setOrientation(Qt::Horizontal);

        gridLayout_8->addWidget(slRStress, 2, 1, 1, 1);

        label_19 = new QLabel(gbStressControls);
        label_19->setObjectName(QString::fromUtf8("label_19"));

        gridLayout_8->addWidget(label_19, 2, 0, 1, 1);

        label_17 = new QLabel(gbStressControls);
        label_17->setObjectName(QString::fromUtf8("label_17"));

        gridLayout_8->addWidget(label_17, 0, 0, 1, 1);

        slSigmaStress = new QSlider(gbStressControls);
        slSigmaStress->setObjectName(QString::fromUtf8("slSigmaStress"));
        slSigmaStress->setMinimum(-1000);
        slSigmaStress->setMaximum(1000);
        slSigmaStress->setOrientation(Qt::Horizontal);

        gridLayout_8->addWidget(slSigmaStress, 1, 1, 1, 1);

        label_20 = new QLabel(gbStressControls);
        label_20->setObjectName(QString::fromUtf8("label_20"));

        gridLayout_8->addWidget(label_20, 3, 0, 1, 1);

        sbS0Stress = new QDoubleSpinBox(gbStressControls);
        sbS0Stress->setObjectName(QString::fromUtf8("sbS0Stress"));
        sbS0Stress->setDecimals(1);
        sbS0Stress->setMinimum(-50.000000000000000);
        sbS0Stress->setMaximum(50.000000000000000);
        sbS0Stress->setSingleStep(0.500000000000000);
        sbS0Stress->setValue(0.000000000000000);

        gridLayout_8->addWidget(sbS0Stress, 0, 2, 1, 1);

        slS0Stress = new QSlider(gbStressControls);
        slS0Stress->setObjectName(QString::fromUtf8("slS0Stress"));
        slS0Stress->setMinimum(-2000);
        slS0Stress->setMaximum(2000);
        slS0Stress->setOrientation(Qt::Horizontal);

        gridLayout_8->addWidget(slS0Stress, 0, 1, 1, 1);

        sbRStress = new QDoubleSpinBox(gbStressControls);
        sbRStress->setObjectName(QString::fromUtf8("sbRStress"));
        sbRStress->setMinimum(-200.000000000000000);
        sbRStress->setMaximum(200.000000000000000);
        sbRStress->setSingleStep(5.000000000000000);

        gridLayout_8->addWidget(sbRStress, 2, 2, 1, 1);

        sbTStress = new QDoubleSpinBox(gbStressControls);
        sbTStress->setObjectName(QString::fromUtf8("sbTStress"));
        sbTStress->setMinimum(-90.000000000000000);
        sbTStress->setMaximum(90.000000000000000);

        gridLayout_8->addWidget(sbTStress, 3, 2, 1, 1);

        slTStress = new QSlider(gbStressControls);
        slTStress->setObjectName(QString::fromUtf8("slTStress"));
        slTStress->setMinimum(-90);
        slTStress->setMaximum(90);
        slTStress->setOrientation(Qt::Horizontal);

        gridLayout_8->addWidget(slTStress, 3, 1, 1, 1);

        cbStressAutoRun = new QCheckBox(gbStressControls);
        cbStressAutoRun->setObjectName(QString::fromUtf8("cbStressAutoRun"));
        cbStressAutoRun->setChecked(true);

        gridLayout_8->addWidget(cbStressAutoRun, 4, 0, 1, 1);

        btnRunStress = new QPushButton(gbStressControls);
        btnRunStress->setObjectName(QString::fromUtf8("btnRunStress"));

        gridLayout_8->addWidget(btnRunStress, 4, 2, 1, 1);


        horizontalLayout_6->addLayout(gridLayout_8);


        horizontalLayout_4->addWidget(gbStressControls);


        horizontalLayout->addWidget(widget);

        gbStressChart = new QGroupBox(Stress);
        gbStressChart->setObjectName(QString::fromUtf8("gbStressChart"));
        sizePolicy.setHeightForWidth(gbStressChart->sizePolicy().hasHeightForWidth());
        gbStressChart->setSizePolicy(sizePolicy);
        horizontalLayout_7 = new QHBoxLayout(gbStressChart);
        horizontalLayout_7->setObjectName(QString::fromUtf8("horizontalLayout_7"));
        verticalLayout_10 = new QVBoxLayout();
        verticalLayout_10->setObjectName(QString::fromUtf8("verticalLayout_10"));
        stressChartContainer = new QWidget(gbStressChart);
        stressChartContainer->setObjectName(QString::fromUtf8("stressChartContainer"));

        verticalLayout_10->addWidget(stressChartContainer);

        gbStressLegend = new QWidget(gbStressChart);
        gbStressLegend->setObjectName(QString::fromUtf8("gbStressLegend"));
        horizontalLayout_8 = new QHBoxLayout(gbStressLegend);
        horizontalLayout_8->setObjectName(QString::fromUtf8("horizontalLayout_8"));
        legendLayoutStress = new QHBoxLayout();
        legendLayoutStress->setObjectName(QString::fromUtf8("legendLayoutStress"));

        horizontalLayout_8->addLayout(legendLayoutStress);


        verticalLayout_10->addWidget(gbStressLegend);

        verticalLayout_10->setStretch(0, 3);
        verticalLayout_10->setStretch(1, 1);

        horizontalLayout_7->addLayout(verticalLayout_10);


        horizontalLayout->addWidget(gbStressChart);

        horizontalLayout->setStretch(0, 1);
        horizontalLayout->setStretch(1, 3);

        verticalLayout_4->addLayout(horizontalLayout);

        verticalLayout_4->setStretch(0, 1);
        verticalLayout_4->setStretch(1, 4);
        tabWidget->addTab(Stress, QString());
        Validation = new QWidget();
        Validation->setObjectName(QString::fromUtf8("Validation"));
        horizontalLayout_9 = new QHBoxLayout(Validation);
        horizontalLayout_9->setObjectName(QString::fromUtf8("horizontalLayout_9"));
        gridLayout_5 = new QGridLayout();
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        gbTestsControls = new QGroupBox(Validation);
        gbTestsControls->setObjectName(QString::fromUtf8("gbTestsControls"));
        horizontalLayout_10 = new QHBoxLayout(gbTestsControls);
        horizontalLayout_10->setObjectName(QString::fromUtf8("horizontalLayout_10"));
        verticalLayout_11 = new QVBoxLayout();
        verticalLayout_11->setObjectName(QString::fromUtf8("verticalLayout_11"));
        btnRunAllTests = new QPushButton(gbTestsControls);
        btnRunAllTests->setObjectName(QString::fromUtf8("btnRunAllTests"));

        verticalLayout_11->addWidget(btnRunAllTests);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_11->addItem(verticalSpacer);


        horizontalLayout_10->addLayout(verticalLayout_11);


        gridLayout_5->addWidget(gbTestsControls, 0, 0, 1, 1);


        horizontalLayout_9->addLayout(gridLayout_5);

        verticalLayout_13 = new QVBoxLayout();
        verticalLayout_13->setObjectName(QString::fromUtf8("verticalLayout_13"));
        twTests = new QTableWidget(Validation);
        twTests->setObjectName(QString::fromUtf8("twTests"));
        twTests->setEditTriggers(QAbstractItemView::NoEditTriggers);
        twTests->setProperty("showDropIndicator", QVariant(false));
        twTests->setAlternatingRowColors(true);
        twTests->setSelectionMode(QAbstractItemView::SingleSelection);
        twTests->setSelectionBehavior(QAbstractItemView::SelectRows);
        twTests->setShowGrid(false);
        twTests->setWordWrap(false);
        twTests->horizontalHeader()->setStretchLastSection(true);

        verticalLayout_13->addWidget(twTests);

        verticalSpacer_2 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_13->addItem(verticalSpacer_2);

        pteLog = new QPlainTextEdit(Validation);
        pteLog->setObjectName(QString::fromUtf8("pteLog"));
        pteLog->setLineWrapMode(QPlainTextEdit::NoWrap);
        pteLog->setReadOnly(true);

        verticalLayout_13->addWidget(pteLog);

        verticalLayout_13->setStretch(0, 3);
        verticalLayout_13->setStretch(2, 1);

        horizontalLayout_9->addLayout(verticalLayout_13);

        tabWidget->addTab(Validation, QString());

        gridLayout_12->addWidget(tabWidget, 0, 0, 1, 1);

        MainWindow->setCentralWidget(centralwidget);

        retranslateUi(MainWindow);

        tabWidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QCoreApplication::translate("MainWindow", "MainWindow", nullptr));
        groupBox_2->setTitle(QCoreApplication::translate("MainWindow", "Instrument", nullptr));
        label_7->setText(QCoreApplication::translate("MainWindow", "K", nullptr));
        label_5->setText(QCoreApplication::translate("MainWindow", "Type", nullptr));
        label_6->setText(QCoreApplication::translate("MainWindow", "T", nullptr));
        cbType->setItemText(0, QCoreApplication::translate("MainWindow", "Call", nullptr));
        cbType->setItemText(1, QCoreApplication::translate("MainWindow", "Put", nullptr));

        btnLoadDefaults->setText(QCoreApplication::translate("MainWindow", "Defaults", nullptr));
        btnSaveProject->setText(QCoreApplication::translate("MainWindow", "Save project", nullptr));
        btnLoadProject->setText(QCoreApplication::translate("MainWindow", "Load project", nullptr));
        btnReset->setText(QCoreApplication::translate("MainWindow", "Reset", nullptr));
        groupBox->setTitle(QCoreApplication::translate("MainWindow", "Market", nullptr));
        label_4->setText(QCoreApplication::translate("MainWindow", "Sigma", nullptr));
        label_3->setText(QCoreApplication::translate("MainWindow", "Q", nullptr));
        label_2->setText(QCoreApplication::translate("MainWindow", "R", nullptr));
        label->setText(QCoreApplication::translate("MainWindow", "S0", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(MarketInstru), QCoreApplication::translate("MainWindow", "Market & Instrument", nullptr));
        btnRunMC->setText(QCoreApplication::translate("MainWindow", "Run MC", nullptr));
        cbLog->setText(QCoreApplication::translate("MainWindow", "Convergence log", nullptr));
        groupBox_3->setTitle(QCoreApplication::translate("MainWindow", "MC Config", nullptr));
        label_8->setText(QCoreApplication::translate("MainWindow", "Number of Steps", nullptr));
        label_9->setText(QCoreApplication::translate("MainWindow", "Tol", nullptr));
        label_11->setText(QCoreApplication::translate("MainWindow", "Target", nullptr));
        label_10->setText(QCoreApplication::translate("MainWindow", "Batch", nullptr));
        label_12->setText(QCoreApplication::translate("MainWindow", "Seed", nullptr));
        groupBox_7->setTitle(QCoreApplication::translate("MainWindow", "Results", nullptr));
        lblNEff->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        lblCIHigh->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        lblPrice->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        lblCILow->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        lblElapsed->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        lblSE->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        label_27->setText(QCoreApplication::translate("MainWindow", "Price", nullptr));
        label_28->setText(QCoreApplication::translate("MainWindow", "SE", nullptr));
        label_29->setText(QCoreApplication::translate("MainWindow", "CI Low", nullptr));
        label_30->setText(QCoreApplication::translate("MainWindow", "CI High", nullptr));
        label_31->setText(QCoreApplication::translate("MainWindow", "Number eff", nullptr));
        label_32->setText(QCoreApplication::translate("MainWindow", "Elapsed", nullptr));
        grpLegendPricing->setTitle(QCoreApplication::translate("MainWindow", "Legend", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(Pricing), QCoreApplication::translate("MainWindow", "Pricing", nullptr));
        groupBox_5->setTitle(QCoreApplication::translate("MainWindow", "Method", nullptr));
        rbBRV->setText(QCoreApplication::translate("MainWindow", "BRV", nullptr));
        rbPW->setText(QCoreApplication::translate("MainWindow", "PW", nullptr));
        rbLRM->setText(QCoreApplication::translate("MainWindow", "LRM", nullptr));
        cbCRN->setText(QCoreApplication::translate("MainWindow", "CRN", nullptr));
        cbAnti->setText(QCoreApplication::translate("MainWindow", "Anti", nullptr));
        groupBox_6->setTitle(QCoreApplication::translate("MainWindow", "Bumps", nullptr));
        label_13->setText(QCoreApplication::translate("MainWindow", "Maturity  (\316\224T)", nullptr));
        label_14->setText(QCoreApplication::translate("MainWindow", "Rate  (\316\224r)", nullptr));
        label_15->setText(QCoreApplication::translate("MainWindow", "Vol (\316\224\317\203)", nullptr));
        label_16->setText(QCoreApplication::translate("MainWindow", "Spot (\316\224S/S)", nullptr));
        btnRunGreeks->setText(QCoreApplication::translate("MainWindow", "Compute Greeks", nullptr));
        groupBox_8->setTitle(QCoreApplication::translate("MainWindow", "Output", nullptr));
        label_25->setText(QCoreApplication::translate("MainWindow", "Rho", nullptr));
        label_26->setText(QCoreApplication::translate("MainWindow", "Theta", nullptr));
        label_22->setText(QCoreApplication::translate("MainWindow", "Delta", nullptr));
        label_23->setText(QCoreApplication::translate("MainWindow", "Vega", nullptr));
        label_24->setText(QCoreApplication::translate("MainWindow", "Gamma", nullptr));
        lblDelta->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        lblVega->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        lblGamma->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        lblRho->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        lblTheta->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(Greeks), QCoreApplication::translate("MainWindow", "Greeks", nullptr));
        gbStressSummary->setTitle(QCoreApplication::translate("MainWindow", "Summary", nullptr));
        label_33->setText(QCoreApplication::translate("MainWindow", "P&L (Taylor)", nullptr));
        lblPLError->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        lblPL->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        lblPriceStressed->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        label_21->setText(QCoreApplication::translate("MainWindow", "P&L", nullptr));
        label_34->setText(QCoreApplication::translate("MainWindow", "MC \342\210\222 Taylor", nullptr));
        lblPLTaylor->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        label_36->setText(QCoreApplication::translate("MainWindow", "Stressed Price", nullptr));
        label_37->setText(QCoreApplication::translate("MainWindow", "Elapsed (ms)", nullptr));
        lblElapsedStressed->setText(QCoreApplication::translate("MainWindow", "0", nullptr));
        gbStressControls->setTitle(QCoreApplication::translate("MainWindow", "Controls", nullptr));
        label_18->setText(QCoreApplication::translate("MainWindow", "Vol (\316\224\317\203, pts)", nullptr));
        label_19->setText(QCoreApplication::translate("MainWindow", "Rate (\316\224r, bps)", nullptr));
        label_17->setText(QCoreApplication::translate("MainWindow", "Spot (\316\224S, %)", nullptr));
        label_20->setText(QCoreApplication::translate("MainWindow", "Maturity (\316\224T, days)", nullptr));
        cbStressAutoRun->setText(QCoreApplication::translate("MainWindow", "AutoRun", nullptr));
        btnRunStress->setText(QCoreApplication::translate("MainWindow", "Reprice", nullptr));
        gbStressChart->setTitle(QCoreApplication::translate("MainWindow", "Chart", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(Stress), QCoreApplication::translate("MainWindow", "Stress", nullptr));
        gbTestsControls->setTitle(QCoreApplication::translate("MainWindow", "Controls", nullptr));
        btnRunAllTests->setText(QCoreApplication::translate("MainWindow", "Run all tests", nullptr));
        pteLog->setPlaceholderText(QString());
        tabWidget->setTabText(tabWidget->indexOf(Validation), QCoreApplication::translate("MainWindow", "Validation", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
