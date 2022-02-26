#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFile>
#include <QFileDialog>
#include <qmessagebox.h>
#include <qmimedata.h>
#include <qcolordialog.h>
#include <vtkRenderer.h>
#include "computenormalsdlg.h"
#include "objectpropertiesdlg.h"
#include "SumatraSettings.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    mSettings = new CSumatraSettings();
    mSettings->ReadSettings();
    ui->sceneWidget->Setup(mSettings);

    createActions();
    createMenus();

    QString message = tr("Default status bar message");
    statusBar()->showMessage(message);

    setWindowTitle(tr("Sumatra (Qt) ") + tr(SUMATRA_VERSION));

    setAcceptDrops(true);
}

MainWindow::~MainWindow()
{
    delete ui;
    if (mObjectPropsDlg)
        delete mObjectPropsDlg;
    if (mSettings)
        delete mSettings;
}

void MainWindow::showOpenFileDialog()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open file"), "",
        "VTK Files (*.vtk) | STL Files (*.stl)");

    // Open file
    QFile file(fileName);
    file.open(QIODevice::ReadOnly);

    // Return on Cancel
    if (!file.exists())
        return;

    openFile(fileName);
}

void MainWindow::openFile(const QString& fileName)
{
    ui->sceneWidget->OpenFile(fileName.toStdString());
    if (mObjectPropsDlg)
        mObjectPropsDlg->FullUpdate();
}

void MainWindow::processComputeNormals()
{
    ComputeNormalsDlg dlg;
    dlg.Set3DScene(ui->sceneWidget->get3DScene());

    if (dlg.exec())
    {
        ui->sceneWidget->get3DScene()->CalculateNormals(dlg.GetSelectedSurface(), dlg.GetReplaceSource(), dlg.GetFlipNormals(), dlg.GetSplitNormals(),
            dlg.GetSplitEdgeAngle());
        forceRendering();
    }
}

void MainWindow::showObjectProperties()
{
    if (ui->sceneWidget->get3DScene()->GetNumberOfSurfaces() < 1)
        return;

    if (!mObjectPropsDlg) 
    {
        C3DScene* t3DScene = ui->sceneWidget->get3DScene();
        mObjectPropsDlg = new ObjectPropertiesDlg(this);
        mObjectPropsDlg->Set3DScene(t3DScene);

        QObject::connect(mObjectPropsDlg, &ObjectPropertiesDlg::valueChanged,
            this, &MainWindow::forceRendering);
    }

    mObjectPropsDlg->PopulateObjectCB();
    mObjectPropsDlg->UpdateAllSceneData();
    mObjectPropsDlg->show();
    mObjectPropsDlg->raise();
    mObjectPropsDlg->activateWindow();
}

void MainWindow::forceRendering()
{
    ui->sceneWidget->ForceRender();
}

void MainWindow::ChooseBackgroundColor()
{
    double col[3];
    ui->sceneWidget->get3DScene()->GetRenderer()->GetBackground(col);

    QColor color = QColorDialog::getColor(QColor(col[0] * 255, col[1] * 255, col[2] * 255), this);
    if (color.isValid())
    {
        ui->sceneWidget->get3DScene()->GetRenderer()->SetBackground(color.red() / 255.0, color.green() / 255.0, color.blue() / 255.0);
        forceRendering();
    }
}

void MainWindow::CutWithPlane()
{
    if (ui->sceneWidget->GetWidgetManipulateObject() == -1)
    {
        //CAnnotateAndManipulateDialog dlg;
        //dlg.Set3Dscene(ui->sceneWidget->Get3DScene());

        //if (dlg.DoModal() == IDOK)
        //{
       //     ui->sceneWidget->SetWidgetManipulateObject(dlg.SelectedObject);
        //    m_wndView.Refresh();
        //}
        ui->sceneWidget->SetWidgetManipulateObject(0);
        CutWithPlaneAct->setChecked(true);
    }
    else
    {
        ui->sceneWidget->SetWidgetManipulateObject(-1);
        //ui->sceneWidget->Refresh();
        CutWithPlaneAct->setChecked(false);
    }
    forceRendering();
}

void MainWindow::createActions()
{
    openAct = new QAction(tr("&Open"), this);
    openAct->setShortcuts(QKeySequence::New);
    openAct->setStatusTip(tr("Open file"));
    connect(openAct, &QAction::triggered, this, &MainWindow::showOpenFileDialog);

    ProcessComputeNormalsAct = new QAction(tr("&Compute normals"), this);
    ProcessComputeNormalsAct->setStatusTip(tr("Compute surface normals"));
    connect(ProcessComputeNormalsAct, &QAction::triggered, this, &MainWindow::processComputeNormals);

    optionsObjectPropAct = new QAction(tr("Object &Properties"), this);
    optionsObjectPropAct->setStatusTip(tr("View and modify object properties"));
    connect(optionsObjectPropAct, &QAction::triggered, this, &MainWindow::showObjectProperties);

    ChooseBackgroundColorAct = new QAction(tr("Choose &Bacground Color"), this);
    ChooseBackgroundColorAct->setStatusTip(tr("Chooose background color"));
    connect(ChooseBackgroundColorAct, &QAction::triggered, this, &MainWindow::ChooseBackgroundColor);

    CutWithPlaneAct = new QAction(tr("Cut With &Plane"), this);
    CutWithPlaneAct->setCheckable(true);
    CutWithPlaneAct->setStatusTip(tr("Cut object with a plane"));
    connect(CutWithPlaneAct, &QAction::triggered, this, &MainWindow::CutWithPlane);
}

void MainWindow::createMenus()
{
    fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(openAct);

    fileMenu = menuBar()->addMenu(tr("&Process"));
    fileMenu->addAction(ProcessComputeNormalsAct);

    fileMenu = menuBar()->addMenu(tr("&Manipulate"));
    fileMenu->addAction(CutWithPlaneAct);

    fileMenu = menuBar()->addMenu(tr("&Options"));
    fileMenu->addAction(optionsObjectPropAct);
    fileMenu->addAction(ChooseBackgroundColorAct);
}

void MainWindow::dragEnterEvent(QDragEnterEvent* event)
{
    event->acceptProposedAction();
}


void MainWindow::dropEvent(QDropEvent* event)
{
    const QMimeData* mimeData = event->mimeData();

    // check for our needed mime type, here a file or a list of files
    if (mimeData->hasUrls())
    {
        QList<QUrl> urlList = mimeData->urls();

        // extract the local paths of the files
        for (int i = 0; i < urlList.size(); i++)
        {
            openFile(urlList.at(i).toLocalFile());
        }
    }    
}
