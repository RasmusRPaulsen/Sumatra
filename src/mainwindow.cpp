#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFile>
#include <QFileDialog>
#include <qmessagebox.h>
#include <qmimedata.h>
#include <qcolordialog.h>
#include <qinputdialog.h>
#include <qstringlist.h>

#include <vtkRenderer.h>

#include <sstream>

#include "computenormalsdlg.h"
#include "objectpropertiesdlg.h"
#include "manipulateobjectdlg.h"
#include "savefiledlg.h"
#include "featureedgesdlg.h"
#include "connectivitydlg.h"
#include "SumatraSettings.h"
#include "createprimitivedlg.h"

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

    statusBar()->showMessage(statusMessage);

    setWindowTitle(tr("Sumatra (Qt) ") + tr(SUMATRA_VERSION) + tr(" (build: ") 
        + tr(__DATE__) + " " + tr(__TIME__) + tr(")"));

    setAcceptDrops(true);

    QObject::connect(ui->sceneWidget, &SceneWidget::updateStatusMessage,
        this, &MainWindow::updateStatusBarMessage);

}

MainWindow::~MainWindow()
{
    delete ui;
    if (mObjectPropsDlg)
        delete mObjectPropsDlg;
    if (mSettings)
        delete mSettings;
}


void MainWindow::showSaveFileDialog()
{
    SaveFileDlg dlg;
    dlg.Set3DScene(ui->sceneWidget->get3DScene());
    if (dlg.exec())
    {
        std::string defname = ui->sceneWidget->get3DScene()->GetSurfaceShortName(dlg.GetSelectedSurface());

        std::string extFilt = "VTK XML File (*.vtp);;VTK File (*.vtk);;STL File (*.stl);;PLY File (*.ply);;All Files (*.*)";
        //std::string extFilt = "Surface File (*.vtk *.vtp *.stl *.txt *.ply);;All Files (*.*)";

        QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
            tr(defname.c_str()), tr(extFilt.c_str()));

        if (fileName != "")
        {
            bool result = ui->sceneWidget->get3DScene()->SaveFile(dlg.GetSelectedSurface(), fileName.toStdString(),
                dlg.getApplyActorTransform(), !dlg.getWriteAsAscii(), dlg.getWriteNormals(), dlg.getWriteScalars());
            if (!result)
                QMessageBox::warning(this, "Warning", "Could not write file");
        }
    }
}

void MainWindow::showAxes()
{
    ui->sceneWidget->get3DScene()->SetAxesVisible(!ui->sceneWidget->get3DScene()->GetAxesVisible());
    showAxesAct->setChecked(ui->sceneWidget->get3DScene()->GetAxesVisible());
    forceRendering();
}

void MainWindow::showScalarBar()
{
    ui->sceneWidget->get3DScene()->SetScalarBarVisible(!ui->sceneWidget->get3DScene()->GetScalarBarVisible());
    forceRendering();
}

void MainWindow::setMarkerSphereValue()
{
    bool ok = false;
    double defVal = ui->sceneWidget->get3DScene()->GetMarkerValue();

    double v = QInputDialog::getDouble(this, tr("Set marker sphere value"),
        tr("Marker sphere value"), defVal, -2147483647, 2147483647, 2, &ok,
        Qt::WindowFlags(), 1);
    if (ok)
    {
        ui->sceneWidget->get3DScene()->SetMarkerValue(v);
        std::ostringstream ost;
        ost << "Marker value: " << ui->sceneWidget->get3DScene()->GetMarkerValue() << 
            " radius: " << ui->sceneWidget->get3DScene()->GetSphereWidgetRadius();
        ui->sceneWidget->get3DScene()->SetStatusText(ost.str(), true);
        updateStatusBarMessage(QString::fromStdString(ost.str()));
        forceRendering();
    }
}

int MainWindow::selectObjectDialog()
{
    if (ui->sceneWidget->get3DScene()->GetNumberOfSurfaces() <= 0)
        return -1;

    QStringList items;
    for (int i = 0; i < ui->sceneWidget->get3DScene()->GetNumberOfSurfaces(); i++)
    {
        std::ostringstream ost;
        ost << std::setfill('0') << std::setw(3) << i + 1 << " : " << ui->sceneWidget->get3DScene()->GetSurfaceShortName(i);
        items << QString::fromStdString(ost.str());
    }
    bool ok;
    QString item = QInputDialog::getItem(this, tr("Select object"),
        tr("Object:"), items, 0, false, &ok);
    if (ok && !item.isEmpty())
    {
        int idx = items.indexOf(item);
        if (idx >= 0)
        {
            // QMessageBox::information(this, "Index returned", QString::number(idx));
            return idx;
        }
    }
    return -1;
}


void MainWindow::visualizeNormals()
{
    int idx = selectObjectDialog();
    if (idx == -1)
        return;

    ui->sceneWidget->get3DScene()->VisualiseNormals(idx, false);
    forceRendering();
}

void MainWindow::visualizeFeatureEdges()
{
    FeatureEdgesDlg dlg;
    dlg.Set3DScene(ui->sceneWidget->get3DScene());

    if (dlg.exec())
    {
        bool boundary = false;
        bool nonManifold = false;
        bool manifold = false;
        bool sharp = false;
        double sharpAngle = 30;
        dlg.getValues(boundary, nonManifold, manifold, sharp, sharpAngle);
        if (!boundary && !nonManifold && !manifold && !sharp)
        {
            QMessageBox::warning(this, "Warning", "At least one edge type should be selected");
        }
        else
        {
            ui->sceneWidget->get3DScene()->ComputeFeatureEdges(dlg.GetSelectedSurface(), boundary, nonManifold,
                manifold, sharp, sharpAngle);
            //ui->sceneWidget->get3DScene()->CalculateNormals(dlg.GetSelectedSurface(), dlg.GetReplaceSource(), dlg.GetFlipNormals(), dlg.GetSplitNormals(),
            //    dlg.GetSplitEdgeAngle());
            forceRendering();
        }
    }
}


void MainWindow::computeConnectivity()
{
    ConnectivityDlg dlg;
    dlg.Set3DScene(ui->sceneWidget->get3DScene());

    if (dlg.exec())
    {
        // void ConnectivityDlg::getValues(bool& replaceSource, int& regionType, double* scalarRange, bool& fullScalarMode, double* point)
        bool replaceSource = false;
        int regionType = 0;
        double scalarRange[2];
        bool fullScalarMode = false;
        double p[3];
        dlg.getValues(replaceSource, regionType, scalarRange, fullScalarMode, p);
        ui->sceneWidget->get3DScene()->Connectivity(dlg.GetSelectedSurface(), replaceSource, regionType, scalarRange, fullScalarMode, p);
        forceRendering();
    }
}

void MainWindow::createPrimitive()
{
    CreatePrimitiveDlg dlg;
    dlg.setPickedPointPosition(ui->sceneWidget->get3DScene()->getLastPickedPoint());
    //dlg.Set3DScene(ui->sceneWidget->get3DScene();

    if (dlg.exec())
    {
        int primitiveType = dlg.getPrimitiveType();
        if (primitiveType == 0)
        {
            double r = 1;
            double phiStart = 0;
            double phiEnd = 180;
            double thetaStart = 0;
            double thetaEnd = 360;
            int phiRes = 50;
            int thetaRes = 50;
            double center[3];
            dlg.getSphereParameters(r, phiStart, phiEnd, thetaStart, thetaEnd, phiRes, thetaRes, center);
            ui->sceneWidget->get3DScene()->CreateSphere(r, thetaStart, thetaEnd, phiStart, phiEnd, thetaRes, phiRes, center);
        }
        else if (primitiveType == 1)
        {
            double size[3];
            double center[3];
            dlg.getCubeParameters(size, center);

            ui->sceneWidget->get3DScene()->CreateCube(size[0], size[1], size[2], center);
        }
        else if (primitiveType == 2)
        {
            double center[3];
            double r = 1;
            double h = 1;
            int res = 50;
            bool capped = false;
            dlg.getCylinderParameters(r, h, res, capped,center);

            ui->sceneWidget->get3DScene()->CreateCylinder(r, h, res, capped, center);
        }
        forceRendering();
    }

}

void MainWindow::exportScene()
{
    std::string extFilt = "Wavefront OBJ (*.obj);; RenderMan RIB (*.rib);; VRML 2.0 (*.wrl);; Single VTP (*.vtp);; All Files (*.*)";

    QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),
        tr(""), tr(extFilt.c_str()));

    if (fileName != "")
    {
        bool result = ui->sceneWidget->ExportSceneToFile(fileName.toStdString());

        if (!result)
            QMessageBox::warning(this, "Warning", "Could not export scene to file");
    }
}

void MainWindow::showOpenFileDialog()
{
    //std::string extFilt = "Surface File (*.vtk *.vtp *.stl *.txt *.ply);;All Files (*.*)";
    std::string extFilt = "Surface Files (*.vtk *.vtp *.stl *.ply *.obj *.wrl);;Raw X Y Z (*.txt *.dat *.pts *.txt *.bnd *.csv);;Landmarks (*.pp);;All Files (*.*)";

    QString fileName = QFileDialog::getOpenFileName(this, tr("Open file"), "", tr(extFilt.c_str()));
    //    "VTK Files (*.vtk) | STL Files (*.stl)");

    if (fileName != "")
    {
        QFile file(fileName);
        file.open(QIODevice::ReadOnly);

        if (!file.exists())
            return;

        openFile(fileName);
    }
}

void MainWindow::openFile(const QString& fileName)
{
    ui->sceneWidget->OpenFile(fileName.toStdString());
    if (mObjectPropsDlg)
        mObjectPropsDlg->FullUpdate();
    updateStatusBarMessage(tr("Read: ") + fileName);
    updateEnabledActions();
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
    updateEnabledActions();
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
        ManipulateObjectDlg dlg;
        dlg.Set3DScene(ui->sceneWidget->get3DScene());

        if (dlg.exec())
        {
            ui->sceneWidget->SetWidgetManipulateObject(dlg.GetSelectedSurface());
            CutWithPlaneAct->setChecked(true);
        }
    }
    else
    {
        ui->sceneWidget->SetWidgetManipulateObject(-1);
        CutWithPlaneAct->setChecked(false);
    }
    forceRendering();
}

void MainWindow::undoManipulate()
{
    ui->sceneWidget->UndoLastOperation();
    forceRendering();
}

void MainWindow::annotateWithSphere()
{
    if (ui->sceneWidget->GetSphereWidgetManipulateObject() == -1)
    {
        ManipulateObjectDlg dlg;
        dlg.Set3DScene(ui->sceneWidget->get3DScene());

        if (dlg.exec())
        {
            ui->sceneWidget->get3DScene()->SetScalarBarVisible(true);
            ui->sceneWidget->SetSphereWidgetManipulateObject(dlg.GetSelectedSurface());
            annotateWithSphereAct->setChecked(true);

            std::ostringstream ost;
            ost << "Marker value (set from menu): " << ui->sceneWidget->get3DScene()->GetMarkerValue()
                << " (size change: keys +/-)";
            ui->sceneWidget->get3DScene()->SetStatusText(ost.str(), true);
            updateStatusBarMessage(QString::fromStdString(ost.str()));
        }
    }
    else
    {
        ui->sceneWidget->SetSphereWidgetManipulateObject(-1);
        annotateWithSphereAct->setChecked(false);
    }
    forceRendering();
}

void MainWindow::updateStatusBarMessage(const QString& str)
{
    statusMessage = str;
    statusBar()->showMessage(statusMessage);
}

void MainWindow::createActions()
{
    openAct = new QAction(tr("&Open"), this);
    //openAct->setShortcuts(QKeySequence::New);
    openAct->setStatusTip(tr("Open file"));
    connect(openAct, &QAction::triggered, this, &MainWindow::showOpenFileDialog);

    saveFileAct = new QAction(tr("&Save"), this);
    //saveFileAct->setShortcuts(QKeySequence::New);
    saveFileAct->setStatusTip(tr("Save file"));
    connect(saveFileAct, &QAction::triggered, this, &MainWindow::showSaveFileDialog);

    exportSceneAct = new QAction(tr("&Export scene"), this);
    exportSceneAct->setStatusTip(tr("Export scene file"));
    connect(exportSceneAct, &QAction::triggered, this, &MainWindow::exportScene);

    ProcessComputeNormalsAct = new QAction(tr("&Compute normals"), this);
    ProcessComputeNormalsAct->setStatusTip(tr("Compute surface normals"));
    connect(ProcessComputeNormalsAct, &QAction::triggered, this, &MainWindow::processComputeNormals);

    visualizeNormalsAct = new QAction(tr("&Visualize normals"), this);
    visualizeNormalsAct->setStatusTip(tr("Visualize normals"));
    connect(visualizeNormalsAct, &QAction::triggered, this, &MainWindow::visualizeNormals);

    visualizeFeatureEdgesAct = new QAction(tr("Visualize &Feature edges"), this);
    visualizeFeatureEdgesAct->setStatusTip(tr("Visualize feature edges"));
    connect(visualizeFeatureEdgesAct, &QAction::triggered, this, &MainWindow::visualizeFeatureEdges);

    computeConnectivityAct = new QAction(tr("Compute &connectivity"), this);
    computeConnectivityAct->setStatusTip(tr("Compute mesh connectivity"));
    connect(computeConnectivityAct, &QAction::triggered, this, &MainWindow::computeConnectivity);

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

    undoManipulateAct = new QAction(tr("&Undo cut"), this);
    undoManipulateAct->setStatusTip(tr("Undo cut with plane"));
    connect(undoManipulateAct, &QAction::triggered, this, &MainWindow::undoManipulate);

    annotateWithSphereAct = new QAction(tr("Annotate with &Sphere"), this);
    annotateWithSphereAct->setCheckable(true);
    annotateWithSphereAct->setStatusTip(tr("Annotate with sphere"));
    connect(annotateWithSphereAct, &QAction::triggered, this, &MainWindow::annotateWithSphere);

    setMarkerSphereValueAct = new QAction(tr("Set Marker Sphere &Value"), this);
    setMarkerSphereValueAct-> setStatusTip(tr("Set sphere marker value"));
    connect(setMarkerSphereValueAct, &QAction::triggered, this, &MainWindow::setMarkerSphereValue);

    showAxesAct = new QAction(tr("Show &Axes"), this);
    showAxesAct->setCheckable(true);
    showAxesAct->setStatusTip(tr("Show coordinate axes"));
    connect(showAxesAct, &QAction::triggered, this, &MainWindow::showAxes);

    showScalarBarAct = new QAction(tr("Show &Scalar Bar"), this);
    showScalarBarAct->setCheckable(true);
    showScalarBarAct->setStatusTip(tr("Show Scalar Bar"));
    connect(showScalarBarAct, &QAction::triggered, this, &MainWindow::showScalarBar);

    createPrimitiveAct = new QAction(tr("Create &Primitive"), this);
    createPrimitiveAct->setStatusTip(tr("Create geometric primitive"));
    connect(createPrimitiveAct, &QAction::triggered, this, &MainWindow::createPrimitive);
}

void MainWindow::createMenus()
{
    fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(openAct);
    fileMenu->addAction(saveFileAct);
    fileMenu->addAction(exportSceneAct);

    fileMenu = menuBar()->addMenu(tr("&View"));
    fileMenu->addAction(showAxesAct);
    fileMenu->addAction(showScalarBarAct);

    fileMenu = menuBar()->addMenu(tr("&Create"));
    fileMenu->addAction(createPrimitiveAct);

    fileMenu = menuBar()->addMenu(tr("&Options"));
    fileMenu->addAction(optionsObjectPropAct);
    fileMenu->addAction(ChooseBackgroundColorAct);

    fileMenu = menuBar()->addMenu(tr("&Manipulate"));
    fileMenu->addAction(CutWithPlaneAct);
    fileMenu->addAction(undoManipulateAct);
    fileMenu->addAction(annotateWithSphereAct);
    fileMenu->addAction(setMarkerSphereValueAct);

    fileMenu = menuBar()->addMenu(tr("&Process"));
    fileMenu->addAction(ProcessComputeNormalsAct);
    fileMenu->addAction(visualizeNormalsAct);
    fileMenu->addAction(visualizeFeatureEdgesAct);
    fileMenu->addAction(computeConnectivityAct);

    updateEnabledActions();
}

void MainWindow::updateEnabledActions()
{
    bool anyObjects = (ui->sceneWidget->get3DScene()->GetNumberOfSurfaces() > 0);
    saveFileAct->setEnabled(anyObjects);
    ProcessComputeNormalsAct->setEnabled(anyObjects);
    visualizeNormalsAct->setEnabled(anyObjects);
    visualizeFeatureEdgesAct->setEnabled(anyObjects);
    computeConnectivityAct->setEnabled(anyObjects);
    CutWithPlaneAct->setEnabled(anyObjects);
    undoManipulateAct->setEnabled(anyObjects);
    annotateWithSphereAct->setEnabled(anyObjects);
    setMarkerSphereValueAct->setEnabled(anyObjects);
    optionsObjectPropAct->setEnabled(anyObjects);
    showScalarBarAct->setChecked(ui->sceneWidget->get3DScene()->GetScalarBarVisible());
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
