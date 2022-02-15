#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFile>
#include <QFileDialog>
#include <qmessagebox.h>
#include <qmimedata.h>
#include "computenormalsdlg.h"
#include "objectpropertiesdlg.h"
#include "SumatraSettings.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    CSumatraSettings mSettings;
    mSettings.ReadSettings();
    ui->sceneWidget->Setup(&mSettings);

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
}

void MainWindow::processComputeNormals()
{
    ComputeNormalsDlg dlg;
    dlg.Setup();
    if (dlg.exec())
    {
        int selected = dlg.getSelectedSurface();
        QMessageBox::information(this, "Hello", tr("Selected surface: ") + QString::number(selected));


        //bool Opt1, Opt2, Opt3;
        //myDialog->GetOptions(Opt1, Opt2, Opt3);
        //DoSomethingWithThoseBooleans(Opt1, Opt2, Opt3);
    }

}

void MainWindow::showObjectProperties()
{
    ObjectPropertiesDlg dlg;
    if (dlg.exec())
    {
    }
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
}

void MainWindow::createMenus()
{
    fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(openAct);

    fileMenu = menuBar()->addMenu(tr("&Process"));
    fileMenu->addAction(ProcessComputeNormalsAct);

    fileMenu = menuBar()->addMenu(tr("&Options"));
    fileMenu->addAction(optionsObjectPropAct);
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
