#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFile>
#include <QFileDialog>
#include <qmessagebox.h>
#include <qmimedata.h>
#include "computenormalsdlg.h"
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

void MainWindow::createActions()
{
    openAct = new QAction(tr("&Open"), this);
    openAct->setShortcuts(QKeySequence::New);
    openAct->setStatusTip(tr("Open file"));
    connect(openAct, &QAction::triggered, this, &MainWindow::showOpenFileDialog);

    ProcessComputeNormalsAct = new QAction(tr("&Compute normals"), this);
    // ProcessComputeNormalsAct->setShortcuts(QKeySequence::New);
    ProcessComputeNormalsAct->setStatusTip(tr("Compute sutface normals"));
    connect(ProcessComputeNormalsAct, &QAction::triggered, this, &MainWindow::processComputeNormals);
}

void MainWindow::createMenus()
{
    fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(openAct);

    fileMenu = menuBar()->addMenu(tr("&Process"));
    fileMenu->addAction(ProcessComputeNormalsAct);
}

void MainWindow::dragEnterEvent(QDragEnterEvent* event)
{
    event->acceptProposedAction();
}


void MainWindow::dropEvent(QDropEvent* event)
{
    //QMessageBox::information(
    //    this, "Hello",
    //    "Test\nbox");

    const QMimeData* mimeData = event->mimeData();

    // check for our needed mime type, here a file or a list of files
    if (mimeData->hasUrls())
    {
        //QStringList pathList;
        QList<QUrl> urlList = mimeData->urls();

        // extract the local paths of the files
        for (int i = 0; i < urlList.size(); i++)
        {
            // pathList.append(urlList.at(i).toLocalFile());
            // call a function to open the files
            openFile(urlList.at(i).toLocalFile());
        }
    }    
}
