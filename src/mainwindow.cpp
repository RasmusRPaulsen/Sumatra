#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFile>
#include <QFileDialog>

#include <vtkDataSetReader.h>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    createActions();
    createMenus();

    QString message = tr("Default status bar message");
    statusBar()->showMessage(message);

    setWindowTitle(tr("Sumatra (Qt)"));
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::showOpenFileDialog()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open file"), "",
        "VTK Files (*.vtk)");

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
    ui->sceneWidget->removeDataSet();

    // Create reader
    vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
    reader->SetFileName(fileName.toStdString().c_str());

    // Read the file
    reader->Update();

    // Add data set to 3D view
    vtkSmartPointer<vtkDataSet> dataSet = reader->GetOutput();
    if (dataSet != nullptr) {
        ui->sceneWidget->addDataSet(reader->GetOutput());
    }
}

void MainWindow::createActions()
{
    openAct = new QAction(tr("&OpenTest"), this);
    openAct->setShortcuts(QKeySequence::New);
    openAct->setStatusTip(tr("Open a file"));
    connect(openAct, &QAction::triggered, this, &MainWindow::showOpenFileDialog);
}

void MainWindow::createMenus()
{
    fileMenu = menuBar()->addMenu(tr("&FileTest"));
    fileMenu->addAction(openAct);
}
