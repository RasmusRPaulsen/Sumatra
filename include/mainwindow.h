#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#define SUMATRA_VERSION "0.0.1"

#include <QMainWindow>
#include <qdragenterevent>
#include "objectpropertiesdlg.h"
#include "SumatraSettings.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void openFile(const QString& fileName);

private slots:
    void showOpenFileDialog();
    void processComputeNormals();
    void showObjectProperties();
    void forceRendering();
    void ChooseBackgroundColor();
    void CutWithPlane();
    void undoManipulate();
    void annotateWithSphere();

protected:
    void dropEvent(QDropEvent* event);
    void dragEnterEvent(QDragEnterEvent* event);

private:
    Ui::MainWindow *ui;

    CSumatraSettings *mSettings = NULL;

private:
    void createActions();
    void createMenus();

    ObjectPropertiesDlg* mObjectPropsDlg = NULL;

    QMenu* fileMenu;
    QMenu* optionsMenu;
    QAction* openAct;
    QAction* ProcessComputeNormalsAct;
    QAction* optionsObjectPropAct;
    QAction* ChooseBackgroundColorAct;
    QAction* CutWithPlaneAct;
    QAction* undoManipulateAct;
    QAction* annotateWithSphereAct;
};

#endif // MAINWINDOW_H
