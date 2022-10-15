#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#define SUMATRA_VERSION "0.2.1"

#include <QMainWindow>
#include <QDragEnterEvent>
#include <QString>
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
    int selectObjectDialog();
    void showOpenFileDialog();
    void processComputeNormals();
    void processDecimate();
    void showObjectProperties();
    void forceRendering();
    void ChooseBackgroundColor();
    void CutWithPlane();
    void undoManipulate();
    void annotateWithSphere();
    void updateStatusBarMessage(const QString& str);
    void showSaveFileDialog();
    void showAxes();
    void showScalarBar();
    void setMarkerSphereValue();
    void visualizeNormals();
    void visualizeFeatureEdges();
    void computeConnectivity();
    void createPrimitive();
    void exportScene();

protected:
    void dropEvent(QDropEvent* event);
    void dragEnterEvent(QDragEnterEvent* event);

private:
    Ui::MainWindow *ui;

    CSumatraSettings *mSettings = NULL;

private:
    void createActions();
    void createMenus();

    void updateEnabledActions();

    ObjectPropertiesDlg* mObjectPropsDlg = NULL;

    // Status bare message
    QString statusMessage = tr("Welcome to the Surface Manipulation and Transformation Toolkit (Sumatra)");

    QMenu* fileMenu;
    QMenu* optionsMenu;
    QAction* openAct;
    QAction* saveFileAct;
    QAction* exportSceneAct;
    QAction* ProcessComputeNormalsAct;
    QAction* ProcessDecimateAct;
    QAction* optionsObjectPropAct;
    QAction* ChooseBackgroundColorAct;
    QAction* CutWithPlaneAct;
    QAction* undoManipulateAct;
    QAction* annotateWithSphereAct;
    QAction* showAxesAct;
    QAction* showScalarBarAct;
    QAction* setMarkerSphereValueAct;
    QAction* visualizeNormalsAct;
    QAction* visualizeFeatureEdgesAct;
    QAction* computeConnectivityAct;
    QAction* createPrimitiveAct;
};

#endif // MAINWINDOW_H
