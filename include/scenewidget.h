#ifndef SCENEWIDGET_H
#define SCENEWIDGET_H

#include <QVTKOpenGLNativeWidget.h>
#include <QtGui>
#include <vtkSmartPointer.h>
#include <vtkInteractorStyleTrackballActor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <string>
#include "3DScene.h"
#include "SumatraSettings.h"

class SceneWidget : public QVTKOpenGLNativeWidget {
    Q_OBJECT
public:
    explicit SceneWidget(QWidget* parent = nullptr);

    void Setup(CSumatraSettings* Settings);

    void OpenFile(const std::string& fname);

    C3DScene* get3DScene();

    // void SetSettings(CSumatraSettings* Settings);

public slots:
    //! Zoom to the extent of the data set in the scene
    // void zoomToExtent();

private:
    virtual ~SceneWidget();

    C3DScene* m3DScene = NULL;

    CSumatraSettings* mSettings = NULL;

    vtkSmartPointer<vtkInteractorStyleTrackballActor> mStyleActor;
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> mStyleCamera;

protected:
    void keyPressEvent(QKeyEvent*);
};

#endif // SCENEWIDGET_H
