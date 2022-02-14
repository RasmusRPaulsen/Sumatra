#ifndef SCENEWIDGET_H
#define SCENEWIDGET_H

#include <QVTKOpenGLNativeWidget.h>
#include <QtGui>
#include <vtkSmartPointer.h>
#include <vtkInteractorStyleTrackballActor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <string>
#include "3DScene.h"

class SceneWidget : public QVTKOpenGLNativeWidget {
    Q_OBJECT
public:
    explicit SceneWidget(QWidget* parent = nullptr);

    void OpenFile(const std::string& fname);

public slots:
    //! Zoom to the extent of the data set in the scene
    // void zoomToExtent();

private:
    virtual ~SceneWidget();

    C3DScene* m3DScene = NULL;

    vtkSmartPointer<vtkInteractorStyleTrackballActor> mStyleActor;
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> mStyleCamera;

protected:
    void keyPressEvent(QKeyEvent*);
};

#endif // SCENEWIDGET_H
