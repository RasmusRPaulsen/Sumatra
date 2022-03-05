#ifndef SCENEWIDGET_H
#define SCENEWIDGET_H

#include <QVTKOpenGLNativeWidget.h>
#include <QtGui>
#include <QString>
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

    void ForceRender();

    C3DScene* get3DScene();

    //! Set that that the widgets should work on the object with the given id
    void SetWidgetManipulateObject(int objId);

    //! Get the object id that the widgets are currently manipulating (-1 in case of none)
    int GetWidgetManipulateObject() const;

    //! Set that that the widgets should work on the object with the given id
    void SetSphereWidgetManipulateObject(int objId);

    //! Get the object id that the widgets are currently manipulating (-1 in case of none)
    int GetSphereWidgetManipulateObject() const;

    //! Space is pressed so if a plane widget is active it should cut
    bool CutWithPlane();

    //! M is pressed so if a plane widget is active it should create a mirror image
    bool MirrorWithPlane();

    //! Space is pressed so if a sphere widget is active it should mark
    bool MarkWithSphere();

    //! Undo last operation
    void UndoLastOperation();

    //! Is undo avaiable
    bool UndoAvailable() const;

public slots:
    //! Zoom to the extent of the data set in the scene
    // void zoomToExtent();

signals:
    void updateStatusMessage(const QString& str);

private:
    virtual ~SceneWidget();

    C3DScene* m3DScene = NULL;

    CSumatraSettings* mSettings = NULL;

    vtkSmartPointer<vtkInteractorStyleTrackballActor> mStyleActor;
    vtkSmartPointer<vtkInteractorStyleTrackballCamera> mStyleCamera;

    // Which object is currently being cut with a plane 
    int m_WidgetManipulateObject = -1;

    int m_SphereWidgetManipulateObject = -1;

    void SpacePressed();

    void changeMarkerSphereSize(bool increase);

protected:
    void keyPressEvent(QKeyEvent*);
};

#endif // SCENEWIDGET_H
