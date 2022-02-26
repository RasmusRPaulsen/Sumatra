#include "scenewidget.h"

#include <qmessagebox.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>

SceneWidget::SceneWidget(QWidget* parent)
    : QVTKOpenGLNativeWidget(parent)
{
    vtkNew<vtkGenericOpenGLRenderWindow> window;
    setRenderWindow(window.Get());
}

SceneWidget::~SceneWidget()
{
    if (m3DScene)
        delete m3DScene;
}

void SceneWidget::ForceRender()
{
    renderWindow()->Render();
}

void SceneWidget::Setup(CSumatraSettings* Settings)
{
    mSettings = Settings;

    m3DScene = new C3DScene(Settings);
    m3DScene->Init();

    mStyleActor = vtkSmartPointer<vtkInteractorStyleTrackballActor>::New();
    mStyleCamera = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();

    renderWindow()->AddRenderer(m3DScene->GetRenderer());
}

void SceneWidget::OpenFile(const std::string& fname)
{
    m3DScene->ReadFile(fname);
    renderWindow()->Render();
}
C3DScene* SceneWidget::get3DScene()
{
    return m3DScene;
}
void SceneWidget::SetWidgetManipulateObject(int objId)
{
    // turn on widget
    if (m_WidgetManipulateObject == -1 && objId >= 0)
    {
        m3DScene->InitialisePlaneWidget(renderWindow()->GetInteractor(), objId);
        m_WidgetManipulateObject = objId;
    }
    else if (m_WidgetManipulateObject != -1)
        // turn off widget if this function is called and the widget is already in action
    {
        m3DScene->DeletePlaneWidget();
        m_WidgetManipulateObject = -1;
    }
}

int SceneWidget::GetWidgetManipulateObject() const
{
    return m_WidgetManipulateObject;
}

void SceneWidget::SetSphereWidgetManipulateObject(int objId)
{
    // turn on widget
    if (m_SphereWidgetManipulateObject == -1 && objId >= 0)
    {
        m3DScene->InitialiseSphereWidget(renderWindow()->GetInteractor(), objId);
        m_SphereWidgetManipulateObject = objId;
        m3DScene->StartMarkWithSphere(m_SphereWidgetManipulateObject);
    }
    else if (m_SphereWidgetManipulateObject != -1)
        // turn off widget if this function is called and the widget is already in action
    {
        m3DScene->DeleteSphereWidget();
        m_SphereWidgetManipulateObject = -1;
        m3DScene->EndMarkWithSphere(m_SphereWidgetManipulateObject);
    }
}

int SceneWidget::GetSphereWidgetManipulateObject() const
{
    return m_SphereWidgetManipulateObject;
}

bool SceneWidget::CutWithPlane()
{
    if (m_WidgetManipulateObject == -1)
        return false;

    m3DScene->CutWithPlane(m_WidgetManipulateObject);
    return true;
}

bool SceneWidget::MirrorWithPlane()
{
    if (m_WidgetManipulateObject == -1)
        return false;

    m3DScene->MirrorWithPlane(m_WidgetManipulateObject);
    return true;
}

bool SceneWidget::MarkWithSphere()
{
    if (m_SphereWidgetManipulateObject == -1)
        return false;

    m3DScene->MarkWithSphere(m_SphereWidgetManipulateObject);
    return true;
}

void SceneWidget::UndoLastOperation()
{
    if (m_WidgetManipulateObject == -1)
        return;

    m3DScene->UndoLastOperation(m_WidgetManipulateObject);
}

bool SceneWidget::UndoAvailable() const
{
    if (m_WidgetManipulateObject == -1)
        return false;

    return m3DScene->GetUndoAvaible(m_WidgetManipulateObject);
}


void SceneWidget::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_A)
    {
        renderWindow()->GetInteractor()->SetInteractorStyle(mStyleActor);
    }
    else if (event->key() == Qt::Key_C)
    {
        renderWindow()->GetInteractor()->SetInteractorStyle(mStyleCamera);
    }
}
