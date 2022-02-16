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

    m3DScene = new C3DScene(mSettings);
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
//
//void SceneWidget::SetSettings(CSumatraSettings* Settings)
//{
//    mSettings = Settings;
//}

// QMessageBox::information(this, "KeyPress", tr("Keypress"));

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
