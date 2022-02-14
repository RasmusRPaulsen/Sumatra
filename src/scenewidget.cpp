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

    m3DScene = new C3DScene();
    m3DScene->Init();

    mStyleActor = vtkSmartPointer<vtkInteractorStyleTrackballActor>::New();
    mStyleCamera = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();

    renderWindow()->AddRenderer(m3DScene->GetRenderer());
}

SceneWidget::~SceneWidget()
{
    if (m3DScene)
        delete m3DScene;
}

void SceneWidget::OpenFile(const std::string& fname)
{
    m3DScene->ReadFile(fname);
    renderWindow()->Render();
}

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
