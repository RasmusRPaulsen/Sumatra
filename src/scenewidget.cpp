#include "scenewidget.h"

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
