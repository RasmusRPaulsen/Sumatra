#include "scenewidget.h"

#include <qmessagebox.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOBJExporter.h>
#include <vtkRIBExporter.h>
#include <vtkVRMLExporter.h>
#include <vtkSingleVTPExporter.h>

#include <sstream>

#include "GeneralUtils.h"

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

bool SceneWidget::ExportSceneToFile(const std::string& fname)
{
    std::string ext = CGeneralUtils::GetExtensionFromFilename(fname);
    std::string prefix = CGeneralUtils::StripExtensionFromFilename(fname);

    bool result = true;
    if (ext == "obj")
    {
        vtkOBJExporter* writer = vtkOBJExporter::New();
        writer->SetRenderWindow(renderWindow());
        writer->SetFilePrefix(prefix.c_str());
        writer->Write();
        writer->Delete();
    }
    else if (ext == "rib")
    {
        vtkRIBExporter* writer = vtkRIBExporter::New();
        writer->SetRenderWindow(renderWindow());
        writer->SetFilePrefix(prefix.c_str());
        writer->Write();
        writer->Delete();
    }
    else if (ext == "wrl")
    {
        vtkVRMLExporter* writer = vtkVRMLExporter::New();
        writer->SetRenderWindow(renderWindow());
        writer->SetFileName(fname.c_str());
        writer->Write();
        writer->Delete();
    }
    else if (ext == "vtp")
    {
        vtkSingleVTPExporter* writer = vtkSingleVTPExporter::New();
        writer->SetRenderWindow(renderWindow());
        writer->SetFilePrefix(prefix.c_str());
        writer->Write();
        writer->Delete();
    }
    else
    {
        result = false;
    }

    return result;
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


void SceneWidget::SpacePressed()
{
    if (!CutWithPlane())
    {
        if (!MarkWithSphere())
        {
            // For strange reasons these methods do not give the correct x,y pos for pointpicking
            //QSize wSize = this->size();
            //QPoint globalPos = QCursor::pos();
            //QPoint localPos = this->mapFromGlobal(globalPos);
            
            int ePos[2];
            renderWindow()->GetInteractor()->GetEventPosition(ePos);

            // std::string text = m3DScene->PickPointAndGetText(localPos.x(), wSize.height() - localPos.y() - 1);
            std::string text = m3DScene->PickPointAndGetText(ePos[0], ePos[1]);

            QString statusMsg = "Pick position (" + QString::number(ePos[0]) + ", " + QString::number(ePos[1]) + ") " +
                QString::fromStdString(text);

            //QString statusMsg = "Globalpos (" + QString::number(globalPos.x()) + ", " + QString::number(globalPos.y()) + ") " +
            //    "LocalPos (" + QString::number(localPos.x()) + ", " + QString::number(localPos.y()) + ") " +
            //    "EventPos (" + QString::number(ePos[0]) + ", " + QString::number(ePos[1]) + ") " +
            //    QString::fromStdString(text);

            emit updateStatusMessage(statusMsg);
        }
    }
    ForceRender();
}

void SceneWidget::changeMarkerSphereSize(bool increase)
{
    double scale = 1.1;
    if (!increase)
        scale = 0.9;

    m3DScene->SetSphereWidgetRadius(m3DScene->GetSphereWidgetRadius() * scale);

    std::ostringstream ost;
    ost << "Marker value: " << m3DScene->GetMarkerValue() <<
        " radius: " << m3DScene->GetSphereWidgetRadius();
    m3DScene->SetStatusText(ost.str(), true);
    emit updateStatusMessage(QString::fromStdString(ost.str()));
    ForceRender();
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
    else if (event->key() == Qt::Key_Space)
    {
        SpacePressed();
    }
    else if (event->key() == Qt::Key_Plus)
    {
        if (m_SphereWidgetManipulateObject != -1)
            changeMarkerSphereSize(true);
    }
    else if (event->key() == Qt::Key_Minus)
    {
        if (m_SphereWidgetManipulateObject != -1)
            changeMarkerSphereSize(false);
    }
}
