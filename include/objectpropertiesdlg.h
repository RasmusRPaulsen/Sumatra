#ifndef OBJECTPROPERTIESDLG_H
#define OBJECTPROPERTIESDLG_H

#include <QDialog>
#include "3DScene.h"

namespace Ui {
class ObjectPropertiesDlg;
}

class ObjectPropertiesDlg : public QDialog
{
    Q_OBJECT

public:
    explicit ObjectPropertiesDlg(QWidget *parent = nullptr);
    ~ObjectPropertiesDlg();

    void Set3DScene(C3DScene* sc);

    void FullUpdate();

    void PopulateObjectCB();

    void UpdateAllSceneData();

public slots:
    void selectionChanged(int);
    void OnRenderingTypePoints(bool);
    void OnRenderingTypeWireframe(bool);
    void OnRenderingTypeSurface(bool);
    void OnRenderingTypeHidden(bool);
    void OnRenderingPointsAsSpheres(bool);
    void OnRenderingLinesAsTubes(bool);
    void OnRenderingFrontFaceCulling(bool);
    void OnRenderingBackFaceCulling(bool);
    void OnRenderingVisibleEdges(bool);
    void OnPointSizeSpin(int);
    void OnLineWidthSpin(int);
    void OnOpacitySpin(double);
    void OnAmbientSpin(double);
    void OnDiffuseSpin(double);
    void OnSpecularSpin(double);
    void OnSpecularPowerSpin(double);
    void OnShowScalarsChk(bool);
    void OnChooseColorBtn();
    void OnDeleteObjectBtn();
    void OnResetTransformBtn();
    void OnRemoveScalarsBtn();
    void OnRemoveNormalsBtn();
    void OnSetMinScalarRange(double);
    void OnSetMaxScalarRange(double);
    void OnSetScalarRangeForAll();

signals:
    void valueChanged();

private:
    Ui::ObjectPropertiesDlg *ui;

    C3DScene *m3DScene = NULL;

};

#endif // OBJECTPROPERTIESDLG_H
