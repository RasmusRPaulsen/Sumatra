#ifndef SMOOTHDLG_H
#define SMOOTHDLG_H

#include <QDialog>
#include "3DScene.h"

namespace Ui {
class SmoothDlg;
}

class SmoothDlg : public QDialog
{
    Q_OBJECT

public:
    explicit SmoothDlg(QWidget *parent = nullptr);
    ~SmoothDlg();
    void Set3DScene(C3DScene* sc);

    void PopulateObjectCB();

    int GetSelectedSurface();

    bool GetReplaceSource();

    // 0: laplacian, 1: Sinc, 2: constrained
    int GetSmoothType();

    int GetIterations();

    double GetRelaxationFactor();

    bool GetBoundarySmoothing();

    bool GetFeatureEdgeSmoothing();

    double GetFeatureAngle();

    bool GetGenerateErrorScalars();


private:
    Ui::SmoothDlg *ui;

    C3DScene* m3DScene = NULL;
};

#endif // SMOOTHDLG_H
