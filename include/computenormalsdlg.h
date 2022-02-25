#ifndef COMPUTENORMALSDLG_H
#define COMPUTENORMALSDLG_H

#include <QDialog>
#include "3DScene.h"

namespace Ui {
class ComputeNormalsDlg;
}

class ComputeNormalsDlg : public QDialog
{
    Q_OBJECT

public:
    explicit ComputeNormalsDlg(QWidget *parent = nullptr);
    ~ComputeNormalsDlg();

    // void Setup();

    void Set3DScene(C3DScene* sc);

    void PopulateObjectCB();

    int GetSelectedSurface();

    bool GetReplaceSource();

    bool GetFlipNormals();

    bool GetSplitNormals();

    double GetSplitEdgeAngle();

private:
    Ui::ComputeNormalsDlg* ui;

    C3DScene* m3DScene = NULL;
};

#endif // COMPUTENORMALSDLG_H
