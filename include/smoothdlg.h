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

private:
    Ui::SmoothDlg *ui;

    C3DScene* m3DScene = NULL;
};

#endif // SMOOTHDLG_H
