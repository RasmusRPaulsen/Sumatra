#ifndef MANIPULATEOBJECTDLG_H
#define MANIPULATEOBJECTDLG_H

#include <QDialog>
#include "3DScene.h"


namespace Ui {
class ManipulateObjectDlg;
}

class ManipulateObjectDlg : public QDialog
{
    Q_OBJECT

public:
    explicit ManipulateObjectDlg(QWidget *parent = nullptr);
    ~ManipulateObjectDlg();

    void Set3DScene(C3DScene* sc);

    void PopulateObjectCB();

    int GetSelectedSurface();

private:
    Ui::ManipulateObjectDlg *ui;

    C3DScene* m3DScene = NULL;
};

#endif // MANIPULATEOBJECTDLG_H
