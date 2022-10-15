#ifndef DECIMATEDLG_H
#define DECIMATEDLG_H

#include <QDialog>
#include "3DScene.h"

namespace Ui {
class DecimateDlg;
}

class DecimateDlg : public QDialog
{
    Q_OBJECT

public:
    explicit DecimateDlg(QWidget *parent = nullptr);
    ~DecimateDlg();

    void Set3DScene(C3DScene* sc);

    void PopulateObjectCB();

    int GetSelectedSurface();

    bool GetReplaceSource();


private:
    Ui::DecimateDlg *ui;

    C3DScene* m3DScene = NULL;
};

#endif // DECIMATEDLG_H
