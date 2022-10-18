#ifndef SUBDIVIDEDLG_H
#define SUBDIVIDEDLG_H

#include <QDialog>
#include "3DScene.h"

namespace Ui {
class SubdivideDlg;
}

class SubdivideDlg : public QDialog
{
    Q_OBJECT

public:
    explicit SubdivideDlg(QWidget *parent = nullptr);
    ~SubdivideDlg();
    void Set3DScene(C3DScene* sc);

    void PopulateObjectCB();

    int GetSelectedSurface();

    bool GetReplaceSource();

    double GetNumberOfSubdivisions();

    // 0: Linear, 1: Butterfly, 2: Loop 
    int GetSubdivisionType();


private:
    Ui::SubdivideDlg *ui;

    C3DScene* m3DScene = NULL;
};

#endif // SUBDIVIDEDLG_H
