#ifndef CONNECTIVITYDLG_H
#define CONNECTIVITYDLG_H

#include <QDialog>
#include "3DScene.h"

namespace Ui {
class ConnectivityDlg;
}

class ConnectivityDlg : public QDialog
{
    Q_OBJECT

public:
    explicit ConnectivityDlg(QWidget *parent = nullptr);
    ~ConnectivityDlg();

    void Set3DScene(C3DScene* sc);

    void setPickedPointPosition(double* p);

    void PopulateObjectCB();

    int GetSelectedSurface();

    void getValues(bool& replaceSource, int& regionType, double* scalarRange, bool& fullScalarMode,
        double* point);

private:
    C3DScene* m3DScene = NULL;

    Ui::ConnectivityDlg *ui;
};

#endif // CONNECTIVITYDLG_H
