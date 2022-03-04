#ifndef SAVEFILEDLG_H
#define SAVEFILEDLG_H

#include <QDialog>

#include "3DScene.h"

namespace Ui {
class SaveFileDlg;
}

class SaveFileDlg : public QDialog
{
    Q_OBJECT

public:
    explicit SaveFileDlg(QWidget *parent = nullptr);
    ~SaveFileDlg();

    void Set3DScene(C3DScene* sc);

    void PopulateObjectCB();

    int GetSelectedSurface();

    bool getApplyActorTransform();

    bool getWriteScalars();

    bool getWriteNormals();

    bool getWriteAsAscii();


private:
    Ui::SaveFileDlg *ui;

    C3DScene* m3DScene = NULL;
};

#endif // SAVEFILEDLG_H
