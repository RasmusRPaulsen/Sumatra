#ifndef COMPUTENORMALSDLG_H
#define COMPUTENORMALSDLG_H

#include <QDialog>

namespace Ui {
class ComputeNormalsDlg;
}

class ComputeNormalsDlg : public QDialog
{
    Q_OBJECT

public:
    explicit ComputeNormalsDlg(QWidget *parent = nullptr);
    ~ComputeNormalsDlg();

    void Setup();

    int getSelectedSurface();

    bool replaceSource();


private:
    Ui::ComputeNormalsDlg *ui;
};

#endif // COMPUTENORMALSDLG_H
