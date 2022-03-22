#ifndef CREATEPRIMITIVEDLG_H
#define CREATEPRIMITIVEDLG_H

#include <QDialog>

namespace Ui {
class CreatePrimitiveDlg;
}

class CreatePrimitiveDlg : public QDialog
{
    Q_OBJECT

public:
    explicit CreatePrimitiveDlg(QWidget *parent = nullptr);
    ~CreatePrimitiveDlg();

private:
    Ui::CreatePrimitiveDlg *ui;
};

#endif // CREATEPRIMITIVEDLG_H
