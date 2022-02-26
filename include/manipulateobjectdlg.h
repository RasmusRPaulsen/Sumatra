#ifndef MANIPULATEOBJECTDLG_H
#define MANIPULATEOBJECTDLG_H

#include <QDialog>

namespace Ui {
class ManipulateObjectDlg;
}

class ManipulateObjectDlg : public QDialog
{
    Q_OBJECT

public:
    explicit ManipulateObjectDlg(QWidget *parent = nullptr);
    ~ManipulateObjectDlg();

private:
    Ui::ManipulateObjectDlg *ui;
};

#endif // MANIPULATEOBJECTDLG_H
