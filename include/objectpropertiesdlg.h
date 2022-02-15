#ifndef OBJECTPROPERTIESDLG_H
#define OBJECTPROPERTIESDLG_H

#include <QDialog>

namespace Ui {
class ObjectPropertiesDlg;
}

class ObjectPropertiesDlg : public QDialog
{
    Q_OBJECT

public:
    explicit ObjectPropertiesDlg(QWidget *parent = nullptr);
    ~ObjectPropertiesDlg();

private:
    Ui::ObjectPropertiesDlg *ui;
};

#endif // OBJECTPROPERTIESDLG_H
