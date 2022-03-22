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

    void setPickedPointPosition(double* p);

    //! 0: Sphere, 1: cube
    int getPrimitiveType();

	void getSphereParameters(double &r, double &phiStart, double &phiEnd, double &thetaStart, double &thetaEnd,
		int &phiRes, int &thetaRes, double *center);

    void getCubeParameters(double *size, double* center);

private:
    Ui::CreatePrimitiveDlg *ui;
};

#endif // CREATEPRIMITIVEDLG_H
