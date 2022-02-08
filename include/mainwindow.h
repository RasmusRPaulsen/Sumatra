#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void showOpenFileDialog();

protected:
    void openFile(const QString& fileName);

private:
    Ui::MainWindow *ui;

private:
    void createActions();
    void createMenus();

    QMenu* fileMenu;
    QAction* openAct;
};

#endif // MAINWINDOW_H
