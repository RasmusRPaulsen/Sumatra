#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#define SUMATRA_VERSION "0.0.1"

#include <QMainWindow>
#include <qdragenterevent>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void openFile(const QString& fileName);

private slots:
    void showOpenFileDialog();
    void processComputeNormals();

protected:
    void dropEvent(QDropEvent* event);
    void dragEnterEvent(QDragEnterEvent* event);

private:
    Ui::MainWindow *ui;

private:
    void createActions();
    void createMenus();

    QMenu* fileMenu;
    QAction* openAct;
    QAction* ProcessComputeNormalsAct;
};

#endif // MAINWINDOW_H
