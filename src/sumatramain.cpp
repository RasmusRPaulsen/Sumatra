#include <QtWidgets/QApplication>
#include <MainWindow.h>
#include <QSurfaceFormat>
#include <QVTKOpenGLNativeWidget.h>


int main(int argc, char *argv[]) {
	QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());

	QApplication app(argc, argv);

	MainWindow mainWindow;
	mainWindow.show();

	return app.exec();
}