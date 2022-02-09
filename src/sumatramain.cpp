#include <QtWidgets/QApplication>
#include <MainWindow.h>
#include <QSurfaceFormat>
#include <QVTKOpenGLNativeWidget.h>
#include <qmessagebox.h>

int main(int argc, char *argv[]) {
	QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());

	QApplication app(argc, argv);

	MainWindow mainWindow;
	mainWindow.show();

	// TODO: Open all files
	if (app.arguments().length() > 1)
	{
		mainWindow.openFile(app.arguments().at(1).toLocal8Bit().constData());
	}
	// QMessageBox::information(NULL, "Arguments:", app.arguments().join(" | "));

	return app.exec();
}