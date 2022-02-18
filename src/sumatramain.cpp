#include <QtWidgets/QApplication>
#include <MainWindow.h>
#include <QSurfaceFormat>
#include <QVTKOpenGLNativeWidget.h>
#include <qmessagebox.h>
#include <qicon.h>

int main(int argc, char *argv[]) {
	QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());

	QApplication app(argc, argv);

	QIcon icon(":/icons/Sumatra.png");
	// QIcon icon(QString("C:\\Users\\rapa\\Documents\\src\\Cpp\\Sumatra\\resources\\icons\\Sumatra.png"));

	MainWindow mainWindow;
	mainWindow.show();
	mainWindow.setWindowIcon(icon);

	// setWindowIcon(QtGui.QIcon('file path'))

	if (app.arguments().length() > 1)
	{
		for (int i = 1; i < app.arguments().length(); i++)
			mainWindow.openFile(app.arguments().at(i));
	}

	return app.exec();
}
