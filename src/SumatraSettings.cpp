#include "SumatraSettings.h"

#include <qdir.h>
#include <qcoreapplication.h>
#include <qmessagebox.h>

CSumatraSettings::CSumatraSettings()
{
}

CSumatraSettings::~CSumatraSettings()
{
}

bool CSumatraSettings::ReadSettings()
{
    // QString dir1 = QDir::currentPath();
        
    QString setdir = QCoreApplication::applicationDirPath();
    // QMessageBox::information(NULL, "Application dir", dir2);


    return true;
}
