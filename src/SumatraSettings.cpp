#include "SumatraSettings.h"

#include <qdir.h>
#include <qcoreapplication.h>
#include <qmessagebox.h>
#include <qfile.h>
#include <QJsonDocument>
#include <qjsonarray.h>
#include <qtextstream.h>

CSumatraSettings::CSumatraSettings()
{
    mColors.clear();
}

CSumatraSettings::~CSumatraSettings()
{
}

bool CSumatraSettings::ReadSettings()
{
    QString setdir = QCoreApplication::applicationDirPath();
    QString setname = setdir + "/SumatraSettings.json";

    QFile loadFile(setname);

    if (!loadFile.open(QIODevice::ReadOnly)) {
        qWarning("Couldn't open settings file.");
        QMessageBox::information(NULL, "Error", "Could not read settings file");
        return false;
    }

    QByteArray saveData = loadFile.readAll();
    QJsonDocument loadDoc(QJsonDocument::fromJson(saveData));

    if (!ParseJSON(loadDoc.object()))
        return false;

    return true;
}

bool CSumatraSettings::ParseJSON(const QJsonObject& json)
{
    if (json.contains("ColorBar") && json["ColorBar"].isDouble())
    {
        int ColorbarType = json["ColorBar"].toInt();
    }
    if (json.contains("PointSize") && json["PointSize"].isDouble())
    {
        int PointSize = json["PointSize"].toInt();
    }
    if (json.contains("BackgroundColor") && json["BackgroundColor"].isArray())
    {
        QJsonArray BackColor = json["BackgroundColor"].toArray();
        if (BackColor.size() != 3)
        {
            QMessageBox::information(NULL, "Warning", "Error reading settings file");
            return false;
        }
        int R = BackColor[0].toInt();
        int G = BackColor[1].toInt();
        int B = BackColor[2].toInt();
        mBackgroundColor.Set(R, G, B);
    }
    if (json.contains("colors") && json["colors"].isArray())
    {
        QJsonArray colors = json["colors"].toArray();
        // int sz = colors.size();
        for (int i = 0; i < colors.size(); i++)
        {
            QJsonArray onecolor = colors[i].toArray();
            int sz2 = onecolor.size();
            if (sz2 != 3)
            {
                QMessageBox::information(NULL, "Warning", "Error reading settings file");
                return false;
            }
                
            int R = onecolor[0].toInt();
            int G = onecolor[1].toInt();
            int B = onecolor[2].toInt();
            vtkColor3ub rgb(R, G, B);
            mColors.push_back(rgb);
        }
    }

    return true;
}
