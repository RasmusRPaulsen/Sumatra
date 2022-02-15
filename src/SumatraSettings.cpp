#include "SumatraSettings.h"

#include <qdir.h>
#include <qcoreapplication.h>
#include <qmessagebox.h>
#include <qfile.h>
#include <QJsonDocument>
#include <qjsonarray.h>

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

    QString setname = setdir + "/SumatraSettings.json";

    QFile loadFile(setname);

    if (!loadFile.open(QIODevice::ReadOnly)) {
        qWarning("Couldn't open setings file.");
        QMessageBox::information(NULL, "Error", "Could not read settings file");
        return false;
    }

    QByteArray saveData = loadFile.readAll();
    QJsonDocument loadDoc(QJsonDocument::fromJson(saveData));

    ParseJSON(loadDoc.object());

    //QTextStream(stdout) << "Loaded save for "
    //    << loadDoc["player"]["name"].toString()
    //    << " using "
    //    << (saveFormat != Json ? "CBOR" : "JSON") << "...\n";
    return true;
}



bool CSumatraSettings::ParseJSON(const QJsonObject& json)
{
    // QMessageBox::information(NULL, "Info", "parsing");

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
        // QString BackgroundColor = json["BackgroundColor"].toString();
        QJsonArray BackgroundColor = json["BackgroundColor"].toArray();
        int R = BackgroundColor[0].toInt();
        int G = BackgroundColor[1].toInt();
        int B = BackgroundColor[2].toInt();
        double alpha = BackgroundColor[3].toDouble();

        // QMessageBox::information(NULL, "Test", BackgroundColor);
        // qWarning
    }
    if (json.contains("colors") && json["colors"].isArray())
    {
        //QStringList keystt = json["colors"].toArray().keys();


        // QString BackgroundColor = json["BackgroundColor"].toString();
        QJsonArray colors = json["colors"].toArray();
        int sz = colors.size();
        for (int i = 0; i < colors.size(); i++)
        {
            auto ttt = colors[i];
            QJsonObject obt = colors[i].toObject();
            QStringList keys = obt.keys();

            QJsonArray onecolor = colors[i].toArray();
            int sz2 = onecolor.size();

            QString tt = QString::number(sz) + QString(", ") + QString::number(sz2);
            if (sz == 0 and sz2 == 0)
                QMessageBox::information(NULL, "Test", tt);

            QString cname = onecolor[0].toString();
            int R = onecolor[1].toInt();
            int G = onecolor[2].toInt();
            int B = onecolor[3].toInt();
            double alpha = onecolor[4].toDouble();
            if (G == 128)
                QMessageBox::information(NULL, "Test", "test");

        }

        //int R = BackgroundColor[0].toInt();
        //int G = BackgroundColor[1].toInt();
        //int B = BackgroundColor[2].toInt();
        //double alpha = BackgroundColor[3].toDouble();

        // QMessageBox::information(NULL, "Test", BackgroundColor);
        // qWarning
    }

    
    //mPlayer.read(json["player"].toObject());

    //if (json.contains("player") && json["player"].isObject())
    //    mPlayer.read(json["player"].toObject());

    //if (json.contains("levels") && json["levels"].isArray()) {
    //    QJsonArray levelArray = json["levels"].toArray();
    //    mLevels.clear();
    //    mLevels.reserve(levelArray.size());
    //    for (int levelIndex = 0; levelIndex < levelArray.size(); ++levelIndex) {
    //        QJsonObject levelObject = levelArray[levelIndex].toObject();
    //        Level level;
    //        level.read(levelObject);
    //        mLevels.append(level);
    //    }
    //}
    return true;
}
