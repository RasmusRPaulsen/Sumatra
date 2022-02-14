#include "SumatraSettings.h"

#include <qdir.h>
#include <qcoreapplication.h>
#include <qmessagebox.h>
#include <qfile.h>
#include <QJsonDocument>

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
        qWarning("Couldn't open save file.");
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
    QMessageBox::information(NULL, "Info", "parsing");

    if (json.contains("ColorBar"))
    {
        QMessageBox::information(NULL, "Test1", QString(json["ColorBar"].toString()));
        // qWarning

    }

    if (json.contains("BackgroundColor") && json["BackgroundColor"].isArray())
    {
        QMessageBox::information(NULL, "Test", QString(json["BackGroundColor"].toString()));
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
