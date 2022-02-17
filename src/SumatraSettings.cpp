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
    mCurrentColor = 0;
    mTest = "Constructed";
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
        mColorbarType = json["ColorBar"].toInt();
    }
    if (json.contains("PointSize") && json["PointSize"].isDouble())
    {
        mPointSize = json["PointSize"].toInt();
    }
    if (json.contains("BackgroundColor") && json["BackgroundColor"].isArray())
    {
        QJsonArray BackColor = json["BackgroundColor"].toArray();
        if (BackColor.size() != 3)
        {
            QMessageBox::information(NULL, "Warning", "Error reading settings file");
            return false;
        }
        mBackgroundColor[0] = BackColor[0].toDouble() / 255.0;
        mBackgroundColor[1] = BackColor[1].toDouble() / 255.0;
        mBackgroundColor[2] = BackColor[2].toDouble() / 255.0;
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

    mTest = mTest + QString(" Parsed");

    return true;
}

void CSumatraSettings::GetNextColor(double* col, bool range1)
{
    if (mColors.size() == 0)
    {
        col[0] = 255;
        col[1] = 255;
        col[2] = 255;
    }
    else
    {
        if (mCurrentColor >= mColors.size())
            mCurrentColor = 0;

        col[0] = mColors[mCurrentColor].GetRed();
        col[1] = mColors[mCurrentColor].GetGreen();
        col[2] = mColors[mCurrentColor].GetBlue();

        mCurrentColor++;
    }

    if (range1)
    {
        col[0] /= 255.0;
        col[1] /= 255.0;
        col[2] /= 255.0;
    }
}

