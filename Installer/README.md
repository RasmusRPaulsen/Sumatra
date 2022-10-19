# Sumatra installer
Installation script for Sumatra

# Howto

- Compile newest version af Sumatra
- Update version info in `config/config.xml` and in `packages/com.dtucompute.sumatra/meta/package.xml`
- Copy `sumatra.exe`, `SumatraSettings.json` and other relevant files to the `packages/com.dtucompute.sumatra/data/` folder
- Start a command prompt in the `packages/com.dtucompute.sumatra/data/` folder and do:
```
C:\Qt\6.2.3\msvc2019_64\bin\windeployqt.exe sumatra.exe
``` 
- Copy all needed VTK DLL files from the `release` folder of the compiled VTK to the `packages/com.dtucompute.sumatra/data/` folder.
- Start a command prompt in the `install` folder and do:
```
C:\Qt\Tools\QtInstallerFramework\4.4\bin\binarycreator.exe -c config/config.xml -p packages SumatraInstaller.exe
```

### hints
(https://doc.qt.io/qtinstallerframework/index.html)
(https://doc.qt.io/qtinstallerframework/operations.html)
(https://www.youtube.com/watch?v=1pKMcwJZay4)
(https://medium.com/swlh/how-to-deploy-your-qt-cross-platform-applications-to-windows-operating-system-by-using-windeployqt-a7cd5663d46e)
(https://stackoverflow.com/questions/27521586/setting-program-files-as-a-default-installation-directory-in-the-qt-installer)
(https://stackoverflow.com/questions/58728532/qt-installer-framework-not-registering-file-type=
(https://github.com/apalomer/image_view/tree/master/installer)

