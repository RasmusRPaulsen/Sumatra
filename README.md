# Sumatra
Surface manipulation and transformation toolkit (Sumatra) - a fast and versatile mesh viewer

This is very much a work in progress and the documentation needs severe upgrading.


## TODO
- Smooth dialog - add vtkConstrainedSmoothingFilter
- ICP
- Compare
- Scalars: smoothing, thresholding, connectivity

## Install Qt

- Download and install Qt in C:\Qt
- Add environment variable QTDIR=C:\Qt\6.2.3\msvc2019_64
- Add path %QTDIR%\bin
- Add path %QTDIR%\lib

## Install CMake

## Install VTK

- Download VTK 9.1.0
- Use CMake
- Enable Qt
- Compile VTK

- Add environment variable VTK_DIR=C:\libs\bin\VTK-9.1-MSVC2022
- Add path %VTK_DIR%\bin\Release

## Install GEL

- Get latest version af GEL from Github
- Compile using CMake
- Add environment variable GEL_SRC_DIR=C:\Users\rapa\Documents\src\libs\GEL

### Videos and inspirations



### Qt creating menus
(https://doc.qt.io/archives/qt-4.8/designer-creating-mainwindows.html)
(https://doc.qt.io/qt-5/qtwidgets-mainwindows-menus-example.html)
(https://www.youtube.com/watch?v=8z191FSk7n4)

### Qt drag and drop
(https://wiki.qt.io/Drag_and_Drop_of_files#Drag_and_Drop_of_files_on_QMainWindows)
(https://doc-snapshots.qt.io/qt6-dev/dnd.html)
(https://flylib.com/books/en/2.18.1/enabling_drag_and_drop.html)

### Qt dialogs
Passing values forth and back:
(https://stackoverflow.com/questions/20948665/passing-a-variable-to-other-dialog-qt)

### Parsing JSON
(https://doc.qt.io/qt-5/qtcore-serialization-savegame-example.html)

### Creating installer
(https://doc.qt.io/qtinstallerframework/index.html)
(https://www.youtube.com/watch?v=1pKMcwJZay4)
(https://medium.com/swlh/how-to-deploy-your-qt-cross-platform-applications-to-windows-operating-system-by-using-windeployqt-a7cd5663d46e)
(https://stackoverflow.com/questions/27521586/setting-program-files-as-a-default-installation-directory-in-the-qt-installer)
(https://stackoverflow.com/questions/58728532/qt-installer-framework-not-registering-file-type=
(https://github.com/apalomer/image_view/tree/master/installer)


### C++ naming conventions 
- from Google: (https://google.github.io/styleguide/cppguide.html#Naming)
- from Qt: (https://wiki.qt.io/Qt_Coding_Style)

## Linux build and installation

- Download and install Qt for Linux (https://doc.qt.io/qt-6/linux.html)
- Get the latest version of VTK (https://gitlab.kitware.com/vtk/vtk/-/blob/master/Documentation/docs/build_instructions/build.md?ref_type=heads)
- Download and install CMake for Linux
- Start CMake for VTK and create a binary output directory for vtk
- Enable build with Qt
- The QT6 include dir should be set to: /home/rapa/Qt/6.9.0/gcc_64/lib/cmake/Qt6
- Go to the vtk bin directory in a terminal and type `make`




### Building Qt application on Linux
(https://vitux.com/compiling-your-first-qt-program-in-ubuntu/)

