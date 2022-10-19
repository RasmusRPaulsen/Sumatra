function Component()
{
  // default constructor
}

Component.prototype.createOperations = function()
{
  component.createOperations();

  if (systemInfo.productType ==  "windows")
  {
    component.addOperation("CreateShortcut",
                   "@TargetDir@/sumatra.exe",
                   "@StartMenuDir@/sumatra.lnk",
                   "iconPath=@TargetDir@/sumatra.ico")
    component.addOperation("CreateShortcut",
                   "@TargetDir@/Sumatra.exe",
                   "@DesktopDir@/sumatra.lnk",
                   "iconPath=@TargetDir@/sumatra.ico")
    component.addOperation("RegisterFileType",
                   "vtp",
                   "@TargetDir@/sumatra.exe \"%1\"",
                   "VTK XML file extension",
                   "text/plain",
                   "@TargetDir@/Sumatra.ico",
				   "ProgId=Sumatra.vtp")
    component.addOperation("RegisterFileType",
                   "stl",
                   "@TargetDir@/sumatra.exe \"%1\"",
                   "STL file extension",
                   "text/plain",
                   "@TargetDir@/Sumatra.ico",
				   "ProgId=Sumatra.stl")
  }
}
