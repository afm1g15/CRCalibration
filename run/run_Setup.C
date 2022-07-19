{
  gROOT->ProcessLine(".L $RHIANA/srcs/Setup/Setup.cpp+");
  gROOT->ProcessLine(".L $RHIANA/srcs/Classes/Plane.cpp+");
  gROOT->ProcessLine(".L $RHIANA/srcs/Classes/Geometry.cpp+");
  gROOT->ProcessLine(".L $RHIANA/srcs/Helpers/PlottingHelpers.cpp+");
  gROOT->ProcessLine(".L $RHIANA/srcs/Helpers/GeometryHelpers.cpp+");
  gROOT->ProcessLine(".L $RHIANA/srcs/Helpers/SliceHelpers.cpp+");
  gROOT->ProcessLine(".L $RHIANA/srcs/Helpers/AnalysisHelpers.cpp+");
  gROOT->ProcessLine(".L $RHIANA/srcs/Helpers/ExternalFunctionHelpers.cpp+");
  gROOT->ProcessLine(".L $RHIANA/srcs/Setup/ReadFiles.cpp+");
  gROOT->ProcessLine(".L $RHIANA/srcs/Classes/ConfigReader.cpp+");
  gROOT->ProcessLine(".L $RHIANA/srcs/Classes/EventProcessor.cpp+");
}
