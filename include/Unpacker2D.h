#ifndef __UNPACKER2D__
#define __UNPACKER2D__

#include <TObject.h>
#include <map>
#include <string>
#include <vector>

#include "Unpacker2.h"

class Unpacker2D : public Unpacker2 {
public:
  void UnpackSingleStep(const char* hldFile, const char* configFile, int numberOfEvents,
                        int refChannelOffset, const char* TOTcalibFile, const char* TDCcalibFile);

  void ParseConfigFile(std::string f, std::string s);
  void DistributeEventsSingleStep(std::string file);
  void BuildEvent(EventIII* e, std::map<UInt_t, std::vector<UInt_t> >* m, std::map<UInt_t, double>* refTimes);
  bool loadTDCcalibFile(const char* calibFile);

};

#endif
