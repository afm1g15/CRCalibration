#include "PDSetup.h"

namespace calib{

  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void GetTime(time_t &rawtime){
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    std::cout << " Local time and date:  " << asctime(timeinfo)                << std::endl;
  } // GetTime
  
  // --------------------------------------------------------------------------------------------------------------------------------------------------
  
  void GetTotalTime(time_t &starttime, time_t &endtime){
    double timediff = difftime(endtime,starttime);

    // Calculate the time in a nice format
    int seconds_per_minute = 60;
    int seconds_per_hour   = 3600;
    int hours              = timediff / seconds_per_hour;
    int minutes            = (timediff - (hours * seconds_per_hour)) / (seconds_per_minute);
    int seconds            = timediff - (hours * seconds_per_hour) - (minutes * seconds_per_minute);

    // Now print
    std::cout << " Total time taken: "    << std::setw(4) << hours << " hours, ";
    std::cout                             << std::setw(4) << minutes << " minutes, ";
    std::cout                             << std::setw(4) << seconds << " seconds." << std::endl;
  } // GetTotalTime

  // --------------------------------------------------------------------------------------------------------------------------------------------------

  void FindReplace(std::string& s,
                   std::string const& toReplace,
                   std::string const& replaceWith,
                   bool verbose) {
    
    if(verbose){
      std::cout << __LINE__ << ". Finding all instances of '" << toReplace << "', within " << s << " and will replace with '" << replaceWith << "'" << std::endl;
    }

    // Find the string of interest within the main string
    std::size_t pos = s.find(toReplace);

    // Setup a boolean to allow for continuous searching,
    // to find all instances of a character
    bool foundNone = (pos == std::string::npos) ? true : false;
    
    // If it doesn't exist, give up
    while(!foundNone){

      if(verbose){
        std::cout << __LINE__ << ". Found an instance of '" << toReplace << "', at position " << pos << " in the string, will replace with '" << replaceWith << "'" << std::endl;
      }
      // Otherwise replace the string of interest with the new string
      s.replace(pos, toReplace.length(), replaceWith);
      
      // Now reassess the logic
      pos = s.find(toReplace.c_str(),pos+replaceWith.length(),s.length());
      if (pos == std::string::npos){
        foundNone = true;
        break;
      }
    } // While foundNone
  } // FindReplace

} // calib (namespace)
