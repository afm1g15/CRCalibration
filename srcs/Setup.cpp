#include "Setup.h"

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

} // calib (namespace)
