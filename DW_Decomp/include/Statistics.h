//
// Created by igor on 13/01/25.
//

#ifndef DW_STATISTICS_H
#define DW_STATISTICS_H

namespace StatisticsNS
{

    class StatisticsData
    {
    public:

        bool        timeLimit  = false;
        int         numIt      = 0;
        double      rootLB     = 0.0;
        double      rootTime   = 0.0;
        double      lowerBound = 0.0;
        double      upperBound = 0.0;
        double      gap        = 0.0;
        double      totalTime  = 0.0;
        std::string date       = "";
        std::string inst       = "";

        StatisticsData()
        {
            std::time_t result = duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
            date = std::string(std::asctime(std::localtime(&result)));
        };

    };

}

#endif //DW_STATISTICS_H
