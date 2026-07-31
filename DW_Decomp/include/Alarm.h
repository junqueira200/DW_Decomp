//
// Created by igor on 10/01/25.
//
#ifndef DW_ALARM_H
#define DW_ALARM_H

#include <csignal>
#include <iostream>
#include <unistd.h>

inline sig_atomic_t alarm_stopG = 0;
inline sig_atomic_t alarmSet    = 0;
// inline unsigned int alarmPeriodG = 5*60; // 5 min

void on_alarm(int signal)
{
    alarm_stopG = 1;
    alarmSet = 0;
    std::cout << "TIME OUT!; alarm_stopG: "<<alarm_stopG<<"\n";
}

void setAlarm(unsigned int alarmPeriod)
{
    //static bool alarmSet = false;
    static bool setSignal = true;
    if(alarmSet == 0)
    {
        //if(setSignal)
        {
            signal(SIGALRM, on_alarm);
            setSignal = false;
        }
        alarm(alarmPeriod);
        alarmSet = 1;
        alarm_stopG = 0;
    }
    else
    {
        std::printf("Error, alarmSet is Already set to true\n");
        throw "ERROR";
    }
}

void setOffAlarm()
{
    alarm(0);
    alarmSet = 0;
    alarm_stopG = 0;
}

#endif // DW_ALARM_H
