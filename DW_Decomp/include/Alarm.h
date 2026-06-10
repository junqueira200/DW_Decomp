//
// Created by igor on 10/01/25.
//
#ifndef DW_ALARM_H
#define DW_ALARM_H

#include <csignal>
#include <iostream>
#include <unistd.h>

inline bool alarm_stopG = false;
inline bool alarmSet = false;
// inline unsigned int alarmPeriodG = 5*60; // 5 min

void on_alarm(int signal)
{
    alarm_stopG = true;
    alarmSet = false;
    std::cout << "TIME OUT!\n\n";
}

void setAlarm(unsigned int alarmPeriod)
{
    //static bool alarmSet = false;

    if(!alarmSet)
    {
        signal(SIGALRM, on_alarm);
        alarm(alarmPeriod);
        alarmSet = true;
        alarm_stopG = false;
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
    alarmSet = false;
    alarm_stopG = false;
}

#endif // DW_ALARM_H
