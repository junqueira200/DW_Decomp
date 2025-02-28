//
// Created by igor on 10/01/25.
//
#ifndef DW_ALARM_H
#define DW_ALARM_H

#include <unistd.h>
#include <csignal>

inline bool alarm_stopG = false;
//inline unsigned int alarmPeriodG = 5*60; // 5 min


void on_alarm(int signal)
{
    alarm_stopG = true;
    std::cout<<"TIME OUT!\n\n";
    exit(-1);
}



void setAlarm(unsigned int alarmPeriod)
{
    static bool alarmSet = false;

    if(!alarmSet)
    {
        signal(SIGALRM, on_alarm);
        alarm(alarmPeriod);
        alarmSet = true;
    }
}

#endif //DW_ALARM_H
