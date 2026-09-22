//
// Created by igor on 10/01/25.
//
#ifndef DW_ALARM_H
#define DW_ALARM_H

#include <csignal>
#include <iostream>
#include <unistd.h>

//inline  sig_atomic_t alarm_stopG;
//inline  sig_atomic_t alarmSet;
// inline unsigned int alarmPeriodG = 5*60; // 5 min

//inline
void on_alarm(int signal);
void setAlarm(unsigned int alarmPeriod);
void setOffAlarm();
bool doStop();

#endif // DW_ALARM_H
