#include "Alarm.h"

volatile sig_atomic_t alarm_stopG = 0;
volatile sig_atomic_t alarmSet    = 0;

void on_alarm(int signal)
{
    alarm_stopG = 1;
    alarmSet = 0;
    std::cout << "TIME OUT!; alarm_stopG: "<<alarm_stopG<<"\n";
}

bool doStop(){return alarm_stopG==1;}

//inline
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
