#ifndef COMUNICATIONCONTROLL_H
#define COMUNICATIONCONTROLL_H

#pragma once

#include <ansi_c.h> 
#include <stdio.h>      
#include <time.h>       
#include <visa.h>       // we need visa (typically found in your VXIPNP\include directory)
#include "TLCCS.h"      // the device driver header

#define MY_INTEGRATION_TIME   0.1            // 100 ms
#define MY_SAMPLE_FILE        "output.txt"  
#define MY_SCAN_COUNT         3              

class comunicationControll
{

public:
	comunicationControll();
	~comunicationControll();

	void error_exit(ViStatus err);
	void waitKeypress();

	
protected:
	

private:
	ViSession   instr = VI_NULL;                 
	FILE*       my_file = NULL;                    
	
public:
	
};
#endif // COMUNICATIONCONTROLL_H