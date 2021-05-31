#ifndef COMUNICATIONCONTROLL_H
#define COMUNICATIONCONTROLL_H

#pragma once

#include <ansi_c.h> 
#include <stdio.h>      
#include <time.h>       
#include <visa.h>       
#include "TLCCS.h"      
#include <vector>


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

	void process(std::vector <std::string> value);

	
protected:
	

private:
	ViSession   instr = VI_NULL;                 
	FILE*       my_file = NULL;                    
	
public:
	ViStatus    err;           
	ViSession   resMgr;              
	ViUInt32    i;                    
	ViUInt32    j;                  
	ViUInt32    cnt;                    
	ViUInt32    status;                    
	ViReal64    data[TLCCS_NUM_PIXELS];          
	ViChar      rscStr[VI_FIND_BUFLEN];          
	ViChar*     rscPtr;                          
	time_t      t;                               
	
};
#endif // COMUNICATIONCONTROLL_H