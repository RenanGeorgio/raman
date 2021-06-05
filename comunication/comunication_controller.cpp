#include "comunication_controller.hpp"



comunicationControll::comunicationControll() {
	err = VI_SUCCESS;       
	resMgr = VI_NULL;              
	i = 0;                
	j = 0;                    
	cnt = 0;                  
	status = 0;                  
}

comunicationControll::~comunicationControll() {

}

void comunicationControll::error_exit(ViStatus err) {
	ViChar ebuf[TLCCS_ERR_DESCR_BUFFER_SIZE];

	tlccs_error_message(instr, err, ebuf); //exibe mensagem de erro
	fprintf(stderr, "ERROR: %s\n", ebuf);

	if (instr != VI_NULL) tlccs_close(instr);

	if (my_file != NULL)  fclose(my_file);

	waitKeypress();

	exit(err);
}

void comunicationControll::waitKeypress() {
	printf("pressione <ENTER> para sair\n");
	while (getchar() == EOF);
}

void comunicationControll::process(char *value[], bool resoucers) {
	// try to open file
	my_file = fopen(MY_SAMPLE_FILE, "w");
	if (my_file == NULL)  return;

	printf("Abrindo aplicacao");

	// Parameter checking / Resource scanning
	if (resoucers)
	{
		// Find resources
		printf("Scanning for CCS instruments ...\n");
		if ((err = viOpenDefaultRM(&resMgr))) error_exit(err);
		if ((err = viFindRsrc(resMgr, TLCCS_FIND_PATTERN, VI_NULL, &cnt, rscStr))) error_exit(err);
		printf("Found %u instrument%s ...\n\n", cnt, (cnt > 1) ? "s" : "");
		rscPtr = rscStr;
		viClose(resMgr);
	}
	else
	{
		// Got resource in command line
		rscPtr = value[0];
	}

	// try to open CCS
	printf("Opening session to '%s' ...\n\n", rscStr);
	err = tlccs_init(rscStr, VI_OFF, VI_OFF, &instr);
	// error handling
	if (err)  error_exit(err);

	// set integration time
	err = tlccs_setIntegrationTime(instr, MY_INTEGRATION_TIME);
	// error handling
	if (err)  error_exit(err);

	// initial scan
	err = tlccs_startScan(instr);
	// error handling
	if (err)  error_exit(err);

	while (i < MY_SCAN_COUNT)
	{
		// request device status
		err = tlccs_getDeviceStatus(instr, status);
		// error handling
		if (err)  error_exit(err);

		// camera is idle -> we can trigger a scan
		if (getStatus & TLCCS_STATUS_SCAN_IDLE)
		{
			// trigger scan
			err = tlccs_startScan(instr);
			// error handling
			if (err)  error_exit(err);
		}

		// camera has data available for transfer
		if (getStatus & TLCCS_STATUS_SCAN_TRANSFER)
		{
			printf("Starting scan %d of %d ...\n\n", i + 1, MY_SCAN_COUNT);

			// trigger scan
			err = tlccs_getScanData(instr, data);
			// error handling
			if (err)  error_exit(err);

			// add seperator
			fprintf(my_file, "----------------- Scan No. %d -----------------\n", i + 1);

			// get time stamp
			t = time(&t);

			// store time stamp to file
			fprintf(my_file, "%s\n", ctime(&t));

			// store data to file
			for (j = 0; j < TLCCS_NUM_PIXELS; j++)
			{
				fprintf(my_file, "Pixel: %4d - Value: %f\n", j + 1, data[j]);
			}

			// one scan is done
			i++;
		}
	}

	// number of scans done

	// close camera
	tlccs_close(instr);
	// close output file
	fclose(my_file);

	waitKeypress();

	// leave main
	//cout << err;

	return;
}