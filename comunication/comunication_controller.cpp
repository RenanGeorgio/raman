#include "comunication_controller.hpp"



comunicationControll::comunicationControll() {
	
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
	printf("Press <ENTER> to exit\n");
	while (getchar() == EOF);
}
