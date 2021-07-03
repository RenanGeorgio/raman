 #include <stdio.h> 
#include <time.h>     
#include "visa.h"   
#include "TLCCS.h"   

//==============================================================================
// Constantes
//===========================================================================   

#define MY_INTEGRATION_TIME   0.1            // 100 ms
#define MY_SAMPLE_FILE        "sample.txt"   // nome do arquivo para salvar os dados
#define MY_SCAN_COUNT         3              // numero de escaneamentos

//===========================================================================
// Globals
//===========================================================================   

ViSession   instr    = VI_NULL;                 // handle do instrumento
FILE*       my_file  = NULL;                    // file handling    

//===========================================================================
// Prototypes
//===========================================================================   

void error_exit(ViStatus err);
void waitKeypress(void);

//==============================================================================
// Main
//==============================================================================
int main (int argc, char *argv[])
{
   ViStatus    err      = VI_SUCCESS;           // error variable
   ViSession   resMgr   = VI_NULL;              // resource manager
   ViUInt32    i        = 0;                    // variavel para o loop
   ViUInt32    j        = 0;                    // variavel de loop secundaria
   ViUInt32    cnt      = 0;                    // contador para dispositivos
   ViUInt32    status   = 0;                    // variavel de status
   ViReal64    data[TLCCS_NUM_PIXELS];          // scan data array
   ViChar      rscStr[VI_FIND_BUFLEN];          // resource string
   ViChar*     rscPtr;                          // pointer para resource string
   time_t      t;                               // time structure

   FILE *filepoint;

   my_file = fopen_s(&filepoint, MY_SAMPLE_FILE, "w");
   //if(my_file == NULL)  return -1;

   printf("Abrindo Thorlabs CCS instrument driver\n");

   // Checando parametros e recursos
   if(argc < 2)
   {
      printf("Procurando por instrumentos CCS ...\n");
      if((err = viOpenDefaultRM(&resMgr))) error_exit(err);
      if((err = viFindRsrc(resMgr, TLCCS_FIND_PATTERN, VI_NULL, &cnt, rscStr))) error_exit(err);
      printf("Encontrado %u instrumento%s ...\n\n", cnt, (cnt>1) ? "s" : "");
      rscPtr = rscStr;
      viClose(resMgr);
   }
   else
   {
      // Recurso recebido na linha de comando
      rscPtr = argv[1];
   }  
   
   // tentando abrir CSS
   printf("Sessão aberta '%s' ...\n\n", rscStr);
   err = tlccs_init(rscStr, VI_OFF, VI_OFF, &instr);
   // error handling
   if(err)  error_exit(err);

   // tempo de integracao
   err = tlccs_setIntegrationTime(instr, MY_INTEGRATION_TIME);
   // error handling
   if(err)  error_exit(err); 

   // scaneamento inicial
   err = tlccs_startScan(instr);
   // error handling
   if(err)  error_exit(err); 

   while( i < MY_SCAN_COUNT )
   {
      // status do dispositivo
      err = tlccs_getDeviceStatus(instr, &status);
      // error handling
      if(err)  error_exit(err); 

      // camera encontrada -> disparar escanemaneto
      if(status & TLCCS_STATUS_SCAN_IDLE)
      {
         // dispara escanemaneto
         err = tlccs_startScan(instr);
         // error handling
         if(err)  error_exit(err); 
      }

      // dispositivo possui dados disponiveis para transferir
      if(status & TLCCS_STATUS_SCAN_TRANSFER)
      {
         printf("Iniciando escaneamento %d de %d ...\n\n", i+1, MY_SCAN_COUNT);
         
         // iniciando escanemnto
         err = tlccs_getScanData(instr, data);
         // error handling
         if(err)  error_exit(err); 

         // separador
         fprintf(my_file, "----------------- Scan Nu. %d -----------------\n", i+1);

         // recebendo tempo
         t = time(&t);

         // arquivando tempo
		 char tmBuff[30];
         fprintf(my_file, "%s\n", ctime_s(tmBuff, sizeof(tmBuff), &t));

         // salvando os dados no arquivo
         for(j = 0; j < TLCCS_NUM_PIXELS; j++)
         {
            fprintf(my_file, "Pixel: %4d - Valor: %f\n", j + 1, data[j]);
         }

         // um lopp de escanemento pronto
         i++;
      }
   }

   // numero total de processos de escanematos finalizados

   // fechando camera
   tlccs_close(instr);
   // close output file
   fclose(my_file);

   waitKeypress();
   
   // saindo da main
   return err;
}


/*---------------------------------------------------------------------------
  Error exit
---------------------------------------------------------------------------*/
void error_exit(ViStatus err)
{
   ViChar ebuf[TLCCS_ERR_DESCR_BUFFER_SIZE];

   // Printando erro
   tlccs_error_message (instr, err, ebuf);
   fprintf(stderr, "ERROR: %s\n", ebuf);
   
   // Fechando o handle do instrumento SE ele estiver aberto
   if(instr != VI_NULL) tlccs_close(instr);
   
   // Fechando stream SE estiver aberto
   if(my_file != NULL)  fclose(my_file); 
   
   // fechando programa
   waitKeypress();
   
   exit (err);
}


/*---------------------------------------------------------------------------
  Print keypress e espere
---------------------------------------------------------------------------*/
void waitKeypress(void)
{
   printf("Pressione <ENTER> e saia\n");
   while(getchar() == EOF);
}
