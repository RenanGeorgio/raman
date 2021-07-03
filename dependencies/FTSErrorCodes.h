////////////////////////////////////////////////////////////////////////
// FTSErrorCodes.h
//	This defines the error codes used by the  
//	Thorlabs Fourier Transform Spectrometer DLL. 
//	
//	Copyright (c) 2019 Thorlabs Sweden, All Rights Reserved
////////////////////////////////////////////////////////////////////////


#ifndef FTSERROR_CODES_H_
#define FTSERROR_CODES_H_

    /* These are possible error codes from the routines
    error codes 1->63 are return codes from the instrument */

/** FTS_SUCCESS indicates successful return from one of the FTS_... function, 
    the remaining FTS_ERROR_.. indicates some kind of error. */
#define FTS_SUCCESS                             0
#define FTS_ERROR_NOT_IMPLEMENTED               0x40
#define FTS_ERROR_NOT_WRITABLE                  0x41
#define FTS_ERROR_RANGE_ERROR                   0x42
#define FTS_ERROR_CMD_IGNORED                   0x43
#define FTS_ERROR_INVALID_CONFIGURATION         0x44
#define FTS_ERROR_COMMUNICATION_FAILURE         0x45
#define FTS_ERROR_TIMEOUT                       0x46
#define FTS_ERROR_SPECTROMETER_FAILURE          0x47
#define FTS_ERROR_SPECTROMETER_BUSY             0x48
#define FTS_ERROR_COMMUNICATION_NOT_OPENED      0x49
#define FTS_ERROR_SPECTROMETER_NOT_FOUND        0x4A
#define FTS_ERROR_IN_SETTINGS                   0x4B
#define FTS_ERROR_SPECTROMETER_UNPLUGGED        0x4C
#define FTS_ERROR_CARRIAGE_RESET                0X4D
#define FTS_ERROR_CARRIAGE_RESET_COMPLETE       0X4E
#define FTS_ERROR_BLOCKING_AIMING_BEAM_OPENED   0X50
#define FTS_ERROR_BLOCKING_AIMING_BEAM_CLOSED   0X51

#define FTS_ERROR_PARAMETER_ERROR               0x100
#define FTS_ERROR_DATA_ERROR                    0x101
#define FTS_ERROR_NOT_POWEROFTWO                0x102
#define FTS_ERROR_DIMENSION_MISMATCH            0x103
#define FTS_ERROR_OUTOFMEMORY                   0x104
#define FTS_ERROR_FILEIO                        0x105
#define FTS_ERROR_FILEERROR                     0x106
#define FTS_ERROR_SPECTRUM_NOT_FOUND            0x107
#define FTS_FAIL                                0x108
#define FTS_ERROR_ACQUISITION_ALREADY_RUNNING   0x109
#define FTS_ERROR_INSTRUMENT_STATUS             0x10A
#define FTS_ERROR_NODATA                        0x10B
#define FTS_ERROR_CONVERSION_ERROR              0x10C
#define FTS_ERROR_CONFIGURATION_ERROR           0x10D
#define FTS_ERROR_FIT_FAILED                    0x10E
#define FTS_ERROR_INVALID_DATA                  0x10F
#define FTS_ERROR_UNIT_MISMATCH                 0x110
#define FTS_ERROR_SPECTRUM_CALCULATION_FAILURE  0x111
#define FTS_ERROR_END_OF_FILE                   0x112
#define FTS_ERROR_SATURATED_INTERFEROGRAM       0x113
#define FTS_ERROR_INTERFEROGRAM_INTENSITY_LOW   0x114
#define FTS_ERROR_WAVELENGTHMETER_IRREGULAR_STRUCTURE    0x115
#define FTS_ERROR_WAVELENGTHMETER_ERROR_TOO_LARGE        0x116
#define FTS_CANCEL                                       0x117
#define FTS_PROGRESS                                     0x118
#define FTS_ANALYSIS_REGION_TOO_SMALL                    0x119
#define FTS_ERROR_INVALID_TIMESTAMP                      0x11A
#define FTS_ERROR_COHERENCE_LENGTH_TOO_LARGE             0x11B
#define FTS_ERROR_WAVELENGTHMETER_OUT_OF_RANGE			 0x11C

#endif