////////////////////////////////////////////////////////////////////////
// FTSData.h
//    This is the shared memory used by the Thorlabs Fourier Transform 
//    Spectrometer DLL. 
//    
//    Copyright (c) 2019 Thorlabs Sweden, All Rights Reserved
////////////////////////////////////////////////////////////////////////


#ifndef FTSDATA_H
#define FTSDATA_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif 

    // ------------------------------------------------------------------------
    // --------------------------- DEFINITIONS --------------------------------
    // ------------------------------------------------------------------------

#define FTS_LIB_VERSION                2.90

#define SERIALNUMBER_LENGTH             32
#define FIRMWARE_REV_LENGTH             30
#define SPECTRUM_COMMENT_LENGTH         128
#define SPECTRUM_NAME_LENGTH            32
#define SPECTRUM_SOURCE_STRLENGTH       24
#define SPECTRUM_OPERATOR_STRLENGTH     24
#define SPECTRUM_MODEL_STRLENGTH        24
#define MAX_GAINLEVELS_PER_SCAN         16
#define MAX_PATH_LENGTH                 512

#define MAX_SPECTROMETER_NUM            8     /** the maximum number of spectrometers that we can manage at any one time */
#define MAX_SPECTROMETER_CHANNELS       2     /** the maximum number of channels on an individual spectrometer that we can manage */
#define MAX_CALIBRATION_COEFFICIENTS    32    /** the maximum number of data points in a spectrometer calibration */
#define MAX_TRACKABLE_PEAKS             2048  /** the maximum number of peaks that we can track in one spectrum */
#define MAX_ZEROFILL_FACTOR             2     /** the maximum allowed value of 'zeroFillFactor' */
#define MAX_TRACEBUFFER_NUM             64
#define MAX_PHASE_POLY_ORDER            3
#define MAX_SENSITIVITY_MODE_NUM        8
#define MAX_RESOLUTION_MODE_NUM         4
#define MAX_DETECTOR_GAIN_LEVEL_NUM     64
#define MAX_SPECTRA_PER_FILE            256
#define MAX_WAVELENGTH_CALIBRATION_COEFFICIENTS    4

/** A SpectrometerIndex is used to enumerate the OSA Devices currently connected to the system. */
typedef unsigned short SpectrometerIndex;

/** A ChannelIndex is used to enumerate the channels on a spectrometer. */
typedef unsigned int ChannelIndex;

/** The FILE_FORMAT_.. flags define the different types of files that we can read or write
*    NOTICE: only the formats CSV and SPF2 can store all data in the spectrum_t
*    In SPC format are some parameters lost, in TXT format are only the arrays stored to file */
#define FILE_FORMAT_SPC        (0)        /** SPC spectrum file format */
#define FILE_FORMAT_CSV        (1)        /** Comma Separated Values */
#define FILE_FORMAT_SPF2       (2)        /** Thorlabs OSA Spectrum File */
#define FILE_FORMAT_TXT        (3)        /** Raw text file, no header */
#define FILE_FORMAT_MAT_V5     (4)        /** Matlab v5 binary data format */
#define FILE_FORMAT_CSV_ZIP    (5)        /** Zipped Comma Separated Values File */
#define FILE_FORMAT_TXT_ZIP    (6)        /** Zipped Raw text file, no header */
#define FILE_FORMAT_JCAMPDX    (7)        /** JCampDX ASCII data file */
#define FILE_FORMAT_HITRAN_XSC (8)        /** HITRAN XCS (Cross section) file */
#define FILE_FORMAT_NUM        (9)

#define ACQUISITION_MODE_DOUBLE_SIDED_F     1    /** only acquisition mode implemented so far */
#define ACQUISITION_MODE_SINGLE_SIDED_F     2    
#define ACQUISITION_MODE_DOUBLE_SIDED_B     3
#define ACQUISITION_MODE_SINGLE_SIDED_B     4
#define ACQUISITION_MODE_DOUBLE_SIDED_FB    5
#define ACQUISITION_MODE_SINGLE_SIDED_FB    6

#define APODISATION_NONE                0
#define APODISATION_NORTON_BEER_WEAK    1
#define APODISATION_NORTON_BEER_MEDIUM  2
#define APODISATION_NORTON_BEER_STRONG  3
#define APODISATION_TRIANGULAR          4
#define APODISATION_COSINE              5
#define APODISATION_HANN                6
#define APODISATION_HAMMING             7    /** a.k.a. Happ-Genzel or Blackman-Harris 2-term */
#define APODISATION_BLACKMANHARRIS3     8
#define APODISATION_BLACKMANHARRIS4     9
#define APODISATION_GAUSSIAN            10
#define APODISATION_TWO_PASS_HANN       11
#define APODISATION_TYPE_NUM            12
#define APODISATION_UNKNOWN             128

#define SMOOTHING_NONE                  0
#define SMOOTHING_BINOMIAL              1
#define SMOOTHING_SAVITSKY_GOLAY        2
#define SMOOTHING_FOURIER               3
#define SMOOTHING_MOVING_AVERAGE        4
#define SMOOTHING_GAUSS                 5
#define SMOOTHING_TYPE_NUM              6

#define PHASE_CORRECTION_NONE           0
#define PHASE_CORRECTION_MERTZ          1
#define PHASE_CORRECTION_FORMAN         2
#define PHASE_CORRECTION_MERTZ_FILE     3    /** not yet implemented in FTS_CalculateSpectrum */
#define PHASE_CORRECTION_FORMAN_FILE    4    /** not yet implemented in FTS_CalculateSpectrum */
#define PHASE_CORRECTION_TYPE_NUM       5
#define PHASE_CORRECTION_UNKNOWN        128

#define DERIVATIVE_FINITE_DIFFERENCE    0
#define DERIVATIVE_SAVITSKY_GOLAY       1
#define DERIVATIVE_TYPE_NUM             2

#define INTERFEROGRAM_PROPERTY_SATURATED         1 /** Flag set to indicate that the interferogram is saturated */
#define INTERFEROGRAM_PROPERTY_NONLINEAR_REGIME  2 /** Flag set to indicate that the intensity of the interferogram is so high that there is a risk of non-linear behavior */

#include "FTSErrorCodes.h"

#define FTS_OPTION_ARRAY_X      256
#define FTS_OPTION_MINMAX_X     512

#define FTS_AIR_MEASURE_TEMP        1
#define FTS_AIR_MEASURE_PRESS       2
#define FTS_AIR_MEASURE_RELHUM      4

#define FTS_BUFFER_INDEX_I      0
#define FTS_BUFFER_INDEX_X      1
#define FTS_BUFFER_INDEX_PHI    2

#define FTS_INSTR_STATUS_OK                         0x00
#define FTS_INSTR_STATUS_VSERROR                    0x01
#define FTS_INSTR_STATUS_REFERENCE_LOW              0x02
#define FTS_INSTR_STATUS_REFERENCE_HIGH             0x04
#define FTS_INSTR_STATUS_REFERENCE_WARMUP           0x08
#define FTS_INSTR_STATUS_REFERENCE_ERROR            0x10
#define FTS_INSTR_STATUS_BIAS_ERROR                 0x20
#define FTS_INSTR_STATUS_DC_SERVO_OVERLOAD          0x40
#define FTS_INSTR_STATUS_DETECTOR_TEMP_OVERLOAD     0x80
#define FTS_INSTR_STATUS_MOTOR_TIMEOUT              0x100
#define FTS_INSTR_STATUS_USB_TRANSFER_ERROR         0x200
#define FTS_INSTR_STATUS_MOTOR_ENDPOINT_ERROR       0x400
#define FTS_INSTR_STATUS_AIMING_BEAM_ON             0x800
#define FTS_INSTR_STATUS_DATA_ERROR                 0x10000

#define FTS_AIMING_BEAM_ON                          (0x01)  /** Flag to indicate that the aiming beam of this device is turned on */
#define FTS_HAS_BLOCKING_AIMING_BEAM                (0x02)  /** Flag to indicate that this device has a blocking, switchable aiming beam */
#define FTS_HAS_NON_BLOCKING_AIMING_BEAM            (0x04)  /** Flag to indicate that this device has a non-blocking aiming beam */


    /** Indicates that the wavelength meter (or coherence length) analysis has failed
    to calculate a value */
#define FTS_WAVELENGTHMETER_ILLEGAL_VALUE        (-1.0)

    // ------------------------------------------------------------------------
    // ----------------------- DATA STRUCTURES --------------------------------
    // ------------------------------------------------------------------------

    typedef unsigned short INSTRUMENT_MODEL;
#define INSTRUMENT_UNKNOWN               0x0000
#define INSTRUMENT_OSA201                0x0001
#define INSTRUMENT_OSA201B               0x0011
#define INSTRUMENT_OSA201C               0x0021
#define INSTRUMENT_OSA202                0x0002
#define INSTRUMENT_OSA202C               0x0022
#define INSTRUMENT_OSA203                0x0003
#define INSTRUMENT_OSA203B               0x0013
#define INSTRUMENT_OSA203C               0x0023
#define INSTRUMENT_OSA205                0x0005
#define INSTRUMENT_OSA205C               0x0025
#define INSTRUMENT_OSA206                0x0006
#define INSTRUMENT_OSA206C               0x0026
#define INSTRUMENT_OSA207                0x0007
#define INSTRUMENT_OSA207C               0x0027
#define INSTRUMENT_VIRTUAL_OSA201        0x0100
#define INSTRUMENT_VIRTUAL_OSA202        0x0200
#define INSTRUMENT_VIRTUAL_OSA203        0x0300
#define INSTRUMENT_VIRTUAL_OSA203B       0x1300
#define INSTRUMENT_VIRTUAL_OSA205        0x0500

    typedef unsigned short TRIGGER_MODE;
#define TRIGGER_NONE            (0)
#define TRIGGER_SOFTWARE        (1)
#define TRIGGER_RISING_FLANK    (2)
#define TRIGGER_FALLING_FLANK   (3)

    typedef unsigned short SPECTRUM_TYPE;
#define SPEC_INTERFEROGRAM    0            /** indicates that this spectrum is an interferogram */
#define SPEC_EMISSION         1            /** emission spectrum, this default for Fourier transformed interferograms */
#define SPEC_TRANSMITTANCE    2            /** indicates that this spectrum is in units of transmittance */
#define SPEC_ABSORBANCE       3            /** indicates that this spectrum is in units of absorbance */
#define SPEC_KUBELKA_MUNK     4            /** indicates that this spectrum is a Kubelka Munk spectrum */
#define SPEC_REFLECTANCE      5            /** indicates that this spectrum is in units of reflectance */
#define SPEC_LOG_1_R          6            /** indicates that this spectrum is log(1 / reflectance) */
#define SPEC_INTERFEROGRAM_ENVELOPE    7   /** indicates that this spectrum is an interferogram envelope */
#define SPEC_PHASE            8            /** indicates that this is a phase spectrum */
#define SPEC_CROSSSECTION     9            /** indicates that this is an absorption cross section (y-axis unit should be Y_UNIT_CM2_MOLEC ). */
#define SPEC_TIME_SERIES      10           /** indicates that this is a time series measurement. The x-axis unit should be seconds. */
#define SPEC_OSCILLOSCOPE     11           /** indicates that this is a oscilloscope mode time series measurement. The x-axis unit should be seconds. */
#define SPEC_UNKNOWN          99           /** unknown type */

    typedef unsigned short INSTRUMENT_TYPE;
#define INSTRUMENT_TYPE_FOURIER_TRANSFORM   (0)
#define INSTRUMENT_TYPE_FABRY_PEROT         (1)
#define INSTRUMENT_TYPE_CCD_GRATING         (2)

    // possible options for spectrum_t->xValueFormat
#define X_VAL_MINMAX    0    /** only min and max values stored, rest are linearly interpolated */
#define X_VAL_ARRAY     1    /** x-values are stored in array 'x' */

    // possible options for spectrum_t->phiValueFormat
#define PHI_VAL_NONE        0    /** No phase data is included in the spectrum */
#define PHI_VAL_ARRAY       1    /** phase values are stored in array 'phi' */
#define PHI_VAL_IMAGINARY   2    /** the 'phi' array is allocated and contains the imaginary component */

    // possible units on the X-axis of a spectrum/an interferogram
    typedef unsigned short X_AXIS_UNIT;
#define UNIT_NO_UNIT    0        /** no unit */
#define X_UNIT_CM       1        /** distance in cm (for interferograms) */
#define X_UNIT_CM_1     2        /** wave numbers in inverse cm (cm^-1) */
#define X_UNIT_THZ      3        /** frequency in Tera Hertz */
#define X_UNIT_NM_VAC   4        /** wavelength in nanometers (in vacuum) */
#define X_UNIT_NM_AIR   5        /** wavelength in nanometers (in air) */
#define X_UNIT_EV       6        /** photon energy in Electron Volts */
#define X_UNIT_INDEX    7        /** x-value is an index */
#define X_UNIT_SECONDS  8        /** time in seconds */
#define X_UNIT_PIXEL    9        /** pixel number on a detector */
#define X_UNIT_UNKNOWN  65534    /** unknown unit */

    // possible units on the Y-axis of a spectrum/an interferogram
    typedef unsigned short Y_AXIS_UNIT;
#define Y_UNIT_AU           32769    /** Arbitrary Units, or Absorbance Units */
#define Y_UNIT_COUNTS       32770    /** AD counts */
#define Y_UNIT_DBM          32771    /** dBm, absolute power */
#define Y_UNIT_DBM_NORM     32772    /** dBm, power density (dBm/x-axis unit, e.g. dBm/nm) */
#define Y_UNIT_MW           32773    /** milli(!!) Watts, absolute power */
#define Y_UNIT_MW_NORM      32774    /** milli(!!) Watts, power density (dBm/x-axis unit, e.g. dBm/nm) */
#define Y_UNIT_PERCENT      32775    /** percent */
#define Y_UNIT_RADIANS      32776    /** angle in radians */
#define Y_UNIT_DEGREES      32777    /** angle in degrees */
#define Y_UNIT_CELSIUS      32778    /** degrees Celsius */
#define Y_UNIT_KELVIN       32779    /** temperature in Kelvin */
#define Y_UNIT_HPA          32780    /** pressure in hectoPascals */
#define Y_UNIT_LOG          32781    /** logarithmic unit (general log unit) */
#define Y_UNIT_DB           32782    /** dB, relative power */
#define Y_UNIT_PP           32783    /** percentage points */
#define Y_UNIT_MIXED        32784    /** mixed units */
#define Y_UNIT_INTENSITY    32785    /** non-calibrated intensity. */
#define Y_UNIT_LOGINTENSITY 32786    /** log non-calibrated intensity (log(Y_UNIT_INTENSITY)) */
#define Y_UNIT_CAL_I        32787    /** calibrated intensity. */
#define Y_UNIT_CAL_LOGI     32788    /** log calibrated intensity (log(Y_UNIT_INTENSITY)) */
#define Y_UNIT_ATMOSPHERE   32789    /** pressure in atmospheres */
#define Y_UNIT_TORR         32790    /** pressure in Torr */
#define Y_UNIT_PSI          32791    /** pressure in pounds per square inch */
#define Y_UNIT_LOG_COUNTS   32792    /** logarithm of counts */
#define Y_UNIT_CM2_MOLEC    32793    /** cross section of cm^2 per molecule */
#define Y_UNIT_MOLEC_CM2    32794    /** column of molecules / cm^2 (number concentration times distance) */
#define Y_UNIT_GRAMS_CM2    32795    /** column of grams / cm^2 (mass concentration times distance) */
#define Y_UNIT_MOLEC_CM3    32796    /** number concentration, molecules / cm^3 */
#define Y_UNIT_GRAMS_CM3    32797    /** mass concentration, grams / cm^3 */
#define Y_UNIT_MOL_DM3      32798    /** number concentration, mol / dm^3 = mol / liter = M */
#define Y_UNIT_TRANSMITTANCE 32799    /** Transmittance unit */
#define Y_UNIT_UNKNOWN      65535    /** unknown unit */

    typedef struct instrument_statistics_t
    {
        /** The total number of time outs for this instrument, 
        *   since the start of the program */
        unsigned int totalTimeoutNum;

        /** The number of consecutive time outs for this instrument.
        *   Reset to zero whenever a successful acquisition is done. */
        unsigned int consecutiveTimeoutNum;

        /** The total number of interferograms that have been
        *   acquired by this specific instrument since the start
        *   of the program. */
        unsigned int totalInterferogramsAcquired;

        /** The total time (in seconds) that the program has spent
        *   waiting for or retrieving interferograms from the instrument
        *   since the start of the program */
        double totalInterferogramAcquisitionTime_s;

        /** The total amount of data acquired from the instrument since
        *   the start of the program, in Mega Bytes (1 000 000 bytes) */
        double totalDataAcquired_MBytes;

        /** The total number of end point resets */
        unsigned int endPointResetCount;

        /** The total number of carriage resets */
        unsigned int carriageResetCount;
    } instrument_statistics_t;
     
    /** An fts_detector_t structure describes the properties of one detector. */
    typedef struct fts_detector_t
    {
        /** The smallest wavenumber (corresponding to the largest wavelength) that this detector can measure */
        float minWaveNr;

        /** The largest wavenumber (corresponding to the smallest wavelength) that this detector can measure */ 
        float maxWaveNr;

        /** Array of wavenumbers at which the calibration was made. These are always in ascending order. */
        float calib_WaveNr[MAX_CALIBRATION_COEFFICIENTS];

        /** Array of calibration coefficients at the wavenumbers specified in calib_WaveNr */
        float calib_power[MAX_CALIBRATION_COEFFICIENTS];

        /** The date that the calibration was made (=(the year)*10000 + (the month)*100 + (the day)*1) */
        unsigned int calib_date;

        /** The number of calibration coefficients, i.e. the length of calib_WaveNr and calib_power */
        unsigned int calib_coeffNum;
        
        /** Describes when the system will issue a non-linearity warning for a too high
        *   input intensity */
        unsigned int non_linear_warning_per_mille;
        
        /** Array containing the available gain levels of the detector.
        *   This array contains the actual gains of the detector. */
        double det_GainLevels[MAX_DETECTOR_GAIN_LEVEL_NUM];

        /** The number of gain levels of the detector, 
        *   i.e. the length of the array det_GainLevels */
        unsigned short det_GainLevelNum;
        
        /** true if the temperature of the detector can be controlled in software */
        bool det_TemperatureControlImplemented;

        /** The current temperature of the detector(s), in degrees C. 
        *   This only makes sense if #det_TemperatureControlImplemented is true. */
        float det_currentTemperature_C;

        /** The number of temperature modes of the detector(s) of the device.
        *   If #det_TemperatureControlImplemented is false then this is zero */
        unsigned short det_temperatureModeNum;

        /** Number of bits on the ADC */
        unsigned short adcBits;

    } fts_detector_t;

    /** An fts_instrument_t structure describes the properties of one instrument.
    *   The FTSData structure holds an array of fts_instrument_t elements, each represents a 
    *   connected instrument and is set up when calling #FTS_InitializeSpectrometers(). */
    typedef struct fts_instrument_t
    {
        /** The serial number of this instrument. Null terminated string */
        char serial[SERIALNUMBER_LENGTH];

        /** The firmware version of this instrument. Null terminated string */
        char firmwareVersion[FIRMWARE_REV_LENGTH];

        /** The model of this instrument. */
        INSTRUMENT_MODEL model;

        /** The properties of the detector(s) */
        fts_detector_t detector[MAX_SPECTROMETER_CHANNELS];

        /** The number of detectors, i.e. number of valid entries in the array 'detector' */
        unsigned short detectorNum;

        /** Length of each package in the data transfer. Internal use only. */
        unsigned short sectorSize;

        /** these are the speeds of the motor at different speed modes (in cm/s) */
        float motor_speed_cm_s[MAX_SENSITIVITY_MODE_NUM];

        /** this is the number of different sensitivity modes of the instrument */
        unsigned short motor_sensitivityModeNum;

        /** Array containing the travel ranges of the motor at different resolution modes (in cm).
        *   This is in units of Optical Path Difference. */
        float motor_travel_cm[MAX_RESOLUTION_MODE_NUM];

        /** acquisitionStart_MinOPD is the Optical Path Difference when the acquisition
        *   begins, i.e. at the first data point in the retrieved interferograms. 
        *   This is usually negative - indicating that the sampling began before
        *   the Zero Path Difference). The OPD at the last point in the interferogram
        *   equals acquisitionStart_MinOPD_cm[modeIndex] + motor_travel_cm[modeIndex]. */
        float acquisitionStart_MinOPD_cm[MAX_RESOLUTION_MODE_NUM];

        /** this is the number of different resolution modes that the instrument has
        *   (and also the length of the array motor_travel_cm) */
        unsigned short motor_resolutionModeNum;

        /** Wavelength of the reference laser, in nm  */
        double referenceWavelength_nm_vac;

        /** Number of samples/reference wavelength in the currently set 
        *   sensitivity mode (NB varies with sensitivity mode) */
        float samplesPerReferenceWavel;

        /** Current temperature inside the instrument, in degrees Celsius. */
        float atm_curTemp_C;

        /** Current atmospheric pressure inside the instrument, in hecto Pascals. */
        float atm_curAtmPress_hPa;

        /** Current relative humidity inside the instrument, in percent. */
        float atm_curRelHum;

        /** Flag describing if the instrument has an temperature, pressure and humidity sensor
        *   inside. This is an OR:ed combination of #FTS_AIR_MEASURE_TEMP, #FTS_AIR_MEASURE_PRESS and #FTS_AIR_MEASURE_RELHUM  */
        unsigned int atmosphericSensor;

        /** The current internal status of the instrument 
        *   (an or:ed combination of the FTS_INSTR_STATUS_.. flags in FTSData.h).
        *   If this equals #FTS_INSTR_STATUS_VSERROR then an internal error has
        *   occurred in the instrument and no acquisitions can be processed.
        *   A more detailed error flag can be found in instrumentDetailedStatus */
        unsigned int instrumentStatus;

        /** The internal error flag of the instrument.
        *   This makes sense only if instrumentStatus equals #FTS_INSTR_STATUS_VSERROR */
        unsigned int instrumentDetailedStatus;

        /** The current status of the communication with this instrument
        *   (can be; \c #FTS_SUCCESS, \c #FTS_ERROR_COMMUNICATION_FAILURE, \c #FTS_ERROR_SPECTROMETER_BUSY, 
        *   \c #FTS_ERROR_COMMUNICATION_NOT_OPENED, \c #FTS_ERROR_SPECTROMETER_NOT_FOUND or \c #FTS_ERROR_INSTRUMENT_STATUS)
        *   This is set by a call to #FTS_CheckSpectrometer */
        unsigned int communicationStatus;

        /** This holds information on the statistics of the communication with
        *   the instrument. Such as speed, quality of the link and how many
        *   interferograms have been acquired */
        instrument_statistics_t statistics;

        /** Array of wavelength correction parameters. Internal use only. */
        double wavelength_correction[MAX_WAVELENGTH_CALIBRATION_COEFFICIENTS];

        /** Stores the time stamp of the last retrieval of an interferogram from this device */
        char timeOfLastRetrieval[24];

        /** The aiming beam properties of this device. This is an or:ed combination of the aiming beam flags;
        *   #FTS_AIMING_BEAM_ON, #FTS_HAS_BLOCKING_AIMING_BEAM, #FTS_HAS_NON_BLOCKING_AIMING_BEAM */
        unsigned int aiming_beam_status;

    } fts_instrument_t;

    /** An acquisition_option_t structure describes the settings to use when retrieving data from an
    *   instrument or when converting an interferogram to a spectrum.
    *   The FTSData structure holds an array of acquisition_option_t elements, each
    *   represents a connected instrument and how to process data from that specific instrument */
    typedef struct acquisition_option_t
    {
        /** level of smoothing to apply to the measured interferogram 
        *   (only first two elements are used, rest are for future use) */
        unsigned short    interferogramSmoothParam[4];

        /** type of smoothing to apply to the measured interferogram. Must be one of the SMOOTHING_.. constants 
        *   defined in FTSData.h */
        unsigned short    interferogramSmoothType;

        /** level of smoothing to apply to the calculated spectrum 
        (only first two elements are used, rest are for future use) */
        unsigned short    spectrumSmoothParam[4];

        /** type of smoothing to apply to the calculated spectrum. Must be one of the SMOOTHING_.. constants 
        *   defined in FTSData.h */
        unsigned short    spectrumSmoothType;

        /** apodisation to make on the interferograms. Must be one of the APODISATION_.. constants 
        *   defined in FTSData.h  */
        unsigned char    apodisationType;

        /** the number of spectra to average */
        unsigned int    averageSpecNum;

        /** if this is true then the spectrum will be updated 
        as a rolling average with period of 'averageSpecNum' */
        bool            rollingAverage;

        /** How many more data points will be added to the spectrum before the FFT, 
        *   number of data points = 2^zeroFillFactor (i.e. zeroFillFactor=1 <=> 2x longer spectrum).
        *   Allowed values are 0 to #MAX_ZEROFILL_FACTOR */
        unsigned char zeroFillFactor;

        /** The type of phase correction to apply. Must be either of #PHASE_CORRECTION_MERTZ, #PHASE_CORRECTION_FORMAN
        *   #PHASE_CORRECTION_MERTZ_FILE or #PHASE_CORRECTION_FORMAN_FILE */
        unsigned char phaseCorrType;

        /** Full file name and path to a file containing the phase information to use in the phase correction.
        *   This is only used if phaseCorrtype equals #PHASE_CORRECTION_MERTZ_FILE or #PHASE_CORRECTION_FORMAN_FILE.
        *   THIS IS CURRENTLY NOT IMPLEMENTED. RESERVED FOR FUTURE USE  */
        char phaseCorrFile[MAX_PATH_LENGTH];

        /** Full file name and path to a file containing the the noise floor for the device in the current mode.
        *   This can be used by the phase correction, if phaseCorrType equals #PHASE_CORRECTION_MERTZ or #PHASE_CORRECTION_FORMAN
        *   to determine the data points which have an intensity above the noise floor. */
        char noise_floor_file[MAX_PATH_LENGTH];

        /** True if a phase spectrum should be calculated while we're calculating the spectrum.  
        *   The calculated phase spectrum will be stored in the spectrum_t::phi array */
        bool calculatePhaseSpectrum;

        /** The collected interferogram will be cropped to the
        *   given distance range after retrieval. Unit in cm
        *   For double sided interferograms should minOPD be negative...
        *   If both are 0 or minOPD_cm >= maxOPD_cm then no cropping will be
        *   performed  */
        float            minOPD_cm;
        float            maxOPD_cm;

        /** The range of wavenumbers to limit the spectrum to. 
        *   Must be within the range of wavenumbers defined for the spectrometer.  */
        float            minWaveNr;
        float            maxWaveNr;

        /** this is the current resolution mode of the instrument
        *   (index into #fts_instrument_t::motor_travel_cm )
        *   Must be less than #fts_instrument_t::motor_resolutionModeNum
        *   This is set by calling #FTS_SetAcquisitionOption_ResolutionMode() */
        unsigned short    motor_curResolutionMode;

        /** this is the current sensitivity mode of the instrument
        *   (index into #fts_instrument_t::motor_speed_cm_s ) 
        *   Must be less than #fts_instrument_t::motor_sensitivityModeNum
        *   This is set by calling #FTS_SetAcquisitionOption_SensitivityMode() */
        unsigned short    motor_curSensitivityMode;

        /** If autoGain is true then the gain levels are automatically determined.
        *   The calculations are done after the collection of an interferogram */
        bool            autoGain;

        /** if singleGain is true then the auto gain will only 
        *   use one gain level throughout the interferogram
        *   ONLY IMPLEMENTED OPTION IS 'TRUE' */
        bool            singleGain;

        /** Number of gain levels that will be used in the scan
        *   ONLY IMPLEMENTED OPTION IS ONE */
        unsigned char    gainLevelNum;

        /** The position in the scan where we will change the gain level.
        *   gainPosition[0] MUST BE ZERO */
        unsigned int     gainPosition[MAX_GAINLEVELS_PER_SCAN];

        /** The level of the gain at each position (index into #fts_instrument_t::det_GainLevels) 
        *   This must be smaller than det_GainLevelNum. Notice that not all gain levels in 
        *   #fts_instrument_t::det_GainLevels are allowed in all sensitivity modes. Set the gain
        *   by calling #FTS_SetAcquisitionOption_SingleGain() to be sure a legal value is set. 
        *   If acquisition_options_t::autoGain is true then this will be updated after each
        *   collected interferogram. */
        unsigned short    gainIndex[MAX_GAINLEVELS_PER_SCAN];

        /** interferogramIntensity_Acceptable defines the acceptable range of intensities for the measured interferogram,
        *   used by the automatic gain routine.
        *   These values must be in the range 0 -> 1.0 (in fractions of the full range of the ADC).
        *   If the interferogram has an intensity within this range then the auto gain will not try to change this */
        float interferogramIntensity_AcceptableLow;
        float interferogramIntensity_AcceptableHigh;

        /** interferogramIntensity_Ideal defines the preferred range of intensities for the measured interferogram,
        *   used by the automatic gain routine.
        *   These values must be in the range 0 -> 1.0 (in fractions of the full range of the ADC).
        *   The range interferogramIntensity_IdealLow to interferogramIntensity_IdealHigh of the full range is what the auto gain routine
        *   will aim for when calculating the next gain level to use. */
        float interferogramIntensity_IdealLow;
        float interferogramIntensity_IdealHigh;

        /** The desired temperature mode of the detector(s). If det_TemperatureControlImplemented
        *   is false in the instrument then this parameter is ignored. */
        unsigned short    det_temperatureMode;

        /** true if we should perform cycle-counting (wavelength metering mode) on the interferograms */
        bool            cycleCounting;

        /** true if we should perform coherence length analysis on the interferograms */
        bool            coherenceLengthAnalysis;

        /** the options for if we should read the data from the 
        *   instrument or use values from the user 
        *   (any combination of #FTS_AIR_MEASURE_TEMP, #FTS_AIR_MEASURE_PRESS and #FTS_AIR_MEASURE_RELHUM) */
        unsigned int    air_measureOption;

        /** Minimum interval between two acquisitions, in milliseconds */
        unsigned int    acqInterval;

        /** This is a description of the parameters used. To be used as a reference
        *   when e.g. storing or loading acquisition_option_t to/from file (see 
        *   #FTS_SaveAcquisitionOptions() or #FTS_LoadAcquisitionOptions() )*/
        char    description[128];

        /** The acquisition mode, should for now equal #ACQUISITION_MODE_DOUBLE_SIDED_F */
        unsigned char    acquisitionMode;

        unsigned short    setResolutionMode;

        unsigned short    setSensitivityMode;

        unsigned short    setDetectorTemperatureMode;
    }acquisition_option_t;

#define spectrum_t_hdr_version 9

    /** An spectrum_t structure contains all data for a collected interferogram or spectrum.
    The data values are stored in the array spectrum_t::I.
    For interferograms are the array spectrum_t::x usually NULL and the x-axis values are described
    by the spectrum_t::x_min and spectrum_t::x_max values.
    */
    typedef struct spectrum_t
    {
        /** This is the size in bytes of the header (note that this is different on 64 and 32 bit machines) */
        unsigned short    hdrsize;

        /** The version of the header */
        unsigned short    hdrversion;

        /** the number of data points in the spectrum_t::I array */
        unsigned int    length;

        /** the unit of the x-axis */
        X_AXIS_UNIT        xAxisUnit;

        /** the unit of the y-axis */
        Y_AXIS_UNIT        yAxisUnit;

        /** the acquisition mode (should be equal to #ACQUISITION_MODE_DOUBLE_SIDED_F ) */
        unsigned char    acquisitionMode;

        /** specifies how the x-axis values are stored, either #X_VAL_MINMAX or #X_VAL_ARRAY, see \ref sec_spectral_x_axis_data. */
        unsigned char    xValueFormat;

        /** the smallest x-value, in wavenumbers. 
        In an interferogram does this define the range for the interferometer */
        float            x_minWnr;

        /** the largest x-value, in wavenumbers. 
        In an interferogram does this define the range for the interferometer */
        float            x_maxWnr;

        /** the smallest x-value */
        float            x_min;

        /** the largest x-value */
        float            x_max;

        /** the smallest value in the spectrum/interferogram array spectrum_t::I */
        float            y_min;

        /** the largest value in the spectrum/interferogram array spectrum_t::I */
        float            y_max;

        /** the gain levels used when collecting this interferogram.
        All non-used values are set to zero. */
        double            gainLevel[MAX_GAINLEVELS_PER_SCAN];

        /** the positions where the gain was set when collecting this interferogram 
        (number is defined by number of non-zero items here) */
        unsigned int    gainPos[MAX_GAINLEVELS_PER_SCAN];

        /** the resolution of the spectrum, always in cm^-1 (ONLY USED FOR SPECTRA) */
        float            resolution;

        /** the type of this spectrum. See the definition of #SPECTRUM_TYPE */
        SPECTRUM_TYPE    type;

        /** the serial number of the  interferometer that collected this spectrum */
        char            interferometerSerial[SERIALNUMBER_LENGTH];

        /** This is the wavelength of the reference laser in the instrument. In nanometers in vacuum */
        double            referenceWavelength_nm_vac;

        /** the spatial distance between two samples in the interferogram, in cm in vacuum */
        double            samplingDistance_cm_vac;

        /** the date when the spectrum was acquired 
        (=(the year)*10000 + (the month)*100 + (the day)*1) */
        unsigned int    date; 

        /** the (local) time of day when the spectrum was acquired 
        (=(the hour)*1000000 + (minutes)*10000 (seconds)*100 + (milli seconds)/10) */
        unsigned int    time; 

        /** the (UTC) time of day when the spectrum was acquired 
        (=(the hour)*1000000 + (minutes)*10000 (seconds)*100 + (milli seconds)/10) */
        unsigned int    gmtTime; 

        /** how many readouts this spectrum is an average of */
        unsigned int    averageNum;

        /** the parameters used for the smoothing applied to this spectrum */
        unsigned short    smoothParam[4]; 

        /** the type of smoothing applied to this spectrum.
        Equals one of the SMOOTHING_... constants in FTSData.h */
        unsigned short    smoothType;

        /** the parameters used for the smoothing applied to the interferogram (before calculating the spectrum) */
        unsigned short    igramSmoothParam[4];

        /** the type of smoothing applied to the interferogram (before calculating the spectrum).
        Equals one of the SMOOTHING_... constants in FTSData.h */
        unsigned short    igramSmoothType;

        /** the type of phase correction applied to this spectrum 
        (must be one of the PHASE_CORRECTION_... constants in FTSData.h) */
        unsigned char    phaseCorrection;

        /** the apodisation performed on this spectrum 
        (must be one of the APODISATION_... constants in FTSData.h) */
        unsigned char    apodisation;

        /** (only for spectra) how many times longer we have made the interferogram 
        before the FFT (1 means nothing, 2 means doubling the length etc) */
        float            zeroFillFactor;

        /** the temperature at the time of the measurement, in degrees C */
        float            air_temp;

        /** the air pressure at the time of the measurement, in hekto Pascal (hPa) */
        float            air_press;

        /** the relative humidity at the time of the measurement, in percent */
        float            air_relHum;

        /** the name of the spectrum     */
        char            name[SPECTRUM_NAME_LENGTH];

        /** null-terminated string with comments... */
        char            comment[SPECTRUM_COMMENT_LENGTH];

        /** this is the number of allocated floats for the array 'I' 
        (should be set! this is used by e.g. #FTS_CopySpectrum() to make sure the arrays are good) */
        unsigned int    __allocatedLengthI;

        /** this is the allocated allocated floats for the array 'x' 
        (should be set! this is used by e.g. #FTS_CopySpectrum() to make sure the arrays are good) */
        unsigned int    __allocatedLengthx;

        /** new in version 2, number of bits on the AD converter */
        unsigned char    adcBits;

        /** new in version 3. The sensitivity mode of the instrument when collecting
        this spectrum_t */
        unsigned short    sensitivityMode;

        /** new in version 3. The resolution mode of the instrument when collecting
        this spectrum_t  */
        unsigned short    resolutionMode;

        /** new in version 3, true if the spectrum is a rolling average */
        bool            rollingAverage;

        /** new in version 3, null-terminated description of the source used */
        char            source[SPECTRUM_SOURCE_STRLENGTH];

        /** new in version 4, the points where the instrument was power calibrated.
        Copied from the #fts_instrument_t struct of the collecting instrument */
        float            calib_WaveNr[MAX_CALIBRATION_COEFFICIENTS];

        /** new in version 4, the power calibration coefficients of the instrument.
        Copied from the #fts_instrument_t struct of the collecting instrument */
        float            calib_Coeff[MAX_CALIBRATION_COEFFICIENTS];

        /** new in version 4, the number of calibration coefficients.
        Copied from the #fts_instrument_t struct of the collecting instrument */
        unsigned char    calib_coeffNum;

        /** new in version 4, the number of samples that were collected on one period
        of the reference laser. */
        unsigned char    samplesPerReference;

        /** new in version 4, The Optical Path Difference when the collection of the 
        interferogram began. This is usually negative (indicating that the sampling
        began before the Zero Path Difference) */
        float            minOPD_cm;

        /** new in version 4,  The Optical Path Difference when the collection of the 
        interferogram stopped. This is usually positive (indicating that the sampling
        stopped after the Zero Path Difference) */
        float            maxOPD_cm;

        /** new in version 4, flag to determine if the temp/pressure was read from the instrument or not.
        This is an or:ed combination of the FTS_AIR_MEASURE_ flags in FTSData.h.
        E.g. if (air_measureOption & #FTS_AIR_MEASURE_PRESS) != 0 then the pressure was measured by the instrument */
        unsigned int    air_measureOption;

        /** new in version 5, coefficients of the phase polynomial as calculated by phase correction. 
        Input is wavenumber, output is radians. Unused parameters are set to zero */
        float            phasePolynomial[MAX_PHASE_POLY_ORDER + 1];

        /** new in version 5, combination of INTERFEROGRAM_PROPERTY_... flags describing attributes of the
        interferogram */
        int                interferogramProperty;

        /** new in version 5, null-terminated string describing the operator of the measurement */
        char            instr_operator[SPECTRUM_OPERATOR_STRLENGTH];

        /** new in version 6, the average value of the retrieved interferogram */
        float            interferogramOffset;

        /** new in version 6, If apodization was performed, then this is the factor with which the amplitude of the (offset removed)
        interferogram was multiplied with to conserve the total power. */
        float            apodizationPowerCompensation;

        /** new in version 6, the model of the acquiring instrument */
        INSTRUMENT_MODEL    instrument_model;

        /** new in version 6, specifies how the phase values are stored 
        (must equal #PHI_VAL_NONE, #PHI_VAL_ARRAY or #PHI_VAL_IMAGINARY, see \ref sec_complex_spectral_data for information n how these are handled.) */
        unsigned char    phiValueFormat;

        /** new in version 6, the allocated length of the phase array. */
        unsigned int    __allocatedLengthPhi;

        /** new in version 7, the type of instrument that created this spectrum. */
        INSTRUMENT_TYPE        instrumentType;

        /** new in version 7, the trigger mode, if any has been used. */
        TRIGGER_MODE        triggerMode;

        /** new in version 7, the integration time, in milli seconds. Relevant for grating spectra */
        double                integrationTime_ms;

        /** new in version 7, a string describing the model of the instrument. This
        is more versatile than the old 'instrument_model' and should be used instead of this. */
        char            instr_model_str[SPECTRUM_MODEL_STRLENGTH];

        /** new in version 8, the allocated length of the logData array, in bytes. */
        unsigned int    __allocatedLengthLog;

        /** new in version 9, the factors applied in the wavelength correction when calculating the spectrum */
        double            wavelength_correction[4];

        /** new in version 9, the (optional) wavenumber shift. These are stored for reference,
        the x-axis values have been saved with their shifted value. If no shift has been applied 
        in wavenumber space then these are all zeros. */
        double            wavenr_shift[4];

        /** new in version 9, the (optional) wavelength shift. These are stored for reference,
        the x-axis values have been saved with their shifted value. If no shift has been applied 
        in wavelength space then these are all zeros. */
        double            wavelength_shift[4];

        /** new in version 9, the temperature of the detector, in degrees C */
        float            detector_temperature_C;

        /** new in version 9, a temperature from an external temperature sensor, in degrees C */
        float            external_temperature_C;

        /** new in version 9, an air pressure from an external pressure sensor, in hekto Pascal (hPa) */
        float            external_pressure_hPa;

        /** new in version 9, a relative humidity from an external humidity sensor, in percent */
        float            external_relHum_percent;

        /** new in version 6, phase or imaginary data (see \ref sec_complex_spectral_data). Set to NULL by default. */
        float*            phi;

        /** the intensity at each x-axis value */
        float*            I;    

        /** the x-scale, this is only used if xValueFormat is X_VAL_ARRAY, otherwise
        it should be set to NULL. */
        float*            x;

        /** new in version 8, pointer to a log for the data. This can contain any type of data and information about the spectrum_t.
        The length of this array is __allocatedLengthLog bytes. This should be set to NULL if not used. */
        void*            logData;
    } spectrum_t;

    /** An acquisition_buffer_t structure is an auxiliary structure containing information
    *   derived from an interferogram / structure but not actually part of the structure. Different routines
    *   in the acquisition and processing of the interferogram as well as the calculation of a spectrum
    *   from an interferogram sets different parameters here. */
    typedef struct acquisition_buffer_t
    {
        /** This is the dominant wave number in the spectrum, as determined by the cycle counting.
        *   equal to #FTS_WAVELENGTHMETER_ILLEGAL_VALUE if the counting failed, or if counting is disabled. */
        double wavelengthMeter_waveNr;

        /** This is the estimated uncertainty in the wavelengthMeter_waveNr as determined by the cycle counting.
        *   This is only defined if acquisition_buffer_t::wavelengthMeter_waveNr does not equal #FTS_WAVELENGTHMETER_ILLEGAL_VALUE */
        double wavelengthMeterError_waveNr;

        /** If the acquisition_buffer_t::wavelengthMeter_waveNr is #FTS_WAVELENGTHMETER_ILLEGAL_VALUE then this describes the error.
        *   Possible values include \c #FTS_ERROR_SATURATED_INTERFEROGRAM, \c #FTS_ERROR_PARAMETER_ERROR,
        *   \c #FTS_ERROR_WAVELENGTHMETER_IRREGULAR_STRUCTURE, \c #FTS_ERROR_INTERFEROGRAM_INTENSITY_LOW or
        *   \c #FTS_ERROR_WAVELENGTHMETER_ERROR_TOO_LARGE */
        unsigned int wavelengthMeter_ErrorFlag;

        /** wavelengthMeter_Date and wavelengthMeter_Time are the date and time of the interferogram last analyzed by the wavelengthMeter analysis */
        unsigned int wavelengthMeter_Date;
        unsigned int wavelengthMeter_Time;

        /** This is the coherence length of the input light, as determined by the coherence
        *   length analysis on the interferogram. Equal to #FTS_WAVELENGTHMETER_ILLEGAL_VALUE if the 
        *   analysis failed, or if the analysis is disabled. */
        double coherenceLength_cm;

        /** This is the estimated uncertainty in the coherenceLength_cm as determined by the coherence length analysis.
        *   This is only defined if acquisition_buffer_t::coherenceLength_cm does not equal #FTS_WAVELENGTHMETER_ILLEGAL_VALUE */
        double coherenceLengthError_cm;

        /** If the coherenceLength_cm is #FTS_WAVELENGTHMETER_ILLEGAL_VALUE then this describes the error */
        unsigned int coherenceLength_ErrorFlag;

        /** wavelengthMeter_Date and wavelengthMeter_Time are the date and time of the interferogram last analyzed by the wavelengthMeter analysis */
        unsigned int coherenceLength_Date;
        unsigned int coherenceLength_Time;

        /** True if the interferogram is saturated */
        bool interferogram_IsSaturated;

        /** True if the interferogram has such a high input intensity that it is in risk of 
        *      being non-linear */
        bool interferogram_InNonLinearRegime;

        /** this is the maximum peak-to-peak signal in the interferogram, in counts */
        double interferogram_MaxPeakToPeakSignal_counts;

        /** This is the Optical Path Difference (in cm) where the peak-to-peak is maximum */
        double interferogram_MaxPeakToPeakPosition_cm;

        /** True of the interferogram is judged to be from a broad-band source */
        bool interferogram_IsBroadBand;

        /** True of the interferogram is judged to be from a monochromatic source */
        bool interferogram_IsMonochromatic;

        /** This is the estimated noise level of the interferogram. Any data point which deviates less
        *   than interferogramNoiseSigma from the mean is probably noise */
        double interferogramNoiseStdev;

        /** This is the coefficients of the phase polynomial as calculated
        by the phase correction. Polynomial input is in wavenumbers, 
        polynomial output is phase in radians. */
        double phasePolynomial_Coeff0;
        double phasePolynomial_Coeff1;
        double phasePolynomial_Coeff2;

        /** the type of phase correction performed. Must equal one of the  PHASE_CORRECTION_ flags in FTSData.h */
        unsigned char phaseCorrectionType;
    }acquisition_buffer_t;

    /** An autostore_option_t structure describes how interferograms/spectra should be stored 
    to disk while we are collecting data. See \ref sec_auto_store_spectrum_t */
    typedef struct autostore_option_t
    {
        /** If true then the automatic writing is enabled */
        bool            autoStoreOn;

        /** Null-terminated string specifying the directory where the data should be stored.
        This may be encoded in UTF8 or ASCII. */
        char            directory[256];

        /** Null-terminated string specifying the common prefix for the stored files.
        This may be encoded in UTF8 or ASCII. */
        char            filePrefix[32];

        /** Null-terminated string specifying the common suffix (not including
        the file ending) for the stored files. This may be encoded in UTF8 or ASCII. */
        char            fileSuffix[32];

        /** The file format that the files should be written in, must be one of the FILE_FORMAT_...
        flags in FTSData.h (but less than FILE_FORMAT_NUM) */
        unsigned char    fileFormat;

        /** If equal to zero then all files will be stored directly in the specified directory.
        If not zero then sub directories will be created in the specified directory, each 
        containing filesPerSubDir stored files */
        unsigned int    filesPerSubDir;

        /** Counter keeping track of the index of the next file to write, for internal use only */
        unsigned int    nextIndex;

        /** Counter keeping track of the index of the next sub directory to write to, for internal use only */
        unsigned int    nextSubDir;

        /** The spectrum / interferogram will be stored with this 
        x-axis unit if possible */
        X_AXIS_UNIT        desiredXAxisUnit;

        /** The spectrum / interferogram will be stored with this 
        y-axis unit if possible */
        Y_AXIS_UNIT        desiredYAxisUnit;

    }autostore_option_t;

    /** The system_settings_t structure describes system wide setting which does not
    *   apply to any specific instrument or spectrum */
    typedef struct system_settings_t
    {
        /** This is what separates the columns when writing CSV-files, default is semi-colon */
        unsigned char    csv_column_separator;

        /** This is what separates the lines when writing CSV-files, default is carriage-return + line-feed */
        char    csv_newline[8];

        /** This is the time from the last acquisition until the system turns of the instrument
        *   zero corresponds to no turning of the instrument. Default is 15 minutes */
        unsigned long    instrument_sleeptime_minutes;

        /** This is the name of the current operator of the instrument. This will be set in 
        *   all the interferograms / spectra that are collected. */
        char            instr_operator[SPECTRUM_OPERATOR_STRLENGTH];

        /** This is the directory to the users personal library. 
        *   This may be encoded in UTF8 or ASCII. */
        char            operator_library_path[MAX_PATH_LENGTH];

        /** If true then the messages from the trace-log will be streamed to the file 'TraceLog_NN.txt' 
        *   (where NN is a unique number) in the directory specified by operator_library_path */
        bool            log_traceLog_messages;

        /** if useGPU is set to true then some of the computations will be offloaded to the GPU
        *   instead of being performed on the CPU. This should be set to false if running on a Windows
        *    version older than Windows 7 or if your GPU is slow */
        bool            useGPU;

    }system_settings_t;

    // ------------------------------------------------------------------------
    // ------------------------------ DATA ------------------------------------
    // ------------------------------------------------------------------------

    /** The FTSData structure contains the properties of the connected instruments, the 
    *   settings for acquiring data from them and the system settings. There is one global
    *   FTSData structure which can be accessed through #FTS_GetFTSData(). */
    typedef struct FTSData
    {
        /** This array of #fts_instrument_t describes the properties of each of 
        *   the connected instruments. Only FTSData::instrument_num elements in
        *   this array are defined. */
        fts_instrument_t                fts_instrument[MAX_SPECTROMETER_NUM];

        /** These are the acquisition options for each of the connected instruments.
        *   Only FTSData::instrument_num elements in this array are defined.  */
        acquisition_option_t            acquisition_options[MAX_SPECTROMETER_NUM];

        /** These are the measured results.
        *   Only FTSData::instrument_num elements in this array are defined. */
        acquisition_buffer_t            data[MAX_SPECTROMETER_NUM];

        /** These are the options for automatic storing of collected interferograms.
        *   Only FTSData::instrument_num elements in this array are defined. 
        *   See \ref sec_auto_store_spectrum_t */
        autostore_option_t              auto_store_interferogram[MAX_SPECTROMETER_NUM];

        /** These are the options for automatic storing of calculated spectra.
        *   Only FTSData::instrument_num elements in this array are defined. 
        *   See \ref sec_auto_store_spectrum_t */
        autostore_option_t              auto_store_spectrum[MAX_SPECTROMETER_NUM];

        /** This is the system settings which applies to all instruments and all channels. */
        system_settings_t               system_settings;

        /** The number of spectrometers that are currently connected to the system...
        *   this is also the number of valid elements in the arrays above */
        unsigned int                    instrument_num;

    } FTSData;

    typedef struct FTSTraceData
    {
        /** This is an array of memory buffers, useful for e.g. copying data to and from
        *   the library. Initially these are initialized of size 100.. 
        *   Access via routines such as #FTS_Trace_Clear() or #FTS_Trace_Allocate() */
        spectrum_t* trace_buffer[MAX_TRACEBUFFER_NUM];
    } FTSTraceData;

    // -------------------------------------------------------------------------------
    // --------------------------- EXPORTED FUNCTIONS --------------------------------
    // -------------------------------------------------------------------------------

    /** Retrieve the address of the global data structure */
    FTSData* __cdecl   FTS_GetFTSData();
    FTSData* __stdcall FTS_GetFTSData_std();


    /** Returns the number of instruments connected to the system .
    To update the list call #FTS_InitializeSpectrometers() */
    unsigned int __cdecl    FTS_GetInstrumentNum();
    unsigned int __stdcall    FTS_GetInstrumentNum_std();

    /** Retrieve the address of instrument data for the given spectrometer 
    *   @param specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers */
    fts_instrument_t*  __cdecl   FTS_GetInstrumentData(SpectrometerIndex specIndex);
    fts_instrument_t*  __stdcall FTS_GetInstrumentData_std(SpectrometerIndex specIndex);

    /** Retrieve the address of the acquisition options for the given spectrometer 
    *   @param specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers */
    acquisition_option_t*  __cdecl   FTS_GetAcquisitionOptions(SpectrometerIndex specIndex);
    acquisition_option_t*  __stdcall FTS_GetAcquisitionOptions_std(SpectrometerIndex specIndex);

    /** Retrieve the address of the acquisition options for the given spectrometer 
    *   @param specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers */
    autostore_option_t * __cdecl   FTS_GetAutoStoreSpectrumOptions(SpectrometerIndex specIndex);
    autostore_option_t * __cdecl   FTS_GetAutoStoreInterferogramOptions(SpectrometerIndex specIndex);
    autostore_option_t * __stdcall FTS_GetAutoStoreSpectrumOptions_std(SpectrometerIndex specIndex);
    autostore_option_t * __stdcall FTS_GetAutoStoreInterferogramOptions_std(SpectrometerIndex specIndex);

    /** Restores the acquisition options for the specified spectrometer to the default values 
    *   @param specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers */
    void __cdecl        FTS_SetAcquisitionOptions_Default(SpectrometerIndex specIndex);
    void __stdcall      FTS_SetAcquisitionOptions_Default_std(SpectrometerIndex specIndex);

    /** Restores the auto store options for both interferograms and spectra to default
    *   values for the specified spectrometer.
    *   @param specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers */
    void __cdecl        FTS_SetAutoStoreOptions_Default(SpectrometerIndex specIndex);
    void __stdcall      FTS_SetAutoStoreOptions_Default_std(SpectrometerIndex specIndex);


    /** \addtogroup grp_acquisition_options Acquisition Options
    *   The FTS_SetAcquisitionOption_... routines sets one or more parameters in the acquisition_options_t structure
    *   that belongs to the specified instrument. 
    *   They all return #FTS_SUCCESS if all parameters were set and #FTS_ERROR_PARAMETER_ERROR if any of the parameters are not valid 
    *  @{
    */

    /** Set the type of, and parameters for, the smoothing that should be applied to the interferograms */
    int __cdecl   FTS_SetAcquisitionOption_InterferogramSmoothing(SpectrometerIndex specIndex, unsigned short iSmoothingType, unsigned short param1, unsigned short param2);
    int __stdcall FTS_SetAcquisitionOption_InterferogramSmoothing_std(SpectrometerIndex specIndex, unsigned short iSmoothingType, unsigned short param1, unsigned short param2);

    /** Set the type of, and parameters for, the smoothing that should be applied to the spectra */
    int __cdecl   FTS_SetAcquisitionOption_SpectrumSmoothing(SpectrometerIndex specIndex, unsigned short sSmoothingType, unsigned short param1, unsigned short param2);
    int __stdcall FTS_SetAcquisitionOption_SpectrumSmoothing_std(SpectrometerIndex specIndex, unsigned short sSmoothingType, unsigned short param1, unsigned short param2);

    /** Set the type of apodization that should be applied to the interferograms before calculating the spectrum. */
    int __cdecl   FTS_SetAcquisitionOption_ApodisationType(SpectrometerIndex specIndex, unsigned char apodisationType);
    int __stdcall FTS_SetAcquisitionOption_ApodisationType_std(SpectrometerIndex specIndex, unsigned char apodisationType);

    /** Set the type of phase correction that should be applied to the interferograms or the spectra. */
    int __cdecl   FTS_SetAcquisitionOption_PhaseCorrectionType(SpectrometerIndex specIndex, unsigned char phaseCorrType);
    int __stdcall FTS_SetAcquisitionOption_PhaseCorrectionType_std(SpectrometerIndex specIndex, unsigned char phaseCorrType);

    /** Enable/Disable averaging the spectra before displaying them. If the number of spectra to average > 1 then the spectra will
    be averaged within the FTSLib and an averaged spectrum be read out using #FTS_GetLastSpectrum. */
    int __cdecl   FTS_SetAcquisitionOption_AverageSpectrum(SpectrometerIndex specIndex, unsigned int averageSpecNum);
    int __stdcall FTS_SetAcquisitionOption_AverageSpectrum_std(SpectrometerIndex specIndex, unsigned int averageSpecNum);

    /** Enable/disable rolling average. */
    int __cdecl   FTS_SetAcquisitionOption_RollingAverage(SpectrometerIndex specIndex, bool rollingAverage);
    int __stdcall FTS_SetAcquisitionOption_RollingAverage_std(SpectrometerIndex specIndex, bool rollingAverage);

    /** Set the lower edge of the wavenumber range to cut the calculated spectra to. 
    This must be higher than the minimum wavenumber of the instrument. */
    int __cdecl   FTS_SetAcquisitionOption_MinWaveNr(SpectrometerIndex specIndex, float value);
    int __stdcall FTS_SetAcquisitionOption_MinWaveNr_std(SpectrometerIndex specIndex, float value);

    /** Set the upper edge of the wavenumber range to cut the calculated spectra to. 
    *   This must be lower than the minimum wavenumber of the instrument. */
    int __cdecl   FTS_SetAcquisitionOption_MaxWaveNr(SpectrometerIndex specIndex, float value);
    int __stdcall FTS_SetAcquisitionOption_MaxWaveNr_std(SpectrometerIndex specIndex, float value);

    /** Sets the zero-filling to apply to the interferograms during calculation of the spectra.
    *   zeroFillFactor = 0 corresponds to no zero filling being applied.
    *   zeroFillFactor = 1 corresponds to filling the interferograms with as many zeros as original values,
    *   zeroFillFactor = 2 corresponds to filling the interferograms with three times as many zeros as original values */
    int __cdecl   FTS_SetAcquisitionOption_ZeroFillFactor(SpectrometerIndex specIndex, unsigned char zeroFillFactor);
    int __stdcall FTS_SetAcquisitionOption_ZeroFillFactor_std(SpectrometerIndex specIndex, unsigned char zeroFillFactor);

    /** Sets the resolution mode of the instrument. The mode number must be 0 to 'motor_resolutionModeNum' (as defined in the #fts_instrument_t).
    *   The provided callback will be called when the operation is complete. The callback may be null. */
    int __cdecl   FTS_SetAcquisitionOption_ResolutionMode(SpectrometerIndex specIndex, unsigned short travelMode, void (__cdecl *pGUIFunction)(unsigned short, unsigned int, unsigned int));
    int __stdcall FTS_SetAcquisitionOption_ResolutionMode_std(SpectrometerIndex specIndex, unsigned short travelMode, void (__stdcall *pGUIFunction)(unsigned short, unsigned int, unsigned int));

    /** Sets the sensitivity mode of the instrument. The mode number must be 0 to 'motor_sensitivityModeNum' (as defined in the #fts_instrument_t).
    *   The provided callback will be called when the operation is complete. The callback may be null. */
    int __cdecl   FTS_SetAcquisitionOption_SensitivityMode(SpectrometerIndex specIndex, unsigned short sensitivityMode, void (__cdecl *pGUIFunction)(unsigned short, unsigned int, unsigned int));
    int __stdcall FTS_SetAcquisitionOption_SensitivityMode_std(SpectrometerIndex specIndex, unsigned short sensitivityMode, void (__stdcall *pGUIFunction)(unsigned short, unsigned int, unsigned int));

    /** Sets the detector temperature mode of the instrument. The mode number must be 0 to 'det_temperatureModeNum' (as defined in the #fts_instrument_t).
    *   Notice: this can only be called if control of the detector temperature is implemented in the device, 
    *   (det_TemperatureControlImplemented is set to true in #fts_instrument_t).
    *   The provided callback will be called when the operation is complete. The callback may be null. */
    int __cdecl   FTS_SetAcquisitionOption_DetectorTemperatureMode(SpectrometerIndex specIndex, unsigned short sensitivityMode, void (__cdecl *pGUIFunction)(unsigned short, unsigned int, unsigned int));
    int __stdcall FTS_SetAcquisitionOption_DetectorTemperatureMode_std(SpectrometerIndex specIndex, unsigned short sensitivityMode, void (__stdcall *pGUIFunction)(unsigned short, unsigned int, unsigned int));


    /** This sets the motor travel range for the given spectrometer.
    *   This sets the mode synchronously and will not return until the mode has been set.
    *   @return #FTS_SUCCESS on successful setting of the parameters
    *   @return #FTS_ERROR_PARAMETER_ERROR if the index or the resolution mode is invalid. */
    int __cdecl   FTS_SetResolutionMode(SpectrometerIndex specIndex, unsigned short rangeMode);
    int __stdcall FTS_SetResolutionMode_std(SpectrometerIndex specIndex, unsigned short rangeMode);

    /** This sets the motor speed mode for the given spectrometer.
    *   This sets the mode synchronously and will not return until the mode has been set.
    *   @return #FTS_SUCCESS on successful setting of the parameters.
    *   @return #FTS_ERROR_PARAMETER_ERROR if the index or the sensitivity mode is invalid. */
    int __cdecl   FTS_SetSensitivityMode(SpectrometerIndex specIndex, unsigned short sensitivityMode);
    int __stdcall FTS_SetSensitivityMode_std(SpectrometerIndex specIndex, unsigned short sensitivityMode);

    /** This sets the detector temperature mode for the given spectrometer.
    *   This sets the mode synchronously and will not return until the mode has been set.
    *   @return #FTS_SUCCESS on successful setting of the parameters
    *   @return #FTS_ERROR_PARAMETER_ERROR if the index or the mode is invalid. */
    int __cdecl   FTS_SetDetectorTemperatureMode(SpectrometerIndex specIndex, unsigned short pdTempMode);
    int __stdcall FTS_SetDetectorTemperatureMode_std(SpectrometerIndex specIndex, unsigned short pdTempMode);

    /** Enable/Disable automatic determination of the gain level to use for the device. */
    int __cdecl   FTS_SetAcquisitionOption_AutoGain(SpectrometerIndex specIndex, bool autoGain);
    int __stdcall FTS_SetAcquisitionOption_AutoGain_std(SpectrometerIndex specIndex, bool autoGain);

    /** Set the gain index to use for the detector(s) of the device. Calling this routine will set auto gain to false */
    int __cdecl   FTS_SetAcquisitionOption_SingleGain(SpectrometerIndex specIndex, unsigned short gainIndex);
    int __stdcall FTS_SetAcquisitionOption_SingleGain_std(SpectrometerIndex specIndex, unsigned short gainIndex);

    /** Enable automatic cycle-counting of measured interferograms. If this is set to true, then
    *   parameters 'wavelengthMeter_waveNr', 'wavelengthMeterError_waveNr','wavelengthMeter_ErrorFlag', wavelengthMeter_Date' and 'wavelengthMeter_Time'
    *   of the acquisition_buffer_t will update for every interferogram that is collected */
    int __cdecl   FTS_SetAcquisitionOption_CycleCounting(SpectrometerIndex specIndex, bool cycleCounting);
    int __stdcall FTS_SetAcquisitionOption_CycleCounting_std(SpectrometerIndex specIndex, bool cycleCounting);

    /** Enable automatic coherence length determination on the measured interferograms. If this is set to true, then
    *   parameters 'coherenceLength_cm', 'coherenceLengthError_cm','coherenceLength_ErrorFlag', coherenceLength_Date' and 'coherenceLength_Time'
    *   of the acquisition_buffer_t will update for every interferogram that is collected */
    int __cdecl   FTS_SetAcquisitionOption_CoherenceAnalysis(SpectrometerIndex specIndex, bool coherenceAnalysis);
    int __stdcall FTS_SetAcquisitionOption_CoherenceAnalysis_std(SpectrometerIndex specIndex, bool coherenceAnalysis);

    /** Set how the environment parameters are determined during the acquisitions, 
    *   see the description of air_measureOption in spectrum_t */
    int __cdecl   FTS_SetAcquisitionOption_AirMeasureOption(SpectrometerIndex specIndex, unsigned int air_measureOption);
    int __stdcall FTS_SetAcquisitionOption_AirMeasureOption_std(SpectrometerIndex specIndex, unsigned int air_measureOption);

    /** Sets the (minimum) time between two acquisitions from the device. 
    *   If this is set to zero then the acquisitions will run as fast as they can. */
    int __cdecl   FTS_SetAcquisitionOption_AcquisitionInterval(SpectrometerIndex specIndex, unsigned int acquisitionInterval_ms);
    int __stdcall FTS_SetAcquisitionOption_AcquisitionInterval_std(SpectrometerIndex specIndex, unsigned int acquisitionInterval_ms);

    ///@}


    /** \addtogroup grp_auto_store Settings for Automatic Data Storage
    *   The FTS_SetAcquisitionOption_AutoStoreInterferogram_... routines sets the options for 
    *   automatic storing of the collected interferograms.
    *   The FTS_SetAcquisitionOption_AutoStoreSpectrum_... routines sets the options for 
    *   automatic storing of the collected interferograms.
    *  @{
    */
    /**  */
    int __cdecl   FTS_SetAcquisitionOption_AutoStoreInterferogram_On(SpectrometerIndex specIndex, bool on);
    int __cdecl   FTS_SetAcquisitionOption_AutoStoreInterferogram_Directory(SpectrometerIndex specIndex, char* directory);
    int __cdecl   FTS_SetAcquisitionOption_AutoStoreInterferogram_FileFormat(SpectrometerIndex specIndex, unsigned char    fileFormat);
    int __cdecl   FTS_SetAcquisitionOption_AutoStoreInterferogram_FilePrefix(SpectrometerIndex specIndex, char* filePrefix);
    int __cdecl   FTS_SetAcquisitionOption_AutoStoreInterferogram_FileSuffix(SpectrometerIndex specIndex, char* fileSuffix);
    int __cdecl   FTS_SetAcquisitionOption_AutoStoreInterferogram_FilesPerSubDir(SpectrometerIndex specIndex, unsigned int filesPerSubDir);
    int __cdecl   FTS_SetAcquisitionOption_AutoStoreInterferogram_XAxisUnit(SpectrometerIndex specIndex, X_AXIS_UNIT unit);
    int __cdecl   FTS_SetAcquisitionOption_AutoStoreInterferogram_YAxisUnit(SpectrometerIndex specIndex, Y_AXIS_UNIT unit);

    int __stdcall FTS_SetAcquisitionOption_AutoStoreInterferogram_On_std(SpectrometerIndex specIndex, bool on);
    int __stdcall FTS_SetAcquisitionOption_AutoStoreInterferogram_Directory_std(SpectrometerIndex specIndex, char* directory);
    int __stdcall FTS_SetAcquisitionOption_AutoStoreInterferogram_FileFormat_std(SpectrometerIndex specIndex, unsigned char    fileFormat);
    int __stdcall FTS_SetAcquisitionOption_AutoStoreInterferogram_FilePrefix_std(SpectrometerIndex specIndex, char* filePrefix);
    int __stdcall FTS_SetAcquisitionOption_AutoStoreInterferogram_FileSuffix_std(SpectrometerIndex specIndex, char* fileSuffix);
    int __stdcall FTS_SetAcquisitionOption_AutoStoreInterferogram_FilesPerSubDir_std(SpectrometerIndex specIndex, unsigned int filesPerSubDir);
    int __stdcall FTS_SetAcquisitionOption_AutoStoreInterferogram_XAxisUnit_std(SpectrometerIndex specIndex, X_AXIS_UNIT unit);
    int __stdcall FTS_SetAcquisitionOption_AutoStoreInterferogram_YAxisUnit_std(SpectrometerIndex specIndex, Y_AXIS_UNIT unit);

    int __cdecl   FTS_SetAcquisitionOption_AutoStoreSpectrum_On(SpectrometerIndex specIndex, bool on);
    int __cdecl   FTS_SetAcquisitionOption_AutoStoreSpectrum_Directory(SpectrometerIndex specIndex, char* directory);
    int __cdecl   FTS_SetAcquisitionOption_AutoStoreSpectrum_FileFormat(SpectrometerIndex specIndex, unsigned char    fileFormat);
    int __cdecl   FTS_SetAcquisitionOption_AutoStoreSpectrum_FilePrefix(SpectrometerIndex specIndex, char* filePrefix);
    int __cdecl   FTS_SetAcquisitionOption_AutoStoreSpectrum_FileSuffix(SpectrometerIndex specIndex, char* fileSuffix);
    int __cdecl   FTS_SetAcquisitionOption_AutoStoreSpectrum_FilesPerSubDir(SpectrometerIndex specIndex, unsigned int filesPerSubDir);
    int __cdecl   FTS_SetAcquisitionOption_AutoStoreSpectrum_XAxisUnit(SpectrometerIndex specIndex, X_AXIS_UNIT unit);
    int __cdecl   FTS_SetAcquisitionOption_AutoStoreSpectrum_YAxisUnit(SpectrometerIndex specIndex, Y_AXIS_UNIT unit);

    int __stdcall FTS_SetAcquisitionOption_AutoStoreSpectrum_On_std(SpectrometerIndex specIndex, bool on);
    int __stdcall FTS_SetAcquisitionOption_AutoStoreSpectrum_Directory_std(SpectrometerIndex specIndex, char* directory);
    int __stdcall FTS_SetAcquisitionOption_AutoStoreSpectrum_FileFormat_std(SpectrometerIndex specIndex, unsigned char    fileFormat);
    int __stdcall FTS_SetAcquisitionOption_AutoStoreSpectrum_FilePrefix_std(SpectrometerIndex specIndex, char* filePrefix);
    int __stdcall FTS_SetAcquisitionOption_AutoStoreSpectrum_FileSuffix_std(SpectrometerIndex specIndex, char* fileSuffix);
    int __stdcall FTS_SetAcquisitionOption_AutoStoreSpectrum_FilesPerSubDir_std(SpectrometerIndex specIndex, unsigned int filesPerSubDir);
    int __stdcall FTS_SetAcquisitionOption_AutoStoreSpectrum_XAxisUnit_std(SpectrometerIndex specIndex, X_AXIS_UNIT unit);
    int __stdcall FTS_SetAcquisitionOption_AutoStoreSpectrum_YAxisUnit_std(SpectrometerIndex specIndex, Y_AXIS_UNIT unit);

    ///@}


    /** \addtogroup grp_acquisition_options Acquisition Options
    *   The FTS_GetAcquisitionOption_... routines retrieves one parameter from the acquisition_options_t structure
    *   that belongs to the specified instrument. 
    *   They all return -1 if specIndex is out of range
    *  @{
    */
    int __cdecl   FTS_GetAcquisitionOption_InterferogramSmoothing(SpectrometerIndex specIndex, unsigned short *param1, unsigned short *param2);
    int __cdecl   FTS_GetAcquisitionOption_SpectrumSmoothing(SpectrometerIndex specIndex, unsigned short *param1, unsigned short *param2);
    int __cdecl   FTS_GetAcquisitionOption_ApodisationType(SpectrometerIndex specIndex);
    int __cdecl   FTS_GetAcquisitionOption_PhaseCorrectionType(SpectrometerIndex specIndex);
    int __cdecl   FTS_GetAcquisitionOption_AverageSpectrum(SpectrometerIndex specIndex);
    bool __cdecl   FTS_GetAcquisitionOption_RollingAverage(SpectrometerIndex specIndex);
    float __cdecl   FTS_GetAcquisitionOption_MinWaveNr(SpectrometerIndex specIndex);
    float __cdecl   FTS_GetAcquisitionOption_MaxWaveNr(SpectrometerIndex specIndex);
    int __cdecl   FTS_GetAcquisitionOption_ZeroFillFactor(SpectrometerIndex specIndex);
    int __cdecl   FTS_GetAcquisitionOption_ResolutionMode(SpectrometerIndex specIndex);
    int __cdecl   FTS_GetAcquisitionOption_SensitivityMode(SpectrometerIndex specIndex);
    int __cdecl   FTS_GetAcquisitionOption_DetectorTemperatureMode(SpectrometerIndex specIndex);

    /** @return true if automatic gain is turned on for the given spectrometer. */
    bool __cdecl   FTS_GetAcquisitionOption_AutoGain(SpectrometerIndex specIndex);

    /** @return the gain <b>index</b> to use for the next acquisition from the spectrometer with the given index */
    int __cdecl   FTS_GetAcquisitionOption_SingleGainLevel(SpectrometerIndex specIndex);
    int __cdecl   FTS_GetAcquisitionOption_CycleCounting(SpectrometerIndex specIndex);
    int __cdecl   FTS_GetAcquisitionOption_CoherenceAnalysis(SpectrometerIndex specIndex);
    int __cdecl   FTS_GetAcquisitionOption_AirMeasureOption(SpectrometerIndex specIndex);
    int __cdecl   FTS_GetAcquisitionOption_AcquisitionInterval(SpectrometerIndex specIndex);

    int __stdcall FTS_GetAcquisitionOption_InterferogramSmoothing_std(SpectrometerIndex specIndex, unsigned short *param1, unsigned short *param2);
    int __stdcall FTS_GetAcquisitionOption_SpectrumSmoothing_std(SpectrometerIndex specIndex, unsigned short *param1, unsigned short *param2);
    int __stdcall FTS_GetAcquisitionOption_ApodisationType_std(SpectrometerIndex specIndex);
    int __stdcall FTS_GetAcquisitionOption_PhaseCorrectionType_std(SpectrometerIndex specIndex);
    int __stdcall FTS_GetAcquisitionOption_AverageSpectrum_std(SpectrometerIndex specIndex);
    bool __stdcall FTS_GetAcquisitionOption_RollingAverage_std(SpectrometerIndex specIndex);
    float __stdcall FTS_GetAcquisitionOption_MinWaveNr_std(SpectrometerIndex specIndex);
    float __stdcall FTS_GetAcquisitionOption_MaxWaveNr_std(SpectrometerIndex specIndex);
    int __stdcall FTS_GetAcquisitionOption_ZeroFillFactor_std(SpectrometerIndex specIndex);
    int __stdcall FTS_GetAcquisitionOption_ResolutionMode_std(SpectrometerIndex specIndex);
    int __stdcall FTS_GetAcquisitionOption_SensitivityMode_std(SpectrometerIndex specIndex);
    int __stdcall  FTS_GetAcquisitionOption_DetectorTemperatureMode_std(SpectrometerIndex specIndex);
    bool __stdcall FTS_GetAcquisitionOption_AutoGain_std(SpectrometerIndex specIndex);
    int __stdcall FTS_GetAcquisitionOption_SingleGainLevel_std(SpectrometerIndex specIndex);
    int __stdcall FTS_GetAcquisitionOption_CycleCounting_std(SpectrometerIndex specIndex);     
    int __stdcall FTS_GetAcquisitionOption_CoherenceAnalysis_std(SpectrometerIndex specIndex);
    int __stdcall FTS_GetAcquisitionOption_AirMeasureOption_std(SpectrometerIndex specIndex);
    int __stdcall FTS_GetAcquisitionOption_AcquisitionInterval_std(SpectrometerIndex specIndex);
    ///@}

    /** \addtogroup grp_instrument_property Getting the Properties of the Instrument
    *   The FTS_GetInstrumentProperty_... routines retrieves the properties of the instrument itself
    *  @{
    */

    /** Read out the serial number of the device with the specified index.
    *   @param specIndex the index of the spectrometer.
    *   @param[out] serial will on return be filled with the serial of the device, this must be a pointer
    *   to a buffer able to hold at least #SERIALNUMBER_LENGTH characters */
    void            __cdecl   FTS_GetInstrumentProperty_Serial(SpectrometerIndex specIndex, char* serial);

    /** Read out the model number of the device */
    INSTRUMENT_MODEL __cdecl   FTS_GetInstrumentProperty_Model(SpectrometerIndex specIndex);

    /** Read out the model of the device as a string.
    *   @param[out] modelStr will on return be filled with the model of the device. This must be a pointer
    *   to a buffer able to hold at least #SPECTRUM_MODEL_STRLENGTH characters */
    void            __cdecl   FTS_GetInstrumentProperty_ModelStr(SpectrometerIndex specIndex, char* modelStr);
    float           __cdecl   FTS_GetInstrumentProperty_MinWavenumber(SpectrometerIndex specIndex);
    float           __cdecl   FTS_GetInstrumentProperty_MaxWavenumber(SpectrometerIndex specIndex);
    unsigned short  __cdecl   FTS_GetInstrumentProperty_ADCBits(SpectrometerIndex specIndex);
    float           __cdecl   FTS_GetInstrumentProperty_CalibrationPower(SpectrometerIndex specIndex, unsigned int index);
    float           __cdecl   FTS_GetInstrumentProperty_CalibrationWavenumber(SpectrometerIndex specIndex, unsigned int index);
    unsigned int    __cdecl   FTS_GetInstrumentProperty_CalibrationCoefficients(SpectrometerIndex specIndex);
    unsigned short  __cdecl   FTS_GetInstrumentProperty_SensitivityModeNum(SpectrometerIndex specIndex);
    float           __cdecl   FTS_GetInstrumentProperty_MotorSpeed(SpectrometerIndex specIndex, unsigned int index);
    unsigned short  __cdecl   FTS_GetInstrumentProperty_ResolutionModeNum(SpectrometerIndex specIndex);
    float           __cdecl   FTS_GetInstrumentProperty_MotorTravel(SpectrometerIndex specIndex, unsigned int index);
    double          __cdecl   FTS_GetInstrumentProperty_ReferenceWavelength(SpectrometerIndex specIndex);
    float           __cdecl   FTS_GetInstrumentProperty_SamplesPerReference(SpectrometerIndex specIndex);
    float           __cdecl   FTS_GetInstrumentProperty_Temperature(SpectrometerIndex specIndex);
    float           __cdecl   FTS_GetInstrumentProperty_AtmosphericPressure(SpectrometerIndex specIndex);
    float           __cdecl   FTS_GetInstrumentProperty_RelativeHumidity(SpectrometerIndex specIndex);
    unsigned int    __cdecl   FTS_GetInstrumentProperty_InstrumentStatus(SpectrometerIndex specIndex);
    unsigned int    __cdecl   FTS_GetInstrumentProperty_DetailedInstrumentStatus(SpectrometerIndex specIndex);
    unsigned int    __cdecl   FTS_GetInstrumentProperty_CommunicationStatus(SpectrometerIndex specIndex);
    float           __cdecl   FTS_GetInstrumentProperty_Resolution_Wnr(SpectrometerIndex specIndex);
    unsigned int    __cdecl   FTS_GetInstrumentProperty_AvailableGainLevels(SpectrometerIndex specIndex, double* GainLevels);
    unsigned int    __cdecl   FTS_GetInstrumentProperty_AvailableGainIndices(SpectrometerIndex specIndex, int* GainLevels);
    bool            __cdecl   FTS_GetInstrumentProperty_HasTemperatureSensor(SpectrometerIndex specIndex);
    bool            __cdecl   FTS_GetInstrumentProperty_HasPressureSensor(SpectrometerIndex specIndex);
    bool            __cdecl   FTS_GetInstrumentProperty_HasHumiditySensor(SpectrometerIndex specIndex);
    bool            __cdecl   FTS_GetInstrumentProperty_HasDetectorTemperatureControl(SpectrometerIndex specIndex);
    float           __cdecl   FTS_GetInstrumentProperty_CurrentDetectorTemperature(SpectrometerIndex specIndex);
    unsigned short  __cdecl   FTS_GetInstrumentProperty_DetectorTemperatureModeNum(SpectrometerIndex specIndex);
    bool            __cdecl   FTS_GetInstrumentProperty_HasBlockingAimingBeam(SpectrometerIndex specIndex);
    bool            __cdecl   FTS_GetInstrumentProperty_IsBlockingAimingBeamOn(SpectrometerIndex specIndex);

    void                __stdcall FTS_GetInstrumentProperty_Serial_std(SpectrometerIndex specIndex, char* serial);
    INSTRUMENT_MODEL    __stdcall FTS_GetInstrumentProperty_Model_std(SpectrometerIndex specIndex);
    void            __stdcall   FTS_GetInstrumentProperty_ModelStr_std(SpectrometerIndex specIndex, char* modelStr);
    float            __stdcall FTS_GetInstrumentProperty_MinWavenumber_std(SpectrometerIndex specIndex);
    float            __stdcall FTS_GetInstrumentProperty_MaxWavenumber_std(SpectrometerIndex specIndex);
    unsigned short    __stdcall FTS_GetInstrumentProperty_ADCBits_std(SpectrometerIndex specIndex);
    float            __stdcall FTS_GetInstrumentProperty_CalibrationPower_std(SpectrometerIndex specIndex, unsigned int index);
    float            __stdcall FTS_GetInstrumentProperty_CalibrationWavenumber_std(SpectrometerIndex specIndex, unsigned int index);
    unsigned int        __stdcall FTS_GetInstrumentProperty_CalibrationCoefficients_std(SpectrometerIndex specIndex);
    unsigned short    __stdcall FTS_GetInstrumentProperty_SensitivityModeNum_std(SpectrometerIndex specIndex);
    float            __stdcall FTS_GetInstrumentProperty_MotorSpeed_std(SpectrometerIndex specIndex, unsigned int index);
    unsigned short    __stdcall FTS_GetInstrumentProperty_ResolutionModeNum_std(SpectrometerIndex specIndex);
    float            __stdcall FTS_GetInstrumentProperty_MotorTravel_std(SpectrometerIndex specIndex, unsigned int index);
    double            __stdcall FTS_GetInstrumentProperty_ReferenceWavelength_std(SpectrometerIndex specIndex);
    float            __stdcall FTS_GetInstrumentProperty_SamplesPerReference_std(SpectrometerIndex specIndex);
    float            __stdcall FTS_GetInstrumentProperty_Temperature_std(SpectrometerIndex specIndex);
    float            __stdcall FTS_GetInstrumentProperty_AtmosphericPressure_std(SpectrometerIndex specIndex);
    float            __stdcall FTS_GetInstrumentProperty_RelativeHumidity_std(SpectrometerIndex specIndex);
    unsigned int        __stdcall FTS_GetInstrumentProperty_InstrumentStatus_std(SpectrometerIndex specIndex);
    unsigned int        __stdcall FTS_GetInstrumentProperty_CommunicationStatus_std(SpectrometerIndex specIndex);
    unsigned int        __stdcall   FTS_GetInstrumentProperty_VSEStatus_std(SpectrometerIndex specIndex);
    float            __stdcall FTS_GetInstrumentProperty_Resolution_Wnr_std(SpectrometerIndex specIndex);
    unsigned int        __stdcall FTS_GetInstrumentProperty_AvailableGainLevels_std(SpectrometerIndex specIndex, double* GainLevels);
    unsigned int        __stdcall FTS_GetInstrumentProperty_AvailableGainIndices_std(SpectrometerIndex specIndex, int* GainLevels);
    bool            __stdcall FTS_GetInstrumentProperty_HasTemperatureSensor_std(SpectrometerIndex specIndex);
    bool            __stdcall FTS_GetInstrumentProperty_HasPressureSensor_std(SpectrometerIndex specIndex);
    bool            __stdcall FTS_GetInstrumentProperty_HasHumiditySensor_std(SpectrometerIndex specIndex);
    bool            __stdcall FTS_GetInstrumentProperty_HasDetectorTemperatureControl_std(SpectrometerIndex specIndex);
    float            __stdcall FTS_GetInstrumentProperty_CurrentDetectorTemperature_std(SpectrometerIndex specIndex);
    unsigned short    __stdcall FTS_GetInstrumentProperty_DetectorTemperatureModeNum_std(SpectrometerIndex specIndex);
    bool            __stdcall   FTS_GetInstrumentProperty_HasBlockingAimingBeam_std(SpectrometerIndex specIndex);
    bool            __stdcall   FTS_GetInstrumentProperty_IsBlockingAimingBeamOn_std(SpectrometerIndex specIndex);
    ///@}

    /** \addtogroup grp_async_acquisition Asynchronous Data Acquisition
    See Also \ref sec_asyncCollection
    *  @{
    */

    /** Retrieve the address of acquired data for the given spectrometer */
    acquisition_buffer_t*  __cdecl   FTS_GetCollectedData(SpectrometerIndex specIndex);
    acquisition_buffer_t*  __stdcall FTS_GetCollectedData_std(SpectrometerIndex specIndex);

    /** Retrieve a copy of the last spectrum acquired asynchronously.
    *   @param[in] specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers.
    *   @param[out] spec will on successful return be filled with the data in the last acquired spectrum.
    *   Notice that this must be allocated before calling this routine with a large enough buffers spectrum_t::I
    *   (to hold the spectral data) and spectrum_t::x (to hold the x-axis data), i.e. spectrum_t::__allocatedLengthI and 
    *   spectrum_t::__allocatedLengthX must both be at least as long as the size of the spectrum to read out.
    *   This may be \c NULL if only the size of the last spectrum should be retrieved.
    *   @param[out] length will on successful return be filled with the length of the last acquired spectrum.
    *   
    *   @return #FTS_SUCCESS if the last acquired spectrum was successfully retrieved
    *   @return #FTS_ERROR_NODATA if no spectrum has yet been calculated asynchronously.
    *   @return #FTS_ERROR_OUTOFMEMORY if the supplied spectrum_t spec was too small to fit in the data.
    *   
    *   To retrieve the size of the last spectrum, set spec to \c NULL and retrieve the value
    *   of the parameter length upon successful return.
    */
    int __cdecl   FTS_GetLastSpectrum(SpectrometerIndex specIndex, spectrum_t* spec, unsigned int* length);
    int __stdcall FTS_GetLastSpectrum_std(SpectrometerIndex specIndex, spectrum_t* spec, unsigned int* length);

    /** \deprecated
    *   Deprecated!  This method is deprecated and should not be used. As of version v2.50 can spectra be stored
    *   at non-even intervals then stored using an x-axis array which cannot be read out here. Use #FTS_GetLastSpectrumArrays instead
    *   Retrieve a copy of the last spectrum acquired asynchronously.
    *   @param[in] specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers.
    *   @param[out] spec will on successful return be filled with the data in the last acquired spectrum.
    *       Notice that this must be allocated before calling this routine with enough length to be able to fill in all
    *       data points. Failure to allocate a large enough array will cause a buffer overrun and crash the program.
    *       This may be \c NULL if only the size of the last spectrum should be retrieved.
    *   @param[out] length will on successful return be filled with the length of the last acquired spectrum.
    *   @param[out] x_min will on successful return be filled with the smallest x-axis value of the spectrum
    *   @param[out] x_max will on successful return be filled with the largest x-axis value of the spectrum
    *   @return #FTS_SUCCESS if the last acquired spectrum was successfully retrieved
    *   @return #FTS_ERROR_NODATA if no spectrum has yet bee retrieved asynchronously.

    *   To retrieve the size of the last spectrum, set spec to \c NULL and retrieve the value
    *   of the parameter length upon successful return.
    *   The spectrum has the x-axis unit of wavenumbers and the y-axis unit of mW.    
    *   The data are sampled at equal intervals from x_min to x_max.
    */
    int __cdecl   FTS_GetLastSpectrum_Array(SpectrometerIndex specIndex, float* spec, unsigned int* length, float* x_min, float* x_max);
    int __stdcall FTS_GetLastSpectrum_Array_std(SpectrometerIndex specIndex, float* spec, unsigned int* length, float* x_min, float* x_max);

    /**  Retrieve a copy of the last spectrum acquired asynchronously, complete with x- and y-axis data.
    *   @param[in] specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers.
    *   @param[out] spec will on successful return be filled with the data in the last acquired spectrum.
    *       Notice that this must be allocated before calling this routine with enough length to be able to fill in all
    *       data points. Failure to allocate a large enough array will cause a buffer overrun and crash the program.
    *       This may be \c NULL if only the size of the last spectrum should be retrieved.
    *   @param[out] xAxisData will on successful return be filled with the x-axis data of the last acquired spectrum.
    *       If the spectrum has its x-axis data stored in X_VAL_ARRAY then this must be of equal length as spec.
    *       If the spectrum has its x-axis data stored in X_VAL_MINMAX then this must have a length of two, and this will only fill
    *       in the smallest and largest x-axis value of the spectrum (the intermediate values can then be calculated by linear interpolation).
    *       Notice that this must be allocated before calling this routine with enough length to be able to fill in all
    *       data points. Failure to allocate a large enough array will cause a buffer overrun and crash the program.
    *       This may be \c NULL if only the size of the last spectrum should be retrieved.
    *   @param[out] phaseData reserved for future use, set to \c NULL. 
    *   @param[out] xValueFormat will on successful return be filled with the x-value format of the last acquired spectrum.
    *       This is set to X_VAL_ARRAY if the x-axis data is stored as an array, or set to X_VAL_MINMAX if the x-axis data is stored
    *       with only the minimum and maximum values.
    *   @param[out] length will on successful return be filled with the length of the last acquired spectrum.
    *   @return #FTS_SUCCESS if the last acquired spectrum was successfully retrieved
    *   @return #FTS_ERROR_NODATA if no spectrum has yet bee retrieved asynchronously.
    *
    *   To retrieve the size and x-axis format of the last spectrum, set all pointers except for length and xValueFormat to \c NULL and retrieve the value
    *   of the parameters length and xValueFormat upon successful return.
    *   The spectrum has the x-axis unit of wavenumbers and the y-axis unit of mW. 
    */
    int __cdecl   FTS_GetLastSpectrumArrays(SpectrometerIndex specIndex, float* spec, float* xAxisData, float* phaseData, int* xValueFormat, unsigned int* length);
    int __stdcall FTS_GetLastSpectrumArrays_std(SpectrometerIndex specIndex, float* spec, float* xAxisData, float* phaseData, int* xValueFormat, unsigned int* length);

    /** This clears the contents of the last acquired spectrum for the specified spectrometer.
    *   This can be useful for e.g. restarting an averaging */
    int __cdecl   FTS_ClearLastSpectrum(SpectrometerIndex specIndex);
    int __stdcall FTS_ClearLastSpectrum_std(SpectrometerIndex specIndex);

    /** Retrieve a copy of the last interferogram acquired asynchronously.
    *   @param[in] specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers.
    *   @param[out] spec will on successful return be filled with the data in the last acquired interferogram.
    *   Notice that this must be allocated before calling this routine with a large enough buffer spectrum_t::I.
    *   This may be \c NULL if only the size of the last interferogram should be retrieved.
    *   @param[out] length will on successful return be filled with the length of the last acquired interferogram.
    *   
    *   @return #FTS_SUCCESS if the last acquired interferogram was successfully retrieved
    *   @return #FTS_ERROR_NODATA if no interferogram has yet been retrieved asynchronously.
    *   @return #FTS_ERROR_OUTOFMEMORY if the supplied spectrum_t spec was too small to fit in the data.
    *   
    *   To retrieve the size of the last interferogram, set spec to \c NULL and retrieve the value
    *   of the parameter length upon successful return.
    */
    int __cdecl   FTS_GetLastInterferogram(SpectrometerIndex specIndex, spectrum_t* spec, unsigned int* length);
    int __stdcall FTS_GetLastInterferogram_std(SpectrometerIndex specIndex, spectrum_t* spec, unsigned int* length);

    /** Retrieve a copy of the last interferogram acquired asynchronously.
    *   @param[in] specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers.
    *   @param[out] spec will on successful return be filled with the data in the last acquired interferogram.
    *   Notice that this must be allocated before calling this routine with enough length to be able to fill in all
    *   data points. Failure to allocate a large enough array will cause a buffer overrun and crash the program.
    *   This may be \c NULL if only the size of the last interferogram should be retrieved.
    *   @param[out] length will on successful return be filled with the length of the last acquired interferogram.
    *   @param[out] x_min will on successful return be filled with the smallest x-axis value of the interferogram
    *   @param[out] x_max will on successful return be filled with the largest x-axis value of the interferogram.
    *   @return #FTS_SUCCESS if the last acquired interferogram was successfully retrieved
    *   @return #FTS_ERROR_NODATA if no interferogram has yet bee retrieved asynchronously.
    *   
    *   To retrieve the size of the last interferogram, set spec to \c NULL and retrieve the value
    *   of the parameter length upon successful return.
    *   The interferogram is sampled at equal distances between x_min and x_max.
    */
    int __cdecl   FTS_GetLastInterferogram_Array(SpectrometerIndex specIndex, float* spec, unsigned int* length, float* x_min, float* x_max);
    int __stdcall FTS_GetLastInterferogram_Array_std(SpectrometerIndex specIndex, float* spec, unsigned int* length, float* x_min, float* x_max);

    /** This reads out a part of the last interferogram acquired asynchronously.
    *   @param[in] specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers.
    *   @param[out] spec will on successful return be filled with a part of the data from the last acquired interferogram.
    *   This must not be \c NULL and must be allocated to be able to hold at least 'length' values.
    *   @param startIndex the first index in the interferogram to retrieve
    *   @param length the number of values to copy out. 
    *   @param[out] valuesUpdated will on successful return be filled with the number of values copied from the last interferogram
    *   to the provided buffer 'spec'. This may be \c NULL if this information is not wanted.
    *   This value may be smaller than 'length' if the end of the interferogram has been reached.
    *   @return #FTS_SUCCESS if the requested part of the last acquired interferogram was successfully retrieved.
    *   @return #FTS_ERROR_NODATA if no interferogram has yet bee retrieved asynchronously.
    *   @return #FTS_ERROR_PARAMETER_ERROR if startIndex is larger than the length of the interferogram or no input parameter was passed.
    *   
    *   To retrieve the size of the last interferogram or the x_min and x_max values, call FTS_GetLastInterferogram_Array with spec set to \c NULL.
    */
    int __cdecl   FTS_GetLastInterferogram_Part(SpectrometerIndex specIndex, float* spec, unsigned int startIndex, unsigned int length, unsigned int* valuesUpdated);
    int __stdcall FTS_GetLastInterferogram_Part_std(SpectrometerIndex specIndex, float* spec, unsigned int startIndex, unsigned int length, unsigned int* valuesUpdated);

    /** This clears the contents of the last acquired interferogram for the specified spectrometer. */
    int __cdecl   FTS_ClearLastInterferogram(SpectrometerIndex specIndex);
    int __stdcall FTS_ClearLastInterferogram_std(SpectrometerIndex specIndex);

    /** Retrieves the last calculated wavenumber in wavelength meter mode
    *   This retrieves the value of acquisition_buffer_t::wavelengthMeter_waveNr in the acquisition_buffer_t associated
    *   with the specified instrument. */
    double __cdecl   FTS_GetLastWavelengthmeter_Wavenr(SpectrometerIndex specIndex);
    double __stdcall FTS_GetLastWavelengthmeter_Wavenr_std(SpectrometerIndex specIndex);

    /** Retrieves the estimated uncertainty in the last calculated wavenumber in wavelength meter mode
    *   This retrieves the value of acquisition_buffer_t::wavelengthMeterError_waveNr in the acquisition_buffer_t associated
    *   with the specified instrument. */
    double __cdecl   FTS_GetLastWavelengthmeterError_Wavenr(SpectrometerIndex specIndex);
    double __stdcall FTS_GetLastWavelengthmeterError_Wavenr_std(SpectrometerIndex specIndex);

    /** Retrieves the value of acquisition_buffer_t::interferogram_IsSaturated for the specified instrument. */
    bool __cdecl   FTS_GetLastInterferogramProperty_IsSaturated(SpectrometerIndex specIndex);
    bool __stdcall FTS_GetLastInterferogramProperty_IsSaturated_std(SpectrometerIndex specIndex);

    /** Retrieves the value of acquisition_buffer_t::interferogram_IsBroadBand for the specified instrument. */
    bool __cdecl   FTS_GetLastInterferogramProperty_IsBroadBand(SpectrometerIndex specIndex);
    bool __stdcall FTS_GetLastInterferogramProperty_IsBroadBand_std(SpectrometerIndex specIndex);

    /** Retrieves the value of acquisition_buffer_t::interferogram_IsMonochromatic for the specified instrument. */
    bool __cdecl   FTS_GetLastInterferogramProperty_IsMonochromatic(SpectrometerIndex specIndex);
    bool __stdcall FTS_GetLastInterferogramProperty_IsMonochromatic_std(SpectrometerIndex specIndex);

    /** Retrieves the value of acquisition_buffer_t::interferogram_MaxPeakToPeakSignal_counts for the specified instrument. */
    double __cdecl   FTS_GetLastInterferogramProperty_MaxPeakToPeakSignal(SpectrometerIndex specIndex);
    double __stdcall FTS_GetLastInterferogramProperty_MaxPeakToPeakSignal_std(SpectrometerIndex specIndex);

    /** Retrieves the value of acquisition_buffer_t::interferogram_MaxPeakToPeakPosition_cm for the specified instrument. */
    double __cdecl   FTS_GetLastInterferogramProperty_MaxPeakToPeakPosition(SpectrometerIndex specIndex);
    double __stdcall FTS_GetLastInterferogramProperty_MaxPeakToPeakPosition_std(SpectrometerIndex specIndex);
    ///@}


    /** \addtogroup grp_memory_management Memory Management
    *  @{
    */

    /** Retrieve the address of the global trace data structure */
    FTSTraceData * __cdecl   FTS_GetFTSTraceData();
    FTSTraceData * __stdcall FTS_GetFTSTraceData_std();

    /** Retrieve the address of a single spectrum_t in the global trace data structure. 
    *   @return the pointer to the spectrum_t structure, or \c NULL if the provided index is >= #MAX_TRACEBUFFER_NUM */
    spectrum_t* __cdecl     FTS_GetTraceSpectrum(unsigned int traceIndex);
    spectrum_t* __stdcall   FTS_GetTraceSpectrum_std(unsigned int traceIndex);

    /** This copies all the data in the spectrum from to the spectrum to.
    *   @param[out] to will on successful return contain the same data as the second spectrum_t.
    *   @param[in] from contains the data to be copied.
    *       NOTE: the spectrum to must be allocated before calling this function 
    *       and the allocated length of the arrays to->I and to->x must be large
    *       enough to hold all the values to be copied
    *   @return #FTS_SUCCESS on success
    *   @return #FTS_ERROR_PARAMETER_ERROR if any of the parameters are not valid
    *   @return #FTS_ERROR_OUTOFMEMORY if the allocated length of the spectrum 'to' is too small     */
    int __cdecl   FTS_CopySpectrum(spectrum_t* to, const spectrum_t* from);
    int __stdcall FTS_CopySpectrum_std(spectrum_t* to, const spectrum_t* from);

    /** This copies all the data in the spectrum one of the trace_buffer:s in #FTSTraceData.
    *   This will, if necessary, allocate the size of the destination trace to make sure the source trace fits.
    *   @param[in] destinationTraceIndex the index of the trace_buffer in #FTSTraceData, must be less than #MAX_TRACEBUFFER_NUM.
    *   @param[in] from contains the data to be copied.
    *   @return #FTS_SUCCESS on success
    *   @return #FTS_ERROR_PARAMETER_ERROR if any of the parameters are not valid */
    int __cdecl   FTS_CopySpectrumToTrace(unsigned int destinationTraceIndex, const spectrum_t* from);
    int __stdcall FTS_CopySpectrumToTrace_std(unsigned int destinationTraceIndex, const spectrum_t* from);

    /** This copies all the properties of the source spectrum to the destination spectrum but no spectral data.
    *   @param[out] to will on successful return contain the same properties as the second spectrum_t.
    *   @param[in] from contains the data to be copied.
    *       NOTE: this will not reallocate the destination spectrum, nor check if the sizes are same.
    *       The buffers to->I, to->x and to->phi will be untoched as well as their allocation info.
    *   @return #FTS_SUCCESS on success
    *   @return #FTS_ERROR_PARAMETER_ERROR if any of the parameters are not valid */
    int __cdecl FTS_CopySpectrumProperties(spectrum_t* to, const spectrum_t* from);
    int __stdcall FTS_CopySpectrumProperties_std(spectrum_t* to, const spectrum_t* from);

    /** This copies data from the provided 'from' buffer to the provided spectrum_t.
    *   @param to will on successful return be filled with the given data.
    *   @param from defines the data to copy to the destination spectrum_t.
    *   @param length the number of data points to copy from 'from' to 'to'.
    *   @param bufferNumber defines the buffer to copy the data to. #FTS_BUFFER_INDEX_I (zero) corresponds to the spectral data buffer 'I', 
    *       #FTS_BUFFER_INDEX_X (one) corresponds to the x-axis data buffer 'x', #FTS_BUFFER_INDEX_PHI (two) corresponds to the phase data buffer 'phi'.
    *   @return #FTS_SUCCESS on success
    *   @return #FTS_ERROR_PARAMETER_ERROR if any of the input parameters is invalid.
    *   @return #FTS_ERROR_OUTOFMEMORY if the destination buffer is not allocated or allocated with too short length. */
    int __cdecl FTS_CopyToSpectrumBuffer(spectrum_t* to, const float* const from, unsigned int length, int bufferNumber);
    int __stdcall FTS_CopyToSpectrumBuffer_std(spectrum_t* to, const float* const from, unsigned int length, int bufferNumber);

    /** This clears one of the trace_buffer:s in #FTSTraceData  
    *   @return #FTS_SUCCESS on success
    *   @return #FTS_ERROR_PARAMETER_ERROR if the bufferIndex is larger than or equal to #MAX_TRACEBUFFER_NUM */
    int __cdecl   FTS_Trace_Clear(unsigned int bufferIndex);
    int __stdcall FTS_Trace_Clear_std(unsigned int bufferIndex);

    /** This allocates a spectrum of the given size and stores it in on of the trace_buffer:s in FTSData
    *   @param bufferIndex The index of the trace_buffer in #FTSTraceData, must be less than #MAX_TRACEBUFFER_NUM.
    *   @param spectrumLength the desired length of the spectral array
    *   @param addX if true then the x-array will be allocated with the same length 
    *       as the spectral array. If false then the x-array will set to \c NULL. 
    *   @return #FTS_SUCCESS if the allocation succeeds. */
    int __cdecl   FTS_Trace_Allocate(unsigned int bufferIndex, unsigned int spectrumLength, bool addX);
    int __stdcall FTS_Trace_Allocate_std(unsigned int bufferIndex, unsigned int spectrumLength, bool addX);

    /** This resizes a spectrum found in one of the trace_buffer:s in #FTSTraceData. 
    *   The contents of the trace buffer arrays will not be kept after the resize.
    *   @param bufferIndex The index of the trace_buffer in #FTSTraceData, must be less than #MAX_TRACEBUFFER_NUM.
    *   @param spectrumLength the desired length of the spectral array
    *   @param addX if true then the x-array will be allocated with the same length 
    *       as the spectral array. If false then the x-array will set to \c NULL. 
    *   @param addPhase if true then the phi-array will be allocated with the same length 
    *       as the spectral array. If false then the phi-array will set to \c NULL. 
    */
    int __cdecl   FTS_Trace_Resize(unsigned int bufferIndex, unsigned int spectrumLength, bool addX, bool addPhase);
    int __stdcall FTS_Trace_Resize_std(unsigned int bufferIndex, unsigned int spectrumLength, bool addX, bool addPhase);

    /** This resizes the 'logData' of a spectrum found in one of the trace_buffer:s in #FTSTraceData. 
    *   The contents of the existing logData (if any) will not be kept after the resize.
    *   @param bufferIndex The index of the trace_buffer in #FTSTraceData, must be less than #MAX_TRACEBUFFER_NUM.
    *   @param logSize the desired size of the logData buffer, in bytes. If zero then the logData will be removed.
    */
    int __cdecl   FTS_Trace_ResizeLogData(unsigned int bufferIndex, unsigned int logSize);
    int __stdcall FTS_Trace_ResizeLogData_std(unsigned int bufferIndex, unsigned int logSize);

    /** Copies out data from a buffer in the provided spectrum_t. 
    *   @param[in] source the spectrum_t to copy the data from.
    *       This must be provided and must be allocated with enough size to be able to hold 'length' values of type 'float'.
    *   @param bufferNumber identifies which buffer should be the source. #FTS_BUFFER_INDEX_I (zero) corresponds to the spectral data buffer 'I', 
    *       #FTS_BUFFER_INDEX_X (one) corresponds to the x-axis data buffer 'x', #FTS_BUFFER_INDEX_PHI (two) corresponds to the phase data buffer 'phi'.
    *   @param[out] destination the array to where the data will be copied. This must be allocated with enough size to be able to copy length values.
    *   @param startIndex the index in the buffer where to start copying.
    *   @param length the number of values to copy out. If this is longer than the number of data available then fewer values will be copied.
    *   @param[in,out] valuesCopied will on successful return be filled with the number of values copied. 

    *   @return #FTS_SUCCESS on successful return.
    *   @return #FTS_ERROR_PARAMETER_ERROR if: the spectrum pointer is invalid,  the bufferNumber is illegal or the startIndex is too large.
    *   @return #FTS_ERROR_NODATA if the buffer with the specified number is not allocated (e.g. bufferNumber is 1 and the buffer 'x' is null)
    */
    int __cdecl FTS_GetSpectrumBuffer(const spectrum_t* source, int bufferNumber, float* destination, unsigned int startIndex, unsigned int length, unsigned int* valuesCopied);
    int __stdcall FTS_GetSpectrumBuffer_std(const spectrum_t* source, int bufferNumber, float* destination, unsigned int startIndex, unsigned int length, unsigned int* valuesCopied);

    /** @return a pointer to the given buffer in the provided spectrum_t.
    *   @return null if either input is invalid.
    *   @param source the spectrum_t to locate a pointer in.
    *   @param bufferNumber identifies which buffer should be the returned. #FTS_BUFFER_INDEX_I (zero) corresponds to the spectral data buffer 'I', 
    *       #FTS_BUFFER_INDEX_X (one) corresponds to the x-axis data buffer 'x', #FTS_BUFFER_INDEX_PHI (two) corresponds to the phase data buffer 'phi'.
    */
    float* __cdecl FTS_GetSpectrumPointer(spectrum_t* source, int bufferNumber);
    float* __stdcall FTS_GetSpectrumPointer_std(spectrum_t* source, int bufferNumber);

    /** Retrieves the value of a parameter in the given spectrum_t using its name, cast to a double.
    *   @param[in] source the spectrum data structure to retrieve parameters from
    *   @param[in] parameterName the name of the parameter to get the value of. This must be a null terminated string, 
    *       exactly matching the name in the spectrum_t, e.g. "gainLevel" or "referenceWavelength_nm".
    *   @param index the index of the parameter, this is only useful if the parameter is an array, e.g. gainLevel.
    *       If the parameter to retrieve is not an array, then the value of this parameter will be ignored.

    *   @return the value of the queried parameter, cast to a double.
    *   @return 0 if the parameter is not found or the parameter cannot be cast to a double (e.g. any of the strings, like "interferometerSerial").

    *   For example, to get all the parameters in the array "smoothParam" this routine needs to be called four times, once with index being 0, once index being 1
    *   once index being 2 and once index being 3. 
    */
    double __cdecl FTS_GetSpectrumParameter_Double(const spectrum_t* source, const char* parameterName, unsigned int index);
    double __stdcall FTS_GetSpectrumParameter_Double_std(const spectrum_t* source, const char* parameterName, unsigned int index);

    /** Sets the value of a parameter in the given spectrum_t using its name, cast to a double.
    *   @param[in] source the spectrum data structure to set parameters in
    *   @param[in] parameterName the name of the parameter to set the value of. This must be a null terminated string, 
    *       exactly matching the name in the spectrum_t, e.g. "gainLevel" or "referenceWavelength_nm".
    *   @param value the new value to set.
    *   @param index the index of the parameter, this is only useful if the parameter is an array, e.g. gainLevel.
    *       If the parameter to retrieve is not an array, then the value of this parameter will be ignored.
    *   @return #FTS_SUCCESS if the value was successfully set.
    *   @return #FTS_ERROR_PARAMETER_ERROR if the name was not found.
    */
    int __cdecl FTS_SetSpectrumParameter_Double(spectrum_t* source, const char* parameterName, double value, unsigned int index);
    int __stdcall FTS_SetSpectrumParameter_Double_std(spectrum_t* source, const char* parameterName, double value, unsigned int index);
   
    /** Retrieves the value of a parameter in the given spectrum_t using its name, cast to an unsigned int. 
    *   for a description of the parameters, see #FTS_GetSpectrumParameter_Double */
    unsigned int __cdecl FTS_GetSpectrumParameter_UInt(const spectrum_t* source, const char* parameterName, unsigned int index);
    unsigned int __stdcall FTS_GetSpectrumParameter_UInt_std(const spectrum_t* source, const char* parameterName, unsigned int index);

    /** Sets the value of a parameter in the given spectrum_t using its name, cast to an unsigned int. 
    *   for a description of the parameters, see #FTS_SetSpectrumParameter_Double
    *   @return #FTS_SUCCESS if the name and index could be successfully found. 
    *   @return #FTS_ERROR_PARAMETER_ERROR if the name or index was not found. */
    int __cdecl FTS_SetSpectrumParameter_UInt(spectrum_t* source, const char* parameterName, unsigned int value, unsigned int index);
    int __stdcall FTS_SetSpectrumParameter_UInt_std(spectrum_t* source, const char* parameterName, unsigned int value, unsigned int index);

    /** Retrieves the value of a string parameter in the given spectrum_t using its name.
    *   for a description of the parameters, see #FTS_GetSpectrumParameter_Double.
    *   If the passed in parameter name is not a string, then null is returned. */
    const char* __cdecl FTS_GetSpectrumParameter_Char(const spectrum_t* source, const char* parameterName, unsigned int index);
    const char*__stdcall FTS_GetSpectrumParameter_Char_std(const spectrum_t* source, const char* parameterName, unsigned int index);

    /** Retrieves the value of a string parameter in the given spectrum_t using its name.
    *   @param source the spectrum to read a string from.
    *   @param parameterName the name of the string to read. This must exactly match the name of the string in the spectrum_t structure (e.g. 'name' or 'source')
    *   @param parameterValue pointer to a char buffer which on successful return will be filled with the desired string.
    *       This buffer must be allocated before calling this function with enough space to hold the entire contents of the desired string.
    *       If the desired string cannot be found then this buffer will be untouched. 
    *   @param index ignored.
    *   @return #FTS_SUCCESS if the name could be successfully found.
    *   @return #FTS_ERROR_PARAMETER_ERROR if the name could not be found or is not a string. */
    int __cdecl     FTS_GetSpectrumParameter_CharValue(const spectrum_t* source, const char* parameterName, char* parameterValue, unsigned int index);
    int __stdcall   FTS_GetSpectrumParameter_CharValue_std(const spectrum_t* source, const char* parameterName, char* parameterValue, unsigned int index);

    /** Sets the value of a string parameter in the given spectrum_t using its name.
    *   @param[in] source the spectrum data structure to set the parameters in
    *   @param[in] parameterName the name of the parameter to get the value of. This must be a null terminated string, 
    *       exactly matching the name in the spectrum_t, e.g. "name" or "comment".
    *   @param[in] newParameterValue the new contents of the string. This must not be longer than the size of the destination char-array.
    *   @return the number of characters written. */
    int __cdecl     FTS_SetSpectrumParameter_Char(spectrum_t* source, const char* parameterName, const char* newParameterValue);
    int __stdcall   FTS_SetSpectrumParameter_Char_std(spectrum_t* source, const char* parameterName, const char* newParameterValue);

    /** Validates the given spectrum_t for common data format errors.
    *   @return #FTS_SUCCESS if the spectrum_t is correct.
    *   @return #FTS_ERROR_DATA_ERROR if the spectrum_t is incorrect.
    *   @param spectrum the spectrum_t to validate, this may be an interferogram or a spectrum.
    *   @param errorMessage if not set to null then this string buffer will be filled with an error message describing what is wrong with the
    *       spectrum_t. If the return value is #FTS_SUCCESS then this will not be touched. 
    *   @param errorMessageLength the length of the errorMessage buffer. */
    int __cdecl     FTS_ValidateSpectrum(const spectrum_t* const spectrum, char* errorMessage, unsigned int errorMessageLength);
    int __stdcall   FTS_ValidateSpectrum_std(const spectrum_t* const spectrum, char* errorMessage, unsigned int errorMessageLength);

    ///@}

    /** \addtogroup grp_error_handling Error Handling
    *  @{
    */

    /** This retrieves a copy of the trace log list
    *   @param[out] buffer the buffer to copy the data to
    *   @param[in] bufferLength the number of chars in 'buffer'. Specifies the maximum number of 
    *   characters to copy. Must be at least 20 characters! */
    int __cdecl   FTS_GetTraceLog(char* buffer, unsigned int bufferLength);
    int __stdcall FTS_GetTraceLog_std(char* buffer, unsigned int bufferLength);

    /** This starts writing the contents of the trace log to file. Every trace log message after this call is appended
    *   to file. The name of the file is 'TraceLog_ID.txt' where is the ID provided in this call. 
    *   This will set #FTS_GetFTSData().system_settings.log_traceLog_messages to true and set the ID of the file.
    *   The file is located at #FTS_GetFTSData().system_settings.operator_library_path.
    *   To stop writing trace log to file, set #FTS_GetFTSData().system_settings.log_traceLog_messages to false.
    *   @param ID the number of the file to write. */
    void __cdecl   FTS_StartWriteTraceLog(unsigned int ID);
    void __stdcall FTS_StartWriteTraceLog_std(unsigned int ID);
    ///@}

#ifdef __cplusplus
}
#endif 

#endif 
