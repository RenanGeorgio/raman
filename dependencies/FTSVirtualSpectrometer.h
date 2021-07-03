////////////////////////////////////////////////////////////////////////
// FTSVirtualSpectrometer.h
//	This is the the methods for simulating virtual
//	Thorlabs Fourier Transform Spectrometers
//	
//	Copyright (c) 2019 Thorlabs Sweden, All Rights Reserved
////////////////////////////////////////////////////////////////////////

#include "FTSData.h"

#ifndef _VIRTUALFTSPECTROMETER_H
#define _VIRTUALFTSPECTROMETER_H

#ifdef __cplusplus
extern "C" {
#endif 

    /** \addtogroup grp_virtual_spectrometer Virtual Spectrometer
    *   @{
    */

    typedef unsigned short VIRTUAL_SOURCE_TYPE;
#define VIRTUAL_SOURCE_MONOCHROMATIC 0x0000
#define VIRTUAL_SOURCE_BROADBAND 0x0001

#define VIRTUAL_SOURCE_MAX_PEAK_NUM 64

    /** virtual_source_parameter_t describes a simulated source for a virtual OSA. */
    typedef struct virtual_source_parameter_t
    {
        /** The number of peaks in the spectrum. Must be strictly smaller than #VIRTUAL_SOURCE_MAX_PEAK_NUM */
        unsigned int    peakNum;

        /** The center wavelength for each peak in the simulated source, in nanometers vacuum.
        *   The number of elements used here is determined by the value of peakNum. */
        double          centerWavelength_nm[VIRTUAL_SOURCE_MAX_PEAK_NUM];

        /** The Full Width at Half Maximum (FWHM) for each peak in the simulated source, in nanometers vacuum.
        *   The number of elements used here is determined by the value of peakNum. */
        double          fwhm_nm[VIRTUAL_SOURCE_MAX_PEAK_NUM];

        /** The amplitude for each peak in the simulated source, in milli watts.
        *   The number of elements used here is determined by the value of peakNum. */
        double          peakAmplitude[VIRTUAL_SOURCE_MAX_PEAK_NUM];
    }virtual_source_parameter_t;

    /** This sets up a virtual Thorlabs Fourier Transform Spectrometers of the given mode
    *   at the given spectrometer index. The type of source that the virtual spectrometer should simulate
    *   is setup by a subsequent call to #FTS_SetupVirtualSource.
    *   @return #FTS_SUCCESS on success */
    int __cdecl FTS_CreateVirtualSpectrometer(SpectrometerIndex specIndex, INSTRUMENT_MODEL model);
    int __stdcall FTS_CreateVirtualSpectrometer_std(SpectrometerIndex specIndex, INSTRUMENT_MODEL model);

    /** Sets up the source properties for a previously created virtual spectrometer.
    *   @param specIndex the index of the previously created virtual spectrometer.
    *   @param type the type of source. Should be either #VIRTUAL_SOURCE_MONOCHROMATIC or #VIRTUAL_SOURCE_BROADBAND.
    *   @param parameters a pointer to a virtual_source_parameter_t describing the type of the source.
    *   @return #FTS_ERROR_PARAMETER_ERROR if the spectrometer at the given index is not a virtual spectrometer.
    *   @return #FTS_SUCCESS on success */
    int __cdecl FTS_SetupVirtualSource(SpectrometerIndex specIndex, VIRTUAL_SOURCE_TYPE type, virtual_source_parameter_t *parameters);
    int __stdcall FTS_SetupVirtualSource_std(SpectrometerIndex specIndex, VIRTUAL_SOURCE_TYPE type, virtual_source_parameter_t *parameters);

    ///@}

#ifdef __cplusplus
}
#endif 

#endif