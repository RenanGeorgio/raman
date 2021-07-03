////////////////////////////////////////////////////////////////////////
// FTSSpectrumAnalysis.h
// These are the Optical Spectrum Analysis functions 
// used by the Thorlabs Fourier Transform Spectrometer DLL. 
// 
// Copyright (c) 2019 Thorlabs Sweden, All Rights Reserved
////////////////////////////////////////////////////////////////////////

#include "FTSData.h"

#ifndef FTSOSA_H
#define FTSOSA_H

#ifdef __cplusplus
extern "C" {
#endif 

    /** \addtogroup grp_osa_analysis OSA Analysis
    *  @{
    */

    // ----------------------------------------------------------------------------------------
    // ---------------------------- WAVELENGTH METER ANALYSIS --------------------------------
    // ----------------------------------------------------------------------------------------

    /** Estimates the wavelength contained in the given interferogram.
    *   Requires an interferogram from a monochromatic light source
    *   @param[in] interferogram the measured interferogram
    *   @param[in] startIdx the data point to start the calculation at (inclusive)
    *   @param[in] stopIdx the data point to stop the calculation at (exclusive)
    *   @param[out] waveNr will on successful return be filled with the calculated wavenumber (in 1/cm).
    *    if the calculation cannot be performed then this will be set to #FTS_WAVELENGTHMETER_ILLEGAL_VALUE
    *   @param[out] waveNrError will on successful return be filled with the estimated error in the calculated wavenumber
    *   @param[out] cycleNum if not NULL, then this will be set to the number of complete cycles counted
    *   @return #FTS_SUCCESS on success
    */
    int __cdecl   FTS_Wavelengthmeter(const spectrum_t* interferogram, unsigned int startIdx, unsigned int stopIdx, double* waveNr, double* waveNrError, int* cycleNum);
    int __stdcall FTS_Wavelengthmeter_std(const spectrum_t* interferogram, unsigned int startIdx, unsigned int stopIdx, double* waveNr, double* waveNrError, int* cycleNum);

    /** A wavelength_meter_result structure is used to handle the result from an Wavelength Metering analysis */ 
    typedef struct wavelength_meter_result
    {
        /** wavenumber_value is the calculated wavenumber of the interferogram, in cm^-1. 
        Set to #FTS_WAVELENGTHMETER_ILLEGAL_VALUE if the analysis has failed. */
        double wavenumber_value;

        /** wavenumber_error is the estimated error in the wavenumber_value, in cm^-1 */
        double wavenumber_error;

        /** cycle_num is the number of cycles counted in the interferogram */
        double cycle_num;

        /** errorMessage describes the reason that the wavelength metering analysis failed, if it did fail */
        char errorMessage[512];
    }wavelength_meter_result;

    /** Estimates the wavelength contained in the given interferogram.
    *   Requires an interferogram from a monochromatic light source
    *   @param[in] interferogram the measured interferogram
    *   @param[in] startIdx the data point to start the calculation at (inclusive)
    *   @param[in] stopIdx the data point to stop the calculation at (exclusive)
    *   @param[out] result will on successful return be filled with the calculated values of the analysis. 
    *   @return #FTS_SUCCESS on success
    *   @return #FTS_ERROR_PARAMETER_ERROR if 'result' is null.
    */
    int __cdecl   FTS_Wavelengthmeter_ext(const spectrum_t* interferogram, unsigned int startIdx, unsigned int stopIdx, wavelength_meter_result* result);
    int __stdcall FTS_Wavelengthmeter_ext_std(const spectrum_t* interferogram, unsigned int startIdx, unsigned int stopIdx, wavelength_meter_result* result);

    /** \deprecated
    *   Deprecated! Estimates the wavelength contained in the given interferogram. 
    *   This method is deprecated and should not be used. Use #FTS_Wavelengthmeter_ext instead!
    *   Requires an interferogram from a monochromatic light source
    *   @param[in] interferogramData the measured interferogram
    *   @param[in] sampleDistance_cm_vac The distance (in centimeters) between two sample points in the interferogram in vacuum
    *   @param[in] startIdx 
    *   @param[in] stopIdx defines the range of data points in the interferogram to calculate over
    *   @param[out] waveNr will on successful return be filled with the calculated wavenumber (in 1/cm).
    *    if the calculation cannot be performed then this will be set to #FTS_WAVELENGTHMETER_ILLEGAL_VALUE
    *   @param[out] waveNrError will on successful return be filled with the estimated error in the calculated wavenumber
    *   @param[out] cycleNum If not NULL, then this will be set to the number of complete cycles counted
    *   @return #FTS_SUCCESS on success
    }
    */
    int __cdecl   FTS_Wavelengthmeter_Array(const float* interferogramData, double sampleDistance_cm_vac, unsigned int startIdx, unsigned int stopIdx, double* waveNr, double* waveNrError, int* cycleNum);
    int __stdcall FTS_Wavelengthmeter_Array_std(const float* interferogramData, double sampleDistance_cm_vac, unsigned int startIdx, unsigned int stopIdx, double* waveNr, double* waveNrError, int* cycleNum);

    // ----------------------------------------------------------------------------------------
    // ---------------------------- COHERENCE LENGTH ANALYSIS --------------------------------
    // ----------------------------------------------------------------------------------------

    /** Estimates the coherence length of the input light source by analyzing the interferogram.
    *   @param[in] interferogramData the measured interferogram
    *   @param[in] sampleDistance_cm_vac The distance (in centimeters) between two sample points in the interferogram
    *   @param[in] startIdx, stopIdx Defines the range of data points in the interferogram to calculate over
    *   @param[out] coherenceLength_cm will on successful return be filled with the calculated coherence length (in cm).
    *   if the calculation cannot be performed then this will be set to #FTS_WAVELENGTHMETER_ILLEGAL_VALUE
    *   @param[out] coherenceLengthError_cm will on successful return be filled with the estimated error in the calculated
    *   coherence length, also this in cm.
    *   @return #FTS_SUCCESS on success
    */
    int __cdecl   FTS_CoherenceLength_Array(const float* interferogramData, double sampleDistance_cm_vac, unsigned int startIdx, unsigned int stopIdx, double* coherenceLength_cm, double* coherenceLengthError_cm);
    int __stdcall FTS_CoherenceLength_Array_std(const float* interferogramData, double sampleDistance_cm_vac, unsigned int startIdx, unsigned int stopIdx, double* coherenceLength_cm, double* coherenceLengthError_cm);

    // ----------------------------------------------------------------------------------------
    // --------------------------------- OSNR ANALYSIS ---------------------------------------
    // ----------------------------------------------------------------------------------------


    /** This performs an OSNR analysis of the given spectrum.
    *   The peaks to analyze must be supplied through the array 'peakCenter' (can e.g. be found by first calling FTS_FindPeaks)
    *   @param[in] spec The spectrum to process
    *   @param[in] peakCenter The position of each of the peaks to analyze
    *       this must be in the same unit as the spectrum x-axis and they must be sorted in order of increasing x-value
    *   @param[out] peakPower_dBm Will on successful return be filled with the power of each of the peaks, in dBm.
    *    If this is NULL then no peak power value will be returned
    *   @param[out] peakNoise_dBm Will on successful return be filled with the noise of each of the peaks, in dBm
    *    If this is NULL then no peak noise value will be returned
    *   @param[out] peakOSNR_dB Will on successful return be filled with the calculated OSNR of each of the peaks, in dB (this must not be NULL)
    *   @param[in] peakNum The number of peaks in 'peakCenter', also the number of values returned in each of the output arrays
    *   @param[in] Br Reference optical bandwidth, in nm
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_OSNR(const spectrum_t* spec, const float* peakCenter, float* peakPower_dBm, float* peakNoise_dBm, float* peakOSNR_dB, unsigned int peakNum, float Br);
    int __stdcall FTS_OSNR_std(const spectrum_t* spec, const float* peakCenter, float* peakPower_dBm, float* peakNoise_dBm, float* peakOSNR_dB, unsigned int peakNum, float Br);

    /** This performs an OSNR analysis of the given spectrum.
    *   The peaks to analyze must be supplied through the array 'peakCenter' 
    *   (can e.g. be found by first calling FTS_FindPeaks)
    *   @param[in] xData can be given in three different forms as specified in \ref sec_x_value_formats.
    *   @param[in] xValueOption determines the options for how the x-data is interpreted, see \ref sec_x_value_formats.
    *   @param[in] yData The spectral data to process
    *   @param[in] length The length of yData (and of xData if xValueOption equals #FTS_OPTION_ARRAY_X)
    *   @param[in] yAxisUnit The unit of the y-axis data, must be one of the Y_UNIT_... constants defined in FTSData.h
    *   @param[in] Bm The optical resolution of the spectrum, in nm
    *   @param[in] peakCenter The position of each of the peaks to analyze
    *       this must be in the same unit as the spectrum x-axis and they must be sorted in order of increasing x-value
    *   @param[out] peakPower_dBm Will on successful return be filled with the power of each of the peaks, in dBm.
    *       If this is NULL then no peak power value will be returned
    *   @param[out] peakNoise_dBm Will on successful return be filled with the noise of each of the peaks, in dBm
    *       If this is NULL then no peak noise value will be returned
    *   @param[out] peakOSNR_dB Will on successful return be filled with the calculated OSNR of each of the peaks, in dB (this must not be NULL)
    *   @param[in] peakNum The number of peaks in 'peakCenter', also the number of values returned in each of the output arrays
    *   @param[in] Br Reference optical bandwidth, in nm
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_OSNR_Array(const float* xData, const float* yData, unsigned int length, int yAxisUnit, unsigned int xValueOption, float Bm, const float* peakCenter, float* peakPower_dBm, float* peakNoise_dBm, float* peakOSNR_dB, unsigned int peakNum, float Br);
    int __stdcall FTS_OSNR_Array_std(const float* xData, const float* yData, unsigned int length, int yAxisUnit, unsigned int xValueOption, float Bm, const float* peakCenter, float* peakPower_dBm, float* peakNoise_dBm, float* peakOSNR_dB, unsigned int peakNum, float Br);

    /** These functions are same as above but also returns the points where the noise were calculated and the calculated noise there 
    *   @param[in] spec The spectrum to analyze
    *   @param[in] peakCenter The position of each of the peaks to analyze
    *       this must be in the same unit as the spectrum x-axis and they must be sorted in order of increasing x-value
    *   @param[out] peakPower_dBm Will on successful return be filled with the power of each of the peaks, in dBm.
    *       If this is NULL then no peak power value will be returned
    *   @param[out] peakNoise_dBm Will on successful return be filled with the noise of each of the peaks, in dBm
    *       If this is NULL then no peak noise value will be returned
    *   @param[out] peakOSNR_dB Will on successful return be filled with the calculated OSNR of each of the peaks, in dB (this must not be NULL)
    *   @param peakNum The number of peaks in 'peakCenter', also the number of values returned in each of the output arrays
    *   @param[in] Br Reference optical bandwidth, in nm
    *   @param[out] noisePt Will on successful return be filled with the positions where the noise was calculated,
    *       must be pre-allocated to hold at least 2*peakNum data points (two points for each peak)
    *   @param[out] noiseLevel Will on successful return be filled with the calculated noise
    *       must be pre-allocated to hold at least 2*peakNum data points (two points for each peak)
    */
    int __cdecl   FTS_OSNR_Ext(const spectrum_t* spec, const float* peakCenter, float* peakPower_dBm, float* peakNoise_dBm, float* peakOSNR_dB, unsigned int peakNum, float Br, float* noisePt, float* noiseLevel);
    int __stdcall FTS_OSNR_Ext_std(const spectrum_t* spec, const float* peakCenter, float* peakPower_dBm, float* peakNoise_dBm, float* peakOSNR_dB, unsigned int peakNum, float Br, float* noisePt, float* noiseLevel);

    /** This performs an OSNR analysis of the given spectrum. These performs the same operations as FTS_OSNR() but also returns the 
    *   points where the noise were calculated and the calculated noise there.
    *   The peaks to analyze must be supplied through the array 'peakCenter' 
    *   (can e.g. be found by first calling FTS_FindPeaks())
    *   @param[in] xData can be given in three different forms as specified in \ref sec_x_value_formats.
    *   @param[in] xValueOption determines the options for how the x-data is interpreted, see \ref sec_x_value_formats.
    *   @param[in] yData The spectral data to process
    *   @param[in] length The length of yData (and of xData if xValueOption equals #FTS_OPTION_ARRAY_X)
    *   @param[in] yAxisUnit The unit of the y-axis data, must be one of the Y_UNIT_... constants defined in FTSData.h
    *   @param[in] Bm The optical resolution of the spectrum, in nm
    *   @param[in] peakCenter The position of each of the peaks to analyze
    *       this must be in the same unit as the spectrum x-axis and they must be sorted in order of increasing x-value
    *   @param[out] peakPower_dBm Will on successful return be filled with the power of each of the peaks, in dBm.
    *       If this is NULL then no peak power value will be returned
    *   @param[out] peakNoise_dBm Will on successful return be filled with the noise of each of the peaks, in dBm
    *       If this is NULL then no peak noise value will be returned
    *   @param[out] peakOSNR_dB Will on successful return be filled with the calculated OSNR of each of the peaks, in dB (this must not be NULL)
    *   @param[out] noisePt Will on successful return be filled with the positions where the noise was calculated,
    *       must be pre-allocated to hold at least 2*peakNum data points (two points for each peak)
    *   @param[out] noiseLevel Will on successful return be filled with the calculated noise
    *       must be pre-allocated to hold at least 2*peakNum data points (two points for each peak)
    *   @param[in] peakNum The number of peaks in 'peakCenter', also the number of values returned in each of the output arrays
    *   @param[in] Br Reference optical bandwidth, in nm
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_OSNR_Ext_Array(const float* xData, const float* yData, unsigned int length, int yAxisUnit, unsigned int xValueOption, float Bm, const float* peakCenter, float* peakPower_dBm, float* peakNoise_dBm, float* peakOSNR_dB, unsigned int peakNum, float Br, float* noisePt, float* noiseLevel);
    int __stdcall FTS_OSNR_Ext_Array_std(const float* xData, const float* yData, unsigned int length, int yAxisUnit, unsigned int xValueOption, float Bm, const float* peakCenter, float* peakPower_dBm, float* peakNoise_dBm, float* peakOSNR_dB, unsigned int peakNum, float Br, float* noisePt, float* noiseLevel);

    // ----------------------------------------------------------------------------------------
    // --------------------------------- OPTICAL POWER ---------------------------------------
    // ----------------------------------------------------------------------------------------

    /** This calculates the total optical power in the given spectrum between the 
    *   points xMin and xMax 
    *   @param[in] spectrum The spectrum to calculate the power of
    *   @param[in] xMin The lower edge of the wavelength/wavenumber/frequency/... region
    *       of which to calculate the power. Must be in the same unit as the
    *       x-axis unit of the spectrum.
    *   @param[in] xMax The upper edge of the wavelength/wavenumber/frequency/... region
    *       of which to calculate the power. Must be in the same unit as the
    *       x-axis unit of the spectrum.
    *   @param[out] result Will on successful return be filled with the calculated power in the same unit as the y-axis unit of the spectrum
    *   @return #FTS_ERROR_PARAMETER_ERROR if the xMin to xMax range does not overlap the spectrum.x_min to spectrum.x_max range
    *   @return #FTS_ERROR_PARAMETER_ERROR if the y - axis unit of the spectrum is not a power calibrated unit (mW, mW_norm, dBm or dBm_norm).
    *   @return #FTS_SUCCESS on successful calculation of the power
    */
    int __cdecl   FTS_OpticalPower(const spectrum_t* spectrum, float xMin, float xMax, double* result);
    int __stdcall FTS_OpticalPower_std(const spectrum_t* spectrum, float xMin, float xMax, double* result);

    /** This calculates the total optical power in the given spectrum between the points xMin and xMax 
    *   @param[in] xData can be given in three different forms as specified in \ref sec_x_value_formats.
    *   @param[in] xValueOption determines the options for how the x-data is interpreted, see \ref sec_x_value_formats.
    *   @param[in] yData The y-value for each data point 
    *   @param[in] length The length of yData (and of xData if xValueOption equals #FTS_OPTION_ARRAY_X)
    *   @param[in] xAxisUnit The unit of the x-axis data, must be one of the X_UNIT_... constants defined in FTSData.h
    *   @param[in] yAxisUnit The unit of the y-axis data, must be one of the Y_UNIT_... constants defined in FTSData.h
    *   @param[in] resolution_wnr The optical resolution of the spectrum, must be given in wavenumbers
    *   @param[in] xMin The lower edge of the wavelength/wavenumber/frequency/... region
    *       of which to calculate the power. Must be in the same unit as the
    *       x-axis unit of the spectrum.
    *   @param[in] xMax The upper edge of the wavelength/wavenumber/frequency/... region
    *       of which to calculate the power. Must be in the same unit as the
    *       x-axis unit of the spectrum.
    *   @param[out] result Will on successful return be filled with the calculated power, in the same unit as the y-axis unit of the spectrum
    *   @return #FTS_ERROR_PARAMETER_ERROR if the xMin to xMax range does not overlap the spectrum.x_min to spectrum.x_max range
    *   @return #FTS_ERROR_PARAMETER_ERROR if the y - axis unit of the spectrum is not a power calibrated unit (mW, mW_norm, dBm or dBm_norm).
    *   @return #FTS_SUCCESS on successful calculation of the power
    */
    int __cdecl   FTS_OpticalPower_Array(const float* xData, const float* yData, unsigned int length, int xAxisUnit, int yAxisUnit, unsigned int xValueOption, float resolution_wnr, float xMin, float xMax, double* result);
    int __stdcall FTS_OpticalPower_Array_std(const float* xData, const float* yData, unsigned int length, int xAxisUnit, int yAxisUnit, unsigned int xValueOption, float resolution_wnr, float xMin, float xMax, double* result);

    typedef struct power_options_t
    {
        /** TxMin The lower edge of the wavelength/wavenumber/frequency/... region
        *       of which to calculate the power. Must be in the same unit as the
        *       x-axis unit of the spectrum. */
        double xMin;

        /** The upper edge of the wavelength/wavenumber/frequency/... region
        *       of which to calculate the power. Must be in the same unit as the
        *       x-axis unit of the spectrum.*/
        double xMax;
    }power_options_t;

    typedef struct power_result_t
    {
        /** The calculated power value */
        double powerValue;

        /** The unit of the calculated power value */
        Y_AXIS_UNIT powerUnit;

        bool isResolutionLimited;
    }power_result_t;

    int __cdecl FTS_OpticalPowerAnalysis(const spectrum_t* const spectrum, const power_options_t* const options, power_result_t* result);
    int __stdcall FTS_OpticalPowerAnalysis_std(const spectrum_t* const spectrum, const power_options_t* const options, power_result_t* result);

    // ----------------------------------------------------------------------------------------
    // --------------------------------- COLOR ANALYSIS ---------------------------------------
    // ----------------------------------------------------------------------------------------

    enum COLOR_SPACE
    {
        COLOR_SPACE_1932_2DEG,
        COLOR_SPACE_1960,
        COLOR_SPACE_1976
    };

    /** This retrieves the chromaticity values for the given spectrum (using 1931 CIE CMF, 2 degrees)
    *   @param[in] spectrum Spectrum to calculate the values for, this must have the x-axis unit of nanometers!
    *   @param[in] threshold All spectral data points with intensity below this value will be ignored,
    *       to reduce the effects of noise
    *   @param[out] x Will on successful return be set to the weight for the X color matching function
    *   @param[out] y Will on successful return be set to the weight for the Y color matching function
    *   @param[out] z Will on successful return be set to the weight for the Z color matching function
    *   @param[out] mainWavelength Will on successful return be set to the main Wavelength in the spectrum
    *   @param[out] colorPurity Will on successful return be set to the calculated color purity, in percent.
    *   @param[out] correlatedColorTemperature Will on successful return be set to the calculated correlated color temperature (CCT) in Kelvin.
    *   @return #FTS_SUCCESS if successful
    *   @return #FTS_ERROR_PARAMETER_ERROR if the input wavelength does not cover the entire range 380 to 700 nm
    *   @return #FTS_ERROR_PARAMETER_ERROR if the unit of the input spectrum is not a calibrated unit (Y_UNIT_DBM, Y_UNIT_MW, Y_UNIT_DBM_NORM, Y_UNIT_MW_NORM, Y_UNIT_CAL_I or Y_UNIT_CAL_LOGI)
    *   @return #FTS_ERROR_NODATA if no data point between 380 and 700 nm is above the threshold
    */
    int __cdecl   FTS_GetColor(const spectrum_t* spectrum, float threshold, float* x, float* y, float* z, float* mainWavelength, float* colorPurity, float* correlatedColorTemperature);
    int __stdcall FTS_GetColor_std(const spectrum_t* spectrum, float threshold, float* x, float* y, float* z, float* mainWavelength, float* colorPurity, float* correlatedColorTemperature);

    typedef struct color_information 
    {
        /** CIE1932_2Deg_x, CIE1932_2Deg_y and CIE1932_2Deg_z are the calculated coordinates
        *   int the CIE 1932 color space using the standard, 2 degrees, observer */
        double CIE1932_2Deg_x;
        double CIE1932_2Deg_y;
        double CIE1932_2Deg_z;

        /** CIE1932_10Deg_x, CIE1932_10Deg_y and CIE1932_10Deg_z are the calculated coordinates
        *   int the CIE 1932 color space using the supplementary, 10 degrees, observer */
        double CIE1932_10Deg_x;
        double CIE1932_10Deg_y;
        double CIE1932_10Deg_z;

        double mainWavelength;

        double colorPurity;

        double correlatedColorTemperature;

        /** returnCode is the return code from the last call to FTS_ColorAnalysis which filled in this struct */
        int returnCode;

        /** If the returnCode is not #FTS_SUCCESS then an error message can be filled in here */
        char userMessage[512];
    }color_information;

    typedef struct color_analysis_parameters
    {
        /** All spectral data points with intensity below the threshold value will be ignored, to reduce the effects of noise */
        float threshold;
    }color_analysis_parameters;

    /** This retrieves the chromaticity values for the given spectrum (using 1931 CIE CMF, 2 degrees and 10 degrees)
    *   @param[in] spectrum Spectrum to calculate the values for, this must have the x-axis unit of nanometers!
    *   @param[in] parameters Defines the options for the analysis.
    *   @param[out] result Will on successful return be filled with the result of the color analysis calculations. 
    *   @return #FTS_SUCCESS if successful
    *   @return #FTS_ERROR_PARAMETER_ERROR if the input wavelength does not cover the entire range 380 to 700 nm
    *   @return #FTS_ERROR_PARAMETER_ERROR if the unit of the input spectrum is not a calibrated unit (Y_UNIT_DBM, Y_UNIT_MW, Y_UNIT_DBM_NORM, Y_UNIT_MW_NORM, Y_UNIT_CAL_I or Y_UNIT_CAL_LOGI)
    *   @return #FTS_ERROR_NODATA if no data point between 380 and 700 nm is above the threshold */
    int __cdecl FTS_ColorAnalysis(const spectrum_t* const spectrum, const color_analysis_parameters* const parameters, color_information* const result);
    int __stdcall FTS_ColorAnalysis_std(const spectrum_t* const spectrum, const color_analysis_parameters* const parameters, color_information* const result);

    /** This retrieves the chromaticity values for the given spectrum (using 1931 CIE CMF, 2 degrees)
    *   @param[in] lambda The wavelength (must be in nanometers) for each data point in 'spectralData'
    *   @param[in] spectralData The value of the spectrum (must be in milli Watts or dBm)
    *   @param[in] dataLength The number of data points in vectors 'lambda' and 'spectralData'
    *   @param[in] xUnit The unit of the 'lambda' array, must be either X_UNIT_NM_VAC or X_UNIT_NM_AIR
    *   @param[in] yUnit The unit of the spectral data, must be either mW or dBm
    *   @param[in] resolution_wnr The resolution of the spectral data, in wavenumbers
    *   @param[in] threshold All spectral data points with intensity below this value will be ignored,
    *       to reduce the effects of noise
    *   @param[out] x Will on successful return be set to the weight for the X color matching function
    *   @param[out] y Will on successful return be set to the weight for the Y color matching function
    *   @param[out] z Will on successful return be set to the weight for the Z color matching function
    *   @param[out] mainWavelength Will on successful return be set to the main Wavelength in the spectrum
    *   @param[out] colorPurity Will on successful return be set to the calculated color purity, in percent.
    *   @param[out] correlatedColorTemperature Will on successful return be set to the calculated correlated color temperature (CCT) in Kelvin.
    *   @return #FTS_SUCCESS if successful
    *   @return #FTS_ERROR_PARAMETER_ERROR if the input wavelength does not cover the entire range 380 to 700 nm
    *   @return #FTS_ERROR_NODATA if no data point between 380 and 700 nm is above the threshold */
    int __cdecl   FTS_GetColor_Array(const float* lambda, const float* spectralData, unsigned int dataLength, X_AXIS_UNIT xUnit, Y_AXIS_UNIT yUnit, float resolution_wnr, float threshold, float* x, float* y, float* z, float* mainWavelength, float* colorPurity, float* correlatedColorTemperature);
    int __stdcall FTS_GetColor_Array_std(const float* lambda, const float* spectralData, unsigned int dataLength, X_AXIS_UNIT xUnit, Y_AXIS_UNIT yUnit, float resolution_wnr, float threshold, float* x, float* y, float* z, float* mainWavelength, float* colorPurity, float* correlatedColorTemperature);

    /** Retrieves the 1932 (2deg) x and y coordinates from the given temperature in Kelvin.
    *   @param[in] temp_Kelvin the temperature in Kelvin, must be between 1000 and 25000 Kelvin.
    *   @param[out] x Will on successful return be set to the x coordinate.
    *   @param[out] y Will on successful return be set to the y coordinate.
    *   @return #FTS_SUCCESS if successful. */
    int __cdecl FTS_ColorAnalysis_Get1932XYFromTemperature(double temp_Kelvin, double* x, double* y);
    int __stdcall FTS_ColorAnalysis_Get1932XYFromTemperature_std(double temp_Kelvin, double* x, double* y);
    
    /** Retrieves the correlated color temperature from a given 1932 (2 deg) x and y coordinate.
    *   @param[in] x The CIE 1932 (2deg) x coordinate.
    *   @param[in] y The CIE 1932 (2deg) y coordinate.
    *   @return The calculated correlated color temperature (CCT) if successful. 
    *   @return -1 if the calculation was not successful. */
    double __cdecl FTS_ColorAnalysis_GetTemperatureFrom1932XY(double x, double y);
    double __stdcall FTS_ColorAnalysis_GetTemperatureFrom1932XY_std(double x, double y);

    /** Converts the two given coordinates from one color space into coordinates in another color space.
    *   @param[in] input1 the first input coordinate.
    *   @param[in] input2 the second input coordinate.
    *   @param[in] inputColorSpace defines how to interpret the input1 and input2 coordinates.
    *       If inputColorSpace is COLOR_SPACE_1932_2DEG then input1 is the x coordinate and input2 is the y coordinate.
    *       If inputColorSpace is COLOR_SPACE_1960 then input1 is the u coordinate and input2 is the v coordinate.
    *       If inputColorSpace is COLOR_SPACE_1976 then input1 is the u' coordinate and input2 is the v' coordinate.
    *   @param[out] output1 Will on successful return be set to the 1960 u coordinate.
    *   @param[out] output2 Will on successful return be set to the 1960 v coordinate.
    *   @param[in] outputColorSpace defines how to interpret the output1 and output2 coordinates
    *       If outputColorSpace is COLOR_SPACE_1932_2DEG then output1 is the x coordinate and output2 is the y coordinate.
    *       If outputColorSpace is COLOR_SPACE_1960 then output1 is the u coordinate and output2 is the v coordinate.
    *       If outputColorSpace is COLOR_SPACE_1976 then output1 is the u' coordinate and output2 is the v' coordinate.
    *   @return #FTS_SUCCESS if successful. */
    int __cdecl   FTS_ColorAnalysis_ConvertColorSpaceCoordinates(float input1, float input2, COLOR_SPACE inputColorSpace, float* output1, float* output2, COLOR_SPACE outputColorSpace);
    int __stdcall FTS_ColorAnalysis_ConvertColorSpaceCoordinates_std(float input1, float input2, COLOR_SPACE inputColorSpace, float* output1, float* output2, COLOR_SPACE outputColorSpace);

    // ----------------------------------------------------------------------------------------
    // ---------------------------------- PEAK FINDING ---------------------------------------
    // ----------------------------------------------------------------------------------------

    /** A peak_t structure is used to handle information on peaks in 
    spectra / interferograms in an efficient way. */
    typedef struct peak_t
    {
        /** The peakPosition is the x-axis location where the highest value is found */
        double  peakPosition;

        /** This is the index where the highest value is found */ 
        unsigned int  peakPositionIndex;

        /** maximum value of the peak, in yAxisUnit:s */
        double  peakLevel;

        /** The center of gravity of the peak, in xAxisUnits */
        double  centroidPosition;

        /** The level at the center of gravity of the peak, in yAxisUnit:s */
        double  centroidLevel;

        /** baseLine is the level at which the base of the peak is found */
        double  baseLine;

        /** width is the width at the specified 'min_peak_height_dB' */
        double  width;

        /** fwhm is the width at 3dB from the peak*/
        double  fwhm;

        /** area is the area below the peak, either from the baseline or from the y=0 depending on the settings */
        double  area;

        /** leftPt is the nearest smaller x-axis position where the intensity has dropped by 'min_peak_height_dB' from the peak level */
        double  leftPt;

        /** rightPt is the nearest larger x-axis position where the intensity has dropped by 'min_peak_height_dB' from the peak level */
        double  rightPt;

        /** leftPt_3dB is the nearest smaller x-axis position where the intensity has dropped by 3dB from the peak level */
        double  leftPt_3dB;

        /** rightPt_3dB is the nearest larger x-axis position where the intensity has dropped by 3dB from the peak level */
        double  rightPt_3dB;

        /** baseLine_rightPt - baseLine_leftPt defines the width of the peak at the baseline level. */
        double  baseLine_leftPt;
        double  baseLine_rightPt;
    }peak_t;


    /** A valley_t structure is used to handle information on valleys in spectra / interferograms in
    an efficient way. */
    typedef struct valley_t
    {
        /** The minimumPosition is the x-axis location where the lowest value is found */
        double minimumPosition;

        /** value at the lowest point, in the yAxisUnit of the input spectrum */
        double minimumLevel;

        /** The center of gravity of the valley, in xAxisUnits */
        double centroidPosition;

        /** The level at the center of gravity of the valley, in yAxisUnits */
        double centroidLevel;

        /** baseLine is the level at which the base of the valley is found */
        double baseLine;

        /** baseLine_leftPt is the smaller x-axis location where the base of the valley is found. */
        double baseLine_leftPt;

        /** baseLine_leftLevel is the value of the input spectrum at the baseLine_leftPt, in xAxisUnits. */
        double baseLine_leftLevel;

        /** baseLine_rightPt is the larger x-axis location where the base of the valley is found. */
        double baseLine_rightPt;

        /** baseLine_rightLevel is the value of the input spectrum at the baseLine_rightPt, in xAxisUnits. */
        double baseLine_rightLevel;

        /** distance from surrounding area to the lowest point */
        double depth;

        /** distance from surrounding area to the lowest point, in dB */
        double depth_dB;

        /** halfWidth is the width of the valley at one-half of the total depth.
        This is equal to the difference between 'rightPt_3dB' and 'leftPt_3dB'. */
        double halfWidth;

        /** width is the width of the valley at 'min_peak_height_dB' above the lowest point.
        This is equal to the difference between 'rightPt' and 'leftPt'. */
        double width;

        /** leftPt is the smaller x-axis location where the width of the valley is 
        'min_peak_height_dB' above the lowest point. */
        double leftPt;

        /** rightPt is the larger x-axis location where the width of the valley is 
        'min_peak_height_dB' above the lowest point. */
        double rightPt;

        /** leftPt_3dB is the nearest smaller x-axis position where the 
        intensity has increased by 3dB from the lowest point. */
        double leftPt_3dB;

        /** leftPt_3dB is the nearest larger x-axis position where the 
        intensity has increased by 3dB from the lowest point. */
        double rightPt_3dB;

    }valley_t;

    /** Creates and returns an array of 'peak_t' structures, which may be passed as argument to FTS_FindPeaks_ext. 
    *   This array must be deleted using FTS_DeleteArray. 
    *   @return the allocated array.
    *   @return null if number is zero. */
    peak_t * __cdecl   FTS_CreateArray_PeakData(unsigned int number);
    peak_t * __stdcall FTS_CreateArray_PeakData_std(unsigned int number);

    /** Creates and returns an array of 'valley_t' structures, which may be passed as argument to FTS_FindValleys_ext. 
    *   This array must be deleted using FTS_DeleteArray. 
    *   @return the allocated array.
    *   @return null if number is zero. */
    valley_t * __cdecl   FTS_CreateArray_ValleyData(unsigned int number);
    valley_t * __stdcall FTS_CreateArray_ValleyData_std(unsigned int number);

    /** Frees up the memory allocated by a previous call to 'FTS_CreateArray_PeakData' or 'FTS_CreateArray_ValleyData'. */
    void __cdecl FTS_DeleteArray(void* peakOrValleyArray);
    void __stdcall FTS_DeleteArray_std(void* peakOrValleyArray);

    /** @return the value of one property of the array of peak_t data.
    *   @param[in] peakArray the array of peak_t structures to retrieve parameters from.
    *   @param indexOfPeak the index of the peak in the array (indices starting at zero)
    *   @param[in] parameterName the name of the parameter to get the value of. This must be a null terminated string, 
    *    exactly matching the name in the peak_t, e.g. "peakPosition" or "fwhm". */
    double __cdecl  FTS_GetPeakProperty(const peak_t* peakArray, unsigned int indexOfPeak, const char* parameterName);
    double __stdcall FTS_GetPeakProperty_std(const peak_t* peakArray, unsigned int indexOfPeak, const char* parameterName);

    /** @return the value of one property of the array of valley_t data.
    *   @param[in] valleyArray the array of valley_t structures to retrieve parameters from.
    *   @param indexOfValley the index of the valley in the array (indices starting at zero)
    *   @param[in] parameterName the name of the parameter to get the value of. This must be a null terminated string, 
    *       exactly matching the name in the valley_t, e.g. "centroidPosition" or "depth". */
    double __cdecl  FTS_GetValleyProperty(const valley_t* valleyArray, unsigned int indexOfValley, const char* parameterName);
    double __stdcall FTS_GetValleyProperty_std(const valley_t* valleyArray, unsigned int indexOfValley, const char* parameterName);

    /** Retrieves all peaks in a collected spectrum/interferogram between the given indices. 
    *   @param[in] spec The spectrum in which to search for peaks
    *   @param[in] startIdx The first index in the data set to use (inclusive)
    *   @param[in] stopIdx The last index in the data set to use (exclusive)
    *   @param[in] threshold Only peaks with intensity above this level will be considered (must be same unit as the y-data)
    *   @param[in] minPeakHeight_dB Only peaks with a height of at least this height over the surrounding baseline will be returned
    *   @param[in] sortOption Determines how the peaks will be sorted
    *       sortOption = 0 -> no sorting
    *       sortOption = 1 -> peaks are sorted in order of increasing center position
    *       sortOption = 2 -> peaks are sorted in order of increasing peak height
    *       sortOption = 3 -> peaks are sorted in order of increasing peak width
    *   @param[out] peakCenter Will on successful return be filled with the center position of each peak found, in the current xAxisUnit of the spectrum.
    *   @param[out] peakWidth Will on successful return be filled with the width of each peak found, in the current xAxisUnit of the spectrum.
    *   @param[out] peakHeight Will on successful return be filled with the height of each peak found, in the current xAxisUnit of the spectrum.
    *   @param[out] peakLeftPt If not NULL, this will on successful return be filled with the first position to the left of the peak
    *       where the intensity has dropped by 'minPeakHeight' dB from the maximum value, in the current xAxisUnit of the spectrum.
    *       This parameter may be NULL.
    *   @param peakRightPt If not NULL, this will on successful return be filled with the first position to the right of the peak
    *       where the intensity has dropped by 'minPeakHeight' dB from the maximum value, in the current xAxisUnit of the spectrum.
    *       This parameter may be NULL.
    *   @param[in] maxPeaks The maximum number of peaks that can be filled into the arrays peakCenter, peakWidth and peakHeight.
    *   @return The number of peaks found. -1 if an error occurs
    */
    int __cdecl   FTS_FindPeaks(spectrum_t* spec, unsigned int startIdx, unsigned int stopIdx, float threshold, float minPeakHeight_dB, unsigned int sortOption, float* peakCenter, float* peakWidth, float* peakHeight, float* peakLeftPt, float* peakRightPt, unsigned int maxPeaks);
    int __stdcall FTS_FindPeaks_std(spectrum_t* spec, unsigned int startIdx, unsigned int stopIdx, float threshold, float minPeakHeight_dB, unsigned int sortOption, float* peakCenter, float* peakWidth, float* peakHeight, float* peakLeftPt, float* peakRightPt, unsigned int maxPeaks);

    /** Retrieves all peaks in a collected spectrum/interferogram between the given indices. 
    *   If you have trouble allocating/deleting or accessing data in struct:s then have a look at the routines #FTS_CreateArray_PeakData, 
    *   #FTS_GetPeakProperty and #FTS_DeleteArray.
    *   @param[in] xData The x-value for each data point, can be given in three different forms as specified in \ref sec_x_value_formats.
    *   @param[in] yData The y-value for each data point 
    *   @param[in] length The length of the yArray (and also the xData array if xValueOption equals #FTS_OPTION_ARRAY_X)
    *   @param[in] xValueOption Determines the options for how the x-data is interpreted, see \ref sec_x_value_formats.
    *   @param[in] yAxisUnit The unit of the y-axis data, must be one of the Y_UNIT_... constants defined in FTSData.h
    *   @param[in] startIdx The first index in the data set to use (inclusive).
    *   @param[in] stopIdx The last index in the data set to use (exclusive).
    *   @param[in] threshold Only peaks with intensity above this level will be considered (must be same unit as the y-data).
    *   @param[in] minPeakHeight_dB Only peaks with a height of at least this height over the surrounding baseline will be returned.
    *   @param[in] sortOption Determines how the peaks will be sorted
    *       sortOption = 0 -> no sorting
    *       sortOption = 1 -> peaks are sorted in order of increasing center position
    *       sortOption = 2 -> peaks are sorted in order of increasing peak height
    *       sortOption = 3 -> peaks are sorted in order of increasing peak width
    *   @param[out] peakArray Will on successful return be filled with info on each peak found.
    *   @param[in] maxPeaks The maximum number of peaks that can be filled into the array peakArray.
    *   @return The number of peaks found. -1 if an error occurs
    */
    int __cdecl   FTS_FindPeaks_ext(const float* xData, const float* yData, unsigned int length, int yAxisUnit, unsigned int xValueOption, unsigned int startIdx, unsigned int stopIdx, float threshold, float minPeakHeight_dB, unsigned int sortOption, peak_t *peakArray, unsigned int maxPeaks);
    int __stdcall FTS_FindPeaks_ext_std(const float* xData, const float* yData, unsigned int length, int yAxisUnit, unsigned int xValueOption, unsigned int startIdx, unsigned int stopIdx, float threshold, float minPeakHeight_dB, unsigned int sortOption, peak_t *peakArray, unsigned int maxPeaks);


    /** Retrieves all peaks in the given data set between the given indices. 
    *   The properties of the found peaks are stored in the given arrays peakCenter, peakWidth and peakHeight.
    *   @param[in] xData The x-value for each data point, can be given in three different forms as specified in \ref sec_x_value_formats.
    *   @param[in] yData The y-value for each data point 
    *   @param[in] length The length of the yArray (and also the xData array if xValueOption equals #FTS_OPTION_ARRAY_X)
    *   @param[in] xValueOption Determines the options for how the x-data is interpreted, see \ref sec_x_value_formats.
    *   @param[in] yAxisUnit The unit of the y-axis data, must be one of the Y_UNIT_... constants defined in FTSData.h
    *   @param[in] startIdx The first index in the data set to use (inclusive)
    *   @param[in] stopIdx The last index in the data set to use (exclusive)
    *   @param[in] threshold Only peaks with intensity above this level will be considered (must be same unit as the y-data)
    *   @param[in] minPeakHeight_dB The minimum height of a peak to be returned, in dB
    *   @param[in] sortOption Determines how the peaks will be sorted
    *       sortOption = 0 -> no sorting
    *       sortOption = 1 -> peaks are sorted in order of increasing center position
    *       sortOption = 2 -> peaks are sorted in order of increasing peak height
    *       sortOption = 3 -> peaks are sorted in order of increasing peak width
    *   @param[out] peakCenter Will on successful return be filled with the center position of each peak found
    *   @param[out] peakWidth Will on successful return be filled with the width of each peak found
    *   @param[out] peakHeight Will on successful return be filled with the height of each peak found
    *   @param[out] peakLeftPt If not NULL, this will on successful return be filled with the first position to the left of the peak
    *       where the intensity has dropped by 'minPeakHeight' dB from the maximum value.
    *       This parameter may be NULL.
    *   @param[out] peakRightPt If not NULL, this will on successful return be filled with the first position to the right of the peak
    *       where the intensity has dropped by 'minPeakHeight' dB from the maximum value
    *       This parameter may be NULL.
    *   @param maxPeaks The maximum number of peaks that can be filled into the arrays peakCenter, peakWidth and peakHeight.
    *   @return The number of peaks found. -1 if an error occurs
    */
    int __cdecl   FTS_FindPeaks_Array(const float* xData, const float* yData, unsigned int length, int yAxisUnit, unsigned int xValueOption, unsigned int startIdx, unsigned int stopIdx, float threshold, float minPeakHeight_dB, unsigned int sortOption, float* peakCenter, float* peakWidth, float* peakHeight, float* peakLeftPt, float* peakRightPt, unsigned int maxPeaks);
    int __stdcall FTS_FindPeaks_Array_std(const float* xData, const float* yData, unsigned int length, int yAxisUnit, unsigned int xValueOption, unsigned int startIdx, unsigned int stopIdx, float threshold, float minPeakHeight_dB, unsigned int sortOption, float* peakCenter, float* peakWidth, float* peakHeight, float* peakLeftPt, float* peakRightPt, unsigned int maxPeaks);

    /** Retrieves all valleys in a collected spectrum/interferogram between the given indices. 
    *   @param[in] spec The spectrum in which to search for valleys
    *   @param[in] startIdx The first index in the data set to use (inclusive)
    *   @param[in] stopIdx The last index in the data set to use (exclusive)
    *   @param[in] threshold Only valleys where the highest part have and intensity above this level will be found
    *   @param[in] minValleyDepth_dB Only valleys which bottom are at least this number of dB below the baseline will be returned
    *   @param[in] sortOption Determines how the valleys will be sorted
    *       sortOption = 0 -> no sorting
    *       sortOption = 1 -> valleys are sorted in order of increasing center position
    *       sortOption = 2 -> valleys are sorted in order of increasing depth
    *       sortOption = 3 -> valleys are sorted in order of increasing width
    *   @param[out] valleyCenter Will on successful return be filled with the center position of each valley found
    *   @param[out] valleyWidth Will on successful return be filled with the width of each valley found
    *   @param[out] valleyDepth Will on successful return be filled with the depth of each valley found
    *   @param[out] valleyLeftPt If not NULL, this will on successful return be filled with the first position to the left of the valley
    *       where the intensity has increased by 'minValleyDepth_dB' dB from the minimum value.
    *       This parameter may be NULL.
    *   @param[out] valleyRightPt If not NULL, this will on successful return be filled with the first position to the right of the valley
    *       where the intensity has increased by 'minValleyDepth_dB' dB from the minimum value
    *       This parameter may be NULL.
    *   @param[in] maxValleys The maximum number of valleys that can be filled into the arrays valleyCenter, valleyWidth and valleyHeight.
    *   @return The number of valleys found. -1 if an error occurs
    */
    int __cdecl   FTS_FindValleys(spectrum_t* spec, unsigned int startIdx, unsigned int stopIdx, float threshold, float minValleyDepth_dB, unsigned int sortOption, float* valleyCenter, float* valleyWidth, float* valleyDepth, float* valleyLeftPt, float* valleyRightPt, unsigned int maxValleys);
    int __stdcall FTS_FindValleys_std(spectrum_t* spec, unsigned int startIdx, unsigned int stopIdx, float threshold, float minValleyDepth_dB, unsigned int sortOption, float* valleyCenter, float* valleyWidth, float* valleyDepth, float* valleyLeftPt, float* valleyRightPt, unsigned int maxValleys);

    /** Retrieves all valleys in a collected spectrum/interferogram between the given indices. 
    *   If you have trouble allocating/deleting or accessing data in struct:s then have a look at the routines #FTS_CreateArray_ValleyData, 
    *   #FTS_GetValleyProperty and #FTS_DeleteArray.
    *   @param[in] xData The x-value for each data point, can be given in three different forms as specified in \ref sec_x_value_formats.
    *   @param[in] yData The y-value for each data point 
    *   @param[in] length The length of the yArray (and also the xData array if xValueOption equals #FTS_OPTION_ARRAY_X)
    *   @param[in] xValueOption Determines the options for how the x-data is interpreted, see \ref sec_x_value_formats.
    *   @param[in] yAxisUnit The unit of the y-axis data, must be one of the Y_UNIT_... constants defined in FTSData.h
    *   @param[in] startIdx The first index in the data set to use (inclusive)
    *   @param[in] stopIdx The last index in the data set to use (exclusive)
    *   @param[in] threshold Only valleys where the highest part have and intensity above this level will be found
    *   @param[in] minValleyDepth_dB Only valleys which bottom are at least this number of dB below the baseline will be returned
    *   @param[in] sortOption Determines how the valleys will be sorted
    *       sortOption = 0 -> no sorting
    *       sortOption = 1 -> valleys are sorted in order of increasing center position
    *       sortOption = 2 -> valleys are sorted in order of increasing depth
    *       sortOption = 3 -> valleys are sorted in order of increasing width
    *   @param[out] valleyArray Will on successful return be filled with the information of each valley found
    *   @param[in] maxValleys The maximum number of valleys that can be filled into the arrays valleyCenter, valleyWidth and valleyHeight.
    *   @return The number of valleys found. -1 if an error occurs
    */
    int __cdecl   FTS_FindValleys_ext(const float* xData, const float* yData, unsigned int length, int yAxisUnit, unsigned int xValueOption, unsigned int startIdx, unsigned int stopIdx, float threshold, float minValleyDepth_dB, unsigned int sortOption, valley_t *valleyArray, unsigned int maxValleys);
    int __stdcall FTS_FindValleys_ext_std(const float* xData, const float* yData, unsigned int length, int yAxisUnit, unsigned int xValueOption, unsigned int startIdx, unsigned int stopIdx, float threshold, float minValleyDepth_dB, unsigned int sortOption, valley_t *valleyArray, unsigned int maxValleys);

    /** Retrieves all valleys in the given data set between the given indices. 
    *   The properties of the found valleys are stored in the given arrays valleyCenter, valleyWidth and valleyHeight.
    *   @param[in] xData The x-value for each data point (this may be NULL, in which case the index of each point will be used)
    *   @param[in] yData The y-value for each data point 
    *   @param[in] length The length of the yArray (and also the xData array if xValueOption equals #FTS_OPTION_ARRAY_X)
    *   @param[in] xValueOption determines the options for how the x-data is interpreted
    *       this is either of the following values;
    *       #FTS_OPTION_ARRAY_X  <-> x-values in array form
    *       #FTS_OPTION_MINMAX_X <-> x-values in min-max form (only first two values in 'xData' are used)
    *   @param[in] yAxisUnit The unit of the y-axis data, must be one of the Y_UNIT_... constants defined in FTSData.h
    *   @param[in] startIdx The first index in the data set to use (inclusive)
    *   @param[in] stopIdx The last index in the data set to use (exclusive)
    *   @param[in] threshold Only valleys where the highest part have and intensity above this level will be found
    *   @param[in] minValleyDepth_dB Only valleys which bottom are at least this number of dB below the baseline will be returned
    *   @param[in] sortOption Determines how the valleys will be sorted
    *       sortOption = 0 -> no sorting
    *       sortOption = 1 -> valleys are sorted in order of increasing center position
    *       sortOption = 2 -> valleys are sorted in order of increasing depth
    *       sortOption = 3 -> valleys are sorted in order of increasing width
    *   @param[out] valleyCenter Will on successful return be filled with the center position of each valley found
    *   @param[out] valleyWidth Will on successful return be filled with the width of each valley found
    *   @param[out] valleyDepth Will on successful return be filled with the depth of each valley found
    *   @param[out] valleyLeftPt If not NULL, this will on successful return be filled with the first position to the left of the valley
    *       where the intensity has increased by 'minValleyDepth_dB' dB from the minimum value.
    *       This parameter may be NULL.
    *   @param[out] valleyRightPt If not NULL, this will on successful return be filled with the first position to the right of the valley
    *       where the intensity has increased by 'minValleyDepth_dB' dB from the minimum value
    *       This parameter may be NULL.
    *   @param[in] maxValleys The maximum number of valleys that can be filled into the arrays valleyCenter, valleyWidth and valleyHeight.
    *   @return The number of valleys found. -1 if an error occurs
    */
    int __cdecl   FTS_FindValleys_Array(const float* xData, const float* yData, unsigned int length, int yAxisUnit, unsigned int xValueOption, unsigned int startIdx, unsigned int stopIdx, float threshold, float minValleyDepth_dB, unsigned int sortOption, float* valleyCenter, float* valleyWidth, float* valleyDepth, float* valleyLeftPt, float* valleyRightPt, unsigned int maxValleys);
    int __stdcall FTS_FindValleys_Array_std(const float* xData, const float* yData, unsigned int length, int yAxisUnit, unsigned int xValueOption, unsigned int startIdx, unsigned int stopIdx, float threshold, float minValleyDepth_dB, unsigned int sortOption, float* valleyCenter, float* valleyWidth, float* valleyDepth, float* valleyLeftPt, float* valleyRightPt, unsigned int maxValleys);

    /** Calculates the (low frequency) amplitude variation of the provided interferogram.
    *   @return #FTS_SUCCESS on successful calculation. */
    int __cdecl FTS_GetAmplitudeVariation(const spectrum_t* const interferogram, unsigned int polynomialOrder, unsigned int nIndicesToIgnore, double* amplitudeVariationInPercent);
    int __stdcall FTS_GetAmplitudeVariation_std(const spectrum_t* const interferogram, unsigned int polynomialOrder, unsigned int nIndicesToIgnore, double* amplitudeVariationInPercent);

    ///@}

#ifdef __cplusplus
}
#endif 

#endif