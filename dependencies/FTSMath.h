////////////////////////////////////////////////////////////////////////
// FTSMath.h
//    This is the math functions used by the Thorlabs Fourier Transform 
//    Spectrometer DLL. 
//    
//    Copyright (c) 2019 Thorlabs Sweden, All Rights Reserved
////////////////////////////////////////////////////////////////////////

#include "FTSData.h"

#ifndef FTSMATH_H
#define FTSMATH_H

#ifdef __cplusplus
extern "C" {
#endif 

    /** The complex_f structure is used for all calculations using complex data */
    typedef struct complex_f{
        float re;
        float im;

        /** Default constructor, the object is initialized to zero */
        complex_f() { re = 0.0F; im = 0.0F; }

        /** Setting the real component only */
        complex_f(float real) { re = real; im = 0.0F; }

        /** Setting the real and imaginary components */
        complex_f(float a, float b) { re = a; im = b; }

    }complex_f;

    enum STITCH_OPTION
    {
        STITCH_KEEP_ORIGINALS = 0,    /* Keep originals means that every data point in the two input spectra will be kept during stitch */
        STITCH_AVERAGE = 1            /* Average means that data points in the overlapping region will be the average of the two input spectra */
    };


    /** The math_operator specifies an operator that can be performed on one or more spectra.
    *   See the routines #FTS_SpectrumBinaryMath(), #FTS_SpectrumScalarMath() and #FTS_SpectrumUnaryMath()
    */
    enum math_operator
    {
        operator_unary_spectrum_sin,
        operator_unary_spectrum_cos,
        operator_unary_spectrum_exp,
        operator_unary_spectrum_log10,
        operator_unary_spectrum_ln,
        operator_unary_spectrum_inverse,
        operator_unary_spectrum_sqrt,
        operator_unary_spectrum_abs,
        operator_unary_spectrum_tan,

        operator_spectrum_addition_linspace,
        operator_spectrum_addition_logspace,
        operator_spectrum_subtraction_linspace,
        operator_spectrum_subtraction_logspace,
        operator_spectrum_multiplication,
        operator_spectrum_division,

        operator_scalar_addition,
        operator_scalar_subtraction,
        operator_scalar_subtraction_scalar_first,
        operator_scalar_multiplication,
        operator_scalar_division,
        operator_scalar_division_scalar_first,
        operator_scalar_power
    };

    // ----------------------------------------------------------------------------------------
    // -------------------------------- UNIT CONVERSTIONS--------------------------------------
    // ----------------------------------------------------------------------------------------

    /** \addtogroup grp_unit_conversions_numbers Unit Conversions - Numbers
    Conversion between different units. Where a conversion between vacuum and air is performed
    then also the atmospheric pressure, temperature and relative humidity needs to be specified.
    *  @{
    */

    double __cdecl   FTS_Convert_WaveNumber_To_nm_air(double waveNr, double atmPressure_hPa, double airTemp_c, double relativeHumididty_perc);
    double __cdecl   FTS_Convert_WaveNumber_To_nm_vac(double waveNr);
    double __cdecl   FTS_Convert_WaveNumber_To_THz(double waveNr);
    double __cdecl   FTS_Convert_WaveNumber_To_eV(double waveNr);

    double __cdecl   FTS_Convert_nm_air_To_WaveNumber(double wavelength, double atmPressure_hPa, double airTemp_c, double relativeHumididty_perc);
    double __cdecl   FTS_Convert_nm_vac_To_WaveNumber(double wavelength);
    double __cdecl   FTS_Convert_nm_air_To_THz(double wavelength, double atmPressure_hPa, double airTemp_c, double relativeHumididty_perc);
    double __cdecl   FTS_Convert_nm_vac_To_THz(double wavelength);
    double __cdecl   FTS_Convert_nm_air_To_eV(double wavelength, double atmPressure_hPa, double airTemp_c, double relativeHumididty_perc);
    double __cdecl   FTS_Convert_nm_vac_To_eV(double wavelength);

    double __cdecl   FTS_Convert_THz_To_WaveNumber(double frequency);
    double __cdecl   FTS_Convert_eV_To_WaveNumber(double energy);

    double __stdcall FTS_Convert_WaveNumber_To_nm_air_std(double waveNr, double atmPressure_hPa, double airTemp_c, double relativeHumididty_perc);
    double __stdcall FTS_Convert_WaveNumber_To_nm_vac_std(double waveNr);
    double __stdcall FTS_Convert_WaveNumber_To_THz_std(double waveNr);
    double __stdcall FTS_Convert_WaveNumber_To_eV_std(double waveNr);

    double __stdcall FTS_Convert_nm_air_To_WaveNumber_std(double wavelength, double atmPressure_hPa, double airTemp_c, double relativeHumididty_perc);
    double __stdcall FTS_Convert_nm_vac_To_WaveNumber_std(double wavelength);
    double __stdcall FTS_Convert_nm_air_To_THz_std(double wavelength, double atmPressure_hPa, double airTemp_c, double relativeHumididty_perc);
    double __stdcall FTS_Convert_nm_vac_To_THz_std(double wavelength);
    double __stdcall FTS_Convert_nm_air_To_eV_std(double wavelength, double atmPressure_hPa, double airTemp_c, double relativeHumididty_perc);
    double __stdcall FTS_Convert_nm_vac_To_eV_std(double wavelength);

    double __stdcall FTS_Convert_THz_To_WaveNumber_std(double frequency);
    double __stdcall FTS_Convert_eV_To_WaveNumber_std(double energy);
    ///@}

    /** \addtogroup grp_unit_conversions_resolutions Unit Conversions - Resolution
    Converting a spectral resolution from wavenumbers to another unit.
    *  @{
    */

    /** @param resolutionWnr the resolution in wavenumbers
    *   @param wavelength_nm_air the wavelength (in nm air) where the resolution should be calculated 
    *   @param atmPressure_hPa the atmospheric pressure (in hPa) at the time of the measurement of the spectrum
    *   @param airTemp_c the air temperature (in Celsius) at the time of the measurement of the spectrum 
    *   @param relativeHumididty_perc the relative humidity (in percent) at the time of the measurement of the spectrum */
    double __cdecl   FTS_ConvertResolution_WaveNumber_To_nm_air(double resolutionWnr, double wavelength_nm_air, double atmPressure_hPa, double airTemp_c, double relativeHumididty_perc);
    double __stdcall FTS_ConvertResolution_WaveNumber_To_nm_air_std(double resolutionWnr, double wavelength_nm_air, double atmPressure_hPa, double airTemp_c, double relativeHumididty_perc);

    /** @param resolutionWnr the resolution in wavenumbers
    *   @param wavelength_nm_vac the wavelength (in nm vacuum) where the resolution should be calculated */
    double __cdecl   FTS_ConvertResolution_WaveNumber_To_nm_vac(double resolutionWnr, double wavelength_nm_vac);
    double __stdcall FTS_ConvertResolution_WaveNumber_To_nm_vac_std(double resolutionWnr, double wavelength_nm_vac);

    /** @param resolutionWnr the resolution in wavenumbers */
    double __cdecl   FTS_ConvertResolution_WaveNumber_To_THz(double resolutionWnr);
    double __stdcall FTS_ConvertResolution_WaveNumber_To_THz_std(double resolutionWnr);

    /** @param resolutionWnr the resolution in wavenumbers */
    double __cdecl   FTS_ConvertResolution_WaveNumber_To_eV(double resolutionWnr);
    double __stdcall FTS_ConvertResolution_WaveNumber_To_eV_std(double resolutionWnr);

    ///@}

    /** \addtogroup grp_unit_conversions_spectra Unit Conversions - Spectra
    *  @{
    */

    /** This converts the given spectrum from the unit it is in to the desired unit
    *   @return #FTS_SUCCESS on success
    *   @return #FTS_ERROR_PARAMETER_ERROR if the supplied spectrum is an interferogram */
    int __cdecl   FTS_ConvertSpectrum_To_nm_air(spectrum_t* spec, bool normalize);
    int __cdecl   FTS_ConvertSpectrum_To_nm_vac(spectrum_t* spec, bool normalize);
    int __cdecl   FTS_ConvertSpectrum_To_WaveNumber(spectrum_t* spec, bool normalize);
    int __cdecl   FTS_ConvertSpectrum_To_THz(spectrum_t* spec, bool normalize);
    int __cdecl   FTS_ConvertSpectrum_To_eV(spectrum_t* spec, bool normalize);

    int __stdcall FTS_ConvertSpectrum_To_nm_air_std(spectrum_t* spec, bool normalize);
    int __stdcall FTS_ConvertSpectrum_To_nm_vac_std(spectrum_t* spec, bool normalize);
    int __stdcall FTS_ConvertSpectrum_To_WaveNumber_std(spectrum_t* spec, bool normalize);
    int __stdcall FTS_ConvertSpectrum_To_THz_std(spectrum_t* spec, bool normalize);
    int __stdcall FTS_ConvertSpectrum_To_eV_std(spectrum_t* spec, bool normalize);

    /** This converts a spectrum with the Y-axis unit #Y_UNIT_MW or #Y_UNIT_MW_NORM to have the Y-axis unit of #Y_UNIT_DBM or #Y_UNIT_DBM_NORM.
    *   @param[in,out] spec the spectrum to convert.
    *   @param[in] normalize if true then the spectrum will on return be in #Y_UNIT_DBM_NORM otherwise #Y_UNIT_DBM
    *   @return #FTS_SUCCESS on success 
    *   @return #FTS_ERROR_PARAMETER_ERROR if the input spectrum does not have the unit of mW or is not a spectrum.*/
    int __cdecl   FTS_ConvertSpectrum_To_dBm(spectrum_t* spec, bool normalize);
    int __stdcall FTS_ConvertSpectrum_To_dBm_std(spectrum_t* spec, bool normalize);

    /** This converts a spectrum with the Y-axis unit #Y_UNIT_MW to have the Y-axis unit of #Y_UNIT_DBM.
    *   @param[in,out] data the spectral data to convert.
    *   @param[in] length the length of the data array.
    *   @return #FTS_SUCCESS on success 
    *   @return #FTS_ERROR_PARAMETER_ERROR if the data is NULL or length is zero.*/
    int __cdecl   FTS_ConvertSpectrum_To_dBm_Array(float* data, unsigned int length);
    int __stdcall FTS_ConvertSpectrum_To_dBm_Array_std(float* data, unsigned int length);

    /** This converts a spectrum with the Y-axis unit #Y_UNIT_DBM or #Y_UNIT_DBM_NORM to have the Y-axis unit of #Y_UNIT_MW or #Y_UNIT_MW_NORM.
    *   @param[in,out] spec the spectrum to convert.
    *   @param[in] normalize if true then the spectrum will on return be in #Y_UNIT_MW_NORM otherwise #Y_UNIT_MW
    *   @return #FTS_SUCCESS on success 
    *   @return #FTS_ERROR_PARAMETER_ERROR if the input spectrum does not have the unit of dBm or is not a spectrum.*/
    int __cdecl   FTS_ConvertSpectrum_To_mW(spectrum_t* spec, bool normalize);
    int __stdcall FTS_ConvertSpectrum_To_mW_std(spectrum_t* spec, bool normalize);

    /** This converts a spectrum with the Y-axis unit #Y_UNIT_DBM to have the Y-axis unit of #Y_UNIT_MW.
    *   @param[in,out] data the spectral data to convert.
    *   @param[in] length the length of the data array.
    *   @return #FTS_SUCCESS on success 
    *   @return #FTS_ERROR_PARAMETER_ERROR if the data is NULL or length is zero.*/
    int __cdecl   FTS_ConvertSpectrum_To_mW_Array(float* data, unsigned int length);
    int __stdcall FTS_ConvertSpectrum_To_mW_Array_std(float* data, unsigned int length);

    /** This converts a spectrum to the desired x and y axis units.
    *   @param[in,out] spec the spectrum to convert
    *   @param[in] xAxisUnit the desired x axis unit.
    *   @param[in] yAxisUnit the desired y axis unit.
    *   @return #FTS_SUCCESS if the operation succeeds.
    *   @return #FTS_ERROR_CONVERSION_ERROR if either the current x axis unit cannot be converted to the
    *       desired xAxisUnit or if the current y axis unit cannot be converted to the desired yAxisUnit */
    int __cdecl   FTS_ConvertSpectrumTo(spectrum_t* spec, X_AXIS_UNIT xAxisUnit, Y_AXIS_UNIT yAxisUnit);
    int __stdcall FTS_ConvertSpectrumTo_std(spectrum_t* spec, X_AXIS_UNIT xAxisUnit, Y_AXIS_UNIT yAxisUnit);

    /** This converts a spectrum from Power Density (y-axis unit of #Y_UNIT_MW_NORM or #Y_UNIT_DBM_NORM) to
    *       Absolute Power (y-axis unit of #Y_UNIT_MW or #Y_UNIT_DBM).
    *   @param[in,out] spec The spectrum to convert.
    *   @return #FTS_SUCCESS if the operation succeeds.
    *   @return #FTS_ERROR_PARAMETER_ERROR if the current y-axis is not #Y_UNIT_MW_NORM or #Y_UNIT_DBM_NORM,
    *       OR if the input is not a spectrum.         */
    int __cdecl   FTS_ConvertSpectrum_To_AbsolutePower(spectrum_t* spec);
    int __stdcall FTS_ConvertSpectrum_To_AbsolutePower_std(spectrum_t* spec);

    /** This converts a spectrum from Power Density (y-axis unit of #Y_UNIT_MW_NORM or #Y_UNIT_DBM_NORM) to
    *       Absolute Power (y-axis unit of #Y_UNIT_MW or #Y_UNIT_DBM)
    *   @param[in] xData the x-axis values. This must be provided if the xAxisUnit is #X_UNIT_NM_AIR or #X_UNIT_NM_VAC. May otherwise be NULL.
    *   @param[in, out] spectrum the intensity data of the spectrum
    *   @param[in] length the length of the spectrum
    *   @param[in] xAxisUnit the unit of the x-axis data, must be one of the X_UNIT_... constants defined in FTSData.h
    *   @param[in] yAxisUnit the unit of the y-axis data, must be one of the Y_UNIT_... constants defined in FTSData.h
    *   @param[in] resolution_wnr the resolution of the spectral data, in wavenumbers
    *   @param[in] air_press_hPa the atmospheric pressure at the time of the collection of this spectrum (in hPa)
    *   @param[in] air_temp_C the atmospheric temperature at the time of the collection of this spectrum (in degrees Celsius)
    *   @param[in] air_relHum_percent the atmospheric relative humidity at the time of the collection of this spectrum (in percent) */
    int __cdecl   FTS_ConvertSpectrum_To_AbsolutePower_Array(float* xData, float* spectrum, unsigned int length, int xAxisUnit, int yAxisUnit, float resolution_wnr, float air_press_hPa, float air_temp_C, float air_relHum_percent);
    int __stdcall FTS_ConvertSpectrum_To_AbsolutePower_Array_std(float* xData, float* spectrum, unsigned int length, int xAxisUnit, int yAxisUnit, float resolution_wnr, float air_press_hPa, float air_temp_C, float air_relHum_percent);

    /** This converts a spectrum from Absolute Power (y-axis unit of #Y_UNIT_MW or #Y_UNIT_DBM) to
    *       Power Density (y-axis unit of #Y_UNIT_MW_NORM or #Y_UNIT_DBM_NORM).
    *   @param[in,out] spec The spectrum to convert.
    *   @return #FTS_SUCCESS if the operation succeeds.
    *   @return #FTS_ERROR_PARAMETER_ERROR if the current y-axis is not #Y_UNIT_MW or #Y_UNIT_DBM,
    *   OR if the input is not a spectrum.         */
    int __cdecl   FTS_ConvertSpectrum_To_PowerDensity(spectrum_t* spec);
    int __stdcall FTS_ConvertSpectrum_To_PowerDensity_std(spectrum_t* spec);

    /** This converts a spectrum from Absolute Power (y-axis unit of #Y_UNIT_MW or #Y_UNIT_DBM) to
    *       Power Density (y-axis unit of #Y_UNIT_MW_NORM or #Y_UNIT_DBM_NORM).
    *   @param[in] xData the x-axis values. This must be provided if the xAxisUnit is #X_UNIT_NM_AIR or #X_UNIT_NM_VAC. May otherwise be NULL.
    *   @param[in, out] spectrum the intensity data of the spectrum
    *   @param[in] length the length of the spectrum
    *   @param[in] xAxisUnit the unit of the x-axis data, must be one of the X_UNIT_... constants defined in FTSData.h
    *   @param[in] yAxisUnit the unit of the y-axis data, must be one of the Y_UNIT_... constants defined in FTSData.h
    *   @param[in] resolution_wnr the resolution of the spectral data, in wavenumbers
    *   @param[in] air_press_hPa the atmospheric pressure at the time of the collection of this spectrum (in hPa)
    *   @param[in] air_temp_C the atmospheric temperature at the time of the collection of this spectrum (in degrees Celsius)
    *   @param[in] air_relHum_percent the atmospheric relative humidity at the time of the collection of this spectrum (in percent) */
    int __cdecl   FTS_ConvertSpectrum_To_PowerDensity_Array(float* xData, float* spectrum, unsigned int length, int xAxisUnit, int yAxisUnit, float resolution_wnr, float air_press_hPa, float air_temp_C, float air_relHum_percent);
    int __stdcall FTS_ConvertSpectrum_To_PowerDensity_Array_std(float* xData, float* spectrum, unsigned int length, int xAxisUnit, int yAxisUnit, float resolution_wnr, float air_press_hPa, float air_temp_C, float air_relHum_percent);

    ///@}

    /** \addtogroup grp_physics_calculations Physics Calculations
    *  @{
    */

    /** This calculates the refractive index of air at the given atmospheric conditions 
    *   @param atmPressure_hPa The atmospheric pressure, in hekto Pascals
    *   @param airTemp_c The temperature of the air, in degrees Celsius
    *   @param relativeHumididty_perc The relative humidity of the air, in percent */
    double __cdecl   FTS_GetRefractiveIndexOfAir(double atmPressure_hPa, double airTemp_c, double relativeHumididty_perc);
    double __stdcall FTS_GetRefractiveIndexOfAir_std(double atmPressure_hPa, double airTemp_c, double relativeHumididty_perc);

    /** This calculates the refractive index of air at the given atmospheric conditions and vacuum wavelength 
    *   @param atmPressure_hPa The atmospheric pressure, in hecto Pascals
    *   @param airTemp_c The temperature of the air, in degrees Celsius
    *   @param relativeHumididty_perc The relative humidity of the air, in percent
    *   @param wavelength_vac_nm The vacuum wavelength, in nanometer */
    double __cdecl   FTS_GetRefractiveIndexOfAir_wavelengthDependent(double atmPressure_hPa, double airTemp_c, double relativeHumididty_perc, double wavelength_vac_nm);
    double __stdcall FTS_GetRefractiveIndexOfAir_wavelengthDependent_std(double atmPressure_hPa, double airTemp_c, double relativeHumididty_perc, double wavelength_vac_nm);

    /** This calculates an ideal black body spectrum for the given temperature.
    *   @param spectrum a spectrum buffer where the calculated black body spectrum will be stored. This must be filled in with
    *       a x-axis data and x-axis unit. The black body spectrum will be filled in in the 'I' array and calculated at the given
    *       x-axis values.
    *   @param temperature_K the temperature to calculate the black body spectrum for, given in Kelvins */
    int __cdecl FTS_CalculateBlackbodySpectrum(spectrum_t* spectrum, double temperature_K);
    int __stdcall FTS_CalculateBlackbodySpectrum_std(spectrum_t* spectrum, double temperature_K);

    ///@}


    // ----------------------------------------------------------------------------------------
    // ----------------------------------- STATISTICS ----------------------------------------
    // ----------------------------------------------------------------------------------------

    /** \addtogroup grp_general_math General Math
    *  @{
    */

    /** Calculates the mean value of the input data 'data' 
    *   @param[in] data The data to calculate the mean of.
    *   @param length The number of data points in data.
    *   @return the mean value */
    float __cdecl   FTS_CalculateMean(const float* data, unsigned int length);
    float __stdcall FTS_CalculateMean_std(const float* data, unsigned int length);

    /** Calculates the minimum and maximum values of the input data 'data'
    *   @param[in] data the data set to find the extreme values of.
    *   @param[in] length the length of the data array.
    *   @param[out] minValue will on successful return be filled with the smallest value in the data array.
    *   @param[out] maxValue will on successful return be filled with the largest value in the data array.
    *   @return #FTS_SUCCESS on success.
    *   @return #FTS_ERROR_PARAMETER_ERROR if either pointer is null or the length is zero. */
    int __cdecl   FTS_CalculateMinMax(const float* data, unsigned int length, float* minValue, float* maxValue);
    int __stdcall FTS_CalculateMinMax_std(const float* data, unsigned int length, float* minValue, float* maxValue);

    // ----------------------------------------------------------------------------------------
    // ----------------------------------- POWER OF TWO ---------------------------------------
    // ----------------------------------------------------------------------------------------

    unsigned int __cdecl   FTS_GetNearestHigherPowerOfTwo(unsigned int value);
    unsigned int __stdcall FTS_GetNearestHigherPowerOfTwo_std(unsigned int value);

    unsigned int __cdecl   FTS_GetNearestHigherPowerOfFour(unsigned int value);
    unsigned int __stdcall FTS_GetNearestHigherPowerOfFour_std(unsigned int value);

    /** Retrieves the SI Prefix for the given order or magnitude. The input order of magnitude must be an
    *   even multiple of three.
    *   e.g. order = -3 returns 'm' and order = +3 returns 'k' */
    char __cdecl   FTS_SI_GetPrefix(unsigned int order);
    char __stdcall FTS_SI_GetPrefix_std(unsigned int order);

    ///@}

    // ----------------------------------------------------------------------------------------
    // ------------------------------ CALCULATE SPECTRUM --------------------------------------
    // ----------------------------------------------------------------------------------------

    /** \addtogroup grp_sync_acquisition Synchronous Data Acquisition
    See Also \ref sec_syncCollection
    *  @{
    */

    /** This will take the given interferogram and process it according to the settings found in g_data
    *       for the given spectrometer. 
    *   @param[in] specIndex the index of the spectrometer. this determines 
    *       the set of parameters that will be used to process the interferogram. See \ref sec_numerating_spectrometers.
    *   @param[in, out] interferogram the interferogram to processes
    *   @return #FTS_SUCCESS on successful return */
    int __cdecl   FTS_ProcessInterferogram(SpectrometerIndex specIndex, spectrum_t* interferogram);
    int __stdcall FTS_ProcessInterferogram_std(SpectrometerIndex specIndex, spectrum_t* interferogram);

    /** This will take the given interferogram and process it according to the settings found in g_data
    *       for the given spectrometer. 
    *   @param[in] specIndex the index of the spectrometer. this determines 
    *       the set of parameters that will be used to process the interferogram. See \ref sec_numerating_spectrometers.
    *   @param[out] acqBuffer will be updated with data from the interferogram. What parameters will be updated
    *       depends on the settings in acqOptions, see the description for acquisition_option_t and acquisition_buffer_t.
    *   @param[in, out] interferogram the interferogram to processes
    *   @return #FTS_SUCCESS on successful return */
    int __cdecl   FTS_ProcessInterferogram_ext1(SpectrometerIndex specIndex, spectrum_t* interferogram, acquisition_buffer_t* acqBuffer);
    int __stdcall FTS_ProcessInterferogram_ext1_std(SpectrometerIndex specIndex, spectrum_t* interferogram, acquisition_buffer_t* acqBuffer);

    /** This will take the given interferogram and process it according to given settings. This performs the
    *       same operation as FTS_ProcessInterferogram() but using only the provided settings and buffers.
    *   @param[in, out] interferogram the interferogram to processes.
    *   @param[in] acqOptions contains the options for what should be done with the interferogram.
    *   @param[out] acqBuffer will be updated with data from the interferogram. What parameters will be updated
    *       depends on the settings in acqOptions, see the description for acquisition_option_t and acquisition_buffer_t.
    *   @return #FTS_SUCCESS on successful return */
    int __cdecl   FTS_ProcessInterferogram_ext(spectrum_t* interferogram, acquisition_option_t* acqOptions, acquisition_buffer_t* acqBuffer);
    int __stdcall FTS_ProcessInterferogram_ext_std(spectrum_t* interferogram, acquisition_option_t* acqOptions, acquisition_buffer_t* acqBuffer);

    /** This will take the given interferogram and process it according to the settings found in g_data. 
    *       The resulting spectrum will (on success) be stored in the supplied buffer spectrum.
    *   @param[in] specIndex the index of the spectrometer that collected the spectrum. this determines 
    *       the set of parameters that will be used to calculate the spectrum. See \ref sec_numerating_spectrometers.
    *   @param[in] interferogram the interferogram to calculate the spectrum from
    *   @param[out] spectrum will on successful return be filled with the calculated spectrum. On non-successful
    *       return is the contents undefined. Notice that this spectrum_t must be allocated with enough memory to hold
    *       the resulting spectrum prior to calling this routine
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_CalculateSpectrum(SpectrometerIndex specIndex, const spectrum_t* interferogram, spectrum_t* spectrum);
    int __stdcall FTS_CalculateSpectrum_std(SpectrometerIndex specIndex, const spectrum_t* interferogram, spectrum_t* spectrum);

    /** This will take the given interferogram and process it according to the settings provided. This performs the same
    *       operations as #FTS_CalculateSpectrum() but uses only the given data structures.
    *       The resulting spectrum will (on success) be stored in the supplied buffer spectrum.
    *   @param[in] interferogram the interferogram to calculate the spectrum from
    *   @param[out] spectrum will on successful return be filled with the calculated spectrum. On non-successful
    *       return is the contents undefined. Notice that this spectrum_t must be allocated with enough memory to hold
    *       the resulting spectrum prior to calling this routine
    *   @param[in] instrumentData the properties of the instrument that collected this data
    *   @param[in] acqOptions the options for calculating the spectrum 
    *   @param[in] acqBuffer if acquisition options specifies Mertz-phase correction then the calculated
    *       phase polynomial will be stored in acqBuffer. This parameter may be NULL.
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_CalculateSpectrum_ext(const spectrum_t* interferogram, spectrum_t* spectrum, const fts_instrument_t* instrumentData, acquisition_option_t* acqOptions, acquisition_buffer_t* acqBuffer);
    int __stdcall FTS_CalculateSpectrum_ext_std(const spectrum_t* interferogram, spectrum_t* spectrum, const fts_instrument_t* instrumentData, acquisition_option_t* acqOptions, acquisition_buffer_t* acqBuffer);

    // ----------------------------------------------------------------------------------------
    // ------------------------------- POWER CALIBRATION --------------------------------------
    // ----------------------------------------------------------------------------------------

    /** This takes the given spectrum and converts the Y-Axis units from being #Y_UNIT_COUNTS to being #Y_UNIT_MW. 
    *   @param[in,out] spec the spectrum to calibrate
    *   @param[in] specIndex the index of the spectrometer that collected the spectrum. This determines 
    *       the set of parameters that will be used to calculate the spectrum. See \ref sec_numerating_spectrometers.  */
    int __cdecl   FTS_PowerCalibrateSpectrum(spectrum_t* spec, SpectrometerIndex specIndex);
    int __stdcall FTS_PowerCalibrateSpectrum_std(spectrum_t* spec, SpectrometerIndex specIndex);

    /** This takes the given spectrum and converts the Y-Axis units from being #Y_UNIT_COUNTS to being #Y_UNIT_MW. 
    *   @param[in,out] spectralData the spectral data to calibrate.
    *   @param[in] wnrMin the wavenumber corresponding to the first data point in the spectrum
    *   @param[in] wnrMax the wavenumber corresponding to the last data point in the spectrum
    *   @param[in] length the length of the buffer 'spectralData'.
    *   @param[in] specIndex the index of the spectrometer that collected the spectrum. This determines 
    *       the set of parameters that will be used to calculate the spectrum. See \ref sec_numerating_spectrometers. 
    */
    int __cdecl   FTS_PowerCalibrateSpectrum_Array(float* spectralData, float wnrMin, float wnrMax, unsigned int length, SpectrometerIndex specIndex);
    int __stdcall FTS_PowerCalibrateSpectrum_Array_std(float* spectralData, float wnrMin, float wnrMax, unsigned int length, SpectrometerIndex specIndex);

    // ----------------------------------------------------------------------------------------
    // --------------------------------- PHASE CORRECTION -------------------------------------
    // ----------------------------------------------------------------------------------------

    /** This will make a phase correction of the given (complex) spectrum.
    *   @param[in] spectralData the spectral data of the spectrum to correct.
    *   @param[out] phaseCorrectedSpec the resulting spectral data after the phase correction. This must be allocated before calling
    *       this routine, with the length (maxIdx - minIdx).
    *   @param[in] spectralDataLength the length of the spectralData array.
    *   @param[in] intSamplingDistance_cm the distance (in Optical Path Difference in cm) between sample points in the    
    *       original interferogram.
    *   @param[in] minIdx, maxIdx The lowest index and highest in spectralData which represents an actual signal. The output phaseCorrectedSpec will be
    *       cut to the range minIdx -> maxIdx and the output will have the length (maxIdx - minIdx).
    *   @param[in] phaseArrayOption If this equals one, then the phase will be read from the supplied
    *       arrays phaseX, phaseY and phaseLength. If this does not equal one then the phase will be calculated from the spectralData and 
    *       stored in the supplied arrays phaseX, phaseY and phaseLength (if not NULL) which then must be preallocated with the length spectralDataLength/2.
    *   @param[in, out] phaseX, phaseY, phaseLength If phaseArrayOption equals one then these arrays will be filled with the calculated phase of the 
    *       spectrum. If phaseArrayOptions is not equal to one then these are the source of the phase information.
    */
    int __cdecl   FTS_PhaseCorrectSpectrum(const complex_f* spectralData, float* phaseCorrectedSpec, unsigned int spectralDataLength, double intSamplingDistance_cm, unsigned int minIdx, unsigned int maxIdx, int phaseArrayOption, double* phaseX, double* phaseY, int* phaseLength);
    int __stdcall FTS_PhaseCorrectSpectrum_std(const complex_f* spectralData, float* phaseCorrectedSpec, unsigned int spectralDataLength, double intSamplingDistance_cm, unsigned int minIdx, unsigned int maxIdx, int phaseArrayOption, double* phaseX, double* phaseY, int* phaseLength);

    /** This will make a phase correction of the given interferogram
    *   @param[in,out] interferogram The interferogram to perform the correction on.
    *   @param[in] phaseResolution The desired resolution of the calculated phase array, in inverse centimeters.
    *   @param[in] minWavenr The smallest wavenumber that is physically present in the interferogram.
    *   @param[in] maxWavenr The largest wavenumber that is physically present in the interferogram.
    *   @param[in,out] phaseX, phaseY, phaseLength
    *   @param[in] phaseArrayOption If this equals one then the phase will be read from the supplied arrays phaseX, phaseY and phaseLenght.
    *       If this does not equal one then the phase will be calculated from the interferogram and the result will be stored in the supplied
    *       arrays phaseX, phaseY and phaseLength. These should then be preallocated with the length of the interferogram/2.
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_PhaseCorrectInterferogram(spectrum_t* interferogram, float phaseResolution, float minWavenr, float maxWavenr, double* phaseX, double* phaseY, int* phaseLength, int phaseArrayOption);
    int __stdcall FTS_PhaseCorrectInterferogram_std(spectrum_t* interferogram, float phaseResolution, float minWavenr, float maxWavenr, double* phaseX, double* phaseY, int* phaseLength, int phaseArrayOption);

    /** This will make a phase correction of the given interferogram
    *   @param[in,out] interferogramData The interferogram to perform the correction on.
    *   @param[in] interferogramLength The length of the interferogram
    *   @param[in] intSamplingDistance_cm The Optical Path Difference between two sample points in the interferogram, in centimeters.
    *   @param[in] phaseResolution The desired resolution of the calculated phase array, in inverse centimeters.
    *   @param[in] minWavenr The smallest wavenumber that is physically present in the interferogram.
    *   @param[in] maxWavenr The largest wavenumber that is physically present in the interferogram.
    *   @param[in,out] phaseX, phaseY, phaseLength
    *   @param[in] phaseArrayOption If this equals one then the phase will be read from the supplied arrays phaseX, phaseY and phaseLenght.
    *       If this does not equal one then the phase will be calculated from the interferogram and the result will be stored in the supplied
    *       arrays phaseX, phaseY and phaseLength. These should then be preallocated with the length of the interferogram/2.
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_PhaseCorrectInterferogram_Array(float* interferogramData, unsigned int interferogramLength, double intSamplingDistance_cm, float phaseResolution, float minWavenr, float maxWavenr, double* phaseX, double* phaseY, int* phaseLength, int phaseArrayOption);
    int __stdcall FTS_PhaseCorrectInterferogram_Array_std(float* interferogramData, unsigned int interferogramLength, double intSamplingDistance_cm, float phaseResolution, float minWavenr, float maxWavenr, double* phaseX, double* phaseY, int* phaseLength, int phaseArrayOption);

    ///@}

    // ----------------------------------------------------------------------------------------
    // -------------------------------- FOURIER ANALYSIS --------------------------------------
    // ----------------------------------------------------------------------------------------

    /** \addtogroup grp_spectrum_FFT Spectrum FFT
    *    \ingroup grp_spectrum_math
    *  @{
    */

    /** This calculates the absolute of the Discrete Fourier Transform of a real-valued input data sequence.
    *       Since the input data series is real is the output series symmetrical around the DC component and only
    *       half of the values are returned. The other half can be calculated from symmetry.
    *   @param[in] inputData the data to calculate the Fourier Transform of.
    *   @param[out] outputData will on successful return be filled with the absolute of the FFT(inputData).
    *       Must be allocated with the length fftLen/2 before calling this function.
    *   @param[in] fftLen the length of inputData. Must be an even power of two (2,4,8,16,32,64,...)
    *   @return #FTS_SUCCESS if all is ok */
    int __cdecl   FTS_FFT_Real_Abs(const float* inputData, float* outputData, unsigned int fftLen);
    int __stdcall FTS_FFT_Real_Abs_std(const float* inputData, float* outputData, unsigned int fftLen);

    /** This calculates the the Discrete Fourier Transform of a real-valued input data sequence.
    *       Since the input data series is real is the output series symmetrical around the DC component and only
    *       half of the values are returned. The other half can be calculated from symmetry.
    *   @param[in] inputData the data to calculate the Fourier Transform of.
    *   @param[out] DFT will on successful return be filled with the calculated FFT(inputData).
    *       Must be allocated with the length fftLen/2 before calling this function.
    *   @param[in] fftLen the length of inputData. Must be an even power of two (2,4,8,16,32,64,...)
    *   @return #FTS_SUCCESS on success. */
    int __cdecl   FTS_FFT_Real(const float* inputData, complex_f* DFT, unsigned int fftLen);
    int __stdcall FTS_FFT_Real_std(const float* inputData, complex_f* DFT, unsigned int fftLen);

    /** Calculates the Discrete Fourier Transform of a complex-valued data series 
    *   @param[in] inputData the data to calculate the Fourier Transform of.
    *   @param[out] DFT will on successful return be filled with the calculated FFT(inputData).
    *       This must be allocated with the length fftLen before calling this function.
    *   @param[in] fftLen the length of inputData. This value must contain only the prime factors 2, 3 or 7,
    *       if this contains any other prime factor then this will return #FTS_ERROR_NOT_POWEROFTWO 
    *   @return #FTS_SUCCESS on success. */
    int __cdecl   FTS_FFT(const complex_f* inputData, complex_f* DFT, unsigned int fftLen);
    int __stdcall FTS_FFT_std(const complex_f* inputData, complex_f* DFT, unsigned int fftLen);

    /** Calculates the Inverse Discrete Fourier Transform of a complex-valued data series 
    *   @param[in] inputData the data to calculate the inverse Fourier Transform of.
    *   @param[out] DFT will on successful return be filled with the calculated IFFT(inputData).
    *       Must be allocated with the length fftLen before calling this function.
    *   @param[in] fftLen the length of inputData. Must be an even power of two (2,4,8,16,32,64,...)
    *   @return #FTS_SUCCESS on success.      */
    int __cdecl   FTS_IFFT(const complex_f* inputData, complex_f* DFT, unsigned int fftLen);
    int __stdcall FTS_IFFT_std(const complex_f* inputData, complex_f* DFT, unsigned int fftLen);

    /** Calculates the Inverse Discrete Fourier Transform of a real-valued data series.
    *   @param[in] inputData the data to calculate the inverse Fourier Transform of.
    *   @param[out] DFT will on successful return be filled with the calculated IFFT(inputData).
    *       Must be allocated with the length fftLen before calling this function.
    *   @param[in] fftLen the length of inputData. Must be an even power of two (2,4,8,16,32,64,...)
    *   @return #FTS_SUCCESS on success.     */
    int __cdecl   FTS_IFFT_Real(const float* inputData, complex_f* DFT, unsigned int fftLen);
    int __stdcall FTS_IFFT_Real_std(const float* inputData, complex_f* DFT, unsigned int fftLen);

    /** @return the nearest allowed length for an FFT calculation. If the queried length is allowed
    *   then the same value is returned */
    unsigned int __cdecl FTS_GetNearestAllowedFFTLength(unsigned int queriedLength);
    unsigned int __stdcall FTS_GetNearestAllowedFFTLength_std(unsigned int queriedLength);

    /** Calculates the Discrete Fourier Transform of a complex-valued spectrum.  
    *   The result of the calculation is a complex-valued spectrum and the provided spectrum
    *   will NOT be reallocated so the argument needs to be allocated with a 'phi' array.  
    *   @param[in, out] spectrum the data to calculate the Fourier Transform of. After successful return will this
    *       also hold the result of the calculation. This may be a complex-valued spectrum (see \ref sec_complex_spectral_data)?
    *   @return #FTS_ERROR_OUTOFMEMORY if the spectrum does not have an 'phi' array allocated, or the size of the 'phi' array is smaller than the length of the spectrum.
    *   @return #FTS_ERROR_NOT_POWEROFTWO if the length of the spectrum is not an allowed length for the FFT calculations.
    *   @return #FTS_SUCCESS on successful calculation of the transform. */
    int __cdecl   FTS_FFTSpectrum(spectrum_t* spectrum);
    int __stdcall FTS_FFTSpectrum_std(spectrum_t* spectrum);

    /** Calculates the Inverse Discrete Fourier Transform of a complex-valued spectrum.  
    *   The result of the calculation is a complex-valued spectrum and the provided spectrum
    *   will NOT be reallocated so the argument needs to be allocated with a 'phi' array.  
    *   @param[in, out] spectrum the data to calculate the Inverse Fourier Transform of. After successful return will this
    *       also hold the result of the calculation. This may be a complex-valued spectrum (see \ref sec_complex_spectral_data)?
    *   @return #FTS_ERROR_OUTOFMEMORY if the spectrum does not have an 'phi' array allocated, or the size of the 'phi' array is smaller than the length of the spectrum.
    *   @return #FTS_ERROR_NOT_POWEROFTWO if the length of the spectrum is not an allowed length for the FFT calculations.
    *   @return #FTS_SUCCESS on successful calculation of the transform. */
    int __cdecl   FTS_IFFTSpectrum(spectrum_t* spectrum);
    int __stdcall FTS_IFFTSpectrum_std(spectrum_t* spectrum);

    /** This performs an FFT-shift of the input data, the data is 
    *   re-organized by putting the zero-frequency component to the 
    *   center of the data. 
    *   This is necessary if one wishes to take the inverse of an already
    *       calculated FT. 
    *   @param[in,out] data the data to reorganize.
    *   @param length the length of the data array. */
    int __cdecl   FTS_FFTShift_Real(float* data, unsigned int length);
    int __stdcall FTS_FFTShift_Real_std(float* data, unsigned int length);

    /** This performs an FFT-shift of the input data, the data is 
    *   re-organized by putting the zero-frequency component to the 
    *   center of the data. 
    *   This is necessary if one wishes to take the inverse of an already
    *       calculated FT. 
    *   @param[in,out] data the data to reorganize.
    *   @param length the length of the data array. */
    int __cdecl   FTS_FFTShift_Comp(complex_f* data, unsigned int length);
    int __stdcall FTS_FFTShift_Comp_std(complex_f* data, unsigned int length);

    ///@}

    /** \addtogroup grp_spectrum_math Spectrum Math
    *  @{
    */

    // ----------------------------------------------------------------------------------------
    // ----------------------------------- COMPLEX MATH ---------------------------------------
    // ----------------------------------------------------------------------------------------

    /** This calculates the absolute value of the complex input data series 
    *   @param[in] compData The complex data to calculate the absolute of.
    *   @param[out] abs Will on successful return be filled with the absolute of the input data.
    *       This needs to be preallocated with the same length as compData.
    *   @param[in] length The length of the input and output arrays. */
    int __cdecl FTS_Complex_Abs(const complex_f* compData, float* abs, unsigned int length);
    int __stdcall FTS_Complex_Abs_std(const complex_f* compData, float* abs, unsigned int length);

    /** This calculates the absolute value of the complex input spectrum 
    *   @param[in, out] complexSpectrum The complex spectrum to calculate the absolute of. This may be a complex-valued spectrum (see \ref sec_complex_spectral_data). */
    int __cdecl FTS_ComplexSpectrum_Abs(spectrum_t* complexSpectrum);
    int __stdcall FTS_ComplexSpectrum_Abs_std(spectrum_t* complexSpectrum);

    /** This calculates the phase (in radians) of the complex input data series 
    *   @param[in] compData The complex data to calculate the phase of.
    *   @param[out] phi Will on successful return be filled with the phase of the input data, in radians.
    *       This needs to be preallocated with the same length as compData.
    *   @param[in] length The length of the input and output arrays. */
    int __cdecl FTS_Complex_Phase(const complex_f* compData, float* phi, unsigned int length);
    int __stdcall FTS_Complex_Phase_std(const complex_f* compData, float* phi, unsigned int length);

    /** This calculates the phase (in radians) of the complex input spectrum 
    *   @param[in, out] complexSpectrum The complex spectrum to calculate the phase of. This may be a complex-valued spectrum (see \ref sec_complex_spectral_data). */
    int __cdecl FTS_ComplexSpectrum_Phase(spectrum_t* complexSpectrum);
    int __stdcall FTS_ComplexSpectrum_Phase_std(spectrum_t* complexSpectrum);

    /** Unwraps an array with angle information given in radians.
    *   @param[in,out] phaseData the phase data to unwrap
    *   @param length the length of the phase information */
    int __cdecl FTS_Unwrap_Array_Radians(float* phaseData, unsigned int length);
    int __stdcall FTS_Unwrap_Array_Radians_std(float* phaseData, unsigned int length);

    /** Unwraps the phase array of the provided spectrum. This will only affect the 'phi' array and will only have any effect
    *   if phiValueFormat equals #PHI_VAL_ARRAY.
    *   @param[in,out] complexSpectrum the spectrum who's phase should be unwrapped.
    *   @return #FTS_SUCCESS if the phase was successfully unwrapped.
    *   @return #FTS_ERROR_NODATA if no phase data exists in the given spectrum. */
    int __cdecl FTS_UnwrapSpectrum(spectrum_t* complexSpectrum);
    int __stdcall FTS_UnwrapSpectrum_std(spectrum_t* complexSpectrum);

    ///@}


    ///@}

    /** \addtogroup grp_spectrum_function_fit Function Fitting
    *    \ingroup grp_spectrum_math
    *  @{
    */

    /** Fits one Gaussian function to the given spectrum in the given range of pixels 
    *   The gaussian is on the form: baseline + A * exp( -(x-mu)^2 / (2 * sigma))
    *   The retrieved parameters 'baseline', 'A', 'mu' and 'sigma' are returned
    *   @param[in] xData can be given in three different forms as specified in \ref sec_x_value_formats.
    *   @param[in] yData The data set to fit the Gaussian to 
    *   @param length The length of the input data series
    *   @param[in] xValueOption determines the options for how the x-data is interpreted, see \ref sec_x_value_formats.
    *   @param startIdx the data point to start the calculation at (inclusive)
    *   @param stopIdx the data point to stop the calculation at (exclusive)
    *   @param[out] mu Will on successful return be filled with the calculated center point 
    *   @param[out] sigma Will on successful return be filled with the calculated width parameter 
    *   @param[out] A Will on successful return be filled with the calculated amplitude parameter
    *   @param[out] baseline Will on successful return be filled with the calculated baseline
    *   @return #FTS_SUCCESS on successful fitting.     */
    int __cdecl   FTS_FitGaussian(const float* xData, const float* yData, unsigned int length, unsigned int startIdx, unsigned int stopIdx, unsigned int xValueOption, double* mu, double* sigma, double* A, double* baseline);    
    int __stdcall FTS_FitGaussian_std(const float* xData, const float* yData, unsigned int length, unsigned int startIdx, unsigned int stopIdx, unsigned int xValueOption, double* mu, double* sigma, double* A, double* baseline);    

    /** Fits one Lorentzian function to the given spectrum in the given range of pixels 
    *   The lorentzian is on the form: baseline + H / ( 1 + 4*(x-x0)^2 / w^2 )
    *   The retrieved parameters 'baseline', 'H', 'x0' and 'w' are returned
    *   @param[in] xData can be given in three different forms as specified in \ref sec_x_value_formats.
    *   @param[in] yData The data set to fit the Lorentzian to 
    *   @param length The length of the input data series
    *   @param[in] xValueOption determines the options for how the x-data is interpreted, see \ref sec_x_value_formats.
    *   @param startIdx the data point to start the calculation at (inclusive)
    *   @param stopIdx the data point to stop the calculation at (exclusive)
    *   @param[out] x0 Will on successful return be filled with the calculated center point 
    *   @param[out] w Will on successful return be filled with the calculated width parameter 
    *   @param[out] H Will on successful return be filled with the calculated amplitude parameter
    *   @param[out] baseline Will on successful return be filled with the calculated baseline
    *   @return #FTS_SUCCESS on successful fitting.     */
    int __cdecl   FTS_FitLorentzian(const float* xData, const float* yData, unsigned int length, unsigned int startIdx, unsigned int stopIdx, unsigned int xValueOption, double* x0, double* w, double* H, double* baseline);    
    int __stdcall FTS_FitLorentzian_std(const float* xData, const float* yData, unsigned int length, unsigned int startIdx, unsigned int stopIdx, unsigned int xValueOption, double* x0, double* w, double* H, double* baseline);    

    /** Fits one Polynomial of order 'polyOrder' to the given spectrum in the given range of pixels 
    *   @param[in] xData can be given in three different forms as specified in \ref sec_x_value_formats.
    *   @param[in] yData The data set to fit the polynomial to 
    *   @param length The length of the input data series
    *   @param[in] xValueOption determines the options for how the x-data is interpreted, see \ref sec_x_value_formats.
    *   @param startIdx the data point to start the calculation at (inclusive)
    *   @param stopIdx the data point to stop the calculation at (exclusive)
    *   @param polyOrder The desired order of the fitted polynomial
    *   @param[out] polynomial Will on successful return be filled with the parameters of the fitted
    *   polynomial. polynomial[0] is the constant and polynomial[N] is the highest order coefficient
    *   NOTE: polynomial must be preallocated with enough length to fit in the calculated polynomial
    *   (i.e. length must be at least polyOrder+1)
    *   @return #FTS_SUCCESS on successful fitting.     */
    int __cdecl   FTS_FitPolynomial(const float* xData, const float* yData, unsigned int length, unsigned int startIdx, unsigned int stopIdx, unsigned int xValueOption, unsigned int polyOrder, double* polynomial);    
    int __stdcall FTS_FitPolynomial_std(const float* xData, const float* yData, unsigned int length, unsigned int startIdx, unsigned int stopIdx, unsigned int xValueOption, unsigned int polyOrder, double* polynomial);    

    /** Fits a baseline polynomial of the desired order to the given spectrum in the desired range 
    *   @param data This is the spectral data to fit the baseline polynomial to
    *   @param length The length of the array 'data'
    *   @param[in] x If not NULL then this specifies the indices in 'data' which will be used
    *       to fit the polynomial. If NULL then these points will be determined automatically
    *       If not NULL then the length of this array must be at least (polyOrder+1)
    *   @param xLength The length of the array 'x'. If 'x' is NULL then this will be ignored
    *   @param polyOrder The desired order of the fitted polynomial
    *   @param[out] polynomial Will on successful return be filled with the parameters of the fitted
    *       polynomial. polynomial[0] is the constant and polynomial[N] is the highest order coefficient
    *   NOTE: polynomial must be preallocated with enough length to fit in the calculated polynomial
    *   (i.e. length must be at least polyOrder+1)
    */
    int __cdecl   FTS_FitBaselinePolynomial(const float* data, unsigned int length, unsigned int* x, unsigned int xLength, unsigned int polyOrder, double* polynomial);
    int __stdcall FTS_FitBaselinePolynomial_std(const float* data, unsigned int length, unsigned int* x, unsigned int xLength, unsigned int polyOrder, double* polynomial);

    ///@}

    /** \addtogroup grp_spectrum_math Spectrum Math
    *  @{
    */
    // ----------------------------------------------------------------------------------------
    // ------------------------------ GENERAL SPECTRUM MATH -----------------------------------
    // ----------------------------------------------------------------------------------------

    /** Performs an apodisation of the given interferogram 
    *   @param[in,out] interferogram The interferogram to perform the apodisation on
    *   @param[in] apodisationType Must be one of the APODISATION_ constants defined in FTSData.h
    *   @return #FTS_SUCCESS upon success
    */
    int __cdecl   FTS_Apodisation(spectrum_t* interferogram, int apodisationType);    
    int __stdcall FTS_Apodisation_std(spectrum_t* interferogram, int apodisationType);    

    /** Performs an apodisation of the given interferogram 
    *   @param[in,out] interferogram The interferogram to perform the apodisation on
    *   @param[in] length The length of the input data
    *   @param[in] apodisationType must be one of the APODISATION_ constants defined in FTSData.h
    *   @return #FTS_SUCCESS upon success
    */
    int __cdecl   FTS_Apodisation_Array(float* interferogram, unsigned int length, int apodisationType);
    int __stdcall FTS_Apodisation_Array_std(float* interferogram, unsigned int length, int apodisationType);

    /** This applies a given level of smoothing to the spectral data
    *   @param[in,out] spectrum The spectrum/interferogram to smooth
    *   @param[in] type Must equal one of the SMOOTHING_... constants defined in FTSData.h
    *   @param[in] param1, param2 The options for the smoothing.
    *   @param[in] startIdx, stopIdx The range of indices in the spectrum between which the smoothing should be performed.
    *   
    *   The interpretation of param1 and param2 depends on the value of 'type'
    *   Gaussian:        param1 is the width of the Gaussian and param2 is the number of passes over the data
    *   Moving average: param1 is the width of the window and param2 is the number of passes over the data
    *   Binomial:        param1 is the width of the window and param2 is the number of passes over the data
    *   Savitsky-Golay:    param1 is the width of the window and param2 is the order of the fitted polynomial
    *   Fourier:        param1 is the fraction of 'length' (in percent, ranging from 0 to 100) to keep in 
    *   the Fourier transform, param2 is the shape of the filter, 0<->boxcar, 1<->parabola, 2<->cosine
    */
    int __cdecl   FTS_Smooth(spectrum_t* spectrum, int type, int param1, int param2, unsigned int startIdx, unsigned int stopIdx);
    int __stdcall FTS_Smooth_std(spectrum_t* spectrum, int type, int param1, int param2, unsigned int startIdx, unsigned int stopIdx);

    /** This applies a given level of smoothing to the spectral data
    *   @param[in,out] data The spectrum/interferogram to smooth
    *   @param[in] type Must equal one of the SMOOTHING_... constants defined in FTSData.h
    *   @param[in] param1, param2 The options for the smoothing.
    *   @param[in] startIdx, stopIdx The range of indices in the spectrum between which the smoothing should be performed.
    *   
    *   The interpretation of param1 and param2 depends on the value of 'type'
    *   Gaussian:        param1 is the width of the Gaussian and param2 is the number of passes over the data
    *   Moving average: param1 is the width of the window and param2 is the number of passes over the data
    *   Binomial:        param1 is the width of the window and param2 is the number of passes over the data
    *   Savitsky-Golay:    param1 is the width of the window and param2 is the order of the fitted polynomial
    *   Fourier:        param1 is the fraction of 'length' (in percent, ranging from 0 to 100) to keep in 
    *   the Fourier transform, param2 is the shape of the filter, 0<->boxcar, 1<->parabola, 2<->cosine
    */
    int __cdecl   FTS_Smooth_Array(float* data, int type, int param1, int param2, unsigned int startIdx, unsigned int stopIdx);
    int __stdcall FTS_Smooth_Array_std(float* data, int type, int param1, int param2, unsigned int startIdx, unsigned int stopIdx);

    /** This calculates the derivative of the given data, which is assumed to be evenly sampled. 
    *   @param[in,out] spectrum The data to calculate the derivative of
    *   @param[in] order The order of the derivative to calculate (if using Savitsky-Golay then this must be <= 5)
    *   @param[in] type Determines the method to use for the calculation
    *   #DERIVATIVE_FINITE_DIFFERENCE <-> finite difference
    *   #DERIVATIVE_SAVITSKY_GOLAY     <-> Savitsky-Golay    */
    int __cdecl   FTS_Derivative(spectrum_t* spectrum, unsigned char order, unsigned char type);
    int __stdcall FTS_Derivative_std(spectrum_t* spectrum, unsigned char order, unsigned char type);

    /** This calculates the derivative of the given data, which is assumed to be 
    *   evenly sampled. 
    *   @param[in,out] yData  The data to calculate the derivative of
    *   @param[in] length The length of the input data series
    *   @param[in] xData can be given in three different forms as specified in \ref sec_x_value_formats.
    *   @param[in] xValueOption determines the options for how the x-data is interpreted, see \ref sec_x_value_formats.
    *   @param[in] order The order of the derivative to calculate (if using Savitsky-Golay then this must be <= 5)
    *   @param[in] type Determines the method to use for the calculation
    *   DERIVATIVE_FINITE_DIFFERENCE <-> finite difference
    *   DERIVATIVE_SAVITSKY_GOLAY     <-> Savitsky-Golay    */
    int __cdecl   FTS_Derivative_Array(const float* xData, float* yData, unsigned int length, unsigned int xValueOption, unsigned char order, unsigned char type);
    int __stdcall FTS_Derivative_Array_std(const float* xData, float* yData, unsigned int length, unsigned int xValueOption, unsigned char order, unsigned char type);

    /** This adds two spectra using the given weighting factor
    *   @param[in,out] spectrum1 Will on successful return be equal to spectrum1 + weight * spectrum2
    *   (note: if weight is negative then the spectra will be subtracted)
    *   @param[in] spectrum2 The second term in the calculation
    *   @param[in] weight The weighting factor for spectrum2. See explanation for spectrum1. 
    *   This will increase the averageNum parameter for spectrum1 and make the 
    *   creation time an average of the times in spectrum1 and spectrum2
    *   @return #FTS_ERROR_DIMENSION_MISMATCH if the lengths of the two spectra are different
    *   @return #FTS_ERROR_DIMENSION_MISMATCH if the x-axis values of the two spectra differ
    *   @return #FTS_SUCCESS on success     */
    int __cdecl   FTS_AddSpectra(spectrum_t* spectrum1, const spectrum_t* spectrum2, float weight);
    int __stdcall FTS_AddSpectra_std(spectrum_t* spectrum1, const spectrum_t* spectrum2, float weight);

    /** This adds two spectra using the given weighting factor
    *   @param[in,out] spectrum1 Will on successful return be equal to spectrum1 + weight * spectrum2
    *   (note: if weight is negative then the spectra will be subtracted)
    *   @param[in] spectrum2 The second term in the calculation
    *   @param[in] weight The weighting factor for spectrum2. See explanation for spectrum1.
    *   @param[in] length The length of the spectra spectrum1 and spectrum2.
    *   This will increase the averageNum parameter for spectrum1 and make the 
    *   creation time an average of the times in spectrum1 and spectrum2
    *   @return #FTS_SUCCESS on success     */
    int __cdecl   FTS_AddSpectra_Array(float* spectrum1, const float* spectrum2, float weight, unsigned int length);
    int __stdcall FTS_AddSpectra_Array_std(float* spectrum1, const float* spectrum2, float weight, unsigned int length);

    /** This adds the second spectrum to the first as an incremental moving average. 
    *   @param[in,out] spectrum1 Will on successful return be equal to either:
    *   <ul>
    *   <li>If N < specToAverageNum:  (spectrum1 * N + spectrum2)/(N+1) where N is the number of averages already in spectrum1.</li>
    *   <li>If N >= specToAverageNum: (spectrum1 * (specToAverageNum-1) + spectrum2)/specToAverageNum where N is the number of averages already in spectrum1.</li>
    *   </ul>
    *   @param[in] spectrum2 Will be added to spectrum1.
    *   @param[in] specToAverageNum The maximum number to average together. See description for spectrum1.
    */
    int __cdecl   FTS_UpdateRollingSpectrumAverage(spectrum_t* spectrum1, const spectrum_t* spectrum2, unsigned int specToAverageNum);
    int __stdcall FTS_UpdateRollingSpectrumAverage_std(spectrum_t* spectrum1, const spectrum_t* spectrum2, unsigned int specToAverageNum);

    /** This multiplies two spectra point wise
    *   On return will spectrum1 be equal to spectrum1 * spectrum2
    *   @param[in,out] spectrum1 The first factor and also the output spectrum
    *   @param[in] spectrum2 The second factor.
    *   @return #FTS_ERROR_DIMENSION_MISMATCH if the lengths of the two spectra are different
    *   @return #FTS_ERROR_DIMENSION_MISMATCH if the x-axis values of the two spectra differ
    *   @return #FTS_SUCCESS on success     */
    int __cdecl   FTS_MultiplySpectra(spectrum_t* spectrum1, const spectrum_t* spectrum2);
    int __stdcall FTS_MultiplySpectra_std(spectrum_t* spectrum1, const spectrum_t* spectrum2);

    /** This multiplies two spectra point wise
    *   On return will spectrum1 be equal to spectrum1 * spectrum2
    *   @param[in,out] spectrum1 The first factor and also the output spectrum
    *   @param[in] spectrum2 The second factor.
    *   @param length The number of valid data points in spectrum1 and spectrum2
    *   @return #FTS_SUCCESS on success     */
    int __cdecl   FTS_MultiplySpectra_Array(float* spectrum1, const float* spectrum2, unsigned int length);
    int __stdcall FTS_MultiplySpectra_Array_std(float* spectrum1, const float* spectrum2, unsigned int length);

    /** This divides two spectra 
    *   On return will spectrum1 be equal to spectrum1 / spectrum2
    *   @param[in,out] spectrum1 The dividend and also the output spectrum
    *   @param[in] spectrum2 The divisor.
    *   @return #FTS_ERROR_DIMENSION_MISMATCH if the lengths of the two spectra are different
    *   @return #FTS_ERROR_DIMENSION_MISMATCH if the x-axis values of the two spectra differ
    *   @return #FTS_SUCCESS on success     */
    int __cdecl   FTS_DivideSpectra(spectrum_t* spectrum1, const spectrum_t* spectrum2);
    int __stdcall FTS_DivideSpectra_std(spectrum_t* spectrum1, const spectrum_t* spectrum2);

    /** This divides two spectra 
    *   On return will spectrum1 be equal to spectrum1 / spectrum2
    *   @param[in,out] spectrum1 The dividend and also the output spectrum
    *   @param[in] spectrum2 The divisor.
    *   @param length The number of valid data points in both spectrum1 and spectrum2.
    *   @return #FTS_SUCCESS on success     */    
    int __cdecl   FTS_DivideSpectra_Array(float* spectrum1, const float* spectrum2, unsigned int length);
    int __stdcall FTS_DivideSpectra_Array_std(float* spectrum1, const float* spectrum2, unsigned int length);

    /** This adds a scalar to a spectrum
    *   @param[in, out] spectrum The spectrum to add the scalar to.
    *   @param scalar The scalar to add.
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_AddScalar(spectrum_t* spectrum, float scalar);
    int __stdcall FTS_AddScalar_std(spectrum_t* spectrum, float scalar);

    /** This adds a scalar to a spectrum
    *   @param[in, out] spectrum The spectrum to add the scalar to.
    *   @param scalar The scalar to add.
    *   @param length The number of valid data points in spectrum.
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_AddScalar_Array(float* spectrum, float scalar, unsigned int length);
    int __stdcall FTS_AddScalar_Array_std(float* spectrum, float scalar, unsigned int length);

    /** This multiplies a spectrum with a scalar
    *   @param[in, out] spectrum The spectrum to multiply with the scalar.
    *   @param scalar The scalar value to multiply the spectrum with.
    *   @return #FTS_SUCCESS on success     */
    int __cdecl   FTS_MultiplyScalar(spectrum_t* spectrum, float scalar);
    int __stdcall FTS_MultiplyScalar_std(spectrum_t* spectrum, float scalar);

    /** This multiplies a spectrum with a scalar
    *   @param[in, out] spectrum The spectrum to multiply with the scalar.
    *   @param scalar The scalar value to multiply the spectrum with.
    *   @param length The number of valid data points in spectrum.
    *   @return #FTS_SUCCESS on success     */
    int __cdecl   FTS_MultiplyScalar_Array(float* spectrum, float scalar, unsigned int length);
    int __stdcall FTS_MultiplyScalar_Array_std(float* spectrum, float scalar, unsigned int length);

    /** This calculates each data point in spectrum1 as the minimum of 
    *   spectrum1 and the same data point in spectrum2.
    *   @param[in,out] spectrum1 Is on successful return the minimum of spectrum1 and spectrum2, point by point
    *   @param[in] spectrum2 The spectrum to compare with. 
    *   @return #FTS_SUCCESS on success 
    *   @return #FTS_ERROR_PARAMETER_ERROR if the lengths differ */
    int __cdecl   FTS_MinSpectrum(spectrum_t* spectrum1, const spectrum_t* spectrum2);
    int __stdcall FTS_MinSpectrum_std(spectrum_t* spectrum1, const spectrum_t* spectrum2);

    /** This calculates each data point in spectrum1 as the minimum of 
    *   spectrum1 and the same data point in spectrum2.
    *   @param[in,out] spectrum1 Is on successful return the minimum of spectrum1 and spectrum2, point by point
    *   @param[in] spectrum2 The spectrum to compare with. 
    *   @param[in] length The length of the data arrays spectrum1 and spectrum2.
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_MinSpectrum_Array(float* spectrum1, const float* spectrum2, unsigned int length);
    int __stdcall FTS_MinSpectrum_Array_std(float* spectrum1, const float* spectrum2, unsigned int length);

    /** This calculates each data point in spectrum1 as the maximum of 
    *   spectrum1 and the same data point in spectrum2.
    *   @param[in,out] spectrum1 Is on successful return the maximum of spectrum1 and spectrum2, point by point
    *   @param[in] spectrum2 The spectrum to compare with. 
    *   @return #FTS_SUCCESS on success 
    *   @return #FTS_ERROR_PARAMETER_ERROR if the lengths differ */
    int __cdecl   FTS_MaxSpectrum(spectrum_t* spectrum1, const spectrum_t* spectrum2);
    int __stdcall FTS_MaxSpectrum_std(spectrum_t* spectrum1, const spectrum_t* spectrum2);

    /** This calculates each data point in spectrum1 as the maximum of 
    *   spectrum1 and the same data point in spectrum2.
    *   @param[in,out] spectrum1 Is on successful return the maximum of spectrum1 and spectrum2, point by point
    *   @param[in] spectrum2 The spectrum to compare with. 
    *   @param[in] length The length of the data arrays spectrum1 and spectrum2.
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_MaxSpectrum_Array(float* spectrum1, const float* spectrum2, unsigned int length);
    int __stdcall FTS_MaxSpectrum_Array_std(float* spectrum1, const float* spectrum2, unsigned int length);

    /** This inverts a spectrum point by point in linear units.
    *   If the spectrum is in a linear y-axis unit then the spectrum will on return be equal to 1/spectrum.
    *   If the spectrum is in a logarithmic y-axis unit then the spectrum will on return be equal to (-spectrum).
    *   @return #FTS_SUCCESS on success     */
    int __cdecl   FTS_InvertSpectrum(spectrum_t* spectrum);
    int __stdcall FTS_InvertSpectrum_std(spectrum_t* spectrum);

    /** This inverts a spectrum point by point.
    *   On return will spectrum be equal to 1/spectrum
    *   It is here assumed that the spectrum is in linear y-axis units.
    *   @return #FTS_SUCCESS on success     */
    int __cdecl   FTS_InvertSpectrum_Array(float* spectrum1, unsigned int length);
    int __stdcall FTS_InvertSpectrum_Array_std(float* spectrum1, unsigned int length);

    /** This zeros-out a spectrum between the two given indices
    *   On return will spectrum be equal to zero at and between the two indices
    *   @return #FTS_SUCCESS on success     */
    int __cdecl   FTS_ZeroOutSpectrum(spectrum_t* spectrum, unsigned int fromIdx, unsigned int toIdx);
    int __stdcall FTS_ZeroOutSpectrum_std(spectrum_t* spectrum, unsigned int fromIdx, unsigned int toIdx);

    /** This zeros-out a spectrum between the two given indices
    *   On return will spectrum be equal to zero at and between the two indices
    *   @return #FTS_SUCCESS on success     */
    int __cdecl   FTS_ZeroOutSpectrum_Array(float* spectrum, unsigned int length, unsigned int fromIdx, unsigned int toIdx);
    int __stdcall FTS_ZeroOutSpectrum_Array_std(float* spectrum, unsigned int length, unsigned int fromIdx, unsigned int toIdx);

    /** This stitches together two spectra
    *   @param[in,out] destination the buffer to hold the output. This must be allocated to hold at least (src1->length + src2->length) elements
    *       for the I array and the x-array.
    *   @param[in] src1, src2 The two spectra to stitch together. These may or may not be overlapping.
    *   @param stitchOption Defines how the stitch will be performed in the overlapping region (if any).
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_StitchSpectra(spectrum_t* destination, const spectrum_t* src1, const spectrum_t* src2, STITCH_OPTION stitchOption);
    int __stdcall FTS_StitchSpectra_std(spectrum_t* destination, const spectrum_t* src1, const spectrum_t* src2, STITCH_OPTION stitchOption);

    /** Calculates the logarithm of the given spectrum 
    *   @param[in,out] spec The spectrum to calculate the logarithm of
    *   @param base The base to use. Must be larger than zero
    *   @return #FTS_SUCCESS upon success
    */
    int __cdecl   FTS_LogSpectrum(spectrum_t* spec, double base);    
    int __stdcall FTS_LogSpectrum_std(spectrum_t* spec, double base);    

    /** Calculates the logarithm of the given spectrum 
    *   @param[in,out] spec The spectrum to calculate the logarithm of
    *   @param[in] length The length of the spectrum.
    *   @param[in] base The base to use. Must be larger than zero
    *   @return #FTS_SUCCESS upon success
    */
    int __cdecl   FTS_LogSpectrum_Array(float* spec, unsigned int length, double base);    
    int __stdcall FTS_LogSpectrum_Array_std(float* spec, unsigned int length, double base);    

    /** Calculates the given spectrum raised to the given base
    *   @param[in,out] spec The spectrum to calculate the power of.
    *   @param base The base to use.
    *   @return #FTS_SUCCESS upon success
    */
    int __cdecl   FTS_PowSpectrum(spectrum_t* spec, double base);    
    int __stdcall FTS_PowSpectrum_std(spectrum_t* spec, double base);    

    /** Calculates the given spectrum raised to the given base
    *   @param[in,out] spec The spectrum to calculate the logarithm of
    *   @param[in] base The base to use. Must be larger than zero
    *   @param[in] length The length of the data arrays spectrum1 and spectrum2.
    *   @return #FTS_SUCCESS upon success
    */
    int __cdecl   FTS_PowSpectrum_Array(float* spec, unsigned int length, double base);    
    int __stdcall FTS_PowSpectrum_Array_std(float* spec, unsigned int length, double base);    

    /** Performs the math_operand operation on the two spectra. 
    *   @param[in,out] spectrum1 This is parameter one to the input and is also the output
    *   @param[in] spectrum2 This is parameter two to the input and will not be modified
    *   @param operand The operation to perform, this must be one of the operator_spectrum_... 
    *       operations as defined in FTSMath.h. 
    *   @return #FTS_SUCCESS if the operation succeeds.
    *   @return #FTS_ERROR_RANGE_ERROR if the two spectra are not of the same length
    *   @return #FTS_ERROR_PARAMETER_ERROR if the operand is of the wrong type */
    int __cdecl   FTS_SpectrumBinaryMath(spectrum_t* spectrum1, const spectrum_t* spectrum2, math_operator operand);
    int __stdcall FTS_SpectrumBinaryMath_std(spectrum_t* spectrum1, const spectrum_t* spectrum2, math_operator operand);

    /** Performs the math_operand operation on spectrum and a scalar 
    *   @param[in,out] spectrum This is the spectrum to operate on. This is input as well as output
    *   @param scalar The scalar to use in the operation.
    *   @param operand The operation to perform, this must be one of the operator_scalar_... 
    *   operations as defined in FTSMath.h.
    *   The following calculations are possible:
    *       <ul>
    *       <li>operator_scalar_addition                =>    spectrum = spectrum + scalar </li>
    *       <li>operator_scalar_subtraction                =>    spectrum = spectrum - scalar </li>
    *       <li>operator_scalar_subtraction_scalar_first=>    spectrum = scalar - spectrum </li>
    *       <li>operator_scalar_multiplication            =>    spectrum = spectrum * scalar </li>
    *       <li>operator_scalar_division                =>    spectrum = spectrum / scalar </li>
    *       <li>operator_scalar_division_scalar_first    =>    spectrum = scalar / spectrum </li>
    *       <li>operator_scalar_power                    =>    spectrum = spectrum ^scalar  </li>
    *       </ul>
    *   @return #FTS_SUCCESS if the operation succeeds.
    *   @return #FTS_ERROR_PARAMETER_ERROR if the operand is of the wrong type */
    int __cdecl   FTS_SpectrumScalarMath(spectrum_t* spectrum, float scalar, math_operator operand);
    int __stdcall FTS_SpectrumScalarMath_std(spectrum_t* spectrum, float scalar, math_operator operand);

    /** Performs the math_operand operation on the single spectrum
    *   @param[in,out] spectrum This is the spectrum to operate on. This is input as well as output
    *   @param operand The operation to perform, this must be one of the operator_unary_... 
    *   operations as defined in FTSMath.h. 
    *   @return #FTS_SUCCESS if the operation succeeds.
    *   @return #FTS_ERROR_PARAMETER_ERROR if the operand is of the wrong type */
    int __cdecl   FTS_SpectrumUnaryMath(spectrum_t* spectrum, math_operator operand);
    int __stdcall FTS_SpectrumUnaryMath_std(spectrum_t* spectrum, math_operator operand);

    /** This integrates the given spectrum between the x-values 'xMin' and 'xMax'
    *   @param[in] spectrum The spectrum to integrate
    *   @param[in] xMin, xMax Defines the range of x-values that will be used in the integration
    *   if these are 0 then the entire data range will be used. It is assumed that the function to
    *   integrate is zero outside the x-data range defined by xData.
    *   @param[out] result Will on successful return be filled with the calculated integral
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_IntegrateSpectrum(const spectrum_t* spectrum, float xMin, float xMax, double* result);
    int __stdcall FTS_IntegrateSpectrum_std(const spectrum_t* spectrum, float xMin, float xMax, double* result);

    /** This integrates the given spectrum (defined by the xData and yData arrays) between the x-values 'xMin' and 'xMax'
    *   @param[in] xData can be given in three different forms as specified in \ref sec_x_value_formats.
    *   @param[in] yData the spectral data to integrate
    *   @param[in] length the length of the yData array (and xData if xValueOption equals #FTS_OPTION_ARRAY_X)
    *   @param[in] xValueOption determines the options for how the x-data is interpreted, see \ref sec_x_value_formats.
    *   @param[in] xLimit1, xLimit2 Defines the range of x-values that will be used in the integration
    *   if these are 0 then the entire data range will be used. It is assumed that the function to
    *   integrate is zero outside the x-data range defined by xData.
    *   @param[out] result Will on successful return be filled with the calculated integral
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_IntegrateSpectrum_Array(const float* xData, const float* yData, unsigned int length, float xLimit1, float xLimit2, unsigned int xValueOption, double* result);
    int __stdcall FTS_IntegrateSpectrum_Array_std(const float* xData, const float* yData, unsigned int length, float xLimit1, float xLimit2, unsigned int xValueOption, double* result);

#define CONVOLVE_SAME (0)
#define CONVOLVE_FAST (1)
#define CONVOLVE_FULL (2)
    /** This performs a convolution between the two input arrays 'data' and 'kernel'. The result is stored in the provided 'outputBuffer'.
    *   @param[in] data the data to convolve
    *   @param[in] dataLength the length of the data
    *   @param[in] kernel the kernel to convolve with
    *   @param[in] kernelLength the length of the kernel.
    *   @param convolutionType the type of the convolution to perform. This can be either:
    *   <ul>
    *   <li>#CONVOLVE_SAME - the central part of the convolution will be returned and outputBuffer will need to be allocated with at least 'dataLength' number of values.</li>
    *   <li>#CONVOLVE_FAST - which will use an FFT to calculate the convolution, this will return the data set as CONVOLVE_SAME, and outputBuffer will need to be allocated with at least 'dataLength' number of values.</li>
    *   <li>#CONVOLVE_FULL - the full convolution will be returned, outputBuffer will need to be allocated with at least 'dataLength+kernelLength+1' number of values.</li>
    *   </ul>
    *   @param[out] outputBuffer will on successful return be filled with the convolution between data and the kernel. This needs to be allocated
    *   before calling this routine, the convolution type determines the necessary length of this buffer.  
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_ConvolveArrays(const float* data, unsigned int dataLength, const float* kernel, unsigned int kernelLength, int convolutionType, float* outputBuffer);
    int __stdcall FTS_ConvolveArrays_std(const float* data, unsigned int dataLength, const float* kernel, unsigned int kernelLength, int convolutionType, float* outputBuffer);

    // ----------------------------------------------------------------------------------------
    // ------------------------------ SPECTRUM PROCESSSING ------------------------------------
    // ----------------------------------------------------------------------------------------

    /** Calculates a transmittance spectrum from a sample and a background spectrum 
    *   The 'sample' and 'background' spectra must both be of the type 'SPEC_EMISSION'
    *   @param[in] sample The spectrum of the sample to calculate the transmission of
    *   @param[in] background The spectrum to use as background spectrum
    *   @param[out] transmittance Will on successful return be filled with the calculated transmission.
    *   @return #FTS_SUCCESS on success    */
    int __cdecl   FTS_CalculateTransmission_2(const spectrum_t* sample, const spectrum_t* background, spectrum_t* transmittance);
    int __stdcall FTS_CalculateTransmission_2_std(const spectrum_t* sample, const spectrum_t* background, spectrum_t* transmittance);

    /** Calculates a transmittance spectrum from an absorbance spectrum 
    *   The 'absorbance' spectrum must be of the type 'SPEC_ABSORBANCE'
    *   @param[in] absorbance An absorbance spectrum to convert to a transmission spectrum
    *   @param[out] transmittance Will on successful return be filled with the calculated transmission.
    *   @return #FTS_SUCCESS on success    */
    int __cdecl   FTS_CalculateTransmission_1(const spectrum_t* absorbance, spectrum_t* transmittance);
    int __stdcall FTS_CalculateTransmission_1_std(const spectrum_t* absorbance, spectrum_t* transmittance);

    /** Calculates a absorbance spectrum from a sample and a background spectrum 
    *   The 'sample' and 'background' spectra must both be of the type 'SPEC_EMISSION'
    *   @param[in] sample The spectrum of the sample to calculate the absorbance of
    *   @param[in] background The spectrum to use as background spectrum
    *   @param[out] absorbance Will on successful return be filled with the calculated absorbance.
    *   @return #FTS_SUCCESS on success    */
    int __cdecl   FTS_CalculateAbsorbance_2(const spectrum_t* sample, const spectrum_t* background, spectrum_t* absorbance);
    int __stdcall FTS_CalculateAbsorbance_2_std(const spectrum_t* sample, const spectrum_t* background, spectrum_t* absorbance);

    /** Calculates a absorbance spectrum from a transmittance spectrum 
    *   The 'transmittance' spectrum must  be of the type 'SPEC_TRANSMITTANCE'
    *   @param[in] transmittance The transmittance spectrum to convert
    *   @param[out] absorbance Will on successful return be filled with the calculated absorbance spectrum.
    *   @return #FTS_SUCCESS on success    */
    int __cdecl   FTS_CalculateAbsorbance_1(const spectrum_t* transmittance, spectrum_t* absorbance);
    int __stdcall FTS_CalculateAbsorbance_1_std(const spectrum_t* transmittance, spectrum_t* absorbance);

    /** Calculates a reflectance spectrum from a sample and a background spectrum 
    *   @param[in] sample The spectrum of the sample to calculate the reflectance of
    *   @param[in] background The spectrum to use as background spectrum
    *   @param[out] reflectance Will on successful return be filled with the calculated reflectance.
    *   The 'sample' and 'background' spectra must both be of the type 'SPEC_EMISSION'
    *   @return #FTS_SUCCESS on success    */
    int __cdecl   FTS_CalculateReflectance_2(const spectrum_t* sample, const spectrum_t* background, spectrum_t* reflectance);
    int __stdcall FTS_CalculateReflectance_2_std(const spectrum_t* sample, const spectrum_t* background, spectrum_t* reflectance);

    /** Calculates a reflectance spectrum from a Log(1/R) spectrum 
    *   The 'log1R' spectrum must  be of the type 'SPEC_LOG_1_R'
    *   @param[in] log1R The Log(1/R) spectrum to convert
    *   @param[out] reflectance Will on successful return be filled with the calculated reflectance spectrum.
    *   @return #FTS_SUCCESS on success    */
    int __cdecl   FTS_CalculateReflectance_1(const spectrum_t* log1R, spectrum_t* reflectance);
    int __stdcall FTS_CalculateReflectance_1_std(const spectrum_t* log1R, spectrum_t* reflectance);

    /** Calculates a Log(1/R) spectrum from a sample and a background spectrum 
    *   The 'sample' and 'background' spectra must both be of the type 'SPEC_EMISSION'
    *   @param[in] sample The spectrum of the sample to calculate the reflectance of
    *   @param[in] background The spectrum to use as background spectrum
    *   @param[out] log1R Will on successful return be filled with the calculated reflectance.
    *   @return #FTS_SUCCESS on success    */
    int __cdecl   FTS_CalculateLog1R_2(const spectrum_t* sample, const spectrum_t* background, spectrum_t* log1R);
    int __stdcall FTS_CalculateLog1R_2_std(const spectrum_t* sample, const spectrum_t* background, spectrum_t* log1R);

    /** Calculates a Log(1/R) spectrum from a reflectance spectrum 
    *   The 'reflectance' spectrum must be of the type 'SPEC_REFLECTANCE'
    *   @param[in] reflectance The reflectance spectrum to convert
    *   @param[out] log1R Will on successful return be filled with the calculated Log(1/R) spectrum.
    *   @return #FTS_SUCCESS on success    */
    int __cdecl   FTS_CalculateLog1R_1(const spectrum_t* reflectance, spectrum_t* log1R);        
    int __stdcall FTS_CalculateLog1R_1_std(const spectrum_t* reflectance, spectrum_t* log1R);        

    /** Calculates the envelope of the given spectrum/interferogram.
    *   @param[in] source The spectrum / interferogram to calculate the envelope of 
    *   @param[out] envelope Will on successful return be filled with the calculated envelope.
    *   If this is not allocated large enough then this will be re-allocated during calculation. 
    *   @return #FTS_SUCCESS if the operation succeeds */
    int __cdecl   FTS_CalculateEnvelope(const spectrum_t* source, spectrum_t* envelope);
    int __stdcall FTS_CalculateEnvelope_std(const spectrum_t* source, spectrum_t* envelope);

    /** Resamples the given spectrum/interferogram with the desired factor
    *   The length of the spectrum/interferogram will be (current length) * factor
    *   @param[in,out] spectrum The spectrum / interferogram to resample.
    *       If the spectrum is to be extended then it will be re-allocated if not enough
    *       memory has been allocated.
    *   @param[in] factor Specifies how much to extend / shorten the spectrum with 
    *       factor > 1.0 results in an extension of the spectrum
    *       factor < 1.0 results in a shortering of the spectrum
    *   @return #FTS_SUCCESS if the operation succeeds
    *   @return #FTS_ERROR_PARAMETER_ERROR if factor <= 0.0 or the spectrum has length zero. */
    int __cdecl   FTS_ResampleSpectrum(spectrum_t* spectrum, double factor);
    int __stdcall FTS_ResampleSpectrum_std(spectrum_t* spectrum, double factor);

    /** Resamples the given spectrum/interferogram by simply picking data points starting at 'startIndex'
    *   and then every 'stepSize' data point. The result will contain the data points: startIndex + N*stepSize.
    *   @param[in,out] spectrum The spectrum / interferogram to resample.
    *   @param[in] factor Specifies how much to extend / shorten the spectrum with 
    *       factor > 1.0 results in an extension of the spectrum
    *       factor < 1.0 results in a shortering of the spectrum
    *   @return #FTS_SUCCESS if the operation succeeds
    *   @return #FTS_ERROR_PARAMETER_ERROR if startIndex < 0 or the spectrum has length zero. */
    int __cdecl   FTS_ResampleSpectrum_Index(spectrum_t* spectrum, int startIndex, int stepSize);
    int __stdcall FTS_ResampleSpectrum_Index_std(spectrum_t* spectrum, int startIndex, int stepSize);

    /** Resamples the given spectrum/interferogram to the provided new x-axis values.
    *   The length of the spectrum/interferogram will be equal newXAxisLength.
    *   @param[in,out] spectrum The spectrum / interferogram to resample.
    *       If the spectrum is to be extended then it will be re-allocated if not enough
    *       memory has been allocated.
    *   @param[in] newXAxisValues specifies the x-axis values to where the spectrum
    *       should be resampled. It is assumed that the spectrum equals zero outside of the
    *       current xMin to xMax range, i.e. if any point in newXAxisValues lies outside of the
    *       current x-axis range of the input spectrum, then the resulting spectrum will be zero there.
    *   @param[in] newXAxisLength the length of the array 'newXAxisValues'.
    *   @return #FTS_SUCCESS if the operation succeeds
    *   @return #FTS_ERROR_PARAMETER_ERROR if newXAxisValues is null or newXAxisLength is zero. */
    int __cdecl   FTS_ResampleSpectrum_xArray(spectrum_t* spectrum, const double* newXAxisValues, unsigned int newXAxisLength);
    int __stdcall FTS_ResampleSpectrum_xArray_std(spectrum_t* spectrum, const double* newXAxisValues, unsigned int newXAxisLength);

    /** Resamples the given spectrum/interferogram to the provided new x-axis values.
    *   The length of the spectrum/interferogram will be equal newXAxisLength.
    *   The output spectrum will have it's x-axis values at equally spaced data points from newXMin to newXMax.
    *   @param[in,out] spectrum The spectrum / interferogram to resample.
    *       If the spectrum is to be extended then it will be re-allocated if not enough
    *       memory has been allocated.
    *   @param newXMin specifies the smallest x-axis value of the resampled spectrum.
    *   @param newXMax specifies the largest x-axis value of the resampled spectrum.
    *   @param newXAxisLength specifies the length of the output spectrum. 
    *   @return #FTS_SUCCESS if the operation succeeds
    *   @return #FTS_ERROR_PARAMETER_ERROR if the spectrum has length zero, newXAxisLength is zero or the interval newXMin to newXMax is not contained
    *       in the range defined by the existing spectrum x-axis data. */
    int __cdecl   FTS_ResampleSpectrum_xMinMax(spectrum_t* spectrum, double newXMin, double newXMax, unsigned int newXAxisLength);
    int __stdcall FTS_ResampleSpectrum_xMinMax_std(spectrum_t* spectrum, double newXMin, double newXMax, unsigned int newXAxisLength);

    ///@}

    int __cdecl FTS_CalculateInterferogram(const spectrum_t* spectrum, spectrum_t* interferogramBuffer, float phaseError_0, float phaseError_1, float phaseError_2);

    /** Calculates the Hilbert-transform of the provided spectrum. The result of the calculation is a complex-valued spectrum and the provided spectrum
    *   will NOT be reallocated so the argument needs to be allocated with a 'phi' array. 
    *   @param[in, out] spectrum the spectrum to calculate the Hilbert transform of. This may be a complex-valued spectrum (see \ref sec_complex_spectral_data)
    *   @param firstIndexToKeep the first index in the calculated spectrum which will _not_ be set to zero.
    *   @param lastIndexToKeep the last index in the calculated spectrum which will _not_ be set to zero.
    *   @return #FTS_ERROR_OUTOFMEMORY if the spectrum does not have an 'phi' array allocated, or the size of the 'phi' array is smaller than the length of the spectrum.
    *   @return #FTS_ERROR_NOT_POWEROFTWO if the length of the spectrum is not an allowed length for the FFT calculations.
    *   @return #FTS_SUCCESS on successful calculation of the transform. */
    int __cdecl FTS_CalculateHilbertTransform(spectrum_t* spectrum, unsigned int firstIndexToKeep, unsigned int lastIndexToKeep);
    int __stdcall FTS_CalculateHilbertTransform_std(spectrum_t* spectrum, unsigned int firstIndexToKeep, unsigned int lastIndexToKeep);

#ifdef __cplusplus
}
#endif 

#endif