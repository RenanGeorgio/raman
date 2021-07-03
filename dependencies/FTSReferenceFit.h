////////////////////////////////////////////////////////////////////////
// FTSReferenceFit.h
//    This is the handling of fitting references to measured spectra used by the
//    Thorlabs Fourier Transform Spectrometer DLL. 
//    
//    Copyright (c) 2019 Thorlabs Sweden, All Rights Reserved
////////////////////////////////////////////////////////////////////////

#include "FTSData.h"
#include "FTSMath.h"

#ifndef FTSREFERENCEFIT_H
#define FTSREFERENCEFIT_H

#ifdef __cplusplus
extern "C" {
#endif 

    enum FitType
    {
        FIT_LINEAR_LEAST_SQUARES = 0,    /** FIT_LINEAR_LEAST_SQUARES type only takes linear parameters into account */
        FIT_LEVENBERG_MARQUARD = 1        /** FIT_LEVENBERG_MARQUARD can also handle non-linear parameters, such as the shift */
    };

    /** A FTSREFERENCEFit structure holds the necessary information for performing a fit. */
    typedef struct FTSREFERENCEFIT
    {
        /** fitRegionLow and fitRegionHigh defines the x-axis range over which the fit is to be done.
        The values are in the unit fitRegionUnit. */
        double fitRegionLow;
        double fitRegionHigh;
        int    fitRegionUnit;

        /** fit is a pointer to the fit object. This can only be accessed through any of the
        routines in FTSReferenceFit.h */
        void * fit;

        /** result is a pointer to the result set. This can only be accessed through any of the
        routines in FTSReferenceFit.h */
        void * result;

        /** userMessage can be filled with a message for the user in case any of the routines, to which this is a parameter, should fail */
        char userMessage[512];
    } FTSREFERENCEFIT;

    enum PARAMETEROPTION
    {
        FIT_FREE = 0,
        FIT_FIX = 1
    };

    typedef struct REFERENCEOPTION
    {
        PARAMETEROPTION        coefficient;
        PARAMETEROPTION        shift;
    } REFERENCEOPTION;

    /** A FTSFITRESULT structure holds the result from a reference fit. */
    typedef struct FTSFITRESULT
    {
        /** coefficient is the scaling coefficient for each of the fitted references */
        double            coefficient[128];

        /** shift is the shift applied to for each of the fitted references */
        double            shift[128];

        /** size is the number of references fitted, this is the number valid values in fitCoefficient */
        unsigned int    size;
    } FTSFITRESULT;

    /** \addtogroup grp_reference_fit Reference Fit
    *   @{
    */

    /** FTS_CreateFit sets up a new #FTSREFERENCEFIT using the provided type of fit. 
    *   @return a newly created #FTSREFERENCEFIT. This must later be deleted using FTS_DeleteFit or a memory leak will occur. 
    *   @return nullptr if the type is invalid. */
    FTSREFERENCEFIT* __cdecl    FTS_CreateFit();
    FTSREFERENCEFIT* __stdcall    FTS_CreateFit_std();

    /** FTS_DeleteFit deletes a #FTSREFERENCEFIT previously created using FTS_CreateFit. 
    *   If the pointer is invalid or was not created using FTS_CreateFit, then undefined behaviors may occur */
    void __cdecl    FTS_DeleteFit(FTSREFERENCEFIT* );
    void __stdcall    FTS_DeleteFit_std(FTSREFERENCEFIT* );

    /** Sets the provided data as a the measured data of the fit, i.e. the data to which the references should be fitted. */
    int    __cdecl        FTS_SetMeasuredDataOfFit(FTSREFERENCEFIT* fit, const spectrum_t* const );
    int    __stdcall    FTS_SetMeasuredDataOfFit_std(FTSREFERENCEFIT* fit, const spectrum_t* const );

    /** Adds a polynomial of the provided order to the fit 
    *   @param[in,out] fit the fit to modify.
    *   @param polynomialOrder the desired order of the added polynomial. Must be >= 0 */
    int __cdecl        FTS_AddPolynomialToFit(FTSREFERENCEFIT* fit, int polynomialOrder);
    int __stdcall    FTS_AddPolynomialToFit_std(FTSREFERENCEFIT* fit, int polynomialOrder);

    /** Adds the provided data as a reference to the fit.
    *   @param[in,out] fit the fit to modify.
    *   @param[in] option the options for how the reference should be handled in the fit.
    *   @param[in] spectrum the reference to add to the fit. */
    int    __cdecl        FTS_AddReferenceToFit(FTSREFERENCEFIT* fit, REFERENCEOPTION* option, const spectrum_t* const spectrum);
    int    __stdcall    FTS_AddReferenceToFit_std(FTSREFERENCEFIT* fit, REFERENCEOPTION* option, const spectrum_t* const spectrum);

    /** Adds the provided data as a reference to the fit.
    *   @param[in, out] fit the fit to update.
    *   @param[in] option the options for how the reference should be handled in the fit.
    *   @param[in] xData can be given in two different forms:
    *       1) as an array of length equal to 'yData', each x-data point is then assumed
    *           to correspond to the yData point at the same index
    *       2) as an array of two values, these then represent the smallest and the largest x-value
    *           for the data in yData between startIdx and stopIdx. The actual x-values will be linearly interpolated 
    *       3) as a NULL pointer The xAxis values are then the index values of the reference data.
    *   @param[in] yData The values of the reference. 
    *   @param length The length of the input data series
    *   @param xAxisUnit The unit for the x-axis data.
    *   @param yAxisUnit The unit for the y-axis data.
    *   @param xValueOption determines the options for how the x-data is interpreted
    *       this is either of the following values;
    *       #FTS_OPTION_ARRAY_X  <-> x-values in array form
    *       #FTS_OPTION_MINMAX_X <-> x-values in min-max form (only first two values in 'xData' are used)
    */
    int    __cdecl        FTS_AddReferenceToFit_Array(FTSREFERENCEFIT* fit, REFERENCEOPTION* option, const float* const xData, const float* const yData, unsigned int length, int xAxisUnit, int yAxisUnit, unsigned int xValueOption);
    int    __stdcall    FTS_AddReferenceToFit_Array_std(FTSREFERENCEFIT* fit, REFERENCEOPTION* option, const float* const xData, const float* const yData, unsigned int length, int xAxisUnit, int yAxisUnit, unsigned int xValueOption);

    /** Sets the x-axis region used to perform the fit */
    int __cdecl        FTS_UpdateFitRange(FTSREFERENCEFIT*, double low, double high, int unit);
    int __stdcall    FTS_UpdateFitRange_std(FTSREFERENCEFIT*, double low, double high, int unit);


    /** Performs the fit.
    *   @return #FTS_SUCCESS if the fit succeeds.
    *   @return #FTS_ERROR_DIMENSION_MISMATCH if the fit failed because one or more of the references doesn't overlap the measured spectrum.
    *   @return #FTS_ERROR_DIMENSION_MISMATCH if the units references doesn't match the unit of the measured spectrum.
    *   @return #FTS_ERROR_CONFIGURATION_ERROR if the fit failed because of errors in the setup
    *   @return #FTS_ERROR_FIT_FAILED if the fit failed for some other reason. 
    *   If the fit does not succeed then the reason why is presented in fit->userMessage.
    *   On success will the provided 'result' parameter be filled with the calculated results
    */
    int __cdecl        FTS_DoFit(FTSREFERENCEFIT* fit, FitType type, FTSFITRESULT* result);
    int __stdcall    FTS_DoFit_std(FTSREFERENCEFIT* fit, FitType type, FTSFITRESULT* result);


    /** Retrieves the scaled (and shifted) reference number 'referenceIdx' from the provided fit.
    *   @param[in] fit the fit to retrieve the result from
    *   @param referenceIdx the index of the reference to retrieve
    *   @param[in,out] destination will on successful return be filled with the scaled and shifted reference. 
    *       This may be null if only the required size is to be requested.
    *   @param[out] length if not null then this will be filled with the required size for the 'destination' spectrum.
    *   @return #FTS_ERROR_PARAMETER_ERROR if any of the parameters are invalid.
    *   @return #FTS_ERROR_OUTOFMEMORY if the destination is too small to be able to hold the result. */
    int __cdecl        FTS_GetFittedReference(FTSREFERENCEFIT* fit, unsigned int referenceIdx, spectrum_t* destination, unsigned int* length);
    int __stdcall    FTS_GetFittedReference_std(FTSREFERENCEFIT* fit, unsigned int referenceIdx, spectrum_t* destination, unsigned int* length);

    typedef struct spectrum_line_t
    {
        /** this is the wavenumber for the line, [cm^-1] */
        double wavenumber;

        /** this is the strength of each line in wavenumbers / (molecule * cm^2) @ 296K */
        double lineStrength;

        /** this is the energy level of the lower state, in cm^-1 */
        double lowerEnergyLevel;

        /** this is the air broadening coefficient (used to calculate the pressure broadening) 
        in wavenumbers/atm @ 296K */
        double airBroadenedHalfWidth;

        /** this is the self broadening coefficient (used to calculate the pressure broadening) 
        in wavenumbers/atm @ 296K */
        double selfBroadenedHalfWidth;

        /** this is the coefficient for the temperature broadening of the line. No unit. */
        double airBroadeningCoefficient;

        /** m_airBroadenedPressureShift is the shift in wavenumber due to pressure broadening 
        in wavenumbers/atm @ 296K */
        double airBroadenedPressureShift;

    } spectrum_line_t;

    typedef struct molecule_state_t
    {
        double    temperature_K;
        double    pressure_atm;
        double    partialPressure_atm;
        double    molecularWeight_g_per_mol;
    } molecule_state_t;

    /** Calculates an absorption cross section from the provided line-by-line description.
    *   The absorption cross section will be calculated using a Voigt profile, calculated on wnrNum equally spaced
    *   x-axis points between wnrMin and wnrMax
    *   @param[in,out] destination will on successful return be filled with the calculated absorption cross section
    *   @param[in] lines points to an array with the spectrum line descriptions
    *   @param lineNum the length of the array pointed to by the parameter lines
    *   @param[in] state Description of the molecule and the environment parameters.
    *   @param wnrMin the smallest wavenumber in the resulting absorption cross section
    *   @param wnrMax the largest wavenumber in the resulting absorption cross section
    *   @param wnrNum the number of data points in the destination.
    *   @param fCancel an optional boolean flag which makes it possible to cancel the calculation. If this is not null and set to true then the operation will cancel. */
    int __cdecl      FTS_CalculateAbsorptionCrossSection(spectrum_t* destination, const spectrum_line_t* const lines, unsigned int lineNum, const molecule_state_t* const state, double wnrMin, double wnrMax, unsigned int wnrNum, bool *fCancel, void (__cdecl *pCallbackFunction)(unsigned int /*identifier*/, int /*status*/, double /*percent_progress*/, void * /*data*/, const char errorMessages[512]));
    int __stdcall    FTS_CalculateAbsorptionCrossSection_std(spectrum_t* destination, const spectrum_line_t* const lines, unsigned int lineNum, const molecule_state_t* const state, double wnrMin, double wnrMax, unsigned int wnrNum, bool *fCancel, void (__stdcall *pCallbackFunction)(unsigned int /*identifier*/, int /*status*/, double /*percent_progress*/, void * /*data*/, const char errorMessages[512]) );

    ///@}

#ifdef __cplusplus
}
#endif 

#endif