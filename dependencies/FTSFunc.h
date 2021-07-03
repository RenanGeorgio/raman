////////////////////////////////////////////////////////////////////////
// FTSFunc.h
//	This is the handling of mathematical functions used by
//	Thorlabs Fourier Transform Spectrometer DLL. 
//	
//	Copyright (c) 2019 Thorlabs Sweden, All Rights Reserved
////////////////////////////////////////////////////////////////////////

#include "FTSData.h"
#include "FTSMath.h"

#ifndef FTSFUNC_H
#define FTSFUNC_H

#ifdef __cplusplus
extern "C" {
#endif 

    /** \addtogroup grp_function_math Function Math
    *   @{
    */

#define FTS_FUNC_GAUSSIAN	(1)
#define FTS_FUNC_LORENTZIAN (2)
#define FTS_FUNC_POLYNOMIAL (3)
#define FTS_FUNC_VOIGT      (4)
#define FTS_FUNC_SPLINE		(5)
#define FTS_FUNC_COMPOSITE  (6)
#define FTS_FUNC_UNDEFINED	(7)

    enum SplineType
    {
        SPLINE_CUBIC_HERMITE,
        SPLINE_QUINTIC_HERMITE
    };

    /** A FunctionParameter structure holds information on a single parameter in a function */
    typedef struct FTSFunctionParameter
    {
        FTSFunctionParameter();
        FTSFunctionParameter(char* setName, double value, bool isFixed, bool linear);

        /** This attempts to set the the value of this FTSFunctionParameter, but will
        *   limit the value to th LowerLimit and UpperLimit values */
        void SetValue(double newValue);

        /** The name of this parameter */
        char name[16];

        /** The value of this parameter */
        double Value;

        /** The standard error of this parameter */
        double Error;

        /** The minimum allowed value for this parameter */
        double LowerLimit;

        /** The maximum allowed value for this parameter */
        double UpperLimit;

        /** this is true if this parameter is to be fixed during the fit process */
        bool fixedValue;

        /** this is set to true if the value of the function is linear with respect to this parameter. 
        *   (don't ever change this) */
        bool linear;

    } FTSFunctionParameter;

    /** A FTSFUNC structure holds a function in the internal FTSLib format. */
    typedef struct FTSFUNC
    {
        /** type is the type of function that this represents. Must equal one of the 
        *   FTS_FUNC_... constants defined in FTSFunc.h */
        unsigned int    type;

        /** func is a pointer to the function. This can only be accessed through any of the
        *   routines in FTSFunc.h */
        void*           func;

        /** This is the list of parameters that this function accepts. This is an
        *   array of pointers to 'FunctionParameter' objects. This must only be allocated/resized/deallocated
        *   through a function call in FTSFunc.h. */
        FTSFunctionParameter**  parameters;

        /** The number of parameters found in the array 'parameters' */
        unsigned int    numberOfParameters;
    } FTSFUNC;

    /** Creates and returns a Gaussian function with the specified amplitude, mean, sigma and baseline level */
    FTSFUNC * __cdecl   FTS_CreateGaussian(double Amplitude, double Mean, double Sigma, double BaseLine);
    FTSFUNC * __stdcall FTS_CreateGaussian_std(double Amplitude, double Mean, double Sigma, double BaseLine);

    /** Creates and returns a Lorentzian function with the desired amplitude, mean, fwhm and baseline level */
    FTSFUNC * __cdecl   FTS_CreateLorentzian(double Amplitude, double Mean, double FWHM, double BaseLine);
    FTSFUNC * __stdcall FTS_CreateLorentzian_std(double Amplitude, double Mean, double FWHM, double BaseLine);

    /** Creates and returns a Voigt function with the specified parameters. 
    *   @param Amplitude the amplitude of the function.
    *   @param Mean the center position of the function.
    *   @param Gaussian_HWHM the Half-Width at Half-Maximum of the Gaussian Component.
    *   @param Lorentzian_HWHM the Half-Width at Half-Maximum of the Lorentzian Component. 
    *   @param BaseLine the baseline leve, this is added to the rest of the Voigt function and can serve as an offset. */
    FTSFUNC * __cdecl   FTS_CreateVoigtFunction(double Amplitude, double Mean, double Gaussian_HWHM, double Lorentzian_HWHM, double BaseLine);
    FTSFUNC * __stdcall FTS_CreateVoigtFunction_std(double Amplitude, double Mean, double Gaussian_HWHM, double Lorentzian_HWHM, double BaseLine);

    /** Creates and returns a Polynomial function with the desired parameters 
    *   @param order the order of the polynomial to create
    *   @param[in] parameters the parameters for the polynomial. There must be (order + 1) elements in this array.
    *       The parameters are stored with the zero:th order coefficient first, followed
    *       by the first order parameter and so forth.  */
    FTSFUNC * __cdecl   FTS_CreatePolynomial(unsigned int order, const double* parameters);
    FTSFUNC * __stdcall FTS_CreatePolynomial_std(unsigned int order, const double* parameters);

    /** Adds the two functions together into one combined function which equals term1 + term2.
    *   @return the newly created function. This needs to be deallocated using #FTS_DeleteFunc once it is no longer needed.
    *   @return null if any of the two input functions is null. */
    FTSFUNC * __cdecl   FTS_AddFunctions(FTSFUNC* term1, FTSFUNC* term2);
    FTSFUNC * __stdcall FTS_AddFunctions_std(FTSFUNC* term1, FTSFUNC* term2);

    /** Creates and returns a Spline function with the desired parameters. The function limits
    *   will automatically be set to the smallest and largest of the xValues. 
    *   @param type The type of spline to create.
    *   @param[in] xValues, yValues The values which the spline will interpolate between. The xValues
    *       must be sorted in increasing order. 
    *   @param valueNum The number of values in the arrays xValues and yValues. 
    *   @return the newly created spline function or NULL if the creation failed. */
    FTSFUNC * __cdecl   FTS_CreateSpline(SplineType type, const float* xValues, const float* yValues, unsigned int valueNum);
    FTSFUNC * __stdcall FTS_CreateSpline_std(SplineType type, const float* xValues, const float* yValues, unsigned int valueNum);

    /** This sets the range of x-value for which the given function is defined.
    *   Any point outside of the definition range has a value of zero. 
    *   @param[in,out] f The function to modify
    *   @param xMin The new smallest x-axis value, must be smaller than xMax
    *   @param xMax The new largest x-axis value, must be larger than xMin. */
    void __cdecl   FTS_SetFunctionLimits(FTSFUNC *f, double xMin, double xMax);
    void __stdcall FTS_SetFunctionLimits_std(FTSFUNC *f, double xMin, double xMax);

    /** Sets the x-axis unit of the given function. No conversion will be done here,
    *   this only sets a flag in the FTSFUNC.
    *   @param[in,out] f The function to modify
    *   @param newUnit The new x-axis unit. */
    void __cdecl   FTS_SetFunctionXAxisUnit(FTSFUNC *f, X_AXIS_UNIT newUnit);
    void __stdcall FTS_SetFunctionXAxisUnit_std(FTSFUNC *f, X_AXIS_UNIT newUnit);

    /** Sets the y-axis unit of the given function. No conversion will be done here,
    *   this only sets a flag in the FTSFUNC
    *   @param[in,out] f The function to modify
    *   @param newUnit The new x-axis unit. */
    void __cdecl   FTS_SetFunctionYAxisUnit(FTSFUNC *f, Y_AXIS_UNIT newUnit);
    void __stdcall FTS_SetFunctionYAxisUnit_std(FTSFUNC *f, Y_AXIS_UNIT newUnit);

    /** Retrieves the minimum x-value for which this function is defined.
    *   @param[in] func The function to query. */
    double __cdecl   FTS_GetXMin(const FTSFUNC *func);
    double __stdcall FTS_GetXMin_std(const FTSFUNC *func);

    /** Retrieves the maximum x-value for which this function is defined
    *   @param[in] func The function to query. */
    double __cdecl   FTS_GetXMax(const FTSFUNC *func);
    double __stdcall FTS_GetXMax_std(const FTSFUNC *func);

    /** Retrieves the x axis unit for this function
    *   @param[in] func The function to query. */
    X_AXIS_UNIT __cdecl   FTS_GetXAxisUnit(const FTSFUNC *func);
    X_AXIS_UNIT __stdcall FTS_GetXAxisUnit_std(const FTSFUNC *func);

    /** Retrieves the y axis unit for this function
    *   @param[in] func The function to query. */
    Y_AXIS_UNIT __cdecl   FTS_GetYAxisUnit(const FTSFUNC *func);
    Y_AXIS_UNIT __stdcall FTS_GetYAxisUnit_std(const FTSFUNC *func);

    /** Retrieves the the value of this function at the desired x value.
    *   @param[in] func The function to query. 
    *   @param xValue The x-axis value where to calculate the functions value. If this is outside of the functions definition range then zero is returned. */
    double __cdecl   FTS_GetValueAt(const FTSFUNC *func, double xValue);
    double __stdcall FTS_GetValueAt_std(const FTSFUNC *func, double xValue);

    /** Performs the math_operand operation on the spectrum and the function
    *   @param[in,out] spectrum this is parameter one to the input and is also the output
    *   @param[in] func this is parameter two to the input and will not be modified
    *   @param[in] operand the operation to perform, this must be one of the operator_spectrum_... 
    *       operations as defined in FTSMath.h. 
    *   @return #FTS_SUCCESS if the operation succeeds.
    *   @return #FTS_ERROR_RANGE_ERROR if the x-axis range of the spectrum does not overlap the	
    *       definition range of the function.
    *   @return #FTS_ERROR_PARAMETER_ERROR if the operand is of the wrong type */
    int __cdecl   FTS_SpectrumFunctionMath(spectrum_t* spectrum, const FTSFUNC *func, math_operator operand);
    int __stdcall FTS_SpectrumFunctionMath_std(spectrum_t* spectrum, const FTSFUNC *func, math_operator operand);

    /** Deletes the function and frees up the memory it uses. The pointer should not be
    *   dereferenced after this call! */
    void __cdecl   FTS_DeleteFunc(FTSFUNC *);
    void __stdcall FTS_DeleteFunc_std(FTSFUNC *);

    typedef struct FTSFunctionFitResult
    {
        // This is set to #FTS_SUCCESS if the fit succeeds, else set to one of the FTS_ERROR_... parameters
        int returnCode;

        // The produced function, in case the fit fails then this is null. 
        //  This needs to be deleted using #FTS_DeleteFunc to avoid memory leaks.
        FTSFUNC*  function;

        /** This is the number of iterations required to produce the fit */
        unsigned int iterationCount;

        /** This is the sum of squared differences between the data and the fitted function in all data points */
        double residualSumOfSquares;

        /** The R^2 (R-Squared) of the fit. */
        double RSquare;

        /** The Degrees of Freedom Adjusted R-Square */
        double AdjustedRSquare;

        /** The degrees of freedom, this equals the number of data points in the input data minus the 
        *   number of fitted coefficients */
        unsigned int residualDegreesOfFreedom;

        /** */
        double RootMeanSquareError;

        /** userMessage can be filled with a message for the user in case the fit fails. */
        char userMessage[512];
    } FTSFunctionFitResult;

    /** Deletes the FTSFunctionFitResult and frees up the memory it uses. 
    *   The pointer should not be dereferenced after this call! */
    void __cdecl   FTS_DeleteResult(FTSFunctionFitResult *);
    void __stdcall FTS_DeleteResult_std(FTSFunctionFitResult *);

    /** Fits a function to the given data set.
    *   @param xData the x-axis component of the data. See \ref sec_x_value_formats for allowed formats.
    *   @param yData the y-axis component of the data. This must be an array of at least length values.
    *   @param startIdx the first index in the xData and yData to fit the data to. Make sure that this is strictly smaller than stopIdx
    *   @param stopIdx the last index in the xData and yData to fit to. Make sure that this is smaller than or equal to (length - 1).
    *   @param xValueOption the option for how to interpret xData, must be either #FTS_OPTION_MINMAX_X or #FTS_OPTION_ARRAY_X, see \ref sec_x_value_formats.
    *   @param functionTemplate describes the function to use, including which parameters should be fitted and which should not.
    *       The set values of the functionPattern will be used as the initialValues for the search for the best fit..
    *   @return a FTSFunctionFitResult which contains the result of the fit. This is a pointer and must be deleted using #FTS_DeleteResult. */
    FTSFunctionFitResult* __cdecl   FTS_FitFunctionToData(const float* xData, const float* yData, unsigned int length, unsigned int startIdx, unsigned int stopIdx, unsigned int xValueOption, FTSFUNC* functionTemplate);
    FTSFunctionFitResult* __stdcall FTS_FitFunctionToData_std(const float* xData, const float* yData, unsigned int length, unsigned int startIdx, unsigned int stopIdx, unsigned int xValueOption, FTSFUNC* functionTemplate);

    ///@}

#ifdef __cplusplus
}
#endif 

#endif