////////////////////////////////////////////////////////////////////////
// FTSpectrometer.h
//	This is the the methods for retrieving data from the 
//	Thorlabs Fourier Transform Spectrometer 
//	
//	Copyright (c) 2019 Thorlabs Sweden, All Rights Reserved
////////////////////////////////////////////////////////////////////////

#ifndef _FTSPECTROMETER_H
#define _FTSPECTROMETER_H

#ifdef __cplusplus
extern "C" {
#endif 

    /** \addtogroup grp_instrument_communication Instrument Communication
    *  @{
    */

    /** This finds all Thorlabs Fourier Transform Spectrometers connected to the 
    *   system and fills in the data about them in the FTSData structure.
    *   This resets all previously existing data in the FTSData structure
    *   This function should be called first to make sure the spectrometers are initialized properly. 
    *   If any spectrometer has been unplugged or plugged in then this routine should be called again 
    *   to update the list of instrument and make sure that the properties for each instrument is set correctly.
    *   This will fill in the fts_instrument_t data structures in the FTSData, if you are using multiple 
    *   spectrometer then please check this data structure to verify the order in which the spectrometers have been indexed.
    *   @return the number of spectrometers found */
    int __cdecl   FTS_InitializeSpectrometers();
    int __stdcall FTS_InitializeSpectrometers_std();

    /** This opens the connection to the spectrometer with the given specIndex
    *   Note that, when using several spectrometers, each spectrometer must be opened.
    *   @param specIndex the index of the spectrometer to open the connection to. See \ref sec_numerating_spectrometers 
    *   @return #FTS_SUCCESS on success
    */
    int __cdecl   FTS_OpenSpectrometer(SpectrometerIndex specIndex);
    int __stdcall FTS_OpenSpectrometer_std(SpectrometerIndex specIndex);

    /** This closes the connection to the spectrometer Note that, when using several spectrometers, each spectrometer must be closed.
    *   @param specIndex the index of the spectrometer to close the connection to. See \ref sec_numerating_spectrometers 
    *   @return #FTS_SUCCESS on success
    */
    int __cdecl   FTS_CloseSpectrometer(SpectrometerIndex specIndex);
    int __stdcall FTS_CloseSpectrometer_std(SpectrometerIndex specIndex);

    /** Sets the spectrometer with the given index to a low power state.
    *   This may mean turning off the motor and/or powering down other components (where supported) */
    int __cdecl   FTS_PowerDownSpectrometer(SpectrometerIndex specIndex);
    int __stdcall FTS_PowerDownSpectrometer_std(SpectrometerIndex specIndex);

    /** This retrieves some information from a connected spectrometer. 
    *   @param specIndex the index of the spectrometer to get info on. See \ref sec_numerating_spectrometers 
    *   @param[out] serial must be a pre-allocated array capable of holding at least #SERIALNUMBER_LENGTH characters
    *   @param[out] minWaveNr will on successful return be filled with the smallest wave number that this instrument can measure.
    *   @param[out] maxWaveNr will on successful return be filled with the largest wave number that this instrument can measure.
    *   @return #FTS_ERROR_SPECTROMETER_NOT_FOUND if no spectrometer exists, e.g. if #FTS_InitializeSpectrometers has not been called.
    *   @return #FTS_ERROR_PARAMETER_ERROR if specIndex is larger than or equal to the number of initialized spectrometer (unless this is zero). 
    *   @return #FTS_SUCCESS otherwise */
    int __cdecl   FTS_GetInfoOnSpectrometer(SpectrometerIndex specIndex, char* serial, float* minWaveNr, float* maxWaveNr);
    int __stdcall FTS_GetInfoOnSpectrometer_std(SpectrometerIndex specIndex, char* serial, float* minWaveNr, float* maxWaveNr);

    /** @param specIndex the index of the spectrometer to query. See \ref sec_numerating_spectrometers 
    *   @return #FTS_SUCCESS if the spectrometer has a built in pressure sensor
    *   @return #FTS_ERROR_NOT_IMPLEMENTED otherwise */
    int __cdecl   FTS_HasPressureSensor(SpectrometerIndex specIndex);
    int __stdcall FTS_HasPressureSensor_std(SpectrometerIndex specIndex);

    /** @param specIndex the index of the spectrometer to query. See \ref sec_numerating_spectrometers 
    *   @return #FTS_SUCCESS if the spectrometer has a built in temperature sensor
    *   @return #FTS_ERROR_NOT_IMPLEMENTED otherwise */
    int __cdecl   FTS_HasTemperatureSensor(SpectrometerIndex specIndex);
    int __stdcall FTS_HasTemperatureSensor_std(SpectrometerIndex specIndex);

    /** @param specIndex the index of the spectrometer to query. See \ref sec_numerating_spectrometers 
    *   @return #FTS_SUCCESS if the spectrometer has a built in humidity sensor
    *   @return #FTS_ERROR_NOT_IMPLEMENTED otherwise */
    int __cdecl   FTS_HasHumiditySensor(SpectrometerIndex specIndex);
    int __stdcall FTS_HasHumiditySensor_std(SpectrometerIndex specIndex);

    /** This retrieves the temperature and/or atmospheric pressure from the sensor in the instrument 
    *   @param specIndex the index of the spectrometer to query. See \ref sec_numerating_spectrometers 
    *   @param paramOption this is an OR:ed combination of #FTS_AIR_MEASURE_PRESS, #FTS_AIR_MEASURE_TEMP and #FTS_AIR_MEASURE_RELHUM
    *       determines which parameters to read
    *   @param[out] parameters will on successful return be filled with the read in values in the order
    *       pressure, temperature and relative humidity. If only one parameter is requested then only that parameter
    *       will be returned.
    *       NOTICE: The relative humidity can only be read from the instrument if #FTS_HasHumiditySensor returns #FTS_SUCCESS
    *   @return #FTS_SUCCESS on successful retrieval of the environmental parameters */
    int __cdecl   FTS_GetEnvironmentParameters(SpectrometerIndex specIndex, int paramOption, float* parameters);
    int __stdcall FTS_GetEnvironmentParameters_std(SpectrometerIndex specIndex, int paramOption, float* parameters);

    /** This the retrieves the (interpolated) power calibration coefficients at the given wavenumbers
    *   @param specIndex the index of the spectrometer to query. See \ref sec_numerating_spectrometers 
    *   @param[out] wavenumberArray an array with the wavenumbers where to retrieve the coefficients
    *   @param[out] cCoeff will on successful return be filled with the calibration coefficients
    *       This must be pre-allocated to be able to hold 'arrayLength' numbers
    *   @param arrayLength the number of data points in 'wavenumberArray' 
    *   @return #FTS_SUCCESS on success
    */
    int __cdecl   FTS_GetPowerCalibrationAt(SpectrometerIndex specIndex, float* wavenumberArray, float* cCoeff, unsigned int arrayLength);
    int __stdcall FTS_GetPowerCalibrationAt_std(SpectrometerIndex specIndex, float* wavenumberArray, float* cCoeff, unsigned int arrayLength);

    /** This will check the status of the spectrometer number 'specIndex'
    *       This automatically fills in the status of the instrument in
    *       #FTSData->fts_instrument[specIndex].communicationStatus
    *   @param specIndex the index of the spectrometer to query. See \ref sec_numerating_spectrometers 
    *   @return #FTS_SUCCESS if all is ok with the spectrometer
    *   @return #FTS_ERROR_PARAMETER_ERROR if the spectrometer has not been initialized (in #FTS_InitializeSpectrometers() )
    *   @return #FTS_ERROR_SPECTROMETER_NOT_FOUND if the spectrometer is missing from the hardware list
    *   @return #FTS_ERROR_COMMUNICATION_NOT_OPENED if the communication with the spectrometer has not been opened using FTS_Open 
    *   @return #FTS_ERROR_SPECTROMETER_FAILURE if the communication with the spectrometer fails for some other reason
    *   @return #FTS_ERROR_INSTRUMENT_STATUS if there is an internal error in the instrument
    *       (in that case is the error code stored in 'instrumentStatus') */
    int __cdecl   FTS_CheckSpectrometer(SpectrometerIndex specIndex);
    int __stdcall FTS_CheckSpectrometer_std(SpectrometerIndex specIndex);

    /** This reads the internal VSE error flag in the instrument. 
    *   @param specIndex the index of the spectrometer to query. See \ref sec_numerating_spectrometers 
    *   @param status will on successful return be filled with the currently set VSE status flag.
    *   @return #FTS_SUCCESS if the reading was successful.
    *   @return #FTS_ERROR_SPECTROMETER_FAILURE if this parameter is not available on the current spectrometer */
    unsigned int __cdecl   FTS_ReadVSEStatus(SpectrometerIndex specIndex, unsigned int* status);
    unsigned int __stdcall FTS_ReadVSEStatus_std(SpectrometerIndex specIndex, unsigned int* status);

    /** This reads the instrument status of the given spectrometer
    *   This automatically fills in the status of the instrument in
    *   #FTS_GetFTSData()->fts_instrument[specIndex].instrumentStatus
    *   (which can be read out using the routine #FTS_GetInstrumentProperty_InstrumentStatus)
    *   @param specIndex the index of the spectrometer to query. See \ref sec_numerating_spectrometers 
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_ReadInstrumentStatus(SpectrometerIndex specIndex);
    int __stdcall FTS_ReadInstrumentStatus_std(SpectrometerIndex specIndex);

    /** This reboots the given spectrometer. The handle must be opened prior to calling this function.
    *   This will automatically close the communication, call #FTS_OpenSpectrometer() to
    *   start the communication again
    *   @param specIndex the index of the spectrometer to query. See \ref sec_numerating_spectrometers 
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_RebootInstrument(SpectrometerIndex specIndex);
    int __stdcall FTS_RebootInstrument_std(SpectrometerIndex specIndex);

    /** @return true if the motor of the provided spectrometer is running.
    *   @param specIndex the index of the spectrometer to query. See \ref sec_numerating_spectrometers */
    bool __cdecl FTS_MotorRunning(SpectrometerIndex specIndex);
    bool __stdcall FTS_MotorRunning_std(SpectrometerIndex specIndex);

    /** This starts or stops the motor of the provided spectrometer.
    *   @param specIndex the index of the spectrometer to query. See \ref sec_numerating_spectrometers.
    *   @param start will start the motor if set to true, otherwise will the motor stop.
    *   @return #FTS_SUCCESS on success. */
    int __cdecl FTS_StartMotor(SpectrometerIndex specIndex, bool start);
    int __stdcall FTS_StartMotor_std(SpectrometerIndex specIndex, bool start);

    /** This tries to reprogram the spectrometer with the given firmware 
    *   @param specIndex the index of the spectrometer to reprogram. See \ref sec_numerating_spectrometers 
    *   @param[out] firmwareFile The full path and filename of the file containing the firmware. This must contain only ASCII characters. 
    */
    int __cdecl   FTS_UpgradeFirmware(SpectrometerIndex specIndex, char* firmwareFile);
    int __stdcall FTS_UpgradeFirmware_std(SpectrometerIndex specIndex, char* firmwareFile);

    ///@}

    /** \addtogroup grp_sync_acquisition Synchronous Data Acquisition
    See Also \ref sec_syncCollection
    *  @{
    */

    /** This retrieves a single interferogram from the spectrometer 
    *   The retrieved interferogram is stored in the supplied buffer 'interferogram'
    *   @param specIndex the index of the spectrometer to get the data from. See \ref sec_numerating_spectrometers 
    *   @param[out] interferogram Will on successful return be filled with the data of the collected
    *   interferogram. This will be the raw interferogram without any processing performed.
    *   The memory for this spectrum must be allocated with enough space prior to calling this function
    *   @return #FTS_SUCCESS on success
    *   @return #FTS_ERROR_TIMEOUT if a timeout occurs
    *   @return #FTS_INSTR_STATUS_REFERENCE_ERROR if there is a problem with the reference, this means that the
    *   retrieved interferogram is faulty and should be ignored
    *   @return #FTS_INSTR_STATUS_VSERROR, #FTS_INSTR_STATUS_BIAS_ERROR or #FTS_INSTR_STATUS_MOTOR_TIMEOUT if an internal error has occurred in the instrument. No data is returned.  */
    int __cdecl   FTS_RetrieveSingleInterferogram(SpectrometerIndex specIndex, spectrum_t* interferogram);
    int __stdcall FTS_RetrieveSingleInterferogram_std(SpectrometerIndex specIndex, spectrum_t* interferogram);

    ///@}

    /** \addtogroup grp_async_acquisition Asynchronous Data Acquisition
    See Also \ref sec_asyncCollection
    *  @{
    */

    /** This starts a continuous acquisition for the given spectrometer with the given callback function.
    *   This will clear the contents of the stored last spectrum and interferogram.
    *   @param[in] specIndex the index of the spectrometer to start acquiring from. See \ref sec_numerating_spectrometers 
    *   @param[in] pGUIFunction a callback which will be called each time an interferogram has been collected
    *   and each time a spectrum has been calculated. See \ref sec_callback_asynch_acquisition for reference to the 
    *   parameters passed to the callback. This may be NULL.
    */
    int __cdecl   FTS_StartContinuousAcquisition(SpectrometerIndex specIndex, void (__cdecl *pGUIFunction)(unsigned short, unsigned int, unsigned int));
    int __stdcall FTS_StartContinuousAcquisition_std(SpectrometerIndex specIndex, void (__stdcall *pGUIFunction)(unsigned short, unsigned int, unsigned int));

    /** This stops the continuous acquisition of data from the given spectrometer 
    *   @param specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers 
    */
    int __cdecl   FTS_StopContinuousAcquisition(SpectrometerIndex specIndex);
    int __stdcall FTS_StopContinuousAcquisition_std(SpectrometerIndex specIndex);

    /** This collects a single interferogram and spectrum asynchronously.
    *       This will clear the contents of the stored last spectrum and interferogram.
    *   @param[in] specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers.
    *   @param[in] pGUIFunction a callback which will be called when the interferogram has been collected
    *       and when the spectrum has been calculated. See \ref sec_callback_asynch_acquisition for reference to the 
    *       parameters passed to the callback. This may be NULL.
    *   @return #FTS_SUCCESS if all is ok.
    *   @return #FTS_ERROR_PARAMETER_ERROR if the specIndex is larger than the index of any connected instrument.
    */
    int __cdecl   FTS_AcquireSingleSpectrum(SpectrometerIndex specIndex, void (__cdecl *pGUIFunction)(unsigned short, unsigned int, unsigned int));
    int __stdcall FTS_AcquireSingleSpectrum_std(SpectrometerIndex specIndex, void (__stdcall *pGUIFunction)(unsigned short, unsigned int, unsigned int));

    /** This collects a given number of interferograms and spectra asynchronously. 
    *   When the given number has been collected then the acquisition will stop.
    *   This will clear the contents of the stored last spectrum and interferogram.
    *   @param specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers.
    *   @param numberToCollect the number of interferograms and spectra to collect.
    *   @param[in] pGUIFunction a callback which will be called when the interferogram has been collected
    *       and when the spectrum has been calculated. See \ref sec_callback_asynch_acquisition for reference to the 
    *       parameters passed to the callback. This may be NULL.
    *   @return #FTS_SUCCESS if all is ok.
    *   @return #FTS_ERROR_PARAMETER_ERROR if the specIndex is larger than the index of any connected instrument.
    */
    int __cdecl   FTS_AcquireSpectra(SpectrometerIndex specIndex, unsigned int numberToCollect, void (__cdecl *pGUIFunction)(unsigned short, unsigned int, unsigned int));
    int __stdcall FTS_AcquireSpectra_std(SpectrometerIndex specIndex, unsigned int numberToCollect, void (__stdcall *pGUIFunction)(unsigned short, unsigned int, unsigned int));

    /** @return true if the spectrometer with the given index is acquiring data.
    *   @return false if the spectrometer with the given index is not acquiring data _or_ the index is invalid. 
    *   Notice this will only return true if this process is collecting data from the given spectrometer, any other process collecting data will _NOT_ be seen. 
    */
    bool __cdecl   FTS_IsAcquiring(SpectrometerIndex specIndex);
    bool __stdcall FTS_IsAcquiring_std(SpectrometerIndex specIndex);

    /** Runs an auto setup procedure asynchronously.
    *   @param[in] specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers.
    *   @param[in] options specifies the options for the auto setup procedure. THIS IS NOT YET IMPLEMENTED AND SHOULD BE NULL.
    *   @param[in] fCancel if not set to NULL then this parameter can be used to cancel the auto setup routine once it is running. If the flag is set to 'true' then
    *       the operation will cancel. The callback will then be called with #FTS_CANCEL.
    *   @param[in] pGUIFunction a callback which will be called when the auto setup routine is finished or has failed. This may be NULL.
    *       The parameters passed to the routine are
    *       1) the index of the spectrometer
    *       2) the current status of the operation. Equal to #FTS_SUCCESS when the operation has finished successfully, #FTS_CANCEL if the operation was canceled,
    *           #FTS_PROGRESS while in progress and #FTS_FAIL if the setup failed.
    *       3) unused
    *       4) a string describing the current state of progress. This may be empty. 
    *   @return #FTS_SUCCESS if all is ok and the job was queued.
    *   @return #FTS_ERROR_PARAMETER_ERROR if any of the input parameters were illegal.
    */
    int __cdecl		FTS_AutoSetup(SpectrometerIndex specIndex, const void* options, bool *fCancel, void (__cdecl *pGUIFunction)(SpectrometerIndex specIndex, unsigned int statusFlag, unsigned int, const char message[512]));
    int __stdcall	FTS_AutoSetup_std(SpectrometerIndex specIndex, const void* options, bool *fCancel, void (__stdcall *pGUIFunction)(SpectrometerIndex specIndex, unsigned int statusFlag, unsigned int, const char message[512]));

    /** #FTS_ClearCallbacks will clear out all references to callbacks passed to any routine providing a callback.
    *   This is useful e.g. when unloading the dll with an acquisition still running */
    void __cdecl    FTS_ClearCallbacks();
    void __stdcall  FTS_ClearCallbacks_std();

    ///@}

    /** This writes the information we have on the spectrometer with the given spectrometer index to the given file 
    *   @param specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers 
    *   @param[in] fileName The full filename and path to the file where to store the instrument information.
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_SaveInstrumentInfo(SpectrometerIndex specIndex, const char* fileName);
    int __stdcall FTS_SaveInstrumentInfo_std(SpectrometerIndex specIndex, const char* fileName);

    /** This writes the currently used acquisition options for the given spectrometers to the given file
    *   @param specIndex the index of the spectrometer. See \ref sec_numerating_spectrometers 
    *   @param[in] fileName The full filename and path to the file where to store the acquisition options.
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_SaveAcquisitionOptions(SpectrometerIndex specIndex, const char* fileName);
    int __stdcall FTS_SaveAcquisitionOptions_std(SpectrometerIndex specIndex, const char* fileName);

    /** This writes the currently used acquisition options for the given spectrometers to the given file
    *   @param[out] option The acquisition options to write to file.
    *   @param[in] fileName The full filename and path to the file where to store the acquisition options.
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_SaveAcquisitionOptions_ext(const acquisition_option_t* option, const char* fileName);
    int __stdcall FTS_SaveAcquisitionOptions_ext_std(const acquisition_option_t* option, const char* fileName);

    /** This writes the information we have on the spectrometer with the given spectrometer index to the given file 
    *   @param specIndex the index of the spectrometer, loads the instrument info to the fts_instrument_t
    *       structure in FTSData belonging to this spectrometer index. See \ref sec_numerating_spectrometers.		
    *   @param[in] fileName The full filename and path to the file where to read the instrument information from.
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_LoadInstrumentInfo(SpectrometerIndex specIndex, const char* fileName);
    int __stdcall FTS_LoadInstrumentInfo_std(SpectrometerIndex specIndex, const char* fileName);

    /** This writes the currently used acquisition options for the given spectrometers to the given file
    *   @param specIndex the index of the spectrometer, loads the acquisition options to the acquisition_option_t 
    *       structure in FTSData belonging to this spectrometer index. See \ref sec_numerating_spectrometers.
    *   @param[in] fileName The full filename and path to the file where to read the acquisition options from.
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_LoadAcquisitionOptions(SpectrometerIndex specIndex, const char* fileName);
    int __stdcall FTS_LoadAcquisitionOptions_std(SpectrometerIndex specIndex, const char* fileName);

    /** This writes the currently used acquisition options for the given spectrometers to the given file
    *   @param[out] option Will on successful return be filled with the read in acquisition options.
    *   @param[in] fileName The full filename and path to the file where to read the acquisition options from.
    *   @return #FTS_SUCCESS on success */
    int __cdecl   FTS_LoadAcquisitionOptions_ext(acquisition_option_t* option, const char* fileName);
    int __stdcall FTS_LoadAcquisitionOptions_ext_std(acquisition_option_t* option, const char* fileName);

    ///@}

#ifdef __cplusplus
}
#endif 

#endif