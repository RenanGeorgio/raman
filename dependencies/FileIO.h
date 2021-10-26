////////////////////////////////////////////////////////////////////////
// FTSFileIo.h
//	This is the file input/output functions used by the 
//	Thorlabs Fourier Transform Spectrometer DLL. 
//	
//	Copyright (c) 2019 Thorlabs Sweden, All Rights Reserved
////////////////////////////////////////////////////////////////////////

#include "FTSData.h"
#include <stdbool.h>

#ifndef FTSFILEIO_H
#define FTSFILEIO_H

#ifdef __cplusplus
extern "C" {
#endif 

    /** \addtogroup grp_file_io File IO
    *   @{
    */

    /** This will check the given file and count the number of
    *   spectra stored in it
    *   @param[in] fileName the full filename and path to the file to investigate (NULL terminated string)
    *       This may be encoded in UTF8 or ASCII.
    *   @return number of spectra in the file, -1 if an error occurs when parsing the file */
    int __cdecl   FTS_CountSpectraInFile(const char* fileName);
    int __stdcall FTS_CountSpectraInFile_std(const char* fileName);

    /** This stores the given spectrum to file in the given data format in the given location
    *   @param[in] spec The spectrum (or interferogram) to write to file
    *   @param[in] fileNameAndPath The full path and file name on the local computer 
    *       where the file should be saved. The directory must exist before calling
    *       this function. This may be encoded in UTF8 or ASCII.
    *   @param[in] fileFormat Must be one of the 'FILE_FORMAT_' flags defined in FTSData.h
    *   @return #FTS_SUCCESS on success 	
    *       NOTICE: The spectrum in memory must not change until this operation returns! */
    int __cdecl   FTS_WriteSpectrum(const spectrum_t* spec, const char* fileNameAndPath, int fileFormat);
    int __stdcall FTS_WriteSpectrum_std(const spectrum_t* spec, const char* fileNameAndPath, int fileFormat);

    /** This stores the given spectrum to file in the given data format in the given location
    *       asynchronously. This will return immediately with #FTS_SUCCESS if the write operation was started.
    *   @param[in] spec The spectrum (or interferogram) to write to file
    *   @param[in] fileNameAndPath The full path and file name on the local computer 
    *       where the file should be saved. The directory must exist before calling
    *       this function. This may be encoded in UTF8 or ASCII.
    *   @param[in] fileFormat Must be one of the 'FILE_FORMAT_' flags defined in FTSData.h
    *   @param[in] identifier An identifier for this call. The callback will be called with this
    *       identifier value as parameter 2.
    *   @param[in] fCancel If this flag is set to true during the write operation, then the writing
    *       is canceled. Otherwise will the operation continue to completion.
    *   @param[in] pCallbackFunction Pointer to a routine that will be called with notifications on
    *       the progress of the writing and possible error messages. See \ref sec_callback_asynch_fileio for reference to the 
    *       parameters passed to the callback. Must not be NULL.Must not be NULL.
    *   @return #FTS_SUCCESS on success.
    *   NOTICE: The spectrum in memory must not change until the callback is called with operation success !
    *   The callback is called with the following parameters
    *   - parameter 1: the identifier of this call. Same as the parameter 'identifier' in the call to FTS_WriteSpectrum_async().
    *   - parameter 2: status flag - set to FTS_SUCCESS if the operation is done, FTS_CANCEL if the operation
    *   has stopped due to cancellation request, FTS_ERROR_FILEERROR if the operation has finished due
    *   to a file error or FTS_PROGRESS if the writing is still in progress. 
    *   - parameter 3: the progress of the writing, in percent.
    */
    int __cdecl   FTS_WriteSpectrum_async(const spectrum_t* spec, const char* fileNameAndPath, int fileFormat, unsigned int identifier, bool *fCancel, void *data, void (__cdecl *pCallbackFunction)(unsigned int /*identifier*/, int /*status*/, double /*percent_progress*/, void * /*data*/, const char errorMessages[512]));
    int __stdcall FTS_WriteSpectrum_async_std(const spectrum_t* spec, const char* fileNameAndPath, int fileFormat, unsigned int identifier, bool *fCancel, void *data, void (__stdcall *pCallbackFunction)(unsigned int /*identifier*/, int /*status*/, double /*percent_progress*/, void * /*data*/, const char errorMessages[512]));

    /** This stores the given array of spectra to file in the given data format in the given location
    *   @param[in] spec An array with the spectra (or interferograms) to write to file
    *   @param[in] specNum The length of the array 'spec' (i.e. the number of spectra to write). This must be smaller than #MAX_SPECTRA_PER_FILE
    *   @param[in] fileNameAndPath The full path and file name on the local computer 
    *       where the file should be saved. The directory must exist before calling
    *       this function. This may be encoded in UTF8 or ASCII.
    *   @param[in] fileFormat Must be a file format that supports multiple spectra in one file, currently only #FILE_FORMAT_SPF2 or #FILE_FORMAT_TXT
    *   @return #FTS_SUCCESS on success
    */
    int __cdecl   FTS_WriteSpectra(spectrum_t** spec, unsigned int specNum, const char* fileNameAndPath, int fileFormat);
    int __stdcall FTS_WriteSpectra_std(spectrum_t** spec, unsigned int specNum, const char* fileNameAndPath, int fileFormat);

    /** This stores the given array of spectra to file in the given data format in the given location
    *       asynchronously. This will return immediately with #FTS_SUCCESS if the write operation was started.
    *   @param[in] spec An array with the spectra (or interferograms) to write to file
    *   @param[in] specNum The length of the array 'spec' (i.e. the number of spectra to write). This must be smaller than #MAX_SPECTRA_PER_FILE
    *   @param[in] fileNameAndPath The full path and file name on the local computer 
    *       where the file should be saved. The directory must exist before calling
    *       this function. This may be encoded in UTF8 or ASCII.
    *   @param[in] fileFormat Must be a file format that supports multiple spectra in one file, currently only #FILE_FORMAT_SPF2 or #FILE_FORMAT_TXT
    *   @param[in] identifier An identifier for this call. The callback will be called with this
    *       identifier value as parameter 2.
    *   @param[in] fCancel If this flag is set to true during the write operation, then the writing
    *       is canceled. Otherwise will the operation continue to completion.
    *   @param[in] pCallbackFunction Pointer to a routine that will be called with notifications on
    *       the progress of the writing and possible error messages. See \ref sec_callback_asynch_fileio for reference to the 
    *       parameters passed to the callback. Must not be NULL.Must not be NULL.
    *   @return #FTS_SUCCESS on success.
    *   NOTICE: The spectra in memory must not change until the callback is called with operation success !
    *   The callback is called with the following parameters
    *   - parameter 1: the identifier of this call. Same as the parameter 'identifier' in the call to FTS_WriteSpectrum_async().
    *   - parameter 2: status flag - set to FTS_SUCCESS if the operation is done, FTS_CANCEL if the operation
    *   has stopped due to cancellation request, FTS_ERROR_FILEERROR if the operation has finished due
    *   to a file error or FTS_PROGRESS if the writing is still in progress. 
    *   - parameter 3: the progress of the writing, in percent. */
    int __cdecl   FTS_WriteSpectra_async(spectrum_t** spec, unsigned int specNum, const char* fileNameAndPath, int fileFormat, unsigned int identifier, bool *fCancel, void *data, void (__cdecl *pCallbackFunction)(unsigned int /*identifier*/, int /*status*/, double /*percent_progress*/, void * /*data*/, const char errorMessages[512]));
    int __stdcall FTS_WriteSpectra_async_std(spectrum_t** spec, unsigned int specNum, const char* fileNameAndPath, int fileFormat, unsigned int identifier, bool *fCancel, void *data, void (__stdcall *pCallbackFunction)(unsigned int /*identifier*/, int /*status*/, double /*percent_progress*/, void * /*data*/, const char errorMessages[512]));

    /** This reads in a spectrum from the given file into the given spectrum
    *   @param[in] fileNameAndPath The full path and file name on the local computer 
    *       from which the spectrum should be read in. This may be encoded in UTF8 or ASCII.
    *   @param[out] spec Will on return be filled with the read in spectrum. If the allocated size of the spectrum
    *       is smaller than the data to read then the spectrum will be re-allocated to be able to hold the spectrum.
    *       NOTICE: It is important that the spectrum is allocated using #FTS_Trace_Allocate since the re-allocation
    *       checks the debug-info attached to the spectrum.
    *   @return #FTS_SUCCESS on success
    */
    int __cdecl   FTS_ReadSpectrum(spectrum_t* spec, const char* fileNameAndPath);
    int __stdcall FTS_ReadSpectrum_std(spectrum_t* spec, const char* fileNameAndPath);

    /** This reads in a spectrum from the given file into the given spectrum
    *   @param[out] spec Will on return be filled with the read in spectrum. If the allocated size of the spectrum
    *       is smaller than the data to read then the spectrum will be re-allocated to be able to hold the spectrum.
    *       It is therefore necessary that the spectrum is allocated either large enough to make sure the spectrum does not need to be reallocated
    *       or is one of the trace_buffer:s allocated using #FTS_Trace_Allocate. Passing in a too small buffer from an externally allocated
    *       source will result in data corruption!
    *   @param[in] fileNameAndPath The full path and file name on the local computer 
    *       from which the spectrum should be read in. This may be encoded in UTF8 or ASCII.
    *   @param[in] spectrumToRead The index (starting at zero) of the spectrum in the file to read. If the file
    *       is in any other format then #FILE_FORMAT_SPF2 then this must be zero. 
    *   NOTICE: It is important that the spectrum is allocated using #FTS_Trace_Allocate since the re-allocation
    *       checks the debug-info attached to the spectrum.
    *   @return #FTS_SUCCESS on success
    */
    int __cdecl   FTS_ReadSpectrum_indexed(spectrum_t* spec, const char* fileNameAndPath, unsigned int spectrumToRead);
    int __stdcall FTS_ReadSpectrum_indexed_std(spectrum_t* spec, const char* fileNameAndPath, unsigned int spectrumToRead);

    /** This reads in a spectrum from the given file into the given spectrum asynchronously. 
    *   This will return immediately with FTS_SUCCESS if the read operation was started. 
    *   @param[out] spec Will when the operation is finished be filled with the read in spectrum. If the allocated size of the spectrum
    *       is smaller than the data to read then the spectrum will be re-allocated to be able to hold the spectrum.
    *   NOTICE: It is important that the spectrum is allocated using #FTS_Trace_Allocate since the re-allocation
    *       checks the debug-info attached to the spectrum.
    *   @param[in] fileNameAndPath The full path and file name on the local computer 
    *       from which the spectrum should be read in. This may be encoded in UTF8 or ASCII.
    *   @param[in] identifier An identifier for this call. The callback will be called with this
    *       identifier value as parameter 2.
    *   @param[in] fCancel If this flag is set to true during the read operation, then the reading
    *       is canceled. Otherwise will the operation continue to completion.
    *   @param[in] pCallbackFunction Pointer to a routine that will be called with notifications on
    *       the progress of the reading and possible error messages. See \ref sec_callback_asynch_fileio for reference to the 
    *       parameters passed to the callback. Must not be NULL.
    *   @return #FTS_SUCCESS on success 	
    *   NOTICE: The spectrum in memory must not change until the callback is called with operation success !
    *   The callback is called with the following parameters
    *   - parameter 1: the identifier of this call. Same as the parameter 'identifier' in the call to FTS_ReadSpectrum_async.
    *   - parameter 2: status flag - set to FTS_SUCCESS if the operation is done, FTS_CANCEL if the operation
    *   has stopped due to cancellation request, FTS_ERROR_FILEERROR if the operation has finished due
    *   to a file error or FTS_PROGRESS if the reading is still in progress. 
    *   - parameter 3: the progress of the reading, in percent.
    */
    int __cdecl   FTS_ReadSpectrum_async(spectrum_t* spec, const char* fileNameAndPath, unsigned int spectrumToRead, unsigned int identifier, bool *fCancel, void *data, void (__cdecl *pCallbackFunction)(unsigned int /*identifier*/, int /*status*/, double /*percent_progress*/, void * /*data*/, const char errorMessages[512]));
    int __stdcall FTS_ReadSpectrum_async_std(spectrum_t* spec, const char* fileNameAndPath, unsigned int spectrumToRead, unsigned int identifier, bool *fCancel, void *data, void (__stdcall *pCallbackFunction)(unsigned int /*identifier*/, int /*status*/, double /*percent_progress*/, void * /*data*/, const char errorMessages[512]));

    /** This reads in an array of spectra from the given file file in the given data format in the given location
    *   @param[out] spec An array with the spectra (or interferograms) to read in
    *   @param[in] specNum The length of the array 'spec' (i.e. the number of spectra to read)
    *   @param[in] fileNameAndPath The full path and file name on the local computer 
    *       where the file can be found. This may be encoded in UTF8 or ASCII.
    *   @param[in] fileFormat Must be a file format that supports multiple spectra in one file, currently only #FILE_FORMAT_SPF2
    *   @param[out] spectraRead Will on successful return be filled with the number of spectra actually read in. Must not be null.
    *   @return #FTS_SUCCESS on success
    */
    int __cdecl   FTS_ReadSpectra(spectrum_t** spec, unsigned int specNum, const char* fileNameAndPath, int fileFormat, unsigned int* spectraRead);
    int __stdcall FTS_ReadSpectra_std(spectrum_t** spec, unsigned int specNum, const char* fileNameAndPath, int fileFormat, unsigned int* spectraRead);

    /** This is same as FTS_ReadSpectra but the read in spectra will be stored in the trace_buffer:s in FTSData.
    *   Up to MAX_TRACEBUFFER_NUM spectra can be read this way */
    int __cdecl   FTS_ReadSpectra_Traces(unsigned int specNum, const char* fileNameAndPath, int fileFormat, unsigned int* spectraRead);
    int __stdcall FTS_ReadSpectra_Traces_std(unsigned int specNum, const char* fileNameAndPath, int fileFormat, unsigned int* spectraRead);

    /** This is same as FTS_ReadSpectra_Traces but the operation will be performed asynchronously.
    *   Up to MAX_TRACEBUFFER_NUM spectra can be read this way */
    int __cdecl   FTS_ReadSpectra_Traces_async(unsigned int specNum, const char* fileNameAndPath, unsigned int* spectraRead, unsigned int identifier, bool *fCancel, void *data, void (__cdecl *pCallbackFunction)(unsigned int /*identifier*/, int /*status*/, double /*percent_progress*/, void * /*data*/, const char errorMessages[512]));
    int __stdcall FTS_ReadSpectra_Traces_async_std(unsigned int specNum, const char* fileNameAndPath, unsigned int* spectraRead, unsigned int identifier, bool *fCancel, void *data, void (__stdcall *pCallbackFunction)(unsigned int /*identifier*/, int /*status*/, double /*percent_progress*/, void * /*data*/, const char errorMessages[512]));

    ///@}

#ifdef __cplusplus
}
#endif 

#endif