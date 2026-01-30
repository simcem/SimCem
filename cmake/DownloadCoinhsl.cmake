# Script to download and extract CoinHSL archive
# This is executed by the download-extract-coinhsl custom target

set(COINHSL_ARCHIVE_URL "https://cloud.bansci.com/public.php/dav/files/iRe55SEGaZB5fZj")
set(COINHSL_ARCHIVE_FILE "${CMAKE_SOURCE_DIR}/coinhsl-archive-2024.05.15.tar")
set(COINHSL_EXTRACT_DIR "${CMAKE_SOURCE_DIR}/submodule/ThirdParty-HSL")
set(COINHSL_EXTRACT_NAME "coinhsl-archive-2024.05.15")

# Check if already extracted
if(NOT EXISTS "${COINHSL_EXTRACT_DIR}/coinhsl")
  # Download the archive
  if(NOT EXISTS "${COINHSL_ARCHIVE_FILE}")
    message(STATUS "Downloading CoinHSL archive from ${COINHSL_ARCHIVE_URL}")
    file(
      DOWNLOAD "${COINHSL_ARCHIVE_URL}"
      "${COINHSL_ARCHIVE_FILE}"
      STATUS download_status
      SHOW_PROGRESS
    )
    
    list(GET download_status 0 download_status_code)
    if(NOT download_status_code EQUAL 0)
      list(GET download_status 1 download_status_msg)
      message(FATAL_ERROR "Failed to download CoinHSL archive: ${download_status_msg}")
    endif()
  endif()

  # Extract the archive
  message(STATUS "Extracting CoinHSL archive to ${COINHSL_EXTRACT_DIR}")
  file(ARCHIVE_EXTRACT INPUT "${COINHSL_ARCHIVE_FILE}" DESTINATION "${COINHSL_EXTRACT_DIR}")
  
  # Rename extracted directory to 'coinhsl'
  if(EXISTS "${COINHSL_EXTRACT_DIR}/${COINHSL_EXTRACT_NAME}")
    message(STATUS "Renaming ${COINHSL_EXTRACT_NAME} to coinhsl")
    file(RENAME "${COINHSL_EXTRACT_DIR}/${COINHSL_EXTRACT_NAME}" "${COINHSL_EXTRACT_DIR}/coinhsl")
  elseif(NOT EXISTS "${COINHSL_EXTRACT_DIR}/coinhsl")
    message(FATAL_ERROR "Could not find extracted directory. Expected either '${COINHSL_EXTRACT_DIR}/${COINHSL_EXTRACT_NAME}' or '${COINHSL_EXTRACT_DIR}/coinhsl'")
  endif()

  message(STATUS "CoinHSL archive extracted and ready")
else()
  message(STATUS "CoinHSL already extracted at ${COINHSL_EXTRACT_DIR}/coinhsl")
endif()
