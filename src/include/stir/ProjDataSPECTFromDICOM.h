/*
    Copyright (C) 2018, Institute of Nuclear Medicine, UCL
    This file is part of STIR.

    This file is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.

    This file is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    See STIR/LICENSE.txt for details
*/
/*!
  \file
  \ingroup projdata
  \brief Declaration of class stir::ProjDataSPECTFromDICOM

  \author Benjamin A. Thomas
*/

#ifndef __stir_ProjDataSPECTFromDICOM_H__
#define __stir_ProjDataSPECTFromDICOM_H__

#include "stir/ExamInfo.h"
#include "stir/ProjDataInMemory.h"
#include "stir/Scanner.h"

#include "boost/shared_array.hpp"

#include <gdcmReader.h>
#include <gdcmStringFilter.h>

#include <string>

START_NAMESPACE_STIR

class Succeeded;

enum class EnergyWindowInfo { LowerThreshold, UpperThreshold, WindowName };

bool is_spect_dicom_file(const char * dicom_filename);
Succeeded get_dicom_tag_info(const gdcm::File &file, const gdcm::Tag tag, std::string &dst);
Succeeded get_energy_window_info(const gdcm::File &file, const EnergyWindowInfo request,  std::string &dst);

std::shared_ptr<ProjDataInMemory> read_spect_dicom(const std::string& filename);

/*!
  \ingroup projdata
  \brief A class which reads/writes SPECT projection data.

*/
class ProjDataSPECTFromDICOM : public ProjDataInMemory {
 public:

//! constructor with only info, but no data
/*!
  \param proj_data_info_ptr object specifying all sizes etc.
    The ProjDataInfo object pointed to will not be modified.
  \param initialise_with_0 specifies if the data should be set to 0.
      If \c false, the data is undefined until you set it yourself.
*/
  ProjDataSPECTFromDICOM (shared_ptr<ExamInfo> const& exam_info_sptr,
                          shared_ptr<ProjDataInfo> const& proj_data_info_ptr,
                          const bool initialise_with_0 = true);

  //! constructor that copies data from another ProjData
  ProjDataSPECTFromDICOM (const ProjData& proj_data);

  //! destructor deallocates all memory the object owns
  virtual ~ProjDataSPECTFromDICOM();

};

END_NAMESPACE_STIR

#endif