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
class ProjDataSPECTFromDICOM;

enum class EnergyWindowInfo { LowerThreshold, UpperThreshold, WindowName };

bool is_spect_dicom_file(const char * dicom_filename);
Succeeded get_dicom_tag_info(const gdcm::File &file, const gdcm::Tag tag, std::string &dst);
Succeeded get_energy_window_info(const gdcm::File &file, const EnergyWindowInfo request,  std::string &dst);

ProjDataSPECTFromDICOM* read_spect_dicom(const std::string& filename);

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
*/
  ProjDataSPECTFromDICOM (shared_ptr<ExamInfo> const& exam_info_sptr,
                          shared_ptr<ProjDataInfo> const& proj_data_info_ptr,
                          const std::string &dicom_filename,
                          StorageOrder o = Segment_View_AxialPos_TangPos,
                          NumericType data_type = NumericType::FLOAT,
                          ByteOrder byte_order = ByteOrder::native,
                          float scale_factor = 1);

  //! constructor that copies data from another ProjData
  ProjDataSPECTFromDICOM (const ProjData& proj_data);

  //! destructor deallocates all memory the object owns
  virtual ~ProjDataSPECTFromDICOM();


  //! Read data into viewgram
  Succeeded read_data(
      const shared_ptr<std::vector<float>> &sino_data,
      uint64_t pos,
      Viewgram<float> &viewgram,
      float scale);

  //! Get & set viewgram
  Viewgram<float> get_viewgram(const int view_num, const int segment_num,const bool make_num_tangential_poss_odd=false) const;
  Succeeded set_viewgram(const Viewgram<float>& v);

  //! Get & set sinogram
  Sinogram<float> get_sinogram(const int ax_pos_num, const int segment_num,const bool make_num_tangential_poss_odd=false) const;
  Succeeded set_sinogram(const Sinogram<float>& s);

  //! Get all sinograms for the given segment
  SegmentBySinogram<float> get_segment_by_sinogram(const int segment_num) const;
  //! Get all viewgrams for the given segment
  SegmentByView<float> get_segment_by_view(const int segment_num) const;

  //! Set all sinograms for the given segment
  Succeeded set_segment(const SegmentBySinogram<float>&);
  //! Set all viewgrams for the given segment
  Succeeded set_segment(const SegmentByView<float>&);

  //! Get the value of bin.
  float get_bin_value(const Bin& this_bin) const;

protected:
  //! the stream with the data
  shared_ptr<std::vector<float>> sino_data = nullptr;
  std::string dicom_filename;

private:
  //! offset of the whole 3d sinogram in the stream
  std::streamoff  offset;


  //!the order in which the segments occur in the stream
  std::vector<int> segment_sequence;

  StorageOrder storage_order;

  NumericType on_disk_data_type;

  ByteOrder on_disk_byte_order;

  // scale_factor is only used when reading data from file. Data are stored in
  // memory as float, with the scale factor multiplied out
  float scale_factor;

  //! Calculate the offset for the given segment
  std::streamoff get_offset_segment(const int segment_num) const;

  //! Calculate offsets for viewgram data
  std::vector<uint64_t> get_offsets(const int view_num, const int segment_num) const;
  //! Calculate offsets for sinogram data
  std::vector<std::streamoff> get_offsets_sino(const int ax_pos_num, const int segment_num) const;

  //! Calculate the offsets for specific bins.
  std::vector<std::streamoff> get_offsets_bin(const Bin) const;

};

END_NAMESPACE_STIR

#include "stir/ProjDataSPECTFromDICOM.inl"

#endif