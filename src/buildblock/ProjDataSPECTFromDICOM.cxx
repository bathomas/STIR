/*!

  \file
  \ingroup projdata
  \brief Implementations for non-inline functions of class stir::ProjDataSPECTFromDICOM

  \author Benjamin A. Thomas
*/
/*
    Copyright (C) 2018, Institute of Nuclear Medicine, UCL
    Copyright (C) 2018, Benjamin A. Thomas
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


#include "stir/ProjDataSPECTFromDICOM.h"
#include "stir/shared_ptr.h"
#include "stir/Succeeded.h"
#include "stir/SegmentByView.h"
#include "stir/Array.h"
#include "stir/Bin.h"
#include "stir/info.h"
#include "stir/ProjDataInfoCylindricalArcCorr.h"

#include <fstream>
#include <sstream>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#ifndef STIR_NO_NAMESPACES
using std::fstream;
using std::iostream;
using std::ios;
using std::string;
#endif

START_NAMESPACE_STIR

bool is_spect_dicom_file(const char * dicom_filename){

  std::unique_ptr<gdcm::Reader> DICOM_reader(new gdcm::Reader);
  DICOM_reader->SetFileName(dicom_filename);

  try {
    if (!DICOM_reader->Read()) {
      return false;
    }
  } catch (const std::string &e){
    error(boost::format("Error reading file %1% as DICOM file") % dicom_filename);
    error(e);
    return false;
  }

  const gdcm::File &file = DICOM_reader->GetFile();

  std::string image_type;
  if (get_dicom_tag_info(file, gdcm::Tag(0x008,0x008), image_type) == Succeeded::no){
    error(boost::format("Error reading image type (0x008,0x008) from %1%") % dicom_filename);
    return false;
  }

  info(boost::format("Image type (0x008,0x008): %1%") % image_type);

  const std::string spect_emission_image_type = "ORIGINAL\\PRIMARY\\TOMO\\EMISSION";

  if (image_type != spect_emission_image_type){
    error(boost::format("Unknown DICOM image type: %1%") % image_type);
    return false;
  }


  return true;
}

Succeeded get_dicom_tag_info(const gdcm::File &file, const gdcm::Tag tag, std::string &dst){

  //Extracts information for a given DICOM tag from a gdcm dataset.
  //Tag contents are returned as a string in dst variable.

  //Tries to read the element associated with the tag. If the read fails, the
  //DataElement should have a ByteValue of NULL.

  try {
    const gdcm::DataSet &ds = file.GetDataSet();

    gdcm::StringFilter sf;
    sf.SetFile(file);

    std::stringstream inStream;
    inStream.exceptions(std::ios::badbit);

    // First try and see if this is a standard tag.
    gdcm::DataElement element = ds.GetDataElement(tag);

    if (element.GetByteValue() != NULL) {
      dst = sf.ToString(tag);
      return Succeeded::yes;
    }

    // Try: RotationInformationSequence     (0054,0052)
    //      DetectorInformationSequence     (0054,0022)
    //      EnergyWindowInformationSequence (0054,0012)
    std::vector<gdcm::Tag> seqs = { gdcm::Tag(0x0054,0x0052), gdcm::Tag(0x0054,0x0022), gdcm::Tag(0x0054,0x0012)};

    for (const auto& t : seqs) {
      const gdcm::DataElement &de = file.GetDataSet().GetDataElement(t);
      const gdcm::SequenceOfItems *sqi = de.GetValueAsSQ();
      const gdcm::Item &item = sqi->GetItem(1);

      element = item.GetDataElement(tag);

      if (element.GetByteValue() != NULL) {
        dst = sf.ToString(element);
        return Succeeded::yes;
      }
    }

  } catch (std::bad_alloc){
    error(boost::format("get_dicom_tag_info: cannot read tag %1%") % tag);
    return Succeeded::no;
  }

  return Succeeded::no;
}

Succeeded get_energy_window_info(const gdcm::File &file, const EnergyWindowInfo request,  std::string &dst){

  if (request == EnergyWindowInfo::WindowName){
    return get_dicom_tag_info(file, gdcm::Tag(0x0054,0x0018), dst);
  }

  try {
    const gdcm::Tag energy_window_info_seq = gdcm::Tag(0x0054,0x0012);
    const gdcm::Tag energy_window_range_seq = gdcm::Tag(0x0054,0x0013);

    const gdcm::Tag lower_energy_window_tag = gdcm::Tag(0x0054,0x0014);
    const gdcm::Tag upper_energy_window_tag = gdcm::Tag(0x0054,0x0015);

    //Get Energy Window Info Sequence
    const gdcm::DataElement &de = file.GetDataSet().GetDataElement(energy_window_info_seq);
    const gdcm::SequenceOfItems *sqi = de.GetValueAsSQ();
    const gdcm::Item &item = sqi->GetItem(1);

    //Get Energy Window Range Sequence
    const gdcm::DataElement &element = item.GetDataElement(energy_window_range_seq);
    const gdcm::SequenceOfItems *sqi2 = element.GetValueAsSQ();
    const gdcm::Item &item2 = sqi2->GetItem(1);

    gdcm::DataElement window_element;

    if (request == EnergyWindowInfo::LowerThreshold)
      window_element = item2.GetDataElement(lower_energy_window_tag);
    else
      window_element = item2.GetDataElement(upper_energy_window_tag);

    if (window_element.GetByteValue() != NULL) {
      gdcm::StringFilter sf;
      sf.SetFile(file);
      dst = sf.ToString(window_element);
      return Succeeded::yes;
    }

  } catch (std::bad_alloc){
    error(boost::format("get_energy_window_info: cannot read energy info"));
    return Succeeded::no;
  }

  return Succeeded::no;
}

std::shared_ptr<ProjDataInMemory> read_spect_dicom(const std::string& filename){

  if (!is_spect_dicom_file(filename.c_str())){
    return nullptr;
  }

  std::unique_ptr<gdcm::Reader> DICOM_reader(new gdcm::Reader);
  DICOM_reader->SetFileName(filename.c_str());

  try {
    if (!DICOM_reader->Read()) {
      return nullptr;
    }
  } catch (const std::string &e){
    error(boost::format("Error reading file %1% as DICOM file") % filename);
    error(e);
    return nullptr;
  }

  const gdcm::File &file = DICOM_reader->GetFile();

  float start_angle = 0.0f;
  int num_of_projections = 0;
  float angular_step = 0.0f;
  int actual_frame_duration = 0; //frame duration in msec
  //int num_of_rotations = 0;
  std::string direction_of_rotation;
  int extent_of_rotation;

  float rotation_radius = 0.0f;

  float lower_en_window_thres = 0.0f;
  float upper_en_window_thres = 0.0f;
  std::string energy_window_name;

  int num_dimensions;
  std::vector<std::string> matrix_labels;
  std::vector<int> matrix_size;
  std::vector<double> pixel_sizes;

  std::string no_of_proj_as_str;
  if (get_dicom_tag_info(file, gdcm::Tag(0x0054,0x0053), no_of_proj_as_str) == stir::Succeeded::yes){
    num_of_projections = std::stoi(no_of_proj_as_str);
    info(boost::format("Number of projections: %1%") % num_of_projections, 2);
  }

  if (get_dicom_tag_info(file, gdcm::Tag(0x0018,0x1140), direction_of_rotation) == stir::Succeeded::yes){
    info(boost::format("Direction of rotation: %1%") % direction_of_rotation, 2);
  }

  std::string start_angle_as_string;
  if (get_dicom_tag_info(file, gdcm::Tag(0x0054,0x0200), start_angle_as_string) == stir::Succeeded::yes){
    start_angle = std::stof(start_angle_as_string);
    info(boost::format("Starting angle: %1%") % start_angle, 2);
  }

  std::string angular_step_as_string;
  if (get_dicom_tag_info(file, gdcm::Tag(0x0018,0x1144), angular_step_as_string) == stir::Succeeded::yes){
    angular_step = std::stof(angular_step_as_string);
    info(boost::format("Angular step: %1%") % angular_step, 2);
  }

  std::string extent_of_rotation_as_string;
  if (get_dicom_tag_info(file, gdcm::Tag(0x0018,0x1143), extent_of_rotation_as_string) == stir::Succeeded::yes){
    extent_of_rotation = std::stoi(extent_of_rotation_as_string);
    info(boost::format("Rotation extent: %1%") % extent_of_rotation, 2);
  }

  std::string radius_as_string;
  if (get_dicom_tag_info(file, gdcm::Tag(0x0018,0x1142), radius_as_string) == stir::Succeeded::yes){
    rotation_radius = std::stof(radius_as_string);
    info(boost::format("Radius: %1%") % radius_as_string, 2);
  }

  std::string actual_frame_duration_as_string;
  if (get_dicom_tag_info(file, gdcm::Tag(0x0018,0x1242), actual_frame_duration_as_string) == stir::Succeeded::yes){
    actual_frame_duration = std::stoi(actual_frame_duration_as_string);
    info(boost::format("Projection frame duration (msec): %1%") % actual_frame_duration, 2);
  }

  num_dimensions = 2;

  matrix_labels.push_back("axial coordinate");
  matrix_labels.push_back("bin coordinate");

  std::string matrix_size_as_string;
  if (get_dicom_tag_info(file, gdcm::Tag(0x0028,0x0010), matrix_size_as_string) == stir::Succeeded::yes){
    matrix_size.push_back(std::stoi(matrix_size_as_string));
    info(boost::format("Matrix size [1]: %1%") % matrix_size_as_string, 2);
  }

  if (get_dicom_tag_info(file, gdcm::Tag(0x0028,0x0011), matrix_size_as_string) == stir::Succeeded::yes){
    matrix_size.push_back(std::stoi(matrix_size_as_string));
    info(boost::format("Matrix size [2]: %1%") % matrix_size_as_string, 2);
  }

  std::string pixel_size_as_string;
  if (get_dicom_tag_info(file, gdcm::Tag(0x0028,0x0030), pixel_size_as_string) == stir::Succeeded::yes){
    info(boost::format("Pixel size: %1%") % pixel_size_as_string, 2);
    size_t  curr, prev = 0;

    while (curr != std::string::npos){
      curr = pixel_size_as_string.find('\\',prev);
      std::string found = pixel_size_as_string.substr(prev, curr-prev);
      pixel_sizes.push_back(std::stof(found));
      prev = curr+1;
    }
  }

  if (get_energy_window_info(file, EnergyWindowInfo::WindowName , energy_window_name) == stir::Succeeded::yes){
    info(boost::format("Energy window: %1%") % energy_window_name, 2);
  }

  std::string lower_window_as_string;
  if (get_energy_window_info(file, EnergyWindowInfo::LowerThreshold , lower_window_as_string) == stir::Succeeded::yes){
    lower_en_window_thres = std::stof(lower_window_as_string);
    info(boost::format("Lower energy window limit: %1%") % lower_en_window_thres, 2);
  }

  std::string upper_window_as_string;
  if (get_energy_window_info(file, EnergyWindowInfo::UpperThreshold , upper_window_as_string) == stir::Succeeded::yes){
    upper_en_window_thres = std::stof(upper_window_as_string);
    info(boost::format("Upper energy window limit: %1%") % upper_en_window_thres, 2);
  }

  const int num_axial_poss = matrix_size.at(1);
  const int num_bins = matrix_size.at(0);
  const int num_views = num_of_projections;

  const double z_spacing_in_cm = pixel_sizes.at(1)/10.;
  const double bin_size_in_cm = pixel_sizes.at(0)/10.;

  VectorWithOffset<float> radii(0, num_views-1);
  //TODO: Ignoring non-circular orbits for now.
  for ( int i = 0 ; i < num_views ; i++ )
    radii[ i ] = static_cast<float>(rotation_radius);

  // somewhat strange values to be compatible with PET
  VectorWithOffset<int> sorted_min_ring_diff(0,0);
  VectorWithOffset<int> sorted_max_ring_diff(0,0);
  VectorWithOffset<int> sorted_num_rings_per_segment(0,0);
  sorted_min_ring_diff[0]=0;
  sorted_max_ring_diff[0]=0;

  sorted_num_rings_per_segment[0]=num_axial_poss;

  // we construct a new scanner object with
  // data from the Interfile header (or the guessed scanner).
  // Initialize the scanner values (most are not used in SPECT reconstruction)

  const int num_rings = sorted_num_rings_per_segment[0];
  const int num_detectors_per_ring = -1;//num_views*2;
  const double average_depth_of_interaction_in_cm = 0;
  const double distance_between_rings_in_cm = z_spacing_in_cm*2; // need to do times 2  such that default z-spacing of reconstruction is z_spacing
  double default_bin_size_in_cm = bin_size_in_cm ;
  const double view_offset_in_degrees = start_angle;
  const int max_num_non_arccorrected_bins = num_bins;
  const int default_num_arccorrected_bins = num_bins;
  const int num_axial_blocks_per_bucket = -1;
  const int num_transaxial_blocks_per_bucket = -1;
  const int num_axial_crystals_per_block = -1;
  const int num_transaxial_crystals_per_block = -1;
  const int num_axial_crystals_per_singles_unit = -1;
  const int num_transaxial_crystals_per_singles_unit = -1;
  const int num_detector_layers = 1;

  shared_ptr<Scanner> scanner_sptr(
      new Scanner(Scanner::User_defined_scanner,
                  "nucmed",
                  num_detectors_per_ring,
                  num_rings,
                  max_num_non_arccorrected_bins,
                  default_num_arccorrected_bins,
                  static_cast<float>(radii[0]),
                  static_cast<float>(average_depth_of_interaction_in_cm*10),
                  static_cast<float>(distance_between_rings_in_cm*10.),
                  static_cast<float>(default_bin_size_in_cm*10),
                  static_cast<float>(view_offset_in_degrees*_PI/180),
                  num_axial_blocks_per_bucket,
                  num_transaxial_blocks_per_bucket,
                  num_axial_crystals_per_block,
                  num_transaxial_crystals_per_block,
                  num_axial_crystals_per_singles_unit,
                  num_transaxial_crystals_per_singles_unit,
                  num_detector_layers)
  );

  shared_ptr<ProjDataInfoCylindricalArcCorr> proj_data_info_sptr(
      new ProjDataInfoCylindricalArcCorr ( scanner_sptr,
                                           float(bin_size_in_cm*10.),
                                           sorted_num_rings_per_segment,
                                           sorted_min_ring_diff,
                                           sorted_max_ring_diff,
                                           num_views,num_bins)
  );

  proj_data_info_sptr->set_ring_radii_for_all_views( radii );

  const float angle_sampling = float (extent_of_rotation)/num_views * float(_PI/180);

  if(direction_of_rotation=="CW")
  {
    proj_data_info_sptr->set_azimuthal_angle_sampling(-angle_sampling);
  }
  else if(direction_of_rotation=="CCW")
  {
    proj_data_info_sptr->set_azimuthal_angle_sampling(angle_sampling);
  }

  shared_ptr<ExamInfo> exam_info_sptr(new ExamInfo());
  exam_info_sptr->imaging_modality = ImagingModality::NM;

  const gdcm::DataElement &de = file.GetDataSet().GetDataElement(gdcm::Tag(0x7fe0,0x0010));
  const gdcm::ByteValue *bv = de.GetByteValue();

  uint64_t len0 = (uint64_t)bv->GetLength()/2;

  Array<1,float> pixel_data_as_float(0, len0);

  uint16_t *ptr = (uint16_t*)bv->GetPointer();

  uint64_t ct = 0;
  while (ct < len0){
    uint16_t val = *ptr;
    pixel_data_as_float[ct] = (float)val;
    ptr++;
    ct++;
  }

  shared_ptr<ProjDataInMemory> proj_data_sptr(new ProjDataInMemory(exam_info_sptr, proj_data_info_sptr));
  proj_data_sptr->fill_from(pixel_data_as_float.begin_all_const());

  return proj_data_sptr;

}

ProjDataSPECTFromDICOM::
~ProjDataSPECTFromDICOM()
{}


END_NAMESPACE_STIR