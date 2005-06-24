//
// $Id$
//
/*
Copyright (C) 2005- $Date$, Hammersmith Imanet Ltd
See STIR/LICENSE.txt for details
*/
/*!
\file
\ingroup utilities
\brief   

  \author Charalampos Tsoumpas
  \author Kris Thielemans
  
	$Date$
	$Revision$
	
	  \par Usage:
	  \code
	   correct_for_scatter [attenuation_image]
						   [no_scatter_viewgram]
						   [scatter_viewgram]
						   [scaled_scatter_filename]
						   [attenuation_threshold]
						   [global_scale_factor]
							
	  Output: Viewgram with name scaled_scatter_filename
              
	  \endcode
	  \param attenuation_threshold defaults to .05 cm^-1	  
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "stir/ProjDataInfoCylindricalNoArcCorr.h"
#include "stir/ProjDataInterfile.h"
#include "stir/utilities.h"
#include "local/stir/Scatter.h"
#include "stir/IndexRange2D.h" 
#include "stir/stream.h"
//#include "stir/Bin.h" 
#ifndef STIR_NO_NAMESPACES
using std::endl;
using std::cout;
using std::cerr;
#endif

/***********************************************************/     

int main(int argc, const char *argv[])                                  
{         
	USING_NAMESPACE_STIR
		using namespace std;
	if (argc< 5 || argc>7)
	{
	   cerr << "Usage:" << argv[0] << "\n"
			<< "\t[attenuation_image]\n"
			<< "\t[no_scatter_projdata]\n"
			<< "\t[scatter_projdata]\n" 
			<< "\t[scale_factor_per_sinogram]\n"
			<< "\t[attenuation_threshold]\n"
			<< "\t[scale_factor_per_sinogram]\n"
			<< "\tattenuation_threshold defaults to 1 cm^-1\n" 
			<< "\tusing defaults to 1 for scaling per sinogram"	;		
		return EXIT_FAILURE;            
	}      
	const float attenuation_threshold = argc>=6 ? atof(argv[5]) : 1 ;
	const int est_scale_factor_per_sino = argc>=7 ? atoi(argv[6]) : 1 ; 
	
	shared_ptr< DiscretisedDensity<3,float> >  	
		density_image_sptr= 
		DiscretisedDensity<3,float>::read_from_file(argv[1]);
	
	warning("\nWARNING: Attenuation image data are supposed to be in units cm^-1\n"
		"\tReference: water has mu .096 cm^-1\n" 
		"\tMax in attenuation image: %g\n" ,
		density_image_sptr->find_max());

	shared_ptr<ProjData> template_proj_data_sptr = ProjData::read_from_file(argv[2]);  
	const ProjDataInfo* proj_data_info_ptr =
		dynamic_cast<ProjDataInfo const *>(
		template_proj_data_sptr->get_proj_data_info_ptr());
	
	if (proj_data_info_ptr==0 || density_image_sptr==0)
		error("Check the input files\n");
	const DiscretisedDensityOnCartesianGrid<3,float>& density_image = 
		dynamic_cast<const DiscretisedDensityOnCartesianGrid<3,float>&  > 
		(*density_image_sptr.get());
	
	string scaled_scatter_filename(argv[4]);    			
	ProjDataInterfile scaled_scatter_proj_data(proj_data_info_ptr->clone(), scaled_scatter_filename);
	
	string att_proj_data_filename("att_proj_data");
	ProjDataInterfile att_proj_data(proj_data_info_ptr->clone(), att_proj_data_filename,ios::out);

	const shared_ptr<ProjData> no_scatter_proj_data_sptr = ProjData::read_from_file(argv[2]);  
	//const ProjDataInfo * projdata_info_ptr = 
    //(*no_scatter_proj_data_sptr).get_proj_data_info_ptr();  
	
	const shared_ptr<ProjData> scatter_proj_data_sptr = ProjData::read_from_file(argv[3]);   
	//const ProjDataInfo * projdata_info_ptr = 
    //(*scatter_proj_data_sptr).get_proj_data_info_ptr();

	estimate_att_viewgram(att_proj_data, density_image);
	Array<2,float> scale_factors;
	if (est_scale_factor_per_sino==1)
	scale_factors =
	scale_factors_per_sinogram(
	    no_scatter_proj_data_sptr, 
		scatter_proj_data_sptr, 
		att_proj_data,
		attenuation_threshold);

	std::cerr << scale_factors;
	scale_scatter_per_sinogram(scaled_scatter_proj_data, 
		scatter_proj_data_sptr, scale_factors) ;
//	inverse_SSRB(scaled_proj_data_3D, scaled_scatter_proj_data);

	return EXIT_SUCCESS;
}                 
