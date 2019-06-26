import os
import logging
from pprint import pprint, pformat

import numpy as np

def main():
  ''' SDE input script generator for a customized structure with fixed charge distribution within oxide
      
  '''
  
  ### General settings
  tag_special_geo = '###SPECIAL_TAG_DO_NOT_REMOVE_poly_solution_regions_geometry###'
  tag_special_doping = '###SPECIAL_TAG_DO_NOT_REMOVE_poly_solution_regions_doping###'
  
  ###################################################
  ### Edit values within this block
  
  
  
  ### Input device parameters
  channel_length    = float(1) # (um)
  channel_thickness = float(2)# (um)
  
  sd_length    = float(1) # (um)
  sd_thickness = float(1) # (um)

  poly_thickness     = float(10) # (um)
  solution_thickness = float(10) # (um)
  
  
  
  ### Input material parameters
  channel_material_name = 'Silicon' # (string)
  channel_material_doping = float(1e13) # (/cm^3) 
  channel_material_type = 'n-type' # (string: ['n-type', 'p-type'])
  
  sd_material_name = 'Silicon' # (string)
  sd_material_doping = float(1e19) # (/cm^3)
  sd_material_type = 'n-type' # (string: ['n-type', 'p-type'])
  
  poly_material_name = 'PolySilicon' # (string)
  solution_material_name = 'PolySilicon' # (string)
  
  ### Input experimental data parameters
  
  data_filepath = 'dipole_data.csv' 
  # The first two rows are headers and skipped
  # The first column is excepted to be position in nanometer (nm)
  # The second columns and onward are excepted to be fixed charge at various times. 
  # The unit for the charge are (/cm^3)

  data_solution_pos = float(0) # (um) position of the boundary of the solution
  data_interface_pos = float(100) # (um) position of the last solution data point 
  data_poly_pos = float(200) # (um) position of the boundary of the polymer

  ### Input meshing parameters
  
  # Channel, source, and drain meshing specifications
  # x-direction is along the channel; y-direction is thickness.
  # All in units of um
  csd_mesh = ''
  csd_mesh += '(sdedr:define-refinement-size "RD_channel" 10 0.1 0.001 1e-3) ; x_max y_max x_min y_min\n'
  csd_mesh += '(sdedr:define-refinement-region "RP_channel" "RD_channel" "channel_region")\n'
  csd_mesh += '(sdedr:define-refinement-size "RD_source" 0.1 0.1 0.001 1e-3) ; x_max y_max x_min y_min\n'
  csd_mesh += '(sdedr:define-refinement-region "RP_source" "RD_source" "source_region")\n'
  csd_mesh += '(sdedr:define-refinement-size "RD_drain" 0.1 0.1 0.001 1e-3) ; x_max y_max x_min y_min\n'
  csd_mesh += '(sdedr:define-refinement-region "RP_drain" "RD_drain" "drain_region")\n'
  csd_mesh += '\n'
  csd_mesh += '(sdedr:define-refinement-size "RD_solution" 0.1 10 0.001 1e-9) ; x_max y_max x_min y_min\n'
  csd_mesh += '(sdedr:define-refinement-material "RP_solution" "RD_solution" "'+solution_material_name+'")\n'
  csd_mesh += '(sdedr:define-refinement-size "RD_poly" 0.1 10 0.001 1e-9) ; x_max y_max x_min y_min\n'
  csd_mesh += '(sdedr:define-refinement-material "RP_poly" "RD_poly" "'+poly_material_name+'")\n'
  csd_mesh += '\n'
  csd_mesh += '(sdedr:define-refinement-function "RD_channel" "MaxLenInt" "channel_region" "source_region" 1e-3 1.1 "DoubleSide" "UseRegionNames")\n'
  csd_mesh += '(sdedr:define-refinement-function "RD_channel" "MaxLenInt" "channel_region" "drain_region" 1e-3 1.1 "DoubleSide" "UseRegionNames")\n'
  csd_mesh += '(sdedr:define-refinement-function "RD_source" "MaxLenInt" "channel_region" "source_region" 1e-3 1.5 "DoubleSide" "UseRegionNames")\n'
  csd_mesh += '(sdedr:define-refinement-function "RD_drain" "MaxLenInt" "channel_region" "drain_region" 1e-3 1.5 "DoubleSide" "UseRegionNames")\n'
  csd_mesh += '\n'
  csd_mesh += '(sdedr:define-refinement-function "RD_channel" "MaxLenInt" "'+channel_material_name+'" "'+poly_material_name+'" 0.5e-3 1.05 "UseMaterialNames")\n'
  
  ### Input numerical parameters
  
  solution_slice_num = 100
  # (um) Number of sub-divided slices within solution
  # If specified with a number, the total charge will be divided equally into these slices
  # If None is specified, the sub-dividied slices will use positions given by input CSV file
  
  poly_slice_num = 100
  # (um) same as solution_slice_num but for polymer layer
  
  output_SDE_filename = 'dipole_sensor_sde.cmd'

  ### End of block for editing
  ###################################################
  
  

  ### PREPROCESS: Check inputs for some consistency issues and prepare inputfile
  if sd_thickness > channel_thickness:
    logging.error('S/D thickness is larger than channel thickness.')
    
  # load experimental datafile
  exp_data = np.genfromtxt(data_filepath, dtype=float, delimiter=',', skip_header=2) 

  total_time_instances = exp_data.shape[1]-1
  
  if total_time_instances < 1:
    logging.error('Unable to proceed. Experiment inputfile has fewer than 2 columns?')
    
  charge_y_pos = exp_data[:,0]*1e-3 # (nm -> um)
  
  if data_solution_pos < data_poly_pos:
    solution_index = np.less_equal(charge_y_pos, data_interface_pos)
    poly_index = np.greater(charge_y_pos, data_interface_pos)
  else:
    poly_index = np.less(charge_y_pos, data_interface_pos)
    solution_index = np.greater_equal(charge_y_pos, data_interface_pos)

  solution_y_pos = charge_y_pos[solution_index] # position of the charge inside solution
  poly_y_pos = charge_y_pos[poly_index] # position of the charge inside poly
    

  if (np.absolute(data_solution_pos - data_interface_pos)*1e-3) < solution_thickness:
    logging.info('The defined thickness of solution layer is LARGER than experimental data provided.')
  else:
    logging.info('The defined thickness of solution layer is SMALLER than experimental data provided. Fixed charge outside of the layer is omitted.')
  if (np.absolute(data_poly_pos - data_interface_pos)*1e-3) < poly_thickness:
    logging.info('The defined thickness of poly layer is LARGER than experimental data provided.')
  else:
    logging.info('The defined thickness of poly layer is SMALLER than experimental data provided. Fixed charge outside of the layer is omitted.')



    
  ### GENERATE: Source + Drain + Channel regions
  # Source is located at left side. x=0 at the leftmost side of the Source
  # Channel/Polymer interface is y=0
  
  sde_str  = ''
  sde_str += ';--- Construct geometry ------------------------------------------\n'
  sde_str += '(sdegeo:set-default-boolean "ABA")\n'
    
  # add channel
  sde_str += cr_rec(0,(sd_length+sd_length+channel_length), \
                    0,channel_thickness, \
                    channel_material_name, 'channel_region') + '\n'
  
  # add source, overwriting the channel region
  sde_str += cr_rec(0,sd_length, \
                    0,sd_thickness, \
                    sd_material_name, 'source_region') + '\n'
  
  # add drain, overwriting the channel region
  sde_str += cr_rec((sd_length+channel_length),(sd_length+sd_length+channel_length), \
                    0,sd_thickness, \
                    sd_material_name, 'drain_region') + '\n'
  sde_str += '\n\n\n'
  
  ### GENERATE: Polymer and Solution regions
  # Create regions only. The fixed charge are added as dopings
  sde_str += tag_special_geo
  
  sde_str_solution_geo = list()
  sde_str_solution_doping = list()
  sde_str_poly_geo = list()
  sde_str_poly_doping = list()
  
  for this_index in range(0, total_time_instances):
    # for each of the time instance from experimental file
    this_solution_conc = exp_data[solution_index,this_index+1]
    this_poly_conc = exp_data[poly_index,this_index+1]

    # arrange data so that the pos = 0 point starts at interface for both solution and poly data
    if data_solution_pos < data_interface_pos:
      solution_y_pos = solution_y_pos - np.max(solution_y_pos)
      poly_y_pos = poly_y_pos - np.min(poly_y_pos)
      
    else:
      solution_y_pos = solution_y_pos - np.min(solution_y_pos)
      poly_y_pos = poly_y_pos - np.max(poly_y_pos)


    
    geo_str, doping_str = create_charged_region(solution_y_pos, this_solution_conc, \
                                                0, (sd_length*2.0+channel_length), \
                                                solution_thickness, (solution_thickness+poly_thickness), \
                                                solution_slice_num, \
                                                solution_material_name, 'sol')

    sde_str_solution_geo.append(geo_str)
    sde_str_solution_doping.append(doping_str)
    
    geo_str, doping_str = create_charged_region(poly_y_pos, this_poly_conc, \
                                                0, (sd_length*2.0+channel_length), \
                                                poly_thickness, 0, \
                                                poly_slice_num, \
                                                poly_material_name, 'poly') 
                                          
    sde_str_poly_geo.append(geo_str)
    sde_str_poly_doping.append(doping_str) 
       
  sde_str += '\n\n\n'
    


  ### DOPING:
  sde_str += ';--- Doping profile --------------------------------------------\n'
  sde_str += '(sdedr:define-constant-profile "channel_doping" "'+ \
             ('BoronActiveConcentration' if channel_material_type == 'p-type' else 'PhosphorusActiveConcentration')+ \
             '" '+np.format_float_scientific(channel_material_doping, precision=5, trim='-')+')\n'
  sde_str += '(sdedr:define-constant-profile-region "Place_channel_doping" "channel_doping" "channel_region")\n'
  
  sde_str += '(sdedr:define-constant-profile "sd_doping" "'+ \
             ('BoronActiveConcentration' if sd_material_type == 'p-type' else 'PhosphorusActiveConcentration')+ \
             '" '+np.format_float_scientific(sd_material_doping, precision=5, trim='-')+')\n'
  sde_str += '(sdedr:define-constant-profile-region "Place_s_doping" "sd_doping" "source_region")\n'
  sde_str += '(sdedr:define-constant-profile-region "Place_d_doping" "sd_doping" "drain_region")\n'
  sde_str += '\n'
  sde_str += tag_special_doping
  sde_str += '\n\n\n'
  
  
  
  ### CONTACT
  sde_str += ';--- Contacts ------------------------------------------\n'
  sde_str += '(sdegeo:define-contact-set "s_contact" 4 (color:rgb 1 1 1) "##")\n'
  sde_str += '(sdegeo:define-contact-set "d_contact" 4 (color:rgb 0 0 0) "##")\n'

  sde_str += '(sdegeo:define-2d-contact (find-edge-id (position '+ \
                  np.format_float_scientific(0, precision=5, trim='-')+' '+ \
                  np.format_float_scientific(sd_thickness/2.0, precision=5, trim='-')+' 0)) "s_contact")\n'
                  
  sde_str += '(sdegeo:define-2d-contact (find-edge-id (position '+ \
                  np.format_float_scientific(sd_length*2.0+channel_length, precision=5, trim='-')+' '+ \
                  np.format_float_scientific(sd_thickness/2.0, precision=5, trim='-')+' 0)) "d_contact")\n'
  
  
  
  ### MESH:
  sde_str += ';--- Mesh refinements -----------------------------------\n'
  sde_str += csd_mesh
  sde_str += '\n\n\n'



  ### FINALIZE: Save files and print some statistics
  
  # Mesh
  sde_str += ';--- Build mesh ------------------------------------------------\n'
  sde_str += '(sde:setrefprops 90 90 0 0)\n'
  sde_str += '(sde:build-mesh "snmesh" "-a -c boxmethod" "n@node@_msh")\n'

  #TODO: write to file
  for this_index in range(0, total_time_instances):

    this_output_filename = 's'+str(this_index)+'_'+output_SDE_filename
    
    this_sde_str = sde_str
    this_sde_str = this_sde_str.replace(tag_special_geo, sde_str_solution_geo[this_index]+sde_str_poly_geo[this_index])
    this_sde_str = this_sde_str.replace(tag_special_doping, sde_str_solution_doping[this_index]+sde_str_poly_doping[this_index])
    
    with open(this_output_filename, 'w+') as f:
      f.write(this_sde_str)
    






def cr_rec(x1,x2,y1,y2,material_name,region_name):
  ''' Create a rectangular region
  '''

  output_str = '(sdegeo:create-rectangle ' + \
               '(position '+' '.join([np.format_float_scientific(x, precision=15, trim='-') for x in [x1,y1,0]])+') ' + \
               '(position '+' '.join([np.format_float_scientific(x, precision=15, trim='-') for x in [x2,y2,0]])+') ' + \
               '"'+material_name +'"' \
               ' "'+region_name+'")'
  
  return output_str
  
  



def create_charged_region(y_pos, this_conc, x1, x2, y1, y2, slice_num, material_name, id_tag):
  ''' Given the concentration profile, thickness of the layer, and material name,
      return two strings containing stacks of slices and doping concentration of these slices.
      y_pos is y1->y2. All in units of um. y1 should always be the interface
  '''
  
  sde_str_geo = ''
  sde_str_doping = ''
  

  # Cut-off the charge profile based on layer thickness (y1-y2), with y_pos going y1->y2
  keep_index = np.less(np.absolute(y_pos), np.absolute(y2-y1))
  
  y_pos = y_pos[keep_index]
  this_conc = this_conc[keep_index]
  
  # layout background material to make sure entire area is covered
  sde_str_geo += cr_rec(x1,x2,-1.0*y1,-1.0*y2,material_name,'reg_'+material_name+'_bg')+'\n'

  # Find total fixed charge concentration
  total_conc = np.absolute(np.trapz(y_pos, this_conc))
  this_doping_type = ('p-type' if np.trapz(y_pos, this_conc) > 0 else 'n-type')
  
  accum_conc = np.empty(y_pos.shape)
  for index, _ in enumerate(y_pos):
    accum_conc[index] = np.absolute(np.trapz(y_pos[0:index], this_conc[0:index]))
  
  y_interp_pos = np.interp(np.linspace(0,total_conc,slice_num+1), accum_conc, y_pos) # obtain interpolated positions
  
  conc_interp = np.interp(y_interp_pos, y_pos, this_conc) # obtain concentration at interpolated positions
  
  conc_interp = conc_interp*1.0e14
  
  for this_index, _ in enumerate(y_interp_pos[1:]): 
    this_region_name = 'reg_'+material_name+'_'+str(this_index)+'_'+id_tag
    # Create geometry
    # Charged regions occupy negative Y axis
    sde_str_geo += cr_rec(x1,x2,-1.0*y1+y_interp_pos[this_index],-1.0*y1+y_interp_pos[this_index+1],material_name,this_region_name)+'\n'

    # Create doping
    sde_str_doping += '(sdedr:define-constant-profile "'+this_region_name+'_doping'+'" "'+ \
                      ('BoronActiveConcentration' if this_doping_type == 'p-type' else 'PhosphorusActiveConcentration')+ \
                      '" '+np.format_float_scientific(np.absolute(conc_interp[this_index]), precision=15, trim='-')+')\n'
                      
    sde_str_doping += '(sdedr:define-constant-profile-region "Place_'+this_region_name+'_doping'+ \
                      '" "'+this_region_name+'_doping'+'" "'+this_region_name+'")\n'
    
  return sde_str_geo, sde_str_doping











if __name__ == '__main__':
  main()
