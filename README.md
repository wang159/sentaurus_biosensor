# Biosensor Structure Generator for Sentaurus Device

This simple script is a SDE inputfile generator for a 2D Channel-Polymer-Solution structure (see figure). The SDE generator takes into account the steep change of charge distribution at the polymer/solution interface, where dipole layers forms. Non-uniform, refined regions and meshes are generated based on the charge distribution to ensure the charge profile is resolved.

![Figure](https://i.imgur.com/8CAAK8P.png)

## Input parameters
### Input device parameters
```
channel_length    = float(1) # (um) Channel length
channel_thickness = float(2) # (um) Channel thickness

sd_length    = float(1) # (um) Source/Drain length
sd_thickness = float(1) # (um) Source/Drain thickness. S/D always starts at Poly-Semiconductor interface

poly_thickness     = float(10) # (um) Polymer thickness
solution_thickness = float(10) # (um) Solution thickness. Ideally, the solution "layer" 
                               #      is semi-infinitly large. In simulation, however, specify a large thickness
                               #      to ensure boundary condition is correct 
```


### Input material parameters
```
channel_material_name = 'Silicon' # (string) Channel material. Make sure this is a material Sentaurus recognizes. See "adding a new material" in Sentaurus manual
channel_material_doping = float(1e13) # (/cm^3) 
channel_material_type = 'n-type' # (string: ['n-type', 'p-type'])

sd_material_name = 'Silicon' # (string)
sd_material_doping = float(1e19) # (/cm^3)
sd_material_type = 'n-type' # (string: ['n-type', 'p-type'])

poly_material_name = 'PolySilicon' # (string)
solution_material_name = 'PolySilicon' # (string)
```
### Input experimental data parameters
```
data_filepath = 'dipole_data.csv' 
# The first two rows are headers and skipped
# The first column is excepted to be position in nanometer (nm)
# The second columns and onward are excepted to be fixed charge at various times. 
# The unit for the charge are (/cm^3)

data_solution_pos = float(0) # (um) position of the boundary of the solution
data_interface_pos = float(100) # (um) position of the last solution data point 
data_poly_pos = float(200) # (um) position of the boundary of the polymer
```
### Input meshing parameters
```
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
```
### Input numerical parameters
```
solution_slice_num = 100
# (um) Number of sub-divided slices within solution
# If specified with a number, the total charge will be divided equally into these slices
# If None is specified, the sub-dividied slices will use positions given by input CSV file

poly_slice_num = 100
# (um) same as solution_slice_num but for polymer layer
  
output_SDE_filename = 'dipole_sensor_sde.cmd'
```
