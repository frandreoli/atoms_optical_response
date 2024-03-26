#############################################################################################################
################## PHYSICAL GLOSSARY ########################################################################
#############################################################################################################
#=
- n_bulk   = refractive index of the embedding (dielectric, lossless) material (for air and vacuum n_bulk = 1)
- Gamma0  = spontaneous (elastic) decay rate of the single atoms in the material. This is set to Gamma0=1, meaning that all the frequencies are expressed in units of Gamma0
- omega0  = atomic resonance frequency between the groun |g> and the excited state |e>
- lambda0 = (2*pi*c)/(omega0*n_bulk) resonant wavelength inside the bulk dielectric material. 
All lengths are calculated in units of lambda0.
=#
#
#
#
#
#
#
#############################################################################################################
################## OVERALL OPTIONS ##########################################################################
#############################################################################################################
#
#
#CODE SETTINGS:
#GEOMETRICAL SETTINGS:
#It defines the geometry of the atomic positions                
const geometry_settings         =      
[
#
"DISORDERED_SPHERE" ;           #1
"DISORDERED_CYLINDER" ;         #2
"DISORDERED_CUBOID" ;           #3
#
"ARRAYS" ;                      #4
#
"METALENS" ;                    #5
#
"CUSTOM_POSITIONS"              #6
#
#Choose here below the number wanted
][5]
#
#
#If YES then the code saves the atomic positions in the file: "atomic_positions.h5"
const pos_save_option           =      ["YES" ; "NO"][1]
#
#
#If YES then the code saves the steady-state, atomic coefficients in the file: "atomic_coeff.h5"
const coeff_save_option         =      ["YES" ; "NO"][1]
#
#
#Possibility of punching defects in the arrays/metalens 
#i.e. removing a fraction defects_fraction of randomly chosen atoms
defects_fraction = 0.0
#
#
#Possibility of randomly shifting the atomic positions from their ordered geometry (array/metalens)
#by a randomly chosen factor (in the 3 dimensions) sampled from a Gaussian with 
#standard deviation small_disorder_std
small_disorder_std = 0.0
#
#
#SYMMETRY OPTION: 
#Default: YES - If YES then the code calculates the atomic positions only in the positive (x>0, y>0) quadrant, 
#and then assumes that the whole physical system (atoms+input light) is symmetric for x->-x and y->-y.
#When applicable, this options allows to highly speed up the code, and consume less memory
const mirror_symmetry_option    =      ["YES" ; "NO"][1]
#
#
#REPETITION NUMBER
#The code will calculate the atomic response a number of times given by n_repetitions
#For disordered geometries, each repetition will sample different atomic positions
#Similarly, if inhom_broad_std>0 (see below), each repetition will sample a different set 
#of random resonance frequencies for the atoms
n_repetitions = 1
#
#
#
#
#
#
#############################################################################################################
################## OVERALL PARAMETERS #######################################################################
#############################################################################################################
#
#
#INITIALIZATION PARAMETERS:
#Adds a name to the simulation, if the Julia file is launched without passing through "Bash_Launcher.sh"
#If the simulation is launched via "Bash_Launcher.sh" then the following name is ignored
const name_simulation = "DEFAULT"
#Maximum RAM available for the computation (in GB)
RAM_GB_max = 450
#
#
#
#
#
#
#
#############################################################################################################
################## PHYSICAL PARAMETERS ######################################################################
#############################################################################################################
#
#
#ATOMIC SYSTEM & ENVIRONMENT:
#Spatial orientation of the atomic dipole-matrix elements
dipoles_polarization       =  [1.0 ; 0.0 ; 0.0] 
#Inelastic decay rate Gamma' of the atoms, in units of Gamma0
const gamma_prime          =  5.75 
#Standard deviation of the Gaussian distribution of inhomogeneous broadening of the atomic resonance frequencies (in units of Gamma0). 
#If =0.0 then no inhomogeneous broadening is added.
const inhom_broad_std      =  0.0
#
#
#INPUT GAUSSIAN BEAM:
#Detuning of the light frequency with respect to the atomic resonance frequency, in units of Gamma0,
#i.e. Delta = (omega_laser - omega0) / Gamma0
#Many values can be set, e.g. [0.0 ; 1.0 ;2.0 ;-10]
#Or as a range, by writing:
#laser_detunings  = range(minimum_value,stop=maximum_value,length=number_of_points)
laser_detunings            =  [0.0] 
#Beam waist of the input Gaussan beam
const w0                   =  4.0
#Direction of the input Gaussian beam. Default: [0.0 ; 0.0 ; 1.0] 
laser_direction            =  [0.0 ; 0.0 ; 1.0] 
#Polarization of the input light. Default: [1.0 ; 0.0 ; 0.0] 
field_polarization         =  [1.0 ; 0.0 ; 0.0] 
#
#
#
#
#
#
#
#############################################################################################################
################## PROBE PARAMETERS #########################################################################
#############################################################################################################
#
#
#Settings for the probe points in the XY plane. 
#If probeXY_option = "NO" then this probe is not computed
probeXY_option        =   ["YES" ; "NO"][1]
#Distance (in the z direction) between the atoms, centered at (x,y,z)=(0,0,0), and the XY probe plane
probeXY_z             =   20
#Number of different probe coordinates in the x direction 
probeXY_points_x      =   250
#Number of different probe coordinates in the y direction 
probeXY_points_y      =   250
#Size of the probe rectangle in the x direction, in the form [min_x ; max_x]
probeXY_range_x       =   [-1 ; 1].*10
#Size of the probe rectangle in the y direction, in the form [min_y ; max_y]
probeXY_range_y       =   [-1 ; 1].*10
#
#
#Settings for the probe points in the YZ plane. 
#If probeYZ_option = "NO" then this probe is not computed
#The meaning of the variables is analogous to above
probeYZ_option        =   ["YES" ; "NO"][1]
probeYZ_x             =   0
probeYZ_points_y      =   200
probeYZ_points_z      =   500
probeYZ_range_y       =   [-1 ; 1].*10
probeYZ_range_z       =   [-1 ; 1].*20
#
#
#Settings for the probe points in the XZ plane. 
#If probeXZ_option = "NO" then this probe is not computed
#The meaning of the variables is analogous to above
probeXZ_option        =   ["YES" ; "NO"][1]
probeXZ_y             =   0
probeXZ_points_x      =   200
probeXZ_points_z      =   500
probeXZ_range_x       =   [-1 ; 1].*10
probeXZ_range_z       =   [-1 ; 1].*20
#
#
#Settings for the probe points in the plane perpendicular to a custom vector probePlane_v1_vec. 
#If probePLANE_option = "NO" then this probe is not computed
#The meaning of the variables is analogous to above
probePLANE_option     =   ["YES" ; "NO"][2]
probePlane_v3_vec     =   [1.0 ; 1.0 ; 1.0]./sqrt(3)
probePlane_v3_value   =   0
probePlane_points_v1  =   50
probePlane_points_v2  =   50
probePlane_range_v1   =   [-1 ; 1].*10
probePlane_range_v2   =   [-1 ; 1].*10
#
#
#Settings for the probe points in a sphere or hemisphere surrounding the atomic cloud. If probeSPHERE_option = "NONE" then this probe is not computedd
probeSPHERE_option    =   ["NONE" ; "FULL_SPHERE" ; "FORWARD_HEMISPHERE" ; "BACKWARD_HEMISPHERE"][1]
probeSphere_radius    =   2
probeSphere_points    =   10000
#
#
#If target_beam_option=="YES" then the code calculates the projection of the output field onto a target Gaussian beam.
#This has the same direction and polarization of the input Gaussian beam, but can 
#have different waist w0_target and different focal point z0_target. 
const target_beam_option    =   ["YES" ; "NO"][1] 
w0_target                   =   2*w0
z0_target                   =   10
#If norm_target_option == "YES" then this target beam is multiplied by (w0/w0_target)*exp(im*k0*z0_target), 
#thus preserving the phase and power of the input beam
normalize_target_option     =   ["YES" ; "NO"][1] 
#
#
#
#
#
#
#
#############################################################################################################
################## METALENS SETTINGS ########################################################################
#############################################################################################################
#
#
#Settings that are relevant only if the geometry is set to METALENS
if geometry_settings == "METALENS"
    #
    #
    #METALENS PARAMETERS:
    #Focal length f
    const focal_length              =    20
    #Total radius of the metalens
    const r_lens                	=    10.5
    #Width of each disk composing the metalens, i.e. r_(j+1) - r_j
    const disks_width               =    7.0/10.5
    #Buffer zone parameter 0<=buffer_smooth<=0.5 (dimensionless fraction)
    const buffer_smooth     	    =    0.2
    #Value of the phase shift at the center of the metalens, i.e. phi_0
    const phase_shift               =    1.0775
    #
    #
    #METALENS OPTIONS:
    #
    #Default: YES - If YES then the code overwrites some of the probe parameters with the default probe parameters for the metalens.
    #Specifically:
    #probeXY_z -> focal distance from the lens
    #probeYZ_x = probeXZ_y = 0.0
    #The default sizes of the probe rectangles are specified in "initialization.jl"
    #The number of probe points remains unchanged
    const default_probe_option      =      ["YES" ; "NO"][1] 
    #
    #Default: YES - If YES then (w0_target, z0_target) are overwritten by the default option.
    #It consists of an ideal Gaussian beam expected after a perfect lens with the given focal_length 
    #Similarly normalize_target_option are overwritten to YES
    const default_target_option     =      ["YES" ; "NO"][1]       
    #
    #
    #
    #The following parameters allow to modify the algorithm by which the atomic metalens is constructed
    #
    #Default: YES - If NO  then xi_z will smoothly vary as function of r, inside any ring (including the buffer zone). If YES then xi_z = xi_z^j is kept fixed in any j-th ring
    const z_fixed_option            =      ["YES" ; "NO"][1] 
    #
    #Default: YES - If NO  then xi_z will smoothly vary as function of r, but only inside the buffer zones
    const z_fixed_buffer_option     =      ["YES" ; "NO"][1] 
    #
    #Default: YES - If YES then the average phase of the i-th ring is given by phi((r_{i-1}+r_{i})/2), otherwise it excludes the buffer zone
    const phase_center_ring_option  =      ["YES" ; "NO"][1] 
    #
end
#
#
#
#
#
#
#
#############################################################################################################
################## ARRAYS SETTINGS ##########################################################################
#############################################################################################################
#
#
#Settings that are relevant only if the geometry is set to ARRAYS
if geometry_settings == "ARRAYS"
    #
    #PHYSICAL SETTINGS OF THE ARRAY:
    #
    #Number of 2D arrays in series
    array_n_layers = 1
    #Lattice constants of the 2D array
    array_xi_x = 0.1
    array_xi_y = 0.1
    #Distance between the arrays
    array_xi_z = 0.3
    #Size of the array
    array_size_x = [-2.0;2.0]
    array_size_y = [-2.0;2.0]
    #
    #
    #OPTIONS:
    #
    #Default: NO - If YES then the laser detunings (saved laser_detunings) are interpreted to be
    #in units of the cooperative decay rate Gamma_coop, rather than the natural linewidth Gamma0, i.e.
    #laser_detunings ./= Gamma_coop
    const array_gamma_coop_option = ["YES" ; "NO"][2]
    #
    #Default: YES - If YES then the laser detunings (saved laser_detunings) are shifted so 
    #that the zero corresponds to the cooperative resonance frequency omega_coop of the array, i.e.
    #laser_detunings .= laser_detunings .- omega_coop 
    const array_omega_coop_option = ["YES" ; "NO"][1]
    #
    #Meaning of the detuning for values of (array_gamma_coop_option,array_omega_coop_option)
    #(NO,  NO ) -> Delta = (omega-omega_0)/Gamma0 (i.e. laser_detunings= 1 means 1*Gamma_0 shift from omega_0)
    #(YES, NO ) -> Delta = (omega-omega_0)/Gamma_coop
    #(NO,  YES) -> Delta = (omega-omega_coop)/Gamma0
    #(YES, YES) -> Delta = (omega-omega_coop)/Gamma_coop (i.e. laser_detunings= 1 means 1*Gamma_coop shift from omega_coop)
    #
end
#
#
#
#
#
#
#
#############################################################################################################
################## DISORDERED POSITIONS SETTINGS ############################################################
#############################################################################################################
#
#
#Settings that are relevant only for DISORDERED_SPHERE, DISORDERED_CYLINDER or DISORDERED_CUBOID
if geometry_settings[1:3]=="DIS"
    #
    #Atomic density, in units of lambda0^3
    dis_atomic_density = 0.2
    #
    #Option for a sphere centered in [0.0; 0.0; 0.0]
    if geometry_settings=="DISORDERED_SPHERE"
        #Radius of the sphere 
        dis_r_sphere = 1.0
    end
    #
    #Options for a cylinder with main axis in the z direction
    if geometry_settings=="DISORDERED_CYLINDER"
        #Radial coordinate of the cylinder 
        dis_r_disk = 2.0
        #Length of the main axis (i.e. thickness of the cylinder)
        dis_z_length = 0.2
    end
    #
    #Options for a cuboid
    if geometry_settings=="DISORDERED_CUBOID"
        #Lengths of the three axes
        dis_x_dim = 2.0
        dis_y_dim = 2.0
        dis_z_dim = 2.0
    end
    #
end
#
#
#
#
#
#
#
#############################################################################################################
################## CUSTOM POSITIONS SETTINGS ################################################################
#############################################################################################################
#
#
#Settings that are relevant only when geometry_settings is set to CUSTOM_POSITIONS
if geometry_settings=="CUSTOM_POSITIONS"
    custom_pos_folder   = "Data_Input"
    custom_pos_file     = "custom_pos.h5"
    custom_pos_variable = "custom_pos"
end