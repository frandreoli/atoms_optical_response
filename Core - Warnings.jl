#############################################################################################################
################## WARNINGS #################################################################################
#############################################################################################################
#
#
println("\n")
#
#
#Converting all directions to real vectors
if abs(conj_norm(imag.(laser_direction)))>ZERO_THRESHOLD
    laser_direction = real.(laser_direction)
    @warn "The direction of the input field, i.e. laser_direction, must be a real vector.\n Neglecting the imaginary parts."
end
if abs(conj_norm(imag.(probePlane_v3_vec)))>ZERO_THRESHOLD
    probePlane_v3_vec = real.(probePlane_v3_vec)
    @warn "The normal direction of the custom probe plane, i.e. probePlane_v3_vec, must be a real vector.\n Neglecting the imaginary parts."
end
#
#
#Checking the normalization of the direction vectors
if probePLANE_option=="YES" && abs(conj_norm(probePlane_v3_vec)-1.0)>ZERO_THRESHOLD
    probePlane_v3_vec=conj_normalize(probePlane_v3_vec)
    @warn "The normal direction of the custom probe plane, i.e. probePlane_v3_vec, must be normalized to 1.\nNormalizing it to 1."
end
if abs(conj_norm(laser_direction)-1.0)>ZERO_THRESHOLD
    laser_direction=conj_normalize(laser_direction)
    @warn "The direction of input beam, i.e. laser_direction, must be normalized to 1.\nNormalizing it to 1."
end
#
#
#Symmetry consistency
if abs(conj_norm(laser_direction.-[0.0 ; 0.0 ; 1.0]))>ZERO_THRESHOLD && mirror_symmetry_option=="YES"
    mirror_symmetry_option = "NO"
    @warn "The input field is not symmetric for x->-x and y->-y. Setting mirror_symmetry_option to NO."
end
#
#
#Transverse field consistency
field_polarization       = conj_normalize(field_polarization)
longitudinal_component   = conj_dot(laser_direction,field_polarization) 
if abs(longitudinal_component) > ZERO_THRESHOLD
    field_polarization   = field_polarization.-(laser_direction.*longitudinal_component)
    field_polarization   = conj_normalize(field_polarization)
    @warn "The field is not a transverse field. Removing from field_polarization its longitudinal component."
end
#
#
#Consistency of the definition of the standard deviation of inhomogenous broadening
#If it is not defined as a number (tested with isa()), the code enters the if and soesn't evaluate the second condition
if !isa(inhom_broad_std,Number) || inhom_broad_std <0 
    @warn "The standard deviation in the Gaussian distribution of inhomogeneous bradening must be a positive number. \nSetting inhom_broad_std=0.0."
    inhom_broad_std = 0.0
end
#
#
#Warns in case the atomic positions should be random, while mirror_symmetry_option=="YES"
if mirror_symmetry_option=="YES" && geometry_settings[1:3]=="DIS"
    @warn "The atomic positions are disordered, but mirror_symmetry_option=='YES'.\nThe atomic positions will be randomly chosen only in the positive (x>=0, y>=0) quadrant, while the other quadrants will be the mirror images."
end
#
#
#Warns in case the atomic positions are given as an input, while mirror_symmetry_option=="YES"
if mirror_symmetry_option=="YES" && geometry_settings=="CUSTOM_POSITIONS"
    @warn "The options mirror_symmetry_option=='YES' and geometry_settings='CUSTOM_POSITIONS' are set.\nThe code will only consider the atomic positions in the positive (x>=0, y>=0) quadrant."
end
#
#
#Consistency of defects_fraction with the rest of options
if !isa(defects_fraction,Number) || defects_fraction<0.0 || defects_fraction>1.0
    @warn "The fraction of defects (i.e. defects_fraction) must be a positive number, lower than unity.\nSetting defects_fraction=0."
    defects_fraction = 0.0
end
if geometry_settings[1:3]=="DIS" && defects_fraction>0.0
    @warn "For disordered geometries, punching random defects is redundant.\nSetting defects_fraction=0."
    defects_fraction = 0.0
end
if mirror_symmetry_option=="YES" && defects_fraction>0.0
    @warn "The mirror symmetry is active and the fraction of (random) defects is non-null.\nThe positions of the defects will be randomized only in the x>0, y>0 quadrant, while they will be symmetric for x->-x and y->-y."
end
#
#
#Consistency of small_disorder_std
if !isa(small_disorder_std,Number) || small_disorder_std<0.0
    @warn "The fraction of defects (i.e. small_disorder_std) must be a positive number.\nSetting small_disorder_std=0."
    small_disorder_std = 0.0
end
if geometry_settings[1:3]=="DIS" && small_disorder_std>0.0
    @warn "For disordered geometries, randomly shifting the atomic positions is redundant.\nSetting small_disorder_std=0."
    small_disorder_std = 0.0
end
if mirror_symmetry_option=="YES" && small_disorder_std>0.0
    @warn "The mirror symmetry is active and the atomic positions are randomly shifted.\nThe atomic positions will be randomly shifte only in the x>0, y>0 quadrant, while they will be symmetric for x->-x and y->-y."
end
#
#
#Consistency of w0
if w0<1
    @warn "For the paraxial approximation to be fully valid, a regime of w0>=1 is preferable."
end
#
#
#Consistency of the ARRAYS settings
if geometry_settings=="ARRAYS"
    #
    lattice_error(xi) = !isa(xi,Number) || xi<0 || abs(imag(xi))>ZERO_THRESHOLD
    #
    if !isa(array_n_layers,Int64) || array_n_layers<1
        @warn "The number of arrays, i.e. array_n_layers, must be a positive integer >0.\nSetting array_n_layers=1."
    end
    #
    if lattice_error(array_xi_x) || lattice_error(array_xi_y) || lattice_error(array_xi_z)
        @error "The lattice constants must be positive, real numbers."
    end
end
#
#
#Metalens consistency
if geometry_settings == "METALENS" 
    #
    #Consistency of the buffer zone definition
    if !isa(buffer_smooth,Number) || (buffer_smooth<0.0 || buffer_smooth>1.0)
        @warn "The buffer zone is ill-defined out of the range 0<=buffer_smooth<=1. Setting buffer_smooth to the closest, valid value."
        !isa(buffer_smooth,Number) || buffer_smooth<0.0 ? buffer_smooth=0.0 : nothing
        buffer_smooth>1.0 ? buffer_smooth=1.0 : nothing
    end
    #Consistency of the disk width
    if !isa(disks_width,Number) || disks_width>r_lens
        @warn "The width of each disk must be a number smaller than the radius of the metalens. Setting disks_width = r_lens."
        disks_width = r_lens
    end
    if disks_width<=0
        @warn "The width of each disk cannot lower or equal zero. Setting disks_width = 0.1*r_lens."
        disks_width = 0.1*r_lens
    end
end
#
#
#Overall randomness in the simulation
no_randomness_option = (geometry_settings[1:3]!="DIS" && abs(inhom_broad_std)<ZERO_THRESHOLD) && abs(defects_fraction)<ZERO_THRESHOLD && abs(small_disorder_std)<ZERO_THRESHOLD
#
#
#Consistency of the repetition number
rep_warning = "The number of repetition n_repetitions must be a positive, integer number.\nSetting n_repetitions=1."
if typeof(n_repetitions)!==Int64
    @warn rep_warning
    n_repetitions = 1
elseif n_repetitions<1
    @warn rep_warning
    n_repetitions = 1
end
if no_randomness_option && n_repetitions>1
    @warn "There are no randomly chosen parameters.\nThere is no reason to repeat the simulation multiple times.\nSetting n_repetitions=1."
    n_repetitions = 1
end
#
#
#Consistency of the symmetry option
if mirror_symmetry_option=="YES" && geometry_settings=="ARRAYS"
    central_pos_x = sum(l_system_x)/2.0
    central_pos_y = sum(l_system_y)/2.0
    if abs(central_pos_x)>ZERO_THRESHOLD || abs(central_pos_y)>ZERO_THRESHOLD
        @warn "The mirror_symmetry_option is incompatible with transversely shifted arrays. \nSetting mirror_symmetry_option=NO."
        mirror_symmetry_option = "NO"
    end
end
#
#
#Consistency of strain option
if strain_option=="CHAIN" && geometry_settings!="CHAIN"
    @warn "The strain_option cannot be set to CHAIN when the geometry_settings is not set to CHAIN. \nSetting strain_option=NONE."
    strain_option = "NONE"
end
#
#
#Checking the normalize_target_option
if geometry_settings == "METALENS"
    if default_target_option=="YES"
        @warn "Due to default_target_option=YES, the option normalize_target_option is changed to YES."
        normalize_target_option = "YES"
    end
end
#
#  
#
#Fixing the options to constants
const mirror_symmetry_option_const        =  mirror_symmetry_option
const target_beam_option_const            =  target_beam_option
const normalize_target_option_const       =  normalize_target_option
if geometry_settings == "METALENS"
    const z_fixed_option_const            =  z_fixed_option
    const z_fixed_buffer_option_const     =  z_fixed_buffer_option
    const phase_center_ring_option_const  =  phase_center_ring_option
end