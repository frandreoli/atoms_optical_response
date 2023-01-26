#############################################################################################################
################## INITIAL OPERATIONS #######################################################################
#############################################################################################################
#Designed for Julia 1.6
using LinearAlgebra, Dates, HDF5, Random, Distributed
include("Core - Functions.jl")
include("Core - Metalens.jl")
include("Core - Evaluation.jl")
Random.seed!()
if nworkers()==1 addprocs() end
BLAS.set_num_threads(nworkers())
const ZERO_THRESHOLD = 10^(-12)
#
#
#
#
#
#
#
#############################################################################################################
################## WARNINGS #################################################################################
#############################################################################################################
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
if probePlane_option=="YES" && abs(conj_norm(probePlane_v3_vec)-1.0)>ZERO_THRESHOLD
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
if inhom_broad_std <0 
    @warn "The standard deviation in the Gaussian distribution of inhomogeneous bradening cannot be negative. \nReplacing inhom_broad_std with |inhom_broad_std|."
    inhom_broad_std = abs(inhom_broad_std)
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
#Metalens consistency
if geometry_settings == "METALENS" 
    #
    #Consistency of the buffer zone definition
    if (buffer_smooth<0.0 || buffer_smooth>1.0)
        @warn "The buffer zone is ill-defined out of the range 0<=buffer_smooth<=1. Choosing the closest, valid value."
        buffer_smooth<0.0 ? buffer_smooth=0.0 : nothing
        buffer_smooth>1.0 ? buffer_smooth=1.0 : nothing
    end
    #Consistency of the disk width
    if disks_width>r_lens
        @warn "The width of each disk cannot be larger than the radius of the metalens. Setting disks_width = r_lens."
        disks_width = r_lens
    end
    if disks_width<=0
        @warn "The width of each disk cannot lower or equal zero. Setting disks_width = 0.1*r_lens."
        disks_width = 0.1*r_lens
    end
end
#
#
#
#
#
#
#
#############################################################################################################
################## DEFINITION OF THE ATOMIC POSITIONS #######################################################
#############################################################################################################
#
#
const time_start=time()
#
#Defines the settings and the atomic positions for an atomic metalens
if geometry_settings == "METALENS"
    (r_atoms, n_atoms, phase_array,lens_disks_r,phase_range_theo, lattice_array) = metalens_creation(r_lens, focal_point, disks_width,buffer_smooth)
    #
    if default_out_beam_option=="YES"
        (w0_target, z0_target) = ideal_beam(w0, k0, focal_point)
        normalize_target_option = "YES"
    end
    #
    #Default settings of the probe points
    if default_probe_option=="YES"
        if probeXY_option == "YES" 
            probe_range_factor = 0.9
            probeXY_z = z0_target
            probeXY_range_x = probeXY_range_y = [-1 ; 1].*(r_lens*probe_range_factor)
        end
        if probeYZ_option == "YES" 
            probe_range_factor = 1.05
            probeYZ_x = 0.0*lambda0
            probeYZ_range_y = [-1 ; 1].*(r_lens*probe_range_factor)
            probeYZ_range_z = [-1 ; 2].*(3.0*z0_target)
        end
        if probeXZ_option == "YES" 
            probe_range_factor = 1.05
            probeXZ_y = 0.0*lambda0
            probeXZ_range_x = [-1 ; 1].*(r_lens*probe_range_factor)
            probeXZ_range_z = [-1 ; 2].*(3.0*z0_target)
        end
    end
end
#
#Defines the settings and the atomic positions for a series of atomic arrays
if geometry_settings == "ARRAYS"
end
#
const time_atomic_pos=time()
#
#
#
#
#
#
#
#############################################################################################################
################## SETTING THE ENVIRONMENT ##################################################################
#############################################################################################################
#
#
#Defines the filename, by adding labels identifying the user choices of parameters for the simulation
length(ARGS)>=1 ? args_checked=ARGS[:] : args_checked=[name_simulation] 
mirror_symmetry_option == "NO" ? file_name="_nAtoms"*string(n_atoms) : file_name="_nAtoms"*string(4*n_atoms)
#
file_name*="_w0"*string(w0/lambda0)[1:min(length(string(w0/lambda0)),3)]
#
nBulk!=1.0         ?  file_name*="_n"*string(nBulk)[1:min(3,length(string(nBulk)))]               : nothing
gamma_prime>0      ?  file_name*="_gPr"*string(gamma_prime)[1:min(5,length(string(gamma_prime)))] : nothing
inhom_broad_std>0  ?  file_name*="_inhom"*string(inhom_broad_std)                                 : nothing
#
#Only if a metalens is computed
if geometry_settings == "METALENS" 
    file_name*="r"*string(r_lens/lambda0)[1:min(length(string(r_lens/lambda0)),3)]
    file_name*="_f"*string(focal_point/lambda0)[1:min(length(string(focal_point/lambda0)),3)]
    file_name*="_widths"*string(disks_width)[1:min(4,length(string(disks_width)))]
    file_name*="_phase"*string(phase_shift)[1:min(5,length(string(phase_shift)))]
    file_name*="_buffer"*string(buffer_smooth)[1:min(4,length(string(buffer_smooth)))]
end
#
mirror_symmetry_option=="YES" ? file_name*="_MIRROR" : nothing
file_name*="_"*string(Dates.today())
file_name=geometry_settings*file_name*"_"*args_checked[1]
#
#
#
#
#
#
#
#############################################################################################################
################## MAIN EVALUATION ##########################################################################
#############################################################################################################
#
#
#Starting the evaluation
println("\nTime: ",now(),"\nStarting evaluation of ", @__FILE__,"\n")
println("Output file name: ",file_name,"\n\n")
final_path_name="Data_Output/"*file_name*"/"
mkpath(final_path_name)
pos_save_option=="YES"  ? h5write_multiple(final_path_name*"pos_atoms", [("r_atoms", r_atoms)])  : nothing
println("Atomic positions created in                 ", time_atomic_pos-time_start)
#
#
#Saving data files with the settings of the simulation
if geometry_settings == "METALENS" 
    n_phase_disks_to_save = length(collect(0.0:disks_width:r_lens))-1
    h5write_multiple(final_path_name*"inputs",        [("inputs",        [lambda0 w0 r_lens focal_point n_phase_disks_to_save nBulk gamma_prime inhom_broad_std] ), ("buffer", buffer_smooth), ("disks_width", disks_width)])
    h5write_multiple(final_path_name*"phase_array",   [("phase_array",   phase_array), ("phase_range_theo", phase_range_theo)                           ])
    h5write_multiple(final_path_name*"lens_disks_r",  [("lens_disks_r",  lens_disks_r                                                                  )])
    h5write_multiple(final_path_name*"lattice_array", [("lattice_array", lattice_array                                                                 )])
end
#
#
#Checking if the RAM estimated for this simulation exceed the threshold set by the user
tot_probe_points = 0
probeXY_option=="YES"     ? tot_probe_points+=probeXY_points_x*probeXY_points_y           : nothing
probeYZ_option=="YES"     ? tot_probe_points+=probeYZ_points_y*probeYZ_points_z           : nothing
probeXZ_option=="YES"     ? tot_probe_points+=probeXZ_points_x*probeXZ_points_z           : nothing
probePlane_option=="YES"  ? tot_probe_points+=probePlane_points_v1*probePlane_points_v2   : nothing
probeSphere_option=="YES" ? tot_probe_points+=probeSphere_points                          : nothing
RAM_GB_estimate = RAM_estimation(TN,tot_probe_points, n_atoms, length(laser_detunings))
if 1.05*RAM_GB_estimate>RAM_GB_max 
    error("Too many atoms: ", n_atoms,"\nThe RAM consumption is estimated as: ", RAM_GB_estimate)
end
#
#
#Main computation
@time CD_main(r_atoms, n_atoms, w0, k0, laser_direction, laser_detunings, dipoles_polarization, field_polarization, w0_target, z0_target )
println("\nEvaluation successfully completed. \nEvaluation time:                            ", time()-time_start,"\n")

