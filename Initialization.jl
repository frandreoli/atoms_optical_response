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
const results_folder_name = "Results"
#
#Uncomment the following line only if the code seems to be leaking RAM memory
#It will be forced to print a line any time the unused memory gets emptied
#These messages will appear in the form "GC: ..."
#GC.enable_logging(true)
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
#Consistency of the repetition number
rep_warning = "The number of repetition n_repetitions must be a positive, integer number.\nSetting n_repetitions=1."
if typeof(n_repetitions)!==Int64
    @warn rep_warning
    n_repetitions = 1
elseif n_repetitions<1
    @warn rep_warning
    n_repetitions = 1
end
if (geometry_settings[1:3]!="DIS" && abs(inhom_broad_std)<ZERO_THRESHOLD)&& n_repetitions>1
    @warn "Neither the atomic positions nor the resonance frequencies are randomly chosen.\nThere is no reason to repeat the simulation multiple times.\nSetting n_repetitions=1."
    n_repetitions = 1
end
#
#
#Consistency of defects_fraction with the rest of options
if geometry_settings[1:3]=="DIS"
    @warn "For disordered geometries, punching random defects is redundant.\nSetting defects_fraction=0."
    defects_fraction = 0.0
end
if !isnumeric(defects_fraction) || defects_fraction<0.0 || defects_fraction>1.0
    @warn "The fraction of defects (i.e. defects_fraction) must be a positive number, lower than unity.\nSetting defects_fraction=0."
    defects_fraction = 0.0
end
if mirror_symmetry_option=="YES" && defects_fraction>0.0
    @warn "The mirror symmetry is active and the fraction of (random) defects is non-null.\nThe positions of the defects will be randomized only in the x>0, y>0 quadrant, while they will be symmetric for x->-x and y->-y."
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
time_start=time()
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
time_atomic_pos=time()
#
#
#Punching defects (for ordered geometries, if specified)
if defects_fraction>0.0
    (r_atoms,n_atoms) = defect_punching(r_atoms,n_atoms, defects_fraction)
end
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
if mirror_symmetry_option == "NO" 
    file_name="_nAtoms"*string(n_atoms)
elseif geometry_settings[1:3]!="DIS"
    positive_points = count((r_atoms[:,1].>0.0).*(r_atoms[:,2].>0.0))
    central_points  = count((abs.(r_atoms[:,1]).<ZERO_THRESHOLD).*(abs.(r_atoms[:,2]).<ZERO_THRESHOLD))
    file_name="_nAtoms"*string(4*positive_points + 2*(n_atoms-(positive_points+central_points)) + central_points)
else
    file_name="_nAtoms"*string(4*n_atoms)
end
#
file_name*="_w0"*string(w0/lambda0)[1:min(length(string(w0/lambda0)),3)]
#
n_bulk!=1.0         ?  file_name*="_n"*string(n_bulk)[1:min(3,length(string(n_bulk)))]               : nothing
gamma_prime>0      ?  file_name*="_gPr"*string(gamma_prime)[1:min(5,length(string(gamma_prime)))] : nothing
inhom_broad_std>0  ?  file_name*="_inhom"*string(inhom_broad_std)                                 : nothing
#
#Only if a metalens is computed
if geometry_settings == "METALENS" 
    file_name*="_r"*string(r_lens/lambda0)[1:min(length(string(r_lens/lambda0)),3)]
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
println("\n\nTime: ",now(),"\nStarting evaluation of ", @__FILE__)
println("Output file name: ",file_name,"\n")
final_path_name="Data_Output/"*file_name*"/"
mkpath(final_path_name)
mkpath(final_path_name*"/"*results_folder_name)
if pos_save_option=="YES" && geometry_settings[1:3]!="DIS"
    #Adding a dummy dimension to uniform the data analysis with the atomic_positions for disordered systems
    h5write_multiple(final_path_name*"atomic_positions", ("r_atoms", add_dimension(r_atoms)) ; open_option="w") 
    println("Atomic positions created in                   ", dig_cut(time_atomic_pos-time_start)," seconds")
end
#
#Saving data files with the settings of the simulation
h5write_multiple(final_path_name*"options", ("pos_save_option", pos_save_option) , ("geometry_settings", geometry_settings) ; open_option="w")
h5write_multiple(final_path_name*"options", ("probeXY_option", probeXY_option) , ("probeYZ_option", probeYZ_option) , ("probeXZ_option", probeXZ_option) , ("probePlane_option", probePlane_option) , ("probeSphere_option", probeSphere_option); open_option="w")
h5write_multiple(final_path_name*"settings", ("lambda0", lambda0) , ("n_bulk",n_bulk) , ("w0", w0) , ("gamma_prime", gamma_prime) , ("inhom_broad_std", inhom_broad_std); open_option="w")
#
#Saving data files with the settings specific of the atomic metalens
if geometry_settings == "METALENS" 
    n_phase_disks_to_save = length(collect(0.0:disks_width:r_lens))-1
    h5write_multiple(final_path_name*"settings_metalens",  ("n_phase_disks_to_save", n_phase_disks_to_save), ("focal_point", focal_point) , ("r_lens",r_lens) , ("buffer", buffer_smooth) , ("disks_width", disks_width) ; open_option="w")
    h5write_multiple(final_path_name*"settings_metalens",  ("phase_array",   phase_array) , ("phase_range_theo", phase_range_theo); open_option="w")
    h5write_multiple(final_path_name*"settings_metalens",  ("lens_disks_r",  lens_disks_r) ; open_option="w")
    h5write_multiple(final_path_name*"settings_metalens",  ("lattice_array", lattice_array) ; open_option="w")
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
#Initializes the data files, or overwrites them if already existing
t_and_r_h5 = h5open(final_path_name*"/"*results_folder_name*"/"*"t_and_r.h5", "w")
close(t_and_r_h5)
#
if (true in (x->x=="YES").([probeXY_option ; probeYZ_option ; probeXZ_option ; probePlane_option ])) || probeSphere_option!="NONE"
    probe_pos_h5=h5open(final_path_name*"/"*results_folder_name*"/"*"probe_positions.h5", "w")
    close(probe_pos_h5)
    probe_field_h5=h5open(final_path_name*"/"*results_folder_name*"/"*"probe_field.h5", "w")
    close(probe_field_h5)
end
#
GC.gc()
#
#Main computation
@time for index_repetition in 1:n_repetitions
    if geometry_settings[1:3]=="DIS"
        #TBA!!!
    end
    #
    GC.gc()
    #
    println("\nStarting the repetition ", index_repetition,"/",n_repetitions,".")
    println("| Current available memory:                   ", dig_cut((Sys.free_memory() / 2^20)/1024), " (GB)")
    @time CD_main(r_atoms, n_atoms, w0, k0, laser_direction, laser_detunings, dipoles_polarization, field_polarization, w0_target, z0_target)
    #
    index_repetition==n_repetitions ? println("\nCore evaluation completed. Performance: ") : nothing
end
#
println("\nEvaluation successfully completed. \nTotal evaluation time:                        ", time()-time_start,"\n")

 