#############################################################################################################
################## INITIAL OPERATIONS #######################################################################
#############################################################################################################
#Designed for Julia 1.6
using LinearAlgebra, Dates, HDF5, Random, Distributed
Random.seed!()
if nworkers()==1 addprocs() end
BLAS.set_num_threads(nworkers())
const ZERO_THRESHOLD = 10^(-12)
const results_folder_name = "Results"
#For large-scale numerics Float32 is preferable compared to the default Float64
const TN = [Float32 ; Float64][1]
#
include("Settings.jl")
include("Core - Functions.jl")
include("Core - Pos - Metalens.jl")
include("Core - Pos - Array.jl")
include("Core - Evaluation.jl")
include("Core - Warnings.jl")
#
#Uncomment the following line only if the code seems to be leaking RAM memory or spending too much time in 
#garbage collection.
#You will force the code to print a line any time the unused memory gets emptied.
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
################## DEFINITION OF THE ATOMIC POSITIONS #######################################################
#############################################################################################################
#
#
#Definition of the wavevector (we recall that all lengths are in units of lambda0)
const k0 = 2.0*pi
#
#
time_start=time()
#
#Defines the settings and the atomic positions for an atomic metalens
if geometry_settings == "METALENS"
    (r_atoms, n_atoms, phase_array,lens_disks_r,phase_range_theo, lattice_array) = metalens_creation(r_lens, focal_length, disks_width,buffer_smooth,phase_shift)
    #
    if default_target_option=="YES"
        (w0_target, z0_target) = ideal_beam(w0, k0, focal_length)
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
            probeYZ_x = 0.0
            probeYZ_range_y = [-1 ; 1].*(r_lens*probe_range_factor)
            probeYZ_range_z = [-1 ; 2].*(3.0*z0_target)
        end
        if probeXZ_option == "YES" 
            probe_range_factor = 1.05
            probeXZ_y = 0.0
            probeXZ_range_x = [-1 ; 1].*(r_lens*probe_range_factor)
            probeXZ_range_z = [-1 ; 2].*(3.0*z0_target)
        end
    end
end
#
#Defines the settings and the atomic positions for a series of atomic arrays
if geometry_settings == "ARRAYS"
    #
    #Calculating the cooperative rates
    time_coop_tic = time()
    (omega_coop, Gamma_coop) = coop_values_function(array_xi_x,array_xi_y,laser_direction[1], laser_direction[2], dipoles_polarization)
    time_coop = time() - time_coop_tic
    #
    #Re-scaling the detuning if the option is on
    if array_gamma_coop_option=="YES"
        laser_detunings ./= Gamma_coop
    end
    #Shifting the detuning if the option is on
    if array_omega_coop_option=="YES"
        laser_detunings .= laser_detunings .- omega_coop 
    end
    #
    #Creation of the atomic positions
    (r_atoms, n_atoms) = arrays_creation(array_n_layers, array_xi_x,array_xi_y,array_xi_z,array_size_x,array_size_y)
    #
end
#
time_atomic_pos=time()
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
    file_name="_nAtoms"*string(Int(round(n_atoms*(1-defects_fraction))))
elseif geometry_settings[1:3]!="DIS"
    positive_points = count((r_atoms[:,1].>0.0).*(r_atoms[:,2].>0.0))
    central_points  = count((abs.(r_atoms[:,1]).<ZERO_THRESHOLD).*(abs.(r_atoms[:,2]).<ZERO_THRESHOLD))
    n_atoms_estimated = 4*positive_points + 2*(n_atoms-(positive_points+central_points)) + central_points
    n_atoms_estimated = Int(round(n_atoms_estimated*(1-defects_fraction)))
    file_name="_nAtoms"*string(n_atoms_estimated)
else
    file_name="_nAtoms"*string(4*Int(round(n_atoms*(1-defects_fraction))))
end
#
if defects_fraction>0.0
    file_name*="_defects"
end
#
if small_disorder_std>0.0
    file_name*="_posShifts"
end
#
file_name*="_w0"*string(w0)[1:min(length(string(w0)),3)]
#
gamma_prime>0      ?  file_name*="_gPr"*string(gamma_prime)[1:min(5,length(string(gamma_prime)))] : nothing
inhom_broad_std>0  ?  file_name*="_inhom"*string(inhom_broad_std)                                 : nothing
#
#Only if a metalens is computed
if geometry_settings == "METALENS" 
    file_name*="_r"*dig_cut(r_lens,3)
    file_name*="_f"*dig_cut(focal_length,3)
    file_name*="_widths"*dig_cut(disks_width,4)
    file_name*="_phase"*dig_cut(phase_shift,5)
    file_name*="_buffer"*dig_cut(buffer_smooth,4)
end
#
#Only if atomic arrays are computed
if geometry_settings == "ARRAYS" 
    file_name*="_nArrays"*string(array_n_layers)
    if array_xi_x==array_xi_y
        file_name*="_xi"*dig_cut(array_xi_x,3)
    else
        file_name*="_xi_x"*dig_cut(array_xi_x,3)
        file_name*="_xi_y"*dig_cut(array_xi_y,3)
    end
    #
    if array_n_layers>1
        file_name*="_xi_z"*dig_cut(array_xi_z,3)
    end
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
#
if geometry_settings == "ARRAYS" 
    println("Cooperative frequency and rate computed in    ", dig_cut(time_coop)," seconds" )
end
#
#Saving the positions, but only if no disorder is present
if pos_save_option=="YES" && no_randomness_option
    #Adding a dummy dimension to uniform the data analysis with the atomic_positions for disordered systems
    h5write_multiple(final_path_name*"atomic_positions", ("r_atoms", add_dimension(r_atoms)) ; open_option="w") 
    println("Atomic positions created in                   ", dig_cut(time_atomic_pos-time_start)," seconds")
end
#
#Saving data files with the settings of the simulation
h5write_multiple(final_path_name*"options", ("pos_save_option", pos_save_option) , ("geometry_settings", geometry_settings) ,("target_beam_option",target_beam_option) ; open_option="w")
h5write_multiple(final_path_name*"options", ("probeXY_option", probeXY_option) , ("probeYZ_option", probeYZ_option) , ("probeXZ_option", probeXZ_option) , ("probePLANE_option", probePLANE_option) , ("probeSPHERE_option", probeSPHERE_option))
h5write_multiple(final_path_name*"options", ("mirror_symmetry_option",mirror_symmetry_option))
h5write_multiple(final_path_name*"settings", ("w0", w0) , ("gamma_prime", gamma_prime) , ("inhom_broad_std", inhom_broad_std); open_option="w")
h5write_multiple(final_path_name*"settings", ("laser_detunings",laser_detunings), ("laser_direction",laser_direction), ("field_polarization",field_polarization) ,("defects_fraction",defects_fraction))
h5write_multiple(final_path_name*"settings", ("n_repetitions",n_repetitions),("probePlane_vec",probePlane_v3_vec))
h5write_multiple(final_path_name*"settings", ("dipoles_polarization",dipoles_polarization))
#
if target_beam_option=="YES"
    h5write_multiple(final_path_name*"settings", ("w0_target",w0_target),("z0_target",z0_target) )
    h5write_multiple(final_path_name*"options", ("normalize_target_option", normalize_target_option) )
end
#
#Saving data files with the settings specific of the atomic metalens
if geometry_settings == "METALENS" 
    n_phase_disks = length(collect(0.0:disks_width:r_lens))-1
    h5write_multiple(final_path_name*"settings_metalens",  ("n_phase_disks", n_phase_disks), ("focal_length", focal_length) , ("r_lens",r_lens) , ("buffer", buffer_smooth) , ("disks_width", disks_width) ; open_option="w")
    h5write_multiple(final_path_name*"settings_metalens",  ("phase_array",   phase_array) , ("phase_range_theo", phase_range_theo))
    h5write_multiple(final_path_name*"settings_metalens",  ("lens_disks_r",  lens_disks_r) , ("phase_shift",phase_shift))
    h5write_multiple(final_path_name*"settings_metalens",  ("lattice_array", lattice_array))
    h5write_multiple(final_path_name*"options_metalens",   ("z_fixed_option", z_fixed_option) , ("z_fixed_buffer_option",z_fixed_buffer_option),("phase_center_ring_option",phase_center_ring_option),("default_probe_option",default_probe_option),("default_target_option",default_target_option) ; open_option="w")
end
#
#Saving data files with the settings specific of the atomic arrays
if geometry_settings == "ARRAYS" 
    h5write_multiple(final_path_name*"settings_arrays", ("array_n_layers", array_n_layers), ("array_xi_x",array_xi_x), ("array_xi_y",array_xi_y), ("array_xi_z",array_xi_z)  ; open_option="w")
    h5write_multiple(final_path_name*"settings_arrays", ("array_size_x", array_size_x), ("array_size_y",array_size_y) )
    h5write_multiple(final_path_name*"coop_arrays",     ("omega_coop", omega_coop), ("Gamma_coop", Gamma_coop)  ; open_option="w")    
    h5write_multiple(final_path_name*"options_arrays",  ("array_gamma_coop_option", array_gamma_coop_option) , ("array_omega_coop_option",array_omega_coop_option) ; open_option="w")
end
#
#Checking if the RAM estimated for this simulation exceed the threshold set by the user
tot_probe_points = 0
probeXY_option=="YES"     ? tot_probe_points+=probeXY_points_x*probeXY_points_y           : nothing
probeYZ_option=="YES"     ? tot_probe_points+=probeYZ_points_y*probeYZ_points_z           : nothing
probeXZ_option=="YES"     ? tot_probe_points+=probeXZ_points_x*probeXZ_points_z           : nothing
probePLANE_option=="YES"  ? tot_probe_points+=probePlane_points_v1*probePlane_points_v2   : nothing
probeSPHERE_option=="YES" ? tot_probe_points+=probeSphere_points                          : nothing
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
if (true in (x->x=="YES").([probeXY_option ; probeYZ_option ; probeXZ_option ; probePLANE_option ])) || probeSPHERE_option!="NONE"
    probe_pos_h5=h5open(final_path_name*"/"*results_folder_name*"/"*"probe_positions.h5", "w")
    close(probe_pos_h5)
    probe_field_h5=h5open(final_path_name*"/"*results_folder_name*"/"*"probe_field.h5", "w")
    close(probe_field_h5)
end
#
GC.gc()
#
#Main computation
performance=@timed for index_repetition in 1:n_repetitions
    #
    if geometry_settings[1:3]=="DIS"
        #TBA!!!
    else
        (r_atoms_here ,n_atoms_here) = (r_atoms[:,:] , n_atoms)
    end
    #
    #Punching defects (only for ordered geometries, and if requested)
    if defects_fraction>0.0
        (r_atoms_here,n_atoms_here) = defect_punching(r_atoms_here,n_atoms_here, defects_fraction)
    end
    #
    #Randomly shifting the positions (only for ordered geometries, and if requested)
    if small_disorder_std>0 
        r_atoms_here = r_atoms_here[:,:].+(randn(MersenneTwister(), Float64, (n_atoms_here,3)).*small_disorder_std)
        n_atoms_here = n_atoms
    end
    #
    #Saving the new atomic positions
    if pos_save_option=="YES" && (!no_randomness_option)
        if index_repetition==1
            h5write_multiple(final_path_name*"atomic_positions", ("r_atoms", add_dimension(r_atoms_here)) ; open_option="w"  )
        else
            h5write_append(final_path_name*"atomic_positions", add_dimension(r_atoms_here),  "r_atoms" )
        end
    end
    #
    GC.gc()
    #
    println("\nStarting the repetition ", index_repetition,"/",n_repetitions,".")
    println("| Current available RAM:                      ", dig_cut((Sys.free_memory() / 2^20)/1024), " (GB)")
    @time CD_main(r_atoms_here, n_atoms_here, w0, k0, laser_direction, laser_detunings, dipoles_polarization, field_polarization, w0_target, z0_target)
    #
    index_repetition==n_repetitions ? println("\nCore evaluation completed. Performance of the core code: ") : nothing
end
println("- Total core-computation time:                ", dig_cut(performance[2]-performance[4])," seconds" )
println("- Total garbage-collection time:              ", dig_cut(performance[4])," seconds" )
#
println("\nEvaluation successfully completed. \n- Total running time:                         ", dig_cut(time()-time_start)," seconds","\n")