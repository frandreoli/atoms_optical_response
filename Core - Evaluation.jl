#############################################################################################################
################## MAIN DEFINITION ##########################################################################
#############################################################################################################
#
#
#Main Function
function CD_main(r_atoms, n_atoms, w0, k0, laser_direction, laser_detunings, dipoles_polarization, field_polarization, w0_target, z0_target,input_field_function, gamma_prime_func )
	#
	#
	#INITIALIZATION AND COUPLED-DIPOLE SOLUTION:
	#
	#Green's function construction and inversion of the couple-dipole equations. 
	#The resulting dipoles are stored in state_coeff[:,:], where the first index defines the detuning of the input light and the second the atom
	#This is the part requiring most of the computational time/memory
	E_field_in    =  input_field_function.(r_atoms[:,1],  r_atoms[:,2],  r_atoms[:,3])
	E_field_in  .*=  conj_dot(dipoles_polarization,field_polarization)
	gamma_prime_array = gamma_prime_func.(r_atoms[:,1],  r_atoms[:,2],  r_atoms[:,3])
	state_coeff   =  CD_inversion(r_atoms, n_atoms, dipoles_polarization, E_field_in, laser_detunings, gamma_prime_array)
	#
	GC.safepoint()
	#
	#If mirror_symmetry_option=="YES" then each atoms counts as well for its mirrored positions. The variable atoms_mult accounts for that 
	if mirror_symmetry_option=="YES"
		atoms_mult = [2*Int(abs(r_atoms[i,1])>ZERO_THRESHOLD)+2*Int(abs(r_atoms[i,2])>ZERO_THRESHOLD) for i in 1:n_atoms]
	else
		atoms_mult = [1 for i in 1:n_atoms]
	end
	#
	n_detunings = length(laser_detunings)
	#
	#
	#CALCULATES THE TRANSMISSION/REFLECTION INTO INPUT AND TARGET MODES:
	#
	#Calculates the projection onto the same Gaussian mode as the input
	time_t_and_r = time()
	t_in = Array{Complex{TN}}(undef,1,n_detunings)
	r_in = Array{Complex{TN}}(undef,1,n_detunings)
	for i in 1:n_detunings
		(t_in[1,i], r_in[1,i]) = CD_t_r_func(E_field_in, w0, 0.0, state_coeff[i,:],atoms_mult, w0)
	end
	h5write_complex_append(final_path_name*"/"*results_folder_name*"/"*"t_and_r", t_in, "t_in")
	h5write_complex_append(final_path_name*"/"*results_folder_name*"/"*"t_and_r", r_in, "r_in")
	#
	#Computes the transmission by projecting onto the target Gaussian beam
	if target_beam_option == "YES"
		#Calculates the target Gaussian mode at the atomic positions
		E_field_target        = gaussian_beam.(r_atoms[:,1],  r_atoms[:,2],  r_atoms[:,3].-z0_target,  w0_target, w0_target, k0, laser_direction[1], laser_direction[2], laser_direction[3])
		E_field_target      .*= conj_dot(dipoles_polarization,field_polarization)
		if normalize_target_option=="YES"
			E_field_target  .*= (w0/w0_target)*exp(2.0im*pi*z0_target)
		end
		#Calculates the projection onto the target Gaussian mode
		t_target = Array{Complex{TN}}(undef,1,n_detunings)
		r_target = Array{Complex{TN}}(undef,1,n_detunings)
		for i in 1:n_detunings
			(t_target[i], r_target[i]) = CD_t_r_func(E_field_target, w0_target, z0_target, state_coeff[i,:],atoms_mult, w0)
		end
		h5write_complex_append(final_path_name*"/"*results_folder_name*"/"*"t_and_r", t_target, "t_target")
		h5write_complex_append(final_path_name*"/"*results_folder_name*"/"*"t_and_r", r_target, "r_target")
	end
	println("| Transmission and reflection computed in     ", dig_cut(time()- time_t_and_r)," seconds")
	#
	#
	#COMPUTES THE FIELD AT THE PROBE POSITIONS:
	#
	#First probe, plane XY
	if probeXY_option=="YES"
		tot_probe_points = probeXY_points_x*probeXY_points_y
		r_probe =  f_probe_PLANE([0.0;0.0;0.1], probeXY_z, probeXY_points_x, probeXY_points_y, probeXY_range_x, probeXY_range_y)
		CD_output_field_wrap("XY", n_detunings, tot_probe_points, n_atoms, r_probe, r_atoms, state_coeff, field_polarization, dipoles_polarization,input_field_function)
	end
	#
	#Second probe, plane YZ:
	if probeYZ_option=="YES"
		tot_probe_points = probeYZ_points_y*probeYZ_points_z
		r_probe =  f_probe_PLANE([1.0;0.0;0.0], probeYZ_x, probeYZ_points_y, probeYZ_points_z, probeYZ_range_y, probeYZ_range_z)
		CD_output_field_wrap("YZ", n_detunings, tot_probe_points, n_atoms, r_probe, r_atoms, state_coeff, field_polarization, dipoles_polarization,input_field_function)
	end
	#
	#Third probe, plane XZ:
	if probeXZ_option=="YES"
		tot_probe_points = probeXZ_points_x*probeXZ_points_z
		r_probe =  f_probe_PLANE([0.0;1.0;0.0], probeXZ_y, probeXZ_points_x, probeXZ_points_z, probeXZ_range_x, probeXZ_range_z)
		CD_output_field_wrap("XZ", n_detunings, tot_probe_points, n_atoms, r_probe, r_atoms, state_coeff, field_polarization, dipoles_polarization,input_field_function)
	end	
	#
	#Fourth probe, custom plane:
	if probePLANE_option=="YES"
		tot_probe_points = probePlane_points_v1*probePlane_points_v2
		r_probe =  f_probe_PLANE(probePlane_v3_vec, probePlane_v3_value, probePlane_points_v1, probePlane_points_v2, probePlane_range_v1, probePlane_range_v2)
		CD_output_field_wrap("PLANE", n_detunings, tot_probe_points, n_atoms, r_probe, r_atoms, state_coeff, field_polarization, dipoles_polarization,input_field_function)
	end	
	#
	#Fifth probe, sphere:
	if probeSPHERE_option!="NONE"
		r_probe = f_probe_SPHERE(probeSphere_radius,probeSphere_points,probeSPHERE_option)
		CD_output_field_wrap("SPHERE", n_detunings, probeSphere_points, n_atoms, r_probe, r_atoms, state_coeff, field_polarization, dipoles_polarization,input_field_function)
	end	
	#
	println("┕ Performance:")
	#
	return state_coeff
end
#
#
#
#
#
#
#
#############################################################################################################
################## COUPLED-DIPOLE FUNCTIONS #################################################################
#############################################################################################################
#
#
#Function to invert the coupled-dipole equations
function CD_inversion(r_atoms, n_atoms,dipoles_polarization,E_field_in, laser_detunings, gamma_prime_array)
	n_detunings = length(laser_detunings)
	CD_green_matrix  =  Array{Complex{TN},2}(undef, n_atoms, n_atoms)
	time_temp=time()
	CD_initialize!(CD_green_matrix, r_atoms, dipoles_polarization, gamma_prime_array)
    println("| Green's matrix initialized in               ", dig_cut(time()-time_temp)," seconds")
	#
	time_temp=time()
	state_coeff = Array{Complex{TN},2}(undef,n_detunings, n_atoms)
	state_coeff[:]
	for i in 1:n_detunings
		state_coeff[i,:]  = (CD_green_matrix-UniformScaling(laser_detunings[i]))\E_field_in
	end
	println("| SM core evaluation finished in              ", dig_cut(time()-time_temp)," seconds")
	return state_coeff
end
#
#
#Function to fill the green tensor for the SM
function CD_initialize!(CD_green_matrix, r_vecs, p, gamma_prime_array)
    na = length(r_vecs[:,1])
	r_vecs_x = r_vecs[:,1].*k0
	r_vecs_y = r_vecs[:,2].*k0
	r_vecs_z = r_vecs[:,3].*k0
	#Total number of steps
	#If mirror_symmetry_option == "YES" then the matrix is not symmetric anymore for those atoms belonging to one of the axes
	mirror_symmetry_option == "YES" ? n_steps = na^2 : n_steps =Int((na^2 + na)/2)
	#The core part is constructed as a single cycle loop to make the parallelization in threads faster
	#The code is structured not to require the loading of external functions, to boost the parallelization efficiency
	Threads.@threads for index in  1:n_steps
		#Constructing two indices from a single index
		if mirror_symmetry_option == "YES"
			i = Int(floor(  (index-1)/na  )) + 1
			j = Int( index-(i-1)*na)
		else
			i = Int(floor( -(1/2) + sqrt(1/4 + 2*(index - 1 ))   )) + 1
			j = Int(index - (i - 1)*i/2)
		end
		#Inhomogeneous broadening option
		inhom_broad_std>0.0 ? inhom_broad_here=randn()*inhom_broad_std : inhom_broad_here=0.0
		#Constructing the values
		x_i=r_vecs_x[i]
		y_i=r_vecs_y[i]
		x_j=r_vecs_x[j]
		y_j=r_vecs_y[j]
		z = r_vecs_z[i]- r_vecs_z[j]
		to_add=0.0+0.0im
		#Value of the Greens's function
		if i==j
			to_add+=-(0.5im*(1+gamma_prime_array[i])) + inhom_broad_here
		else
			x = x_i-x_j
			y = y_i-y_j
			r =sqrt((x^2)+(y^2)+(z^2))
			cos_theta_square = ((x*p[1])+(y*p[2])+(p[3]*z))^2
			to_add+=-(3*exp(im*r)/(4*r^3)*((r^2 + im*r - 1)*1+(3 - r^2 - 3im*r)*cos_theta_square/r^2))
		end
		#
		#The following steps apply only if mirror_symmetry_option == "YES"
		if mirror_symmetry_option == "YES"
			#Second step
			if abs(x_j)>ZERO_THRESHOLD
				x = x_i+x_j
				y = y_i-y_j
				r =sqrt((x^2)+(y^2)+(z^2))
				cos_theta_square = ((x*p[1])+(y*p[2])+(p[3]*z))^2
				to_add+=-(3*exp(im*r)/(4*r^3)*((r^2 + im*r - 1)*1+(3 - r^2 - 3im*r)*cos_theta_square/r^2))
			end
			#Third step
			if abs(y_j)>ZERO_THRESHOLD
				x = x_i-x_j
				y = y_i+y_j
				r =sqrt((x^2)+(y^2)+(z^2))
				cos_theta_square = ((x*p[1])+(y*p[2])+(p[3]*z))^2
				to_add+=-(3*exp(im*r)/(4*r^3)*((r^2 + im*r - 1)*1+(3 - r^2 - 3im*r)*cos_theta_square/r^2))
			end
			#Fourth step
			if abs(x_j)>ZERO_THRESHOLD && abs(y_j)>ZERO_THRESHOLD
				x = x_i+x_j
				y = y_i+y_j
				r =sqrt((x^2)+(y^2)+(z^2))
				cos_theta_square = ((x*p[1])+(y*p[2])+(p[3]*z))^2
				to_add+=-(3*exp(im*r)/(4*r^3)*((r^2 + im*r - 1)*1+(3 - r^2 - 3im*r)*cos_theta_square/r^2))
			end
		end
		#
		#If mirror_symmetry_option == "YES" then the matrix is not symmetric anymore, due to those atoms belonging to one of the axes
		if mirror_symmetry_option == "YES"
			CD_green_matrix[i, j] = Complex{TN}(to_add) 
		else
			CD_green_matrix[i, j] = CD_green_matrix[j, i] = Complex{TN}(to_add)
		end
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
################## INPUT-OUTPUT FUNCTIONS ###################################################################
#############################################################################################################
#
#
#Function to evaluate the SM transmission coefficient
function CD_t_r_func(E_field_target, w0_target, z0_target, state_coeff,atoms_mult, w0_in)
	E_field_target_conj = conj.(E_field_target)
	t_alpha_term        = (im*(3/(4*pi^2))*(1.0/(w0_in^2)))
	t0                  = (2*pi*w0_in*w0_target)/(pi*(w0_in^2+w0_target^2)+1.0im*z0_target)
    t                   = t0 + (t_alpha_term*sum(E_field_target_conj.*state_coeff.*atoms_mult))
	r                   =      (t_alpha_term*sum(E_field_target.*     state_coeff.*atoms_mult))
    return (t, r)
end
#
#
#Wrapper of the function to evaluate the field at the probe points
function CD_output_field_wrap(name, n_detunings, tot_probe_points, n_atoms, r_probe, r_atoms,  state_coeff, field_polarization, dipoles_polarization,input_field_function)
	time_temp = time()
	E_field_out_probe = Array{Complex{Float64}}(undef, 1,n_detunings, tot_probe_points, 3)
	#Computing the field for the 3 polarizations and for each detuning
	#It stores the final result in the matrix E_field_out_probe, whose first index is the repetition,
	#the second is the detuning, the third is the position of the probe, 
	#and the fourth is the field polarization (in the basis x,y,z)
	for det_index in 1:n_detunings
		E_field_out_probe[1,det_index, :,:] = CD_output_field_func(tot_probe_points, n_atoms, r_probe, r_atoms, state_coeff[det_index,:], field_polarization, dipoles_polarization,input_field_function)
	end
	#
	#Saving the data
	h5write_append(final_path_name*"/"*results_folder_name*"/"*"probe_positions",  r_probe, "probe_pos_"*name)
	h5write_complex_append(final_path_name*"/"*results_folder_name*"/"*"probe_field", E_field_out_probe, "probe_field_"*name)
	#
	length(name)>=2 ? space_add=" "^(6-(length(name)-2)) : space_add=" "^6
	println("| Evaluation of the "*name*" probe finished in"*space_add, dig_cut(time()-time_temp)," seconds")
end
#
#
#Function to evaluate the field at the probe points
function CD_output_field_func(tot_probe_points, n_atoms, r_probe, r_atoms, state_coeff_here, field_polarization, dipoles_polarization,input_field_function)
	G_matrix_probe =  Array{Complex{TN},2}(undef, tot_probe_points, n_atoms)
	E_field_in_at_probe = input_field_function.(r_probe[:,1], r_probe[:,2], r_probe[:,3])
	unit_vec_xyz =([1.0 ; 0.0 ; 0.0],[0.0 ; 1.0 ; 0.0],[0.0 ; 0.0 ; 1.0])
	E_field_out_probe = Array{Complex{TN},2}(undef, tot_probe_points, 3)
	for pol_index = 1:3
		pol_vec = unit_vec_xyz[pol_index]
		CD_initialize_probe!(G_matrix_probe, r_probe, r_atoms, pol_vec, dipoles_polarization)
		E_field_out_probe[:,pol_index] = CD_output_field_core(G_matrix_probe, E_field_in_at_probe.*field_polarization[pol_index], state_coeff_here)
	end
	return E_field_out_probe
end
#
#
#Core function to evaluate the field at the probe points
function CD_output_field_core(CD_green_matrix_probe, E_field_in_probe, state_coeff)
	n_probes=length(E_field_in_probe[:,1])
	E_tot_probe=Array{Complex{Float64}}(undef, n_probes)
	Threads.@threads for i in 1:n_probes
		E_tot_probe[i] = E_field_in_probe[i] + sum(CD_green_matrix_probe[i,:].*state_coeff[:])
	end
    return E_tot_probe
end
#
#
#Function to fill the green tensor for the output
function CD_initialize_probe!(CD_green_matrix_probe, r_vecs1, r_vecs2, p1,p2)
	p_dot     = dot(p1, p2)
    na1       = length(r_vecs1[:,1])
	na2       = length(r_vecs2[:,1])
	r_vecs1_x = r_vecs1[:,1].*k0
	r_vecs1_y = r_vecs1[:,2].*k0
	r_vecs1_z = r_vecs1[:,3].*k0
	r_vecs2_x = r_vecs2[:,1].*k0
	r_vecs2_y = r_vecs2[:,2].*k0
	r_vecs2_z = r_vecs2[:,3].*k0
	n_tot=na1*na2
	#
	#The core part is constructed as a single cycle loop to make the parallelization in threads faster 
	#The code is structured not to require the loading of external functions, to boost the parallelization efficiency
	#Both statement are based on empirical observation
	Threads.@threads for index in 1:n_tot
		#Constructing two indices from a single index
		i = Int(floor(  (index-1)/na2  )) + 1
		j = Int( index-(i-1)*na2)
		#Constructing the distance variables
		z = r_vecs1_z[i]-r_vecs2_z[j]
		x_i=r_vecs1_x[i]
		y_i=r_vecs1_y[i]
		x_j=r_vecs2_x[j]
		y_j=r_vecs2_y[j]
		to_add=0.0+0.0im
		#First step
		x = x_i-x_j
		y = y_i-y_j
		r =sqrt((x^2)+(y^2)+(z^2))
		cos_theta_square=(p1[1]*x+p1[2]*y+p1[3]*z)*(p2[1]*x+p2[2]*y+p2[3]*z)
		to_add+=(3*exp(im*r)/(4*r^3)*((r^2 + im*r - 1)*p_dot+(3 - r^2 - 3im*r)*cos_theta_square/r^2))
		#
		#The following steps apply only if mirror_symmetry_option == "YES"
		if mirror_symmetry_option == "YES"
			#Second step
			if abs(x_j)>ZERO_THRESHOLD
				x = x_i+x_j
				y = y_i-y_j
				r =sqrt((x^2)+(y^2)+(z^2))
				cos_theta_square=(p1[1]*x+p1[2]*y+p1[3]*z)*(p2[1]*x+p2[2]*y+p2[3]*z)
				to_add+=(3*exp(im*r)/(4*r^3)*((r^2 + im*r - 1)*p_dot+(3 - r^2 - 3im*r)*cos_theta_square/r^2))
			end
			#Third step
			if abs(y_j)>ZERO_THRESHOLD
				x = x_i-x_j
				y = y_i+y_j
				r =sqrt((x^2)+(y^2)+(z^2))
				cos_theta_square=(p1[1]*x+p1[2]*y+p1[3]*z)*(p2[1]*x+p2[2]*y+p2[3]*z)
				to_add+=(3*exp(im*r)/(4*r^3)*((r^2 + im*r - 1)*p_dot+(3 - r^2 - 3im*r)*cos_theta_square/r^2))
			end
			#Fourth step
			if abs(x_j)>ZERO_THRESHOLD && abs(y_j)>ZERO_THRESHOLD
				x = x_i+x_j
				y = y_i+y_j
				r =sqrt((x^2)+(y^2)+(z^2))
				cos_theta_square=(p1[1]*x+p1[2]*y+p1[3]*z)*(p2[1]*x+p2[2]*y+p2[3]*z)
				to_add+=(3*exp(im*r)/(4*r^3)*((r^2 + im*r - 1)*p_dot+(3 - r^2 - 3im*r)*cos_theta_square/r^2))
			end
		end
		#Converts the values into TN type to save RAM space
		CD_green_matrix_probe[i, j] = Complex{TN}(to_add)
	end
end