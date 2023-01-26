#############################################################################################################
################## METALENS FUNCTIONS #######################################################################
#############################################################################################################
#
#
#Calculate the waist and focal point of the theoretical, ideal Gaussian beam after a lens with focal_point
#given an initial waist w0 and focal point z0_in
function ideal_beam(w0, k0, focal_point, z0_in=0.0)
	zR_in = (k0*(w0^2))/2.0
	magn  = focal_point/sqrt((z0_in-focal_point)^2+zR_in^2)
	w0_f  = w0*magn
	z0_f  = focal_point+(z0_in-focal_point)*magn^2
	return (w0_f, z0_f)
end
#
#
#
#
#
#
#
#############################################################################################################
################## METALENS ATOMIC POSITIONS ################################################################
#############################################################################################################
#
#
#FUNCTIONS TO CONSTRUCT THE METALENS:
#
#
#Find element of "array" closest to "value"
function nearest_choice(array, value)
    findmin(abs.(array.-value))[2]
end
#
#Main function to construct the atomic metalens
function metalens_creation(r_lens, focal_point, disks_width, buffer_zone)
	#The variable "/PhaseData" in "phase_constants_data.h5" contains a 4xM array.
	#This stores M vectors with 4 elements each, i.e. xi_x^i, xi_y^i, xi_z(xi_x^i,xi_y^i) and phi(xi_x^i,xi_y^i)
	#When the algorithm needs a value of phi_j in the j-th ring, it finds the closest match phi_j=phi(xi_x^i,xi_y^i)
	#And then it infers the corresponding lattice constants
    phase_full_data=h5read("Data_Input/phase_constants_data"*".h5", "/PhaseData")
	phase_full_data[4,:]=shift_phase.(phase_full_data[4,:])
    phi_func(r1,r2) = mod(k0*(focal_point-sqrt(focal_point^2+((r1+r2)/2)^2))+phase_shift,2*pi)
	#Calculates the radii of the rings
    r_range       = collect(0.0:disks_width:r_lens)
	#Calculates the number of rings
	n_phase_rings = length(r_range)-1
	#Calculates the widths of the buffer zones
	buffer_range  = vcat(0.0,r_range[2:end-1].+(r_range[2]*buffer_zone))
	#Decides if including the buffer zones or not in the evaluation of the average phase needed in the ring
	if phase_center_ring_option=="NO"
		phase_range = phi_func.(buffer_range, r_range[2:end])
	else
		phase_range = phi_func.(r_range[1:end-1], r_range[2:end])
	end
	#Constructs the arrays and initializes the variables
    r_atoms       = Array{Float64}(undef,0,3)
	r_atoms_old   = Array{Float64}(undef,0,3)
	r_atoms_new   = Array{Float64}(undef,0,3)
    phase_array   = Array{Float64}(undef,n_phase_rings)
	lattice_array = Array{Float64}(undef,n_phase_rings,3)
	max_lattice_const = 1.0
	old_x_max = 0.0
	old_y_max = 0.0
    for i in 1:n_phase_rings
		#Infers the value of the lattie constants given the required phase in the disk
        current_phase_index = nearest_choice(phase_full_data[4,:], phase_range[i])
        phase_array[i]      = phase_full_data[4,   current_phase_index]
        lattice_array[i,:]  = phase_full_data[1:3, current_phase_index]
		#Checks if (old_x_max, old_y_max) are undefined 
		if (old_x_max, old_y_max) == (0.0,0.0)
			old_x_max = lattice_array[i,1]/2
			old_y_max = lattice_array[i,2]/2
		end
		#Starts constructing the rings
		if lattice_array[i,1]<max_lattice_const && lattice_array[i,2]<max_lattice_const
			#Calculates the atomic positions inside each ring
			(r_atoms_new, old_x_max, old_y_max) = metalens_creation_core(r_lens, lattice_array[i,:],r_range[i+1],r_range[i])
			#The following applies only if there where atoms in the ring calculated at the previous step (i.e. the step cannot be the first one)
			#Otherwise see the "else" below
			if i>1 && length(r_atoms_old[:,1])>0
				#Computes the variation of lattice constants between the previous and the current ring
				a_x_variation = abs(lattice_array[i,1]-lattice_array[i-1,1])/lattice_array[i,1]
				a_y_variation = abs(lattice_array[i,2]-lattice_array[i-1,2])/lattice_array[i,2]
				#Creates the buffer zones (only if buffer_zone>0.0)
				#If both xi_x and xi_y (here a_x and a_y) have varied substantially, no buffer zone is drawn
				if !(a_x_variation>0.01 && a_y_variation>0.01) && buffer_zone>0.0
					(r_atoms_new, old_x_max, old_y_max) = metalens_creation_core(r_lens, lattice_array[i,:],r_range[i+1],buffer_range[i])
					#Decides if xi_x or xi_y has varied more and defines 1 as the one varying more and 2 as that varying less
					if a_x_variation > a_y_variation
						#Changing a_x while a_y remains fixed
						indices_order=[1;2]
					else
						#Changing a_y while a_x remains fixed
						indices_order=[2;1]
					end
					#
					#Finds the lower positions of the current ring 
					(up_nodes_1_start, up_nodes_2_start)   = get_lower(r_atoms_new[:,indices_order[1]],r_atoms_new[:,indices_order[2]])
					#Finds the upper positions of the previous ring
					(low_nodes_1, low_nodes_2) = get_upper(r_atoms_old[:,indices_order[1]],r_atoms_old[:,indices_order[2]])
					#The following code will create atoms connecting the lower positions of the current ring with the upper positions of the previous
					max_low = maximum(low_nodes_1)+lattice_array[i-1,indices_order[1]]*0.75
					select_up = up_nodes_1_start.<max_low
					up_nodes_1 = up_nodes_1_start[select_up]
					up_nodes_2 = up_nodes_2_start[select_up]
					select_up_axis = (up_nodes_1_start.>max_low).*(up_nodes_2_start.>lattice_array[i,indices_order[2]])
					up_nodes_1_to_axis = up_nodes_1_start[select_up_axis]
					up_nodes_2_to_axis = up_nodes_2_start[select_up_axis]
					#
					#Decides if connecting from the lower to the upper ring or viceversa
					if length(up_nodes_1)<=length(low_nodes_1)
						(up_nodes_1_buff, up_nodes_2_buff,low_nodes_1_buff, low_nodes_2_buff) = (up_nodes_1, up_nodes_2, low_nodes_1, low_nodes_2)
						index_up_down = 1
					else
						(up_nodes_1_buff, up_nodes_2_buff,low_nodes_1_buff, low_nodes_2_buff) = (low_nodes_1, low_nodes_2, up_nodes_1, up_nodes_2)
						index_up_down = 2
					end
					#
					#Calculates the atomic positions in the buffer zones. 
					#First for those elements connecting the rings and then for those atoms going straight to the axes (black boxes of Fig.A4)
					(r_atoms_buffer_x, r_atoms_buffer_y)           = (metalens_creation_buffer(up_nodes_1_buff, up_nodes_2_buff,low_nodes_1_buff, low_nodes_2_buff,lattice_array[i,indices_order[2]] ,index_up_down))[indices_order]
					(r_atoms_buffer_x_axis, r_atoms_buffer_y_axis) = (metalens_creation_buffer(up_nodes_1_to_axis, up_nodes_2_to_axis,up_nodes_1_to_axis, [-lattice_array[i,indices_order[2]]/2 for ii in 1:length(up_nodes_2_to_axis)],lattice_array[i,indices_order[2]],1))[indices_order]
					#
					#Decides if the z-component of the positions of the atoms in the buffer zone is fixed or smoothly varies
					if z_fixed_buffer_option=="YES"
						r_atoms_buffer = metalens_creation_buffer_z(r_atoms_buffer_x, r_atoms_buffer_y, lattice_array[i,3])
					else
						r_atoms_buffer = metalens_creation_buffer_z(r_atoms_buffer_x, r_atoms_buffer_y, 0.0)
					end
					#Adds the element to the list of atomic positions
					r_atoms = vcat(r_atoms,r_atoms_buffer)
					#
					#Decides if the z-component of the positions of the atoms in the buffer zone is fixed or smoothly varies. 
					#This applies to those atoms connecting directly to the axes (black boxes of Fig.A4)
					if z_fixed_buffer_option=="YES"
						r_atoms_buffer_axis = metalens_creation_buffer_z(r_atoms_buffer_x_axis, r_atoms_buffer_y_axis, lattice_array[i,3])#!!!!!!!!!!!!!!!!!
					else
						r_atoms_buffer_axis = metalens_creation_buffer_z(r_atoms_buffer_x_axis, r_atoms_buffer_y_axis, 0.0)
					end
					#Adds the element to the list of atomic positions
					#This applies to those atoms connecting directly to the axes (black boxes of Fig.A4)
					r_atoms = vcat(r_atoms,r_atoms_buffer_axis)
				else
					(r_atoms_new, old_x_max, old_y_max) = metalens_creation_core(r_lens, lattice_array[i,:],r_range[i+1],r_range[i]+0.0*(r_range[2]*buffer_zone/3))
				end
			end
			r_atoms = vcat(r_atoms,r_atoms_new)
			r_atoms_old=r_atoms_new[:,:]
		else
			#It applies if it is the first ring created or if the previous ring was empty
			r_atoms_old=Array{Float64}(undef,0,3)
			(old_x_max, old_y_max) == (0.0,0.0)
		end
    end
	#
	#Adds the atomic positions in the other radiants in case mirror_symmetry_option=="NO"
	if mirror_symmetry_option=="NO"
		mirror_selection = (r_atoms[:,1].>0.0).*(r_atoms[:,2].>0.0)
		r_atoms_temp = r_atoms[mirror_selection,:]
		r_atoms_temp[:,1].*=-1
		r_atoms = vcat(r_atoms,r_atoms_temp)
		r_atoms_temp[:,2].*=-1
		r_atoms = vcat(r_atoms,r_atoms_temp)
		r_atoms_temp[:,1].*=-1
		r_atoms = vcat(r_atoms,r_atoms_temp)
	end
	#
    return (r_atoms[:,:],length(r_atoms[:,1]),phase_array,r_range,collect(phase_range), lattice_array)
end
#
#Overall shift of the phase to range [0,2pi] rather than [-pi,pi]
function shift_phase(phi)
	if phi<=0.0
		return phi+2*pi
	else
		return phi
	end
end
#
#Functions to construct a lattice
function boundaries_choice(nAtoms)
    if nAtoms%2==0
        return [-nAtoms/2; nAtoms/2-1]
    else
        return [-(nAtoms-1)/2; (nAtoms-1)/2]
    end
end
#
#Core part of the ring creation. 
#It creates an array with given lattice constants and then selects only the atoms with radius within the ring
function metalens_creation_core(r_lens, lattice_constants, r_max,r_min)
    a_x=lattice_constants[1]*lambda0
    a_y=lattice_constants[2]*lambda0
    a_z=lattice_constants[3]*lambda0
    naX = Int(floor(2*r_lens/a_x))+1
    naY = Int(floor(2*r_lens/a_y))+1
	naX%2!=0 ? naX+=1 : nothing
	naY%2!=0 ? naY+=1 : nothing
    boundX=boundaries_choice(naX)
    boundY=boundaries_choice(naY)
    xOption=(naX+1)%2
    yOption=(naY+1)%2
    lattice_array_x=[(i + xOption/2)*a_x for i in boundX[1]:boundX[2] for j in boundY[1]:boundY[2]]
    lattice_array_y=[(j + yOption/2)*a_y for i in boundX[1]:boundX[2] for j in boundY[1]:boundY[2]]
	#
    lattice_norms=((xx,yy)->sqrt(xx^2+yy^2)).(lattice_array_x,lattice_array_y)
	#Selects the atoms inside the chosen ring
    selected_points=(lattice_norms.<r_max).*(lattice_norms.>=r_min).*(lattice_array_x.>=0.0).*(lattice_array_y.>=0.0)
    lattice_array_x=lattice_array_x[selected_points]
    lattice_array_y=lattice_array_y[selected_points]
    n_points_selected=length(lattice_array_x)
	#
	#Same as metalens_creation_buffer_z() below, but for the whole ring, rather than just the buffer zones
	if z_fixed_option=="YES"
		return (hcat(repeat(lattice_array_x,3),repeat(lattice_array_y,3),[z*a_z for z in [-1 ; 0 ; 1 ].+0.0 for i in 1:n_points_selected ]), maximum(lattice_array_x),maximum(lattice_array_y) )
	else
		phi_func(rr) = mod(k0*(focal_point-sqrt(focal_point^2+rr^2))+phase_shift,2*pi)
		z_func(rr)   = (2*pi-mod(phi_func(rr),pi) )/(6*pi)
		return ( hcat(repeat(lattice_array_x,3),repeat(lattice_array_y,3),[z*(z_func(sqrt(lattice_array_x[i]^2+lattice_array_y[i]^2))) for z in [-1 ; 0 ; 1 ].+0.0 for i in 1:n_points_selected ]) , maximum(lattice_array_x),maximum(lattice_array_y) )
	end
	#
	#
end
#
#Uppest positions of a series of columns (i.e. a collection of the higher point of each column)
#This also applies to rows, where now uppest means the ones placed most on the right side
function get_upper(pos_array_1,pos_array_2)
	index_sorted=sort(collect(1:length(pos_array_1)), by=(xx->pos_array_1[xx]))
	pos_array_x=pos_array_1[index_sorted]
	pos_array_y=pos_array_2[index_sorted]
	current_max=pos_array_y[1]
	current_x=pos_array_x[1]
	up_x=[current_x]
	up_y=[current_max]
	current_index=1
	for i in 2:length(pos_array_x)
		x_here = pos_array_x[i]
		y_here = pos_array_y[i]
		if x_here>current_x
			up_x=vcat(up_x,x_here)
			up_y=vcat(up_y,y_here)
			current_x=x_here
			current_max = y_here
			current_index+=1
		else
			if y_here>current_max
				up_x[current_index]=x_here
				up_y[current_index]=y_here
				current_max = y_here
			end
		end
	end
	return (up_x,up_y)
end
#
#Same as before but for lowest/most on the left positions
function get_lower(pos_array_1,pos_array_2)
	index_sorted=sort(collect(1:length(pos_array_1)), by=(xx->pos_array_1[xx]))
	pos_array_x=pos_array_1[index_sorted]
	pos_array_y=pos_array_2[index_sorted]
	current_min=pos_array_y[1]
	current_x=pos_array_x[1]
	down_x=[current_x]
	down_y=[current_min]
	current_index=1
	for i in 2:length(pos_array_x)
		x_here = pos_array_x[i]
		y_here = pos_array_y[i]
		if x_here>current_x
			down_x=vcat(down_x,x_here)
			down_y=vcat(down_y,y_here)
			current_x=x_here
			current_min = y_here
			current_index+=1
		else
			if y_here<current_min
				down_x[current_index]=x_here
				down_y[current_index]=y_here
				current_min = y_here
			end
		end
	end
	return (down_x,down_y)
end
#
#
#Core functions to create the buffer zones
function metalens_creation_buffer(up_nodes_1, up_nodes_2,low_nodes_1, low_nodes_2, a_2, index_up_down)
	index_low = 1
	index_low_max = length(low_nodes_2)
	r_atoms_buffer_x = Array{Float64}(undef,0)
	r_atoms_buffer_y = Array{Float64}(undef,0)
	for i in 1:length(up_nodes_1)
		if index_low<=index_low_max
			p1_x=up_nodes_1[i]
			p1_y=up_nodes_2[i]
			p2_x=low_nodes_1[index_low]
			p2_y=low_nodes_2[index_low]
			final_index_low = index_low
			for j in (index_low+1):index_low_max
				p2_x_test = low_nodes_1[j]
				p2_y_test = low_nodes_2[j]
				if abs(p1_x-p2_x_test)<abs(p1_x-p2_x)
					final_index_low = j
					p2_x = p2_x_test
					p2_y = p2_y_test
				end
			end
			if (p1_y>p2_y,p1_y<p2_y )[index_up_down]
				index_low = final_index_low
				(x_points_to_add, y_points_to_add) = line_connect(p1_x, p1_y,p2_x,p2_y, a_2)
				r_atoms_buffer_x = vcat(r_atoms_buffer_x,x_points_to_add)
				r_atoms_buffer_y = vcat(r_atoms_buffer_y,y_points_to_add)
				index_low+=1
			end
		end
	end
	return (r_atoms_buffer_x, r_atoms_buffer_y)
end
#
function line_connect(p1_x, p1_y,p2_x,p2_y, a_2)
	y_min = min(p1_y, p2_y)+a_2*0.9
	y_max = max(p1_y, p2_y)-a_2*0.9
	y_points = collect(y_min:a_2:y_max)
	x_points = @. (p1_x*y_points - p2_x*y_points + p2_x*p1_y - p1_x*p2_y)/(p1_y - p2_y)
	return (x_points,y_points)
end
#
#Function to create the z coordinate of the points inside the buffer zones
#If a_z_sharp is different from zero, then the code assign the fixed value a_z_sharp to the z-coordinate of all the points
#Otherwise it connects them in space by calculating xi_z(phi_lens(r)). This latter case applies only if z_fixed_buffer_option == "NO"
function metalens_creation_buffer_z(r_atoms_buffer_x, r_atoms_buffer_y, a_z_sharp)
	n_points_selected = length(r_atoms_buffer_x)
	if a_z_sharp>0.0
		return hcat(repeat(r_atoms_buffer_x,3),repeat(r_atoms_buffer_y,3),[z*a_z_sharp for z in [-1 ; 0 ; 1 ].+0.0 for i in 1:n_points_selected ])
	else
		phi_func(rr) = mod(k0*(focal_point-sqrt(focal_point^2+rr^2))+phase_shift,2*pi)
		z_func(rr)   = (2*pi-mod(phi_func(rr),pi) )/(6*pi)
		return hcat(repeat(r_atoms_buffer_x,3),repeat(r_atoms_buffer_y,3),[z*(z_func(sqrt(r_atoms_buffer_x[i]^2+r_atoms_buffer_y[i]^2))) for z in [-1 ; 0 ; 1 ].+0.0 for i in 1:n_points_selected ])
	end
end
