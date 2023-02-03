#############################################################################################################
################## FUNCTIONS TO RANDOMLY SAMPLE POINTS ######################################################
#############################################################################################################
#
#
#Function to uniformly sample in a cylinder with main axis oriented along z
function uniform_sampling_cylinder!(r_atoms,n_atoms, rLength, zLength)
  length(r_atoms[:,1])!=n_atoms ? error("Incompatible number of atoms in positions sampling.") : nothing
  length(r_atoms[1,:])!=3       ? error("Wrong vector space dimension in positions sampling.") : nothing
  phasesRandom=2.0*pi*rand(n_atoms)
  radiusRandom=sqrt.(rand(n_atoms))
  #
  #See: http://mathworld.wolfram.com/DiskPointPicking.html
  xPoints=rLength*radiusRandom.*cos.(phasesRandom)
  yPoints=rLength*radiusRandom.*sin.(phasesRandom)
  zPoints=rand(n_atoms).*(zLength).-(zLength/2)
  for i in 1:n_atoms
      r_atoms[i,:]=[xPoints[i];yPoints[i];zPoints[i]]
  end
end
#
#Function to uniformly sample in a 3D ball/sphere
function uniform_sampling_sphere!(r_atoms,n_atoms, r_sphere, type_sampling="BALL",hemisphere="FULL_SPHERE")
  length(r_atoms[:,1])!=n_atoms ? error("Incompatible number of atoms in positions sampling.") : nothing
  length(r_atoms[1,:])!=3       ? error("Wrong vector space dimension in positions sampling.") : nothing
  #
  #See: 
  #http://mathworld.wolfram.com/SpherePointPicking.html 
  #http://mathworld.wolfram.com/DiskPointPicking.html 
  #http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
  #
  thetaRandom=(2.0*pi).*rand(n_atoms)
  uRandom=2.0.*rand(n_atoms).-1.0
  #
  #Deciding if spanning the full solid angle or half (in the forward/backward direction)
  if hemisphere!="FULL_SPHERE" && hemisphere!="FORWARD_HEMISPHERE" && hemisphere!="BACKWARD_HEMISPHERE"
    error("Invalid definition of hemisphere in uniform_sampling_sphere!()")
  end
  hemisphere=="FORWARD_HEMISPHERE"  ? uRandom=abs.(uRandom)  : nothing
  hemisphere=="BACKWARD_HEMISPHERE" ? uRandom=-abs.(uRandom) : nothing
  #
  #Deciding if sampling a ball or a sphere
  if type_sampling!="BALL" && type_sampling!="SPHERE" 
    error("Invalid definition of type_sampling in uniform_sampling_sphere!()")
  end
  if type_sampling=="BALL"
    rRandom=(rand(n_atoms)).^(1/3)
  else
    rRandom=[1.0 for i in 1:n_atoms]
  end
  #
  #Calculating the points
  xPoints=r_sphere.*rRandom.*sqrt.(1.0.-(uRandom).^2).*cos.(thetaRandom)
  yPoints=r_sphere.*rRandom.*sqrt.(1.0.-(uRandom).^2).*sin.(thetaRandom)
  zPoints=r_sphere.*rRandom.*uRandom
  #
  #Storing them into an array
  for i in 1:n_atoms
      r_atoms[i,:]=[xPoints[i];yPoints[i];zPoints[i]]
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
################## FUNCTIONS TO CALCULATE THE PROBE POINTS ##################################################
#############################################################################################################
#
#
#FUNCTIONS TO CALCULATE THE POSITIONS WHERE THE TOTAL FIELD MUST BE COMPUTED:
#
#Probe points computed on the plane perpendicular to v3_vec
function f_probe_PLANE(v3_vec, v3_value, points_v1, points_v2, range_v1, range_v2)
	rot_matrix = inv(rot_matrix_function(v3_vec))
	v1Range=range(range_v1[1],stop=range_v1[2],length=points_v1)
	v2Range=range(range_v2[1],stop=range_v2[2],length=points_v2)
	r_probe=Array{Float64}(undef,points_v1,points_v2,3)
	for i_v1 in 1:points_v1
		for i_v2 in 1:points_v2
      if v3_vec!=[1.0;0.0;0.0]
			  r_probe[i_v1,i_v2,:]=rot_matrix*[v1Range[i_v1] ; v2Range[i_v2] ; v3_value]
      else
        r_probe[i_v1,i_v2,:]=rot_matrix*[v2Range[i_v2] ; v1Range[i_v1] ; v3_value]
      end
		end
	end
	return reshape(r_probe,points_v1*points_v2,3)
end
#
#
#Probe points computed on a sphere/hemisphere surrounding the atoms
function f_probe_SPHERE(radius,points,probeSphere_option)
  r_probe=Array{Float64}(undef,points,3)
  r_probe[:,:].=0.0
  uniform_sampling_sphere!(r_probe,points, radius, "SPHERE",probeSphere_option)
  return r_probe
end
#
#
#
#
#
#
#
#############################################################################################################
################## LINEAR ALGEBRA FUNCTIONS #################################################################
#############################################################################################################
#
#The following perform linear algebra operations on complex vectors
#They are used only if the field/atomic polarizations are complex
function conj_dot(arr_a,arr_b)
	dot(arr_a,arr_b)
end
#
function conj_norm(arr_a)
	sqrt(conj_dot(arr_a,arr_a))
end
#
function conj_normalize(arr_a)
	arr_a./conj_norm(arr_a)
end
#
#Function to rotate the vector vec_start to [0;0;1] in 3D
#inv(rot_matrix_function(vec_start)) rotates the vector [0;0;1] to vec_start
#It is the analogous of RotationMatrix[{{x, y, z}, {0, 0, 1}}] in Wolfram Mathematica
function rot_matrix_function(vec_start) 
  vec_start = conj_normalize(vec_start)
  if norm(vec_start.-[0.0 ; 0.0 ; 1.0])<=ZERO_THRESHOLD
    matrix = I
  else
    (x,y,z) = Tuple(vec_start)
    matrix = [
    ((y^2 + z + z^2)/(1 + z))   (-((x*y)/(1 + z)))           -x ;
    (-((x*y)/(1 + z)))           ((1 - y^2 + z)/(1 + z))     -y ;
    ((x - x*z^2)/(x^2 + y^2))    ((y - y*z^2)/(x^2 + y^2))   ((z - z^3)/(x^2 + y^2))
    ]
  end
  return matrix
end
#
#
#
#
#
#
#
#############################################################################################################
################## PHYSICAL FUNCTIONS DEFINITIONS ###########################################################
#############################################################################################################
#
#
#Functions to calculate input Gaussian beam at position r_vec. 
#laser_direction is the propagation direction (the default is z) and the beam has waists w0_x and w0_y in the x and y directions.
#In this code we always set w0_x=w0_y
#
function zr(w0, lambda)
  pi*w0^2/lambda
end
#Beam width w(z)
function w(z, w0, lambda)
  w0*sqrt(1 + (z/zr(w0, lambda))^2)
end
#Inverse of Radius of curvature of beam R(z)
function invr(z, w0, lambda)
  z/(z^2 + zr(w0, lambda)^2)
end
#
function gaussian_beam(x_vec_start,y_vec_start, z_vec_start, w0_x, w0_y, k, laser_direction_x, laser_direction_y, laser_direction_z)
  laser_direction_vec = [laser_direction_x ; laser_direction_y ; laser_direction_z]
  #Accounting for Gaussian beams travelling in different directions
  if norm(laser_direction_vec.-[0.0;0.0;1.0])>ZERO_THRESHOLD
    #
    #Rotation matrix that transforms laser_direction_vec -> [0 ; 0 ; 1]
    rot_matrix=rot_matrix_function(laser_direction_vec) 
    #
    (x_vec, y_vec, z_vec) = Tuple(rot_matrix*[x_vec_start ; y_vec_start ; z_vec_start])
  else
    (x_vec, y_vec, z_vec) = (x_vec_start , y_vec_start , z_vec_start)
  end
  #
  lambda = 2*pi/k
  w_x = w(z_vec, w0_x, lambda)
  w_y = w(z_vec, w0_y, lambda)
  norm_factor = sqrt(w0_x*w0_y/(w_x*w_y))
  real_arg    = -(x_vec^2/w_x^2 + y_vec^2/w_y^2)
  imag_arg    = k*z_vec + k*(x_vec^2*invr(z_vec, w0_x, lambda) + y_vec^2*invr(z_vec, w0_y, lambda))/2 - atan(z_vec/zr(w0_x, lambda))/2 - atan(z_vec/zr(w0_y, lambda))/2
  return norm_factor*exp(real_arg + im*imag_arg)
end
#
#
#
#
#
#
#
#############################################################################################################
################## SAVING FUNCTIONS #########################################################################
#############################################################################################################
#
#
#FUNCTIONS TO SAVE DATA FILES IN HDF5 FORMAT:
#
#Function to save any variable. 
#"cw" -> Do not overwrite existing data
#"w"  -> Overwrite existing data 
function h5write_basic(file_name,data, name_variable="" ; open_option="cw")
  file_h5=h5open(file_name*".h5", open_option)
  haskey(file_h5, name_variable) ? delete_object(file_h5, name_variable) : nothing
  file_h5[name_variable]=data
  close(file_h5)
end
#
#Function to save an array of complex number.
#In the HDF5 file the variable "_re" ("_im") contain an array with the the real (imaginary) parts
function h5write_complex(file_name,data, name_variable="" ; open_option="cw")
  length(name_variable)>0 ? add_name=name_variable*"_" : add_name = ""
  file_h5=h5open(file_name*".h5", open_option)
  haskey(file_h5, name_variable*"re") ? delete_object(file_h5, name_variable*"re") : nothing
  haskey(file_h5, name_variable*"im") ? delete_object(file_h5, name_variable*"im") : nothing
  file_h5[name_variable*"re"]=real.(data)
  file_h5[name_variable*"im"]=imag.(data)
  close(file_h5)
end
#
#Function to save multiples new variables into a file. 
#If the file already exists it preserves the data, otherwise it creates the file
#If the variable already exist in the file, it overwrites it
#data_array is an array of tuples whose first element is the name of the variable
#while the second is the variable to save. 
#The option open_option="cw" only overwrites the chosen variable in the .h5 file,
#while setting open_option="w" first deletes all elements stored in the file.
function h5write_multiple(file_name,data_array... ; open_option="cw")
  file_h5=h5open(file_name*".h5", open_option)
  for index in 1:length(data_array)
    name_variable = data_array[index][1]
    variable_data = data_array[index][2]
    haskey(file_h5, name_variable) ? delete_object(file_h5, name_variable) : nothing
    file_h5[name_variable]=variable_data
  end
  close(file_h5)
end
#
#
#Functions to add new elements to an array already stored
function h5write_append(file_name,data, name_variable="")
  open_option="cw"
  file_h5=h5open(file_name*".h5", open_option)
  if haskey(file_h5, name_variable)
    old_data = read(file_h5[name_variable])
    delete_object(file_h5, name_variable)
    try
      file_h5[name_variable]=vcat(old_data, data)
    catch error_writing
      close(file_h5)
      error(error_writing)
    end
  else
    try
      file_h5[name_variable]=data
    catch error_writing
      close(file_h5)
      error(error_writing)
    end
  end
  close(file_h5)
end
#
function h5write_complex_append(file_name,data, name_variable="")
  length(name_variable)>0 ? add_name="_" : add_name = ""
  for domain_f in ((real,"re") , (imag, "im"))
    h5write_append(file_name,(domain_f[1]).(data), name_variable*add_name*domain_f[2]) 
  end
end
#
#
#
#
#
#
#
#
#
#############################################################################################################
################## OTHER FUNCTIONS ##########################################################################
#############################################################################################################
#
#
#Function to estimate the RAM used by the code
function RAM_estimation(TN,tot_probe_points, n_atoms, n_detunings)
  bit_size = 1.0/8
  TN===Float32 ? sizeTN = 32*bit_size : sizeTN = 64*bit_size
  sizeTN_complex    = 2.0*sizeTN 
  size_matrix_G     = sizeTN_complex*n_atoms^2
  size_matrix_probe = sizeTN_complex*n_atoms*tot_probe_points
  size_coeff_array  = n_atoms*n_detunings
  (max(size_matrix_probe, size_matrix_G)+size_matrix_G+size_coeff_array)/(1024^3)
end
#
#Function to clear all methods of a function.
#The argument must be passed as func=:function_name
function clear_function_all(func)
  for method_to_clear in methods(func)
    Base.delete_method(method_to_clear)
  end
end
#
#Adding a dummy dimension to a matrix
function add_dimension(matrix)
  matrix_temp = Array{Float64}(undef, 1, length(matrix[:,1]), length(matrix[1,:]) )
	matrix_temp[1,:,:] = matrix[:,:]
  matrix_temp
end
#
#Cuts the output number to print. 
#It only works for numbers whose integer part has less than target_digits digits
function dig_cut(x_number, target_digits=15)
  string_x_number=string(x_number)
  x_digits = length(string_x_number)
  if x_digits>target_digits
    string_x_number=string_x_number[1:target_digits]
  elseif x_digits<target_digits
    if !occursin(".",string_x_number) 
      string_x_number*="."
      x_digits+=1
    end
    string_x_number*="0"^(target_digits-x_digits)
  end
  string_x_number
end
#
#
#
#
#
#
#
#############################################################################################################
################## PUNCH HOLES ##############################################################################
#############################################################################################################
#
#
#Function to punch defects in an atomic array/metalens
function defect_punching(r_atoms,n_atoms, defects_fraction)
  n_atoms_new = Int(round(n_atoms*(1-defects_fraction)))
  r_atoms_new = r_atoms[sort((shuffle(1:length(r_atoms[:,1])))[1:n_atoms_new]),:]
  (r_atoms_new, n_atoms_new)
end
#
#
#
#
#
#
#
#############################################################################################################
################## DISORDERED POSITIONS  ####################################################################
#############################################################################################################
#
#
#Function to define the atomic positions by uniformly sampling inside a sphere
function dis_sphere_creation(r_sphere,atomic_density)
  sphere_volume = (4/3)*pi*r_sphere^3
  n_atoms = Int(round(sphere_volume*atomic_density))
  r_atoms = Array{Float64}(undef, n_atoms, 3)
  uniform_sampling_sphere!(r_atoms,n_atoms, r_sphere, "BALL","FULL_SPHERE")
  return (r_atoms,n_atoms)
end
#
#
#Function to define the atomic positions by uniformly sampling inside a cylinder with main axis in the z direction
function dis_cyl_creation(r_disk,z_length,atomic_density)
  cylinder_volume = pi*r_disk^2*z_length
  n_atoms = Int(round(cylinder_volume*atomic_density))
  r_atoms = Array{Float64}(undef, n_atoms, 3)
  uniform_sampling_cylinder!(r_atoms,n_atoms, rLength, zLength)
  return (r_atoms,n_atoms)
end
#
#
#Function to define the atomic positions by uniformly sampling inside a cylinder with main axis in the z direction
function dis_cuboid_creation(x_dim,y_dim,z_dim,atomic_density)
  cuboid_volume = x_dim*y_dim*z_dim
  n_atoms = Int(round(cuboid_volume*atomic_density))
  r_atoms = Array{Float64}(undef, n_atoms, 3)
  for ii in 1:3
    r_atoms[:,ii].=(rand(Float64,n_atoms).*2.0 .- 1.0)*([x_dim ; y_dim ; z_dim][ii]/2)
  end
  return (r_atoms,n_atoms)
end
