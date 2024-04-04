#############################################################################################################
################## COOPERATIVE PROPERTIES OF ARRAY ##########################################################
#############################################################################################################
#
#
#Initialization
using SpecialFunctions: erfi
#
#
#Functions for the regularization of the self-energy
function g0(theta)
    -(
    (erfi(theta/sqrt(2))-1.0im) / (exp(theta^2/2))
    - (-0.5+theta^2) / (sqrt(pi/2)*theta^3)
    ) /2
end
#In the following the g_x and g_y are intended in units of k0
function Lambda_func(g_x,g_y)
    sqrt(complex(1- g_x^2 - g_y^2 ))
end
#
function f1_func(theta,Lambda)
    ( -1.0im+erfi(theta*Lambda/sqrt(2)) )/Lambda
end
#
function I0_core(theta,Lambda,f1_val, g_x, g_y, d_vec)
    solution = zeros(Complex{Float64}, 3,3)
    solution[1,1] = (1-g_x^2)*f1_val
    solution[2,2] = (1-g_y^2)*f1_val
    solution[2,1] = solution[1,2] = -g_x*g_y*f1_val
    solution[3,3] = (1-Lambda^2)*f1_val + exp(theta^2*Lambda^2/2)*sqrt(2/pi)/theta
    return dot(d_vec,-solution*d_vec)
end
#
function I0_func_generic(a,b,theta,xi_x,xi_y,k_in_x,k_in_y, d_vec)
    k_tot_x = a/xi_x - k_in_x/(2*pi)
    k_tot_y = b/xi_y - k_in_y/(2*pi)
    Lambda=Lambda_func(k_tot_x, k_tot_y)
    f1_val = f1_func(theta,Lambda)
    return I0_core(theta,Lambda,f1_val, k_tot_x, k_tot_y, d_vec)
end
#
#Function to sum up to convergence in a smart way
function smart_sum(func,args...)
    #Settings for convergence
    max_steps_conv=2000
    convergence_countdown_max=3
    convergence_threshold=10.0^(-8)
    #Initialization
    solution=old_solution=0.0
    convergence_countdown=0
    #Smart sum
    for i in 1:max_steps_conv
        to_add=func(i,args...)
        isnan(to_add) ? nothing : solution+=to_add
        if abs(solution)!=0.0
            abs(solution-old_solution)/abs(solution) < convergence_threshold ? convergence_countdown+=1 : convergence_countdown=0
        else
            solution-old_solution==0 ? convergence_countdown+=1 : convergence_countdown=0
        end
        convergence_countdown>=convergence_countdown_max ?  break  : nothing
        old_solution=solution
        i==max_steps_conv && length(args)==0 ? error("Reached convergence threshold in smart_sum()") : nothing
    end
    return solution
end
#
#Function to evaluate the self energy
function coop_value_theta(theta,xi_x,xi_y,k_in_x,k_in_y, d_vec)
    I0_func(a,b) = I0_func_generic(a,b,theta,xi_x,xi_y,k_in_x,k_in_y, d_vec)
    solution=0
    solution+=I0_func(0,0)
    #
    solution+=smart_sum((a,b) -> I0_func(a,b) , 0)
    solution+=smart_sum((a,b) -> I0_func(-a,b) , 0)
    #
    solution+=smart_sum((b,a) -> I0_func(a,b) , 0)
    solution+=smart_sum((b,a) -> I0_func(a,-b) , 0)
    #
    solution+=smart_sum( b_main ->  smart_sum( (a,b) -> I0_func(a,b)  , b_main)  )
    solution+=smart_sum( b_main ->  smart_sum( (a,b) -> I0_func(-a,b)  , b_main)  )
    solution+=smart_sum( b_main ->  smart_sum( (a,b) -> I0_func(a,-b)  , b_main)  )
    solution+=smart_sum( b_main ->  smart_sum( (a,b) -> I0_func(-a,-b)  , b_main)  )
    #
    (3/(8*pi*xi_x*xi_y))*exp(-theta^2/2)*solution-g0(theta) + 0.5im
end
#
#Accurate convergence function (it returns omega_coop and Gamma_coop)
function coop_values_function(xi_x,xi_y,k_in_x,k_in_y, d_vec)
    error_threshold=10.0^(-2)
    theta_now=(xi_x+xi_y)/(10)#0.04
    theta_min=0.0001
    decrement_ratio=1.5#1.2
    omega_coop_old=0
    coop_values_complex = coop_value_theta(theta_now,xi_x,xi_y,k_in_x,k_in_y, d_vec) 
    omega_coop_new= -real(coop_values_complex)
    #
    while theta_now>theta_min
        theta_now /= decrement_ratio
        coop_values_complex = coop_value_theta(theta_now,xi_x,xi_y,k_in_x,k_in_y, d_vec)
        omega_coop_new = -real(coop_values_complex)
        abs((omega_coop_new-omega_coop_old))/abs(omega_coop_new) < error_threshold ? break : nothing
        omega_coop_old = omega_coop_new
        theta_now<=theta_min ? error("theta not convergent") : nothing
    end
    (omega_coop_new, 2*imag(coop_values_complex))
end
#
#
#
#
#
#
#
#############################################################################################################
################## ARRAY CONSTRUCTION #######################################################################
#############################################################################################################
#
#
#Functions to construct a lattice
function array_boundaries(nAtoms)
    if nAtoms%2==0
        return [-nAtoms/2; nAtoms/2-1]
    else
        return [-(nAtoms-1)/2; (nAtoms-1)/2]
    end
end
#
#Function to create the atomic array
function arrays_creation(m_planes, xi_x,xi_y,xi_z,array_size_x,array_size_y)
    #
    #Lattice properties
    l_system_x = abs(sum(array_size_x.*[-1;1]))
    l_system_y = abs(sum(array_size_y.*[-1;1]))
	naX = Int(floor(l_system_x/xi_x))+1
    naY = Int(floor(l_system_y/xi_x))+1
    #
    #Shift to complain to mirror symmetry
	boundX=array_boundaries(naX)
	boundY=array_boundaries(naY)
	boundZ=array_boundaries(m_planes)
	xOption=(naX+1)%2
	yOption=(naY+1)%2
	zOption=(m_planes+1)%2
    #
    #Lattice creation
	lattice_array=[[(i + xOption/2)*xi_x; (j + yOption/2)*xi_y ; (k + zOption/2)*xi_z] for j in boundY[1]:boundY[2] for i in boundX[1]:boundX[2] for k in boundZ[1]:boundZ[2]]
	#
    if mirror_symmetry_option_const == "YES"
        selected_positions = (x->(x[1]>=0.0 && x[2]>=0.0)).(lattice_array)
        lattice_array = lattice_array[selected_positions]
    end
    #
    #Saving the positions
    r_atoms = Array{Float64}(undef, length(lattice_array), 3)
	n_atoms = length(lattice_array)
    central_pos_x = sum(array_size_x)/2.0
    central_pos_y = sum(array_size_y)/2.0
	for i in 1:n_atoms
		r_atoms[i,:]= lattice_array[i] .+ [central_pos_x; central_pos_y; 0.0]
	end
    #
	return (r_atoms, n_atoms)
end