#############################################################################################################
################## COOPERATIVE PROPERTIES ARRAY##############################################################
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
#Accurate convergence function
function coop_values_function(xi_x,xi_y,k_in_x,k_in_y, d_vec)
    error_threshold=10.0^(-2)
    theta_now=0.04
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

xi_x=0.3
xi_y=xi_x*2
k_x=1/sqrt(3)
k_y=2/sqrt(3)

println(coop_values_function(xi_x,xi_y,k_x,k_y, [1;0;0]))