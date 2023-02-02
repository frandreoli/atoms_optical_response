  
#############################################################################################################
################## INITIALIZATION ###########################################################################
#############################################################################################################
#
#
using SpecialFunctions: erfi
#
#
#
#
#
#
#
#############################################################################################################
################## COOPERATIVE PROPERTIES ###################################################################
#############################################################################################################
#
#=
#Functions for the regularization of the self-energy
function g0(theta)
    -(
    (erfi(theta/sqrt(2))-1.0im) / (exp(theta^2/2))
    - (-0.5+theta^2) / (sqrt(pi/2)*theta^3)
    ) /2
end
#
function I0_lambda(a,b,xi_x,xi_y)
    sqrt(complex(1- (a/xi_x)^2 - (b/xi_y)^2 ))
end
#
function I0_func(a,b,theta,xi_x,xi_y)
    lambda=I0_lambda(a,b,xi_x,xi_y)
    -pi*(-1.0im+erfi(theta*lambda/sqrt(2)))/lambda
end
#
#Function to sum up to convergence in a smart way
function smart_sum(func,args...)
    #Settings for convergence
    max_steps_conv=2000
    convergence_countdown_max=3
    convergence_threshold=10.0^(-20)
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
function delta_coop_theta(theta,xi_x,xi_y)
    solution=0
    solution+=I0_func(0,0,theta,xi_x,xi_y)
    solution+=4*smart_sum( b_main ->  smart_sum( (a,b) -> (1-(a/xi_x)^2)*I0_func(a,b,theta,xi_x,xi_y)  , b_main)  )
    solution+=2*smart_sum((a,b) -> (1-(a/xi_x)^2)*I0_func(a,b,theta,xi_x,xi_y) , 0)
    solution+=2*smart_sum((b,a) -> (1-(a/xi_x)^2)*I0_func(a,b,theta,xi_x,xi_y) , 0)
    (3/(8*pi^2*xi_x*xi_y))*exp(-theta^2/2)*solution-g0(theta)
end
#
#Accurate convergence function
function omega_coop_function(xi_x,xi_y)
    error_threshold=10.0^(-2)
    theta_now=0.04
    theta_min=0.005
    decrement_ratio=1.2
    delta_coop_old=0
    delta_coop_new=real(delta_coop_theta(theta_now,xi_x,xi_y))
    #
    while theta_now>theta_min
        theta_now /= decrement_ratio
        delta_coop_new = real(delta_coop_theta(theta_now,xi_x,xi_y))
        abs((delta_coop_new-delta_coop_old))/abs(delta_coop_new) < error_threshold ? break : nothing
        delta_coop_old = delta_coop_new
        theta_now<=theta_min ? error("theta not convergent") : nothing
    end
    -delta_coop_new
end
#
#Function to evaluate the Gamma_coop, when ξ<1
function gamma_coop_function(xi_x,xi_y)
    3/(4 *pi* xi_x*xi_y)
end
=#
#