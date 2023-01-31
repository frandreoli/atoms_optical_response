### atoms_optical_response

# Introduction

The code performs numerical simulations of the steady-state, optical response of large-scale systems of quantum, two-level emitters (atoms, color centers, etc.).
Various options are provided on the density, geometry and physical properties of the light scatterers.

The original code has been developed as a part of a PhD project by Francesco Andreoli, under the supervision of Prof. Dr. Darrick Chang.

## Code specifications

The core code is written in Julia (tested on version Julia 1.8), while the output data are saved in the HDF5 format. 
A Mathematica notebook is provided to read the data and generate simple plots, with an intuitive, user-friendly interface. This latter is currently available only for the "Metalens" option (see below).
The simulation is specifically optimized for multi-threaded computation, and the code was extensively used on a multi-core (>32 threads) cluster with ~400GB of RAM to simulate systems of $\sim 5\times 10^5$ atoms.




## Physical background and motivation

In this project, we implement a code to simulate the collective behavior of a mesoscopic ensemble of quantum emitters. This task is crucial to predict new phenomena that can occur in actual experiments, where cooperative effects can drastically alter the optical response [[1](Dicke1954CoherenceProcesses)–[16](Asenjo-Garcia2017ExponentialArrays)]. 

A single, two-level emitter is characterized by a dipole matrux element $d$ associated with the transition $\omega_0$ from the ground to the excited state. A collection of $N$ interacting atomic dipoles $\mathbf{d}$  $_j$
illuminated by an external field $\bm{E}_{\rm in}$(r_j) can be described by the coupled-dipole formalism, where each dipole is driven by both the input field and the field re-scattered by the other dipoles. With good approximation, the dipoles align in the same direction of the field polarization (say x ̂), leading to (14,16)

$$ ... $$

where G ̿(r_j-r_k,ω_0 ) is the electromagnetic Green’s tensor in vacuum describing the scattering pattern of each dipole, while α_0 is the atomic polarizability. This problem is equivalent to inverting the algebraic system ∑_j▒〖M_jk d_j 〗=E_in (r_j ), characterized by the N×N complex, non-Hermitian matrix M_jk=δ_jk/(α_0 ϵ_0)-(1-δ_jk)μ_0 ω_0^2  x ̂⋅G ̿(r_j-r_k,ω_0 )⋅x ̂. Similarly, the total field can be reconstructed at N_probe probe points in space, as E_tot (r_k )=E_in (r_j )+∑_j▒〖L_jk  d_j 〗, where we define the 3N_probe×N matrix  L_jk=μ_0 ω_0^2  G ̿(r_k-r_j,ω_0 )⋅x ̂. Hereafter, we discuss the design of a code which we used to simulate up to N∼7×〖10〗^5 atoms, roughly two orders of magnitude higher than comparable works  (5,10–12,14,17–23).


cite Dicke [1](Dicke1954CoherenceProcesses)


$$SE = \frac{\sigma}{\sqrt{n}}$$
 



# Initializing the code

## Overall settings

## Specific settings



### Disordered geometries settings

### Array settings

### Atomic metalens settings

### Custom geometry settings


# Run a simulation

# Data outputs



# References 



<a id="Dicke1954CoherenceProcesses">[1]</a> 
Dicke RH. 
*Coherence in Spontaneous Radiation Processes.* 
[Physical Review 93(1):99–110](https://link.aps.org/doi/10.1103/PhysRev.93.99), (1954) 

<a id="Asenjo-Garcia2017ExponentialArrays">[16]</a> 
Asenjo-Garcia A, Moreno-Cardoner M, Albrecht A, Kimble HJ, Chang DE. 
*Exponential Improvement in Photon Storage Fidelities Using Subradiance and “Selective Radiance” in Atomic Arrays.*
[Phys Rev X 7(3):31024](https://link.aps.org/doi/10.1103/PhysRevX.7.031024) (2017)



