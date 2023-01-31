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

In this project, we implement a code to simulate the collective behavior of a mesoscopic ensemble of quantum emitters. This task is crucial to predict new phenomena that can occur in actual experiments, where cooperative effects can drastically alter the optical response [[1](Dicke1954)–[16](Asenjo-Garcia2017)]. 

We consider two-level emitters with dipole matrix elements $d\_0$, associated to the transition $\omega_0$ from the ground to the excited state. This latter has a spontaneous emission rate $\Gamma\_0$, which defines the elastic scattering rate of the atoms. A supplementary rate $\Gamma'$ (in units of $\Gamma\_0$ accounts for the possibility of inelastic scattering, which might be due to the system environment (such as the excitation of phononic modes). 
A collection of $N$ emitters is characterized by their dipoles $\mathbf{d}\_j$. When illuminated by an external field $\mathbf{E}\_{\text{in}}(\mathbf{r}\_j)$, the optical response can be described by the coupled-dipole formalism, where each dipole is driven by both the input field and the field re-scattered by the other dipoles. This is valid as long as the external drive has a low intensity (with Rabi frequecy $\Omega_0\ll \Gamma_0$), so that non-linear effects are negligible. In these simulations, we assume that the dipoles have a well defined direction $\mathbf{d}\_j = d\_j \hat{\mathbf{d}}\_0$. When dealing with atoms in free space, this is a valid approximation provided that the input field is $\hat{\mathbf{x}}$-polarized and that it retains more intense that the scattere field at the atomic positions. 
On the contrary, a constrain on the dipole orientations is imposed on solid-state emitters by the geometry of the environment (such as the surrounding diamond lattice for color centers in diamond). Finally, we assume that the emitters are in a bulk, dielectric material of (possitive, real) index $n\_{\text{bulk}}$, meaning that $k_0=2\pi/\lambda_0 = n\_{\text{bulk}}\omega_0/c$ denotes the resonant wavevector. Under these conditions, the steady-state response is described by [14,[16](Asenjo-Garcia2017)]

 $$ d\_j = \alpha\_0\epsilon\_0 \left(\hat{\mathbf{d}}\_0\cdot \mathbf{E}\_{\text{in}}(\mathbf{r}\_j) + \mu\_0 \omega\_0^2 \displaystyle \sum\_{k\neq j}^N \hat{\mathbf{d}}\_0 \cdot \mathbf{G}(\mathbf{r}\_j-\mathbf{r}\_k,\omega\_0)\cdot \hat{\mathbf{d}}\_0 d_k \right), $$

where $\alpha\_0=-3\pi \epsilon\_0 /\[(\Delta + i(1+\Gamma')/2)k\_0^3\]$ is the atomic polarizability, with $\Delta=(\omega - \omega_0)/\Gamma_0$ representing the dimensionless detuning between the frequency $\omega$ of the incident light and the atomic resonance $\omega_0$, in units of the natural linewidth $\Gamma_0$. At the same time, $\mathbf{G}(\mathbf{r}\_j-\mathbf{r}\_k,\omega\_0)$ is the electromagnetic Green’s tensor. This problem is equivalent to inverting the algebraic system $\sum\_j \mathcal{M}\_{jk} d\_j = \hat{\mathbf{d}}\_0\cdot \mathbf{E}\_{\text{in}}(\mathbf{r}\_j)$, characterized by the $N\times N$ complex, non-Hermitian matrix $$. Similarly, the total field can be reconstructed at N_probe probe points in space, as $$, where we define the $3N\_{\text{probe}}$ matrix $$. 


With this code, it is possible to calculate the total field at different probe positions, given many geometric choices of the atomic positions (either ordered or with random distributions) and a Gaussian input beam of waist $w\_0$ and focal point $\mathbf{r}=0$, with custom direction and polarization. Moreover, the code calculates the transmission and reflection by projecting the scattered field onto the same Gaussian mode as the input [,[14](Andreoli2021)]. On top of that, the user has the possibility to calculate as well the projection and reflection onto a different Gaussian mode, with custom waist $w\_{\text{target}}$ and custom focal point. At the same time, the computation can be performed for several values of the detuning $\Delta$, and given any choice of the dipole orientation $\hat{\mathbf{d}}\_0$. Finally, when simulating lattice geometries, the user can test scenarios where a fixed fractions of random defects (i.e. missing atoms at lattice positions) is present. For any random setting, the computation can be repeated multiple times, each sampling different configurations.




This code has been used to simulate systems of up to $\sim 5\times 10^5$ emitters, roughly two orders of magnitude larger than comparable works [5,10–12,14,17–23].


## Key code 

# Code guide
## Initializing the code

### Overall settings

### Specific settings



#### Disordered geometries settings

#### Array settings

#### Atomic metalens settings

#### Custom geometry settings


## Run the simulation

## Data outputs



# References 



<a id="Dicke1954">[1]</a> 
Dicke RH. 
*Coherence in Spontaneous Radiation Processes.* 
[Physical Review 93(1):99–110](https://link.aps.org/doi/10.1103/PhysRev.93.99), (1954) 

<a id="Gross1982">[2]</a> 
Gross M, Haroche S. *Superradiance: An essay on the theory of collective spontaneous emission.* 
Physics Reports. 1982 Dec 1;93(5):301–96. 

<a id="Asenjo-Garcia2017">[16]</a> 
Asenjo-Garcia A, Moreno-Cardoner M, Albrecht A, Kimble HJ, Chang DE. 
*Exponential Improvement in Photon Storage Fidelities Using Subradiance and “Selective Radiance” in Atomic Arrays.*
[Phys Rev X 7(3):31024](https://link.aps.org/doi/10.1103/PhysRevX.7.031024) (2017)



