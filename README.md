# 1) Introduction

The code performs numerical simulations of the steady-state, optical response of large-scale systems of quantum, two-level emitters (atoms, color centers, etc.).
Various options are provided on the density, geometry and physical properties of the light scatterers.

The original code has been developed as a part of a PhD project by Francesco Andreoli, under the supervision of Prof. Dr. Darrick Chang.

## 1.1) Code specifications

The core code is written in Julia (initially written for Julia 1.6 and well tested on the current version Julia 1.8), while the output data are saved in the HDF5 format. 
A Mathematica notebook is provided to read the data and generate simple plots, with an intuitive, user-friendly interface. This latter is currently available only for the "Metalens" option (see below).
The simulation is specifically optimized for multi-threaded computation, and the code was extensively used on a multi-core (>32 threads) cluster with ~400GB of RAM to simulate systems of up to $\sim 7\times 10^5$ atoms.




## 1.2) Physical background and motivation

In this project, we implement a code to simulate the collective behavior of a mesoscopic ensemble of quantum emitters. This task is crucial to predict new phenomena that can occur in actual experiments, where cooperative effects can drastically alter the optical response [[1](Dicke1954)–[16](Asenjo-Garcia2017)]. 

We consider two-level emitters with dipole matrix elements $d\_0$, associated to the transition $\omega_0$ from the ground to the excited state. This latter has a spontaneous emission rate $\Gamma\_0$, which defines the elastic scattering rate of the atoms. A supplementary rate $\Gamma'$ (in units of $\Gamma\_0$ accounts for the possibility of inelastic scattering, which might be due to the system environment (such as the excitation of phononic modes). 
A collection of $N$ emitters is characterized by their dipoles $\mathbf{d}\_j$. When illuminated by an external field $\mathbf{E}\_{\text{in}}(\mathbf{r}\_j)$, the optical response can be described by the coupled-dipole formalism, where each dipole is driven by both the input field and the field re-scattered by the other dipoles. This is valid as long as the external drive has a low intensity (with Rabi frequecy $\Omega_0\ll \Gamma_0$), so that non-linear effects are negligible. In these simulations, we assume that the dipoles have a well defined direction $\mathbf{d}\_j = d\_j \hat{\mathbf{d}}\_0$. When dealing with atoms in free space, this is a valid approximation provided that the input field is $\hat{\mathbf{x}}$-polarized and that it retains more intense that the scattere field at the atomic positions. 
On the contrary, a constrain on the dipole orientations is imposed on solid-state emitters by the geometry of the environment (such as the surrounding diamond lattice for color centers in diamond). Finally, we assume that the emitters are in a bulk, dielectric material of (possitive, real) index $n\_{\text{bulk}}$, meaning that $k_0=2\pi/\lambda_0 = n\_{\text{bulk}}\omega_0/c$ denotes the resonant wavevector. Under these conditions, the steady-state response is described by [14,[16](Asenjo-Garcia2017)]

 $$ d\_j = \alpha\_0\epsilon\_0 \left(\hat{\mathbf{d}}\_0\cdot \mathbf{E}\_{\text{in}}(\mathbf{r}\_j) + \mu\_0 \omega\_0^2 \displaystyle \sum\_{k\neq j}^N \hat{\mathbf{d}}\_0 \cdot \mathbf{G}(\mathbf{r}\_j-\mathbf{r}\_k,\omega\_0)\cdot \hat{\mathbf{d}}\_0 d_k \right), $$

where $\alpha\_0=-3\pi \epsilon\_0 /\[(\Delta + i(1+\Gamma')/2)k\_0^3\]$ is the atomic polarizability, with $\Delta=(\omega - \omega_0)/\Gamma_0$ representing the dimensionless detuning between the frequency $\omega$ of the incident light and the atomic resonance $\omega_0$, in units of the natural linewidth $\Gamma_0$. At the same time, $\mathbf{G}(\mathbf{r}\_j-\mathbf{r}\_k,\omega\_0)$ is the electromagnetic Green’s tensor. This problem is equivalent to inverting the algebraic system $\sum\_j \mathcal{M}\_{jk} d\_j = \hat{\mathbf{d}}\_0\cdot \mathbf{E}\_{\text{in}}(\mathbf{r}\_j)$, characterized by the $N\times N$ complex, non-Hermitian matrix $\mathcal{M}\_{jk}=\delta\_{jk}/(\alpha\_0\epsilon\_0)-(1-\delta\_{jk})\mu\_0\omega\_0^2 \hat{\mathbf{d}}\_0 \cdot \mathbf{G}(\mathbf{r}\_j-\mathbf{r}\_k,\omega\_0)\cdot \hat{\mathbf{d}}\_0$. 
Similarly, the total field can be reconstructed at $N\_{\text{probe}}$ probe points $\textbf{r}\_{k}$ in space, as $\mathbf{E}(\mathbf{r}\_{k})=\mathbf{E}\_{\text{in}}(\mathbf{r}\_{k})+\sum\_j\mathcal{L}\_{jk}d\_j$, where we define the $3N\_{\text{probe}}\times N$ matrix $\mathcal{L}\_{jk}=\mu\_0\omega\_0^2\mathbf{G}(\mathbf{r}\_j-\mathbf{r}\_k,\omega\_0)\cdot  \hat{\mathbf{d}}\_0 $. 


This code allows to calculate the total field at different probe positions, given many geometric choices of the atomic positions (either ordered or with random distributions). The input field consists of a Gaussian input beam of waist $w\_0$ and focal point $\mathbf{r}=0$, with custom direction and polarization. Moreover, the code calculates the transmission and reflection by projecting the scattered field onto the same Gaussian mode as the input [,[14](Andreoli2021)]. On top of that, the user has the possibility to calculate the projection and reflection onto a different Gaussian mode as well, with custom waist $w\_{\text{target}}$ and custom focal point (but same direction and polarization as the input). At the same time, the computation can be performed for several values of the detuning $\Delta$, and given any (fixed) choice of the dipole orientation $\hat{\mathbf{d}}\_0$. Finally, when simulating lattice geometries, the user can test scenarios where a fixed fractions of defects (i.e. missing atoms at random lattice positions) is present. For any random setting that is enabled, the computation can be repeated multiple times, each sampling different configurations. 

The number of emitters that we were able to simulate with this code is roughly two orders of magnitude larger than comparable works [5,10–12,14,17–23].






## 1.3) Computational insights

The core of our simulations is the inversion of $\mathcal{M}\_{jk}$, implemented with the backslash operator "\\" (see ["right division operator"](https://docs.julialang.org/en/v1/base/math/)). This performs adaptive algorithms based on the structure of random matrix, so that, in the worst-case scenario, the number of elementary operations scales as $\sim N^3$. This process is speeded up through the openBLAS library for linear algebra, which can evaluate the problem in parallel over up to 32 cores. Asymptotically, this task dominates the time complexity of the simulation, overcoming the other sources of time consumption, such as the creation of both $\mathcal{M}\_{jk}$ and $\mathcal{L}\_{jk}$, scaling as $\sim(N+3N_{\text{probe}} )N$. Nonetheless, the asymptotic scaling does not assure that these contributions are negligible in finite computations, due to both large pre-factors and different durations of the elementary operations. We empirically noticed that this is not the case for finite values of $N$, and expecially when this is comparable to $N_{\text{probe}}$. Although we privileged operations performed in a vectorized fashion, this does not straightforwardly apply to the creation of the matrix $\mathcal{M}\_{jk}$, since the diagonal elements would exhibit infinite values corresponding to the dipole self-energy, which must be removed (16). To perform this task efficiently, we designed a (sigle) loop cycle to run in parallel over several threads. To this aim, we designed the loop with only elementary operations avoiding the possible bottleneck of multiple threads recalling the same complicated function.

Similar considerations can be drawn regarding the memory consumption. The main allocations of RAM are associated to the creation of both the complex-diagonal, dense matrix $\mathcal{M}\_{jk}$ and the probe matrices $\mathcal{J}\_{jk}$. By properly arranging the algorithm, we avoid the allocation of unnecessay memory at the same time, flushing the RAM when desirable. For example, we invert $\mathcal{M}\_{jk}$ at the beginning of the core part and we promptly clear the memory, only conserving the solutions $d\_j$. The probe matrix $\mathcal{J}\_{jk}$ is constructed only afterwards, and if more probe geometries are selected each matrix $\mathcal{J}\_{jk}$ is used and flushed before allocating the next one. Empirically, we observed that the backslash operator "\" allocates $\mathcal{M}\_{jk}$ twice, when inverting it (observed in Julia 1.6). 

We





When the algorithm inverts the M_jk it needs to allocate it twice, so the total memory consumption is given by ∼(2N+3N_probe )N complex Float. 

Nonetheless, one can properly order the operations, and flush the memory when necessary: in particular, we invert M_jk at the beginning and we promptly clear the memory, only conserving the solutions d_j. Since 〖3N〗_probe≲N, the upper bound of memory allocation becomes 〖2×N〗^2. 

By defining the matrix as Complex{Float32} (64 bit) rather than the custom Complex{Float64} (128 bit). We numerically checked that we were operating with enough precision.

### 1.3.1) Physical simplifications
To further simplify the computational problem, we give the user the option to take advantage of some common symmetries in the pysical system. In particular, a typical problem consists of studying the cooperative properties of ordered atomic lattices (26–28). Often, such a system can be arranged to be symmetric for $x\to -x$ and $y\to -y$, without much loss of generality in the physical conclusions. This implies that each dipole $d\_j$ is equal to those at the mirrored positions. The actual degrees of freedom are thus given by the number of atoms satisfying $x\_j\geq 0$ and $y\_j\geq 0$ (roughly $\sim N/4$). The coupled-dipole equations can be then simplified by accounting only for these atoms, and considering as if each of them scattered light from the mirrored positions as well. 




# 2) Code guide
## 2.1) Initializing the code

### 2.1.1) Overall settings

### 2.1.2) Specific settings



#### 2.1.2.1) Disordered geometries settings

#### 2.1.2.2) Array settings

#### 2.1.2.3) Atomic metalens settings

#### 2.1.2.4) Custom geometry settings


## 2.2) Run the simulation

## 2.3) Data outputs



# 3) References 



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



