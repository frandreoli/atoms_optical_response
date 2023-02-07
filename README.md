# 1 - Introduction

The code performs numerical simulations of the steady-state, optical response of large-scale systems of quantum, two-level emitters (atoms, color centers, etc.).
Various options are provided on the density, geometry and physical properties of the light scatterers.

The original code has been developed as a part of a PhD project by Francesco Andreoli, under the supervision of Prof. Dr. Darrick Chang. 
Some parts of this code have been used in \[[1](#Andreoli2021)-[3](#Andreoli2023b)\].

## 1.1 - Code specifications

The core code is written in [Julia](https://julialang.org/) (initially written for Julia 1.6 and well tested on the current version Julia 1.8), while the output data are saved in the HDF5 format. 
A Mathematica notebook is provided to read the data and generate simple plots, with an intuitive, user-friendly interface. This latter is currently available only for the "METALENS" option (see [Sec. 2.1.2](#211---overall-settings)).
The simulation is specifically optimized for multi-threaded computation, and the code was extensively used on a multi-core (>32 threads) cluster with ~400GB of RAM to simulate systems of up to $\sim 7\times 10^5$ atoms.




## 1.2 - Physical background and motivation

In this project, we implement a code to simulate the collective behavior of a mesoscopic ensemble of quantum emitters. This task is crucial to predict new phenomena that can occur in actual experiments, where cooperative effects can drastically alter the optical response \[[4](#Dicke1954CoherenceProcesses)–[21](#Manzoni2018OptimizationArrays)\]. 

We consider two-level emitters with dipole matrix elements $d\_0$, associated to the transition $\omega_0$ from the ground to the excited state. This latter has a spontaneous emission rate $\Gamma\_0$, which defines the elastic scattering rate of the atoms. A supplementary rate $\Gamma'$ (in units of $\Gamma\_0$ accounts for the possibility of inelastic scattering, which might be due to the system environment (such as the excitation of phononic modes). 
A collection of $N$ emitters is characterized by their dipoles $\mathbf{d}\_j$. When illuminated by an external field $\mathbf{E}\_{\text{in}}(\mathbf{r}\_j)$, the optical response can be described by the coupled-dipole formalism, where each dipole is driven by both the input field and the field re-scattered by the other dipoles. This is valid as long as the external drive has a low intensity (with Rabi frequecy $\Omega_0\ll \Gamma_0$), so that non-linear effects are negligible. In these simulations, we assume that the dipoles have a well defined direction $\mathbf{d}\_j = d\_j \hat{\mathbf{d}}\_0$. When dealing with atoms in free space, this is a valid approximation provided that the input field is $\hat{\mathbf{x}}$-polarized and that it retains more intense that the scattere field at the atomic positions. 
On the contrary, a constrain on the dipole orientations is imposed on solid-state emitters by the geometry of the environment (such as the surrounding diamond lattice for color centers in diamond). Finally, we assume that the emitters are in a bulk, dielectric material of (possitive, real) index $n\_{\text{bulk}}$, meaning that $k_0=2\pi/\lambda_0 = n\_{\text{bulk}}\omega_0/c$ denotes the resonant wavevector. Under these conditions, the steady-state response is described by \[[6](#AsenjoGarcia2017ExponentialArrays)\]

 $$ d\_j = \alpha\_0\epsilon\_0 \left(\hat{\mathbf{d}}\_0\cdot \mathbf{E}\_{\text{in}}(\mathbf{r}\_j) + \mu\_0 \omega\_0^2 \displaystyle \sum\_{k\neq j}^N \hat{\mathbf{d}}\_0 \cdot \mathbf{G}(\mathbf{r}\_j-\mathbf{r}\_k,\omega\_0)\cdot \hat{\mathbf{d}}\_0 d_k \right), $$

where $\alpha\_0=-3\pi \epsilon\_0 /\[(\Delta + i(1+\Gamma')/2)k\_0^3\]$ is the atomic polarizability, with $\Delta=(\omega - \omega_0)/\Gamma_0$ representing the dimensionless detuning between the frequency $\omega$ of the incident light and the atomic resonance $\omega_0$, in units of the natural linewidth $\Gamma_0$. At the same time, $\mathbf{G}(\mathbf{r}\_j-\mathbf{r}\_k,\omega\_0)$ is the electromagnetic Green’s tensor. This problem is equivalent to inverting the algebraic system $\sum\_j \mathcal{M}\_{jk} d\_j = \hat{\mathbf{d}}\_0\cdot \mathbf{E}\_{\text{in}}(\mathbf{r}\_j)$, characterized by the $N\times N$ complex, non-Hermitian matrix $\mathcal{M}\_{jk}=\delta\_{jk}/(\alpha\_0\epsilon\_0)-(1-\delta\_{jk})\mu\_0\omega\_0^2 \hat{\mathbf{d}}\_0 \cdot \mathbf{G}(\mathbf{r}\_j-\mathbf{r}\_k,\omega\_0)\cdot \hat{\mathbf{d}}\_0$. 
Similarly, the total field can be reconstructed at $N\_{\text{probe}}$ probe points $\textbf{r}\_{k}$ in space, as $\mathbf{E}(\mathbf{r}\_{k})=\mathbf{E}\_{\text{in}}(\mathbf{r}\_{k})+\sum\_j\mathcal{L}\_{jk}d\_j$, where we define the $3N\_{\text{probe}}\times N$ matrix $\mathcal{L}\_{jk}=\mu\_0\omega\_0^2\mathbf{G}(\mathbf{r}\_j-\mathbf{r}\_k,\omega\_0)\cdot  \hat{\mathbf{d}}\_0 $. 


This code allows to calculate the total field at different probe positions, given many geometric choices of the atomic positions (either ordered or with random distributions). The input field consists of a Gaussian input beam of waist $w\_0$ and focal point $\mathbf{r}=0$, with custom direction and polarization. Moreover, the code calculates the transmission and reflection by projecting the scattered field onto the same Gaussian mode as the input \[[1](#Andreoli2021), [21](#Manzoni2018OptimizationArrays)\]. On top of that, the user has the possibility to calculate the projection and reflection onto a different Gaussian mode as well, with custom waist $w\_{\text{target}}$ and custom focal point (but same direction and polarization as the input). At the same time, the computation can be performed for several values of the detuning $\Delta$, and given any (fixed) choice of the dipole orientation $\hat{\mathbf{d}}\_0$. Finally, when simulating lattice geometries, the user can test scenarios where a fixed fractions of defects (i.e. missing atoms at random lattice positions) is present. For any random setting that is enabled, the computation can be repeated multiple times, each sampling different configurations. 


The possibility of simulating systems with a large number $N\gg 1$ of atom-like emitters is a thrilling idea, especially if $N$ becomes comparable to the typical values of actual experiments. At the same time, the quest for large-scale simulations has also theoretical roots. For example, this is an essential requirement to investigate the bulk properties of dense atomic media, where the effects of the boundary must be negligible. To obtain so, one can illuminate the atomic cloud with a finite light beam whose waist $w\_0$ is smaller than the size of the ensemble $l$. Due to the paraxial approximation, however, this imposes the constraint $\lambda\_0 \ll w\_0\ll l$, which converts to large $N$ for dense ensembles. The number of emitters that we were able to simulate with this code is roughly two orders of magnitude larger than comparable works \[[7](#Javanainen2014ShiftsSample)-[17](#Dobbertin2020CollectiveNanocavities)\].






## 1.3 - Computational insights

The core of our simulations is the inversion of $\mathcal{M}\_{jk}$, implemented with the backslash operator "\\" (see [*left division operator*](https://docs.julialang.org/en/v1/base/math/)). This performs adaptive algorithms based on the structure of the matrix, so that, in the worst-case scenario, the number of elementary operations scales as $\sim N^3$. This process is speeded up through the openBLAS library for linear algebra, which can evaluate the problem in parallel over up to 32 cores. Asymptotically, this task dominates the time complexity of the simulation, overcoming the other sources of time consumption, such as the creation of both $\mathcal{M}\_{jk}$ and $\mathcal{L}\_{jk}$, which scale as $\sim(N+3N_{\text{probe}} )N$. Nonetheless, the asymptotic scaling does not assure that these contributions are negligible in finite computations, due to both large pre-factors and different durations of the elementary operations. We empirically noticed that this is often not the case for reasonable values of $N$, expecially when this is comparable to $\sim N_{\text{probe}}$. Although we privileged operations performed in a vectorized fashion, this does not straightforwardly apply to the creation of the matrix $\mathcal{M}\_{jk}$, since the diagonal elements would exhibit infinite values corresponding to the dipole self-energy, which must be removed \[[6](#AsenjoGarcia2017ExponentialArrays)\]. To perform this task efficiently, we implemented a (sigle) loop cycle to run in parallel over several threads. To this aim, we designed the loop with only elementary operations, avoiding the possible bottleneck of multiple threads calling the same complicated functions.

Similar considerations can be drawn regarding the memory consumption. The main allocations of RAM are associated to the creation of both the complex-diagonal, dense $N\times N$ matrix $\mathcal{M}\_{jk}$ and the $N\_{\text{probe}}\times N$ probe matrices $\mathcal{L}\_{jk}$. By properly arranging the algorithm, we avoid the allocation of unnecessay memory at the same time, flushing the RAM when desirable. For example, we fill in and invert $\mathcal{M}\_{jk}$ at the beginning of the core part and we promptly clear the memory, only conserving the $N$ solutions $d\_j$. The probe matrix $\mathcal{L}\_{jk}$ is constructed only afterwards, and if many probe geometries are selected, thn each matrix $\mathcal{L}\_{jk}$ is filled in, used and flushed before allocating the next one. Empirically, we observed that the backslash operator "\\" allocates $\mathcal{M}\_{jk}$ twice, when inverting it (observed in Julia 1.6), so that we roughly estimate the maximum allocated RAM by $\sim \max(2N,3N\_{\text{probe}} )N$ complex Floats. Finally, we drastically reduce this number by defining these matrices as Complex{Float32} (64 bit) rather than the custom Complex{Float64} (128 bit). We numerically checked that we were operating with enough precision.

### 1.3.1 - Physical simplifications
To further simplify the computational problem, we give the user the option to take advantage of some common symmetries in the pysical system. In particular, a typical problem consists of studying the cooperative properties of ordered atomic lattices (26–28). Often, such a system can be arranged to be symmetric for $x\to -x$ and $y\to -y$, without much loss of generality in the physical conclusions. This implies that each dipole $d\_j$ is equal to those at the mirrored positions. The actual degrees of freedom are thus given by the number of atoms satisfying $x\_j\geq 0$ and $y\_j\geq 0$ (roughly $\sim N/4$). The coupled-dipole equations can be then simplified by accounting only for these atoms, and considering as if each of them scattered light from the mirrored positions as well. 




# 2 - Documentation
In this section, further information will be provided on the use the code, detailing the options and settings available to the user. The core simulation requires an installed version of [Julia](https://julialang.org/) higher than 1.6.


## 2.1 - Initialization
Firts, we describe the settings and options available for the user, which can be found and edited in the file ***Settings.jl***. The options are often provided in the form:

```Julia
option = ["option_1" ; "option_2"][i]
```

and the user can change the integer number $i$ to choose the $i$-th option from the list. 


### 2.1.1 - Overall settings
Here, we describe the generic settings that can be implemented in the code.


- `geometry_settings` \
The option defines the geometry of the positions of the quantum emitters. It can be set to three *disordered* geometries `"DISORDERED_SPHERE"`, `"DISORDERED_CYLINDER"` and `"DISORDERED_CUBOID"`, whose choice will randomly sample the positions from a uniform distribution inside the selected shape. When setting it to `"ARRAYS"` the emitters are arranged on a 3D lattice composed of many arrays in a row, whose number, size, lattice constants and distance can be later set. The choice of `"METALENS"` arranges the emitter position to form a three-layer atomic metalens (further information can be found in [[3](#Andreoli2023b)]). Finally, by selecting `"CUSTOM_POSITIONS"` the user can feed the simulation with its own set of emitter positions (in the 3D space), which must be formatted as a $N\times 3$ matrix and saved in a *.h5* file (HDF5 format) whose name can be later selected.

- `pos_save_option` \
When this option is set to `YES` the simulation will save the atomic positions as a `n_repetitions` $\times N\times 3$ tensor named *r\_atoms*, in the file *atomic\_positions.h5*. The value and meaning of `n_repetitions` will be defined below.

- `defects_fraction`\
Given an ordered geometry (array, metalens or custom) of the atomic positions, the user can remove a fraction $0\leq$ `defects_fraction` $\leq 1$ of atoms, randomly chosen, to simulate the presence of defects in the geometrical construction.

- `small_disorder_std`\
Given an ordered geometry (array, metalens or custom), the user can shift the atomic positions from their standard value by randomly sampling (in each of the three dimensions) a shift from a Gaussian distribution with standard deviation `small_disorder_std`. If `small_disorder_std = 0.0`, then no shift will be applied.

- `mirror_symmetry_option`\
When it is set to `YES`, the code will assume that the atomic positions are symmetric for $x\to -x$ and $y\to -y$, as described in [Sec. 1.3.1](#131---physical-simplifications). 

- `n_repetitions`\
For disordered geometries, or settings that include random sampling, it is often convenient to repeat the simulation many times, each sampling a different configuration. The value of `n_repetitions` defines the number of times the optical response will be calculated.

- `name_simulation`\
This string allows to define a specific label for the simulation, to distinguish it from other runs.

- `RAM_GB_max`\
The maximum amount of RAM that the user would prefer to allocate. The code will try to (roughly) estimate the RAM that it will allocate and abort the simulation if it exceeds `RAM_GB_max`.




### 2.1.2 - Physical settings
Here, we describe the settings related to the physical system. We recall that all lengths are intended in units of $\lambda\_0$ and all rates/frequencies in units of $\Gamma\_0$.

- `dipoles_polarization`\
It defines the dipole matrix element $\hat{\mathbf{d}}\_0$ of the quantum emitters. This must be a $3$-fold unit vector representing the $x$, $y$ and $z$ coordinates, in the form `[x ; y ; z]`.

- `gamma_prime`\
It defines the inelastic scattering rate $\Gamma'$, which quantifies the energy losses from the standpoint of the optical linear response.


- `inhom_broad_std`\
The user can opt for the atoms to have their resonant frequency $\omega\_0$ randomly shifted, by sampling (for each atom) from a Gaussian distribution of mean value $\omega\_0$ and standard deviation `inhom_broad_std` (in units of $\Gamma\_0$). This addresses the possibility of inhomogeneous broadening. If `inhom_broad_std = 0.0`, then all atoms have the same resonant frequency $\omega\_0$.

#### 2.1.2.1 - Input Gaussian beam

- `laser_detunings`\
Artay that defines the dimensionless detuning $\Delta=(\omega-\omega\_0)/\Gamma\_0$ between the frequency of the input beam $\omega$ and the resonance frequency $\omega\_0$, in units of $\Gamma_0$. Many values can be added, in the form `laser_detunings = [value_1 ; value_2 ; value_3 ...]`, and the results will be calculated for each of them

- `w0`\
Beam waist of the input Gaussian beam. For the paraxial approximation to be fully valid, one must have $w\_0\gtrsim 1$ (all in units of $\lambda\_0$).

- `laser_direction`\
Direction of propagation of the input Gaussian beam (unit vector with the same structure of `dipoles_polarization`).

- `field_polarization`\
Polarization of the input Gaussian beam (unit vector with the same structure of `dipoles_polarization`).


### 2.1.3 - Probe settings
Here, we list the settings related to the probe points where the total field (input and scattered) will be calculated.

- `probeXY_option`\
When it is set to `probeXY_option="YES"`, the total field will be calculated on a number `probeXY_points_x*probeXY_points_y` of points on the $(x,y)$ plane at a distance $z=$`probeXY_z`. These points will form a grid of dimensions `probeXY_range_x=[x_min ; x_max]` and `probeXY_range_y=[y_min ; y_max]`. 

- `probeYZ_option` and `probeXZ_option`\
Same as `probeXY_option`, but related to the $(y,z)$ and $(x,z)$ plane.

- `probePLANE_option`\
This option gives the possibility of calculating the probe points on a rectangle on a custom plane, which is orthogonal to the unit vector `probePlane_v3_vec`, and placed at a distance `probePlane_v3_value`.

- `probeSPHERE_option`\
This option allows to calculate the total field on a number `probeSphere_points` of probe points randomly chosen over a sphere of radius `probeSphere_radius`, around the atomic ensemble, by setting it to `FULL_SPHERE`. Otherwise, the values `FORWARD_HEMISPHERE` and `BACKWARD_HEMISPHERE` constrain these probe points on only, respectively, the $z>0$ and $z<0$ hemispheres. Finally, the choice of `NONE` will neglect this option.


### 2.1.3.1 - Target transmission mode
The user has the possibility of calculating the transmission and reflection coefficients not only in the same optical mode as the input Gaussian beam (which is automatic), but also by projecting onto a custom, target mode. This can be accomplished by setting `target_beam_option="YES"`, and then defining the target mode by its waist `w0_target` and focal distance `z0_target`. This target mode will retain the same direction and polarization of the input mode. Finally, in some cases (such as for the `METALENS` geometry), one might want to normalize that mode to carry the same power as the input field (to ensure energy conservation), and the same phase. This is accomplished by setting `normalize_target_option="YES"`.



### 2.1.4 - Specific settings
Some geometrical configurations of the positions of the emitters need further, specific settings, which are described in this section.


#### 2.1.4.1 - Disordered geometries settings
When `geometry_settings="DISORDERED_SPHERE"` or `geometry_settings="DISORDERED_CYLINDER"` or `geometry_settings="DISORDERED_CUBOID"`, then the atomic positions are randomly sampled from a uniform distribution, inside the desired shape, centered at $\mathbf{r}=0$. 

- `dis_atomic_density`\
The number of atoms to be sampled is define through the atomic density (in units of $\lambda\_0^3$, which can be defined with this variable.

- `dis_r_sphere`\
When `geometry_settings="DISORDERED_SPHERE"`, then the shape of the atomic cloud is a 3D ball of radius `dis_r_sphere` (in units of $\lambda\_0$).

- `dis_r_disk` and `dis_z_length`\
When `geometry_settings="DISORDERED_CYLINDER"`, then the shape of the atomic cloud is a cylinder with axis along the $\hat{\mathbf{z}}$ direction, and whose radius and thickness are respectively given by `dis_r_disk` and `dis_z_length` (in units of $\lambda\_0$).

- `dis_x_dim`, `dis_y_dim` and `dis_z_dim`\
When `geometry_settings="DISORDERED_CUBOID"`, then the shape of the atomic cloud is a cuboid with sizes `dis_x_dim`, `dis_y_dim` and `dis_z_dim`.




#### 2.1.4.2 - Array settings
When `geometry_settings="ARRAYS"`, the system arrange the atomic position to into a series of finite 2D, rectangular arrays, extending in the $(x,y)$-plane and separated by a fixed distance in the $\hat{\mathbf{z}}$ direction. The specific settings for this geometry are listed here below.

- `array_n_layers`\
This variable identifies the number of 2D arrays that are placed in series. It must be an integer number with `array_n_layers`>0.

- `array_xi_x` and `array_xi_y`\
These variables define the lattice constants (in units of $\lambda\_0$) of the 2D, rectangular arrays, along the $\hat{\mathbf{x}}$ and $\hat{\mathbf{y}}$ directions.

- `array_xi_z`\
If `array_n_layers`>1, then `array_xi_z` identifies the longitudinal distance between nearest-neighbour 2D arrays (in units of $\lambda\_0$).

- `array_size_x` and `array_size_y`\
These two variables identify the extension of the finite 2D arrays in the $\hat{\mathbf{x}}$ and $\hat{\mathbf{y}}$ directions. They must be written in the following form: `array_size_i = [i_min, i_max]` for $i=x,y$.

- `array_gamma_coop_option` and `array_omega_coop_option`\
These two options allow the user to re-define the array `laser_detunings` with scales dictated by the cooperative properties of a single 2D array, which are more suitable in many computations. Specifically, we consider the cooperative resonance $\omega\_{\text{coop}}$ and decay rate $\Gamma\_{\text{coop}}$ associated to the excitation by the input laser of the 2D atomic spin wave with wavevector $\mathbf{k}\_{\rm in}^{xy}$ in the $(x,y)$-plane \[[6](#AsenjoGarcia2017ExponentialArrays), [19](#Shahmoon2017CooperativeArrays)\]. This process also accounts for the alignment of the atomic dipoles $\hat{\mathbf{d}}\_0$. When `array_gamma_coop_option="YES"`, then the values `laser_detunings` will be scaled by a factor $\Gamma\_{\text{coop}}$. Similarly, when `array_omega_coop_option = "YES"`, then `laser_detunings` will be shifted so that the zero corresponds to $\omega\_{\text{coop}}$.


#### 2.1.4.3 - Atomic metalens settings
When `geometry_settings="METALENS"`, the system arrange the atomic position to act as a metalens, following the recipe of \[[3](#Andreoli2023b)\]. This metalens extends in the $(x,y)$ plane, and it is meant to focus a Gaussian beam at normal incidence (i.e. travelling in the $\hat{\mathbf{z}}$ direction). 

- `focal_length`\
Its values defines the focal length $f$ of the lens. We recall that the actual distance between the lens and the focal spot is close to $f$. 

- `r_lens`\
This parameter defines the radius of the atomic metalens (in units of $\lambda_0$).

- `disks_width`\
It corresponds to the width $r\_{j+1}-r\_j$ of each disk composing the atomic metalens (in units of $\lambda_0$).

- `buffer_smooth`\
It describes the fraction of `disks_width` that is covered by a buffer zone, as detailed in [3](#Andreoli2023b).

- `phase_shift`\
Its value defines the overal phase shift impressed by the lens, i.e. the phase shift $\phi\_0$ at the radial center $(x,y)=(0,0)$.

- `default_probe_option`\
When it is set to `YES`, then some of the probe-point parameters will be overwritten by default settings. For example, the value of `probeXY_z` is set to correspond to the focal plane, while also imposing `probeYZ_x = probeXZ_y = 0.0`.

- `default_target_option`\
By setting this option to `YES` the user has the possibility of overwriting the specifications of the target Gaussian beam (as defined in [Sec. 2.1.3.1](#2131---target-transmission-mode) with those of the theoretical mode where an ideal thin lens with the same focal length would transmit.




#### 2.1.4.4 - Custom geometry settings
When `geometry_settings="CUSTOM_POSITIONS"`, then code will simulate the optical response of a system of quantum emitters whose positions in space are chosen by the user. These must be saved as a $N\times 3$ real matrix in a file .h5 with the HDF5 encoding. The folder where the file is located can be specified via `custom_pos_folder`, while the name of the file and the name of the variable with the relevant matrix are identified by respectively `custom_pos_file` and `custom_pos_variable`.




## 2.2 - Run the simulation
To run the simulation, the user must first define the correct settings in the file ***Settings.jl***, following the guidelines of [Sec. 2.1](#21---initialization). Afterwards, the simulation can be performed by simply running the file ***Launcher.jl***, within Julia (the code has been written and tested for versions 1.6-1.8). Attention must be paid to have already installed the package [HDF5](https://juliaio.github.io/HDF5.jl/stable/).

### 2.2.1 - Linux environment
If the simulation is run on a Linux operating system, then a handy wrapper is provided in the file ***Shell_Launcher.sh***. This can be run from the command line by typing `bash Shell_Launcher.sh`. This launcher provides the possibility of controlling the number of cores (used for operations of linear algebra) and threads (used to fill in the matrices with the data) that the simulation will make use of. The command-line outputs of the simulation will be stored in a specific file ".out" in the directory ***Data_Output/***.

Through this launcher, the user can also choose to monitor the RAM and CPUs used by the main simulation by setting `yesORnoMonitor="y"`. This option will launch an additional code (contained in ***Resources\_monitor/Resources\_monitor.jl***), which is aimed to keep track of the usage of these computational resources. This supplementary tool will produce two plots (in png format) with either the RAM or CPU used as a function of time. It requires the Julia package [Plots](https://docs.juliaplots.org/stable/) and it has been tested on the operating system *CentOS Linux 7 (Core)*.



## 2.3 - Data outputs
The results of the simulation are stored in different output .h5 files, that we list here below. The data are have a HDF5 encoding, which allows to store many variables in the same file. 

### 2.3.1 - Generic parameters and options

- ***options.h5***\
This file stores the values of the options defined by the user. Specifically, it contains the variables `pos_save_option`, `geometry_settings`, `target_beam_option`, `probeXY_option`, `probeYZ_option`, `probeXZ_option`, `probePLANE_option`, `probeSPHERE_option` and `mirror_symmetry_option`. If `target_beam_option=="YES"`, then the option `normalize_target_option` is stored too.

- ***settings.h5***\
This files stores the parameters of the simulations: `w0`, `gamma_prime`, `inhom_broad_std`, `laser_detunings`, `laser_direction`, `field_polarization`, `defects_fraction`, `n_repetitions`, `probePlane_vec` and `dipoles_polarization`. If `target_beam_option=="YES"`, then the variables `w0_target` and `z0_target` are stored too.

- ***atomic_positions.h5***\
If `pos_save_option=="YES"`, then the atomic positions are stored in this file, under the variable name `r_atoms`. This consist of a `n_repetitions` $\times N\times 3$ tensor, where the first index defines the repetition, the second is a label for each atom and the third identifies the spatial coordinate (x, y or z). When random configurations are enabled (such as disordered geometries) each repetition is associated to a different set of atomic positions.

### 2.3.2 - Specific parameters and options
Specific settings are also saved, related to each different chosen geometry.

#### 2.3.2.1 - Disordered geometries
For disordered geometries, the specific initializing parameters are stored in the file ***settings_disorder.h5***, which always contains the variable `dis_atomic_density`. If `geometry_settings == "DISORDERED_SPHERE"` then also the variable `dis_r_sphere` is saved. If `geometry_settings == "DISORDERED_CYLINDER"` then the two variables `dis_r_disk` and `dis_z_length` are stored. Finally, if `geometry_settings == "DISORDERED_CUBOID"` the file contains the variables: `dis_x_dim`, `dis_y_dim` and `dis_z_dim`.


#### 2.3.2.2 - Array geometry

#### 2.3.2.3 - Atomic metalens


### 2.3.3 - Simulation results


``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``,



 ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``, ``,


# 3 - References 


<a id="Andreoli2021">[1]</a> 
Andreoli F, Gullans MJ, High AA, Browaeys A, Chang DE, 
*Maximum Refractive Index of an Atomic Medium*, 
[Physical Review X 11, 011026](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.11.011026) (2021)

 <a id="Andreoli2023a">[2]</a> 
Andreoli F, Windt B, Grava S, Andolina GM, Gullans MJ, High AA, Chang DE, 
*The maximum refractive index of an atomic crystal - from quantum optics to quantum chemistry*, 
paper in preparation (2023)

<a id="Andreoli2023b">[3]</a> 
Andreoli F, High AA, Chang DE, 
*Metalens formed by structured, sub-wavelength atomic arrays*, 
paper in preparation (2023)

<a id="Dicke1954CoherenceProcesses">[4]</a>
Dicke R H,
*Coherence in Spontaneous Radiation Processes*,
[Physical Review 93, 99-110](https://dx.doi.org/10.1103/PhysRev.93.99) (1954)

<a id="Gross1982SuperradianceEmission">[5]</a>
Gross M, Haroche S,
*Superradiance: An essay on the theory of collective spontaneous emission*,
[Physics Reports 93, 301-396](https://dx.doi.org/10.1016/0370-1573(82)90102-8) (1982)

<a id="AsenjoGarcia2017ExponentialArrays">[6]</a>
Asenjo-Garcia A, Moreno-Cardoner M, Albrecht A, Kimble H J, Chang D E,
*Exponential Improvement in Photon Storage Fidelities Using Subradiance and "Selective Radiance" in Atomic Arrays*,
[Physical Review X 7, 31024](https://dx.doi.org/10.1103/PhysRevX.7.031024) (2017)

<a id="Javanainen2014ShiftsSample">[7]</a>
Javanainen Juha, Ruostekoski Janne, Li Yi, Yoo Sung Mi,
*Shifts of a resonance line in a dense atomic sample*,
[Physical Review Letters 112, ](https://dx.doi.org/10.1103/PhysRevLett.112.113603) (2014)

<a id="Schilder2016PolaritonicAtoms">[8]</a>
Schilder N J, Sauvan C, Hugonin J-P, Jennewein S, Sortais Y R P, Browaeys A, Greffet J-J,
*Polaritonic modes in a dense cloud of cold atoms*,
[Physical Review A 93, 063835](https://dx.doi.org/10.1103/PhysRevA.93.063835) (2016)

<a id="Schilder2017HomogenizationScatterers">[9]</a>
Schilder N J, Sauvan C, Sortais Y R P, Browaeys A, Greffet J-J,
*Homogenization of an ensemble of interacting resonant scatterers*,
[Physical Review A 96, 013825](https://dx.doi.org/10.1103/PhysRevA.96.013825) (2017)

<a id="Schilder2020NearResonantAtoms">[10]</a>
Schilder N J, Sauvan C, Sortais Y R P, Browaeys A, Greffet J-J,
*Near-Resonant Light Scattering by a Subwavelength Ensemble of Identical Atoms*,
[Physical Review Letters 124, 073403](https://dx.doi.org/10.1103/PhysRevLett.124.073403) (2020)

<a id="Chomaz2012AbsorptionAnalysis">[11]</a>
Chomaz L, Corman L, Yefsah T, Desbuquois R, Dalibard J,
*Absorption imaging of a quasi-two-dimensional gas: a multiple scattering analysis*,
[New Journal of Physics 14, 055001](https://dx.doi.org/10.1088/1367-2630/14/5/055001) (2012)

<a id="Javanainen2016LightOptics">[12]</a>
Javanainen Juha, Ruostekoski Janne,
*Light propagation beyond the mean-field theory of standard optics*,
[Optics Express 24, 993](https://dx.doi.org/10.1364/oe.24.000993) (2016)

<a id="Zhu2016LightMedia">[13]</a>
Zhu Bihui, Cooper John, Ye Jun, Rey Ana Maria,
*Light scattering from dense cold atomic media*,
[Physical Review A 94, 023612](https://dx.doi.org/10.1103/PhysRevA.94.023612) (2016)

<a id="Jennewein2016PropagationAtoms">[14]</a>
Jennewein S, Sortais Y RP, Greffet J J, Browaeys A,
*Propagation of light through small clouds of cold interacting atoms*,
[Physical Review A 94, ](https://dx.doi.org/10.1103/PhysRevA.94.053828) (2016)

<a id="Corman2017TransmissionAtoms">[15]</a>
Corman L, Ville J L, Saint-Jalm R, Aidelsburger M, Bienaim'e T, Nascimb`ene S, Dalibard J, Beugnon J,
*Transmission of near-resonant light through a dense slab of cold atoms*,
[Physical Review A 96, 53629](https://dx.doi.org/10.1103/PhysRevA.96.053629) (2017)

<a id="Jenkins2016CollectiveExperiment">[16]</a>
Jenkins S D, Ruostekoski J, Javanainen J, Jennewein S, Bourgain R, Pellegrino J, Sortais Y R P, Browaeys A,
*Collective resonance fluorescence in small and dense atom clouds: Comparison between theory and experiment*,
[Physical Review A 94, 023842](https://dx.doi.org/10.1103/PhysRevA.94.023842) (2016)

<a id="Dobbertin2020CollectiveNanocavities">[17]</a>
Dobbertin Helge, Low Robert, Scheel Stefan,
*Collective dipole-dipole interactions in planar nanocavities*,
[Physical Review A 102, 031701](https://dx.doi.org/10.1103/PhysRevA.102.031701) (2020)

<a id="Bettles2016EnhancedArray">[18]</a>
Bettles Robert J, Gardiner Simon A, Adams Charles S,
*Enhanced Optical Cross Section via Collective Coupling of Atomic Dipoles in a 2D Array*,
[Physical Review Letters 116, 103602](https://dx.doi.org/10.1103/PhysRevLett.116.103602) (2016)

<a id="Shahmoon2017CooperativeArrays">[19]</a>
Shahmoon Ephraim, Wild Dominik S, Lukin Mikhail D, Yelin Susanne F,
*Cooperative Resonances in Light Scattering from Two-Dimensional Atomic Arrays*,
[Physical Review Letters 118, 113601](https://dx.doi.org/10.1103/PhysRevLett.118.113601) (2017)

<a id="Rui2020ALayer">[20]</a>
Rui Jun, Wei David, Rubio-Abadal Antonio, Hollerith Simon, Zeiher Johannes, Stamper-Kurn Dan M, Gross Christian, Bloch Immanuel,
*A subradiant optical mirror formed by a single structured atomic layer*,
[Nature 583, 369-374](https://dx.doi.org/10.1038/s41586-020-2463-x) (2020)

<a id="Manzoni2018OptimizationArrays">[21]</a>
Manzoni M T, Moreno-Cardoner M, Asenjo-Garcia A, Porto J V, Gorshkov A V, Chang D E,
*Optimization of photon storage fidelity in ordered atomic arrays*,
[New Journal of Physics 20, 83048](https://dx.doi.org/10.1088/1367-2630/aadb74) (2018)





