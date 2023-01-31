### atoms_optical_response

# Introduction

The code performs numerical simulations of the steady-state, optical response of large-scale systems of quantum, two-level emitters (atoms, color centers, etc.).
Various options are provided on the density, geometry and physical properties of the light scatterers.

The original code has been developed as a part of a PhD project by Francesco Andreoli, under the supervision of Prof. Dr. Darrick Chang.

## Code specifications

The core code is written in Julia (tested on version Julia 1.8), while the output data are saved in the HDF5 format. 
A Mathematica notebook is provided to read the data and generate simple plots, with an intuitive, user-friendly interface. This latter is currently available only for the "Metalens" option (see below).
The simulation is specifically optimized for multi-threaded computation, and the code was extensively used on a multi-core (>32 threads) cluster with ~400GB of RAM to simulate systems of $\sim 5\times 10^5$ atoms.

## Physical background




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




---
output:
  md_document:
    variant: markdown_github
bibliography: bibliography.bib
---


