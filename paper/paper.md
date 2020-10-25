---
title: 'rustbca: A High-Performance Binary-Collision-Approximation Code for Ion-Material Interactions'
tags:
  - Rust
  - plasma material interactions
  - binary collision approximation
  - ion solid interactions
  - sputtering
  - reflection
  - implantation
authors:
  - name: Jon. T Drobny
    orcid: 0000-0002-9733-6058
    affiliation: 1
  - name: Davide Curreli
    affiliation: 1
affiliations:
  - name: Department of Nuclear, Plasma, and Radiological Engineering, University of Illinois at Urbana-Champaign
    index: 1
date: 24 October 2020
bibliography: paper.bib
---

# Summary

Ion-material interactions are of vital importance in industrial applications, the study and design of nuclear fusion devices, engineering survivable spacecraft, and more. In particular, plasma-material interactions, including the phenomena of sputtering, reflection, and implantation are dominated by ion-material interactions. These phenomena are difficult to model analytically, and many such models rely on empirical or semi-empirical formulas, such as the Yamamura sputtering yield formula[1], or the Thomas reflection coefficient[2]. However, such models are inherently limited, and of little applicability to complex geometry, multi-component surfaces, or for coupling to plasma or material dynamics codes. Since ion-material interactions span a range of energies from sub-eV to GeV and beyond, n-body approaches such as molecular dynamics are computationally infeasible for most applications. Instead, approximations to the full n-body problem can be used. Most commonly, the Binary Collision Approximation (BCA), a set of simplifying assumptions to the full n-body ion-material interaction are used. rustbca is a high-performance, general purpose ion-material interactions BCA code, built for performance, flexibility, and ease of use. rustbca includes 2D, inhomogeneous, arbitrary composition geometry, electronic stopping formulations for low energy (up to 25 keV/nucleon) and high energy (up to 1 GeV/nucleon), Kr-C, ZBL, Moliere, and Lenz-Jensen screened coulomb interatomic potentials, Lennard-Jones and Morse attractive-repulsive potentials, the unique capability of using multiple interatomic potentials in a single simulation, choice of Gaussian quadrature and the MAGIC algorithm for determining the scattering angle, full trajectory tracking of ions and material atoms, human-readable configuration file, and full 6D output of all particles that leave the simulation (via sputtering or reflection).

# Binary Collision Approximation

rustbca is an amorphous-material BCA code. The TRIM family of codes, which includes Tridyn, SDTrimSP, F-TRIDYN, and SRIM, are also amorphous-material BCA codes, and historically the most popular implementation of the BCA. Based on the number of citations, SRIM is the most popular amorphous-material BCA code, followed by the original TRIM, Tridyn, and SDTrimSP. Crystalline-material BCA codes have also been developed, such as MARLOWE and ACAT, but are not as widely used. Assumptions used in the amorphous-material BCA can be summarized as follows:

* Particles in the code, ions and material atoms both, are "superparticles" that represent many real ions or atoms each.
* Energetic particles interact with stationary atoms in the material through elastic, binary collisions
* Collisions occur at mean-free-path lengths
* Particle trajectories are approximated by the classical asymptotic trajectories
* Electronic interactions are separated from the nuclear, elastic interactions
* Local electronic energy losses occur at each collision
* Nonlocal electronic energy losses occur along each segment of the asymptotic trajectories
* Material atoms are mobile and transfer momentum following all collisions
* Particles are stopped when their energy drops below a threshold, cutoff energy
* Particles that leave a surface experience reflection by or refraction through a locally planar surface binding potential

# Examples

rustbca includes multiple example input files, under the examples/ folder on the directory, as well as discussion of each on the rustbca github wiki page.

# Statement of Need

Many BCA codes have been developed. Many of these, including SRIM, Tridyn, F-TRIDYN, SDTrimSP and its derivatives, FTRIM and VFTRIM, are based on the original TRIM code. However, each has limitations that prevent widespread adoption across a broad range of applications. In particular, SRIM, which is closed-source, suffers from relatively poor performance and significant anomalies in sputtered atom angular distributions and light ion sputtering yields. Tridyn and F-TRIDYN, which are not open source, are limited to low energy, screened-coulomb potentials, mono-angular ion beams, atomically flat and atomically rough surfaces respectively, and are single-threaded. SDTrimSP, although significantly more developed than the preceding codes, is valid only for low energy and not open source. As far as the authors are aware, there is no widely-accepted open-source BCA code. Additionally, those that are available are not well suited to a wide range of problems, including direct coupling of BCA codes to particle and subsurface dynamics codes, such as those performed for ITER divertor simulations. rustbca has been developed to fill that gap and expand upon the feature set included in currently available BCA codes. Features unique to rustbca include the ability to handle attractive-repulsive interatomic potentials, including  the ability to use multiple interatomic potentials in one simulation, handling high-energy incident ions and 2D geometry, large file input of incident particles to facilitate coupling to other codes via HDF5, a human-readable, TOML configuration file, modern programming techniques, robust error-handling, and multi-threading. We hope that giving users direct access to a user-friendly, flexible, high-performance, open-source BCA will encourage and enable exciting research avenues in industrial plasmas, nuclear fusion, astrophysics, and more.
