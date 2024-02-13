module rustbca

    use, intrinsic :: iso_c_binding

    interface

        subroutine reflect_single_ion_c(num_species_target, ux, uy, uz, E1, &
            Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, Eb2, n2) bind(c)

            !Runs a single ion BCA trajectory with no recoils
            !Args:
            !   num_species_target (integer): number of species in target
            !   ux (real(c_double)): x-direction of incident ion, x-direction of reflected ion
            !   uy (real(c_double)): y-direction of incident ion, y-direction of reflected ion
            !   uz (real(c_double)): z-direction of incident ion, z-direction of reflected ion
            !   E1 (real(c_double)): initial energy of incident ion in eV
            !   Z1 (real(c_double)): atomic number of incident ion
            !   m1 (real(c_double)): atomic mass of incident ion in eV
            !   Ec1 (real(c_double)): cutoff energy of incident ion in eV
            !   Es1 (real(c_double)): surface binding energy of incident ion in eV
            !   Z2 (real(c_double), dimension(:)): list of atomic numbers of target speciesd
            !   m2 (real(c_double), dimension(:)): list of atomic masses of target species in amu
            !   Ec2 (real(c_double), dimension(:)): list of cutoff energies of target species in eV
            !   Es2 (real(c_double), dimension(:)): list of surface binding energies of target species in eV
            !   Eb2 (real(c_double), dimension(:)): list of bulk binding energies of target species in eV
            !   n2 (real(c_double), dimension(:)): list of number densities of target species in 1/angstrom^3

            use, intrinsic :: iso_c_binding
            real(c_double), intent(inout) :: ux, uy, uz, E1
            real(c_double), intent(in) :: Z1, m1, Ec1, Es1
            integer(c_int), intent(in) :: num_species_target
            real(c_double), intent(in), dimension(*) :: Z2, m2, Ec2, Es2, Eb2, n2

        end subroutine reflect_single_ion_c

        subroutine rotate_given_surface_normal(nx, ny, nz, ux, uy, uz) bind(c)

          use, intrinsic :: iso_c_binding
          real(c_double), intent(in) :: nx, ny, nz
          real(c_double), intent(inout) :: ux, uy, uz

        end subroutine rotate_given_surface_normal

        subroutine rotate_back(nx, ny, nz, ux, uy, uz) bind(c)

          use, intrinsic :: iso_c_binding
          real(c_double), intent(in) :: nx, ny, nz
          real(c_double), intent(inout) :: ux, uy, uz

        end subroutine rotate_back

        function compound_bca_list_fortran(num_incident_ions, track_recoils, ux, uy, uz, E1, &
            Z1, m1, Ec1, Es1, &
            num_species_target, Z2, m2, Ec2, Es2, Eb2, n2, &
            num_emitted_particles) bind(c) result(output)


            !Runs a homogeneous, flat, compound target BCA with an arbitrary list of ions.
            !Args:
            !   num_incident_ion (integer(c_int)): number of incident ions
            !   track_recoils (logical(c_bool)): whether to generate recoils (disable to turn off sputtering)
            !   ux (real(c_double), dimension(:)): x-direction of incident ion, x-direction of reflected ion
            !   uy (real(c_double), dimension(:)): y-direction of incident ion, y-direction of reflected ion
            !   uz (real(c_double), dimension(:)): z-direction of incident ion, z-direction of reflected ion
            !   E1 (real(c_double), dimension(:)): initial energy of incident ion in eV
            !   Z1 (real(c_double), dimension(:)): atomic number of incident ion
            !   m1 (real(c_double), dimension(:)): atomic mass of incident ion in eV
            !   Ec1 (real(c_double), dimension(:)): cutoff energy of incident ion in eV
            !   num_species_target(integer(c_int)): number of species in target
            !   Es1 (real(c_double), dimension(:)): surface binding energy of incident ion in eV
            !   Z2 (real(c_double), dimension(:)): list of atomic numbers of target speciesd
            !   m2 (real(c_double), dimension(:)): list of atomic masses of target species in amu
            !   Ec2 (real(c_double), dimension(:)): list of cutoff energies of target species in eV
            !   Es2 (real(c_double), dimension(:)): list of surface binding energies of target species in eV
            !   Eb2 (real(c_double), dimension(:)): list of bulk binding energies of target species in eV
            !   n2 (real(c_double), dimension(:)): list of number densities of target species in 1/angstrom^3
            !   num_emitted_particles (integer(c_int), intent(out)): NOTE THAT THIS IS INTENT(OUT) number of emitted particles in output
            !Returns:
            !   output (type(c_ptr)): a c pointer to a 2D array of size (num_emitted_particles, 6) that consists of Z, m, E, ux, uy, uz

            use, intrinsic :: iso_c_binding
            logical(c_bool), intent(in) :: track_recoils
            integer(c_int), intent(in) :: num_incident_ions, num_species_target
            integer(c_int), intent(out) :: num_emitted_particles
            real(c_double), intent(in), dimension(*) :: ux, uy, uz, E1, Z1, m1, Ec1, Es1
            real(c_double), intent(in), dimension(*) :: Z2, m2, Ec2, Es2, EB2, n2
            type(c_ptr) :: output
        end function compound_bca_list_fortran

    end interface

    contains

        subroutine transform_to_local_angle(ux, uy, uz, alpha)

            !Rotates a vector in 2D
            !Args:
            !   ux (real(c_double)): x-direction
            !   uy (real(c_double)): y-direction
            !   uz (real(c_double)): z-direction
            !   alpha (real(c_double)): local surface angle measured counter-clockwise from x-axis in radians

            real(8), intent(inout) :: ux, uy, uz
            real(8), intent(in) :: alpha

            ux = ux*cos(alpha) - uy*sin(alpha)
            uy = ux*sin(alpha) + uy*cos(alpha)

        end subroutine transform_to_local_angle

end module rustbca
