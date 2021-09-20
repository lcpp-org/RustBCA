module rustbca

    use, intrinsic :: iso_c_binding
    
    interface

        subroutine reflect_single_ion_c(num_species_target, ux, uy, uz, E1, &
            Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, Eb2, n2) bind(c)

            !Runs a single ion BCA trajectory with no recoils
            !Args:
            !   num_species_target (integer): number of species in target      
            !   ux (real): x-direction of incident ion, x-direction of reflected ion
            !   uy (real): y-direction of incident ion, y-direction of reflected ion
            !   uz (real): z-direction of incident ion, z-direction of reflected ion
            !   E1 (real): initial energy of incident ion in eV
            !   Z1 (real): atomic number of incident ion
            !   m1 (real): atomic mass of incident ion in eV
            !   Ec1 (real): cutoff energy of incident ion in eV
            !   Es1 (real): surface binding energy of incident ion in eV
            !   Z2 (array(real)): list of atomic numbers of target speciesd
            !   m2 (array(real)): list of atomic masses of target species in amu
            !   Ec2 (array(real)): list of cutoff energies of target species in eV
            !   Es2 (array(real)): list of surface binding energies of target species in eV
            !   Eb2 (array(real)): list of bulk binding energies of target species in eV
            !   n2 (array(real)): list of number densities of target species in 1/angstrom^3

            use, intrinsic :: iso_c_binding
            real(c_double), intent(inout) :: ux, uy, uz, E1
            real(c_double), intent(in) :: Z1, m1, Ec1, Es1
            integer(c_int), intent(in) :: num_species_target
            real(c_double), intent(in), dimension(*) :: Z2, m2, Ec2, Es2, Eb2, n2

        end subroutine reflect_single_ion_c

        function compound_bca_list_fortran(num_incident_ions, track_recoils, ux, uy, uz, E1, &
            Z1, m1, Ec1, Es1, &
            num_species_target, Z2, m2, Ec2, Es2, Eb2, n2, &
            num_emitted_particles) bind(c) result(output)

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
            !   ux (real): x-direction
            !   uy (real): y-direction
            !   uz (real): z-direction
            !   alpha (real): local surface angle measured counter-clockwise from x-axis in radians

            real(8), intent(inout) :: ux, uy, uz
            real(8), intent(in) :: alpha

            ux = ux*cos(alpha) - uy*sin(alpha)
            uy = ux*sin(alpha) + uy*cos(alpha)

        end subroutine transform_to_local_angle

end module rustbca