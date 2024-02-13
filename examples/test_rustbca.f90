
program test_rustbca

    use rustbca
    use, intrinsic :: iso_c_binding

    integer :: N_ions
    real(c_double), allocatable, dimension(:) :: ux, uy, uz, E, Z1, m1, Ec1, Es1
    integer(c_int) :: num_species_target, num_emitted_particles
    real(c_double), target :: Z2(2), m2(2), Ec2(2), Es2(2), Eb2(2), n2(2)
    real(c_double) :: ux1, uy1, uz1, E1
    type(c_ptr) :: bca_output_c
    real(c_double), pointer, dimension(:,:) :: bca_output_f
    real :: start, stop
    logical(c_bool) :: track_recoils

    !Initial ion conditions
    N_ions = 1000
    allocate(ux(N_ions), uy(N_ions), uz(N_ions), E(N_ions), Z1(N_ions), m1(N_ions), Ec1(N_ions), Es1(N_ions))
    ux(:) = 0.999
    uy(:) = sqrt(1.0 - 0.999*0.999)
    uz(:) = 0.0
    E(:) = 1000.0_8

    !Hydrogen
    Z1(:) = 1.0_8
    m1(:) = 1.008_8
    Ec1(:) = 1.0_8
    Es1(:) = 1.5_8

    !Titanium Hydride
    num_species_target = 2
    Z2(1) = 22.0_8
    m2(1) = 47.867_8
    Ec2(1) = 4.84_8
    Es2(1) = 4.84_8
    Eb2(1) = 3.0_8
    n2(1) = 0.04527_8

    Z2(2) = 1.0_8
    m2(2) = 1.008_8
    Ec2(2) = 1.5_8
    Es2(2) = 1.5_8
    Eb2(2) = 0.0_8
    n2(2) = 0.09054_8

    track_recoils = .false.

    call cpu_time(start)
    bca_output_c = compound_bca_list_fortran(N_ions, track_recoils, ux, uy, uz, E, &
        Z1, m1, Ec1, Es1, &
        num_species_target, Z2, m2, Ec2, Es2, Eb2, n2, &
        num_emitted_particles)
    call c_f_pointer(bca_output_c, bca_output_f, [num_emitted_particles, 6])

    call cpu_time(stop)

    write(*,*) "Elapsed time in seconds per ion per eV: ", (stop - start)/N_ions/1000.0

    deallocate(bca_output_f)

    call cpu_time(start)
    do i = 0, N_ions
        !Test reflect_single_ion routine
        ux1 = 0.999
        uy1 = sqrt(1.0 - 0.999*0.999)
        uz1 = 0.0
        E1 = 1000.0
        call reflect_single_ion_c(num_species_target, ux1, uy1, uz1, E1, Z1(1), m1(1), Ec1(1), Es1(1), Z2, m2, Ec2, Es2, Eb2, n2)
    end do
    call cpu_time(stop)
    write(*,*) "Elapsed time in ions per eV per s: ", (stop - start)/N_ions/1000.0

    deallocate(ux, uy, uz, E, Z1, m1, Ec1, Es1)

end program test_rustbca
