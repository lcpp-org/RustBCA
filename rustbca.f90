module rustbca
    
    interface

        subroutine reflect_single_ion_c(num_species_target, ux, uy, uz, E1, &
            Z1, m1, Ec1, Es1, Z2, m2, Ec2, Es2, Eb2, n2) bind(c)

            use, intrinsic :: iso_c_binding
            real(c_double), intent(inout) :: ux, uy, uz, E1
            real(c_double), intent(in) :: Z1, m1, Ec1, Es1
            integer(c_int), intent(in) :: num_species_target
            real(c_double), intent(in), dimension(*) :: Z2, m2, Ec2, Es2, Eb2, n2

        end subroutine reflect_single_ion_c
        
    end interface

end module rustbca