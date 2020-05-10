subroutine build_impedance_matrix(w_l,w_x,w_y,w_z,k,a,Zmn)
    use KernelFunctions
    implicit NONE
    double precision, INTENT(IN), dimension(:) :: w_l,w_x,w_y,w_z
    double precision, INTENT(IN) :: k,a
    double complex, INTENT(INOUT), dimension(:,:) :: Zmn
    double complex :: Amn,Phimn
    integer :: i_n, i_m, n
    n = size(w_l)-2
    call set_kernel_parameters(w_l,w_x,w_y,w_z,k,a)
    do i_m = 1,n
        do i_n = 1,n
            call calc_A_m_n(i_m+1,i_n+1,Amn)
            call calc_Phi_m_n(i_m+1,i_n+1,Phimn)
            Zmn(i_m,i_n) = (k**2)*Amn-Phimn
        end do
    end do
end subroutine build_impedance_matrix

subroutine build_sources(w_l,w_x,w_y,w_z,k,a,i_source,V)
    use KernelFunctions
    implicit NONE
    double precision, INTENT(IN), dimension(:) :: w_l, w_x,w_y,w_z
    double precision, INTENT(IN) :: k,a
    integer, INTENT(IN) :: i_source
    double complex, INTENT(INOUT), dimension(:) :: V
    double precision :: Vm
    integer :: i_m, n
    n = size(w_z)-2
    call set_kernel_parameters(w_l,w_x,w_y,w_z,k,a)
    do i_m = 1,n
        call calc_V_m(i_m, i_source, Vm)
        V(i_m) = Vm
    end do
end subroutine build_sources

