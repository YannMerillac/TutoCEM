module KernelFunctions
    implicit NONE
    double precision, PARAMETER :: Pi = 3.1415927
    double complex, PARAMETER :: j = (0,1)
    double precision, private :: k,a
    double precision, dimension(:), allocatable, private :: ctrl_l
    double precision, private :: z_m,z_n
    integer, private :: m,n,i_source
    ! integration parameters
    integer, PARAMETER :: limit=50
    double precision , PARAMETER :: epsabs=1.e-8
    double precision , PARAMETER :: epsrel=0.
    integer, dimension(limit) :: iwork1
    double precision, dimension(4*limit) :: work1
    integer, dimension(limit) :: iwork2
    double precision, dimension(4*limit) :: work2
    
contains
subroutine set_m_n(m_in,n_in)
    integer, INTENT(IN) :: m_in,n_in
    m = m_in
    n = n_in
end subroutine set_m_n

subroutine set_kernel_parameters(ctrl_l_in,k_in,a_in)
    double precision, INTENT(IN) :: k_in,a_in
    double precision, INTENT(IN), dimension(:) :: ctrl_l_in
    k = k_in
    a = a_in
    ctrl_l = ctrl_l_in
end subroutine set_kernel_parameters

function G(t1,t2)
    double precision :: r
    double complex :: G
    double precision, INTENT(IN) :: t1,t2
    r = sqrt((t1-t2)**2+a**2)
    G = exp(-k*r*j)/(4.*Pi*r)
end function G

subroutine get_ctrl_pts(i,l1,l2,l3)
    integer, INTENT(IN) :: i
    double precision, INTENT(INOUT) :: l1,l2,l3
    integer :: n
    n = size(ctrl_l)
    l2 = ctrl_l(i)
    if (i==1) then
        l1 = ctrl_l(1)
        l3 = ctrl_l(2)
    else if (i==n) then
        l1 = ctrl_l(n-1)
        l3 = ctrl_l(n)
    else
        l1 = ctrl_l(i-1)
        l3 = ctrl_l(i+1)
    end if
end subroutine get_ctrl_pts

function tri_f(t,i)
    double precision :: tri_f,z1,z2,z3
    double precision, INTENT(IN) :: t
    integer, INTENT(IN) :: i
    call get_ctrl_pts(i,z1,z2,z3)
    if (t<=z1 .or. t>=z3) then
        tri_f = 0.
    else if (t<=z2) then
        tri_f = (t-z1)/(z2-z1)
    else
        tri_f = (z3-t)/(z3-z2)
    end if
end function tri_f

function tri_df(t,i)
    double precision :: tri_df,z1,z2,z3
    double precision, INTENT(IN) :: t
    integer, INTENT(IN) :: i
    call get_ctrl_pts(i,z1,z2,z3)
    if (t<=z1 .or. t>=z3) then
        tri_df = 0.
    else if (t<=z2) then
        tri_df = 1./(z2-z1)
    else
        tri_df = -1/(z3-z2)
    end if
end function tri_df

function Ia_1_re(zp)
    double precision, INTENT(IN) :: zp
    double complex :: Ia_1
    double precision :: Ia_1_re
    Ia_1 = G(z_m,zp)*tri_f(zp,n)
    Ia_1_re = dreal(Ia_1)
end function Ia_1_re

function Ia_1_im(zp)
    double precision, INTENT(IN) :: zp
    double complex :: Ia_1
    double precision :: Ia_1_im
    Ia_1 = G(z_m,zp)*tri_f(zp,n)
    Ia_1_im = dimag(Ia_1)
end function Ia_1_im

function Ia_re(z)
    double precision, INTENT(IN) :: z
    double precision :: Ia_re,z1,z2,z3,output
    double precision :: abserr
    integer :: neval, ier, last
    z_m = z
    call get_ctrl_pts(n,z1,z2,z3)
    call dqags(Ia_1_re,z1,z3,epsabs,epsrel,output,abserr,neval,ier,limit,4*limit,last,iwork1,work1)
    Ia_re = output*tri_f(z_m,m)
end function Ia_re

function Ia_im(z)
    double precision, INTENT(IN) :: z
    double precision :: Ia_im,z1,z2,z3,output
    double precision :: abserr
    integer :: neval, ier, last
    z_m = z
    call get_ctrl_pts(n,z1,z2,z3)
    call dqags(Ia_1_im,z1,z3,epsabs,epsrel,output,abserr,neval,ier,limit,4*limit,last,iwork1,work1)
    Ia_im = output*tri_f(z_m,m)
end function Ia_im

function Iphi_1_re(zp)
    double precision, INTENT(IN) :: zp
    double complex :: Iphi_1
    double precision :: Iphi_1_re
    Iphi_1 = G(z_m,zp)*tri_df(zp,n)
    Iphi_1_re = dreal(Iphi_1)
end function Iphi_1_re

function Iphi_1_im(zp)
    double precision, INTENT(IN) :: zp
    double complex :: Iphi_1
    double precision :: Iphi_1_im
    Iphi_1 = G(z_m,zp)*tri_df(zp,n)
    Iphi_1_im = dimag(Iphi_1)
end function Iphi_1_im

function Iphi_re(z)
    double precision, INTENT(IN) :: z
    double precision :: Iphi_re,z1,z2,z3,output
    double precision :: abserr
    integer :: neval, ier, last
    z_m = z
    call get_ctrl_pts(n,z1,z2,z3)
    call dqags(Iphi_1_re,z1,z3,epsabs,epsrel,output,abserr,neval,ier,limit,4*limit,last,iwork1,work1)
    Iphi_re = output*tri_df(z_m,m)
end function Iphi_re

function Iv(z)
    double precision, INTENT(IN) :: z
    double precision :: Iv,z1,z2,dz
    z1 = ctrl_l(i_source)
    z2 = ctrl_l(i_source+1)
    dz = z2-z1
    if (z>=z1 .and. z<=z2) then
        Iv=(1./dz)*tri_f(z,m)
    else
        Iv=0.
    end if
end function Iv

function Iphi_im(z)
    double precision, INTENT(IN) :: z
    double precision :: Iphi_im,z1,z2,z3,output
    double precision :: abserr
    integer :: neval, ier, last
    z_m = z
    call get_ctrl_pts(n,z1,z2,z3)
    call dqags(Iphi_1_im,z1,z3,epsabs,epsrel,output,abserr,neval,ier,limit,4*limit,last,iwork1,work1)
    Iphi_im = output*tri_df(z_m,m)
end function Iphi_im

subroutine calc_A_m_n(m_in,n_in,A_m_n)
    integer, INTENT(IN) :: m_in,n_in
    double complex, INTENT(INOUT) :: A_m_n
    double precision :: zm1,zm2,zm3, Amn_re, Amn_im
    double precision :: abserr
    integer :: neval, ier, last
    m = m_in
    n = n_in
    ! integrate
    call get_ctrl_pts(m,zm1,zm2,zm3)
    call dqags(Ia_re,zm1,zm3,epsabs,epsrel,Amn_re,abserr,neval,ier,limit,4*limit,last,iwork2,work2)
    call dqags(Ia_im,zm1,zm3,epsabs,epsrel,Amn_im,abserr,neval,ier,limit,4*limit,last,iwork2,work2)
    A_m_n = Amn_re+j*Amn_im
end subroutine calc_A_m_n

subroutine calc_Phi_m_n(m_in,n_in,A_m_n)
    integer, INTENT(IN) :: m_in,n_in
    double complex, INTENT(INOUT) :: A_m_n
    double precision :: zm1,zm2,zm3, Amn_re, Amn_im
    double precision :: abserr
    integer :: neval, ier, last
    m = m_in
    n = n_in
    ! integrate
    call get_ctrl_pts(m,zm1,zm2,zm3)
    call dqags(Iphi_re,zm1,zm3,epsabs,epsrel,Amn_re,abserr,neval,ier,limit,4*limit,last,iwork2,work2)
    call dqags(Iphi_im,zm1,zm3,epsabs,epsrel,Amn_im,abserr,neval,ier,limit,4*limit,last,iwork2,work2)
    A_m_n = Amn_re+j*Amn_im
end subroutine calc_Phi_m_n

subroutine calc_V_m(m_in, i_source_in, Vm)
    integer, INTENT(IN) :: m_in, i_source_in
    double precision, INTENT(INOUT) :: Vm
    double precision :: abserr, zm1, zm2, zm3
    integer :: neval, ier, last
    m = m_in
    i_source = i_source_in
    call get_ctrl_pts(m,zm1,zm2,zm3)
    call dqags(Iv,zm1,zm3,epsabs,epsrel,Vm,abserr,neval,ier,limit,4*limit,last,iwork2,work2)
end subroutine calc_V_m

end module KernelFunctions
