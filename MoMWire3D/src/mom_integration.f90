module KernelFunctions
    implicit NONE
    double precision, PARAMETER :: Pi = 3.1415927
    double complex, PARAMETER :: j = (0,1)
    double precision, private :: k,a
    double precision, dimension(:), allocatable, private :: ctrl_l, ctrl_x, ctrl_y, ctrl_z
    double precision, private :: l_m
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

subroutine set_kernel_parameters(ctrl_l_in,ctrl_x_in,ctrl_y_in,ctrl_z_in,k_in,a_in)
    double precision, INTENT(IN) :: k_in,a_in
    double precision, INTENT(IN), dimension(:) :: ctrl_l_in, ctrl_x_in, ctrl_y_in, ctrl_z_in
    k = k_in
    a = a_in
    ctrl_l = ctrl_l_in
    ctrl_x = ctrl_x_in
    ctrl_y = ctrl_y_in
    ctrl_z = ctrl_z_in
end subroutine set_kernel_parameters

function interval(t)
    double precision, INTENT(IN) :: t
    integer :: i
    integer :: interval
    do i=1,size(ctrl_l)-1
        if (t>=ctrl_l(i) .and. t<ctrl_l(i+1)) then
            interval=i
            exit
        end if
    end do
end function interval

subroutine get_r(t,rx,ry,rz)
    double precision, INTENT(IN) :: t
    double precision, INTENT(OUT) :: rx, ry, rz
    double precision :: alpha,t1,t2
    integer :: i_pos
    i_pos = interval(t)
    t1 = ctrl_l(i_pos)
    t2 = ctrl_l(i_pos+1)
    alpha = (t-t1)/(t2-t1)
    rx = alpha*ctrl_x(i_pos+1)+(1-alpha)*ctrl_x(i_pos)
    ry = alpha*ctrl_y(i_pos+1)+(1-alpha)*ctrl_y(i_pos)
    rz = alpha*ctrl_z(i_pos+1)+(1-alpha)*ctrl_z(i_pos)
end subroutine get_r

subroutine get_t(t,tx,ty,tz)
    double precision, INTENT(IN) :: t
    double precision, INTENT(OUT) :: tx, ty, tz
    double precision :: dl
    integer :: i_pos
    i_pos = interval(t)
    dl = ctrl_l(i_pos+1) - ctrl_l(i_pos)
    tx = (ctrl_x(i_pos+1)-ctrl_x(i_pos))/dl
    ty = (ctrl_y(i_pos+1)-ctrl_y(i_pos))/dl
    tz = (ctrl_z(i_pos+1)-ctrl_z(i_pos))/dl
end subroutine get_t

function G(t1,t2)
    double precision :: r,x1,y1,z1,x2,y2,z2
    double complex :: G
    double precision, INTENT(IN) :: t1,t2
    call get_r(t1,x1,y1,z1)
    call get_r(t2,x2,y2,z2)
    r = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2+a**2)
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

function Ia_1_re_x(lp)
    double precision, INTENT(IN) :: lp
    double complex :: Ia_1
    double precision :: Ia_1_re_x,txp,typ,tzp
    call get_t(lp,txp,typ,tzp)
    Ia_1 = G(l_m,lp)*tri_f(lp,n)*txp
    Ia_1_re_x = dreal(Ia_1)
end function Ia_1_re_x

function Ia_1_re_y(lp)
    double precision, INTENT(IN) :: lp
    double complex :: Ia_1
    double precision :: Ia_1_re_y,txp,typ,tzp
    call get_t(lp,txp,typ,tzp)
    Ia_1 = G(l_m,lp)*tri_f(lp,n)*typ
    Ia_1_re_y = dreal(Ia_1)
end function Ia_1_re_y

function Ia_1_re_z(lp)
    double precision, INTENT(IN) :: lp
    double complex :: Ia_1
    double precision :: Ia_1_re_z,txp,typ,tzp
    call get_t(lp,txp,typ,tzp)
    Ia_1 = G(l_m,lp)*tri_f(lp,n)*tzp
    Ia_1_re_z = dreal(Ia_1)
end function Ia_1_re_z

function Ia_1_im_x(lp)
    double precision, INTENT(IN) :: lp
    double complex :: Ia_1
    double precision :: Ia_1_im_x,txp,typ,tzp
    call get_t(lp,txp,typ,tzp)
    Ia_1 = G(l_m,lp)*tri_f(lp,n)*txp
    Ia_1_im_x = dimag(Ia_1)
end function Ia_1_im_x

function Ia_1_im_y(lp)
    double precision, INTENT(IN) :: lp
    double complex :: Ia_1
    double precision :: Ia_1_im_y,txp,typ,tzp
    call get_t(lp,txp,typ,tzp)
    Ia_1 = G(l_m,lp)*tri_f(lp,n)*typ
    Ia_1_im_y = dimag(Ia_1)
end function Ia_1_im_y

function Ia_1_im_z(lp)
    double precision, INTENT(IN) :: lp
    double complex :: Ia_1
    double precision :: Ia_1_im_z,txp,typ,tzp
    call get_t(lp,txp,typ,tzp)
    Ia_1 = G(l_m,lp)*tri_f(lp,n)*tzp
    Ia_1_im_z = dimag(Ia_1)
end function Ia_1_im_z

function Ia_re(l)
    double precision, INTENT(IN) :: l
    double precision :: Ia_re, a_re_x, a_re_y, a_re_z, l1, l2, l3, tx, ty, tz
    double precision :: abserr
    integer :: neval, ier, last
    l_m = l
    call get_t(l,tx,ty,tz)
    call get_ctrl_pts(n,l1,l2,l3)
    call dqags(Ia_1_re_x,l1,l3,epsabs,epsrel,a_re_x,abserr,neval,ier,limit,4*limit,last,iwork1,work1)
    call dqags(Ia_1_re_y,l1,l3,epsabs,epsrel,a_re_y,abserr,neval,ier,limit,4*limit,last,iwork1,work1)
    call dqags(Ia_1_re_z,l1,l3,epsabs,epsrel,a_re_z,abserr,neval,ier,limit,4*limit,last,iwork1,work1)
    Ia_re = tri_f(l_m,m)*(a_re_x*tx+a_re_y*ty+a_re_z*tz)
end function Ia_re

function Ia_im(l)
    double precision, INTENT(IN) :: l
    double precision :: Ia_im, a_im_x, a_im_y, a_im_z, l1, l2, l3, tx,ty,tz
    double precision :: abserr
    integer :: neval, ier, last
    l_m = l
    call get_t(l,tx,ty,tz)
    call get_ctrl_pts(n,l1,l2,l3)
    call dqags(Ia_1_im_x,l1,l3,epsabs,epsrel,a_im_x,abserr,neval,ier,limit,4*limit,last,iwork1,work1)
    call dqags(Ia_1_im_y,l1,l3,epsabs,epsrel,a_im_y,abserr,neval,ier,limit,4*limit,last,iwork1,work1)
    call dqags(Ia_1_im_z,l1,l3,epsabs,epsrel,a_im_z,abserr,neval,ier,limit,4*limit,last,iwork1,work1)
    Ia_im = tri_f(l_m,m)*(a_im_x*tx+a_im_y*ty+a_im_z*tz)
end function Ia_im

function Iphi_1_re(lp)
    double precision, INTENT(IN) :: lp
    double complex :: Iphi_1
    double precision :: Iphi_1_re
    Iphi_1 = G(l_m,lp)*tri_df(lp,n)
    Iphi_1_re = dreal(Iphi_1)
end function Iphi_1_re

function Iphi_1_im(lp)
    double precision, INTENT(IN) :: lp
    double complex :: Iphi_1
    double precision :: Iphi_1_im
    Iphi_1 = G(l_m,lp)*tri_df(lp,n)
    Iphi_1_im = dimag(Iphi_1)
end function Iphi_1_im

function Iphi_re(l)
    double precision, INTENT(IN) :: l
    double precision :: Iphi_re,l1,l2,l3,output
    double precision :: abserr
    integer :: neval, ier, last
    l_m = l
    call get_ctrl_pts(n,l1,l2,l3)
    call dqags(Iphi_1_re,l1,l3,epsabs,epsrel,output,abserr,neval,ier,limit,4*limit,last,iwork1,work1)
    Iphi_re = output*tri_df(l_m,m)
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

function Iphi_im(l)
    double precision, INTENT(IN) :: l
    double precision :: Iphi_im,l1,l2,l3,output
    double precision :: abserr
    integer :: neval, ier, last
    l_m = l
    call get_ctrl_pts(n,l1,l2,l3)
    call dqags(Iphi_1_im,l1,l3,epsabs,epsrel,output,abserr,neval,ier,limit,4*limit,last,iwork1,work1)
    Iphi_im = output*tri_df(l_m,m)
end function Iphi_im

subroutine calc_A_m_n(m_in,n_in,A_m_n)
    integer, INTENT(IN) :: m_in,n_in
    double complex, INTENT(INOUT) :: A_m_n
    double precision :: lm1,lm2,lm3, Amn_re, Amn_im
    double precision :: abserr
    integer :: neval, ier, last
    m = m_in
    n = n_in
    ! integrate
    call get_ctrl_pts(m,lm1,lm2,lm3)
    call dqags(Ia_re,lm1,lm3,epsabs,epsrel,Amn_re,abserr,neval,ier,limit,4*limit,last,iwork2,work2)
    call dqags(Ia_im,lm1,lm3,epsabs,epsrel,Amn_im,abserr,neval,ier,limit,4*limit,last,iwork2,work2)
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
