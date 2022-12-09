module global_variable

    implicit none
    !------------------
    ! パラメーター
    !------------------
    
    !smac
    integer, parameter :: width = 100 
    integer, parameter :: height = 50 
    integer, parameter :: frames = 5000 
    real(8), parameter :: Lx = 8.0d0
    real(8), parameter :: Ly = 4.0d0
    real(8), parameter :: h = Lx / dble(width)
    real(8), parameter :: Re = 150.0d0
    real(8), parameter :: dt = 0.02d0
    real(8), parameter :: u_inflow = 1.0d0
    logical, parameter :: is_restart = .false.

    !bicgstab
    real(8), parameter :: eps = 1e-2

    !ibm
    integer, parameter :: Nx = 21, Ny = 21
    integer, parameter :: N_prism = 2*Nx + 2*Ny - 4
    real(8), parameter :: D = 1.0d0
    real(8), parameter :: prism_left = 5.5d0
    real(8), parameter :: prism_right = prism_left + D
    real(8), parameter :: prism_under = 5.5d0
    real(8), parameter :: prism_above = prism_under + D
    integer, parameter :: i_left =  int(prism_left /h) + 1
    integer, parameter :: i_right = int(prism_right/h) + 1
    integer, parameter :: j_under = int(prism_under/h) + 1
    integer, parameter :: j_above = int(prism_above/h) + 1

    real(8), parameter :: center_x = 2.0d0 
    real(8), parameter :: center_y = 2.0d0 
    real(8), parameter :: radius = 0.25d0
    integer, parameter :: N_circle =20 

    ! 円の時と角柱の時で書き換えること
    integer, parameter :: force_num = N_circle


    !-------------------
    ! 変数
    !-------------------

    !smac関係
    real(8), save ::       u(0:width+1, 0:height+1) = 0.0d0 !u(n)
    real(8), save ::   u_old(0:width+1, 0:height+1) = 0.0d0 !u(n-1)
    real(8), save ::   u_new(0:width+1, 0:height+1) = 0.0d0 !u(n+1)
    real(8), save ::       v(0:width+1, 0:height+1) = 0.0d0 !v(n)
    real(8), save ::   v_old(0:width+1, 0:height+1) = 0.0d0 !v(n-1)
    real(8), save ::   v_new(0:width+1, 0:height+1) = 0.0d0 !v(n+1)
    real(8), save ::       p(0:width+1, 0:height+1) = 0.0d0 !p(n)
    real(8), save ::     phi(0:width+1, 0:height+1) = 0.0d0 !ポアソン方程式の解

    real(8), save ::   u_sol(0:width+1, 0:height+1) = 0.0d0
    real(8), save ::   v_sol(0:width+1, 0:height+1) = 0.0d0
    
    !ibm関係(defined at staggerd grid)
    real(8), save :: u_tilde(1:width, 1:height) = 0.0d0
    real(8), save :: v_tilde(1:width, 1:height) = 0.0d0
    real(8), save :: fx(1:width, 1:height) = 0.0d0
    real(8), save :: fy(1:width, 1:height) = 0.0d0

    !ibm関係(defined at prism surface)
    real(8), save :: x_f(1:force_num) = 0.0d0
    real(8), save :: y_f(1:force_num) = 0.0d0
    real(8), save :: u_tilde_in_surface(1:force_num) = 0.0d0
    real(8), save :: v_tilde_in_surface(1:force_num) = 0.0d0
    real(8), save :: fx_in_surface(1:force_num) = 0.0d0
    real(8), save :: fy_in_surface(1:force_num) = 0.0d0

    !bicgstab
    real(8), save :: a1(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: a2(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: a3(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: a4(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: a5(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: b_crank_x(0:width+1, 0:height+1) = 0.0d0

    real(8), save :: c1(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: c2(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: c3(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: c4(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: c5(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: b_crank_y(0:width+1, 0:height+1) = 0.0d0

    real(8), save :: e1(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: e2(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: e3(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: e4(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: e5(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: b_poisson(0:width+1, 0:height+1) = 0.0d0

    real(8), save :: r(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: r0(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: stab1(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: stab2(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: stab3(0:width+1, 0:height+1) = 0.0d0
    real(8), save :: stab4(0:width+1, 0:height+1) = 0.0d0


    !------------------------
    ! アウトプット
    !------------------------
    character(len=256), parameter :: file_name = "output.dat"
    integer, parameter :: file_number = 10 !(適当w)

    character(len=256), parameter :: file_name_save_output = "save_output.dat"
    integer, parameter :: file_number_save_output =  17 !(適当w)

    character(len=256), parameter :: file_name_save_input = "save_input.dat"
    integer, parameter :: file_number_save_input =  20 !(適当w)

    character(len=256), parameter :: file_name_point = "output_point.dat"
    integer, parameter :: file_number_point = 23 
end module global_variable
