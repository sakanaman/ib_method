module intf_mod

    interface
        !-----------------------
        !離散化に関するサブルーチン
        !-----------------------

        !圧力項のx方向の離散化
        subroutine discrete_pressure_x(i, j, p , width, height, dx, dpdx)
            integer, intent(in) :: i, j, width, height
            real(8), intent(in) :: dx
            real(8), intent(in) :: p(0 : width+1, 0 : height+1) !形状明示配列
            real(8), intent(out) ::  dpdx
        end subroutine discrete_pressure_x

        !圧力項のy方向の離散化
        subroutine discrete_pressure_y(i, j, p , width, height, dy, dpdy)
            integer, intent(in) :: i, j, width, height
            real(8), intent(in) :: dy
            real(8), intent(in) :: p(0 : width+1, 0 : height+1) !形状明示配列
            real(8), intent(out) :: dpdy
        end subroutine discrete_pressure_y

        !対流項のx方向への離散化
        subroutine discrete_advection_x(i, j, u, v, dx, dy, width, height, ad_x)
            integer, intent(in) :: i, j, width, height
            real(8), intent(in) :: dx, dy
            real(8), intent(in) :: u(0 : width+1, 0 : height+1), v(0 : width+1, 0 : height+1)
            real(8), intent(out) :: ad_x
        end subroutine discrete_advection_x

        !対流項のy方向への離散化
        subroutine discrete_advection_y(i, j, u, v, dx, dy, width, height, ad_y)
            integer, intent(in) :: i, j, width, height
            real(8), intent(in) :: dx, dy
            real(8), intent(in) :: u(0 : width+1, 0 : height+1), v(0 : width+1, 0 : height+1)
            real(8), intent(out) :: ad_y
        end subroutine discrete_advection_y

        !粘性項のx方向への離散化
        subroutine discrete_viscosity_x(i, j, u, v, dx, dy, Re, width, height, visc_x)
            integer, intent(in) :: i, j, width, height
            real(8), intent(in) :: u(0 : width+1, 0 : height+1),v(0 : width+1, 0 : height+1)
            real(8), intent(in) :: dx, dy, Re
            real(8), intent(out) :: visc_x            
        end subroutine discrete_viscosity_x

        !粘性項のy方向への離散化
        subroutine discrete_viscosity_y(i, j, u, v, dx, dy, Re, width, height, visc_y)

            integer, intent(in) :: i, j, width, height
            real(8), intent(in) :: u(0 : width+1, 0 : height+1), v(0 : width+1, 0 : height+1)
            real(8), intent(in) :: dx, dy, Re
            real(8), intent(out) :: visc_y

        end subroutine discrete_viscosity_y


        !---------------------------
        !楕円型偏微分方程式のソルバー
        !---------------------------

        !SOR法で解く(for 楕円型偏微分方程式)
        subroutine solve_sor(x, a1, a2, a3, a4, a5, b ,eps, alpha, width, height)
            integer, intent(in) :: width,height
            real(8), intent(in) :: eps, alpha
            real(8), intent(in) :: a1(0:width+1, 0:height+1)
            real(8), intent(in) :: a2(0:width+1, 0:height+1)
            real(8), intent(in) :: a3(0:width+1, 0:height+1)
            real(8), intent(in) :: a4(0:width+1, 0:height+1)
            real(8), intent(in) :: a5(0:width+1, 0:height+1)
            real(8),intent(in) :: b(0:width+1, 0:height+1)
            real(8), intent(out) :: x(0:width+1, 0:height+1)
        end subroutine solve_sor

        !Bi-CGStab法で解く(for 楕円型偏微分方程式)
        subroutine solve_bicgstab(x, a1, a2, a3, a4, a5, b ,eps, width, height, r, r0, p, e, y, v)
            integer, intent(in) :: width,height
            real(8), intent(in) :: eps
            ! a1~a5,bは調整済みとする(境界条件など考慮している形になっているということ)
            real(8), intent(in) :: a1(0:width+1, 0:height+1)
            real(8), intent(in) :: a2(0:width+1, 0:height+1)
            real(8), intent(in) :: a3(0:width+1, 0:height+1)
            real(8), intent(in) :: a4(0:width+1, 0:height+1)
            real(8), intent(in) :: a5(0:width+1, 0:height+1)
            real(8),intent(in) :: b(0:width+1, 0:height+1)
            real(8), intent(out) :: x(0:width+1, 0:height+1)

            real(8),intent(inout) :: r(0:width+1, 0:height+1)
            real(8),intent(inout) :: r0(0:width+1, 0:height+1)
            real(8),intent(inout) :: p(0:width+1, 0:height+1)
            real(8),intent(inout) :: e(0:width+1, 0:height+1)
            real(8),intent(inout) :: y(0:width+1, 0:height+1)
            real(8),intent(inout) :: v(0:width+1, 0:height+1)
        end subroutine solve_bicgstab

        !-------------------------
        !境界条件に関するサブルーチン
        !-------------------------

        subroutine set_upper_boundary_slip(width ,height, u, v)
            integer, intent(in) :: width, height
            real(8), intent(out) :: u(0:width+1, 0:height+1)
            real(8), intent(out) :: v(0:width+1, 0:height+1)
        end subroutine set_upper_boundary_slip

        subroutine set_upper_boundary_slip_for_vnew(width, height, v_new)
            integer, intent(in) :: width, height
            real(8), intent(out) :: v_new(0:width+1, 0:height+1)
        end subroutine set_upper_boundary_slip_for_vnew

        subroutine set_under_boundary_slip(width, height, u, v)
            integer, intent(in) :: width, height
            real(8), intent(inout) :: u(0:width+1, 0:height+1)
            real(8), intent(inout) :: v(0:width+1, 0:height+1)
        end subroutine set_under_boundary_slip

        subroutine set_under_boundary_slip_vnew(width, height, v_new)
            integer, intent(in) :: width, height
            real(8), intent(out) :: v_new(0:width+1, 0:height+1)
        end subroutine set_under_boundary_slip_vnew


        subroutine set_inflow_boundary(width, height, u_inflow, u, v)
            integer, intent(in) :: width, height
            real(8), intent(in) :: u_inflow !流入速度
            real(8), intent(out) :: u(0:width+1, 0:height+1)
            real(8), intent(out) :: v(0:width+1, 0:height+1)
        end subroutine set_inflow_boundary

        subroutine set_inflow_boundary_for_unew(width, height, u_inflow, u_new)
            integer, intent(in) :: width, height
            real(8), intent(in) :: u_inflow !流入速度
            real(8), intent(out) :: u_new(0:width+1, 0:height+1)
        end subroutine set_inflow_boundary_for_unew

        subroutine set_outflow_boundary(width, height, dt, dx, u, v, u_new, v_new)
            integer, intent(in) :: width, height
            real(8), intent(in) :: dt, dx
            real(8), intent(in) :: u(0:width+1, 0:height+1)
            real(8), intent(in) :: v(0:width+1, 0:height+1)
            real(8), intent(out) :: u_new(0:width+1, 0:height+1)
            real(8), intent(out) :: v_new(0:width+1, 0:height+1)
        end subroutine set_outflow_boundary

        subroutine set_phi_boundary(width, height, phi)
            integer, intent(in) :: width, height
            real(8), intent(out) :: phi(0:width+1, 0:height+1)
        end subroutine set_phi_boundary

        !----------------------
        ! 初期条件の設定
        !----------------------
        subroutine intial_condition(width, height, u_inflow, u, v, u_old, v_old, p)
            integer, intent(in) :: width, height
            real(8), intent(in) :: u_inflow
            real(8), intent(out) :: u(0: width+1, 0:height+1)
            real(8), intent(out) :: v(0: width+1, 0:height+1)
            real(8), intent(out) :: u_old(0: width+1, 0:height+1) ! for アダムス・バッシュフォース
            real(8), intent(out) :: v_old(0: width+1, 0:height+1) ! for アダムス・バッシュフォース
            real(8), intent(out) :: p(0: width+1, 0:height+1)
        end subroutine intial_condition


        subroutine restart_initialize(width, height, u_inflow, u, v, u_old, v_old, p, file_name, file_number)
            integer, intent(in) :: width, height
            real(8), intent(in) :: u_inflow
            character(len=256), intent(in) :: file_name
            integer, intent(in) :: file_number

            real(8), intent(out) :: u(0: width+1, 0:height+1)
            real(8), intent(out) :: v(0: width+1, 0:height+1)
            real(8), intent(out) :: u_old(0: width+1, 0:height+1) ! for アダムス・バッシュフォース
            real(8), intent(out) :: v_old(0: width+1, 0:height+1) ! for アダムス・バッシュフォース(以下、AB法と呼ぶ)
            real(8), intent(out) :: p(0: width+1, 0:height+1)

        end subroutine restart_initialize


        !-----------------------------
        ! 埋め込み境界法に関するサブルーチン
        !-----------------------------

        subroutine set_force_points_in_prism(x_left, x_right, y_under, y_above, Nx, Ny, x_f, y_f)
            real(8), intent(in) :: x_left, x_right  ! 角柱左側x座標, 右側x座標
            real(8), intent(in) :: y_under, y_above  ! 角柱下側y座標, 上側y座標
            integer, intent(in) :: Nx, Ny ! 各軸方向ごとの体積力点を置く個数
            real(8), intent(out) :: x_f(1 : Nx*2 + Ny*2 - 4)
            real(8), intent(out) :: y_f(1 : Nx*2 + Ny*2 - 4)
        end subroutine set_force_points_in_prism

        subroutine set_force_points_in_circle(center_x, center_y, radius, N_circle, x_f, y_f, &
                                                x_min, x_max, y_min, y_max)
            real(8), intent(in) :: center_x, center_y
            real(8), intent(in) :: radius
            integer, intent(in) :: N_circle
            real(8), intent(out) :: x_f(1 : N_circle), y_f(1 : N_circle)
            real(8), intent(out) :: x_min, x_max, y_min, y_max
        end subroutine set_force_points_in_circle

        subroutine delta_h_1d(x,x0,h,res)
            real(8), intent(in) :: x, x0, h
            real(8), intent(out) :: res
        end subroutine delta_h_1d


        !-----------------------------
        ! 楕円型偏微分方程式の調整
        !-----------------------------
        subroutine crank_nicolson_left_side_x(width, height, h, Re, dt, &
                                     a1, a2, a3, a4, a5)
            integer, intent(in) :: width, height
            real(8), intent(in) :: h, Re, dt
            real(8), intent(out) :: a1(0:width+1, 0:height+1)
            real(8), intent(out) :: a2(0:width+1, 0:height+1)
            real(8), intent(out) :: a3(0:width+1, 0:height+1)
            real(8), intent(out) :: a4(0:width+1, 0:height+1)
            real(8), intent(out) :: a5(0:width+1, 0:height+1)
        end subroutine crank_nicolson_left_side_x

        subroutine crank_nicolson_left_side_y(width, height, h, Re, dt, &
                                    c1, c2, c3, c4, c5)
            integer, intent(in) :: width, height
            real(8), intent(in) :: h, Re, dt
            real(8), intent(out) :: c1(0:width+1, 0:height+1)
            real(8), intent(out) :: c2(0:width+1, 0:height+1)
            real(8), intent(out) :: c3(0:width+1, 0:height+1)
            real(8), intent(out) :: c4(0:width+1, 0:height+1)
            real(8), intent(out) :: c5(0:width+1, 0:height+1)
        end subroutine crank_nicolson_left_side_y

        subroutine poisson_left_side(width, height, h, e1, e2, e3, e4, e5)
            integer, intent(in) :: width, height
            real(8), intent(in) :: h
            real(8), intent(out) :: e1(0:width+1, 0:height+1)
            real(8), intent(out) :: e2(0:width+1, 0:height+1)
            real(8), intent(out) :: e3(0:width+1, 0:height+1)
            real(8), intent(out) :: e4(0:width+1, 0:height+1)
            real(8), intent(out) :: e5(0:width+1, 0:height+1)
        end subroutine poisson_left_side


        !------------------
        !　ユーティリティ
        !------------------
        subroutine calc_u_point(i, j, h, res_x, res_y)
            integer,intent(in) :: i, j
            real(8), intent(in) :: h
            real(8), intent(out) :: res_x, res_y
        end subroutine calc_u_point

        subroutine calc_v_point(i, j, h, res_x, res_y)
            integer,intent(in) :: i, j
            real(8), intent(in) :: h
            real(8), intent(out) :: res_x, res_y
        end subroutine calc_v_point

        !---------------------
        ! 記録サブルーチン
        !---------------------
        subroutine initialize_output_file
        end subroutine initialize_output_file

        subroutine log_velocity
        end subroutine log_velocity

        subroutine initialize_outputfile_for_point
        end subroutine initialize_outputfile_for_point

        subroutine log_velocity_point(x, y)
            real(8), intent(in) :: x, y
        end subroutine log_velocity_point


        subroutine savedata(width, height, u_inflow, u, v, u_old, v_old, p, file_number, file_name)
            integer, intent(in) :: width, height
            real(8), intent(in) :: u_inflow
            character(len=256), intent(in) :: file_name
            integer, intent(in) :: file_number
            real(8), intent(out) :: u(0: width+1, 0:height+1)
            real(8), intent(out) :: v(0: width+1, 0:height+1)
            real(8), intent(out) :: u_old(0: width+1, 0:height+1) ! for アダムス・バッシュフォース
            real(8), intent(out) :: v_old(0: width+1, 0:height+1) ! for アダムス・バッシュフォース(以下、AB法と呼ぶ)
            real(8), intent(out) :: p(0: width+1, 0:height+1)
        end subroutine savedata

        subroutine close_file
        end subroutine close_file

    end interface

end module intf_mod
