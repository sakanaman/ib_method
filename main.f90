program main
    use intf_mod
    use global_variable
    implicit none

    !-------------------
    ! 補助的な変数
    !-------------------
    real(8) :: aabb_left
    real(8) :: aabb_right
    real(8) :: aabb_under
    real(8) :: aabb_above

    integer :: aabb_left_i
    integer :: aabb_right_i
    integer :: aabb_under_j
    integer :: aabb_above_j

    real(8) :: pressure_term_x
    real(8) :: pressure_term_y

    real(8) :: advection_term_x
    real(8) :: advection_term_y

    real(8) :: advection_term_x_old
    real(8) :: advection_term_y_old

    real(8) :: viscosity_term_x
    real(8) :: viscosity_term_y

    real(8) :: delta_h_x, delta_h_y
    real(8) :: time = 0.0d0
    integer :: i, j, k, iter

    integer(8) :: force_index_x, force_index_y

    real(8) :: u_point_x, u_point_y
    real(8) :: v_point_x, v_point_y

    real(8) :: modify_phi_term

    call initialize_outputfile_for_point
    call initialize_output_file

    ! 初期速度場、圧力場の設定
    call intial_condition(width, height, u_inflow, u, v, u_old, v_old, p)


    ! 体積力の定義点を設定
    ! call set_force_points_in_prism(prism_left, prism_right, prism_under, prism_above, Nx, Ny, x_f, y_f)
    call set_force_points_in_circle(center_x, center_y, radius, N_circle, x_f, y_f, &
                                    aabb_left, aabb_right, aabb_under, aabb_above)
    
    ! 円一個の時のみ有効。アドホックだがコメントアウトでの調整をせよ。
    do j = 1, height
    do i = 1, width
        call calc_u_point(i, j, h, u_point_x, u_point_y)
        call calc_v_point(i, j, h, v_point_x, v_point_y)
        if((u_point_x - center_x)**2.0d0 + (u_point_y - center_y)**2.0d0 <= radius**2.0d0) then
            u(i, j) = 0.0d0
        endif
        if((v_point_x - center_x)**2.0d0 + (v_point_y - center_y)**2.0d0 <= radius**2.0d0) then
            v(i, j) = 0.0d0
        endif
    enddo
    enddo

    write(*, *) "radius = ", radius
    write(*, *) "Re = ", re
    write(*, *) "frames = ", frames
    write(*, *) "center_x = ", center_x
    write(*, *) "center_y = ", center_y
    write(*, *) "force_num = ", force_num
    write(*, *) "aabb_left = ", aabb_left
    write(*, *) "aabb_right = ", aabb_right
    write(*, *) "aabb_under = ", aabb_under
    write(*, *) "aabb_above = ", aabb_above

    aabb_left_i = int(aabb_left / h) + 1
    aabb_right_i = int(aabb_right / h) + 1
    aabb_under_j = int(aabb_under / h) + 1
    aabb_above_j = int(aabb_above / h) + 1

    !ポアソン方程式の左辺の係数調整
    call poisson_left_side(width, height, h, e1, e2, e3, e4, e5)

    !クランクニコルソンの左辺項の設定(x)
    call crank_nicolson_left_side_x(width, height, h, Re, dt, &
                                        a1, a2, a3, a4, a5)

    !クランクニコルソンの左辺項の設定(y)
    call crank_nicolson_left_side_y(width, height, h, Re, dt, &
                                         c1, c2, c3, c4, c5)


    if (is_restart .eqv. .true.) then
        call restart_initialize(width, height, u_inflow, u, v, u_old, v_old, p, file_name_save_input, file_number_save_input) 
    endif


    !-----------------------------
    ! Start SMAC...
    !-----------------------------
    do iter = 1, frames
        write(*, *) "iteration = ", iter
        time = time + dt

        !-------------------------------------
        ! u_new, v_newにおいて境界での値を設定
        !-------------------------------------
        call set_upper_boundary_slip_for_vnew(width, height, v_new)
        call set_under_boundary_slip_vnew(width, height, v_new)
        call set_inflow_boundary_for_unew(width, height, u_inflow, u_new)
        call set_outflow_boundary(width, height, dt, h, u, v, u_new, v_new)


        ! TODO: 実装が形状依存なのでこれを直す。
        !-------------------------------
        !　u_tilde, v_tildeを計算
        !-------------------------------
        
        do j = aabb_under_j - 3, aabb_above_j + 3
        do i = aabb_left_i - 3, aabb_right_i + 3
            !x
            call discrete_advection_x(i, j, u, v, h, h, width, height, advection_term_x)
            call discrete_advection_x(i, j, u_old, v_old, h, h, width, height, advection_term_x_old)
            call discrete_viscosity_x(i, j, u, v, h, h, Re, width, height, viscosity_term_x)
            call discrete_pressure_x(i, j, p, width, height, h, pressure_term_x)

            !y
            call discrete_advection_y(i, j, u, v, h, h, width, height, advection_term_y)
            call discrete_advection_y(i, j, u_old, v_old, h, h, width, height, advection_term_y_old)
            call discrete_viscosity_y(i, j, u, v, h, h, Re, width, height, viscosity_term_y)
            call discrete_pressure_y(i, j, p, width, height, h, pressure_term_y)

            u_tilde(i, j) = u(i, j) + dt &
                          * (-pressure_term_x - 0.5d0 * (3.0d0 * advection_term_x - advection_term_x_old) + viscosity_term_x)
            v_tilde(i, j) = v(i, j) + dt &
                          * (-pressure_term_y - 0.5d0 * (3.0d0 * advection_term_y - advection_term_y_old) + viscosity_term_y)
        enddo
        enddo


        ! do j = j_under-3, j_under+3
        ! do i = i_left-3, i_right+3
        !     !x
        !     call discrete_advection_x(i, j, u, v, h, h, width, height, advection_term_x)
        !     call discrete_advection_x(i, j, u_old, v_old, h, h, width, height, advection_term_x_old)
        !     call discrete_viscosity_x(i, j, u, v, h, h, Re, width, height, viscosity_term_x)
        !     call discrete_pressure_x(i, j, p, width, height, h, pressure_term_x)

        !     !y
        !     call discrete_advection_y(i, j, u, v, h, h, width, height, advection_term_y)
        !     call discrete_advection_y(i, j, u_old, v_old, h, h, width, height, advection_term_y_old)
        !     call discrete_viscosity_y(i, j, u, v, h, h, Re, width, height, viscosity_term_y)
        !     call discrete_pressure_y(i, j, p, width, height, h, pressure_term_y)

        !     u_tilde(i, j) = u(i, j) + dt &
        !                   * (-pressure_term_x - 0.5d0 * (3.0d0 * advection_term_x - advection_term_x_old) + viscosity_term_x)
        !     v_tilde(i, j) = v(i, j) + dt &
        !                   * (-pressure_term_y - 0.5d0 * (3.0d0 * advection_term_y - advection_term_y_old) + viscosity_term_y)
        ! enddo
        ! enddo

        ! do j = j_under+4, j_above-4
        ! do i = i_left-3, i_left+3
        !     !x
        !     call discrete_advection_x(i, j, u, v, h, h, width, height, advection_term_x)
        !     call discrete_advection_x(i, j, u_old, v_old, h, h, width, height, advection_term_x_old)
        !     call discrete_viscosity_x(i, j, u, v, h, h, Re, width, height, viscosity_term_x)
        !     call discrete_pressure_x(i, j, p, width, height, h, pressure_term_x)

        !     !y
        !     call discrete_advection_y(i, j, u, v, h, h, width, height, advection_term_y)
        !     call discrete_advection_y(i, j, u_old, v_old, h, h, width, height, advection_term_y_old)
        !     call discrete_viscosity_y(i, j, u, v, h, h, Re, width, height, viscosity_term_y)
        !     call discrete_pressure_y(i, j, p, width, height, h, pressure_term_y)

        !     u_tilde(i, j) = u(i, j) + dt &
        !                   * (-pressure_term_x - 0.5d0 * (3.0d0 * advection_term_x - advection_term_x_old) + viscosity_term_x)
        !     v_tilde(i, j) = v(i, j) + dt &
        !                   * (-pressure_term_y - 0.5d0 * (3.0d0 * advection_term_y - advection_term_y_old) + viscosity_term_y)
        ! enddo
        ! enddo

        ! do j = j_under+4, j_above-4
        ! do i = i_right-3, i_right+3
        !     !x
        !     call discrete_advection_x(i, j, u, v, h, h, width, height, advection_term_x)
        !     call discrete_advection_x(i, j, u_old, v_old, h, h, width, height, advection_term_x_old)
        !     call discrete_viscosity_x(i, j, u, v, h, h, Re, width, height, viscosity_term_x)
        !     call discrete_pressure_x(i, j, p, width, height, h, pressure_term_x)

        !     !y
        !     call discrete_advection_y(i, j, u, v, h, h, width, height, advection_term_y)
        !     call discrete_advection_y(i, j, u_old, v_old, h, h, width, height, advection_term_y_old)
        !     call discrete_viscosity_y(i, j, u, v, h, h, Re, width, height, viscosity_term_y)
        !     call discrete_pressure_y(i, j, p, width, height, h, pressure_term_y)

        !     u_tilde(i, j) = u(i, j) + dt &
        !                   * (-pressure_term_x - 0.5d0 * (3.0d0 * advection_term_x - advection_term_x_old) + viscosity_term_x)
        !     v_tilde(i, j) = v(i, j) + dt &
        !                   * (-pressure_term_y - 0.5d0 * (3.0d0 * advection_term_y - advection_term_y_old) + viscosity_term_y)
        ! enddo
        ! enddo

        ! do j = j_above-3, j_above+3
        ! do i = i_left-3, i_right+3
        !     !x
        !     call discrete_advection_x(i, j, u, v, h, h, width, height, advection_term_x)
        !     call discrete_advection_x(i, j, u_old, v_old, h, h, width, height, advection_term_x_old)
        !     call discrete_viscosity_x(i, j, u, v, h, h, Re, width, height, viscosity_term_x)
        !     call discrete_pressure_x(i, j, p, width, height, h, pressure_term_x)

        !     !y
        !     call discrete_advection_y(i, j, u, v, h, h, width, height, advection_term_y)
        !     call discrete_advection_y(i, j, u_old, v_old, h, h, width, height, advection_term_y_old)
        !     call discrete_viscosity_y(i, j, u, v, h, h, Re, width, height, viscosity_term_y)
        !     call discrete_pressure_y(i, j, p, width, height, h, pressure_term_y)

        !     u_tilde(i, j) = u(i, j) + dt &
        !                   * (-pressure_term_x - 0.5d0 * (3.0d0 * advection_term_x - advection_term_x_old) + viscosity_term_x)
        !     v_tilde(i, j) = v(i, j) + dt &
        !                   * (-pressure_term_y - 0.5d0 * (3.0d0 * advection_term_y - advection_term_y_old) + viscosity_term_y)
        ! enddo
        ! enddo


        !-----------------------------
        ! 角柱表面上のu_tilde, v_tilde
        !-----------------------------
        do k = 1, force_num
            force_index_x = int(x_f(k)/h) + 1
            force_index_y = int(y_f(k)/h) + 1

            u_tilde_in_surface(k) = 0.0d0
            v_tilde_in_surface(k) = 0.0d0

            do j = force_index_y-3, force_index_y+3
            do i = force_index_x-3, force_index_x+3
                call calc_u_point(i, j, h, u_point_x, u_point_y)
                call calc_v_point(i, j, h, v_point_x, v_point_y)

                call delta_h_1d(u_point_x, x_f(k), h, delta_h_x)
                call delta_h_1d(u_point_y, y_f(k), h, delta_h_y)
                u_tilde_in_surface(k) = u_tilde_in_surface(k) + u_tilde(i, j) * delta_h_x * delta_h_y * h * h

                call delta_h_1d(v_point_x, x_f(k), h, delta_h_x)
                call delta_h_1d(v_point_y, y_f(k), h, delta_h_y)
                v_tilde_in_surface(k) = v_tilde_in_surface(k) + v_tilde(i, j) * delta_h_x * delta_h_y * h * h
            enddo
            enddo

        enddo


        !-------------------------
        ! 角柱表面上の体積力
        !-------------------------

        fx_in_surface = -1.0d0 * u_tilde_in_surface / dt
        fy_in_surface = -1.0d0 * v_tilde_in_surface / dt


        !------------------------------
        ! 各点における体積力を求める
        !------------------------------
        fx = 0.0d0
        fy = 0.0d0
        do k = 1, force_num
            force_index_x = int(x_f(k)/h) + 1
            force_index_y = int(y_f(k)/h) + 1
            do j = force_index_y-3, force_index_y+3
            do i = force_index_x-3, force_index_x+3
                call calc_u_point(i, j, h, u_point_x, u_point_y)
                call calc_v_point(i, j, h, v_point_x, v_point_y)

                call delta_h_1d(u_point_x, x_f(k), h, delta_h_x)
                call delta_h_1d(u_point_y, y_f(k), h, delta_h_y)
                fx(i, j) = fx(i, j) + fx_in_surface(k) * delta_h_x * delta_h_y * h * h

                call delta_h_1d(v_point_x, x_f(k), h, delta_h_x)
                call delta_h_1d(v_point_y, y_f(k), h, delta_h_y)
                fy(i, j) = fy(i, j) + fy_in_surface(k) * delta_h_x * delta_h_y * h * h
            enddo
            enddo
        enddo


        !-----------------------------
        ! クランク・ニコルソン(x)
        !-----------------------------
        !右辺項の設定(x)
        do j = 1, height
        do i = 1, width-1
            call discrete_advection_x(i, j, u, v, h, h, width, height, advection_term_x)
            call discrete_advection_x(i, j, u_old, v_old, h, h, width, height, advection_term_x_old)
            call discrete_viscosity_x(i, j, u, v, h, h, Re, width, height, viscosity_term_x)
            call discrete_pressure_x(i, j, p, width, height, h, pressure_term_x)

            b_crank_x(i, j) = u(i, j) + dt &
                            * (0.5d0 * viscosity_term_x - pressure_term_x &
                            - 0.5d0*(3.0d0 * advection_term_x - advection_term_x_old) + fx(i, j))

            if(i == 1) then
                b_crank_x(i, j) = b_crank_x(i, j) + dt / (2.0d0 * Re) * 1.0d0/(h*h) * u_new(i-1, j)
            endif

            if(i == width-1) then
                b_crank_x(i, j) = b_crank_x(i, j) + dt / (2.0d0 * Re) * 1.0d0/(h*h) * u_new(i+1, j)
            endif

        enddo
        enddo

        ! solve
        call solve_bicgstab(u_sol(0:width, 0:height+1), &
                               a1(0:width, 0:height+1), &
                               a2(0:width, 0:height+1), &
                               a3(0:width, 0:height+1), &
                               a4(0:width, 0:height+1), &
                               a5(0:width, 0:height+1), &
                        b_crank_x(0:width, 0:height+1), &
                                eps, width-1, height, &
                                r(0:width, 0:height+1), &
                               r0(0:width, 0:height+1), &
                            stab1(0:width, 0:height+1), &
                            stab2(0:width, 0:height+1), &
                            stab3(0:width, 0:height+1), &
                            stab4(0:width, 0:height+1))

        u_new(1:width-1, 1:height) = u_sol(1:width-1, 1:height)


        !-----------------------------
        ! クランク・ニコルソン(y)
        !-----------------------------
        do j = 1, height-1
        do i = 1, width
            call discrete_advection_y(i, j, u, v, h, h, width, height, advection_term_y)
            call discrete_advection_y(i, j, u_old, v_old, h, h, width, height, advection_term_y_old)
            call discrete_viscosity_y(i, j, u, v, h, h, Re, width, height, viscosity_term_y)
            call discrete_pressure_y(i, j, p, width, height, h, pressure_term_y)

            b_crank_y(i, j) = v(i, j) + dt &
                            * (0.5d0 * viscosity_term_y - pressure_term_y &
                            - 0.5d0 * (3.0d0*advection_term_y - advection_term_y_old) + fy(i, j))

            if(j == 1) then
                b_crank_y(i, j) = b_crank_y(i, j) + dt / (2.0d0 * Re) * 1.0d0/(h*h) * v_new(i, j-1)
            endif
            
            if(j == height - 1) then
                b_crank_y(i, j) = b_crank_y(i, j) + dt / (2.0d0 * Re) * 1.0d0/(h*h) * v_new(i, j+1)
            endif
        enddo
        enddo
        

        call solve_bicgstab(v_sol(0:width+1, 0:height), &
                               c1(0:width+1, 0:height), &
                               c2(0:width+1, 0:height), &
                               c3(0:width+1, 0:height), &
                               c4(0:width+1, 0:height), &
                               c5(0:width+1, 0:height), &
                        b_crank_y(0:width+1, 0:height), &
                               eps, width, height-1, &
                               r(0:width+1, 0:height), &
                              r0(0:width+1, 0:height), &
                            stab1(0:width+1, 0:height), &
                            stab2(0:width+1, 0:height), &
                            stab3(0:width+1, 0:height), &
                            stab4(0:width+1, 0:height) )

        v_new(1:width, 1:height-1) = v_sol(1:width, 1:height-1)



        !-----------------------------
        ! ポアソン方程式を解く
        !-----------------------------
        do j = 1, height
        do i = 1, width
            b_poisson(i, j) = 1.0d0/(dt * h) * (-u_new(i-1, j) + u_new(i, j) - v_new(i, j-1) + v_new(i, j))
        enddo
        enddo


        call solve_bicgstab(phi, e1, e2, e3, e4, e5, b_poisson, eps, width, height, &
                            r, r0, stab1, stab2, stab3, stab4)

        call set_phi_boundary(width, height, phi)


        !-------------------------------
        ! 速度場、圧力場の補正
        !-------------------------------
        do j = 1, height
        do i = 1, width - 1
            u_new(i, j) = u_new(i, j) - dt * (-phi(i, j) + phi(i+1, j)) / h
        enddo
        enddo

        do j = 1, height-1
        do i = 1, width
            v_new(i, j) = v_new(i, j) - dt * (-phi(i, j) + phi(i, j+1)) / h
        enddo
        enddo

        do j = 1, height
        do i = 1, width
            call discrete_viscosity_x(i, j, phi, phi, h, h, Re, width, height, modify_phi_term)

            p(i, j) = p(i, j) + phi(i, j) - 0.5d0 * dt * modify_phi_term
        enddo
        enddo


        call set_inflow_boundary(width, height, u_inflow, u_new, v_new)
        call set_under_boundary_slip(width, height, u_new, v_new)
        call set_upper_boundary_slip(width, height, u_new, v_new)

        u_old = u
        v_old = v
        u = u_new
        v = v_new
        u_new = 0.0d0
        v_new = 0.0d0

        if (iter > frames-1000) then
            call log_velocity
        endif

        if (iter == frames) then
            call savedata(width, height, u_inflow, u, v, u_old, v_old, p, file_number_save_output, file_name_save_output)
        endif

        call log_velocity_point(5.0d0, 2.0d0)

    enddo


    call close_file

end program main
