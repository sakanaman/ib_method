!================================
! 角柱一個(カルマン渦の実験)
!================================
subroutine set_force_points_in_prism(x_left, x_right, y_under, y_above, Nx, Ny, x_f, y_f)
    implicit none
    real(8), intent(in) :: x_left, x_right  ! 角柱左側x座標, 右側x座標
    real(8), intent(in) :: y_under, y_above  ! 角柱下側y座標, 上側y座標
    integer, intent(in) :: Nx, Ny ! 各軸方向ごとの体積力点を置く個数
    real(8), intent(out) :: x_f(1 : Nx*2 + Ny*2 - 4)
    real(8), intent(out) :: y_f(1 : Nx*2 + Ny*2 - 4)

    real(8) :: ds_x, ds_y
    integer :: i, j, offset

    !各辺重複込みでNx, Ny個の点を置く前提なので植木算的に間隔はNx-1個, Ny-1個になる
    ds_x = (x_right - x_left) / dble(Nx - 1)
    ds_y = (y_above - y_under) / dble(Ny - 1)

    

    !------------------
    ! まずは下側に
    !------------------
    do i = 1, Nx-1
        x_f(i) = x_left + ds_x * i
        y_f(i) = y_under
    enddo

    !------------------
    ! 次は右側
    !------------------
    offset = Nx - 1
    do j = 1, Ny-1
        x_f(j + offset) = x_right
        y_f(j + offset) = y_under + ds_y * j
    enddo

    !-----------------
    ! 上側
    !-----------------
    offset = Nx + Ny - 2
    do i = 1, Nx-1
        x_f(i + offset) = x_right - ds_x * i
        y_f(i + offset) = y_above
    enddo

    !------------------
    ! 左側
    !------------------
    offset = 2 * Nx + Ny - 3
    do j = 1, Ny - 1
        x_f(j + offset) = x_left
        y_f(j + offset) = y_above - ds_y * j
    enddo

end subroutine set_force_points_in_prism


subroutine set_force_points_in_circle(center_x, center_y, radius, N_circle, x_f, y_f, &
                                      x_min, x_max, y_min, y_max)
    implicit none
    real(8), intent(in) :: center_x, center_y
    real(8), intent(in) :: radius
    integer, intent(in) :: N_circle
    real(8), intent(out) :: x_f(1 : N_circle), y_f(1 : N_circle)
    real(8), intent(out) :: x_min, x_max, y_min, y_max

    integer :: i, j
    real(8), parameter :: pi = 3.141592653589793238462643383279502884d0
    real(8) :: delta_theta

    delta_theta = 2.0d0 * pi / dble(N_circle)

    x_min = 100000.0d0
    x_max = -100000.0d0
    y_min = 100000.0d0
    y_max = -100000.0d0


    do i = 1, N_circle
        x_f(i) = center_x + radius * cos(delta_theta * dble(i))
        y_f(i) = center_y + radius * sin(delta_theta * dble(i))

        if (x_f(i) < x_min) then
            x_min = x_f(i)
        endif

        if (x_f(i) > x_max) then
            x_max = x_f(i)
        endif

        if (y_f(i) < y_min) then
            y_min = y_f(i)
        endif

        if (y_f(i) > y_max) then
            y_max = y_f(i)
        endif
    enddo
end subroutine set_force_points_in_circle

subroutine delta_h_1d(x,x0,h,res)
    implicit none
    real(8), intent(in) :: x, x0, h
    real(8), intent(out) :: res

    real(8) :: r


    r = abs(x - x0)/h

    if(0.5 <= r .and. r <= 1.5) then
        res = 1.0d0/h &
            * 1.0d0/6.0d0 * (5.0d0 - 3.0d0 * r - sqrt(-3.0d0 * (1.0d0 - r) ** 2.0d0 + 1.0d0))
    else if (r < 0.5) then
        res = 1.0d0/h &
            * 1.0d0/3.0d0 * (1.0d0 + sqrt(-3.0d0 * r ** 2.0d0 + 1.0d0))
    else
        res = 0.0d0
    endif
end subroutine delta_h_1d
