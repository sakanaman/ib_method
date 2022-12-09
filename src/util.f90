subroutine calc_u_point(i, j, h, res_x, res_y)
    implicit none
    integer,intent(in) :: i, j
    real(8), intent(in) :: h
    real(8), intent(out) :: res_x, res_y
    
    res_x =              dble(i) * h
    res_y = -0.5d0 * h + dble(j) * h 
end subroutine calc_u_point


subroutine calc_v_point(i, j, h, res_x, res_y)
    implicit none
    integer,intent(in) :: i, j
    real(8), intent(in) :: h
    real(8), intent(out) :: res_x, res_y

    res_x = -0.5d0 * h + dble(i) * h
    res_y =              dble(j) * h
end subroutine calc_v_point