subroutine initialize_output_file
    use global_variable
    implicit none

    open(file=file_name, unit=file_number, status="replace")

    write(file_number, '(a)', advance="no") "frames: "
    write(file_number, '(i0)') frames    

    write(file_number, '(a)', advance="no") "width: "
    write(file_number, '(i0)') width

    write(file_number, '(a)', advance="no") "height: "
    write(file_number, '(i0)') height

    write(file_number, '(a)', advance="no") "Lx: "
    write(file_number, '(f0.15)') Lx

    write(file_number, '(a)', advance="no") "Ly: "
    write(file_number, '(f0.15)') Ly

    write(file_number, '(a)', advance="no") "prism_left: "
    write(file_number, '(f0.15)') prism_left

    write(file_number, '(a)', advance="no") "prism_right: "
    write(file_number, '(f0.15)') prism_right

    write(file_number, '(a)', advance="no") "prism_under: "
    write(file_number, '(f0.15)') prism_under

    write(file_number, '(a)', advance="no") "prism_above: "
    write(file_number, '(f0.15)') prism_above

    write(file_number, '(a)', advance="no") "dt: "
    write(file_number, '(f0.15)') dt

end subroutine initialize_output_file

subroutine log_velocity
    use global_variable
    implicit none

    integer :: i, j
    
    do j = height, 1, -1
        do i = 1, width-1
            write(file_number, '(f0.15, " ")', advance="no") u(i, j)
        enddo
        write(file_number, '(f0.15, " ")') u(width, j)
    enddo

    do j = height, 1, -1 
        do i = 1, width-1
            write(file_number, '(f0.15, " ")', advance="no") v(i, j)
        enddo
        write(file_number, '(f0.15, " ")') v(width, j)
    enddo
end subroutine log_velocity


subroutine initialize_outputfile_for_point
    use global_variable
    implicit none
    
    open(unit=file_number_point, file=file_name_point, status="replace")

end subroutine initialize_outputfile_for_point

subroutine log_velocity_point(x, y)
    use global_variable
    implicit none
    real(8), intent(in) :: x, y

    integer :: i_x, j_y

    i_x = int(x/h) + 1
    j_y = int(y/h) + 1

    write(file_number_point, '(f0.15)') v(i_x, j_y)
    
end subroutine log_velocity_point

subroutine savedata(width, height, u_inflow, u, v, u_old, v_old, p, file_number, file_name)
    implicit none
    
    integer, intent(in) :: width, height
    real(8), intent(in) :: u_inflow
    character(len=256), intent(in) :: file_name
    integer, intent(in) :: file_number

    real(8), intent(out) :: u(0: width+1, 0:height+1)
    real(8), intent(out) :: v(0: width+1, 0:height+1)
    real(8), intent(out) :: u_old(0: width+1, 0:height+1) ! for アダムス・バッシュフォース
    real(8), intent(out) :: v_old(0: width+1, 0:height+1) ! for アダムス・バッシュフォース(以下、AB法と呼ぶ)
    real(8), intent(out) :: p(0: width+1, 0:height+1)

    integer :: i, j

    open(file_number, file=file_name, status="replace")

    do j = 0, height+1
        do i = 0, width
            write(file_number, '(f0.15, " ")', advance="no") u(i, j)
        enddo
        write(file_number, '(f0.15, " ")') u(width+1, j)
    enddo

    do j = 0, height+1
        do i = 0, width
            write(file_number, '(f0.15, " ")', advance="no") v(i, j)
        enddo
        write(file_number, '(f0.15, " ")') v(width+1, j)
    enddo

    do j = 0, height+1
        do i = 0, width
            write(file_number, '(f0.15, " ")', advance="no") u_old(i, j)
        enddo
        write(file_number, '(f0.15, " ")') u_old(width+1, j)
    enddo

    do j = 0, height+1
        do i = 0, width
            write(file_number, '(f0.15, " ")', advance="no") v_old(i, j)
        enddo
        write(file_number, '(f0.15, " ")') v_old(width+1, j)
    enddo

    do j = 0, height+1
        do i = 0, width
            write(file_number, '(f0.15, " ")', advance="no") p(i, j)
        enddo
        write(file_number, '(f0.15, " ")') p(width+1, j)
    enddo

    close(file_number)
end subroutine savedata


subroutine close_file
    use global_variable
    implicit none
    close(file_number_point)
    close(file_number)
end subroutine close_file