!クランクニコルソン法で陰的に計算する際の連立方程式の左辺項の調整
! (前提: 右辺項の離散化、u_pの流入/流出境界の値は設定済みとする)
subroutine crank_nicolson_left_side_x(width, height, h, Re, dt, &
                                     a1, a2, a3, a4, a5)
    implicit none
    integer, intent(in) :: width, height
    real(8), intent(in) :: h, Re, dt
    real(8), intent(out) :: a1(0:width+1, 0:height+1)
    real(8), intent(out) :: a2(0:width+1, 0:height+1)
    real(8), intent(out) :: a3(0:width+1, 0:height+1)
    real(8), intent(out) :: a4(0:width+1, 0:height+1)
    real(8), intent(out) :: a5(0:width+1, 0:height+1)

    integer :: i, j


    do j = 1, height
    do i = 1, width-1
        a1(i, j) = 1.0d0 + dt / (2.0d0*Re) * 4.0d0/(h*h)
        a2(i, j) = - dt / (2.0d0 * Re) * 1.0d0/(h*h)
        a3(i, j) = - dt / (2.0d0 * Re) * 1.0d0/(h*h)
        a4(i, j) = - dt / (2.0d0 * Re) * 1.0d0/(h*h)
        a5(i, j) = - dt / (2.0d0 * Re) * 1.0d0/(h*h)

        !左側
        if(i == 1) then
            a2(i, j) = 0.0d0
        endif

        !右側
        if(i == width - 1) then
            a3(i, j) = 0.0d0
        endif

        !下側
        if(j == 1) then
            a1(i, j) = a1(i, j) + a4(i, j)
            a4(i, j) = 0.0d0
        endif

        !上側
        if(j == height) then
            a1(i, j) = a1(i, j) + a5(i, j)
            a5(i, j) = 0.0d0
        endif


    enddo
    enddo
end subroutine crank_nicolson_left_side_x


!クランクニコルソン法で陰的に計算する際の連立方程式の左辺項の調整
! (前提: 右辺項の離散化、v_pの流入/流出境界の値は設定済みとする)
subroutine crank_nicolson_left_side_y(width, height, h, Re, dt, &
                                    c1, c2, c3, c4, c5)
    implicit none
    integer, intent(in) :: width, height
    real(8), intent(in) :: h, Re, dt
    real(8), intent(out) :: c1(0:width+1, 0:height+1)
    real(8), intent(out) :: c2(0:width+1, 0:height+1)
    real(8), intent(out) :: c3(0:width+1, 0:height+1)
    real(8), intent(out) :: c4(0:width+1, 0:height+1)
    real(8), intent(out) :: c5(0:width+1, 0:height+1)

    integer :: i, j

    do j = 1, height - 1
    do i = 1, width
        c1(i, j) = 1.0d0 + dt / (2.0d0*Re) * 4.0d0/(h*h)
        c2(i, j) = - dt / (2.0d0 * Re) * 1.0d0/(h*h)
        c3(i, j) = - dt / (2.0d0 * Re) * 1.0d0/(h*h)
        c4(i, j) = - dt / (2.0d0 * Re) * 1.0d0/(h*h)
        c5(i, j) = - dt / (2.0d0 * Re) * 1.0d0/(h*h)

        if(i == 1) then
            c1(i, j) = c1(i, j) + c2(i, j)
            c2(i, j) = 0.0d0
        endif

        if(i == width) then
            c1(i, j) = c1(i, j) + c3(i, j)
            c3(i, j) = 0.0d0
        endif

        if(j == 1) then
            c4(i, j) = 0.0d0
        endif

        if(j == height-1) then
            c5(i, j) = 0.0d0
        endif

    enddo
    enddo
end subroutine crank_nicolson_left_side_y


subroutine poisson_left_side(width, height, h, e1, e2, e3, e4, e5)
    implicit none
    integer, intent(in) :: width, height
    real(8), intent(in) :: h
    real(8), intent(out) :: e1(0:width+1, 0:height+1)
    real(8), intent(out) :: e2(0:width+1, 0:height+1)
    real(8), intent(out) :: e3(0:width+1, 0:height+1)
    real(8), intent(out) :: e4(0:width+1, 0:height+1)
    real(8), intent(out) :: e5(0:width+1, 0:height+1)

    integer :: i, j

    do j = 1, height
    do i = 1, width
        e1(i, j) = -4.0d0/(h*h)
        e2(i, j) = 1.0d0/(h * h)
        e3(i, j) = 1.0d0/(h * h)
        e4(i, j) = 1.0d0/(h * h)
        e5(i, j) = 1.0d0/(h * h)

        if(i == 1) then
            e1(i, j) = e1(i, j) + e2(i, j)
            e2(i, j) = 0.0d0
        endif

        if(i == width) then
            e1(i, j) = e1(i, j) + e3(i, j)
            e3(i, j) = 0.0d0
        endif

        if(j == 1) then
            e1(i, j) = e1(i, j) + e4(i, j)
            e4(i, j) = 0.0d0
        endif

        if(j == height) then
            e1(i, j) = e1(i, j) + e5(i, j)
            e5(i, j) = 0.0d0
        endif
    enddo
    enddo
end subroutine poisson_left_side