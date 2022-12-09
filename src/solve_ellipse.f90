subroutine solve_sor(x, a1, a2, a3, a4, a5, b ,eps, alpha, width, height)
    implicit none
    integer, intent(in) :: width,height
    real(8), intent(in) :: eps, alpha
    real(8), intent(in) :: a1(0:width+1, 0:height+1)
    real(8), intent(in) :: a2(0:width+1, 0:height+1)
    real(8), intent(in) :: a3(0:width+1, 0:height+1)
    real(8), intent(in) :: a4(0:width+1, 0:height+1)
    real(8), intent(in) :: a5(0:width+1, 0:height+1)
    real(8),intent(in) :: b(0:width+1, 0:height+1)
    real(8), intent(out) :: x(0:width+1, 0:height+1)


    real(8) :: x_dash
    real(8) :: x_new
    real(8) :: error = 0.0d0
    integer :: i,j

    do
        error = 0.0d0
        do j = 1, height
            do i = 1, width
                x_dash = (b(i,j) &
                            - ( a2(i, j) * x(i-1, j) + a3(i, j) * x(i+1, j) + a4(i, j) * x(i, j-1) + a5(i, j) * x(i, j+1) )) &
                        / a1(i, j)
                x_new = (1.0d0 - alpha) * x(i,j) + alpha * x_dash
                if( error < abs(x_dash - x(i,j)) * abs(a1(i,j))) then
                    error = abs(x_dash - x(i,j)) * abs(a1(i,j))
                endif
                x(i, j) = x_new
            end do
        end do
        if(error < eps) then
            exit !収束したのでループを抜ける(L∞ノルム)
        endif
    end do
end subroutine solve_sor


subroutine solve_bicgstab(x, a1, a2, a3, a4, a5, b ,eps, width, height, r, r0, p, e, y, v)
    !$ use omp_lib
    implicit none
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

    real(8) :: C1,C2,C3, alpha, beta, ev, vv, rr
    integer :: i, j


    !$ call omp_set_num_threads(8)

    !------------------------------------
    ! 初期ベクトルの設定(x = 0, r = b - Ax = b, r0 = r, p = r)
    !------------------------------------
    ! x = 0.0d0
    ! r = b
    ! r0 = r
    ! p = r
    
    !$OMP parallel private(i, j)
    !$OMP do
    do j = 0, height+1
    do i = 0, width+1
        x(i, j) = 0.0d0
        r(i, j) = b(i, j)
        r0(i, j) = r(i, j)
        p(i, j) = r(i, j)
    enddo
    enddo
    !$OMP enddo
    !$OMP end parallel


    !------------------------------------
    ! C1 = (r0, r0)
    !------------------------------------
    C1 = 0.0d0
    !$OMP parallel private(i, j)
    !$OMP do reduction(+ : C1)
    do j = 1, height
    do i = 1, width
        C1 = C1 + r(i, j) ** 2.0d0
    enddo
    enddo
    !$OMP enddo
    !$OMP end parallel

    !------------------------------------
    ! 奇跡的に何もしないで収束してるとき
    !------------------------------------
    if(sqrt(C1) < eps) then
        return
    endif

    !------------------------------------
    ! ループ開始
    !------------------------------------
    do
        !------------------------------------
        ! y = A*p
        !------------------------------------
        !$OMP parallel private(i, j)
        !$OMP do
        do j = 1, height
        do i = 1, width
            y(i, j) = a1(i,j)*p(i,j) + a2(i,j)*p(i-1,j) + a3(i,j)*p(i+1,j) &
                                     + a4(i,j)*p(i,j-1) + a5(i,j)*p(i,j+1)
        enddo
        enddo
        !$OMP enddo
        !$OMP end parallel

        !------------------------------------
        ! C2 = (r0, y)
        !------------------------------------
        C2 = 0.0d0
        !$OMP parallel private(i, j)
        !$OMP do reduction(+ : C2)
        do j = 1, height
        do i = 1, width
            C2 = C2 + r0(i, j) * y(i, j)
        enddo
        enddo
        !$OMP enddo
        !$OMP end parallel

        !------------------------------------
        ! alpha = C1/C2
        !------------------------------------
        alpha = C1/C2

        !------------------------------------
        ! e = r_k - alpha * y
        !------------------------------------        
        !e = r - alpha * y

        !$OMP parallel private(i, j)
        !$OMP do
        do j = 1, height
        do i = 1, width
            e(i, j) = r(i, j) - alpha * y(i, j)
        enddo
        enddo
        !$OMP enddo
        !$OMP end parallel        

        !------------------------------------
        ! v = A*e
        !------------------------------------
        !$OMP parallel private(i, j)
        !$OMP do
        do j = 1, height
        do i = 1, width
            v(i, j) = a1(i,j)*e(i,j) + a2(i,j)*e(i-1,j) + a3(i,j)*e(i+1,j)&
                                        + a4(i,j)*e(i,j-1) + a5(i,j)*e(i,j+1)
        enddo
        enddo
        !$OMP enddo
        !$OMP end parallel        

        !------------------------------------
        ! ev = (e, v)
        !------------------------------------
        ev = 0.0d0
        !$OMP parallel private(i, j)
        !$OMP do reduction(+ : ev)
        do j = 1, height
        do i = 1, width
            ev = ev + e(i,j) * v(i, j)
        enddo
        enddo
        !$OMP enddo
        !$OMP end parallel

        !------------------------------------
        ! vv = (v, v)
        !------------------------------------
        vv = 0.0d0
        !$OMP parallel private(i, j)
        !$OMP do reduction(+ : vv)
        do j = 1, height
        do i = 1, width
            vv = vv + v(i,j) ** 2.0d0
        enddo
        enddo
        !$OMP enddo
        !$OMP end parallel


        !------------------------------------
        ! C3 = (e, v)/(v, v)
        !------------------------------------
        C3 = ev/vv



        !------------------------------------
        ! x_{k+1} = x_{k} + alpha_k * p_k + C3*e
        !------------------------------------
        ! x = x + alpha * p + C3 * e
        !$OMP parallel private(i, j)
        !$OMP do
        do j = 1, height
        do i = 1, width
            x(i,j) = x(i, j) + alpha * p(i, j) + C3 * e(i, j)
        enddo
        enddo
        !$OMP enddo
        !$OMP end parallel



        !------------------------------------
        ! r_{k+1} = e - C3 * v
        !------------------------------------
        ! r = e - C3 * v

        !$OMP parallel private(i, j)
        !$OMP do
        do j = 1, height
        do i = 1, width
            r(i, j) = e(i, j) - C3 * v(i, j)
        enddo
        enddo
        !$OMP enddo
        !$OMP end parallel


        !------------------------------------
        ! rr = (r, r)     <-- norm(r)
        !------------------------------------
        rr = 0.0d0
        !$OMP parallel private(i, j)
        !$OMP do reduction(+ : rr)
        do j = 1, height
        do i = 1, width
            rr = rr + r(i,j) ** 2.0d0
        enddo
        enddo
        !$OMP enddo
        !$OMP end parallel


        !------------------------------------
        ! 収束判定
        !------------------------------------
        if(sqrt(rr) < eps) then
            exit
        endif



        !------------------------------------
        ! C1 = (r0, r_{k+1})
        !------------------------------------
        C1 = 0.0d0
        !$OMP parallel private(i, j)
        !$OMP do reduction(+ : C1)
        do j = 1, height
        do i = 1, width
            C1 = C1 + r0(i, j) * r(i, j)
        enddo
        enddo
        !$OMP enddo
        !$OMP end parallel


        !------------------------------------
        ! beta_{k} = C1/(C2 * C3)
        !------------------------------------
        beta = C1/(C2 * C3)


        !------------------------------------
        ! p_{k+1} = r_{k+1} + beta_{k+1}*(p_{k} - C3*y)
        !------------------------------------
        ! p = r + beta * (p - C3*y)
        !$OMP parallel private(i, j)
        !$OMP do
        do j = 1, height
        do i = 1, width
            p(i, j) = r(i, j) + beta * (p(i, j) - C3 * y(i, j))
        enddo
        enddo
        !$OMP enddo
        !$OMP end parallel

    enddo


end subroutine solve_bicgstab


!-------- Gauss-Seidel --------
! a1 * x(i, j) + a2 * x(i-1, j) + a3 * x(i+1, j) + a4 * x(i, j-1) + a5 * x(i, j+1) = b
!   --> [ b - { a2 * x(i-1, j) + a3 * x(i+1, j) + a4 * x(i, j-1) + a5 * x(i, j+1) } ] / a1 で x(i, j)を更新する
!
!
!--------- SOR ---------
! x_dash(i, j) = [ b - { a2 * x(i-1, j) + a3 * x(i+1, j) + a4 * x(i, j-1) + a5 * x(i, j+1) } ] / a1
!  ---> new_x(i,j) = (1 - alpha) * x(i,j) + alpha * x_dash(i,j)
!
!