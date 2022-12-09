!==============================注意============================
! 下のサブルーチンの使い方を間違えないようにこの注意書きを見ること！
!
! 時間更新の手順は、
! u_old = u(-1), u = u(0), u_new = ?  <--- 初期条件
!         ↓
! u_old = u(-1), u = u(0), u_new = u(1)
!         ↓
! u_old = u(0),  u = u(1), u_new = u(2)
! といった感じで進んでいくことを想定(vも同様)している
! また、unew,vnewは予測速度up, vpの役割も兼ねる
!=============================================================


!上側境界(すべり境界)の設定(楕円型偏微分方程式の係数調整はここではしないので注意)
! u_old, u, v_old, vに対して適用
subroutine set_upper_boundary_slip(width ,height, u, v)
    implicit none
    integer, intent(in) :: width, height
    real(8), intent(out) :: u(0:width+1, 0:height+1)
    real(8), intent(out) :: v(0:width+1, 0:height+1)

    integer :: i

    !-----------------------
    ! 垂直方向の速度勾配 = 0
    !-----------------------
    do i = 1, width-1
        u(i, height+1) = u(i, height)
    enddo
end subroutine set_upper_boundary_slip


!上側境界(すべり境界)の設定(楕円型偏微分方程式の係数調整はここではしないので注意)
! v_newに対して適用
subroutine set_upper_boundary_slip_for_vnew(width, height, v_new)
    implicit none
    integer, intent(in) :: width, height
    real(8), intent(out) :: v_new(0:width+1, 0:height+1)

    integer :: i, j

    !-----------------------
    !  非透過条件
    !-----------------------
    do i = 1, width
        v_new(i, height) = 0.0d0
    enddo
end subroutine set_upper_boundary_slip_for_vnew



!下側境界(すべり境界)の設定(楕円型偏微分方程式の係数調整はここではしないので注意)
!u_old, u, v_old, vに対して適用
subroutine set_under_boundary_slip(width, height, u, v)
    implicit none
    integer, intent(in) :: width, height
    real(8), intent(inout) :: u(0:width+1, 0:height+1)
    real(8), intent(inout) :: v(0:width+1, 0:height+1)
    integer :: i


    !----------------------
    ! 垂直方向の速度勾配 = 0
    !----------------------
    do i = 1, width-1
        u(i, 0) = u(i, 1)
    enddo

end subroutine set_under_boundary_slip


!下側境界(すべり境界)の設定(楕円型偏微分方程式の係数調整はここではしないので注意)
!v_newに対して適用
subroutine set_under_boundary_slip_vnew(width, height, v_new)
    implicit none
    integer, intent(in) :: width, height
    real(8), intent(out) :: v_new(0: width+1, 0: height+1)

    integer :: i
    !----------------------
    ! 非透過条件
    !----------------------
    do i = 1, width
        v_new(i, 0) = 0.0d0
    enddo
end subroutine set_under_boundary_slip_vnew



!流入境界(u_old, uに対して適用)
subroutine set_inflow_boundary(width, height, u_inflow, u, v)
    implicit none
    integer, intent(in) :: width, height
    real(8), intent(in) :: u_inflow !流入速度
    real(8), intent(out) :: u(0:width+1, 0:height+1)
    real(8), intent(out) :: v(0:width+1, 0:height+1)

    integer :: j

    !---------------
    ! 流入速度vの設定
    !---------------
    do j = 1, height-1
        v(0, j) = -v(1, j)
    enddo
end subroutine set_inflow_boundary


!u_newに対して適用
subroutine set_inflow_boundary_for_unew(width, height, u_inflow, u_new)
    implicit none
    integer, intent(in) :: width, height
    real(8), intent(in) :: u_inflow !流入速度
    real(8), intent(out) :: u_new(0:width+1, 0:height+1)

    integer :: j

    !---------------
    ! 流入速度uの設定
    !---------------
    do j = 1, height
        u_new(0, j) = u_inflow
    enddo
end subroutine set_inflow_boundary_for_unew


!流出境界(u -> u_newの流出境界を更新, uの情報は全て既知となっているように実装すること)
subroutine set_outflow_boundary(width, height, dt, dx, u, v, u_new, v_new)
    implicit none
    integer, intent(in) :: width, height
    real(8), intent(in) :: dt, dx
    real(8), intent(in) :: u(0:width+1, 0:height+1)
    real(8), intent(in) :: v(0:width+1, 0:height+1)
    real(8), intent(out) :: u_new(0:width+1, 0:height+1)
    real(8), intent(out) :: v_new(0:width+1, 0:height+1)

    real(8) :: um
    integer :: j

    !-------------------------
    ! 主流方向の流出平均速度の計算
    !-------------------------
    um = 0.0d0
    do j = 1, height
        um = um + u(width, j)
    enddo
    um = um / dble(height)


    !-----------------------------
    ! 流出速度uの計算(with 移流方程式)
    !-----------------------------
    do j = 1, height
        u_new(width, j) =  u(width, j)    - dt * um * (u(width, j)   - u(width-1, j)) / dx
    enddo

    !-----------------------------
    ! 流出速度vの計算(with 移流方程式)
    !-----------------------------
    do j = 1, height-1
        v_new(width+1, j) = v(width+1, j) - dt * um * (v(width+1, j) - v(width, j))   / dx
    enddo
end subroutine set_outflow_boundary


subroutine set_phi_boundary(width, height, phi)
    implicit none
    integer, intent(in) :: width, height
    real(8), intent(out) :: phi(0:width+1, 0:height+1)

    integer :: i, j

    do i = 1, width
        phi(i, 0) = phi(i, 1)
        phi(i, height+1) = phi(i, height)
    enddo

    do j = 1, height
        phi(0, j) = phi(1, j)
        phi(width+1, j) = phi(width, j)
    enddo
end subroutine set_phi_boundary

