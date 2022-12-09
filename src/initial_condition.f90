!==========================注意===========================
! ここでは、IBMの初期化はしない(ibm.f90を参照)
! また、パラメーター設定や配列のアローケーションは主プログラムで行う
! そして、いかなる変数も宣言時にゼロ初期化されている前提で進める。
!========================================================

subroutine intial_condition(width, height, u_inflow, u, v, u_old, v_old, p)
    implicit none

    integer, intent(in) :: width, height
    real(8), intent(in) :: u_inflow

    real(8), intent(out) :: u(0: width+1, 0:height+1)
    real(8), intent(out) :: v(0: width+1, 0:height+1)
    real(8), intent(out) :: u_old(0: width+1, 0:height+1) ! for アダムス・バッシュフォース
    real(8), intent(out) :: v_old(0: width+1, 0:height+1) ! for アダムス・バッシュフォース(以下、AB法と呼ぶ)
    real(8), intent(out) :: p(0: width+1, 0:height+1)

    integer :: i, j

    !----------------------
    ! 内部の主流方向uの初期化
    !----------------------
    do j = 1, height
    do i = 1, width-1
        u(i, j) = u_inflow
    enddo
    enddo

    !------------------------
    ! 流出境界の主流方向uの初期化
    !------------------------
    do j = 1, height
        u(width, j) = u_inflow
    enddo

    !---------------------------------
    ! 内部、境界の主流垂直方向速度vの初期化
    ! (0なので必要ないが、明示的に)
    !---------------------------------
    v = 0.0d0


    !-----------------------------------------------------
    ! u_old, v_oldの初期化
    ! (u,vと同じ値にすれば１回目のAB法はオイラー陽解法になる..はず)
    !-----------------------------------------------------
    u_old = u
    v_old = v

    !----------------------
    ! 圧力の初期化(0スタート)
    ! (明示的にここでも初期化する)
    !----------------------
    p = 0.0d0

end subroutine intial_condition


subroutine restart_initialize(width, height, u_inflow, u, v, u_old, v_old, p, file_name, file_number)
     
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

    open(file_number, file = file_name, status = "old")

    ! save fileの仕様: U -> V -> U_old -> V_old -> P の順番で入ってること

    read(file_number, *) u
    read(file_number, *) v
    read(file_number, *) u_old
    read(file_number, *) v_old
    read(file_number, *) p
    
    close(file_number) 
end subroutine restart_initialize