subroutine discrete_pressure_x(i, j, p , width, height, dx, dpdx)

    implicit none

    integer, intent(in) :: i, j, width, height
    real(8), intent(in) :: dx
    real(8), intent(in) :: p(0:width+1, 0:height+1) !形状明示配列
    real(8), intent(out) ::  dpdx

    dpdx = (-p(i, j) + p(i + 1, j))/dx

end subroutine discrete_pressure_x


subroutine discrete_pressure_y(i, j, p , width, height, dy, dpdy)

    implicit none

    integer, intent(in) :: i, j, width, height
    real(8), intent(in) :: dy
    real(8), intent(in) :: p(0 : width+1, 0 : height+1) !形状明示配列
    real(8), intent(out) :: dpdy

    dpdy = (-p(i, j) + p(i, j + 1))/dy

end subroutine discrete_pressure_y


!対流項のx方向への離散化
subroutine discrete_advection_x(i, j, u, v, dx, dy, width, height, ad_x)
    
    implicit none
    
    integer, intent(in) :: i, j, width, height
    real(8), intent(in) :: dx, dy
    real(8), intent(in) :: u(0 : width+1, 0 : height+1), v(0 : width+1, 0 : height+1)
    real(8), intent(out) :: ad_x

    real(8) :: tmp1, tmp2, tmp3

    tmp1 = 1.0d0/dx * ( ((u(i, j) + u(i+1, j)) * 0.5d0) ** 2.0d0 - ((u(i-1, j) + u(i, j)) * 0.5d0) ** 2.0d0)

    tmp2 = 1.0d0/dy * ( (u(i, j) + u(i, j+1)) * 0.5d0 * (v(i, j) + v(i+1, j)) * 0.5d0)
    tmp3 = 1.0d0/dy * (- (u(i, j-1) + u(i, j)) * 0.5d0 * (v(i, j-1) + v(i+1, j-1)) * 0.5d0)

    ad_x = tmp1 + tmp2 + tmp3
    
    
end subroutine discrete_advection_x


!対流項のy方向への離散化
subroutine discrete_advection_y(i, j, u, v, dx, dy, width, height, ad_y)
    
    implicit none
    
    integer, intent(in) :: i, j, width, height
    real(8), intent(in) :: dx, dy
    real(8), intent(in) :: u(0 : width+1, 0 : height+1), v(0 : width+1, 0 : height+1)
    real(8), intent(out) :: ad_y

    real(8) :: tmp1, tmp2, tmp3, tmp4

    tmp1 = 1.0d0/dx *  ((u(i, j) + u(i, j+1)) * 0.5d0 * (v(i, j) + v(i+1, j)) * 0.5d0)
    tmp2 = 1.0d0/dx *  (- (u(i-1, j) + u(i-1, j+1)) * 0.5d0 * (v(i-1, j) + v(i, j)) * 0.5d0)

    tmp3 = 1.0d0/dy *  (((v(i, j) + v(i, j+1)) * 0.5d0) ** 2.0d0 - ((v(i, j-1) + v(i, j)) * 0.5d0) ** 2.0d0)


    ad_y = tmp1 + tmp2 + tmp3
    
end subroutine discrete_advection_y



!粘性項のx方向への離散化
subroutine discrete_viscosity_x(i, j, u, v, dx, dy, Re, width, height, visc_x)
    
    implicit none
    integer, intent(in) :: i, j, width, height
    real(8), intent(in) :: u(0 : width+1, 0 : height+1) ,v(0 : width+1, 0 : height+1)
    real(8), intent(in) :: dx, dy, Re
    real(8), intent(out) :: visc_x

    real(8) :: tmp1, tmp2

    tmp1 =  1.0d0/(dx * dx) * (u(i-1, j) - 2.0d0 * u(i, j) + u(i+1, j))
    tmp2 =  1.0d0/(dy * dy) * (u(i, j-1) - 2.0d0 * u(i, j) + u(i, j+1))

    visc_x = 1.0d0/Re * (tmp1 + tmp2)

    
end subroutine discrete_viscosity_x



!粘性項のy方向への離散化
subroutine discrete_viscosity_y(i, j, u, v, dx, dy, Re, width, height, visc_y)
    
    implicit none
    
    integer, intent(in) :: i, j, width, height
    real(8), intent(in) :: u(0 : width+1, 0 : height+1), v(0 : width+1, 0 : height+1)
    real(8), intent(in) :: dx, dy, Re
    real(8), intent(out) :: visc_y

    real(8) :: tmp1, tmp2

    tmp1 = 1.0d0/(dx * dx) * (v(i-1, j) - 2.0d0 * v(i, j) + v(i + 1, j))
    tmp2 = 1.0d0/(dy * dy) * (v(i, j-1) - 2.0d0 * v(i, j) + v(i, j + 1))

    visc_y = 1.0d0/Re * (tmp1 + tmp2)

end subroutine discrete_viscosity_y