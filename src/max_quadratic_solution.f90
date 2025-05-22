! Esta rutina calcula el tiempo en que deven salir las partículas para evitar que salgan muy pegadas
  
! Subrutina para calcular la solución mayor
subroutine max_quadratic_solution(a, b, c, max_root, has_solution)
  implicit none
  double precision, intent(in) :: a, b, c
  double precision, intent(out) :: max_root
  logical, intent(out) :: has_solution
  double precision :: discriminant, root1, root2
  
  ! Calcular el discriminante
  discriminant = b**2 - 4.0d0 * a * c
  
  ! Verificar si las soluciones son reales
  if (discriminant < 0.0d0) then
     has_solution = .false.
     max_root = 0.0d0
     return
  endif
  
  ! Calcular las dos raíces
  root1 = (-b + sqrt(discriminant)) / (2.0d0 * a)
  root2 = (-b - sqrt(discriminant)) / (2.0d0 * a)
  
  ! Devolver la raíz mayor
  max_root = max(root1, root2)
  has_solution = .true.
  
end subroutine max_quadratic_solution


