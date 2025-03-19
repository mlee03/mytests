program this_program

  implicit none
  
  integer, allocatable :: this_allocatable
  integer, pointer :: this_pointer => NULL()
  
  call this_subroutine(this_allocatable)
  call this_subroutine(this_pointer)

contains

  subroutine this_subroutine(this_argument)
    implicit none
    integer, intent(in), optional :: this_argument
    write(*,*) "ITS PRESENCE", present(this_argument)
  end subroutine this_subroutine
  
end program this_program
  
