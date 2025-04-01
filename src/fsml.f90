module fsml
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, fsml!"
  end subroutine say_hello
end module fsml
