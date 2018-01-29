module shift
use variables
IMPLICIT none

contains

  subroutine shift_format(iout,fname,library_name)
    integer,intent(in) :: iout
    character(*),intent(in) :: fname, library_name
    integer :: ios
    integer :: counter ! idk what this is but it seems to increment on each line

    character(*),parameter :: isDefault = 'isDefault="false"'
    character(*),parameter :: isUsed = 'isUsed="true"'
    character(*),parameter :: nline = '<Parameter docString="" '

    101 format(a) ! plain-text format

    open(unit=iout,file=fname,status='replace',action='write',iostat=ios)
    if (ios .ne. 0) then
      write(*,'(a,i6,2a)') 'FATAL -- error opening unit -- ', iout, ' -- ', fname
      write(*,'(a,i6)') 'ios = ', ios
    endif

    counter = 0

    write(iout,'(3a)') '<ParameterList name="', library_name, '">'
    write(iout,*) nline, 'id = ', counter, '" ', isDefault, ' ', isUsed, ' ',&
      'name="num groups" type="int" value="', ngroup, '"/>'
    counter = counter + 1

    
  endsubroutine shift_format
endmodule shift
