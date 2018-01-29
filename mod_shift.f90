module shift
use variables
IMPLICIT none

contains

  subroutine shift_format(iout,fname,library_name)
    integer,intent(in) :: iout
    character(*),intent(in) :: fname, library_name
    integer :: ios
    integer :: counter ! idk what this is but it seems to increment on each line
    character(10) :: tmpChar
    integer :: i,g,gprime

    character(*),parameter :: isDefault = 'isDefault="false"'
    character(*),parameter :: isUsed = 'isUsed="true"'
    character(*),parameter :: nline = '<Parameter docString="" '

    101 format(a) ! plain-text format
    201 format(i10) ! 10-character integer for tmpChar

    open(unit=iout,file=fname,status='replace',action='write',iostat=ios)
    if (ios .ne. 0) then
      write(*,'(a,i6,2a)') 'FATAL -- error opening unit -- ', iout, ' -- ', fname
      write(*,'(a,i6)') 'ios = ', ios
    endif

    counter = 0
    write(tmpChar,201) counter
    tmpChar = trim(adjustl(tmpChar))
    write(iout,'(3a)') '<ParameterList name="', trim(library_name), '">'
    write(iout,'(9a)',advance='no') nline, 'id="', trim(tmpChar), '" ', &
      isDefault, ' ', isUsed
    write(tmpChar,201) ngroup
    tmpChar = trim(adjustl(tmpChar))
    write(iout,*) 'name="num groups" type="int" value="', trim(tmpChar), '"/>'
    
    counter = counter + 1
    write(tmpChar,201) counter
    tmpChar = trim(adjustl(tmpChar))
    write(iout,'(9a)',advance='no') nline, 'id="', trim(tmpChar), '" ', &
      isDefault, ' ', isUsed
    write(tmpChar,201) maxord
    tmpChar = trim(adjustl(tmpChar))
    write(iout,*) 'name="pn order" type="int" value="', trim(tmpChar), '"/>'
    
    counter = counter + 1
    write(tmpChar,201) counter
    tmpChar = trim(adjustl(tmpChar))
    write(iout,'(9a)',advance='no') nline, 'id="', trim(tmpChar), '" ', &
      isDefault, ' ', isUsed
    write(iout,'(a)',advance='no') 'name="material names" type="Array(string)" value="{'
    do i = 1,niso-1
      write(iout,'(2a)',advance='no') trim(adjustl(hisonm(i))), ','
    enddo
    write(iout,'(2a)') trim(adjustl(hisonm(niso))), '}"/>'

    do i = 1,niso
      counter = counter + 1
      write(tmpChar,201) counter
      tmpChar = trim(adjustl(tmpChar))
      write(iout,'(3a)',advance='no') '<ParameterList id="', trim(tmpChar), '"'
      write(iout,*) 'name="', trim(adjustl(hisonm(i))), '">'

      counter = counter + 1
      write(tmpChar,201) counter
      tmpChar = trim(adjustl(tmpChar))
      write(iout,'(9a)',advance='no') nline, 'id="', trim(tmpChar), '" ', &
        isDefault, ' ', isUsed
      write(iout,101,advance='no') 'name="total" type="Array(double)" value="{'
      do g = 1,ngroup-1
        write(iout,'(e12.6,a)',advance='no') xs(i)%p0trans(g), ','
      enddo
      write(iout,'(e12.6,a)') xs(i)%p0trans(ngroup), '}"/>'

      counter = counter + 1
      write(tmpChar,201) counter
      tmpChar = trim(adjustl(tmpChar))
      write(iout,'(9a)',advance='no') nline, 'id="', trim(tmpChar), '" ', &
        isDefault, ' ', isUsed
      write(iout,101,advance='no') 'name="scattering" type="Array(double)" value="{'
      do g = 1,ngroup-1
        write(iout,'(e12.6,a)',advance='no') xs(i)%p0trans(g), ','
      enddo
      write(iout,'(e12.6,a)') xs(i)%p0trans(ngroup), '}"/>'
      
      if (any(sfis(i,:) > 0.0)) then
        ! is fissionable
        counter = counter + 1
        write(tmpChar,201) counter
        tmpChar = trim(adjustl(tmpChar))
        write(iout,'(9a)',advance='no') nline, 'id="', trim(tmpChar), '" ', &
          isDefault, ' ', isUsed
        write(iout,101,advance='no') 'name="fission" type="Array(double)" value="{'
        do g = 1,ngroup-1
          write(iout,'(e12.6,a)',advance='no') xs(i)%sigf(g), ','
        enddo
        write(iout,'(e12.6,a)') xs(i)%sigf(ngroup), '}"/>'

        counter = counter + 1
        write(tmpChar,201) counter
        tmpChar = trim(adjustl(tmpChar))
        write(iout,'(9a)',advance='no') nline, 'id="', trim(tmpChar), '" ', &
          isDefault, ' ', isUsed
        write(iout,101,advance='no') 'name="nu_fission" type="Array(double)" value="{'
        do g = 1,ngroup-1
          write(iout,'(e12.6,a)',advance='no') xs(i)%nuf(g) * xs(i)%sigf(g), ','
        enddo
        write(iout,'(e12.6,a)') xs(i)%nuf(ngroup) * xs(i)%sigf(ngroup), '}"/>'

        counter = counter + 1
        write(tmpChar,201) counter
        tmpChar = trim(adjustl(tmpChar))
        write(iout,'(9a)',advance='no') nline, 'id="', trim(tmpChar), '" ', &
          isDefault, ' ', isUsed
        write(iout,101,advance='no') 'name="chi" type="Array(double)" value="{'
        do g = 1,ngroup-1
          write(iout,'(e12.6,a)',advance='no') xs(i)%chi(g), ','
        enddo
        write(iout,'(e12.6,a)') xs(i)%chi(ngroup), '}"/>'
      endif

      counter = counter + 1
      write(tmpChar,201) counter
      tmpChar = trim(adjustl(tmpChar))
      write(iout,'(9a)',advance='no') nline, 'id="', trim(tmpChar), '" ', &
        isDefault, ' ', isUsed
      write(iout,101,advance='no') 'name="sigma_s0" type="TwoDArray(double)"'
      write(tmpChar,201) ngroup
      tmpChar = trim(adjustl(tmpChar))
      write(iout,'(5a)',advance='no') 'value="', trim(tmpChar), 'x', trim(tmpChar), ':{'
      do g = 1,ngroup
        do gprime = 1,ngroup-1
          write(iout,'(e12.6,a)',advance='no') 1.0d0, ','
        enddo
      enddo
      write(iout,'(e12.6,a)') 1.0d0, '}"/>'

      write(iout,101) '</ParameterList>'
    enddo ! niso

    write(iout,101) '</ParameterList>'

  endsubroutine shift_format
endmodule shift
