module shift
use variables
use spectrum_calc
IMPLICIT none

contains

  subroutine shift_format(iout,fname,library_name)
    integer,intent(in) :: iout
    character(*),intent(in) :: fname, library_name
    integer :: ios
    integer :: counter ! idk what this is but it seems to increment on each line
    character(10) :: tmpChar
    integer :: i,g,gprime

    real(8),dimension(:),allocatable :: shift_total, shift_scattering, &
      shift_fission, shift_nu_fission, shift_chi
    real(8),dimension(:,:),allocatable :: shift_sigma_s0

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
  
    ! do homogenization
    if (.not. allocated(chi_tilde)) then
      call spectrum_solve()
    endif
    allocate(shift_total(ngroup))
    allocate(shift_scattering(ngroup))
    allocate(shift_fission(ngroup))
    allocate(shift_nu_fission(ngroup))
    allocate(shift_chi(ngroup))
    allocate(shift_sigma_s0(ngroup,ngroup))
    shift_chi = chi_tilde
    do i = 1,niso
      shift_total(:) = shift_total(:) + xs(i)%p0trans(:)
      do g = 1,ngroup
        shift_scattering(g) = shift_scattering(g) + sum(xs(i)%mpact_scat(g,:))
      enddo
      shift_fission(:) = shift_fission(:) + xs(i)%sigf(:)
      shift_nu_fission(:) = shift_nu_fission(:) + xs(i)%nuf(:) * xs(i)%sigf(:)
      shift_sigma_s0(:,:) = shift_sigma_s0(:,:) + xs(i)%mpact_scat(:,:)
    enddo

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
    !do i = 1,niso-1
    !  write(iout,'(2a)',advance='no') trim(adjustl(hisonm(i))), ','
    !enddo
    write(iout,'(a)',advance='no') 'homog'
    write(iout,'(2a)') trim(adjustl(hisonm(niso))), '}"/>'

    counter = counter + 1
    write(tmpChar,201) counter
    tmpChar = trim(adjustl(tmpChar))
    write(iout,'(3a)',advance='no') '<ParameterList id="', trim(tmpChar), '"'
    write(iout,*) 'name="', 'homog', '">'

    counter = counter + 1
    write(tmpChar,201) counter
    tmpChar = trim(adjustl(tmpChar))
    write(iout,'(9a)',advance='no') nline, 'id="', trim(tmpChar), '" ', &
      isDefault, ' ', isUsed
    write(iout,101,advance='no') 'name="total" type="Array(double)" value="{'
    do g = 1,ngroup-1
      write(iout,'(e12.6,a)',advance='no') shift_total(g), ','
    enddo
    write(iout,'(e12.6,a)') shift_total(ngroup), '}"/>'

    counter = counter + 1
    write(tmpChar,201) counter
    tmpChar = trim(adjustl(tmpChar))
    write(iout,'(9a)',advance='no') nline, 'id="', trim(tmpChar), '" ', &
      isDefault, ' ', isUsed
    write(iout,101,advance='no') 'name="scattering" type="Array(double)" value="{'
    do g = 1,ngroup-1
      write(iout,'(e12.6,a)',advance='no') shift_scattering(g), ','
    enddo
    write(iout,'(e12.6,a)') shift_scattering(ngroup), '}"/>'
    
    if (any(shift_chi(:) > 0.0d0)) then
      ! is fissionable
      counter = counter + 1
      write(tmpChar,201) counter
      tmpChar = trim(adjustl(tmpChar))
      write(iout,'(9a)',advance='no') nline, 'id="', trim(tmpChar), '" ', &
        isDefault, ' ', isUsed
      write(iout,101,advance='no') 'name="fission" type="Array(double)" value="{'
      do g = 1,ngroup-1
        write(iout,'(e12.6,a)',advance='no') shift_fission(g), ','
      enddo
      write(iout,'(e12.6,a)') shift_fission(ngroup), '}"/>'

      counter = counter + 1
      write(tmpChar,201) counter
      tmpChar = trim(adjustl(tmpChar))
      write(iout,'(9a)',advance='no') nline, 'id="', trim(tmpChar), '" ', &
        isDefault, ' ', isUsed
      write(iout,101,advance='no') 'name="nu_fission" type="Array(double)" value="{'
      do g = 1,ngroup-1
        write(iout,'(e12.6,a)',advance='no') shift_nu_fission(g), ','
      enddo
      write(iout,'(e12.6,a)') shift_nu_fission(ngroup), '}"/>'

      counter = counter + 1
      write(tmpChar,201) counter
      tmpChar = trim(adjustl(tmpChar))
      write(iout,'(9a)',advance='no') nline, 'id="', trim(tmpChar), '" ', &
        isDefault, ' ', isUsed
      write(iout,101,advance='no') 'name="chi" type="Array(double)" value="{'
      do g = 1,ngroup-1
        write(iout,'(e12.6,a)',advance='no') shift_chi(g), ','
      enddo
      write(iout,'(e12.6,a)') shift_chi(ngroup), '}"/>'
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
      do gprime = 1,ngroup
        if ((g == ngroup) .and. (gprime == ngroup)) then
          exit
        endif
        write(iout,'(e12.6,a)',advance='no') shift_sigma_s0(gprime,g), ','
      enddo
    enddo
    write(iout,'(e12.6,a)') shift_sigma_s0(ngroup,ngroup), '}"/>'

    write(iout,101) '</ParameterList>'

    write(iout,101) '</ParameterList>'

  endsubroutine shift_format
endmodule shift
