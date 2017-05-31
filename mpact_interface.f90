module mpact_interface
use variables
contains

subroutine mpact_format(iout,fname,library_name)
IMPLICIT NONE
integer,intent(in) :: iout
character(80),intent(in) :: fname, library_name

integer :: i, g, j, ios = 0

101 format (a) ! plain-text descriptor

open(unit = iout, file = fname, status = 'replace', action = 'write', iostat = ios)
if (ios .ne. 0) then
	write(*,'(a,i6,2a)') 'FATAL -- error opening unit -- ', iout, ' -- ', fname
	write(*,'(a,i6)') 'ios = ', ios
endif

! Line 1: Title
write(iout,101) library_name

! Line 2: number_of_gropus number_of_xsSets [number_of_delayed_groups]
! TO-DO: (niso == xs sets) ???
write(iout,'(2i6)') ngroup, niso

! Line 3: Upper energy bound in [ev]
! group1_ubound group2_ubound ... groupNG_ubound
call write_list(iout,emax,ngroup)

write(iout,*)
do i = 1,niso
	! XSMACRO xs_name scat_order
	write(iout,'(a,x,a,x,i1)') 'XSMACRO', hisonm(i), 0

	! abs1 nfs1 kfs1 chi1
	do g = 1,ngroup
		write(iout,'(4(e12.6,x))') xs(i)%sigabs(g), xs(i)%nuf(g) * xs(i)%sigf(g), xs(i)%kappa * xs(i)%sigf(g), xs(i)%chi(g)
	enddo
	
	! scat0_1->1 scat0_2->1 ... scat0_NG->1
	! scat0_1->2 scat0_2->2 ... scat0_NG->2
	! ...
	! scat0_1->NG scat0_2->NG ... scat0_NG->NG
	! scat1_1->1 scat1_2->1 ... scat1_NG->1
	! scat1_1->2 scat1_2->2 ... scat1_NG->2
	! ...
	! scat1_1->NG scat1_2->NG ... scat1_NG->NG
	do g = 1,ngroup
		do j = 1,ngroup
			write(iout,'(e12.6,x)', advance = 'no') xs(i)%scat(j,g)
		enddo
		write(iout,*)
	enddo
	
	write(iout,*)
enddo


endsubroutine mpact_format

subroutine write_list(iout,vector,length)
IMPLICIT NONE
integer,intent(in) :: iout, length
real(4),dimension(:) :: vector

integer :: i

do i = 1,length
	write(iout,'(e12.6,x)',advance = 'no') vector(i)
enddo
write(iout,*)

endsubroutine write_list

endmodule mpact_interface