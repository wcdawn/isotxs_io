module mpact_interface
use variables
use spectrum_calc
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
		write(iout,'(4(e12.6,x))') xs(i)%mpact_abs(g), xs(i)%nuf(g) * xs(i)%sigf(g), xs(i)%kappa * xs(i)%sigf(g), xs(i)%chi(g)
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
			write(iout,'(e12.6,x)', advance = 'no') xs(i)%mpact_scat(j,g)
		enddo
		write(iout,*)
	enddo
	
	write(iout,*)
	
enddo
close(unit = iout)


endsubroutine mpact_format

subroutine mpact_homogenize(iout,fname,library_name)
IMPLICIT NONE
integer,intent(in) :: iout
character(80),intent(in) :: fname, library_name

integer :: i, g, j
integer :: ios = 0

101 format (a) ! plain-text descriptor

open(unit = iout, file = fname, status = 'replace', action = 'write', iostat = ios)
if (ios .ne. 0) then
	write(*,'(a,i6,2a)') 'FATAL -- error opening unit -- ', iout, ' -- ', fname
	write(*,'(a,i6)') 'ios = ', ios
endif

! Line 1: Title
write(iout,101) library_name

! Line 2: number_of_gropus number_of_xsSets [number_of_delayed_groups]
write(iout,'(2i6)') ngroup, 1

! Line 3: Upper energy bound in [ev]
! group1_ubound group2_ubound ... groupNG_ubound
call write_list(iout,emax,ngroup)
write(iout,*)

! XSMACRO xs_name scat_order
write(iout,'(a,x,a,x,i1)') 'XSMACRO', 'homog', 0

if (.not. allocated(chi_tilde)) then
	call spectrum_solve()
endif

! perform homogenization
! the only one not performed here is chi_tilde becasue this requires the spectrum solution
allocate(mpact_abs(ngroup))
allocate(mpact_nusigf(ngroup))
allocate(mpact_kappasigf(ngroup))
allocate(mpact_scat(ngroup,ngroup))
mpact_abs(:) = 0.0d0
mpact_nusigf(:) = 0.0d0
mpact_kappasigf(:) = 0.0d0
mpact_scat(:,:) = 0.0d0
do i = 1,niso
	do g = 1,ngroup
		mpact_abs(g) = mpact_abs(g) + xs(i)%mpact_abs(g)
		mpact_nusigf(g) = mpact_nusigf(g) + (xs(i)%nuf(g) * xs(i)%sigf(g))
		mpact_kappasigf(g) = mpact_kappasigf(g) + (xs(i)%kappa * xs(i)%sigf(g))
		do j = 1,ngroup
			mpact_scat(g,j) = mpact_scat(g,j) + xs(i)%mpact_scat(g,j)
		enddo
	enddo
enddo
! zero out fission if non-fissile
if (vector_sum(chi_tilde,ngroup) .le. 1e-3) then
	mpact_kappasigf(:) = 0.0d0
	mpact_nusigf(:) = 0.0d0
endif
do g = 1,ngroup
	! abs1 nfs1 kfs1 chi1
	write(iout,'(4(e12.6,x))') mpact_abs(g), mpact_nusigf(g), mpact_kappasigf(g), chi_tilde(g)
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
		write(iout,'(e12.6,x)', advance = 'no') mpact_scat(j,g)
	enddo
	write(iout,*)
enddo
close(unit = iout)

endsubroutine mpact_homogenize

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