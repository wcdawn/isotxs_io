program isotxsio
IMPLICIT NONE

! General control variables
integer :: i, j, r, c
integer :: ifl, iout, ios
character(80) :: fname, fname_out
logical :: lfixstr

! FILE IDENTIFICATION
character(8) :: hname
character(8),dimension(2) :: huse
integer :: ivers
! FILE CONTROL   (1D RECORD)
integer :: ngroup, niso, maxup, maxdn, maxord, ichist, nscmax, nsblok
! FILE DATA   (2D RECORD)
character(8),dimension(12) :: hsetid
character(8),dimension(:),allocatable :: hisonm
real(4),dimension(:),allocatable :: chi, vel, emax
real(4) :: emin
integer,dimension(:),allocatable :: loca

! Formats
! CHARACTER
101 format(a)  ! plain-text descriptor
! INTEGER
201 format(i6) ! integer*6
! RECORDS
900 format(/, 'FILE IDENTIFICATION', &
           /, a,  ' HNAME         HOLLERITH FILE NAME - ISOTXS', &
           /, a,  ' HUSE(1)       HOLLERITH USER IDENTIFICATION', &
           /, a,  ' HUSE(2)       HOLLERITH USER IDENTIFICATION', &
           /, i8, ' IVERS         FILE VERSION NUMBER')
901 format(/, 'FILE CONTROL   (1D RECORD)', &
           /, i6, ' NGROUP        NUMBER OF ENERGY GROUPS IN FILE', &
           /, i6, ' NISO          NUMBER OF ISOTOPES IN FILE', &
           /, i6, ' MAXUP         MAXIMUM NUMBER OF UPSCATTER GROUPS', &
           /, i6, ' MAXDN         MAXIMUM NUMBER OF DOWNSCATTER GROUPS', &
           /, i6, ' MAXORD        MAXIMUM SCATTERING ORDER', &
           /, i6, ' ICHIST        FILE-WIDE FISSION SPECTRUM FLAG', &
           /, i6, ' NSCMAX        MAXIMUM NUMBER OF BLOCKS OF SCATTERING DATA', &
           /, i6, ' NSBLOK        SUBBLOCKING CONTROL FOR SCATTER MATRICES')

ifl = 11
fname = 'ISOTXS.20'
iout = 21
fname_out = 'ascii.out'
lfixstr = .false.

!------------------------------------------------------------------------------!
! OPEN FILES
!------------------------------------------------------------------------------!
! ISOTXS in
! default open with little_endian. May be closed and reopened later
open(ifl, file = fname, status = 'old', form = 'unformatted', action = 'read', iostat = ios)
call checkopen(ifl,fname,ios)

! ASCII out
open(unit = iout, file = fname_out, status = 'replace', action = 'write', iostat = ios)
call checkopen(iout,fname_out,ios)

!------------------------------------------------------------------------------!
! READ FILE IDENTIFICATION
!------------------------------------------------------------------------------!
read(ifl) hname, (huse(i), i = 1,2), ivers
if ((ivers .lt. 0) .or. (ivers .gt. 20)) then
	! need big-endian open
	write(*,101) 'closing ISOTXS and opening with BIG_ENDIAN'
	close(ifl)
	open(ifl, file = fname, status = 'old', form = 'unformatted', action = 'read', convert = 'BIG_ENDIAN', iostat = ios)
	read(ifl) hname, (huse(i), i = 1,2), ivers	
endif
	

! Check if characters need fixing
if (hname .eq. 'TOSI  SX') then
	lfixstr = .true.
	write(*,101) 'fixing strings'
endif
! fix characters if necessary
if (lfixstr) then
	call fixstr(hname)
	call fixstr(huse(1))
	call fixstr(huse(2))
endif

write(iout,900) hname, (huse(i), i = 1,2), ivers

! make sure this is an ISOTXS file
if (hname .ne. 'ISOTXS  ') then
	write(*,101) 'FATAL -- invalid file name'
	write(*,101) hname
	stop
endif

!------------------------------------------------------------------------------!
! READ FILE CONTROL   (1D RECORD)
!------------------------------------------------------------------------------!
read(ifl,iostat = ios) ngroup, niso, maxup, maxdn, maxord, ichist, nscmax, nsblok
write(iout,901) ngroup, niso, maxup, maxdn, maxord, ichist, nscmax, nsblok
if (ios .ne. 0) then
	write(*,101) 'FATAL -- error reading FILE CONTROL   (1D RECORD)'
	write(*,'(a,i6)') 'return code = ', ios
	stop
endif

write(*,101) 'allocating memory'
allocate(hisonm(niso))
allocate(chi(ngroup))
allocate(vel(ngroup))
allocate(emax(ngroup))
allocate(loca(niso))

!------------------------------------------------------------------------------!
! READ FILE DATA   (2D RECORD)
!------------------------------------------------------------------------------!
if (ichist .eq. 1) then
	read(ifl,iostat = ios) (hsetid(i), i = 1,12), (hisonm(i), i = 1,niso), (chi(j), j = 1,ngroup), (vel(j), j = 1,ngroup),&
	                       (emax(j), j = 1,ngroup), emin, (loca(i), i = 1,niso)
else
	read(ifl,iostat = ios) (hsetid(i), i = 1,12), (hisonm(i), i = 1,niso), (vel(j), j = 1,ngroup), &
	                       (emax(j), j = 1,ngroup), emin, (loca(i), i = 1,niso)
endif

if (lfixstr) then
	do i = 1,12
		call fixstr(hsetid(i))
	enddo
	do i = 1,niso
		call fixstr(hisonm(i))
	enddo
endif

write(iout,'(/,a)') 'FILE DATA   (2D RECORD)'
write(iout,101) 'HSETID(I)     HOLLERITH IDENTIFICATION OF FILE'
write(iout,'(12a)') (hsetid(i), i = 1,12) 
write(iout,101) 'HISONM(I)     HOLLERITH ISOTOPE LABEL FOR ISOTOPE I'
i = 0
do r = 1,((niso / 5) + 1)
	do c = 1,5
		i = i + 1
		if (i .gt. niso) exit
		write(iout,'(i3,x,a,3x)',advance = 'no') i, hisonm(i)
	enddo
	write(iout,*)
enddo
if (ichist .eq. 1) then
	write(iout,101) 'CHI(J)        FILE-WIDE FISSION SPECTRUM'
	call write_fivetable_real(chi,ngroup,iout)
endif
write(iout,101) 'VEL(J)        MEAN NEUTRON VELOCITY IN GROUP J (CM/SEC)'
call write_fivetable_real(vel,ngroup,iout)
write(iout,101) 'EMAX(J)       MAXIMUM ENERGY BOUND OF GROUP J (EV)'
call write_fivetable_real(emax,ngroup,iout)
write(iout,101) 'EMIN          MINIMUM ENERGY BOUND OF SET (EV)'
write(iout,'(e12.6)') emin
write(iout,101) 'LOCA(I)       NUMBER OF RECORDS TO BE SKIPPED TO READ DATA FOR ISOTOPE I'
call write_fivetable_int(loca,ngroup,iout)

close(ifl)
contains


subroutine fixstr(str)
IMPLICIT NONE
character(8) :: str, str0
str0=str
str(1:1)=str0(4:4) 
str(2:2)=str0(3:3) 
str(3:3)=str0(2:2) 
str(4:4)=str0(1:1) 
str(5:5)=str0(8:8) 
str(6:6)=str0(7:7) 
str(7:7)=str0(6:6) 
str(8:8)=str0(5:5)
endsubroutine fixstr

subroutine checkopen(ifl,fname,ios)
IMPLICIT NONE
integer,intent(in) :: ifl, ios
character(80),intent(in) :: fname
if (ios .ne. 0) then
	write(*,'(a,i2,a,a)') 'FATAL -- error opening unit -- ', ifl, ' -- ', fname
	write(*,'(a,i3)') 'ios = ', ios
	stop
endif
endsubroutine checkopen

subroutine write_fivetable_real(var,length,iout)
IMPLICIT NONE
real(4),dimension(:),intent(in) :: var
integer,intent(in) :: length, iout
integer :: r, c, i
i = 0
do r = 1,((length / 5) + 1)
	do c = 1,5
		i = i + 1
		if (i .gt. length) exit
		write(iout,'(i3,x,e12.6,3x)',advance = 'no') i, var(i)
	enddo
write(iout,*)
enddo
endsubroutine write_fivetable_real

subroutine write_fivetable_int(var,length,iout)
IMPLICIT NONE
integer,dimension(:),intent(in) :: var
integer,intent(in) :: length, iout
integer :: r, c, i
i = 0
do r = 1,((length / 5) + 1)
	do c = 1,5
		i = i + 1
		if (i .gt. length) exit
		write(iout,'(i3,x,i6,3x)',advance = 'no') i, var(i)
	enddo
write(iout,*)
enddo
endsubroutine write_fivetable_int

endprogram isotxsio