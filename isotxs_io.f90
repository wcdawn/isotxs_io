program isotxsio
use text_io
use variables
use spectrum_calc
use mpact_interface
IMPLICIT NONE

! General control variables
integer :: i, j, k, r, c
integer :: ifl, ios = 0, mpact_iout
integer :: ichiread, ifisread, ialfread, inpread, in2nread, indread, intread ! reading checks
character(80) :: fname, mpact_fname, mpact_library_name
logical :: lfixstr, lascii, lspectrum, lmpact, lmpact_homog

! Formats
! CHARACTER
101 format(a)  ! plain-text descriptor
ifl = 11
! fname = '16.4_Fuel.ISOTXS_complete'
! fname = 'ISOTXS.20'
fname = 'ISOTXS.soft_reflector'
! fname = 'ISOTXS.u235'
lfixstr   = .false.
lascii    = .true.
lspectrum = .false.
lmpact    = .true.
lmpact_homog = .true.

!------------------------------------------------------------------------------!
! OPEN FILES
!------------------------------------------------------------------------------!
! ISOTXS in
! default open with little_endian. May be closed and reopened later
open(ifl, file = fname, status = 'old', form = 'unformatted', action = 'read', convert = 'LITTLE_ENDIAN', iostat = ios)
call checkopen(ifl,fname,ios)

!------------------------------------------------------------------------------!
! READ FILE IDENTIFICATION
!------------------------------------------------------------------------------!
read(ifl) hname, huse(1:2), ivers
if ((ivers .lt. 0) .or. (ivers .gt. 20)) then
	! need big-endian open
	write(*,101) 'closing ISOTXS and opening with BIG_ENDIAN'
	close(ifl)
	open(ifl, file = fname, status = 'old', form = 'unformatted', action = 'read', convert = 'BIG_ENDIAN', iostat = ios)
	read(ifl) hname, huse(1:2), ivers	
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
if (ios .ne. 0) then
	write(*,101) 'FATAL -- error reading FILE CONTROL   (1D RECORD)'
	write(*,'(a,i6)') 'return code = ', ios
	stop
endif

! TO-DO: If this is a problem, need to add ANOTHER loop to scattering block processing
if (nsblok .ne. 1) then
	write(*,101) 'FATAL -- nsblok must .eq. 1'
	write(*,'(a,i3)') 'nsblok = ', nsblok
	stop
endif

write(*,101) 'allocating memory'
call allocate_memory(niso,ngroup,nscmax,ichist)

! READ PRINCIPAL CROSS SECTIONS   (5D RECORD)
! TO-DO: Is this ever going to be important?
! CD    ISTRPD        NUMBER OF COORDINATE DIRECTIONS FOR WHICH          -
! CD                     COORDINATE DEPENDENT TRANSPORT CROSS SECTIONS   -
! CD                     ARE GIVEN. IF ISTRPD=0, NO COORDINATE DEPENDENT -
! CD                     TRANSPORT CROSS SECTIONS ARE GIVEN.             -
! CD    STRPD(J,I)    COORDINATE DIRECTION I TRANSPORT CROSS SECTION     -
! CD                               (PRESENT IF ISTRPD.GT.0)              -
! allocate(strpd(niso,ngroup,maxval(istrpd)))


!------------------------------------------------------------------------------!
! READ FILE DATA   (2D RECORD)
!------------------------------------------------------------------------------!
if (ichist .eq. 1) then
	ichiread = ngroup
else
	ichiread = 0
endif
read(ifl,iostat = ios) hsetid(1:12), hisonm(1:niso), chi(1:ichiread), vel(1:ngroup), emax(1:ngroup), emin, loca(1:niso)
if (ios .ne. 0) then
	write(*,101) 'FATAL -- error reading 2D RECORD'
	write(*,'(a,i6)') 'error code ', ios
	stop
endif

if (lfixstr) then
	do i = 1,12
		call fixstr(hsetid(i))
	enddo
	do i = 1,niso
		call fixstr(hisonm(i))
	enddo
endif

!------------------------------------------------------------------------------!
! READ FILE-WIDE CHI DATA   (3D RECORD)
! TO-DO: This is untested. No file from the test suite has this.
!------------------------------------------------------------------------------!
if (ichist .gt. 1) then
	write(*,101) 'WARNING -- This is untested. No file from the test suite has this.'
	read(ifl,iostat = ios) chi_fw(1:ichist,1:ngroup), isspec(1:ngroup)
	if (ios .ne. 0) then
		write(*,101) 'FATAL -- error reading 3D RECORD'
		write(*,'(a,i6)') 'error code ', ios
		stop
	endif
endif

! *************(REPEAT FOR ALL ISOTOPES)                   
! *         ISOTOPE CONTROL AND GROUP                      
! *                        INDEPENDENT DATA    ALWAYS      
! *         PRINCIPAL CROSS SECTIONS           ALWAYS      
! *         ISOTOPE CHI DATA                   ICHI.GT.1   
! *                                                        
! *  **********(REPEAT TO NSCMAX SCATTERING BLOCKS)        
! *  *  *******(REPEAT FROM 1 TO NSBLOK)                   
! *  *  *   SCATTERING SUB-BLOCK               LORD(N).GT.0
! *************                                            
do i = 1,niso
	!------------------------------------------------------------------------------!
	! READ ISOTOPE CONTROL AND GROUP INDEPENDENT DATA   (4D RECORD)
	!------------------------------------------------------------------------------!
	read(ifl,iostat = ios) habsid(i), hident(i), hmat(i), amass(i), efiss(i), ecapt(i), temp(i), sigpot(i), adens(i), &
	                       kbr(i), ichi(i), ifis(i), ialf(i), inp(i), in2n(i), ind(i), int(i), ltot(i), ltrn(i), &
	                       istrpd(i), idsct(i,1:nscmax), lord(i,1:nscmax), jband(i,1:ngroup,1:nscmax), &
	                       ijj(i,1:ngroup,1:nscmax)
	if (ios .ne. 0) then
		write(*,101) 'FATAL -- error reading 4D RECORD'
		write(*,'(a,i3,x,a)') 'isotope ', i, habsid(i)
		write(*,'(a,i6)') 'error code ', ios
		stop
	endif
	if (istrpd(i) .ne. 0) then
		write(*,101) 'FATAL -- istrpd(i) .ne. 0'
		write(*,101) 'This mode is not supported. See 4D RECORD'
		write(*,101) 'CD    ISTRPD        NUMBER OF COORDINATE DIRECTIONS FOR WHICH          -'
		write(*,101) 'CD                     COORDINATE DEPENDENT TRANSPORT CROSS SECTIONS   -'
		write(*,101) 'CD                     ARE GIVEN. IF ISTRPD=0, NO COORDINATE DEPENDENT -'
		write(*,101) 'CD                     TRANSPORT CROSS SECTIONS ARE GIVEN.             -'
		write(*,'(a,i3,x,a)') 'isotope ', i, habsid(i)
		stop
	endif
	if (lfixstr) then
		call fixstr(habsid(i))
		call fixstr(hident(i))
		call fixstr(hmat(i))
	endif
	
	!------------------------------------------------------------------------------!
	! READ PRINCIPAL CROSS SECTIONS   (5D RECORD)
	!------------------------------------------------------------------------------!
	ifisread = 0
	ichiread = 0
	ialfread = 0
	inpread  = 0
	in2nread = 0
	indread  = 0
	intread  = 0
	if (ifis(i) .gt. 0) ifisread = ngroup
	if (ichi(i) .eq. 1) ichiread = ngroup
	if (ialf(i) .gt. 0) ialfread = ngroup
	if (inp(i)  .gt. 0) inpread  = ngroup
	if (in2n(i) .gt. 0) in2nread = ngroup
	if (ind(i)  .gt. 0) indread  = ngroup
	if (int(i)  .gt. 0) intread  = ngroup
	read(ifl,iostat = ios) strpl(i,1:ngroup,1:ltrn(i)), stotpl(i,1:ngroup,1:ltot(i)), sngam(i,1:ngroup), sfis(i,1:ifisread), &
	                       snutot(i,1:ifisread), chiso(i,1:ichiread), snalf(i,1:ialfread), snp(i,1:inpread), sn2n(i,1:in2nread), &
	                       snd(i,1:indread), snt(i,1:intread)
	if (ios .ne. 0) then
		write(*,101) 'FATAL -- error reading 5D RECORD'
		write(*,'(a,i3,x,a)') 'isotope ', i, habsid(i)
		write(*,'(a,i6)') 'error code ', ios
		stop
	endif
	
	!------------------------------------------------------------------------------!
	! READ ISOTOPE CHI DATA   (6D RECORD)
	! TO-DO: This is untested. No file from the test suite has this.
	!------------------------------------------------------------------------------!
	if (ichi(i) .gt. 1) then
		write(*,101) 'WARNING -- This is untested. No file from the test suite has this.'
		read(ifl,iostat = ios) chiiso(i,1:ichi(i),1:ngroup), isopec(i,1:ngroup)
		if (ios .ne. 0) then
			write(*,101) 'FATAL -- error reading 6D RECORD'
			write(*,'(a,i3,x,a)') 'isotope ', i, habsid(i)
			write(*,'(a,i6)') 'error code ', ios
			stop
		endif
	endif
	
	!------------------------------------------------------------------------------!
	! READ SCATTERING SUB-BLOCK   (7D RECORD)
	!------------------------------------------------------------------------------!
	do j = 1,nscmax
		if (lord(i,j) .gt. 0) then
			kmax = 0
			do k = 1,ngroup
				kmax = kmax + jband(i,k,j)
			enddo
			read(ifl,iostat = ios) scat(i,j,1:kmax,1:lord(i,j))
			if (ios .ne. 0) then
				write(*,101) 'FATAL -- error reading 7D RECORD'
				write(*,'(a,i3,x,a)') 'isotope ', i, habsid(i)
				write(*,'(a,i3)') 'scattering matrix # ', j
				write(*,'(a,i6)') 'error code ', ios
				stop
			endif
		endif
	enddo	
enddo

if (lascii) then
	write(*,101) 'writing ascii output'
	call ascii_out()
endif

write(*,101) 'building xs structure'
call xs_structure()

if (lspectrum) then
	write(*,101) 'spectral calc'
	call spectrum_solve()
endif

if (lmpact) then
	write(*,101) 'writing to MPACT user format'
	mpact_fname = 'soft_reflector.xsl'
	mpact_library_name = 'soft reflector MPACT'
	! mpact_fname = 'u235.xsl'
	! mpact_library_name = 'u235 single isotope test'
	mpact_iout = 22
	call mpact_format(mpact_iout,mpact_fname,mpact_library_name)
endif

if (lmpact_homog) then
	write(*,101) 'homogenizing and writing to MPACT user format'
	mpact_fname = 'soft_reflector_homog.xsl'
	mpact_library_name = 'soft reflector MPACT homog'
	mpact_iout = 23
	call mpact_homogenize(mpact_iout,mpact_fname,mpact_library_name)
endif

close(ifl)
endprogram isotxsio

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