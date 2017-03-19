module text_io
contains

subroutine ascii_out(hname,huse,ivers, &
                     ngroup,niso,maxup,maxdn,maxord,ichist,nscmax,nsblok, &
                     hsetid,hisonm,chi,vel,emax,emin,loca, &
                     chi_fw, isspec, &
                     habsid,hident,hmat,amass,efiss,ecapt,temp,sigpot,adens,kbr,ichi,ifis,ialf,inp,in2n,ind,int,ltot,ltrn,istrpd, &
                     idsct,lord,jband,ijj, &
                     strpl,stotpl,strpd,sngam,sfis,snutot,chiso,snalf,snp,sn2n,snd,snt, &
                     chiiso,isopec, &
                     kmax,scat)
IMPLICIT NONE
integer :: i, j, k, iout, ios = 0
integer :: scat_point
character(80) :: fname_out

! FILE IDENTIFICATION
character(8),intent(in) :: hname
character(8),dimension(2),intent(in) :: huse
integer,intent(in) :: ivers
! FILE CONTROL   (1D RECORD)
integer,intent(in) :: ngroup, niso, maxup, maxdn, maxord, ichist, nscmax, nsblok
! FILE DATA   (2D RECORD)
character(8),dimension(12),intent(in) :: hsetid
character(8),dimension(:),allocatable,intent(in) :: hisonm
real(4),dimension(:),allocatable,intent(in) :: chi, vel, emax
real(4),intent(in) :: emin
integer,dimension(:),allocatable,intent(in) :: loca
! FILE-WIDE CHI DATA   (3D RECORD)
real(4),dimension(:,:),allocatable,intent(in) :: chi_fw
integer,dimension(:),allocatable,intent(in) :: isspec
! ISOTOPE CONTROL AND GROUP INDEPENDENT DATA   (4D RECORD)
character(8),dimension(:),allocatable,intent(in) :: habsid, hident, hmat
real(4),dimension(:),allocatable,intent(in) :: amass, efiss, ecapt, temp, sigpot, adens
integer,dimension(:),allocatable,intent(in) :: kbr,  ichi, ifis, ialf, inp, in2n, ind, int, ltot, ltrn, istrpd
integer,dimension(:,:),allocatable,intent(in) :: idsct, lord
integer,dimension(:,:,:),allocatable,intent(in) :: jband, ijj
! PRINCIPAL CROSS SECTIONS   (5D RECORD)
real(4),dimension(:,:,:),allocatable,intent(in) :: strpl, stotpl, strpd
real(4),dimension(:,:),allocatable,intent(in) :: sngam, sfis, snutot, chiso, snalf, snp, sn2n, snd, snt
! ISOTOPE CHI DATA (6D RECORD)
real(4),dimension(:,:,:),allocatable,intent(in) :: chiiso
integer,dimension(:,:),allocatable,intent(in) :: isopec
! SCATTERING SUB-BLOCK   (7D RECORD)
integer,intent(in) :: kmax
real(4),dimension(:,:,:,:),allocatable,intent(in) :: scat

iout = 21
fname_out = 'ascii.out'
!------------------------------------------------------------------------------!
! OPEN FILES
!------------------------------------------------------------------------------!
! ASCII out
open(unit = iout, file = fname_out, status = 'replace', action = 'write', iostat = ios)
call checkopen(iout,fname_out,ios)

! Formats
! CHARACTER
101 format(a)  ! plain-text descriptor
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
902 format('ISOTOPE CONTROL AND GROUP INDEPENDENT DATA   (4D RECORD)', &
           /, a,     ' HABSID        HOLLERITH ABSOLUTE ISOTOPE LABEL', &
           /, a,     ' HIDENT        IDENTIFIER OF LIBRARY FROM WHICH BASIC DATA CAME', &
           /, a,     ' HMAT          ISOTOPE IDENTIFICATION', &
           /, f6.2, ' AMASS         GRAM ATOMIC WEIGHT', &
           /, e6.1, ' EFISS         TOTAL THERMAL ENERGY YIELD/FISSION (W.SEC/FISS)', &
           /, e6.1, ' ECAPT         TOTAL THERMAL ENERGY YIELD/CAPTURE (W.SEC/CAPT)', &
           /, e6.1, ' TEMP          ISOTOPE TEMPERATURE (DEGREES KELVIN)', &
           /, e6.1, ' SIGPOT        AVERAGE EFFECTIVE POTENTIAL SCATTERING IN RESONANCE RANGE (BARNS/ATOM)', &
           /, e6.1, ' ADENS         DENSITY OF ISOTOPE IN MIXTURE IN WHICH ISOTOPE CROSS SECTIONS WERE GENERATED (A/BARN-CM)', &
           /, i6,    ' KBR           ISOTOPE CLASSIFICATION (SEE DOCUMENTATION)', &
           /, i6,    ' ICHI          ISOTOPE FISSION SPECTRUM FLAG', &
           /, i6,    ' IFIS          (N,F) CROSS SECTION FLAG ', &
           /, i6,    ' IALF          (N,ALPHA) CROSS SECTION FLAG', &
           /, i6,    ' INP           (N,P) CROSS SECTION FLAG', &
           /, i6,    ' IN2N          (N,2N) CROSS SECTION FLAG', &
           /, i6,    ' IND           (N,D) CROSS SECTION FLAG', &
           /, i6,    ' INT           (N,T) CROSS SECTION FLAG', &
           /, i6,    ' LTOT          NUMBER OF MOMENTS OF TOTAL CROSS SECTION PROVIDED', &
           /, i6,    ' LTRN          NUMBER OF MOMENTS OF TRANSPORT CROSS SECTION', &
           /, i6,    ' ISTRPD        NUMBER OF COORDINATE DIRECTIONS ... (MUST .EQ. 0)')

! FILE IDENTIFICATION
write(iout,900) hname, huse(1:2), ivers
! FILE CONTROL   (1D RECORD)
write(iout,901) ngroup, niso, maxup, maxdn, maxord, ichist, nscmax, nsblok
! FILE DATA   (2D RECORD)
write(iout,'(/,a)') 'FILE DATA   (2D RECORD)'
write(iout,101) 'HSETID(I)     HOLLERITH IDENTIFICATION OF FILE'
write(iout,'(12a)') hsetid(1:12)
write(iout,101) 'HISONM(I)     HOLLERITH ISOTOPE LABEL FOR ISOTOPE I'
call write_fivetable_char(hisonm,niso,iout)
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
! FILE-WIDE CHI DATA   (3D RECORD)
if (ichist .gt. 1) then
	write(iout,'(/,a)') 'FILE-WIDE CHI DATA   (3D RECORD)'
	write(iout,101) 'CHI(K,J)      FRACTION OF NEUTRONS EMITTED INTO GROUP J AS A RESULT OF FISSION IN ANY GROUP,USING SPECTRUM K'
	call write_realtable(chi_fw,ichist,ngroup,iout)
	write(iout,101) 'ISSPEC(I)     ISSPEC(I)=K IMPLIES THAT SPECTRUM K IS USED TO CALCULATE EMISSION SPECTRUM FROM FISSION IN GROUP I'
	call write_fivetable_int(isspec,ngroup,iout)
endif
write(iout,101)
write(iout,101) '********************************************************************************'
write(iout,101) 'INDIVIDUAL ISOTOPIC DATA'
write(iout,101) '********************************************************************************'
do i = 1,niso
	write(iout,'(/,a,i3)') 'ISOTOPE ', i
	! ISOTOPE CONTROL AND GROUP INDEPENDENT DATA   (4D RECORD)
	write(iout,902) habsid(i), hident(i), hmat(i), amass(i), efiss(i), ecapt(i), temp(i), sigpot(i), adens(i), &
                    kbr(i), ichi(i), ifis(i), ialf(i), inp(i), in2n(i), ind(i), int(i), ltot(i), ltrn(i), &
                    istrpd(i)
	write(iout,101) 'IDSCT(N)      SCATTERING MATRIX TYPE IDENTIFICATION FOR'
	write(iout,101) 'SCATTERING BLOCK N.  SIGNIFICANT ONLY IF LORD(N).GT.0'
	write(iout,101) 'IDSCT(N)=000 + NN, TOTAL SCATTERING, (SUM OF'
	write(iout,101) '    ELASTIC,INELASTIC, AND N,2N SCATTERING'
	write(iout,101) '    MATRIX TERMS).'
	write(iout,101) '        =100 + NN, ELASTIC SCATTERING'
	write(iout,101) '        =200 + NN, INELASTIC SCATTERING'
	write(iout,101) '        =300 + NN, (N,2N) SCATTERING'
	write(iout,101) 'WHERE NN IS THE LEGENDRE EXPANSION INDEX OF THE'
	write(iout,101) 'FIRST MATRIX IN BLOCK N'
	call write_fivetable_int(idsct(i,:),nscmax,iout)
	write(iout,101) 'LORD(N)       NUMBER OF SCATTERING ORDERS IN BLOCK N.'
	call write_fivetable_int(lord(i,:),nscmax,iout)
	write(iout,101) 'JBAND(J,N)    NUMBER OF GROUPS THAT SCATTER INTO GROUP J, INCLUDING SELF-SCATTER, IN SCATTERING BLOCK N'
	call write_inttable(jband(i,:,:),ngroup,nscmax,iout)
	write(iout,101) 'IJJ(J,N)      POSITION OF IN-GROUP SCATTERING CROSS SECTION IN SCATTERING DATA FOR GROUP J, SCATTERING BLOCK N'
	! PRINCIPAL CROSS SECTIONS   (5D RECORD)
	write(iout,101) 'PRINCIPAL CROSS SECTIONS   (5D RECORD)'
	write(iout,101) 'STRPL(J,L)    PL WEIGHTED TRANSPORT CROSS SECTION '
	call write_realtable(strpl(i,:,:),ngroup,ltrn(i),iout)
	write(iout,101) 'STOTPL(J,L)   PL WEIGHTED TOTAL CROSS SECTION'
	call write_realtable(stotpl(i,:,:),ngroup,ltot(i),iout)
	write(iout,101) 'SNGAM        SFIS         SNUTOT       CHISO        SNALF        SNP          SN2N         SND          SNT'
	do j = 1,ngroup
		write(iout,'(9(e12.6,x))') sngam(i,j), sfis(i,j), snutot(i,j), chiso(i,j), snalf(i,j), snp(i,j), sn2n(i,j), &
		                           snd(i,j), snt(i,j)
	enddo
	! ISOTOPE CHI DATA   (6D RECORD)
	! not tested
	! SCATTERING SUB-BLOCK   (7D RECORD)
	write(iout,101) 'SCATTERING SUB-BLOCK   (7D RECORD)'
	do j = 1,nscmax
		if (lord(i,j) .gt. 0) then
			write(iout,'(a,i3,a,i3)') 'block id ', idsct(i,j), ' order ', lord(i,j)
			call write_fivetable_real(scat(i,j,:,lord(i,j)),kmax,iout)
			scat_point = 0
			do k = 1,ngroup
				write(iout,*) 'scattering into group ', k, ' from position ', scat_point + 1, ' through ', scat_point + jband(i,k,j)
				call write_fivetable_real(scat(i,j,scat_point + 1:scat_point + jband(i,k,j),lord(i,j)),jband(i,k,j),iout)
				scat_point = scat_point + jband(i,k,j)
			enddo
		endif
	enddo
enddo
	












endsubroutine ascii_out

subroutine checkopen(ifl,fname,ios)
IMPLICIT NONE
integer,intent(in) :: ifl, ios
character(80),intent(in) :: fname
if (ios .ne. 0) then
	write(*,'(a,i6,a,a)') 'FATAL -- error opening unit -- ', ifl, ' -- ', fname
	write(*,'(a,i6)') 'ios = ', ios
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
if ((i .ne. length) .and. (r .ne. 5)) write(iout,*)
enddo
endsubroutine write_fivetable_real

subroutine write_fivetable_char(var,length,iout)
IMPLICIT NONE
character(8),dimension(:),intent(in) :: var
integer,intent(in) :: length, iout
integer :: r, c, i
i = 0
do r = 1,((length / 5) + 1)
	do c = 1,5
		i = i + 1
		if (i .gt. length) exit
		write(iout,'(i3,x,a,3x)',advance = 'no') i, var(i)
	enddo
if ((i .ne. length) .and. (r .ne. 5)) write(iout,*)
enddo
endsubroutine write_fivetable_char

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
if ((i .ne. length) .and. (r .ne. 5)) write(iout,*)
enddo
endsubroutine write_fivetable_int

subroutine write_inttable(var,dim1,dim2,iout)
IMPLICIT NONE
integer,dimension(:,:),intent(in) :: var
integer,intent(in) :: dim1, dim2, iout
integer :: i, j
do i = 1,dim1
	do j = 1,dim2
		write(iout,'(i3,x)',advance = 'no') var(i,j)
	enddo
	write(iout,*)
enddo
endsubroutine write_inttable

subroutine write_realtable(var,dim1,dim2,iout)
IMPLICIT NONE
real(4),dimension(:,:),intent(in) :: var
integer,intent(in) :: dim1, dim2, iout
integer :: i, j
do i = 1,dim1
	do j = 1,dim2
		write(iout,'(e12.6,x)',advance = 'no') var(i,j)
	enddo
	write(iout,*)
enddo
endsubroutine write_realtable

endmodule text_io