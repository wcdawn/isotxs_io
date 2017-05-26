module variables
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
! FILE-WIDE CHI DATA   (3D RECORD)
real(4),dimension(:,:),allocatable :: chi_fw
integer,dimension(:),allocatable :: isspec
! ISOTOPE CONTROL AND GROUP INDEPENDENT DATA   (4D RECORD)
character(8),dimension(:),allocatable :: habsid, hident, hmat
real(4),dimension(:),allocatable :: amass, efiss, ecapt, temp, sigpot, adens
integer,dimension(:),allocatable :: kbr,  ichi, ifis, ialf, inp, in2n, ind, int, ltot, ltrn, istrpd
integer,dimension(:,:),allocatable :: idsct, lord
integer,dimension(:,:,:),allocatable :: jband, ijj
! PRINCIPAL CROSS SECTIONS   (5D RECORD)
real(4),dimension(:,:,:),allocatable :: strpl, stotpl, strpd
real(4),dimension(:,:),allocatable :: sngam, sfis, snutot, chiso, snalf, snp, sn2n, snd, snt
! ISOTOPE CHI DATA (6D RECORD)
real(4),dimension(:,:,:),allocatable :: chiiso
integer,dimension(:,:),allocatable :: isopec
! SCATTERING SUB-BLOCK   (7D RECORD)
integer :: kmax
real(4),dimension(:,:,:,:),allocatable :: scat
! xs_structure
type xs_library
	real(8),allocatable,dimension(:,:) :: sigtr, sigtot, scat, n2n
	real(8),allocatable,dimension(:) :: sigtotp0, signg, sigf, nuf, chi, sigalf, sigp, sign2n, sigd, sigt
endtype
type(xs_library),allocatable,dimension(:) :: xs

contains

subroutine allocate_memory(niso,ngroup,nscmax,ichist)
IMPLICIT NONE
integer,intent(in) :: niso,ngroup,nscmax,ichist
! FILE DATA   (2D RECORD)
allocate(hisonm(niso))
if (ichist .eq. 1) then
	allocate(chi(ngroup))
endif
allocate(vel(ngroup))
allocate(emax(ngroup))
allocate(loca(niso))
! FILE-WIDE CHI DATA   (3D RECORD)
if (ichist .gt. 1) then
	allocate(chi_fw(ichist,ngroup))
	allocate(isspec(ngroup))
endif
! ISOTOPE CONTROL AND GROUP INDEPENDENT DATA   (4D RECORD)
allocate(habsid(niso))
allocate(hident(niso))
allocate(hmat(niso))
allocate(amass(niso))
allocate(efiss(niso))
allocate(ecapt(niso))
allocate(temp(niso))
allocate(sigpot(niso))
allocate(adens(niso))
allocate(kbr(niso))
allocate(ichi(niso))
allocate(ifis(niso))
allocate(ialf(niso))
allocate(inp(niso))
allocate(in2n(niso))
allocate(ind(niso))
allocate(int(niso))
allocate(ltot(niso))
allocate(ltrn(niso))
allocate(istrpd(niso))
allocate(idsct(niso,nscmax))
allocate(lord(niso,nscmax))
allocate(jband(niso,ngroup,nscmax))
allocate(ijj(niso,ngroup,nscmax))
! PRINCIPAL CROSS SECTIONS   (5D RECORD)
allocate(strpl(niso,ngroup,nscmax))
allocate(stotpl(niso,ngroup,nscmax))
allocate(sngam(niso,ngroup))
allocate(sfis(niso,ngroup))
allocate(snutot(niso,ngroup))
allocate(snalf(niso,ngroup))
allocate(chiso(niso,ngroup))
allocate(snp(niso,ngroup))
allocate(sn2n(niso,ngroup))
allocate(snd(niso,ngroup))
allocate(snt(niso,ngroup))
! ISOTOPE CHI DATA   (6D RECORD)
allocate(chiiso(niso,ichist,ngroup))
allocate(isopec(niso,ngroup))
! SCATTERING SUB-BLOCK   (7D RECORD)
! TO-DO : ADDRESS THE SIZE OF THIS (IN DIMENSION 4)
allocate(scat(niso,nscmax,ngroup * ngroup,3))
endsubroutine allocate_memory

subroutine xs_structure
IMPLICIT NONE
integer :: i,j,k
integer :: point, group_offset, group_start, group_end

allocate(xs(niso))
do i = 1,niso
	allocate(xs(i)%sigtr(ngroup,nscmax))
	allocate(xs(i)%sigtot(ngroup,nscmax))
	allocate(xs(i)%sigtotp0(ngroup))
	allocate(xs(i)%signg(ngroup))
	allocate(xs(i)%sigf(ngroup))
	allocate(xs(i)%nuf(ngroup))
	allocate(xs(i)%chi(ngroup))
	allocate(xs(i)%sigalf(ngroup))
	allocate(xs(i)%sigp(ngroup))
	allocate(xs(i)%sign2n(ngroup))
	allocate(xs(i)%sigd(ngroup))
	allocate(xs(i)%sigt(ngroup))
	allocate(xs(i)%scat(ngroup,ngroup))
	allocate(xs(i)%n2n(ngroup,ngroup))
enddo

do i = 1,niso
	if (abs(adens(i) - 1.0d0) .lt. 1.0d-2) then
		write(*,'(a)') 'WARNING -- adens(i) .eq. 1.0d0'
		write(*,'(a)') 'are you working with macro or micro xs?'
		write(*,'(a,i3)') 'isotope ', i
		write(*,'(a,e12.6)') 'adens ', adens(i)
	endif
	
	xs(i)%sigtr(:,:)  = strpl(i,:,:)  * adens(i)
	xs(i)%sigtot(:,:) = stotpl(i,:,:) * adens(i)
	xs(i)%sigtotp0(:) = stotpl(i,:,1) * adens(i)
	xs(i)%signg(:)    = sngam(i,:)    * adens(i)
	xs(i)%sigf(:)     = sfis(i,:)     * adens(i)
	xs(i)%nuf(:)      = snutot(i,:)             
	xs(i)%chi(:)      = chiso(i,:)              
	xs(i)%sigalf(:)   = snalf(i,:)    * adens(i)
	xs(i)%sigp(:)     = snp(i,:)      * adens(i)
	xs(i)%sign2n(:)   = sn2n(i,:)     * adens(i)
	xs(i)%sigd(:)     = snd(i,:)      * adens(i)
	xs(i)%sigt(:)     = snt(i,:)      * adens(i)
	
	xs(i)%scat(:,:) = 0.0d0
	xs(i)%n2n(:,:)  = 0.0d0
	
	do j = 1,nscmax
		point = 0
		do k = 1,ngroup
			group_offset = ijj(i,k,j) - 1
			group_end    = k - group_offset
			group_start  = group_end - jband(i,k,j) + 1
			if (idsct(i,j) .ge. 300) then
				! n2n
				if ((idsct(i,j) - 300 .ne. 0)) then
					point = point + jband(i,k,j)
					cycle
				endif
				! write(*,*) group_start, group_end, (group_end - group_start)
				! write(*,*) point + 1, point + jband(i,k,j), ((point + jband(i,k,j)) - (point + 1))
				xs(i)%n2n(k,group_start:group_end) = scat(i,j,point + 1:point + jband(i,k,j),lord(i,j))
			elseif (idsct(i,j) .ge. 200) then
				! inelastic
				if ((idsct(i,j) - 200 .ne. 0)) then
					point = point + jband(i,k,j)
					cycle
				endif
				xs(i)%scat(k,group_start:group_end) = xs(i)%scat(k,group_start:group_end) + &
				scat(i,j,point + 1:point + jband(i,k,j),lord(i,j))
			elseif (idsct(i,j) .ge. 100) then
				! elastic
				if ((idsct(i,j) - 100 .ne. 0)) then
					point = point + jband(i,k,j)
					cycle
				endif
				xs(i)%scat(k,group_start:group_end) = xs(i)%scat(k,group_start:group_end) + &
				scat(i,j,point + 1:point + jband(i,k,j),lord(i,j))
			else
				! write(*,'(a)') 'FATAL -- unsupported idsct value'
				! write(*,'(a,i3)') 'isotope', i
				! write(*,'(a,i3)') 'idsct', idsct(i,j)
				! stop
			endif
			point = point + jband(i,k,j)
		enddo
	enddo
	write(20,*)
	do j = 1,ngroup
		write(20,'(2e13.6)') xs(i)%scat(:,j)
	enddo
enddo

endsubroutine xs_structure

endmodule variables