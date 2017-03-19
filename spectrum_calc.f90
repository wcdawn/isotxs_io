module spectrum_calc
use variables
! type xs_library
	! real(8),allocatable,dimension(:,:) :: sigtr, sigtot
	! real(8),allocatable,dimension(:) :: signg, sigf, nuf, chi, sigalf, sigp, sign2n, sigd, sigt
! endtype

contains

subroutine spectrum_solve
IMPLICIT NONE
! type(xs_library),allocatable,dimension(:),intent(in) :: xs
! integer,intent(in) :: niso, ngroup
integer :: i, j, g
! real(8),allocatable,dimension(:,:) :: sigtr, sigtot
real(8),allocatable,dimension(:)   :: signg, nusigf, chi, sigalf, sigp, sign2n, sigd, sigt
real(8),allocatable,dimension(:)   :: phi, source, xs_total
real(8) :: k, tol, scat_sum

allocate(signg(ngroup))
allocate(nusigf(ngroup))
allocate(sigalf(ngroup))
allocate(sigp(ngroup))
allocate(sign2n(ngroup))
allocate(sigd(ngroup))
allocate(sigt(ngroup))
allocate(phi(ngroup))
allocate(source(ngroup))
allocate(xs_total(ngroup))

! homogenize
! sum of macro xs
! sigtr(:)  = 0.0d0
! sigtot = 0.0d0
signg  = 0.0d0
nusigf = 0.0d0
sigalf = 0.0d0
sigp   = 0.0d0
sign2n = 0.0d0
sigd   = 0.0d0
sigt   = 0.0d0
do j = 1,ngroup
	do i = 1,niso
		signg(j)  = signg(j)  + xs(i)%signg(j)
		nusigf(j) = nusigf(j) + (xs(i)%sigf(j) * xs(i)%nuf(j))
		sigalf(j) = sigalf(j) + xs(i)%sigalf(j)
		sigp(j)   = sigp(j)   + xs(i)%sigp(j)
		sign2n(j) = sign2n(j) + xs(i)%sign2n(j)
		sigd(j)   = sigd(j)   + xs(i)%sigd(j)
		sigt(j)   = sigt(j)   + xs(i)%sigt(j)
	enddo
enddo

! iterate
! intialize
do j = 1,ngroup
	scat_sum = 0.0d0
	do g = 1,ngroup
		scat_sum = scat_sum + 1.0d0
	enddo
	xs_total(j) = signg(j) + sigalf(j) + sigp(j) + sign2n(j) + scat_sum
enddo
phi = 1.0d0
k   = 1.0d0


endsubroutine spectrum_solve

endmodule spectrum_calc