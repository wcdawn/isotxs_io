module spectrum_calc
contains

subroutine spectrum_solve(xs,niso,ngroup)
IMPLICIT NONE
type xs_library
	real(8),allocatable,dimension(:,:) :: sigtr, sigtot
	real(8),allocatable,dimension(:) :: signg, sigf, nuf, chi, sigalf, sigp, sign2n, sigd, sigt
endtype
type(xs_library) xs(niso)
integer,intent(in) :: niso, ngroup
! intent(in) :: xs

! homogenize

! iterate

endsubroutine spectrum_solve

endmodule spectrum_calc