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
integer :: i, j, g, gprime
integer :: iteration
! real(8),allocatable,dimension(:,:) :: sigtr, sigtot
real(8),allocatable,dimension(:)   :: signg, nusigf, chi, sigalf, sigp, sign2n, sigd, sigt
real(8),allocatable,dimension(:)   :: phi, phi_old, source, xs_total, chi_tilde
real(8),allocatable,dimension(:,:) :: scatter, n2n
real(8) :: k, tol, converge, scat_sum, fiss_sum, scat_source_sum, lambda, lambda_old, numerator, denominator
real(8) :: chi_top, chi_bot
real(8) :: chi_tilde_sum

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
allocate(scatter(ngroup,ngroup))
allocate(n2n(ngroup,ngroup))
allocate(chi_tilde(ngroup))

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
		scatter(j,:) = scatter(j,:) + xs(i)%scat(j,:)
		n2n(j,:) = n2n(j,:) + xs(i)%n2n(j,:)
	enddo
enddo
! iterate
! intialize
do j = 1,ngroup
	scat_sum = 0.0d0
	do g = 1,ngroup
		scat_sum = scat_sum + scatter(g,j)
	enddo
	xs_total(j) = signg(j) + sigalf(j) + sigp(j) + sign2n(j) + scat_sum
enddo
phi = 1.0d0
phi_old = phi
k   = 1.0d0
lambda = 1.0d0
iteration = 0
tol = 1.0d-5
converge = 1.0d0
do while (converge .gt. tol)
	phi_old = phi
	do g = 1,ngroup
		fiss_sum = 0.0d0
		do gprime = 1,ngroup
			fiss_sum = fiss_sum + nusigf(gprime) * phi(gprime)
		enddo
		scat_source_sum = 0.0d0
		do gprime = 1,ngroup
			scat_source_sum = scat_source_sum + (2.0d0 * n2n(g,gprime) + scatter(g,gprime)) * phi(gprime)
		enddo
		chi_top = 0.0d0
		chi_bot = 0.0d0
		do i = 1,niso
			chi_top = chi_top + xs(i)%chi(g) * xs(i)%nuf(g) * xs(i)%sigf(g) * phi(g)
		enddo
		chi_bot = nusigf(g) * phi(g)
		chi_tilde(g) = chi_top / chi_bot
		chi_tilde_sum = 0.0d0
		do gprime = 1,ngroup
			chi_tilde_sum = chi_tilde_sum + chi_tilde(gprime)
		enddo
		! write(*,'(a,f12.10)') 'chi_tilde_sum = ', chi_tilde_sum
		do gprime = 1,ngroup
			chi_tilde(gprime) = chi_tilde(gprime) / chi_tilde_sum
		enddo
		chi_tilde_sum = 0.0d0
		do gprime = 1,ngroup
			chi_tilde_sum = chi_tilde_sum + chi_tilde(gprime)
		enddo
		write(*,'(a,f12.10)') 'chi_tilde_sum = ', chi_tilde_sum
		source(g) = (chi_tilde(g) / lambda) * fiss_sum + scat_source_sum
		phi(g) = source(g) / xs_total(g)
	enddo
	numerator = 0.0d0
	denominator = 0.0d0
	do gprime = 1,ngroup
		numerator = numerator + nusigf(gprime) * phi(gprime)
		denominator = denominator + nusigf(gprime) * phi_old(gprime)
	enddo
	write(20,*) numerator, denominator
	lambda_old = lambda
	lambda = lambda * (numerator / denominator)
	iteration = iteration + 1
	write(*,*) iteration, lambda
	! converge = 1.0d-10
	converge = abs(lambda - lambda_old) / lambda
	if (iteration .eq. 101) stop
enddo
write(*,'(a,i3)') 'iteration ', iteration
write(*,'(a,f12.10)') 'k ', k
write(*,'(a,f12.10)') 'lambda ', lambda

endsubroutine spectrum_solve

endmodule spectrum_calc