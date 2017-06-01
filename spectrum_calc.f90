module spectrum_calc
use variables

contains

subroutine spectrum_solve
IMPLICIT NONE
integer :: i, j, g, gprime
integer :: iteration
! real(8),allocatable,dimension(:,:) :: sigtr, sigtot
real(8),allocatable,dimension(:)   :: signg, sigf, nusigf, sigalf, sigp, sign2n, sigd, sigt, mpact_absnusigf
real(8),allocatable,dimension(:)   :: phi, phi_old, source, xs_total, chi_tilde
real(8),allocatable,dimension(:,:) :: scatter, n2n
real(8) :: tol, converge, scat_sum, fiss_sum, scat_source_sum, lambda, lambda_old, numerator, denominator
real(8) :: chi_top, chi_bot
real(8) :: chi_tilde_sum, iso_fiss_sum

allocate(signg(ngroup))
allocate(sigf(ngroup))
allocate(nusigf(ngroup))
allocate(sigalf(ngroup))
allocate(sigp(ngroup))
allocate(sign2n(ngroup))
allocate(sigd(ngroup))
allocate(sigt(ngroup))
allocate(scatter(ngroup,ngroup))
allocate(n2n(ngroup,ngroup))
allocate(xs_total(ngroup))
allocate(phi(ngroup))
allocate(phi_old(ngroup))
allocate(chi_tilde(ngroup))
allocate(source(ngroup))

! homogenize
! sum of macro xs
signg(:)  = 0.0d0
sigf(:)   = 0.0d0
nusigf(:) = 0.0d0
sigalf(:) = 0.0d0
sigp(:)   = 0.0d0
sign2n(:) = 0.0d0
sigd(:)   = 0.0d0
sigt(:)   = 0.0d0
scatter(:,:) = 0.0d0
n2n(:,:) = 0.0d0
do g = 1,ngroup
	do i = 1,niso
		signg(g)  = signg(g)  + xs(i)%signg(g)
		sigf(g)   = sigf(g)   + xs(i)%sigf(g)
		nusigf(g) = nusigf(g) + (xs(i)%nuf(g) * xs(i)%sigf(g))
		sigalf(g) = sigalf(g) + xs(i)%sigalf(g)
		sigp(g)   = sigp(g)   + xs(i)%sigp(g)
		sign2n(g) = sign2n(g) + xs(i)%sign2n(g)
		sigd(g)   = sigd(g)   + xs(i)%sigd(g)
		sigt(g)   = sigt(g)   + xs(i)%sigt(g)
		scatter(g,:) = scatter(g,:) + xs(i)%scat(g,:)
		n2n(g,:) = n2n(g,:) + xs(i)%n2n(g,:)
	enddo
enddo

! iterate
! intialize
do g = 1,ngroup
	scat_sum = 0.0d0
	do gprime = 1,ngroup
		scat_sum = scat_sum + scatter(g,gprime)
	enddo
	xs_total(g) = signg(g) + sigalf(g) + sigp(g) + sigf(g) + sign2n(g) + sigd(g) + sigt(g) + scat_sum
enddo
phi = 1.0d0
lambda = 1.0d0
iteration = 0
tol = 1.0d-7
converge = 1.0d0
do while (converge .gt. tol)
	phi_old = phi
	lambda_old = lambda
	
	
	! chi collapse denominator
	chi_bot = 0.0d0
	do gprime = 1,ngroup
		chi_bot = chi_bot + nusigf(gprime) * phi_old(gprime)
	enddo
	
	! fission sum in source term
	fiss_sum = 0.0d0
	do gprime = 1,ngroup
		fiss_sum = fiss_sum + nusigf(gprime) * phi_old(gprime)
	enddo
	
	
	do g = 1,ngroup
		chi_top = 0.0d0
		do i = 1,niso
			! isotopic fission sum in chi collapse numerator
			iso_fiss_sum = 0.0d0
			do gprime = 1,ngroup
				iso_fiss_sum = iso_fiss_sum + xs(i)%nuf(gprime) * xs(i)%sigf(gprime) * phi(gprime)
			enddo
			chi_top = chi_top + iso_fiss_sum * xs(i)%chi(g)
		enddo
		chi_tilde(g) = chi_top / chi_bot
	enddo
	chi_tilde_sum = vector_sum(chi_tilde,ngroup)
	
	
	
	
	do g = 1,ngroup
		! satter sum in source term
		scat_source_sum = 0.0d0
		do gprime = 1,ngroup
			scat_source_sum = scat_source_sum + (2.0d0 * n2n(gprime,g) + scatter(gprime,g)) * phi(gprime)
		enddo
		source(g) = (chi_tilde(g) / lambda) * fiss_sum + scat_source_sum
		phi(g) = source(g) / xs_total(g)
	enddo
	numerator = 0.0d0
	denominator = 0.0d0
	do gprime = 1,ngroup
		numerator = numerator + nusigf(gprime) * phi(gprime)
		denominator = denominator + nusigf(gprime) * phi_old(gprime)
	enddo
	lambda = lambda_old * (numerator / denominator)
	iteration = iteration + 1
	write(*,'(i3,x,f12.10,x,f12.10)') iteration, lambda, chi_tilde_sum
	converge = abs(lambda - lambda_old) / lambda
	if (iteration .eq. 1001) stop 'MAX ITERATION'
enddo
write(*,'(a,i3)') 'iteration ', iteration
write(*,'(a,f12.10)') 'lambda ', lambda
write(*,'(a)') '---PHI---'
do g = 1,ngroup
	write(*,'(i3,x,e12.6)') g, phi(g)
enddo

endsubroutine spectrum_solve

function vector_sum(var,length)
IMPLICIT NONE
	real(8) :: vector_sum
	real(8),dimension(:),allocatable,intent(in) :: var
	integer,intent(in) :: length
	integer :: i
	vector_sum = 0.0d0
	do i = 1,length
		vector_sum = vector_sum + var(i)
	enddo
endfunction vector_sum

endmodule spectrum_calc