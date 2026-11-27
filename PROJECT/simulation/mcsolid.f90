! *****************************************************************************
! Monte Carlo code to simulate at NVT collectivity
! a system of N atoms of Ar inside a cubic box, interacting
! through a truncated Lennard-Jones potential (rc = min(5, L/2) sigma)
! plus analytic tail corrections.
! Metropolis algorithm: one trial displacement of a random atom per cycle.
!
! PROJECT version (solid Ar), built from the P3 code mcnvt.f90:
!   - initial configuration = perfect fcc lattice (N = 4 n^3) at density rho*
!     (no leap-conf.data needed)
!   - output files labelled by temperature (T scan at fixed rho*)
!   - position array allocated with the actual N
!   - structure factor S(k) on the k grid of the box, Bragg peaks S(G) of
!     the fcc lattice and translational order parameter Phi (111 peaks)
!
! Input : mc-solid.data
! Output: energy-mc_TX.XX.dat   (cycle, U*/N)
!         pressure-mc_TX.XX.dat (cycle, P*)
!         rdf-mc_TX.XX.dat      (r*, g(r))
!         sk-mc_TX.XX.dat       (k*, S(k), number of k vectors in the shell)
!         sk-lattice.dat        (same for the perfect fcc lattice, T* = 0)
!         bragg-mc_TX.XX.dat    (h2+k2+l2, G*, S(G)/N, sigma, error of mean, nvec)
!         order-mc_TX.XX.dat    (cycle, Phi)
!         mc-summary.dat        (one line per run: parameters + results)
!         newconf.data          (last configuration in A)
! *****************************************************************************

program mcsolid
    implicit none
    integer, parameter :: dp = kind(0.d0)
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
    integer, parameter :: nhis = 200           ! bins of the g(r) histogram
    integer, parameter :: nblocks = 10         ! blocks for the error of the mean

    ! ##################################################
    ! 1. Defining dimensions and declaring variables
    ! ##################################################
    integer :: ncycl, nsamp, nequil, natoms, iseed, nadjust
    integer :: icycl, is, l, k, nseed
    integer(8) :: nacc, natt
    real(dp) :: tref, sigma, epsil, deltax, rho, rc
    real(dp) :: boxlength, vol, utail, ptail, epot, accrate
    real(dp), allocatable :: r(:,:)

    ! g(r)
    real(dp) :: g(nhis), delg, rlo, rhi, shell

    ! sampling
    integer :: nsamples
    real(dp), allocatable :: useries(:), pseries(:)
    real(dp) :: mean_u, sig_u, err_u, mean_p, sig_p, err_p

    ! S(k), Bragg peaks S(G) and order parameter Phi
    !   k = dk * (nx,ny,nz), dk = 2 pi / L, half space only (S(-k) = S(k))
    !   kn2(ik) = nx^2+ny^2+nz^2          -> shell index of S(k)
    !   km2(ik) = h^2+k^2+l^2 if k = G_hkl -> Bragg family index (0 otherwise)
    integer :: nsampk, ncell, nmax, n2max, maxm2, nk, nsk, ik
    real(dp) :: kmax, dk, phi
    integer, allocatable :: kv(:,:), kn2(:), km2(:), nvec_n2(:), nvec_m2(:)
    real(dp), allocatable :: sk_now(:), skacc(:), sfam_now(:), sfam_series(:,:), phiseries(:)
    real(dp) :: mean_phi, sig_phi, err_phi, mean_s, sig_s, err_s

    integer, allocatable :: seed(:)
    character(len=64) :: fener, fpres, frdf, fsk, fbragg, forder
    logical :: exists

    ! ##################################################
    ! 2. Reading data and computing related quantities
    ! ##################################################
    open(1, file='mc-solid.data', status='old')
    read(1,*) ncycl, nsamp, nequil
    read(1,*) natoms
    read(1,*) tref
    read(1,*) sigma, epsil
    read(1,*) deltax
    read(1,*) rho
    read(1,*) iseed
    read(1,*) nsampk, kmax                  ! S(k) every nsampk cycles, for |k| <= kmax (1/sigma)
    close(1)

    ! Reproducible random numbers
    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed = iseed + 37 * [(k - 1, k = 1, nseed)]
    call random_seed(put=seed)

    if (tref <= 0.0_dp) stop 'ERROR: T* must be > 0'
    if (rho  <= 0.0_dp) stop 'ERROR: rho* must be > 0'
    if (mod(ncycl, nsampk) /= 0) stop 'ERROR: ncycl must be a multiple of nsampk'

    ! ##################################################
    ! 3-4. Initial configuration: perfect fcc lattice, reduced units
    ! ##################################################
    allocate(r(3, natoms))
    call fcc_lattice(natoms, rho, r, boxlength, ncell)

    ! ##################################################
    ! 4b. k vectors compatible with the periodic box, |k| <= kmax
    ! ##################################################
    dk = 2.0_dp * pi / boxlength
    nmax = int(kmax / dk)
    n2max = nmax**2
    maxm2 = n2max / ncell**2
    allocate(kv(3, (2*nmax + 1)**3), kn2((2*nmax + 1)**3), km2((2*nmax + 1)**3))
    call kvectors(nmax, ncell, size(kn2), kv, kn2, km2, nk)

    allocate(nvec_n2(n2max), nvec_m2(maxm2))
    nvec_n2 = 0
    nvec_m2 = 0
    do ik = 1, nk
        nvec_n2(kn2(ik)) = nvec_n2(kn2(ik)) + 1
        if (km2(ik) > 0) nvec_m2(km2(ik)) = nvec_m2(km2(ik)) + 1
    end do
    if (maxm2 < 3) stop 'ERROR: kmax too small to include the (111) Bragg peak'
    if (nvec_m2(3) /= 4) stop 'ERROR: expected 4 (111) vectors in the half space'

    allocate(sk_now(n2max), skacc(n2max), sfam_now(maxm2))
    allocate(sfam_series(ncycl / nsampk, maxm2), phiseries(ncycl / nsampk))
    skacc = 0.0_dp
    nsk = 0
    write(*,'(a,f8.4,a,i4,a,f8.4,a,i8)') ' S(k): dk* =', dk, '   nmax =', nmax, &
                                         '   kmax* =', nmax * dk, '   k vectors =', nk

    ! S(k) of the perfect lattice (reference + check: S(G)/N = 1, Phi = 1)
    call sample_sk(natoms, r, nmax, dk, nk, kv, kn2, km2, n2max, nvec_n2, maxm2, nvec_m2, &
                   sk_now, sfam_now, phi)
    write(*,'(a,f10.6,a,f10.6,a,f10.6)') ' perfect fcc: Phi =', phi, '   S(111)/N =', &
                                         sfam_now(3) / natoms, '   S(200)/N =', sfam_now(4) / natoms
    call write_sk('sk-lattice.dat', n2max, nvec_n2, dk, sk_now, 1)

    vol = boxlength**3
    rc  = min(5.0_dp, 0.5_dp * boxlength)
    utail = (8.0_dp/3.0_dp) * pi * rho * ((1.0_dp/3.0_dp) / rc**9 - 1.0_dp / rc**3)          ! per atom
    ptail = (16.0_dp/3.0_dp) * pi * rho**2 * ((2.0_dp/3.0_dp) / rc**9 - 1.0_dp / rc**3)

    delg = 0.5_dp * boxlength / real(nhis, dp)
    g = 0.0_dp
    allocate(useries(ncycl / nsamp), pseries(ncycl / nsamp))
    nsamples = 0

    nadjust = 100 * natoms                  ! adjust deltax every 100 sweeps (equilibration only)

    call energy(natoms, r, boxlength, rc, epot)

    write(*,'(a,f8.4,a,f8.4,a,f8.4,a,f8.4)') ' T* =', tref, '   rho* =', rho, &
                                             '   L* =', boxlength, '   rc* =', rc
    ! static energy of the perfect lattice (check against the analytic fcc lattice sum)
    write(*,'(a,f12.6,a,f12.6)') ' U*/N perfect fcc (truncated) =', epot / real(natoms, dp), &
                                 '   + tail =', epot / real(natoms, dp) + utail
    write(*,'(a,i12,a,i12,a,i8)') ' nequil =', nequil, '   ncycl =', ncycl, '   nsamp =', nsamp

    ! Opening files to write results
    write(fener, '(a,f4.2,a)') 'energy-mc_T',   tref, '.dat'
    write(fpres, '(a,f4.2,a)') 'pressure-mc_T', tref, '.dat'
    write(frdf,  '(a,f4.2,a)') 'rdf-mc_T',      tref, '.dat'
    write(fsk,   '(a,f4.2,a)') 'sk-mc_T',       tref, '.dat'
    write(fbragg,'(a,f4.2,a)') 'bragg-mc_T',    tref, '.dat'
    write(forder,'(a,f4.2,a)') 'order-mc_T',    tref, '.dat'
    open(3, file=trim(fener), status='replace')
    open(4, file=trim(fpres), status='replace')
    open(7, file=trim(forder), status='replace')

    ! ################################################
    ! 5. Start the loop to generate new configurations
    ! ################################################
    nacc = 0
    natt = 0
    do icycl = 1, nequil + ncycl
        call mcmove(natoms, r, boxlength, rc, tref, deltax, epot, nacc)
        natt = natt + 1

        if (icycl <= nequil) then
            ! 5.1 Equilibration: tune deltax towards 50 % acceptance, no sampling
            if (mod(icycl, nadjust) == 0) then
                accrate = real(nacc, dp) / real(natt, dp)
                deltax = deltax * min(max(accrate / 0.5_dp, 0.5_dp), 1.5_dp)
                deltax = min(max(deltax, 0.01_dp), 0.5_dp * boxlength)
                nacc = 0
                natt = 0
            end if
            if (icycl == nequil) then
                nacc = 0
                natt = 0
            end if
        else
            ! 5.2 Production
            if (mod(icycl - nequil, nsamp) == 0) then
                call sample(natoms, r, boxlength, rc, rho, tref, utail, ptail, delg, nhis, g, &
                            epot, icycl - nequil, nsamples, useries, pseries)
            end if
            if (mod(icycl - nequil, nsampk) == 0) then
                call sample_sk(natoms, r, nmax, dk, nk, kv, kn2, km2, n2max, nvec_n2, maxm2, &
                               nvec_m2, sk_now, sfam_now, phi)
                nsk = nsk + 1
                skacc = skacc + sk_now
                sfam_series(nsk, :) = sfam_now / real(natoms, dp)     ! S(G)/N per family
                phiseries(nsk) = phi
                write(7,*) icycl - nequil, phi
            end if
        end if
    end do
    close(3)
    close(4)
    close(7)

    accrate = real(nacc, dp) / real(natt, dp)

    ! ##################
    ! 6. Result analysis
    ! ##################
    ! 6.1 Averages, standard deviations and block errors of the mean
    call stats(nsamples, useries, nblocks, mean_u, sig_u, err_u)
    call stats(nsamples, pseries, nblocks, mean_p, sig_p, err_p)
    call stats(nsk, phiseries, nblocks, mean_phi, sig_phi, err_phi)

    write(*,'(a,f10.5,a,f8.4)') ' deltax* =', deltax, '   acceptance =', accrate
    write(*,'(a,i10,a,i10)')    ' samples =', nsamples, '   S(k) samples =', nsk
    write(*,'(a,3f12.5)')       ' U*/N, sigma_U, err =', mean_u, sig_u, err_u
    write(*,'(a,3f12.5)')       ' P*  , sigma_P, err =', mean_p, sig_p, err_p
    write(*,'(a,3f12.5)')       ' Phi , sigma  , err =', mean_phi, sig_phi, err_phi

    inquire(file='mc-summary.dat', exist=exists)
    open(9, file='mc-summary.dat', status='unknown', position='append')
    if (.not. exists) then
        write(9,'(a)') '#  T*     rho*      N        nequil         ncycl   nsamp  deltax*  acc     rc*     ' // &
                       'U*/N        sigma_U     err_U       P*          sigma_P     err_P       ' // &
                       'Phi         sigma_Phi   err_Phi'
    end if
    write(9,'(2f8.4,i6,2i14,i8,3f8.4,9f12.6)') tref, rho, natoms, nequil, ncycl, nsamp, deltax, &
                                               accrate, rc, mean_u, sig_u, err_u, mean_p, sig_p, err_p, &
                                               mean_phi, sig_phi, err_phi
    close(9)

    ! 6.1b S(k) averaged over the run, and Bragg peaks S(G)/N per {hkl} family
    call write_sk(trim(fsk), n2max, nvec_n2, dk, skacc, nsk)
    open(8, file=trim(fbragg), status='replace')
    write(8,'(a)') '# h2+k2+l2   G*            S(G)/N        sigma         err_mean      nvec'
    do ik = 1, maxm2
        if (nvec_m2(ik) == 0) cycle
        call stats(nsk, sfam_series(1:nsk, ik), nblocks, mean_s, sig_s, err_s)
        write(8,'(i8,4es14.6,i8)') ik, dk * ncell * sqrt(real(ik, dp)), mean_s, sig_s, err_s, nvec_m2(ik)
    end do
    close(8)

    ! 6.2 Normalize g(r): pairs counted twice -> divide by N * rho * shell volume
    open(5, file=trim(frdf), status='replace')
    do k = 1, nhis
        rlo = real(k - 1, dp) * delg
        rhi = real(k, dp) * delg
        shell = (4.0_dp/3.0_dp) * pi * (rhi**3 - rlo**3)
        write(5,'(2f14.6)') 0.5_dp * (rlo + rhi), &
                            g(k) / (real(nsamples, dp) * real(natoms, dp) * rho * shell)
    end do
    close(5)

    ! 6.3 Saving last configuration in A (no velocities)
    open(11, file='newconf.data', status='unknown')
    do is = 1, natoms
        write(11,*) (r(l,is)*sigma, l=1,3)
    end do
    write(11,*) boxlength*sigma
    close(11)

end program mcsolid

! *********************************************************
! subroutine fcc_lattice: perfect fcc crystal in a cubic box (reduced units)
!   fcc = 4 interpenetrating simple cubic lattices of side a, shifted by
!   the basis vectors (0,0,0), (1/2,1/2,0), (1/2,0,1/2), (0,1/2,1/2) * a
!   N = 4 ncell^3 atoms, 4 atoms per cubic cell -> rho* = 4 / a^3
! *********************************************************
subroutine fcc_lattice(natoms, rho, r, boxlength, ncell)
    implicit none
    integer, parameter :: dp = kind(0.d0)
    integer, intent(in) :: natoms
    real(dp), intent(in) :: rho
    real(dp), intent(out) :: r(3,natoms), boxlength
    integer, intent(out) :: ncell

    real(dp), parameter :: basis(3,4) = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                                 0.5_dp, 0.5_dp, 0.0_dp, &
                                                 0.5_dp, 0.0_dp, 0.5_dp, &
                                                 0.0_dp, 0.5_dp, 0.5_dp], [3,4])
    real(dp) :: a
    integer :: ix, iy, iz, ib, i

    ! the box must contain an integer number of cubic cells, otherwise the
    ! crystal does not match itself through the periodic boundaries
    ncell = nint((real(natoms, dp) / 4.0_dp)**(1.0_dp/3.0_dp))
    if (4 * ncell**3 /= natoms) stop 'ERROR: fcc lattice needs N = 4 n^3 (108, 256, 500, 864, ...)'

    a = (4.0_dp / rho)**(1.0_dp/3.0_dp)     ! lattice constant (sigma)
    boxlength = real(ncell, dp) * a

    i = 0
    do ix = 0, ncell - 1
        do iy = 0, ncell - 1
            do iz = 0, ncell - 1
                do ib = 1, 4
                    i = i + 1
                    r(:,i) = a * (real([ix, iy, iz], dp) + basis(:,ib))
                end do
            end do
        end do
    end do

    write(*,'(a,i4,a,f9.5,a,f9.5)') ' fcc lattice: ncell =', ncell, '   a* =', a, &
                                    '   r_nn* = a/sqrt(2) =', a / sqrt(2.0_dp)
end subroutine fcc_lattice

! *********************************************************
! subroutine kvectors: list of wave vectors k = (2 pi / L) (nx,ny,nz)
!   - only the half space (nz>0, or nz=0 & ny>0, or nz=ny=0 & nx>0),
!     since S(-k) = S(k) for real positions
!   - |n|^2 <= nmax^2 (sphere), k = 0 excluded
!   - Bragg vectors of the fcc crystal: G = (2 pi / a)(h,k,l) with h,k,l all odd
!     or all even; since L = ncell*a, n = ncell*(h,k,l). km2 = h^2+k^2+l^2
! *********************************************************
subroutine kvectors(nmax, ncell, nkdim, kv, kn2, km2, nk)
    implicit none
    integer, intent(in) :: nmax, ncell, nkdim
    integer, intent(out) :: kv(3,nkdim), kn2(nkdim), km2(nkdim), nk

    integer :: nx, ny, nz, n2, hkl(3)

    nk = 0
    do nz = 0, nmax
        do ny = -nmax, nmax
            do nx = -nmax, nmax
                n2 = nx**2 + ny**2 + nz**2
                if (n2 == 0 .or. n2 > nmax**2) cycle
                if (nz == 0 .and. (ny < 0 .or. (ny == 0 .and. nx < 0))) cycle
                nk = nk + 1
                if (nk > nkdim) stop 'ERROR: kvectors, nkdim too small'
                kv(:,nk) = [nx, ny, nz]
                kn2(nk) = n2
                km2(nk) = 0
                if (all(mod(kv(:,nk), ncell) == 0)) then
                    hkl = kv(:,nk) / ncell
                    if (all(modulo(hkl, 2) == modulo(hkl(1), 2))) km2(nk) = sum(hkl**2)
                end if
            end do
        end do
    end do
end subroutine kvectors

! *********************************************************
! subroutine sample_sk: structure factor of the current configuration
!   rho(k) = sum_j exp(i k.r_j),   S(k) = |rho(k)|^2 / N
!   sk_now(n2)  = S averaged over the k vectors of the shell |n|^2 = n2
!   sfam_now(m) = S averaged over the Bragg vectors of family h2+k2+l2 = m
!   phi         = |rho(G)| / N averaged over the 4 (111) vectors (order parameter)
!   exp(i k.r) = ex(nx) * ey(ny) * ez(nz), with the 1D factors computed once per atom
! *********************************************************
subroutine sample_sk(natoms, r, nmax, dk, nk, kv, kn2, km2, n2max, nvec_n2, maxm2, nvec_m2, &
                     sk_now, sfam_now, phi)
    implicit none
    integer, parameter :: dp = kind(0.d0)
    integer, intent(in) :: natoms, nmax, nk, n2max, maxm2
    integer, intent(in) :: kv(3,nk), kn2(nk), km2(nk), nvec_n2(n2max), nvec_m2(maxm2)
    real(dp), intent(in) :: r(3,natoms), dk
    real(dp), intent(out) :: sk_now(n2max), sfam_now(maxm2), phi

    complex(dp), allocatable :: ex(:,:), ey(:,:), ez(:,:)
    complex(dp) :: rhok
    real(dp) :: s
    integer :: n, ik

    allocate(ex(natoms,-nmax:nmax), ey(natoms,-nmax:nmax), ez(natoms,-nmax:nmax))
    do n = -nmax, nmax
        ex(:,n) = cmplx(cos(n * dk * r(1,:)), sin(n * dk * r(1,:)), dp)
        ey(:,n) = cmplx(cos(n * dk * r(2,:)), sin(n * dk * r(2,:)), dp)
        ez(:,n) = cmplx(cos(n * dk * r(3,:)), sin(n * dk * r(3,:)), dp)
    end do

    sk_now = 0.0_dp
    sfam_now = 0.0_dp
    phi = 0.0_dp
    do ik = 1, nk
        rhok = sum(ex(:,kv(1,ik)) * ey(:,kv(2,ik)) * ez(:,kv(3,ik)))
        s = (real(rhok, dp)**2 + aimag(rhok)**2) / real(natoms, dp)
        sk_now(kn2(ik)) = sk_now(kn2(ik)) + s
        if (km2(ik) > 0) sfam_now(km2(ik)) = sfam_now(km2(ik)) + s
        if (km2(ik) == 3) phi = phi + abs(rhok) / real(natoms, dp)
    end do

    where (nvec_n2 > 0) sk_now = sk_now / real(nvec_n2, dp)
    where (nvec_m2 > 0) sfam_now = sfam_now / real(nvec_m2, dp)
    phi = phi / real(nvec_m2(3), dp)
end subroutine sample_sk

! *********************************************************
! subroutine write_sk: k*, S(k), number of k vectors in the shell
!   skacc = sum over nsk samples of the shell averages
! *********************************************************
subroutine write_sk(fname, n2max, nvec_n2, dk, skacc, nsk)
    implicit none
    integer, parameter :: dp = kind(0.d0)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: n2max, nvec_n2(n2max), nsk
    real(dp), intent(in) :: dk, skacc(n2max)

    integer :: n2

    if (nsk < 1) stop 'ERROR: write_sk with no S(k) samples'
    open(12, file=fname, status='replace')
    write(12,'(a)') '#  k*            S(k)          nvec'
    do n2 = 1, n2max
        if (nvec_n2(n2) == 0) cycle
        write(12,'(2es14.6,i8)') dk * sqrt(real(n2, dp)), skacc(n2) / real(nsk, dp), nvec_n2(n2)
    end do
    close(12)
end subroutine write_sk

! *********************************************************
! subroutine energy (was forces): total potential energy, no accelerations
! *********************************************************
subroutine energy(natoms, r, boxlength, rc, epot)
    implicit none
    integer, parameter :: dp = kind(0.d0)
    integer, intent(in) :: natoms
    real(dp), intent(in) :: r(3,natoms), boxlength, rc
    real(dp), intent(out) :: epot

    real(dp) :: pot, vir, rr
    integer :: is, js

    epot = 0.0_dp
    do is = 1, natoms - 1
        do js = is + 1, natoms
            call lj(r(:,is), r(:,js), boxlength, rc, pot, vir, rr)
            epot = epot + pot
        end do
    end do
end subroutine energy

! *********************************************************
! subroutine eneri: energy of atom is placed at ri with all the others
! *********************************************************
subroutine eneri(is, ri, natoms, r, boxlength, rc, ei)
    implicit none
    integer, parameter :: dp = kind(0.d0)
    integer, intent(in) :: is, natoms
    real(dp), intent(in) :: ri(3), r(3,natoms), boxlength, rc
    real(dp), intent(out) :: ei

    real(dp) :: pot, vir, rr
    integer :: js

    ei = 0.0_dp
    do js = 1, natoms
        if (js == is) cycle
        call lj(ri, r(:,js), boxlength, rc, pot, vir, rr)
        ei = ei + pot
    end do
end subroutine eneri

! *********************************************************
! subroutine Lennard-Jones (simplified: energy and virial, no accelerations)
!   pot = 4 (1/r^12 - 1/r^6),   vir = r_ij . F_ij = 24 (2/r^12 - 1/r^6)
!   rr  = minimum-image distance (used for g(r))
! *********************************************************
subroutine lj(ri, rj, boxlength, rc, pot, vir, rr)
    implicit none
    integer, parameter :: dp = kind(0.d0)
    real(dp), intent(in) :: ri(3), rj(3), boxlength, rc
    real(dp), intent(out) :: pot, vir, rr

    real(dp) :: rij(3), rr2, ynvrr2, ynvrr6, ynvrr12
    integer :: l

    rr2 = 0.0_dp
    pot = 0.0_dp
    vir = 0.0_dp
    do l = 1, 3
        rij(l) = rj(l) - ri(l)
        rij(l) = rij(l) - boxlength * anint(rij(l) / boxlength)
        rr2 = rr2 + rij(l) * rij(l)
    end do
    rr = sqrt(rr2)

    if (rr < rc) then
        ynvrr2 = 1.0_dp / rr2
        ynvrr6 = ynvrr2 * ynvrr2 * ynvrr2
        ynvrr12 = ynvrr6 * ynvrr6
        pot = 4.0_dp * (ynvrr12 - ynvrr6)
        vir = 24.0_dp * (2.0_dp * ynvrr12 - ynvrr6)
    end if
end subroutine lj

! *********************************************************
! subroutine mcmove (replaces velpos): one Metropolis trial displacement
! *********************************************************
subroutine mcmove(natoms, r, boxlength, rc, tref, deltax, epot, nacc)
    implicit none
    integer, parameter :: dp = kind(0.d0)
    integer, intent(in) :: natoms
    real(dp), intent(inout) :: r(3,natoms), epot
    real(dp), intent(in) :: boxlength, rc, tref, deltax
    integer(8), intent(inout) :: nacc

    integer :: is, l
    real(dp) :: rnd, xi(3), rnew(3), eold, enew, du

    ! choose an atom at random
    call random_number(rnd)
    is = min(int(rnd * natoms) + 1, natoms)

    call eneri(is, r(:,is), natoms, r, boxlength, rc, eold)

    ! random displacement in [-deltax/2, deltax/2] and periodic boundary conditions
    call random_number(xi)
    do l = 1, 3
        rnew(l) = r(l,is) + deltax * (xi(l) - 0.5_dp)
        if (rnew(l) < 0.0_dp) rnew(l) = rnew(l) + boxlength
        if (rnew(l) >= boxlength) rnew(l) = rnew(l) - boxlength
    end do

    call eneri(is, rnew, natoms, r, boxlength, rc, enew)
    du = enew - eold

    ! Metropolis criterion with beta* = 1/T*
    ! (if rejected, the old configuration is kept and counted again)
    call random_number(rnd)
    if (du <= 0.0_dp .or. rnd < exp(-du / tref)) then
        r(:,is) = rnew
        epot = epot + du
        nacc = nacc + 1
    end if
end subroutine mcmove

! *********************************************************
! subroutine sample: U, P (virial) and g(r) histogram of the current configuration
! *********************************************************
subroutine sample(natoms, r, boxlength, rc, rho, tref, utail, ptail, delg, nhis, g, &
                  epot, isamp, nsamples, useries, pseries)
    implicit none
    integer, parameter :: dp = kind(0.d0)
    integer, intent(in) :: natoms, nhis, isamp
    real(dp), intent(in) :: r(3,natoms), boxlength, rc, rho, tref, utail, ptail, delg
    real(dp), intent(inout) :: g(nhis), epot, useries(*), pseries(*)
    integer, intent(inout) :: nsamples

    integer :: is, js, ig
    real(dp) :: pot, vir, rr, utot, virial, un, p

    utot = 0.0_dp
    virial = 0.0_dp
    do is = 1, natoms - 1
        do js = is + 1, natoms
            call lj(r(:,is), r(:,js), boxlength, rc, pot, vir, rr)
            utot = utot + pot
            virial = virial + vir
            if (rr < 0.5_dp * boxlength) then
                ig = int(rr / delg) + 1
                if (ig <= nhis) g(ig) = g(ig) + 2.0_dp
            end if
        end do
    end do

    ! check of the running energy kept by mcmove
    if (abs(utot - epot) > 1.0e-6_dp * max(1.0_dp, abs(utot))) then
        write(*,*) 'WARNING: energy drift', isamp, epot, utot
    end if
    epot = utot

    un = utot / real(natoms, dp) + utail
    p  = rho * tref + virial / (3.0_dp * boxlength**3) + ptail

    nsamples = nsamples + 1
    useries(nsamples) = un
    pseries(nsamples) = p
    write(3,*) isamp, un
    write(4,*) isamp, p
end subroutine sample

! *********************************************************
! subroutine pbc: put all atoms inside [0, L)
! *********************************************************
subroutine pbc(natoms, r, boxlength)
    implicit none
    integer, parameter :: dp = kind(0.d0)
    integer, intent(in) :: natoms
    real(dp), intent(inout) :: r(3,natoms)
    real(dp), intent(in) :: boxlength

    r(:,1:natoms) = r(:,1:natoms) - boxlength * floor(r(:,1:natoms) / boxlength)
end subroutine pbc

! *********************************************************
! subroutine stats: mean, standard deviation of the samples and
! error of the mean from nb block averages
! *********************************************************
subroutine stats(n, x, nb, mean, sig, err)
    implicit none
    integer, parameter :: dp = kind(0.d0)
    integer, intent(in) :: n, nb
    real(dp), intent(in) :: x(n)
    real(dp), intent(out) :: mean, sig, err

    integer :: b, lb
    real(dp) :: bm, s2

    mean = sum(x) / real(n, dp)
    sig  = sqrt(sum((x - mean)**2) / real(n - 1, dp))

    lb = n / nb
    if (lb < 1) then
        err = -1.0_dp
        return
    end if
    s2 = 0.0_dp
    do b = 1, nb
        bm = sum(x((b-1)*lb + 1 : b*lb)) / real(lb, dp)
        s2 = s2 + (bm - mean)**2
    end do
    err = sqrt(s2 / real(nb * (nb - 1), dp))
end subroutine stats
