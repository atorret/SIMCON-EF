! *****************************************************************************
! Monte Carlo code to simulate at NVT collectivity
! a system of N atoms of Ar inside a cubic box, interacting
! through a truncated Lennard-Jones potential (rc = min(5, L/2) sigma)
! plus analytic tail corrections.
! Metropolis algorithm: one trial displacement of a random atom per cycle.
!
! Built from leapfroglj.f following the proposed modifications:
!   - no velocities / accelerations
!   - forces + lj  -> energy-only routines (plus virial for the pressure)
!   - velpos       -> mcmove
!   - every nsamp cycles -> sample (U, P, g(r))
!
! Input : mc-lj.data, leap-conf.data (positions in A; velocities are skipped)
! Output: energy-mc_rhoX.XX.dat   (cycle, U*/N)
!         pressure-mc_rhoX.XX.dat (cycle, P*)
!         rdf-mc_rhoX.XX.dat      (r*, g(r))
!         mc-summary.dat          (one line per run: parameters + results)
!         newconf.data            (last configuration in A)
! *****************************************************************************

program mcnvt
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
    real(dp) :: tref, sigma, epsil, deltax, rho, rc, rhoconf
    real(dp) :: boxlength, vol, fac, utail, ptail, epot, accrate
    real(dp) :: r(3,1000), vdummy(3)

    ! g(r)
    real(dp) :: g(nhis), delg, rlo, rhi, shell

    ! sampling
    integer :: nsamples
    real(dp), allocatable :: useries(:), pseries(:)
    real(dp) :: mean_u, sig_u, err_u, mean_p, sig_p, err_p

    integer, allocatable :: seed(:)
    character(len=64) :: fener, fpres, frdf
    logical :: exists

    ! ##################################################
    ! 2. Reading data and computing related quantities
    ! ##################################################
    open(1, file='mc-lj.data', status='old')
    read(1,*) ncycl, nsamp, nequil
    read(1,*) natoms
    read(1,*) tref
    read(1,*) sigma, epsil
    read(1,*) deltax
    read(1,*) rho
    read(1,*) iseed
    close(1)

    ! Reproducible random numbers
    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed = iseed + 37 * [(k - 1, k = 1, nseed)]
    call random_seed(put=seed)

    ! ##################################################
    ! 3. Reading initial configuration (positions) in A
    ! ##################################################
    open(2, file='leap-conf.data', status='old')
    do is = 1, natoms
        read(2,*) (r(l,is), l=1,3)
        read(2,*) (vdummy(l), l=1,3)        ! velocities are not used in MC
    end do
    read(2,*) boxlength
    close(2)

    ! ##################################################
    ! 4. Change to reduced units and set the density
    ! ##################################################
    call reduced(natoms, r, boxlength, sigma)

    rhoconf = real(natoms, dp) / boxlength**3
    if (rho <= 0.0_dp) rho = rhoconf        ! keep the density of leap-conf.data
    fac = (rhoconf / rho)**(1.0_dp/3.0_dp)  ! scale box and positions to the target rho
    boxlength = boxlength * fac
    r(:,1:natoms) = r(:,1:natoms) * fac
    call pbc(natoms, r, boxlength)

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
    write(*,'(a,i12,a,i12,a,i8)') ' nequil =', nequil, '   ncycl =', ncycl, '   nsamp =', nsamp

    ! Opening files to write results
    write(fener, '(a,f4.2,a)') 'energy-mc_rho',   rho, '.dat'
    write(fpres, '(a,f4.2,a)') 'pressure-mc_rho', rho, '.dat'
    write(frdf,  '(a,f4.2,a)') 'rdf-mc_rho',      rho, '.dat'
    open(3, file=trim(fener), status='replace')
    open(4, file=trim(fpres), status='replace')

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
        end if
    end do
    close(3)
    close(4)

    accrate = real(nacc, dp) / real(natt, dp)

    ! ##################
    ! 6. Result analysis
    ! ##################
    ! 6.1 Averages, standard deviations and block errors of the mean
    call stats(nsamples, useries, nblocks, mean_u, sig_u, err_u)
    call stats(nsamples, pseries, nblocks, mean_p, sig_p, err_p)

    write(*,'(a,f10.5,a,f8.4)') ' deltax* =', deltax, '   acceptance =', accrate
    write(*,'(a,i10)')          ' samples =', nsamples
    write(*,'(a,3f12.5)')       ' U*/N, sigma_U, err =', mean_u, sig_u, err_u
    write(*,'(a,3f12.5)')       ' P*  , sigma_P, err =', mean_p, sig_p, err_p

    inquire(file='mc-summary.dat', exist=exists)
    open(9, file='mc-summary.dat', status='unknown', position='append')
    if (.not. exists) then
        write(9,'(a)') '#  rho*     T*      N        nequil         ncycl   nsamp  deltax*  acc     rc*     ' // &
                       'U*/N        sigma_U     err_U       P*          sigma_P     err_P'
    end if
    write(9,'(2f8.4,i6,2i14,i8,3f8.4,6f12.6)') rho, tref, natoms, nequil, ncycl, nsamp, deltax, &
                                               accrate, rc, mean_u, sig_u, err_u, mean_p, sig_p, err_p
    close(9)

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

end program mcnvt

! *********************************************************
! subroutine reduced (simplified: no velocities)
! *********************************************************
subroutine reduced(natoms, r, boxlength, sigma)
    implicit none
    integer, parameter :: dp = kind(0.d0)
    integer, intent(in) :: natoms
    real(dp), intent(inout) :: r(3,1000), boxlength
    real(dp), intent(in) :: sigma

    boxlength = boxlength / sigma
    r(:,1:natoms) = r(:,1:natoms) / sigma
end subroutine reduced

! *********************************************************
! subroutine energy (was forces): total potential energy, no accelerations
! *********************************************************
subroutine energy(natoms, r, boxlength, rc, epot)
    implicit none
    integer, parameter :: dp = kind(0.d0)
    integer, intent(in) :: natoms
    real(dp), intent(in) :: r(3,1000), boxlength, rc
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
    real(dp), intent(in) :: ri(3), r(3,1000), boxlength, rc
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
    real(dp), intent(inout) :: r(3,1000), epot
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
    real(dp), intent(in) :: r(3,1000), boxlength, rc, rho, tref, utail, ptail, delg
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
    real(dp), intent(inout) :: r(3,1000)
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
