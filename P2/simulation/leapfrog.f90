! *****************************************************************************
! Molecular Dynamics code to simulate at NVE collectivity 
! a system of N atoms of Ar inside a cubic box, interacting
! through a truncated Lennard-Jones force field at 2.5 sigma. 
! Leap-frog Verlet integration algorithm.
! *****************************************************************************

program leapfroglj
    implicit none
    integer, parameter :: dp = kind(0.d0)

    ! 0. Constants
    integer, parameter :: nhis = 500
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

    ! ##################################################
    ! 1. Defining dimensions and declaring ALL variables
    ! ##################################################

    ! 1.1 Setup variables
    integer :: nconf, natoms, nf, is, l, i
    real(dp) :: mass, sigma, epsil, deltat, rc
    real(dp) :: boxlength, uvel, epot, ecin, temp, etot
    real(dp) :: r(3,1000), vinf(3,1000), accel(3,1000)
    
    ! 1.2 g(r) variables
    real(dp) :: g(nhis), rr
    real(dp) :: delg, vb, nid, rho
   
    ! 1.3 Pressure variables
    real(dp) :: virial, p 

    ! 1.4 VACF and MSD variables
    integer, parameter :: tmax = 1001
    integer, parameter :: it0 = 50
    integer, parameter :: t0max = 25
    integer :: nsamp, ntel, t0, tt0, delt, t
    real(dp) :: dtime, real_time

    real(dp) :: ntime(tmax), vacf(tmax), r2t(tmax)
    integer :: time0(t0max)
    real(dp) :: r0(3, 1000, t0max), v0(3, 1000, t0max)
    real(dp) :: runw(3, 1000)

    ! 1.5 Spectrum parameters
    real(dp) :: time_array(tmax)
    integer, parameter :: num_w = 500
    real(dp) :: omega_array(0:num_w)
    real(dp) :: spectrum(0:num_w)
    integer :: w_idx

    ! ################################################
    ! 2. Reading data and computing related quantities
    ! ################################################
    open(1, file='leap-lj.data', status='old')
    read(1,*) nconf
    read(1,*) natoms
    read(1,*) mass
    read(1,*) sigma, epsil
    read(1,*) deltat 
    close(1)
    
    nf = 3 * natoms - 3 
    rc = 2.5_dp

    ! ######################################################################
    ! 3. Reading initial configuration (positions, velocities) in A and A/ps
    ! ######################################################################
    ! 3.1 Reading the initial configuration from the file 'leap-conf.data'
    open(2, file='leap-conf.data', status='old')
    do is = 1, natoms
        read(2,*) (r(l,is), l=1,3)
        read(2,*) (vinf(l,is), l=1,3)
    end do        
    read(2,*) boxlength
    close(2)

    ! 3.2 Opening files to write results
    open(3, file='energy-leap.dat', status='unknown')
    open(4, file='temp-leap.dat', status='unknown')
    open(5, file='rdf_output.dat', status='unknown') ! For g(r) 
    open(6, file='pressure.dat', status='unknown')   ! For pressures

    ! #######################################################################################
    ! 4. Change to reduced units and initialize variables for g(r) and diffusion calculations
    ! #######################################################################################
    call reduced(natoms, r, vinf, boxlength, deltat, epsil, sigma, mass, uvel)

    ! 4.1 RDF implementation
    delg = boxlength / (2.0_dp * real(nhis, dp))
    do i = 1, nhis
        g(i) = 0.0_dp
    end do

    rho = real(natoms, dp) / (boxlength**3) ! Density for g(r) and pressure calculations

    ! 4.2 Initialize counters for diffusion calculations
    nsamp = 1
    dtime = deltat * nsamp
    ntel = 0
    t0 = 0
    do i = 1, tmax
        ntime(i) = 0.0_dp
        vacf(i) = 0.0_dp
        r2t(i) = 0.0_dp
    end do
    
    ! 4.3 Initialize positions and velocities at t=0 for diffusion calculations
    do is = 1, natoms
        do l = 1, 3
            runw(l,is) = r(l,is)
        end do
    end do

    ! ################################################
    ! 5. Start the loop to generate new configurations
    ! ################################################ 
    do i = 1, nconf
        call forces(natoms, r, boxlength, accel, rc, epot, g, delg, virial)
        call velpos(natoms, vinf, accel, deltat, r, runw, nf, ecin, temp, boxlength)

        etot = ecin + epot
        write(3,*) i*deltat, etot 
        write(4,*) i*deltat, temp

        ! 5.1 Compute only the pression in rc (not taking into account the tail correction) includint tail correction
        p = rho * temp + virial / (3.0_dp * (boxlength**3))

        p = p + (16.0_dp * pi * rho**2) * ((2.0_dp/3.0_dp)*(1.0_dp/rc**9) - (1.0_dp/rc**3)) / (3.0_dp*(boxlength**3))
        write(6,*) i*deltat, p

        ! 5.2 Update the counters and store positions and velocities for VACF and MSD calculations (code proposed)
        ntel = ntel + 1
        if (mod(ntel, it0) == 0) then
            t0 = t0 + 1
            tt0 = mod(t0 - 1, t0max) + 1
            time0(tt0) = ntel
            do is = 1, natoms
                do l = 1, 3     ! Positions and velocities in 3D
                    r0(l, is, tt0) = runw(l, is)
                    v0(l, is, tt0) = vinf(l, is)
                end do
            end do
        end if
        do t = 1, min(t0, t0max)
            delt = ntel - time0(t) + 1
            if (delt <= tmax) then
                ntime(delt) = ntime(delt) + 1.0_dp                
                do is = 1, natoms
                    vacf(delt) = vacf(delt) + (vinf(1,is)*v0(1,is,t) + vinf(2,is)*v0(2,is,t) + vinf(3,is)*v0(3,is,t))
                    r2t(delt) = r2t(delt) + (runw(1,is)-r0(1,is,t))**2 + (runw(2,is)-r0(2,is,t))**2 + (runw(3,is)-r0(3,is,t))**2
                end do
            end if
        end do
    end do

    ! ##################
    ! 6. Result analysis
    ! ##################

    ! 6.1 Normalize g(r)
    do i = 1, nhis
        rr = delg * (real(i, dp) + 0.5_dp) 
        vb = ((real(i+1, dp))**3 - (real(i, dp))**3) * delg**3 
        nid = (4.0_dp / 3.0_dp) * pi * vb
        g(i) = g(i) / (real(natoms, dp) * real(nconf, dp) * nid) 
        write(5, *) rr, g(i)
    end do
    
    ! 6.2 Normalize VACF and MSD and write to file
    open(7, file='diffusion.dat', status='unknown')
    do i = 1, tmax      ! Normalization of VACF and MSD (code proposed)
        time_array(i) = dtime * real(i - 1, dp)

        if (ntime(i) > 0.0_dp) then
            real_time = dtime * real(i - 1, dp)
            vacf(i) = vacf(i) / (real(natoms, dp) * ntime(i))
            r2t(i) = r2t(i) / (real(natoms, dp) * ntime(i))
            
            write(7, '(3F15.6)') real_time, vacf(i), r2t(i)
        end if
    end do

    ! 6.3 Fourier Transform of VACF to obtain the spectrum
    call simpson(200, vacf, time_array, spectrum)

    open(8, file='spectrum.dat', status='unknown')
    do w_idx = 0, num_w
        omega_array(w_idx) = real(w_idx, dp)
        write(8, *) omega_array(w_idx), spectrum(w_idx)
    end do
    
    close(3)
    close(4)
    close(5)
    close(6)
    close(7)
    close(8)

    ! 6.4 Saving last configuration in A and A/ps 
    open(11, file='newconf.data', status='unknown')
    do is = 1, natoms
        write(11,*) (r(l,is)*sigma, l=1,3)
        write(11,*) (vinf(l,is)*uvel, l=1,3)
    end do        
    write(11,*) boxlength*sigma
    close(11)

end program leapfroglj

! *********************************************************
! subroutine reduced
! *********************************************************
subroutine reduced(natoms, r, vinf, boxlength, deltat, epsil, sigma, mass, uvel)
    implicit none
    integer, parameter :: dp = kind(0.d0)
    
    integer, intent(in) :: natoms
    real(dp), intent(inout) :: r(3,1000), vinf(3,1000)
    real(dp), intent(inout) :: boxlength, deltat
    real(dp), intent(in) :: epsil, sigma, mass
    real(dp), intent(out) :: uvel
    
    real(dp) :: rgas, utime
    integer :: is, l

    rgas = 8.314472673_dp ! J/(mol*K) 
    utime = sigma * sqrt(mass/epsil) * sqrt(10.0_dp/rgas)
    uvel = sigma / utime  ! unit of velocity, expressed in A/ps

    boxlength = boxlength / sigma
    deltat = deltat / utime
    
    do is = 1, natoms
        do l = 1, 3
            r(l,is) = r(l,is) / sigma
            vinf(l,is) = vinf(l,is) / uvel 
        end do
    end do
end subroutine reduced

! *********************************************************
! subroutine forces
! *********************************************************
subroutine forces(natoms, r, boxlength, accel, rc, epot, g, delg, virial)
    implicit none
    integer, parameter :: dp = kind(0.d0)

    integer, intent(in) :: natoms
    real(dp), intent(in) :: r(3,1000), boxlength, rc
    real(dp), intent(out) :: accel(3,1000), epot
    
    real(dp) :: pot
    integer :: is, js

    ! Modificado: Se declaran explícitamente los intent para estas variables
    real(dp), intent(inout) :: g(500) 
    real(dp), intent(in) :: delg
    real(dp), intent(out) :: virial

    ! Inicialización
    accel = 0.0_dp 
    epot = 0.0_dp 
    virial = 0.0_dp

    ! Atom-atom interactions
    do is = 1, natoms - 1
        do js = is + 1, natoms
            call lj(is, js, r, boxlength, accel, rc, pot, g, delg, virial)
            epot = epot + pot
        end do
    end do

end subroutine forces

! *********************************************************
! subroutine Lennard-Jones
! *********************************************************
subroutine lj(is, js, r, boxlength, accel, rc, pot, g, delg, virial)
    implicit none
    integer, parameter :: dp = kind(0.d0)

    integer, intent(in) :: is, js
    real(dp), intent(in) :: r(3,1000), boxlength, rc, delg
    real(dp), intent(inout) :: accel(3,1000)
    real(dp), intent(out) :: pot
    
    real(dp) :: rij(3), rr2, rr, rijl
    real(dp) :: ynvrr2, ynvrr6, ynvrr12, forcedist
    integer :: l

    integer :: ig
    real(dp), intent(inout) :: g(500)
    real(dp), intent(inout) :: virial

    rr2 = 0.0_dp
    pot = 0.0_dp
    
    do l = 1, 3
        rijl = r(l,js) - r(l,is)
        rij(l) = rijl - boxlength * anint(rijl/boxlength)
        rr2 = rr2 + rij(l)*rij(l)
    end do
    rr = sqrt(rr2)
    
    ! Cálculo para RDF (Añadida validación para que ig no supere 500)
    if (rr < boxlength / 2.0_dp) then 
        ig = int(rr / delg) + 1
        if (ig <= 500) then
            g(ig) = g(ig) + 2.0_dp
        end if
    end if 

    if (rr < rc) then
        ynvrr2 = 1.0_dp / rr2
        ynvrr6 = ynvrr2 * ynvrr2 * ynvrr2
        ynvrr12 = ynvrr6 * ynvrr6
        forcedist = 24.0_dp * (2.0_dp*ynvrr12 - ynvrr6) * ynvrr2
        pot = 4.0_dp * (ynvrr12 - ynvrr6)  
        virial = virial + forcedist*rr2
        
        do l = 1, 3
            accel(l,is) = accel(l,is) - forcedist*rij(l)
            accel(l,js) = accel(l,js) + forcedist*rij(l)
        end do
    end if
end subroutine lj

! *********************************************************
! subroutine velpos
! *********************************************************
! Calculating velocity at instants t+delta/2 and t 
! and getting temporal position at instant t + deltat.
subroutine velpos(natoms, vinf, accel, deltat, r, runw, nf, ecin, temp, boxlength)
    implicit none
    integer, parameter :: dp = kind(0.d0)

    integer, intent(in) :: natoms, nf
    real(dp), intent(inout) :: vinf(3,1000), r(3,1000), runw(3,1000)
    real(dp), intent(in) :: accel(3,1000), deltat, boxlength
    real(dp), intent(out) :: ecin, temp
    
    real(dp) :: v2, vsup, v
    integer :: is, l

    ecin = 0.0_dp
    
    do is = 1, natoms
        v2 = 0.0_dp
        do l = 1, 3
            vsup = vinf(l,is) + accel(l,is)*deltat
            r(l,is) = r(l,is) + vsup*deltat
            runw(l,is) = runw(l,is) + vsup*deltat
            v = (vsup + vinf(l,is)) / 2.0_dp
            v2 = v2 + v*v
            vinf(l,is) = vsup
        end do
        ecin = ecin + 0.5_dp * v2
    end do
    
    temp = 2.0_dp * ecin / real(nf, dp)

    ! Applying periodic boundary conditions
    do is = 1, natoms
        do l = 1, 3
            if (r(l,is) < 0.0_dp) r(l,is) = r(l,is) + boxlength
            if (r(l,is) > boxlength) r(l,is) = r(l,is) - boxlength
        end do
    end do
end subroutine velpos

! ***********************************************************
! subroutine simpson (Fourier Transform of VACF)
! Modified from practice 1 to include the cutoff of the array
! ***********************************************************
subroutine simpson(n, C, t, result)
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = KIND(0.d0)
    INTEGER, INTENT(IN) :: n        ! Added n for cutoff of the array
    REAL(dp), INTENT(IN) :: C(n), t(n)
    REAL(dp), INTENT(OUT) :: result(0:500)
    
    INTEGER :: i, w
    REAL(dp) :: h, s
    REAL(dp), ALLOCATABLE :: f(:)

    h = t(2) - t(1)
    ALLOCATE(f(n))

    DO w = 0, 500
        DO i = 1, n
            f(i) = C(i) * COS(REAL(w, dp) * t(i)) 
        END DO

        s = 0.0_dp
        DO i = 2, n - 1, 2
            s = s + h / 3.0_dp * (f(i-1) + 4.0_dp * f(i) + f(i+1))
        END DO
        result(w) = s
    END DO

    DEALLOCATE(f)
end subroutine simpson