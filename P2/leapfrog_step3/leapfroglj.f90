!*****************************************************************************
!          Molecular Dynamics code to simulate at NVE ensemble
!          a system of N atoms of Ar inside a cubic box, interacting
!          through a truncated Lennard-Jones force field at 2.5 sigma.
!          Leap-frog Verlet integration algorithm.
!*****************************************************************************
program leapfroglj
   implicit double precision(a-h,o-z)
   double precision mass, nid
 
   ! 1. Defining dimensions
   dimension r(3,1000), vinf(3,1000), accel(3,1000), runw(3,1000), rold(3,1000)
   dimension c(0:10000), r2(0:10000)
   parameter (nhis = 1000)
   dimension g(nhis)
   dimension t(0:10000), D(10000), c_spectrum(0:1000)
   ! Multiple time origins for VACF and MSD (lecture SAMPLE algorithm).
   ! init=0 initialize, init=1 accumulate a sample, init=2 normalize.
   integer :: ncor, it0, ntmax, t0max, ntel, t0
   integer, allocatable :: time0(:), ntime(:)
   double precision, allocatable :: vx0(:,:,:), x0(:,:,:), vacf(:), r2t(:)
 
   INCLUDE '../../chdir_to_code.inc'
 
   ! 2. Reading data and computing related quantities
   open(1, file='../data/leap-lj.data', status='old')
   read(1,*) nconf
   read(1,*) natoms
   read(1,*) mass
   read(1,*) sigma, epsil
   read(1,*) deltat
   close(1)
 
   nf = 3*natoms - 3  ! number of degrees of freedom
   rc = 2.5d0         ! cutoff radius in units of sigma
 
   ! 3. Reading initial configuration (positions, velocities) in A and A/ps
   open(2, file='../data/leap-conf.data', status='old')
   do is = 1, natoms
      read(2,*) (r(l,is), l=1,3)
      read(2,*) (vinf(l,is), l=1,3)
   end do
   read(2,*) boxlength
   close(2)
 
   ! Opening files to write results
   open(3, file='energy-leap.dat', status='unknown')
   open(4, file='temp-leap.dat', status='unknown')
   open(8, file='thermo-leap.dat', status='unknown')
 
   ! 4. Change to reduced units
   call reduced(natoms, r, vinf, boxlength, deltat, epsil, sigma, &
                mass, uvel)

   ! runw holds positions without periodic wrapping, needed for the MSD.
   runw = r

   ! 5. Start the loop to generate new configurations
   pi = 4.d0*datan(1.d0)
   delg = boxlength/(2.d0*dfloat(nhis))
   rho = dfloat(natoms)/boxlength**3
   etail = (8.d0*pi/3.d0)*rho*dfloat(natoms)* &
           (1.d0/(3.d0*rc**9) - 1.d0/rc**3)
   ptail = (16.d0*pi/3.d0)*rho**2* &
           (2.d0/(3.d0*rc**9) - 1.d0/rc**3)
   do j = 1, nhis
      g(j) = 0.d0
   end do
   ekin_sum = 0.d0
   epot_sum = 0.d0
   pkin_sum = 0.d0
   pvir_sum = 0.d0

   ! Correlation window of about 5 reduced time units: the VACF has
   ! decayed and the MSD is already linear. ncor is a multiple of 6 so
   ! the Simpson 3/8 rule has a complete set of panels. A new time
   ! origin is taken every it0 samples; t0max keeps every origin still
   ! inside the window. ntmax is tmax in the lecture listing (1-based,
   ! with C(0) stored in vacf(1)).
   it0 = 5
   ncor = nint(5.d0/deltat)
   ncor = ncor - mod(ncor, 6)
   if (ncor < 6) ncor = 6
   if (ncor > nconf) ncor = nconf - mod(nconf, 6)
   ntmax = ncor + 2
   t0max = ncor/it0 + 3
   allocate(time0(t0max), ntime(ntmax))
   allocate(vx0(3, natoms, t0max), x0(3, natoms, t0max))
   allocate(vacf(ntmax), r2t(ntmax))
   call sample(0)
 
   write(8,'(a)') '# t etot pkin pvir ptail ptot'
 
   print *, 'Iniciando simulacion MD (', nconf, ' pasos)...'
   do i = 1, nconf
      if (mod(i, 1000) == 0) then
         print *, '  -> Completado paso ', i, ' de ', nconf
      end if
      call forces(natoms, r, boxlength, accel, rc, epot, vir, nhis, g, delg)
      rold = r
      call velpos(natoms, vinf, accel, deltat, r, nf, ecin, temp, boxlength)

      ! Unfold this step (SAMPLE does not apply periodic wrapping to r2t).
      do j = 1, natoms
         do l = 1, 3
            dr = r(l,j) - rold(l,j)
            runw(l,j) = runw(l,j) + dr - boxlength*dnint(dr/boxlength)
         end do
      end do
      call sample(1)
      
      ! Conserved NVE energy (no tail). Thermodynamic totals include it.
      etot = ecin + epot
      etot_corr = etot + etail
      pkin = 2.d0*ecin/(3.d0*boxlength**3)
      pvir = vir/(3.d0*boxlength**3)
      ptot = pkin + pvir + ptail
      ekin_sum = ekin_sum + ecin
      epot_sum = epot_sum + epot
      pkin_sum = pkin_sum + pkin
      pvir_sum = pvir_sum + pvir
      write(3,*) i*deltat, etot
      write(4,*) i*deltat, temp
      write(8,*) i*deltat, etot_corr, pkin, pvir, ptail, ptot
   end do
   close(3)
   close(4)
   close(8)
 
   ekin_avg = ekin_sum/dfloat(nconf)
   epot_avg = epot_sum/dfloat(nconf)
   etot_avg = ekin_avg + epot_avg + etail
   pkin_avg = pkin_sum/dfloat(nconf)
   pvir_avg = pvir_sum/dfloat(nconf)
   ptot_avg = pkin_avg + pvir_avg + ptail
 
   open(9, file='averages-leap.dat', status='unknown')
   write(9,*) 'nconf', nconf
   write(9,*) 'ekin', ekin_avg
   write(9,*) 'epot', epot_avg
   write(9,*) 'etail', etail
   write(9,*) 'etot', etot_avg
   write(9,*) 'pkin', pkin_avg
   write(9,*) 'pvir', pvir_avg
   write(9,*) 'ptail', ptail
   write(9,*) 'ptot', ptot_avg
   close(9)
   print *, 'Mean values in reduced units'
   print *, '  Ekin  =', ekin_avg
   print *, '  Epot  =', epot_avg
   print *, '  Etail =', etail
   print *, '  Etot  =', etot_avg
   print *, '  Pkin  =', pkin_avg
   print *, '  Pvir  =', pvir_avg
   print *, '  Ptail =', ptail
   print *, '  Ptot  =', ptot_avg
 
   ! 6. Radial distribution function in reduced units.
   open(15, file='g-leap.dat', status='unknown')
   do j = 1, nhis
      rr = delg*(dfloat(j) - 0.5d0)
      vb = (dfloat(j)**3 - dfloat(j-1)**3)*delg**3
      nid = (4.d0/3.d0)*pi*vb*rho
      g(j) = g(j)/(dfloat(nconf)*dfloat(natoms)*nid)
      write(15,*) rr, g(j)
   end do
   close(15)

   ! 7. Normalize VACF and MSD (SAMPLE init=2) and copy to 0-based arrays.
   call sample(2)
   do i = 0, ncor
      c(i) = vacf(i+1)
      r2(i) = r2t(i+1)
      t(i) = i*deltat
   end do
   print *, 'VACF and MSD averaged over time origins'
   print *, '  origins at t = 0     =', ntime(1)
   print *, '  origins at t = tmax  =', ntime(ncor+1)
   print *, '  C(0) =', c(0)

   open(16, file='c-leap.dat', status='unknown')
   open(17, file='r2-leap.dat', status='unknown')
   do i = 0, ncor
      write(16,*) t(i), c(i)
      write(17,*) t(i), r2(i)
   end do
   close(16)
   close(17)

   ! 8. Spectral density of the averaged VACF, Simpson 3/8 rule.
   !    Frequencies w = 0, 1, ..., 100. Only those are written: the
   !    previous loop stopped at 100 and the rest of the file was unset.
   call spectrum(c(0:ncor), t(0:ncor), c_spectrum)
   open(18, file='c_spectrum-leap.dat', status='unknown')
   do i = 0, 100
      write(18,*) dfloat(i), c_spectrum(i)
   end do
   close(18)

   ! 9. Instantaneous diffusion coefficient from the mean square displacement
   open(19, file='D-leap.dat', status='unknown')
   do i = 1, ncor
      D(i) = r2(i)/(6.d0*t(i))
      write(19,*) t(i), D(i)
   end do
   close(19)

   ! 10. Least-squares linear fit r2(t) = 6*D*t + b over the diffusive regime.
   !     The first quarter of the window (ballistic regime) is discarded.
   !     yint is used instead of "intercept": implicit typing would make it integer.
   ifit = ncor/4
   npts = ncor - ifit + 1
   sx = 0.d0
   sy = 0.d0
   sxx = 0.d0
   sxy = 0.d0
   do i = ifit, ncor
      sx = sx + t(i)
      sy = sy + r2(i)
      sxx = sxx + t(i)*t(i)
      sxy = sxy + t(i)*r2(i)
   end do
   slope = (dfloat(npts)*sxy - sx*sy)/(dfloat(npts)*sxx - sx*sx)
   yint = (sy - slope*sx)/dfloat(npts)
   d_einstein = slope/6.d0
   tmin = t(ifit)
   tmax = t(ncor)

   ! Green-Kubo: D = (1/3) * integral of C(t) from 0 to the correlation window.
   d_gk = deltat*(sum(c(0:ncor)) - 0.5d0*c(0) - 0.5d0*c(ncor))/3.d0

   open(20, file='diffusion-leap.dat', status='unknown')
   write(20,'(a)') '# D_einstein slope intercept tmin tmax D_greenkubo'
   write(20,*) d_einstein, slope, yint, tmin, tmax, d_gk
   close(20)
   print *, 'Diffusion coefficient (reduced units)'
   print *, '  D (Einstein fit) =', d_einstein
   print *, '  D (Green-Kubo)   =', d_gk
 
   ! 7. Saving last configuration in A and A/ps
   open(11, file='newconf.data', status='unknown')
   do is = 1, natoms
      write(11,*) (r(l,is)*sigma, l=1,3)
      write(11,*) (vinf(l,is)*uvel, l=1,3)
   end do
   write(11,*) boxlength*sigma
   close(11)
 
   stop

 contains

 !*********************************************************
 !     SAMPLE: VACF and MSD with multiple time origins
 !     init=0 initialize counters
 !     init=1 accumulate a sample of data
 !     init=2 average the results
 !     Positions in x0 / runw are unwrapped (PBC not assumed).
 !*********************************************************

 subroutine sample(init)
   implicit none
   integer, intent(in) :: init
   integer :: i, j, l, tt0, torg, delt
   double precision :: dx

   if (init == 0) then
      ntel = 0
      t0 = 0
      do i = 1, ntmax
         ntime(i) = 0
         vacf(i) = 0.d0
         r2t(i) = 0.d0
      end do
   else if (init == 1) then
      ntel = ntel + 1
      if (mod(ntel, it0) == 0) then
         t0 = t0 + 1
         tt0 = mod(t0-1, t0max) + 1
         time0(tt0) = ntel
         do j = 1, natoms
            do l = 1, 3
               x0(l,j,tt0) = runw(l,j)
               vx0(l,j,tt0) = vinf(l,j)
            end do
         end do
      end if
      do torg = 1, min(t0, t0max)
         delt = ntel - time0(torg) + 1
         if (delt < ntmax) then
            ntime(delt) = ntime(delt) + 1
            do j = 1, natoms
               do l = 1, 3
                  vacf(delt) = vacf(delt) + vinf(l,j)*vx0(l,j,torg)
                  dx = runw(l,j) - x0(l,j,torg)
                  r2t(delt) = r2t(delt) + dx*dx
               end do
            end do
         end if
      end do
   else if (init == 2) then
      do i = 1, ntmax
         if (ntime(i) > 0) then
            vacf(i) = vacf(i)/(dfloat(natoms)*dfloat(ntime(i)))
            r2t(i) = r2t(i)/(dfloat(natoms)*dfloat(ntime(i)))
         end if
      end do
   end if
 end subroutine sample

 !*********************************************************
 !     spectral density, composite Simpson 3/8 rule
 !     S(w) = integral C(t) cos(w t) dt, w = 0, 1, ..., 100
 !*********************************************************

 subroutine spectrum(C, t, result)
   implicit none
   double precision, dimension(:), intent(in) :: C, t
   double precision, dimension(0:), intent(out) :: result
   integer :: i, w, n, panel
   double precision :: h, s
   double precision, allocatable :: f(:)

   n = int(0.25*size(t))
   h = t(2) - t(1)
   allocate(f(n))
   result = 0.d0

   do w = 0, 100
      do i = 1, n
         f(i) = C(i)*cos(dble(w)*t(i))
      end do

      ! Composite Simpson 3/8: panels come in groups of three.
      s = f(1) + f(n)
      do i = 2, n - 1
         if (mod(i - 1, 3) == 0) then
            s = s + 2.d0*f(i)
         else
            s = s + 3.d0*f(i)
         end if
      end do
      result(w) = 3.d0*h/8.d0*s
   end do

   deallocate(f)
 end subroutine spectrum

 end program leapfroglj
 
 !*********************************************************
 !              subroutine reduced
 !*********************************************************
 
 subroutine reduced(natoms, r, vinf, boxlength, deltat, &
                    epsil, sigma, mass, uvel)
   implicit double precision(a-h,o-z)
   double precision mass
   dimension r(3,1000), vinf(3,1000)
 
   rgas = 8.314472673d0  ! J/(mol*K)
   utime = sigma*dsqrt(mass/epsil)*dsqrt(10.d0/rgas)
   uvel = sigma/utime    ! unit of velocity, expressed in A/ps
 
   boxlength = boxlength/sigma
   deltat = deltat/utime
   do is = 1, natoms
      do l = 1, 3
         r(l,is) = r(l,is)/sigma
         vinf(l,is) = vinf(l,is)/uvel
      end do
   end do
 
   return
 end subroutine reduced
 
 !*********************************************************
 !              subroutine forces
 !*********************************************************
 
 subroutine forces(natoms, r, boxlength, accel, rc, epot, vir, nhis, g, delg)
   implicit double precision(a-h,o-z)
   dimension r(3,1000), accel(3,1000)
   dimension g(nhis)
 
   do is = 1, natoms
      do l = 1, 3
         accel(l,is) = 0.d0  ! sets accelerations to 0
      end do
   end do
   epot = 0.d0
   vir = 0.d0
 
   ! atom-atom interactions
   do is = 1, natoms-1
      do js = is+1, natoms
         call lj(is, js, r, boxlength, accel, rc, pot, virij, nhis, g, delg)
         epot = epot + pot
         vir = vir + virij
      end do
   end do
 
   return
 end subroutine forces
 
 !*********************************************************
 !              subroutine Lennard-Jones
 !*********************************************************
 
 subroutine lj(is, js, r, boxlength, accel, rc, pot, virij, nhis, g, delg)
   implicit double precision(a-h,o-z)
   dimension r(3,1000), accel(3,1000)
   dimension rij(3)
   dimension g(nhis)
 
   rr2 = 0.d0
   pot = 0.d0
   virij = 0.d0
 
   do l = 1, 3
      rijl = r(l,js) - r(l,is)
      rij(l) = rijl - boxlength*dnint(rijl/boxlength)
      rr2 = rr2 + rij(l)*rij(l)
   end do
 
   rr = dsqrt(rr2)
 
   if (rr .lt. rc) then
      ynvrr2 = 1.d0/rr2
      ynvrr6 = ynvrr2*ynvrr2*ynvrr2
      ynvrr12 = ynvrr6*ynvrr6
      forcedist = 24.d0*(2.d0*ynvrr12 - ynvrr6)*ynvrr2
      pot = 4.d0*(ynvrr12 - ynvrr6)
      do l = 1, 3
         accel(l,is) = accel(l,is) - forcedist*rij(l)
         accel(l,js) = accel(l,js) + forcedist*rij(l)
      end do
      virij = forcedist*rr2
   end if
 
   ! Histogram of pair distances up to L/2 (minimum image).
   if (rr .lt. boxlength/2.d0) then
      ig = int(rr/delg) + 1
      if (ig .ge. 1 .and. ig .le. nhis) g(ig) = g(ig) + 2.d0
   end if
 
   return
 end subroutine lj
 
 !*********************************************************
 !              subroutine velpos
 !*********************************************************
 
 subroutine velpos(natoms, vinf, accel, deltat, r, nf, ecin, temp, &
                   boxlength)
   implicit double precision(a-h,o-z)
   dimension vinf(3,1000), accel(3,1000), r(3,1000)
 
   ecin = 0.d0
   do is = 1, natoms
      v2 = 0.d0
      do l = 1, 3
         vsup = vinf(l,is) + accel(l,is)*deltat
         r(l,is) = r(l,is) + vsup*deltat
         v = (vsup + vinf(l,is))/2.d0
         v2 = v2 + v*v
         vinf(l,is) = vsup
      end do
      ecin = ecin + 0.5d0*v2
   end do
   temp = 2.d0*ecin/dfloat(nf)
 
   ! Applying periodic boundary conditions
   do is = 1, natoms
      do l = 1, 3
         if (r(l,is) .lt. 0.d0) r(l,is) = r(l,is) + boxlength
         if (r(l,is) .gt. boxlength) r(l,is) = r(l,is) - boxlength
      end do
   end do
 
   return
 end subroutine velpos