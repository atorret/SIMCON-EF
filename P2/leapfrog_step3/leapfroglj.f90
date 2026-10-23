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
   dimension r(3,1000), vinf(3,1000), accel(3,1000), runw(3,1000)
   parameter (nhis = 1000)
   dimension g(nhis)
   ! VACF and MSD accumulated over multiple time origins (circular buffer)
   integer, allocatable :: time0(:), ntime(:)
   double precision, allocatable :: x0(:,:,:), v0(:,:,:), vacf(:), r2t(:)

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
   sumt = 0.d0
   sumt2 = 0.d0
   sume = 0.d0
   sume2 = 0.d0

   ! Window of ~5 time units (VACF decayed, MSD linear), ncor a multiple of 6.
   ! A new time origin every it0 steps; n0max origins fit in the window.
   it0 = 5
   ncor = 6*(min(nint(5.d0/deltat), nconf)/6)
   n0max = ncor/it0 + 1
   allocate(time0(n0max), x0(3,natoms,n0max), v0(3,natoms,n0max))
   allocate(ntime(0:ncor), vacf(0:ncor), r2t(0:ncor))
   ntime = 0
   vacf = 0.d0
   r2t = 0.d0
   n0 = 0

   write(8,'(a)') '# t etot pkin pvir ptail ptot'

   print *, 'Iniciando simulacion MD (', nconf, ' pasos)...'
   do i = 1, nconf
      if (mod(i, 1000) == 0) then
         print *, '  -> Completado paso ', i, ' de ', nconf
      end if
      call forces(natoms, r, boxlength, accel, rc, epot, vir, nhis, g, delg)
      call velpos(natoms, vinf, accel, deltat, r, nf, ecin, temp, boxlength)

      ! Leap-frog displacement is vinf*deltat; added without PBC.
      runw(:,1:natoms) = runw(:,1:natoms) + vinf(:,1:natoms)*deltat

      ! Store a new time origin, overwriting the oldest one.
      if (mod(i, it0) == 0) then
         n0 = n0 + 1
         k = mod(n0-1, n0max) + 1
         time0(k) = i
         x0(:,:,k) = runw(:,1:natoms)
         v0(:,:,k) = vinf(:,1:natoms)
      end if
      ! Accumulate v(t).v(0) and [r(t)-r(0)]^2 for every stored origin.
      do k = 1, min(n0, n0max)
         it = i - time0(k)
         if (it <= ncor) then
            ntime(it) = ntime(it) + 1
            vacf(it) = vacf(it) + sum(vinf(:,1:natoms)*v0(:,:,k))
            r2t(it) = r2t(it) + sum((runw(:,1:natoms) - x0(:,:,k))**2)
         end if
      end do

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
      sumt = sumt + temp
      sumt2 = sumt2 + temp*temp
      sume = sume + etot
      sume2 = sume2 + etot*etot
   end do
   close(3)
   close(4)
   close(8)
 
   ! Mean and standard deviation of T and E_tot over the run
   tmean = sumt/nconf
   emean = sume/nconf
   tsd = dsqrt(sumt2/nconf - tmean**2)
   esd = dsqrt(sume2/nconf - emean**2)
   print *, 'steps  = ', nconf
   print *, '<T>    = ', tmean, '  sigma_T = ', tsd
   print *, '<Etot> = ', emean, '  sigma_E = ', esd
 
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

   ! 7. VACF and MSD averaged over atoms and time origins.
   vacf = vacf/(dfloat(natoms)*ntime)
   r2t = r2t/(dfloat(natoms)*ntime)
   open(16, file='c-leap.dat', status='unknown')
   open(17, file='r2-leap.dat', status='unknown')
   open(19, file='D-leap.dat', status='unknown')
   do it = 0, ncor
      write(16,*) it*deltat, vacf(it)
      write(17,*) it*deltat, r2t(it)
      if (it > 0) write(19,*) it*deltat, r2t(it)/(6.d0*it*deltat)
   end do
   close(16)
   close(17)
   close(19)
   print *, '  C(0) =', vacf(0)

   ! 8. Spectral density S(w) = integral C(t) cos(w t) dt, w = 0, ..., 100.
   !    Composite Simpson 3/8 rule on the first quarter of the window
   !    (C(t) ~ 0 beyond), with a multiple of 3 intervals.
   nsdw = 3*(ncor/12)
   open(18, file='c_spectrum-leap.dat', status='unknown')
   do iw = 0, 100
      s = 0.d0
      do it = 0, nsdw
         if (it == 0 .or. it == nsdw) then
            wgt = 1.d0
         else if (mod(it, 3) == 0) then
            wgt = 2.d0
         else
            wgt = 3.d0
         end if
         s = s + wgt*vacf(it)*dcos(iw*it*deltat)
      end do
      write(18,*) dfloat(iw), 3.d0*deltat/8.d0*s
   end do
   close(18)

   ! 9. Einstein: least-squares fit r2(t) = 6*D*t + b over the diffusive
   !    regime (the first quarter of the window, ballistic, is discarded).
   !    yint is used instead of "intercept": implicit typing would make it integer.
   ifit = ncor/4
   npts = ncor - ifit + 1
   sx = 0.d0
   sy = 0.d0
   sxx = 0.d0
   sxy = 0.d0
   do it = ifit, ncor
      tt = it*deltat
      sx = sx + tt
      sy = sy + r2t(it)
      sxx = sxx + tt*tt
      sxy = sxy + tt*r2t(it)
   end do
   slope = (dfloat(npts)*sxy - sx*sy)/(dfloat(npts)*sxx - sx*sx)
   yint = (sy - slope*sx)/dfloat(npts)
   d_einstein = slope/6.d0

   open(20, file='diffusion-leap.dat', status='unknown')
   write(20,'(a)') '# D_einstein slope intercept tmin tmax'
   write(20,*) d_einstein, slope, yint, ifit*deltat, ncor*deltat
   close(20)
   print *, '  D (Einstein fit, reduced units) =', d_einstein

   ! 10. Saving last configuration in A and A/ps
   open(11, file='newconf.data', status='unknown')
   do is = 1, natoms
      write(11,*) (r(l,is)*sigma, l=1,3)
      write(11,*) (vinf(l,is)*uvel, l=1,3)
   end do
   write(11,*) boxlength*sigma
   close(11)

   stop
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