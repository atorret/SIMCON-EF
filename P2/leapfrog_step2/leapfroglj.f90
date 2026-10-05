!*****************************************************************************
!          Molecular Dynamics code to simulate at NVE ensemble
!          a system of N atoms of Ar inside a cubic box, interacting
!          through a truncated Lennard-Jones force field at 2.5 sigma.
!          Leap-frog Verlet integration algorithm.
!
!          Shared inputs (do not copy into other versions):
!            ../data/leap-lj.data    simulation parameters
!            ../data/leap-conf.data  initial configuration
!          A new version is a sibling directory of this one and opens
!          those same files. Energy, temperature, pressure
!          (kinetic, virial and tail), the radial distribution
!          function g(r), and the final configuration are written
!          in this directory.
!
!*****************************************************************************
program leapfroglj
  implicit double precision(a-h,o-z)
  double precision mass, nid

  ! 1. Defining dimensions
  dimension r(3,1000), vinf(3,1000), accel(3,1000)
  parameter (nhis = 1000)
  dimension g(nhis)

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

  ! 5. Start the loop to generate new configurations
  !    delg is the RDF bin width. Bins cover r in [0, L/2].
  !    Tail corrections assume g(r) = 1 for r > rc (reduced units):
  !      E_tail = (8*pi/3)*rho*N*[1/(3*rc^9) - 1/rc^3]
  !      P_tail = (16*pi/3)*rho^2*[2/(3*rc^9) - 1/rc^3]
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

  write(8,'(a)') '# t ekin epot etail etot pkin pvir ptail ptot'

  print *, 'Iniciando simulacion MD (', nconf, ' pasos)...'
  do i = 1, nconf
     if (mod(i, 1000) == 0) then
        print *, '  -> Completado paso ', i, ' de ', nconf
     end if
     call forces(natoms, r, boxlength, accel, rc, epot, vir, nhis, g, delg)
     
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
     write(8,*) i*deltat, ecin, epot, etail, etot_corr, pkin, pvir, ptail, ptot
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
  !    Bin j is the shell [(j-1)*delg, j*delg). Its center is
  !    (j-0.5)*delg and its volume is (4/3)*pi*(r_out^3-r_in^3).
  !    nid (ideal-gas occupancy of that shell) is real: with the
  !    implicit typing, a name starting with n would be integer.
  open(15, file='g-leap.dat', status='unknown')
  do j = 1, nhis
     rr = delg*(dfloat(j) - 0.5d0)
     vb = (dfloat(j)**3 - dfloat(j-1)**3)*delg**3
     nid = (4.d0/3.d0)*pi*vb*rho
     g(j) = g(j)/(dfloat(nconf)*dfloat(natoms)*nid)
     write(15,*) rr, g(j)
  end do
  close(15)
  close(5)

  ! 7. Saving last configuration in A and A/ps
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
!*********************************************************
!              subroutine reduced
!*********************************************************
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
!*********************************************************
!              subroutine forces
!*********************************************************
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
!*********************************************************
!              subroutine Lennard-Jones
!*********************************************************
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
  ! ig = 1 is the shell [0, delg). Fortran arrays are 1-based.
  if (rr .lt. boxlength/2.d0) then
     ig = int(rr/delg) + 1
     if (ig .ge. 1 .and. ig .le. nhis) g(ig) = g(ig) + 2.d0
  end if

  return
end subroutine lj

!*********************************************************
!*********************************************************
!              subroutine velpos
!*********************************************************
!*********************************************************

! Calculating velocity at instants t+delta/2 and t
! and position at instant t + deltat.

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
        if (r(l,is) .lt. 0) r(l,is) = r(l,is) + boxlength
        if (r(l,is) .gt. boxlength) r(l,is) = r(l,is) - boxlength
     end do
  end do

  return
end subroutine velpos
