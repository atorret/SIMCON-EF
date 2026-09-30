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
!          those same files. Energy, temperature, the radial
!          distribution function g(r), and the final configuration
!          are written in this directory.
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

  ! 4. Change to reduced units
  call reduced(natoms, r, vinf, boxlength, deltat, epsil, sigma, &
               mass, uvel)

  ! 5. Start the loop to generate new configurations
  !    delg is the RDF bin width. Bins cover r in [0, L/2].
  pi = 4.d0*datan(1.d0)
  delg = boxlength/(2.d0*dfloat(nhis))
  do j = 1, nhis
     g(j) = 0.d0
  end do

  do i = 1, nconf
     call forces(natoms, r, boxlength, accel, rc, epot, nhis, g, delg)
     call velpos(natoms, vinf, accel, deltat, r, nf, ecin, temp, &
                 boxlength)
     etot = ecin + epot
     write(3,*) i*deltat, etot
     write(4,*) i*deltat, temp
  end do
  close(3)
  close(4)

  ! 6. Radial distribution function in reduced units.
  !    Bin j is the shell [(j-1)*delg, j*delg). Its center is
  !    (j-0.5)*delg and its volume is (4/3)*pi*(r_out^3-r_in^3).
  !    nid (ideal-gas occupancy of that shell) is real: with the
  !    implicit typing, a name starting with n would be integer.
  rho = dfloat(natoms)/boxlength**3
  open(5, file='g-leap.dat', status='unknown')
  do j = 1, nhis
     rr = delg*(dfloat(j) - 0.5d0)
     vb = (dfloat(j)**3 - dfloat(j-1)**3)*delg**3
     nid = (4.d0/3.d0)*pi*vb*rho
     g(j) = g(j)/(dfloat(nconf)*dfloat(natoms)*nid)
     write(5,*) rr, g(j)
  end do
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

subroutine forces(natoms, r, boxlength, accel, rc, epot, nhis, g, delg)
  implicit double precision(a-h,o-z)
  dimension r(3,1000), accel(3,1000)
  dimension g(nhis)

  do is = 1, natoms
     do l = 1, 3
        accel(l,is) = 0.d0  ! sets accelerations to 0
     end do
  end do
  epot = 0.d0

  ! atom-atom interactions
  do is = 1, natoms-1
     do js = is+1, natoms
        call lj(is, js, r, boxlength, accel, rc, pot, nhis, g, delg)
        epot = epot + pot
     end do
  end do

  return
end subroutine forces

!*********************************************************
!*********************************************************
!              subroutine Lennard-Jones
!*********************************************************
!*********************************************************

subroutine lj(is, js, r, boxlength, accel, rc, pot, nhis, g, delg)
  implicit double precision(a-h,o-z)
  dimension r(3,1000), accel(3,1000)
  dimension rij(3)
  dimension g(nhis)

  rr2 = 0.d0
  pot = 0.d0

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
