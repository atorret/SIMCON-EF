!*****************************************************************************
!          Monte Carlo code to simulate at NVT ensemble
!          a system of N atoms of Ar inside a cubic box, interacting
!          through a truncated Lennard-Jones potential at 2.5 sigma.
!          Metropolis algorithm.
!*****************************************************************************
program mcnvt
  implicit double precision(a-h,o-z)
  double precision nid

  ! 1. Defining dimensions
  dimension r(3,1000)
  parameter (nhis = 500)
  dimension g(nhis)

  INCLUDE '../chdir_to_code.inc'

  ! 2. Reading data and computing related quantities
  open(1, file='mc-lj.data', status='old')
  read(1,*) ncycl, nsamp
  read(1,*) natoms
  read(1,*) tref
  read(1,*) sigma, epsil
  read(1,*) deltax
  close(1)

  rc = 2.5d0  ! cutoff radius in units of sigma

  ! 3. Reading initial configuration (positions) in A
  open(2, file='mc-conf.data', status='old')
  do is = 1, natoms
     read(2,*) (r(l,is), l=1,3)
  end do
  read(2,*) boxlength
  close(2)

  ! Opening files to write results
  open(3, file='energy-mc.dat', status='unknown')
  open(4, file='pressure-mc.dat', status='unknown')

  ! 4. Change to reduced units
  call reduced(natoms, r, boxlength, deltax, sigma)

  temp = tref/epsil
  beta = 1.d0/temp
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
  macct = 0
  mnonacct = 0
  nsamples = 0
  sumu = 0.d0
  sumu2 = 0.d0
  sump = 0.d0
  sump2 = 0.d0

  ! 5. Start the loop to generate new configurations
  do icycl = 1, ncycl
     call mcmove(natoms, r, beta, boxlength, rc, deltax, macct, mnonacct)
     if (mod(icycl,nsamp) .eq. 0) then
        call sample(natoms, r, boxlength, rc, rho, temp, etail, ptail, &
                    upot, press, nhis, g, delg)
        upa = upot/dfloat(natoms)
        write(3,*) icycl, upa
        write(4,*) icycl, press
        nsamples = nsamples + 1
        sumu = sumu + upa
        sumu2 = sumu2 + upa*upa
        sump = sump + press
        sump2 = sump2 + press*press
     end if
  end do
  close(3)
  close(4)

  umean = sumu/nsamples
  pmean = sump/nsamples
  usd = dsqrt(sumu2/nsamples - umean**2)
  psd = dsqrt(sump2/nsamples - pmean**2)
  acc = dfloat(macct)/dfloat(macct + mnonacct)

  print *, 'T*          = ', temp
  print *, 'rho*        = ', rho
  print *, 'rc*         = ', rc
  print *, 'delta*      = ', deltax
  print *, 'ncycl       = ', ncycl, '  samples = ', nsamples
  print *, 'acceptance  = ', acc
  print *, '<U*>/N      = ', umean, '  sigma_U = ', usd
  print *, '<P*>        = ', pmean, '  sigma_P = ', psd
  print *, 'U*tail/N    = ', etail/dfloat(natoms), '  P*tail = ', ptail

  ! 6. Radial distribution function in reduced units
  open(12, file='gdr-mc.dat', status='unknown')
  do j = 1, nhis
     rr = delg*(dfloat(j) - 0.5d0)
     vb = (dfloat(j)**3 - dfloat(j-1)**3)*delg**3
     nid = (4.d0/3.d0)*pi*vb*rho
     g(j) = g(j)/(dfloat(nsamples)*dfloat(natoms)*nid)
     write(12,*) rr, g(j)
  end do
  close(12)

  ! 7. Saving last configuration in A
  open(11, file='newconf.data', status='unknown')
  do is = 1, natoms
     write(11,*) (r(l,is)*sigma, l=1,3)
  end do
  write(11,*) boxlength*sigma
  close(11)

  stop
end program mcnvt

!*********************************************************
!              subroutine reduced
!*********************************************************

subroutine reduced(natoms, r, boxlength, deltax, sigma)
  implicit double precision(a-h,o-z)
  dimension r(3,1000)

  boxlength = boxlength/sigma
  deltax = deltax/sigma
  do is = 1, natoms
     do l = 1, 3
        r(l,is) = r(l,is)/sigma
     end do
  end do

  return
end subroutine reduced

!*********************************************************
!              subroutine sample
!*********************************************************

! U* and P* (virial route) of the configuration, with tail corrections.
subroutine sample(natoms, r, boxlength, rc, rho, temp, etail, ptail, &
                  upot, press, nhis, g, delg)
  implicit double precision(a-h,o-z)
  dimension r(3,1000)
  dimension g(nhis)

  upot = 0.d0
  vir = 0.d0

  ! atom-atom interactions
  do is = 1, natoms-1
     do js = is+1, natoms
        call lj(is, js, r, boxlength, rc, pot, virij, nhis, g, delg)
        upot = upot + pot
        vir = vir + virij
     end do
  end do

  upot = upot + etail
  press = rho*temp + vir/(3.d0*boxlength**3) + ptail

  return
end subroutine sample

!*********************************************************
!              subroutine Lennard-Jones
!*********************************************************

subroutine lj(is, js, r, boxlength, rc, pot, virij, nhis, g, delg)
  implicit double precision(a-h,o-z)
  dimension r(3,1000)
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
!              subroutine ener
!*********************************************************

subroutine ener(natoms, io, rio, enerio, r, boxlength, rc)
  implicit double precision(a-h,o-z)
  dimension r(3,1000)
  dimension rio(3)

  enerio = 0.d0
  do js = 1, natoms
     if (js .ne. io) then
        call ljmc(io, js, rio, r, boxlength, rc, pot)
        enerio = enerio + pot
     end if
  end do

  return
end subroutine ener

!*********************************************************
!              subroutine Lennard-Jones (MC)
!*********************************************************

subroutine ljmc(io, js, rio, r, boxlength, rc, pot)
  implicit double precision(a-h,o-z)
  dimension r(3,1000)
  dimension rij(3), rio(3)

  rr2 = 0.d0
  pot = 0.d0
  do l = 1, 3
     rijl = r(l,js) - rio(l)
     rij(l) = rijl - boxlength*dnint(rijl/boxlength)
     rr2 = rr2 + rij(l)*rij(l)
  end do

  rr = dsqrt(rr2)

  if (rr .lt. rc) then
     ynvrr2 = 1.d0/rr2
     ynvrr6 = ynvrr2*ynvrr2*ynvrr2
     ynvrr12 = ynvrr6*ynvrr6
     pot = 4.d0*(ynvrr12 - ynvrr6)
  end if

  return
end subroutine ljmc

!*********************************************************
!              subroutine mcmove
!*********************************************************

subroutine mcmove(natoms, r, beta, boxlength, rc, deltax, macct, mnonacct)
  implicit double precision(a-h,o-z)
  dimension r(3,1000), rold(3), rnew(3)

  ! select a particle at random
  roll = rand()
  io = int(roll*natoms) + 1
  do l = 1, 3
     rold(l) = r(l,io)
  end do

  call ener(natoms, io, rold, enold, r, boxlength, rc)
  do l = 1, 3
     rnew(l) = rold(l) + (rand() - 0.5)*deltax
  end do

  call ener(natoms, io, rnew, ennew, r, boxlength, rc)
  roll = rand()
  if (roll .lt. exp(-beta*(ennew - enold))) then
     do l = 1, 3
        r(l,io) = rnew(l)
     end do
     macct = macct + 1
  else
     mnonacct = mnonacct + 1
  end if

  ! Applying periodic boundary conditions
  do l = 1, 3
     if (r(l,io) .lt. 0) r(l,io) = r(l,io) + boxlength
     if (r(l,io) .gt. boxlength) r(l,io) = r(l,io) - boxlength
  end do

  return
end subroutine mcmove
