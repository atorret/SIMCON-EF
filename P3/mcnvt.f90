!*****************************************************************************
!          Monte Carlo code to simulate at NVT ensemble
!          a system of N atoms of Ar inside a cubic box, interacting
!          through a truncated Lennard-Jones potential at 2.5 sigma.
!          Metropolis algorithm with single-particle trial moves.
!          Built from leapfroglj.f90 (Practice 2).
!
!          Inputs (in this directory):
!            mc-lj.data    simulation parameters (last line: rho*)
!            mc-conf.data  only read if rho* <= 0 in mc-lj.data
!                          (positions in A, box in A). If rho* > 0 the
!                          initial FCC lattice (natoms = 4 n^3) is built
!                          inside the program for that density.
!          The first 20% of the cycles are equilibration (not averaged).
!          Outputs:
!            energy-mc.dat    cycle, U*/N   (every nsamp moves)
!            pressure-mc.dat  cycle, P*     (every nsamp moves)
!            gdr-mc.dat       r*, g(r)      (averaged over the samples)
!            newconf.data     last configuration, same format as mc-conf.data
!
!*****************************************************************************
program mcnvt
  implicit double precision(a-h,o-z)
  double precision, parameter :: pi = 3.14159265358979d0

  integer, parameter :: nhis = 500  ! number of g(r) bins

  ! 1. Defining dimensions
  dimension r(3,1000), g(nhis), bas(3,4)
  data bas /0.d0,0.d0,0.d0, 0.5d0,0.5d0,0.d0, 0.5d0,0.d0,0.5d0, 0.d0,0.5d0,0.5d0/

  INCLUDE '../chdir_to_code.inc'

  ! 2. Reading data and computing related quantities
  open(1, file='mc-lj.data', status='old')
  read(1,*) ncycl, nsamp
  read(1,*) natoms
  read(1,*) tref
  read(1,*) sigma, epsil
  read(1,*) deltax
  read(1,*) rhoin
  close(1)

  rc = 2.5d0  ! cutoff radius in units of sigma

  neq = ncycl/5  ! equilibration cycles, not averaged

  ! 3. Initial configuration (positions) in A
  if (rhoin .gt. 0.d0) then
     ! FCC lattice with 4*nc**3 atoms at density rhoin (reduced)
     nc = nint((dfloat(natoms)/4.d0)**(1.d0/3.d0))
     if (4*nc**3 .ne. natoms) stop 'natoms must be 4*n**3 for an FCC start'
     boxlength = (dfloat(natoms)/rhoin)**(1.d0/3.d0)*sigma
     a = boxlength/nc
     is = 0
     do i = 0, nc-1
        do j = 0, nc-1
           do k = 0, nc-1
              do ib = 1, 4
                 is = is + 1
                 r(1,is) = (i + 0.25d0 + bas(1,ib))*a
                 r(2,is) = (j + 0.25d0 + bas(2,ib))*a
                 r(3,is) = (k + 0.25d0 + bas(3,ib))*a
              end do
           end do
        end do
     end do
  else
     open(2, file='mc-conf.data', status='old')
     do is = 1, natoms
        read(2,*) (r(l,is), l=1,3)
     end do
     read(2,*) boxlength
     close(2)
  end if

  ! Opening files to write results
  open(3, file='energy-mc.dat', status='unknown')
  open(4, file='pressure-mc.dat', status='unknown')

  ! 4. Change to reduced units
  call reduced(natoms, r, boxlength, deltax, sigma)

  temp = tref/epsil   ! T*
  beta = 1.d0/temp    ! beta* = 1/T*
  vol = boxlength**3
  rho = dfloat(natoms)/vol

  ! Tail corrections (truncated LJ potential, g(r)=1 beyond rc)
  rci3 = 1.d0/rc**3
  rci9 = rci3**3
  utail = 8.d0/3.d0*pi*dfloat(natoms)*rho*(rci9/3.d0 - rci3)
  ptail = 16.d0/3.d0*pi*rho*rho*(2.d0/3.d0*rci9 - rci3)

  macct = 0
  mnonacct = 0
  nsamples = 0
  sumu = 0.d0
  sumu2 = 0.d0
  sump = 0.d0
  sump2 = 0.d0

  ! g(r) histogram from 0 to boxlength/2
  delg = boxlength/(2.d0*nhis)
  do ig = 1, nhis
     g(ig) = 0.d0
  end do

  ! 5. Start the loop to generate new configurations
  do icycl = 1, ncycl
     call mcmove(natoms, r, beta, boxlength, rc, deltax, macct, mnonacct)
     if (icycl .gt. neq .and. mod(icycl,nsamp) .eq. 0) then
        call sample(natoms, r, boxlength, rc, rho, temp, utail, ptail, &
                    upot, press, nhis, g, delg)
        upa = upot/dfloat(natoms)  ! interaction energy per atom
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

  ! Mean and standard deviation of U*/N and P* over the samples
  umean = sumu/nsamples
  pmean = sump/nsamples
  usd = dsqrt(max(0.d0, sumu2/nsamples - umean**2))
  psd = dsqrt(max(0.d0, sump2/nsamples - pmean**2))
  acc = dfloat(macct)/dfloat(macct + mnonacct)

  print *, 'T*          = ', temp
  print *, 'rho*        = ', rho
  print *, 'rc*         = ', rc
  print *, 'delta*      = ', deltax
  print *, 'ncycl       = ', ncycl, '  samples = ', nsamples
  print *, 'acceptance  = ', acc
  print *, '<U*>/N      = ', umean, '  sigma_U = ', usd
  print *, '<P*>        = ', pmean, '  sigma_P = ', psd
  print *, 'U*tail/N    = ', utail/dfloat(natoms), '  P*tail = ', ptail

  ! Normalising g(r): counts per atom and sample / ideal-gas counts in the shell
  open(12, file='gdr-mc.dat', status='unknown')
  do ig = 1, nhis
     rbin = delg*(dfloat(ig) - 0.5d0)
     vshell = 4.d0/3.d0*pi*(dfloat(ig)**3 - dfloat(ig-1)**3)*delg**3
     g(ig) = g(ig)/(dfloat(nsamples)*dfloat(natoms)*rho*vshell)
     write(12,*) rbin, g(ig)
  end do
  close(12)

  ! 6. Saving last configuration in A
  open(11, file='newconf.data', status='unknown')
  do is = 1, natoms
     write(11,*) (r(l,is)*sigma, l=1,3)
  end do
  write(11,*) boxlength*sigma
  close(11)

  stop
end program mcnvt

!*********************************************************
!*********************************************************
!              subroutine reduced
!*********************************************************
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
!*********************************************************
!              subroutine sample
!*********************************************************
!*********************************************************

! Total interaction energy (U*) and pressure (P*, virial route) of the
! current configuration, both including tail corrections.
! Also accumulates the g(r) histogram.

subroutine sample(natoms, r, boxlength, rc, rho, temp, utail, ptail, &
                  upot, press, nhis, g, delg)
  implicit double precision(a-h,o-z)
  dimension r(3,1000)
  dimension g(nhis)

  upot = 0.d0
  vir = 0.d0

  ! atom-atom interactions
  do is = 1, natoms-1
     do js = is+1, natoms
        call lj(is, js, r, boxlength, rc, pot, rf, nhis, g, delg)
        upot = upot + pot
        vir = vir + rf
     end do
  end do

  upot = upot + utail
  press = rho*temp + vir/(3.d0*boxlength**3) + ptail

  return
end subroutine sample

!*********************************************************
!*********************************************************
!              subroutine Lennard-Jones For Monte Carlo
!*********************************************************
!*********************************************************

! Interaction between atom js and an atom io placed at rio.
! Returns the pair energy (pot) and the virial term r_ij . F_ij (rf).

subroutine ljmc(js, rio, r, boxlength, rc, pot, rf)
  implicit double precision(a-h,o-z)
  dimension r(3,1000)
  dimension rij(3), rio(3)

  rr2 = 0.d0
  pot = 0.d0
  rf = 0.d0
  do l = 1, 3
     rijl = r(l,js) - rio(l)
     rij(l) = rijl - boxlength*dnint(rijl/boxlength)
     rr2 = rr2 + rij(l)*rij(l)
  end do

  if (rr2 .lt. rc*rc) then  ! same as rr < rc, without the square root
     ynvrr2 = 1.d0/rr2
     ynvrr6 = ynvrr2*ynvrr2*ynvrr2
     ynvrr12 = ynvrr6*ynvrr6
     pot = 4.d0*(ynvrr12 - ynvrr6)
     rf = 24.d0*(2.d0*ynvrr12 - ynvrr6)
  end if

  return
end subroutine ljmc

!*********************************************************
!*********************************************************
!              subroutine Lennard-Jones (sampling)
!*********************************************************
!*********************************************************

! Interaction between atoms is and js of the configuration.
! Returns the pair energy (pot) and the virial term r_ij . F_ij (virij),
! and adds the pair to the g(r) histogram (bins of width delg).

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

   if (rr.lt.boxlength/2.d0) then
      ig = int(rr/delg) + 1   ! bins 1..nhis
      g(ig) = g(ig) + 2.d0    ! pair counted for both atoms

      if (rr.lt.rc) then
         ynvrr1 = 1.d0/rr
         ynvrr2 = ynvrr1*ynvrr1
         ynvrr6 = ynvrr2*ynvrr2*ynvrr2
         ynvrr12 = ynvrr6*ynvrr6
         forcedist = 24.d0*(2.d0*ynvrr12 - ynvrr6)*ynvrr1  ! |F_ij|
         pot = 4.d0*(ynvrr12 - ynvrr6)
         virij = forcedist*rr                             ! r_ij . F_ij
      end if
   end if

   return
end subroutine lj

!*********************************************************
!*********************************************************
!              subroutine ener
!*********************************************************
!*********************************************************

! Interaction energy of atom io, placed at rio, with all the other atoms.

subroutine ener(natoms, io, rio, enerio, r, boxlength, rc)
  implicit double precision(a-h,o-z)
  dimension r(3,1000)
  dimension rio(3)

  enerio = 0.d0
  do js = 1, natoms
     if (js .ne. io) then
        call ljmc(js, rio, r, boxlength, rc, pot, rf)
        enerio = enerio + pot
     end if
  end do

  return
end subroutine ener

!*********************************************************
!*********************************************************
!              subroutine mcmove
!*********************************************************
!*********************************************************

! One Metropolis trial move: displace a random atom by up to deltax/2
! in each direction and accept with probability min(1, exp(-beta*dU)).

subroutine mcmove(natoms, r, beta, boxlength, rc, deltax, macct, mnonacct)
  implicit double precision(a-h,o-z)
  dimension r(3,1000), rold(3), rnew(3)

  ! select a particle at random
  call random_number(roll)
  io = min(int(roll*natoms) + 1, natoms)
  do l = 1, 3
     rold(l) = r(l,io)
  end do

  ! energy of the old configuration
  call ener(natoms, io, rold, enold, r, boxlength, rc)

  ! trial position and energy of the new configuration
  do l = 1, 3
     call random_number(roll)
     rnew(l) = rold(l) + (roll - 0.5d0)*deltax
  end do
  call ener(natoms, io, rnew, ennew, r, boxlength, rc)

  ! Metropolis acceptance criterion
  call random_number(roll)
  if (roll .lt. dexp(-beta*(ennew - enold))) then
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
