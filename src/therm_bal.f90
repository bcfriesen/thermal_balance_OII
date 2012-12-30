program main
implicit none

!-----------------------VARIABLE DECLARATIONS-------------------------------
integer, parameter :: dp = kind(1.0d0)
! M: all the coefficients that will multiply the N's
! A: all the spontaneous transition rates
! N: before calling DGESV (the LAPACK linear solver): all zeros
!    (except one)
!    after calling DGESV: all the fractional populations of levels
!    1 through 5
! q: all the collision rates
real(kind=dp) :: M(5,5), A(5,5), N(5), q(5,5)

! i: loop index
! pivot: something DGESV needs to solve this matrix equation (not
!        relevant to the user)
! ok: DGESV error status
integer :: i, pivot(5), ok

! N_e: electron number density
! T_e: electron temperature
real(kind=dp) :: N_e, T_e

! A_0: a constant for calculating the q's
real(kind=dp), parameter :: A_0 = 8.629d-06

! omega: degeneracies of the levels
integer :: omega(5)

! k: Boltzmann's constant (erg/K)
real(kind=dp), parameter :: k_B = 1.381d-16

! E: energies between the levels (erg)
real(kind=dp) :: E(5,5)

! upsilon: velocity-averaged collision strengths
real(kind=dp) :: upsilon(5,5)

! S: spin of each level
! L: orbital angular momentum of each level
! J: total angular momentum of each level
real(kind=dp) :: S(5), L(5), J(5)
!---------------------------------------------------------------------------


!----------------------------MISC. PARAMETERS-------------------------------
! assume electron number density of 10^2 cm^{-3}
N_e = 1.0d+02

! Enter spin and angular momentum stuff
S(1) = 1.5d+00
S(2) = 0.5d+00
S(3) = 0.5d+00
S(4) = 0.5d+00
S(5) = 0.5d+00
L(1) = 0.0d+00
L(2) = 2.0d+00
L(3) = 2.0d+00
L(4) = 1.0d+00
L(5) = 1.0d+00
J(1) = 1.5d+00
J(2) = 2.5d+00
J(3) = 2.5d+00
J(4) = 1.5d+00
J(5) = 1.5d+00

! degeneracy = 2J+1
do i = 1, 5
  omega(i) = (2.0d+00 * J(i)) + 1.0d+00
end do

!----------------------------LEVEL ENERGIES---------------------------------
! Looked these up on NIST. Units: erg

E(1,1) = 0.0d+00
E(1,2) = 5.319d-12
E(1,3) = 5.323d-12
E(1,4) = 8.028d-12
E(1,5) = 8.028d-12
      
E(2,1) = E(1,2)
E(2,2) = 0.0d+00
E(2,3) = 3.979d-15
E(2,4) = 2.709d-12
E(2,5) = 2.710d-12

E(3,1) = E(1,3)
E(3,2) = E(2,3)
E(3,3) = 0.0d+00
E(3,4) = 2.705d-12
E(3,5) = 2.706d-12

E(4,1) = E(1,4)
E(4,2) = E(2,4)
E(4,3) = E(3,4)
E(4,4) = 0.0d+00
E(4,5) = 3.952d-16

E(5,1) = E(1,5)
E(5,2) = E(2,5)
E(5,3) = E(3,5)
E(5,4) = E(4,5)
E(5,5) = 0.0d+00
!---------------------------------------------------------------------------


!--------------------------COLLISION STRENGTHS------------------------------
! Found these in (AGN)^2 Table 3.7.  Note that we have to use
! Eq. 3.21 in (AGN)^2 for 4S <--> 2D and 4S <--> 2P

upsilon(1,1) = 0.00d+0
upsilon(1,2) = (((2.0d+0 * J(2)) + 1.0d+0) / ((2.0d+0 * S(2)) + &
               1.0d+0) * ((2.0d+0 * L(2)) + 1.0d+0)) * 1.34d0
upsilon(1,3) = (((2.0d+0 * J(3)) + 1.0d+0) / ((2.0d+0 * S(3)) + &
               1.0d+0) * ((2.0d+0 * L(3)) + 1.0d+0)) * 1.34d0
upsilon(1,4) = (((2.0d+0 * J(4)) + 1.0d+0) / ((2.0d+0 * S(4)) + &
               1.0d+0) * ((2.0d+0 * L(4)) + 1.0d+0)) * 0.40d0
upsilon(1,5) = (((2.0d+0 * J(5)) + 1.0d+0) / ((2.0d+0 * S(5)) + &
               1.0d+0) * ((2.0d+0 * L(5)) + 1.0d+0)) * 0.40d0

upsilon(2,1) = upsilon(1,2)
upsilon(2,2) = 0.0d+0
upsilon(2,3) = 1.17d+0
! Can't find this value anywhere so I'm using the same as 2D_{5/2} <--> 2P_{1/2}
upsilon(2,4) = 0.33d0
upsilon(2,5) = 0.33d0

upsilon(3,1) = upsilon(1,3)
upsilon(3,2) = upsilon(2,3)
upsilon(3,3) = 0.0d+0
upsilon(3,4) = 0.82d0
upsilon(3,5) = 0.28d0

upsilon(4,1) = upsilon(1,4)
upsilon(4,2) = upsilon(2,4)
upsilon(4,3) = upsilon(3,4)
upsilon(4,4) = 0.0d+0
upsilon(4,5) = 0.157d0
   
upsilon(5,1) = upsilon(1,5)
upsilon(5,2) = upsilon(2,5)
upsilon(5,3) = upsilon(3,5)
upsilon(5,4) = upsilon(4,5)
upsilon(5,5) = 0.0d+0
!---------------------------------------------------------------------------


!----------------------------RATE DATA--------------------------------------
! Now we enter the A values in by hand.  These are taken from NIST.

! 2D_{5/2} --> 4S_{3/2}
A(2,1) = 2.86d-5
! 2D_{3/2} --> 4S_{3/2}
A(3,1) = 1.59d-4
! 2D_{3/2} --> 2D_{5/2}
A(3,2) = 1.30d-7
! 2P_{3/2} --> 4S_{3/2}
A(4,1) = 5.22d-2
! 2P_{3/2} --> 2D_{5/2}
A(4,2) = 9.07d-2
! 2P_{3/2} --> 2D_{3/2}
A(4,3) = 1.49d-2
! 2P_{1/2} --> 4S_{3/2}
A(5,1) = 2.12d-2
! 2P_{1/2} --> 2D_{5/2}
A(5,2) = 5.19d-2
! 2P_{1/2} --> 2D_{3/2}
A(5,3) = 7.74d-2
! 2P_{1/2} --> 2P_{3/2}
A(5,4) = 1.41d-10
!---------------------------------------------------------------------------

! Start at a low temperature
T_e = 3.0d+2

! Write all this stuff to therm_bal.out
open (unit=11, file='therm_bal_OII.out')

! Start the temperature loop!
do

! Now we calculute the q values using Eqs. 3.18 and 3.20 in (AGN)^2
  q(2,1) = (8.629d-6 / dsqrt(T_e)) * (upsilon(1,2) / omega(2))
  q(1,2) = (omega(2) / omega(1)) * q(2,1) * dexp(E(2,1) / (k_B * T_e))

  q(3,1) = (8.629d-6 / dsqrt(T_e)) * (upsilon(1,3) / omega(3))
  q(1,3) = (omega(3) / omega(1)) * q(3,1) * dexp(E(3,1) / (k_B * T_e))

  q(3,2) = (8.629d-6 / dsqrt(T_e)) * (upsilon(2,3) / omega(3))
  q(2,3) = (omega(3) / omega(2)) * q(3,2) * dexp(E(3,2) / (k_B * T_e))

  q(4,1) = (8.629d-6 / dsqrt(T_e)) * (upsilon(1,4) / omega(4))
  q(1,4) = (omega(4) / omega(1)) * q(4,1) * dexp(E(4,1) / (k_B * T_e))

  q(4,2) = (8.629d-6 / dsqrt(T_e)) * (upsilon(2,4) / omega(4))
  q(2,4) = (omega(4) / omega(2)) * q(4,2) * dexp(E(4,2) / (k_B * T_e))

  q(4,3) = (8.629d-6 / dsqrt(T_e)) * (upsilon(3,4) / omega(4))
  q(3,4) = (omega(4) / omega(3)) * q(4,3) * dexp(E(4,3) / (k_B * T_e))

  q(5,1) = (8.629d-6 / dsqrt(T_e)) * (upsilon(1,5) / omega(5))
  q(1,5) = (omega(5) / omega(1)) * q(5,1) * dexp(E(5,1) / (k_B * T_e))

  q(5,2) = (8.629d-6 / dsqrt(T_e)) * (upsilon(2,5) / omega(5))
  q(2,5) = (omega(5) / omega(2)) * q(5,2) * dexp(E(5,2) / (k_B * T_e))

  q(5,3) = (8.629d-6 / dsqrt(T_e)) * (upsilon(3,5) / omega(5))
  q(3,5) = (omega(5) / omega(3)) * q(5,3) * dexp(E(5,3) / (k_B * T_e))

  q(5,4) = (8.629d-6 / dsqrt(T_e)) * (upsilon(4,5) / omega(5))
  q(4,5) = (omega(5) / omega(4)) * q(5,4) * dexp(E(5,4) / (k_B * T_e))
!---------------------------------------------------------------------------

!---------------------------SET UP MATRICES---------------------------------
! level 1 equilibrium
  M(1,1) = N_e * (q(1,2) + q(1,3) + q(1,4) + q(1,5))
  M(1,2) = -((N_e * q(2,1)) + A(2,1))
  M(1,3) = -((N_e * q(3,1)) + A(3,1))
  M(1,4) = -((N_e * q(4,1)) + A(4,1))
  M(1,5) = -((N_e * q(5,1)) + A(5,1))

! level 2 equilibrium
  M(2,1) = -(N_e * q(1,2))
  M(2,2) = (N_e * (q(2,1))) + A(2,1)
  M(2,3) = -A(3,2)
  M(2,4) = -A(4,2)
  M(2,5) = -A(5,2)

! level 3 equilibrium
  M(3,1) = -(N_e * q(1,3))
  M(3,2) = 0.0d0
  M(3,3) = (N_e * (q(3,1))) + A(3,2) + A(3,1)
  M(3,4) = -A(4,3)
  M(3,5) = -A(5,3)

! level 4 equilbrium
  M(4,1) = -(N_e * q(1,4))
  M(4,2) = 0.0d0
  M(4,3) = 0.0d0
  M(4,4) = (N_e * (q(4,1))) + A(4,3) + A(4,2) + A(4,1)
  M(4,5) = -A(5,4)

! We scrapped level 5 stuff in order to provide this constraint
! equation, which says that the sum of the level populations
! must be the total population of O II at temperature T.
  M(5,1) = 1.0d+0
  M(5,2) = 1.0d+0
  M(5,3) = 1.0d+0
  M(5,4) = 1.0d+0
  M(5,5) = 1.0d+0
      
! Before we call DGESV, this will be the RHS of the matrix
! equation, which is all zeros (due to equilibrium), except for
! row 5, which is N(O II) from the constraint equation on the
! level populations.  After DGESV runs, this matrix will
! contain the fractional populations of the first 5 levels.
  N(1) = 0.0d+0
  N(2) = 0.0d+0
  N(3) = 0.0d+0
  N(4) = 0.0d+0
  N(5) = 1.0d+0
!-----------------------------------------------------------------------


!------------------------SOLVE THESE EQUATIONS!-------------------------
  call DGESV(5, 1, M, 5, pivot, N, 5, ok)
!-----------------------------------------------------------------------


! Write out the results
  write(11,'(6es12.4e2)') T_e, N(1), N(2), N(3), N(4), N(5)

! Do uniform increments in log(T_e)
  T_e = dlog(T_e)
  T_e = T_e + 1.0d-2
  T_e = dexp(T_e)

! Quit at T_e = 10^6 K
  if (T_e .gt. 1.0d+6) exit
end do

close(11)
write(*,*) 'SUCCESS!!!!'

end program main
!-------------------------------------------------------------------------------
