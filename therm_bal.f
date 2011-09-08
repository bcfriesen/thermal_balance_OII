c     therm_bal: Calculates fractional populations of the first five
c     levels of O II, taking into account only collisional excitations/
c     deexcitations, and spontaneous emissions

      program main
      implicit none
      
c-----------------------VARIABLE DECLARATIONS-------------------------------
c     M: all the coefficients that will multiply the N's
c     A: all the spontaneous transition rates
c     N: before calling DGESV (the LAPACK linear solver): all zeros
c        (except one)
c        after calling DGESV: all the fractional populations of levels
c        1 through 5
c     q: all the collision rates
      real*8 :: M(5,5), A(5,5), N(5), q(5,5)

c     i: loop index
c     pivot: something DGESV needs to solve this matrix equation (not
c            relevant to the user)
c     ok: DGESV error status
      integer :: i, pivot(5), ok

c     N_e: electron number density
c     T_e: electron temperature
      real*8 :: N_e, T_e

c     A_0: a constant for calculating the q's
      real*8, parameter :: A_0 = 8.629d-06

c     omega: degeneracies of the levels
      integer :: omega(5)

c     k: Boltzmann's constant (cgs)
      real*8, parameter :: k = 1.381d-16

c     E: energies between the levels
      real*8 :: E(5,5)

c     upsilon: velocity-averaged collision strengths
      real*8 :: upsilon(5,5)

c     S: spin of each level
c     L: orbital angular momentum of each level
c     J: total angular momentum of each level
      real*8 :: S(5), L(5), J(5)

c     N_OII: the temperature-dependent number density of O II
      real*8 :: N_OII
c---------------------------------------------------------------------------


c----------------------------MISC. PARAMETERS-------------------------------
c     assume electron number density of 10^2 cm^{-3}      
      N_e = 1.0d+02

c     Enter spin and angular momentum stuff
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

c     degeneracy = 2J+1
      do i = 1, 5
         omega(i) = (2.0d+00 * J(i)) + 1.0d+00
      end do

c----------------------------LEVEL ENERGIES---------------------------------
c     Looked these up on NIST

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
c---------------------------------------------------------------------------


c--------------------------COLLISION STRENGTHS------------------------------
c     Found these in (AGN)^2 Table 3.7.  Note that we have to use
c     Eq. 3.21 in (AGN)^2 for 4S <--> 2D and 4S <--> 2P

      upsilon(1,1) = 0.00d+0
      upsilon(1,2) = (((2.0d+0 * J(2)) + 1.0d+0) / ((2.0d+0 * S(2)) + 
     1               1.0d+0) * ((2.0d+0 * L(2)) + 1.0d+0)) * 1.34d0
      upsilon(1,3) = (((2.0d+0 * J(3)) + 1.0d+0) / ((2.0d+0 * S(3)) + 
     1               1.0d+0) * ((2.0d+0 * L(3)) + 1.0d+0)) * 1.34d0
      upsilon(1,4) = (((2.0d+0 * J(4)) + 1.0d+0) / ((2.0d+0 * S(4)) + 
     1               1.0d+0) * ((2.0d+0 * L(4)) + 1.0d+0)) * 0.40d0
      upsilon(1,5) = (((2.0d+0 * J(5)) + 1.0d+0) / ((2.0d+0 * S(5)) +
     1               1.0d+0) * ((2.0d+0 * L(5)) + 1.0d+0)) * 0.40d0

      upsilon(2,1) = upsilon(1,2)
      upsilon(2,2) = 0.0d+0
      upsilon(2,3) = 1.17d+0
c     Can't find this value anywhere so I'm using the same as 2D_{5/2} <-->
c     2P_{1/2}
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
c---------------------------------------------------------------------------


c----------------------------RATE DATA--------------------------------------
c     Now we enter the A values in by hand.  These are taken from NIST.

c     2D_{5/2} --> 4S_{3/2}
      A(2,1) = 2.86d-5
c     2D_{3/2} --> 4S_{3/2}
      A(3,1) = 1.59d-4
c     2D_{3/2} --> 2D_{5/2}
      A(3,2) = 1.30d-7
c     2P_{3/2} --> 4S_{3/2}
      A(4,1) = 5.22d-2
c     2P_{3/2} --> 2D_{5/2}
      A(4,2) = 9.07d-2
c     2P_{3/2} --> 2D_{3/2}
      A(4,3) = 1.49d-2
c     2P_{1/2} --> 4S_{3/2}
      A(5,1) = 2.12d-2
c     2P_{1/2} --> 2D_{5/2}
      A(5,2) = 5.19d-2
c     2P_{1/2} --> 2D_{3/2}
      A(5,3) = 7.74d-2
c     2P_{1/2} --> 2P_{3/2}
      A(5,4) = 1.41d-10
c---------------------------------------------------------------------------

c     Start at a low temperature
      T_e = 3.0d3

c     Write all this stuff to therm_bal.out
      open (unit=11, file='therm_bal_OII.out')

c     Start the loop!
      do

c        Now we calculute the q values using Eqs. 3.18 and 3.20 in (AGN)^2

         q(2,1) = 0.0d0
         q(1,2) = 0.0d0

         q(3,1) = 0.0d0
         q(1,3) = 0.0d0

         q(3,2) = 0.0d0
         q(2,3) = 0.0d0

         q(4,1) = 0.0d0
         q(1,4) = 0.0d0

         q(4,2) = 0.0d0
         q(2,4) = 0.0d0

         q(4,3) = 0.0d0
         q(3,4) = 0.0d0

         q(5,1) = 0.0d0
         q(1,5) = 0.0d0

         q(5,2) = 0.0d0
         q(2,5) = 0.0d0

         q(5,3) = 0.0d0
         q(3,5) = 0.0d0

         q(5,4) = 0.0d0
         q(4,5) = 0.0d0
c---------------------------------------------------------------------------

c---------------------------SET UP MATRICES---------------------------------
c        level 1 equilibrium
         M(1,1) = N_e * (q(1,2) + q(1,3) + q(1,4) + q(1,5))
         M(1,2) = -((N_e * q(2,1)) + A(2,1))
         M(1,3) = -((N_e * q(3,1)) + A(3,1))
         M(1,4) = -((N_e * q(4,1)) + A(4,1))
         M(1,5) = -((N_e * q(5,1)) + A(5,1))

c        level 2 equilibrium
         M(2,1) = -(N_e * q(1,2))
         M(2,2) = (N_e * (q(2,1)))
     1            + A(2,1)
         M(2,3) = -A(3,2)
         M(2,4) = -A(4,2)
         M(2,5) = -A(5,2)

c        level 3 equilibrium
         M(3,1) = -(N_e * q(1,3))
         M(3,2) = 0.0d0
         M(3,3) = (N_e * (q(3,1)))
     1            + A(3,2) + A(3,1)
         M(3,4) = -A(4,3)
         M(3,5) = -A(5,3)

c        level 4 equilbrium
         M(4,1) = -(N_e * q(1,4))
         M(4,2) = 0.0d0
         M(4,3) = 0.0d0
         M(4,4) = (N_e * (q(4,1)))
     1            + A(4,3) + A(4,2) + A(4,1)
         M(4,5) = -A(5,4)

c        We scrapped level 5 stuff in order to provide this constraint
c        equation, which says that the sum of the level populations
c        must be the total population of O II at temperature T.
         M(5,1) = 1.0d+0
         M(5,2) = 1.0d+0
         M(5,3) = 1.0d+0
         M(5,4) = 1.0d+0
         M(5,5) = 1.0d+0
      
c        Before we call DGESV, this will be the RHS of the matrix
c        equation, which is all zeros (due to equilibrium), except for
c        row 5, which is N(O II) from the constraint equation on the
c        level populations.  After DGESV runs, this matrix will
c        contain the fractional populations of the first 5 levels.
         N(1) = 0.0d+0
         N(2) = 0.0d+0
         N(3) = 0.0d+0
         N(4) = 0.0d+0
         N(5) = 1.0d+0
c-----------------------------------------------------------------------


c------------------------SOLVE THESE EQUATIONS!-------------------------
         call DGESV(5, 1, M, 5, pivot, N, 5, ok)
c-----------------------------------------------------------------------


c        Write out the results
         write(11,'(6es12.4e2)') T_e, N(1), N(2), N(3), N(4), N(5)

c        Do uniform increments in log(T_e)
         T_e = log(T_e)
         T_e = T_e + 1.0d-2
         T_e = exp(T_e)

c        Quit at T_e = 10^6 K
         if (T_e .gt. 2.0d+4) exit
      end do

      close(11)
      write(*,*) 'SUCCESS!!!!'

      end program main
!-------------------------------------------------------------------------------




!-------------------------------------------------------------------------------
      real*8 function N_OII(T)
      implicit none
      real*8 :: T
      N_OII = -(5.3d4 * (T - 8.6d4) * (T - 2.9d4) * (T + 5.0d4)) /
     1        ((T - 5.7d5) * (T - 8.6d4) * (T - 2.9d4) * (T + 5.1d4))
      end function N_OII
!-------------------------------------------------------------------------------
