! ***********************************************************************
!
!   Copyright (C) 2010-2024  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

module fingering_modes

   use const_def, only: dp
   use num_lib
   use math_lib
   use utils_lib
   use auto_diff

   implicit none

   private

   public :: eval_fastest_fingering

contains

   subroutine eval_fastest_fingering(Pr, tau, R_0, lam_hat, l2_hat, ierr)

      real(dp), intent(in)  :: Pr
      real(dp), intent(in)  :: tau
      type(auto_diff_real_1var_order1), intent(in)  :: R_0
      type(auto_diff_real_1var_order1), intent(out) :: lam_hat
      type(auto_diff_real_1var_order1), intent(out) :: l2_hat
      integer, intent(out)  :: ierr

      type(auto_diff_real_1var_order1) :: a_0
      type(auto_diff_real_1var_order1) :: a_1
      type(auto_diff_real_1var_order1) :: a_2
      type(auto_diff_real_1var_order1) :: a_3
      type(auto_diff_real_1var_order1) :: q
      type(auto_diff_real_1var_order1) :: r
      type(auto_diff_real_1var_order1) :: snq
      type(auto_diff_real_1var_order1) :: tlam

      ! Find the growth rate lam_hat and wavenumber-squared l2_hat of
      ! the maximally-growing fingering mode, by directly solving the
      ! cubic problem (eqn. 19 of Brown, Garaud, & Stellmach, ApJ
      ! 768:34, 2013). The cubic is recast in terms of the reduced
      ! growth rate tlam = lam_hat/l2_hat, and the single positive
      ! root is used to set up lam_hat


      ! Initialize auto_diff variables
      a_0%val = 0d0
      a_0%d1val1 = 0d0

      a_1%val = 0d0
      a_1%d1val1 = 0d0

      a_2%val = 0d0
      a_2%d1val1 = 0d0

      a_3%val = 0d0
      a_3%d1val1 = 0d0

      q%val = 0d0
      q%d1val1 = 0d0

      r%val = 0d0
      r%d1val1 = 0d0

      snq%val = 0d0
      snq%d1val1 = 0d0

      tlam%val = 0d0
      tlam%d1val1 = 0d0

      ! Set up cubic coefficients

      a_3 = Pr - R_0 - Pr*R_0 + tau
      a_2 = -2*(R_0-1)*(Pr + tau + pr*tau)
      a_1 = pr + tau - 4*Pr*(R_0-1)*tau - (1+Pr)*R_0*pow2(tau)
      a_0 = -2*Pr*tau*(R_0*tau - 1)

      ! Determine q and r

      q = a_1/(3*a_3) - pow2(a_2/(3*a_3))
      r = a_1*a_2/(6*pow2(a_3)) - a_0/(2*a_3) - pow3(a_2/(3*a_3))

      ! Sanity check (ensures that the cubic has three real roots)

      if (pow2(r%val) + pow3(q%val) > 0._dp) then
         write(*, *) 'invalid cubic in eval_fastest_fingering'
         ierr = -1
         return
      end if

      ! Calculate the root

      if (q%val < 0._dp) then
         snq = sqrt(-q)
         tlam = 2*snq*cos(acos(r/pow3(snq))/3) - a_2/(3*a_3)
      else
         tlam = -a_2/(3*a_3)
      endif

      ! Evaluate outputs

      l2_hat = sqrt(l4_hat_(Pr, tau, R_0, tlam))
      lam_hat = l2_hat*tlam

      ! Finish

      ierr = 0

      return

   end subroutine eval_fastest_fingering

   !****

   function l4_hat_(Pr, tau, R_0, tlam)

      real(dp), intent(in) :: Pr
      real(dp), intent(in) :: tau
      type(auto_diff_real_1var_order1), intent(in) :: R_0
      type(auto_diff_real_1var_order1), intent(in) :: tlam
      type(auto_diff_real_1var_order1)             :: l4_hat_

      ! Evaluate the horizontal wavenumber^4 given the reduced growth
      ! rate tlam

      l4_hat_ = Pr*(1d0 + tlam - R_0*(tau + tlam))/(R_0*(1d0 + tlam)*(Pr + tlam)*(tau + tlam))

      ! Finish

      return

   end function l4_hat_

end module fingering_modes
