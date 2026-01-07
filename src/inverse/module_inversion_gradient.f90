!
! © 2025. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by
! Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by
! Triad National Security, LLC, and the U.S. Department of Energy/National
! Nuclear Security Administration. The Government is granted for itself and
! others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!


module inversion_gradient

    use mod_parameters
    use gradient
    use inversion_regularization

    implicit none

contains

    subroutine compute_gradient

        ! Initialization
        call zero_gradient

        ! Compute model parameter gradients
        yn_misfit_only = .false.
        yn_save_adjsrc = .true.
        yn_reconstruct = .true.
        ! yn_save_entire = .true.
        ! yn_encoded_synthetic_stage = .true.
        call compute_gradient_shots
        ! yn_encoded_synthetic_stage = .false.
        call mpibarrier

        ! Process and regularize gradients
        call process_gradient
        call regularize_gradient
        call output_gradient
        call mpibarrier

        ! Print progress
        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Gradient computation completed. ')
        end if

    end subroutine compute_gradient

end module inversion_gradient
