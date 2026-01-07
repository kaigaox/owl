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


module inversion_regularization

    use mod_parameters
    use regularization

#ifdef _dim2_
#define model_dimension dimension(:, :)
#endif

#ifdef _dim3_
#define model_dimension dimension(:, :, :)
#endif

    implicit none

contains

    !
    !> Add L2-norm regularization term to gradient for a single parameter
    !
    subroutine add_l2reg_single_parameter(reg, model, grad, name)  !

        real, model_dimension, intent(in) :: model !, admm_model
        real, model_dimension, intent(inout) :: reg, grad
        character(len=*), intent(in) :: name
        character(len=1024) :: file_mask

        real :: reg_coef, reg_scale
        logical :: const_reg
        real, allocatable, model_dimension :: grad_mask

        ! Update dual variable
        reg = model - reg

        ! Strength of regularization
        call readpar_xlogical(file_parameter, 'const_reg', const_reg, .false., iter*1.0)
        if (const_reg) then

            call readpar_xfloat(file_parameter, 'reg_lambda_'//tidy(name), reg_coef, 0.0, iter*1.0)

        else

            call readpar_xfloat(file_parameter, 'reg_scale_'//tidy(name), reg_scale, 0.2, iter*1.0)
            if (rankid == 0) then
                call warn(date_time_compact()//' Regularization scale for '//tidy(name)//' = '//num2str(reg_scale, '(es)'))
            end if

            ! Compute regularization coefficient for each parameter
            if (maxval(abs(reg)) == 0) then
                reg_coef = 0.0
            else
                reg_coef = mean(grad, 2)/mean(reg, 2)*reg_scale
            end if

        end if

        if (rankid == 0) then
            call warn(date_time_compact()//' Regularization lambda for '//tidy(name)//' = '//num2str(reg_coef, '(es)'))
        end if

        ! modify gradients
        grad = grad + reg_coef*reg

        ! Mask gradient again
        if (any(gradient_processing == 'mask')) then
            call readpar_xstring(file_parameter, 'grad_mask', file_mask, file_mask, iter*1.0)
            call prepare_model_single_parameter(grad_mask, 'mask', file_mask, update=.false.)
            grad = grad*grad_mask
        end if

    end subroutine

    !
    !> L2 regularization
    !
    subroutine add_l2reg

        integer :: i

        do i = 1, nmodel
            if (yn_regularize_model .or. yn_regularize_source) then
                call add_l2reg_single_parameter(model_reg(i)%array, model_m(i)%array, &
                    model_grad(i)%array, model_name(i))
            end if
        end do

        call mpibarrier

        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Regularization finished. ')
        end if

    end subroutine add_l2reg

    !
    !> Regularize gradient
    !
    subroutine regularize_gradient

        if (yn_regularize_model .or. yn_regularize_source) then
            call add_l2reg
        end if

    end subroutine regularize_gradient

    !
    !> Add regularization
    !
    subroutine compute_regularization

        ! Regularize model parameters
        if (yn_regularize_model) then
            call model_regularization
        end if

        !        ! Regularize source parameters
        !        if (yn_regularize_source) then
        !            call source_regularization
        !        end if

    end subroutine

end module
