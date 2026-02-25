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


module mod_model

    use libflit
    use mod_parameters

    implicit none

    integer :: nx0, ny0, nz0
    real :: dx0, dy0, dz0
    real :: ox0, oy0, oz0

    integer :: nx, ny, nz
    real :: dx, dy, dz
    real :: ox, oy, oz

    real :: dz_max

    real :: xmin = -float_large
    real :: ymin = -float_large
    real :: zmin = -float_large
    real :: xmax = +float_large
    real :: ymax = +float_large
    real :: zmax = +float_large

    real :: shot_xbeg, shot_ybeg, shot_zbeg
    real :: shot_xend, shot_yend, shot_zend
    integer :: shot_nxbeg, shot_nybeg, shot_nzbeg
    integer :: shot_nxend, shot_nyend, shot_nzend
    integer :: shot_nx, shot_ny, shot_nz

    integer :: pml

#ifdef _dim2_
    type(meta_array2_real), allocatable, dimension(:) :: model_m
    type(meta_array2_real), allocatable, dimension(:) :: model_m_backup
    type(meta_array2_real), allocatable, dimension(:) :: model_aux
    type(meta_array2_real), allocatable, dimension(:) :: model_grad
    type(meta_array2_real), allocatable, dimension(:) :: model_srch
    type(meta_array2_real), allocatable, dimension(:) :: model_reg
#endif

#ifdef _dim3_
    type(meta_array3_real), allocatable, dimension(:) :: model_m
    type(meta_array3_real), allocatable, dimension(:) :: model_m_backup
    type(meta_array3_real), allocatable, dimension(:) :: model_aux
    type(meta_array3_real), allocatable, dimension(:) :: model_grad
    type(meta_array3_real), allocatable, dimension(:) :: model_srch
    type(meta_array3_real), allocatable, dimension(:) :: model_reg
#endif

contains

    subroutine set_regular_space

        ! Original model dimensions
        call readpar_int(file_parameter, 'nx', nx0, 0, required=.true.)
        call readpar_int(file_parameter, 'ny', ny0, 1)
        call readpar_int(file_parameter, 'nz', nz0, 0, required=.true.)
        call readpar_float(file_parameter, 'dx', dx0, 1.0, required=.true.)
        call readpar_float(file_parameter, 'dy', dy0, 1.0)
        call readpar_float(file_parameter, 'dz', dz0, 1.0, required=.true.)
        call readpar_float(file_parameter, 'ox', ox0, 0.0)
        call readpar_float(file_parameter, 'oy', oy0, 0.0)
        call readpar_float(file_parameter, 'oz', oz0, 0.0)

        ! Target model dimensions
        call readpar_int(file_parameter, 'nnx', nx, nx0)
        call readpar_int(file_parameter, 'nny', ny, ny0)
        call readpar_int(file_parameter, 'nnz', nz, nz0)
        call readpar_float(file_parameter, 'ddx', dx, dx0)
        call readpar_float(file_parameter, 'ddy', dy, dy0)
        call readpar_float(file_parameter, 'ddz', dz, dz0)
        call readpar_float(file_parameter, 'oox', ox, ox0)
        call readpar_float(file_parameter, 'ooy', oy, oy0)
        call readpar_float(file_parameter, 'ooz', oz, oz0)

        ! Target model ranges
        call readpar_float(file_parameter, 'xmin', xmin, ox)
        call readpar_float(file_parameter, 'xmax', xmax, ox + (nx - 1)*dx)
        call assert(xmax >= xmin, ' <set_regular_space> Error: xmax must >= xmin')
        call readpar_float(file_parameter, 'ymin', ymin, oy)
        call readpar_float(file_parameter, 'ymax', ymax, oy + (ny - 1)*dy)
        call assert(ymax >= ymin, ' <set_regular_space> Error: ymax must >= ymin')
        call readpar_float(file_parameter, 'zmin', zmin, oz)
        call readpar_float(file_parameter, 'zmax', zmax, oz + (nz - 1)*dz)
        call assert(zmax >= zmin, ' <set_regular_space> Error: zmax must >= zmin')

        ! Target model origins
        ox = xmin
        oy = ymin
        oz = zmin
        nx = ceiling((xmax - xmin)/dx) + 1
        ny = ceiling((ymax - ymin)/dy) + 1
        nz = ceiling((zmax - zmin)/dz) + 1

        ! Check if requires model interpolation
        if (nx /= nx0 .or. ny /= ny0 .or. nz /= nz0 .or. &
                dx /= dx0 .or. dy /= dy0 .or. dz /= dz0 .or. &
                ox /= ox0 .or. oy /= oy0 .or. oz /= oz0) then
            require_model_interp = .true.
        end if

        ! Check non-negativeness
        call assert(nx >= 1 .and. ny >= 1 .and. nz >= 1, ' <set_regular_space> Error: n# must >= 1')
        call assert(dx > 0 .and. dy > 0 .and. dz > 0, ' <set_regular_space> Error: d# must > 0')

        ! Print original and target model dimensions
        if (rankid == 0) then

#ifdef _dim2_
            call warn(' ')
            call warn(' Original model dimensions: ')
            call warn('      nx, nz = '//num2str(nx0)//', '//num2str(nz0))
            call warn('      dx, dz = '//num2str(dx0, '(es)')//', '//num2str(dz0, '(es)'))
            call warn('      ox, oz = '//num2str(ox0, '(es)')//', '//num2str(oz0, '(es)'))
            call warn('      ex, ez = ' &
                //num2str(ox0 + (nx0 - 1)*dx0, '(es)')//', ' &
                //num2str(oz0 + (nz0 - 1)*dz0, '(es)'))
            call warn(' ')
            call warn(' Target model dimensions: ')
            call warn('      nx, nz = '//num2str(nx)//', '//num2str(nz))
            call warn('      dx, dz = '//num2str(dx, '(es)')//', '//num2str(dz, '(es)'))
            call warn('      ox, oz = '//num2str(ox, '(es)')//', '//num2str(oz, '(es)'))
            call warn('      ex, ez = ' &
                //num2str(ox + (nx - 1)*dx, '(es)')//', ' &
                //num2str(oz + (nz - 1)*dz, '(es)'))
            call warn(' ')
#endif

#ifdef _dim3_
            call warn(' ')
            call warn(' Original model dimensions: ')
            call warn('      nx, ny, nz = '//num2str(nx0)//', '//num2str(ny0)//', '//num2str(nz0))
            call warn('      dx, dy, dz = '//num2str(dx0, '(es)')//', '//num2str(dy0, '(es)')//', '//num2str(dz0, '(es)'))
            call warn('      ox, oy, oz = '//num2str(ox0, '(es)')//', '//num2str(oy0, '(es)')//', '//num2str(oz0, '(es)'))
            call warn('      ex, ey, ez = ' &
                //num2str(ox0 + (nx0 - 1)*dx0, '(es)')//', ' &
                //num2str(oy0 + (ny0 - 1)*dy0, '(es)')//', ' &
                //num2str(oz0 + (nz0 - 1)*dz0, '(es)'))
            call warn(' ')
            call warn(' Target model dimensions: ')
            call warn('      nx, ny, nz = '//num2str(nx)//', '//num2str(ny)//', '//num2str(nz))
            call warn('      dx, dy, dz = '//num2str(dx, '(es)')//', '//num2str(dy, '(es)')//', '//num2str(dz, '(es)'))
            call warn('      ox, oy, oz = '//num2str(ox, '(es)')//', '//num2str(oy, '(es)')//', '//num2str(oz, '(es)'))
            call warn('      ex, ey, ez = ' &
                //num2str(ox + (nx - 1)*dx, '(es)')//', ' &
                //num2str(oy + (ny - 1)*dy, '(es)')//', ' &
                //num2str(oz + (nz - 1)*dz, '(es)'))
            call warn(' ')
#endif

        end if

        ! Number of PML layers
        call readpar_int(file_parameter, 'npml', pml, 15)

        ! For free-surface modeling, dz_max constrains the maximum grid size along z
        if (yn_free_surface) then
            call readpar_float(file_parameter, 'dz_max', dz_max, 1.5*dz)
        else
            dz_max = dz
        end if

    end subroutine


#ifdef _dim2_

    !
    !> Get model or source parameters from meta arrays by name
    !
    function get_model(name, default_value) result(v)

        character(len=*), intent(in) :: name
        real, intent(in), optional :: default_value
        real, allocatable, dimension(:, :) :: v

        select case (name)

            case ('mt', 'stf')

                if (any(name == model_name)) then
                    v = get_meta_array_core(model_m, name)
                else if (any(name == model_name_aux)) then
                    v = get_meta_array_core(model_aux, name)
                else
                    v = zeros(nc_mt, ns)
                end if

            case default

                if (any(name == model_name)) then
                    v = get_meta_array_core(model_m, name)
                    v = v(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend)
                else if (any(name == model_name_aux)) then
                    v = get_meta_array_core(model_aux, name)
                    v = v(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend)
                else
                    if (.not. present(default_value)) then
                        if (rankid_group == 0) then
                            write(error_unit, *) ' <get_model> Warning: Default_value for '//tidy(name) &
                                //' not given; set to zero. '
                        end if
                    end if
                    v = zeros(shot_nz, shot_nx)
                    if (present(default_value)) then
                        v = v + default_value
                    end if
                end if

        end select

    end function

    !
    !> Get model or source parameters for inversion
    !
    subroutine prepare_model_single_parameter(w, name, file_w, const, source, update)

        real, allocatable, dimension(:, :), intent(inout) :: w
        character(len=*), intent(in) :: name
        character(len=*), intent(in), optional :: file_w
        real, intent(in), optional :: const
        real, dimension(:, :), intent(in), optional :: source
        logical, intent(in), optional :: update

        logical :: update_this_model
        character(len=1024) :: file_update
        real, allocatable, dimension(:, :) :: wt

        select case (name)

            case ('mt', 'stf')

                w = zeros(nc_mt, ns)

                ! If read in or assign const value
                if (present(file_w) .and. file_w /= '') then
                    call input_array(w, file_w)
                else
                    if (present(const)) then
                        w = const
                    end if
                    if (file_w == '') then
                        if (rankid == 0) then
                            call warn(' Warning: Moodel '//tidy(name)//' is empty. ')
                        end if
                    end if
                end if

            case default

                w = zeros(nz, nx)

                ! If read in or assign const value
                if (present(file_w) .and. file_w /= '') then

                    if (require_model_interp) then
                        ! If it is necessary to resample

                        call alloc_array(wt, [1, nz0, 1, nx0])
                        call input_array(wt, file_w)
                        w = interp(wt, [nz0, nx0], [dz0, dx0], [oz0, ox0], &
                            [nz, nx], [dz, dx], [oz, ox], ['linear', 'linear'])
                        deallocate (wt)

                    else
                        ! If resampling is not required
                        call input_array(w, file_w)

                    end if

                else
                    if (present(const)) then
                        w = const
                    end if
                    if (file_w == '') then
                        if (rankid == 0) then
                            call warn(' Warning: Moodel '//tidy(name)//' is empty. ')
                        end if
                    end if
                end if

        end select

        ! If the source is given
        if (present(source)) then
            w = source
        end if

        ! If the inversion starts from certain iteration other than 1
        if (present(update)) then
            if (.not.update) then
                update_this_model = .false.
            else
                update_this_model = .true.
            end if
        else
            update_this_model = .true.
        end if

        if (resume_from_iter > 1 .and. update_this_model) then
            file_update = tidy(dir_working)//'/iteration_'// &
                num2str(resume_from_iter - 1)//'/model/updated_'//tidy(name)//'.bin'
            if (file_exists(file_update)) then
                call input_array(w, file_update)
            else
                if (rankid == 0) then
                    call warn(' Error: Updated model '//tidy(name)//' not found. Exiting. ')
                end if
                stop
            end if
        end if

    end subroutine

#endif


#ifdef _dim3_

    function get_model(name, default_value) result(v)

        character(len=*), intent(in) :: name
        real, intent(in), optional :: default_value
        real, allocatable, dimension(:, :, :) :: v

        select case (name)

            case ('mt', 'stf')

                if (any(name == model_name)) then
                    v = get_meta_array_core(model_m, name)
                else if (any(name == model_name_aux)) then
                    v = get_meta_array_core(model_aux, name)
                else
                    v = zeros(nc_mt, ns, 1)
                end if

            case default

                if (any(name == model_name)) then
                    v = get_meta_array_core(model_m, name)
                    v = v(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)
                else if (any(name == model_name_aux)) then
                    v = get_meta_array_core(model_aux, name)
                    v = v(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend)
                else
                    if (.not. present(default_value)) then
                        if (rankid_group == 0) then
                            write(error_unit, *) ' <get_model> Warning: default_value for '//tidy(name) &
                                //' is not given; setting to zero. '
                        end if
                    end if
                    v = zeros(shot_nz, shot_ny, shot_nx)
                    if (present(default_value)) then
                        v = v + default_value
                    end if
                end if

        end select

    end function

    subroutine prepare_model_single_parameter(w, name, file_w, const, source, update)

        real, allocatable, dimension(:, :, :), intent(inout) :: w
        character(len=*), intent(in) :: name
        character(len=*), intent(in), optional :: file_w
        real, intent(in), optional :: const
        real, dimension(:, :, :), intent(in), optional :: source
        logical, intent(in), optional :: update

        logical :: update_this_model
        character(len=1024) :: file_update
        real, allocatable, dimension(:, :, :) :: wt

        integer, dimension(1:3) :: mdim

        if (rankid == 0) then

            select case (name)

                case ('mt', 'stf')

                    w = zeros(nc_mt, ns, 1)

                    ! If read in or assign const value
                    if (present(file_w) .and. file_w /= '') then
                        call input_array(w, file_w)
                    else
                        if (present(const)) then
                            w = const
                        end if
                        if (file_w == '') then
                            if (rankid == 0) then
                                call warn(' Warning: Moodel '//tidy(name)//' is empty. ')
                            end if
                        end if
                    end if


                case default

                    w = zeros(nz, ny, nx)

                    ! If read in or assign const value
                    if (present(file_w) .and. file_w /= '') then

                        if (require_model_interp) then
                            ! If it is necessary to resample

                            call alloc_array(wt, [1, nz0, 1, ny0, 1, nx0])
                            call input_array(wt, file_w)
                            w = interp(wt, [nz0, ny0, nx0], [dz0, dy0, dx0], [oz0, oy0, ox0], &
                                [nz, ny, nx], [dz, dy, dx], [oz, oy, ox], &
                                ['linear', 'linear', 'linear'])
                            deallocate (wt)

                        else
                            ! If resampling is not required
                            call input_array(w, file_w)

                        end if

                    else
                        if (present(const)) then
                            w = const
                        end if
                        if (file_w == '') then
                            if (rankid == 0) then
                                call warn(' Warning: Parameter '//tidy(name)//' is empty. ')
                            end if
                        end if
                    end if

            end select

            ! if assign a given source
            if (present(source)) then
                w = source
            end if

            ! if the inversion starts from certain iteration other than 1
            if (present(update)) then
                if (.not.update) then
                    update_this_model = .false.
                else
                    update_this_model = .true.
                end if
            else
                update_this_model = .true.
            end if

            if (resume_from_iter > 1 .and. update_this_model) then
                file_update = tidy(dir_working)//'/iteration_'// &
                    num2str(resume_from_iter - 1)//'/model/updated_'//tidy(name)//'.bin'
                if (file_exists(file_update)) then
                    call input_array(w, file_update)
                else
                    if (rankid == 0) then
                        call warn(' Error: Updated parameter model '//tidy(name)//' not found. Exiting. ')
                    end if
                    stop
                end if
            end if

            mdim = [size(w, 1), size(w, 2), size(w, 3)]

        end if

        if (nrank > 1) then
            call mpibarrier
            call bcast_array(mdim)
            call mpibarrier
            if (rankid /= 0) then
                call alloc_array(w, [1, mdim(1), 1, mdim(2), 1, mdim(3)])
            end if
            call mpibarrier
            call bcast_array(w)
        end if

    end subroutine

#endif

end module
