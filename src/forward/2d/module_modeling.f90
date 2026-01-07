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


module mod_forward_modeling

    use mod_parameters
    use mod_data_processing
    use mod_utility
    use mod_model
    use acoustic_iso_2d
    use elastic_vhtiort_2d
    use elastic_tti_2d

    implicit none

contains

    subroutine prepare_model

        integer :: i
        character(len=1024) :: fname

        allocate(model_m(1:nmodel))

        do i = 1, nmodel
            call readpar_string(file_parameter, 'file_'//tidy(model_name(i)), fname, '')
            call prepare_model_single_parameter(model_m(i)%array, model_name(i), fname)
            model_m(i)%name = model_name(i)
        end do

    end subroutine prepare_model

    subroutine forward_modeling

        integer :: i

        type(wave_solver_acoustic_iso_2d) :: solver_acoustic_iso
        type(wave_solver_elastic_vhtiort_2d) :: solver_elastic_vhtiort
        type(wave_solver_elastic_tti_2d) :: solver_elastic_tti

        ! Show the medium parameter statistics
        if (rankid == 0) then
            do i = 1, nmodel
                call plot_histogram(model_m(i)%array, &
                    label=date_time_compact()//' '//tidy(model_name(i))//' distribution ')
            end do
        end if
        call mpibarrier

        do ishot = shot_in_group(groupid, 1), shot_in_group(groupid, 2)

            call set_adaptive_range(gmtr(ishot))

            select case (which_medium)

                case ('acoustic-iso')
                    solver_acoustic_iso%nx = shot_nx
                    solver_acoustic_iso%nz = shot_nz
                    solver_acoustic_iso%dx = dx
                    solver_acoustic_iso%dz = dz
                    solver_acoustic_iso%ox = shot_xbeg
                    solver_acoustic_iso%oz = shot_zbeg
                    solver_acoustic_iso%dt = dt
                    solver_acoustic_iso%tmax = tmax
                    solver_acoustic_iso%data_dt = data_dt
                    solver_acoustic_iso%data_tmax = data_tmax
                    solver_acoustic_iso%dir_synthetic = tidy(dir_synthetic)
                    if (sum(snaps) > 0) then
                        solver_acoustic_iso%dir_snapshot = tidy(dir_snapshot)
                        solver_acoustic_iso%snaps = regspace(snaps(1), snaps(2), snaps(3))
                    end if
                    solver_acoustic_iso%pml = pml
                    solver_acoustic_iso%free_surface = yn_free_surface
                    solver_acoustic_iso%verbose = verbose
                    solver_acoustic_iso%gmtr = gmtr(ishot)
                    solver_acoustic_iso%vp = get_model('vp')
                    solver_acoustic_iso%rho = get_model('rho', 1.0)

                    call solver_acoustic_iso%forward

                case ('elastic-iso', 'elastic-vhtiort')
                    solver_elastic_vhtiort%nx = shot_nx
                    solver_elastic_vhtiort%nz = shot_nz
                    solver_elastic_vhtiort%dx = dx
                    solver_elastic_vhtiort%dz = dz
                    solver_elastic_vhtiort%ox = shot_xbeg
                    solver_elastic_vhtiort%oz = shot_zbeg
                    solver_elastic_vhtiort%dt = dt
                    solver_elastic_vhtiort%tmax = tmax
                    solver_elastic_vhtiort%data_dt = data_dt
                    solver_elastic_vhtiort%data_tmax = data_tmax
                    solver_elastic_vhtiort%dir_synthetic = tidy(dir_synthetic)
                    if (sum(snaps) > 0) then
                        solver_elastic_vhtiort%dir_snapshot = tidy(dir_snapshot)
                        solver_elastic_vhtiort%snaps = regspace(snaps(1), snaps(2), snaps(3))
                    end if
                    solver_elastic_vhtiort%pml = pml
                    solver_elastic_vhtiort%free_surface_dz_refine = free_surface_dz_refine
                    solver_elastic_vhtiort%dz_max = dz_max
                    solver_elastic_vhtiort%free_surface = yn_free_surface
                    solver_elastic_vhtiort%verbose = verbose
                    solver_elastic_vhtiort%gmtr = gmtr(ishot)
                    solver_elastic_vhtiort%compx = yn_compx
                    solver_elastic_vhtiort%compz = yn_compz

                    select case (aniso_param)

                        case ('iso')
                            solver_elastic_vhtiort%anisotropy_type = 'iso'
                            solver_elastic_vhtiort%vp = get_model('vp')
                            solver_elastic_vhtiort%vs = get_model('vs')
                            solver_elastic_vhtiort%rho = get_model('rho', 1.0)

                        case ('thomsen')
                            solver_elastic_vhtiort%anisotropy_type = 'thomsen'
                            solver_elastic_vhtiort%vp = get_model('vp')
                            solver_elastic_vhtiort%vs = get_model('vs')
                            solver_elastic_vhtiort%tieps = get_model('epsilon', 0.0)
                            solver_elastic_vhtiort%tidel = get_model('delta', 0.0)
                            solver_elastic_vhtiort%tithe = get_model('theta', 0.0)
                            solver_elastic_vhtiort%rho = get_model('rho', 1.0)

                        case ('a-t')
                            solver_elastic_vhtiort%anisotropy_type = 'a-t'
                            solver_elastic_vhtiort%vp = get_model('vp')
                            solver_elastic_vhtiort%vs = get_model('vs')
                            solver_elastic_vhtiort%tieps = get_model('epsilon', 0.0)
                            solver_elastic_vhtiort%tieta = get_model('eta', 0.0)
                            solver_elastic_vhtiort%tithe = get_model('theta', 0.0)
                            solver_elastic_vhtiort%rho = get_model('rho', 1.0)

                        case ('cij')
                            solver_elastic_vhtiort%anisotropy_type = 'cij'
                            solver_elastic_vhtiort%c11 = get_model('c11')
                            solver_elastic_vhtiort%c13 = get_model('c13')
                            solver_elastic_vhtiort%c33 = get_model('c33')
                            solver_elastic_vhtiort%c55 = get_model('c55')
                            solver_elastic_vhtiort%rho = get_model('rho', 1.0)

                    end select

                    call solver_elastic_vhtiort%forward

                case ('elastic-tti')
                    solver_elastic_tti%nx = shot_nx
                    solver_elastic_tti%nz = shot_nz
                    solver_elastic_tti%dx = dx
                    solver_elastic_tti%dz = dz
                    solver_elastic_tti%ox = shot_xbeg
                    solver_elastic_tti%oz = shot_zbeg
                    solver_elastic_tti%dt = dt
                    solver_elastic_tti%tmax = tmax
                    solver_elastic_tti%data_dt = data_dt
                    solver_elastic_tti%data_tmax = data_tmax
                    solver_elastic_tti%dir_synthetic = tidy(dir_synthetic)
                    if (sum(snaps) > 0) then
                        solver_elastic_tti%dir_snapshot = tidy(dir_snapshot)
                        solver_elastic_tti%snaps = regspace(snaps(1), snaps(2), snaps(3))
                    end if
                    solver_elastic_tti%pml = pml
                    solver_elastic_tti%free_surface_dz_refine = free_surface_dz_refine
                    solver_elastic_tti%dz_max = dz_max
                    solver_elastic_tti%free_surface = yn_free_surface
                    solver_elastic_tti%file_topo = tidy(file_topo)
                    solver_elastic_tti%topo_interp = tidy(topo_interp)
                    solver_elastic_tti%measure_source_depth_from_surface = measure_source_depth_from_surface
                    solver_elastic_tti%measure_receiver_depth_from_surface = measure_receiver_depth_from_surface
                    solver_elastic_tti%source_vertical_to_surface = source_vertical_to_surface
                    solver_elastic_tti%receiver_vertical_to_surface = receiver_vertical_to_surface
                    solver_elastic_tti%verbose = verbose
                    solver_elastic_tti%gmtr = gmtr(ishot)
                    solver_elastic_tti%compx = yn_compx
                    solver_elastic_tti%compz = yn_compz
                    solver_elastic_tti%save_mesh = yn_save_mesh

                    select case (aniso_param)

                        case ('iso')
                            solver_elastic_tti%anisotropy_type = 'iso'
                            solver_elastic_tti%vp = get_model('vp')
                            solver_elastic_tti%vs = get_model('vs')
                            solver_elastic_tti%rho = get_model('rho', 1.0)

                        case ('thomsen')
                            solver_elastic_tti%anisotropy_type = 'thomsen'
                            solver_elastic_tti%vp = get_model('vp')
                            solver_elastic_tti%vs = get_model('vs')
                            solver_elastic_tti%tieps = get_model('epsilon', 0.0)
                            solver_elastic_tti%tidel = get_model('delta', 0.0)
                            solver_elastic_tti%tithe = get_model('theta', 0.0)
                            solver_elastic_tti%rho = get_model('rho', 1.0)

                        case ('a-t')
                            solver_elastic_tti%anisotropy_type = 'a-t'
                            solver_elastic_tti%vp = get_model('vp')
                            solver_elastic_tti%vs = get_model('vs')
                            solver_elastic_tti%tieps = get_model('epsilon', 0.0)
                            solver_elastic_tti%tieta = get_model('eta', 0.0)
                            solver_elastic_tti%tithe = get_model('theta', 0.0)
                            solver_elastic_tti%rho = get_model('rho', 1.0)

                        case ('cij')
                            solver_elastic_tti%anisotropy_type = 'cij'
                            solver_elastic_tti%c11 = get_model('c11')
                            solver_elastic_tti%c13 = get_model('c13')
                            solver_elastic_tti%c15 = get_model('c15')
                            solver_elastic_tti%c33 = get_model('c33')
                            solver_elastic_tti%c35 = get_model('c35')
                            solver_elastic_tti%c55 = get_model('c55')
                            solver_elastic_tti%rho = get_model('rho', 1.0)

                    end select

                    call solver_elastic_tti%forward

            end select

            call warn(date_time_compact()//' Shot '//num2str(set_srcid(ishot))//' forward modeling completed.')

            call process_synthetic(ishot)

        end do

        call mpibarrier

    end subroutine forward_modeling

end module mod_forward_modeling
