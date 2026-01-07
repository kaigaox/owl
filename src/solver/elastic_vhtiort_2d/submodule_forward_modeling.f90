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


submodule(elastic_vhtiort_2d) elastic_vhtiort_2d_forward_modeling

    use libflit
    use elastic_vhtiort_2d_vars
    use elastic_vhtiort_2d_boundary_saving
    use elastic_vhtiort_2d_wavefield
    use elastic_vhtiort_2d_cfspml

    implicit none

contains

    module subroutine compute_forward(this)

        class(wave_solver_elastic_vhtiort_2d), intent(inout) :: this

        integer :: l, ir, irx, irz, rgx, rgz, rhx, rhz
        integer :: i, j, t

        call prepare_modeling(this)
        call compute_cfspml_damping_coef
        call alloc_forward_wavefield

        if (yn_compx) then
            call this%seis_vx%init(nt=nt, dt=dt, nr=sgmtr%nr)
        end if
        if (yn_compz) then
            call this%seis_vz%init(nt=nt, dt=dt, nr=sgmtr%nr)
        end if

        if (yn_reconstruct) then
            call prepare_boundary_saving
            call open_boundary_saving
        end if

        ! start modeling
        l = 1
        do t = 1, nt

            ! save boundary wavefield
            if (yn_reconstruct) then
                call save_boundary_wavefield(t)
                if (t == nt) then
                    call output_final_step_wavefield
                end if
            end if

            if (yn_free_surface) then
                ! This is to ensure near surface receivers get correct amplitude
                !$omp parallel do private(j) schedule(auto)
                do j = 1, fdhalf
                    vx(:, 1 - j) = vx(:, 1 + j)
                    vz(:, 2 - j) = vz(:, 1 + j)
                end do
                !$omp end parallel do
            end if

            !$omp parallel do private(ir, irx, irz, rgx, rgz, rhx, rhz) schedule(auto)
            do ir = 1, sgmtr%nr
                if (sgmtr%recr(ir)%weight /= 0) then

                    rgx = sgmtr%recr(ir)%gx
                    rgz = sgmtr%recr(ir)%gz
                    rhx = sgmtr%recr(ir)%hx
                    rhz = sgmtr%recr(ir)%hz

                    if (yn_compx) then
                        do irz = -nkw, nkw
                            do irx = -nkw, nkw
                                this%seis_vx%trace(ir)%data(t) = &
                                    this%seis_vx%trace(ir)%data(t) + &
                                    vx(rhx + irx, rgz + irz) &
                                    *sgmtr%recr(ir)%interp_hx(irx) &
                                    *sgmtr%recr(ir)%interp_iz(irz) &
                                    *sgmtr%recr(ir)%weight
                            end do
                        end do
                    end if

                    if (yn_compz) then
                        do irz = -nkw, nkw
                            do irx = -nkw, nkw
                                this%seis_vz%trace(ir)%data(t) = &
                                    this%seis_vz%trace(ir)%data(t) + &
                                    vz(rgx + irx, rhz + irz) &
                                    *sgmtr%recr(ir)%interp_ix(irx) &
                                    *sgmtr%recr(ir)%interp_hz(irz) &
                                    *sgmtr%recr(ir)%weight
                            end do
                        end do
                    end if

                end if
            end do
            !$omp end parallel do

            ! Output wavefield snapshot if necessary
            if (np /= 0 .and. l <= np) then
                if (t - 1 == nint(snaps(l)/dt)) then

                    call alloc_array(snapvx, [1, nx, 1, nz], pad=pml)
                    call alloc_array(snapvz, [1, nx, 1, nz], pad=pml)

                    if (yn_free_surface) then

                        !$omp parallel do private(i, j) collapse(2)
                        do j = -pml + 1, nz + pml
                            do i = 1, nx
                                snapvx(i, j) = 0.5*sum(vx(i:i + 1, j))
                                snapvz(i, j) = 0.5*sum(vz(i, j:j + 1))
                            end do
                        end do
                        !$omp end parallel do

                        call map_irregular_to_regular(snapvx, this, [1, this%nx, 1, this%nz])
                        call map_irregular_to_regular(snapvz, this, [1, this%nx, 1, this%nz])

                        call output_array(snapvx, tidy(dir_snapshot)//'/shot_' &
                            //num2str(sgmtr%id) &
                            //'_forward_wavefield_x_' &
                            //num2str(l)//'.bin', transp=.true.)
                        call output_array(snapvz, tidy(dir_snapshot)//'/shot_' &
                            //num2str(sgmtr%id) &
                            //'_forward_wavefield_z_' &
                            //num2str(l)//'.bin', transp=.true.)

                    else

                        !$omp parallel do private(i, j) collapse(2)
                        do j = 1, nz
                            do i = 1, nx
                                snapvx(i, j) = 0.5*sum(vx(i:i + 1, j))
                                snapvz(i, j) = 0.5*sum(vz(i, j:j + 1))
                            end do
                        end do
                        !$omp end parallel do

                        call output_array(snapvx(1:nx, 1:nz), tidy(dir_snapshot)//'/shot_' &
                            //num2str(sgmtr%id) &
                            //'_forward_wavefield_x_' &
                            //num2str(l)//'.bin', transp=.true.)
                        call output_array(snapvz(1:nx, 1:nz), tidy(dir_snapshot)//'/shot_' &
                            //num2str(sgmtr%id) &
                            //'_forward_wavefield_z_' &
                            //num2str(l)//'.bin', transp=.true.)

                    end if

                    l = l + 1
                end if
            end if

            ! Wavefield update
            if (yn_free_surface) then
                call update_wavefield_free_surface(dt, &
                    stressxx, stresszz, stressxz, vx, vz, &
                    memory_pdxvx, memory_pdzvx, memory_pdxvz, memory_pdzvz, &
                    memory_pdxxx, memory_pdzxz, memory_pdxxz, memory_pdzzz)
            else
                call update_wavefield(dt, &
                    stressxx, stresszz, stressxz, vx, vz, &
                    memory_pdxvx, memory_pdzvx, memory_pdxvz, memory_pdzvz, &
                    memory_pdxxx, memory_pdzxz, memory_pdxxz, memory_pdzzz)
            end if

            ! Source term
            call add_source(t)

            if (verbose .and. (mod(t, max(nint(nt/10.0), 1)) == 0 .or. t == 1 .or. t == nt)) then
                call warn(date_time_compact()//' >> Shot '//num2str(sgmtr%id) &
                    //' forward modeling step '//num2str(t)//' of '//num2str(nt))
                if (any(isnan(vx)) .or. any(isnan(vz))) then
                    call warn(date_time_compact()//' >> Vx, Vz contain NaN!')
                    stop
                else
                    call warn(date_time_compact()//' >> Vx, Vz value range = ')
                    call warn(date_time_compact()//'      '//num2str(minval(vx), '(es)') &
                        //' ~ '//num2str(maxval(vx), '(es)'))
                    call warn(date_time_compact()//'      '//num2str(minval(vz), '(es)') &
                        //' ~ '//num2str(maxval(vz), '(es)'))
                end if
            end if

        end do

        ! Close boundary saving
        if (yn_reconstruct) then
            call close_boundary_saving
        end if

        ! Output data
        if (yn_compx) then
            call this%seis_vx%resamp(nnt=data_nt, ddt=data_dt)
            call this%seis_vx%output(tidy(dir_synthetic)//'/shot_'//num2str(sgmtr%id)//'_seismogram_x.su')
        end if

        if (yn_compz) then
            call this%seis_vz%resamp(nnt=data_nt, ddt=data_dt)
            call this%seis_vz%output(tidy(dir_synthetic)//'/shot_'//num2str(sgmtr%id)//'_seismogram_z.su')
        end if

        call mpibarrier_group

    end subroutine compute_forward

end submodule
