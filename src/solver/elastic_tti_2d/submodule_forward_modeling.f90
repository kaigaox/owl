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


submodule(elastic_tti_2d) elastic_tti_2d_forward_modeling

    use libflit
    use elastic_tti_2d_vars
    use elastic_tti_2d_boundary_saving
    use elastic_tti_2d_wavefield
    use elastic_tti_2d_cfspml

    implicit none

contains

    module subroutine compute_forward(this)

        class(wave_solver_elastic_tti_2d), intent(inout) :: this

        integer :: l, ir, irx, irz, rgx, rgz, rhx, rhz
        integer :: i, j, t
        real :: amp1, amp2
        real, allocatable, dimension(:) :: rec_slopex
        real, allocatable, dimension(:, :) :: datax, dataz

        call prepare_modeling(this)
        call compute_cfspml_damping_coef
        call alloc_forward_wavefield

        call this%seis_vx%init(nt=nt, dt=dt, nr=sgmtr%nr)
        call this%seis_vz%init(nt=nt, dt=dt, nr=sgmtr%nr)

        ! Prepare for boundary saving
        if (yn_reconstruct) then
            call prepare_boundary_saving
            call open_boundary_saving
        end if

        l = 1
        do t = 1, nt

            ! Save boundary wavefield
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
                    vx_hxiz(:, 1 - j) = vx_hxiz(:, 1 + j)
                    vz_hxiz(:, 1 - j) = vz_hxiz(:, 1 + j)
                    vx_ixhz(:, 2 - j) = vx_ixhz(:, 1 + j)
                    vz_ixhz(:, 2 - j) = vz_ixhz(:, 1 + j)
                end do
                !$omp end parallel do
            end if

            !$omp parallel do private(ir, irx, irz, rgx, rgz, rhx, rhz, amp1, amp2) schedule(auto)
            do ir = 1, sgmtr%nr
                if (sgmtr%recr(ir)%weight /= 0) then

                    rgx = sgmtr%recr(ir)%gx
                    rgz = sgmtr%recr(ir)%gz
                    rhx = sgmtr%recr(ir)%hx
                    rhz = sgmtr%recr(ir)%hz

                    do irz = -nkw, nkw
                        do irx = -nkw, nkw
                            amp1 = sgmtr%recr(ir)%interp_hx(irx) &
                                *sgmtr%recr(ir)%interp_iz(irz) &
                                *sgmtr%recr(ir)%weight*0.5
                            amp2 = sgmtr%recr(ir)%interp_ix(irx) &
                                *sgmtr%recr(ir)%interp_hz(irz) &
                                *sgmtr%recr(ir)%weight*0.5
                            this%seis_vx%trace(ir)%data(t) = &
                                this%seis_vx%trace(ir)%data(t) &
                                + vx_hxiz(rhx + irx, rgz + irz)*amp1 &
                                + vx_ixhz(rgx + irx, rhz + irz)*amp2
                        end do
                    end do

                    do irz = -nkw, nkw
                        do irx = -nkw, nkw
                            amp1 = sgmtr%recr(ir)%interp_hx(irx) &
                                *sgmtr%recr(ir)%interp_iz(irz) &
                                *sgmtr%recr(ir)%weight*0.5
                            amp2 = sgmtr%recr(ir)%interp_ix(irx) &
                                *sgmtr%recr(ir)%interp_hz(irz) &
                                *sgmtr%recr(ir)%weight*0.5
                            this%seis_vz%trace(ir)%data(t) = &
                                this%seis_vz%trace(ir)%data(t) &
                                + vz_hxiz(rhx + irx, rgz + irz)*amp1 &
                                + vz_ixhz(rgx + irx, rhz + irz)*amp2
                        end do
                    end do

                end if

            end do
            !$omp end parallel do

            ! Record wavefield snapshot if necessary
            if (np /= 0 .and. l <= np) then
                if (t - 1 == nint(snaps(l)/dt)) then

                    call alloc_array(snapvx, [1, nx, 1, nz], pad=pml)
                    call alloc_array(snapvz, [1, nx, 1, nz], pad=pml)

                    if (yn_free_surface) then

                        !$omp parallel do private(i, j) collapse(2)
                        do j = -pml + 1, nz + pml
                            do i = -pml + 1, nx + pml
                                snapvx(i, j) = 0.5*(0.5*sum(vx_hxiz(i:i + 1, j)) + 0.5*sum(vx_ixhz(i, j:j + 1)))
                                snapvz(i, j) = 0.5*(0.5*sum(vz_hxiz(i:i + 1, j)) + 0.5*sum(vz_ixhz(i, j:j + 1)))
                            end do
                        end do
                        !$omp end parallel do

                        open(3, file=tidy(dir_snapshot)//'/shot_' &
                            //num2str(sgmtr%id) &
                            //'_forward_wavefield_' &
                            //num2str(l)//'.txt')
                        do i = -pml + 1, nx + pml
                            do j = 1, nz + pml
                                write(3, '(4es)') (i - 1)*dx + ox, zz_i(i, j), snapvx(i, j), snapvz(i, j)
                            end do
                        end do
                        close(3)

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
                        do j = -pml + 1, nz + pml - 1
                            do i = -pml + 1, nx + pml - 1
                                snapvx(i, j) = 0.5*(0.5*sum(vx_hxiz(i:i + 1, j)) + 0.5*sum(vx_ixhz(i, j:j + 1)))
                                snapvz(i, j) = 0.5*(0.5*sum(vz_hxiz(i:i + 1, j)) + 0.5*sum(vz_ixhz(i, j:j + 1)))
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

            ! Update wavefield
            if (yn_free_surface) then
                call update_wavefield_free_surface(dt, &
                    stressxx_ixiz, stresszz_ixiz, stressxz_ixiz, &
                    stressxx_hxhz, stresszz_hxhz, stressxz_hxhz, &
                    vx_ixhz, vz_ixhz, &
                    vx_hxiz, vz_hxiz, &
                    memory_pdxvx_hxhz, &
                    memory_pdxvz_hxhz, &
                    memory_pdxvx_ixiz, &
                    memory_pdxvz_ixiz, &
                    memory_pdxxx_hxiz, &
                    memory_pdxxz_hxiz, &
                    memory_pdxxx_ixhz, &
                    memory_pdxxz_ixhz, &
                    memory_pdzvx_hxhz, &
                    memory_pdzvz_hxhz, &
                    memory_pdzvx_ixiz, &
                    memory_pdzvz_ixiz, &
                    memory_pdzxx_hxiz, &
                    memory_pdzxz_hxiz, &
                    memory_pdzzz_hxiz, &
                    memory_pdzxx_ixhz, &
                    memory_pdzxz_ixhz, &
                    memory_pdzzz_ixhz)
            else
                call update_wavefield(dt, &
                    stressxx_ixiz, stresszz_ixiz, stressxz_ixiz, &
                    stressxx_hxhz, stresszz_hxhz, stressxz_hxhz, &
                    vx_ixhz, vz_ixhz, &
                    vx_hxiz, vz_hxiz, &
                    memory_pdxvx_hxhz, &
                    memory_pdxvz_hxhz, &
                    memory_pdxvx_ixiz, &
                    memory_pdxvz_ixiz, &
                    memory_pdxxx_hxiz, &
                    memory_pdxxz_hxiz, &
                    memory_pdxxx_ixhz, &
                    memory_pdxxz_ixhz, &
                    memory_pdzvx_hxhz, &
                    memory_pdzvz_hxhz, &
                    memory_pdzvx_ixiz, &
                    memory_pdzvz_ixiz, &
                    memory_pdzxz_hxiz, &
                    memory_pdzzz_hxiz, &
                    memory_pdzxz_ixhz, &
                    memory_pdzzz_ixhz)
            end if

            ! Source term
            call add_source(t)

            if (verbose .and. (mod(t, max(nint(nt/10.0), 1)) == 0 .or. t == 1 .or. t == nt)) then
                call warn(date_time_compact()//' >> Shot '//num2str(sgmtr%id) &
                    //' forward modeling step '//num2str(t)//' of '//num2str(nt))
                if (any(isnan(vx_ixhz)) .or. any(isnan(vz_ixhz))) then
                    call warn(date_time_compact()//' >> Vx, Vz contain NaN!')
                    stop
                else
                    call warn(date_time_compact()//' >> Vx, Vz value range = ')
                    call warn(date_time_compact()//'      '//num2str(min(minval(vx_hxiz), minval(vx_ixhz)), '(es)') &
                        //' ~ '//num2str(max(maxval(vx_hxiz), maxval(vx_ixhz)), '(es)'))
                    call warn(date_time_compact()//'      '//num2str(min(minval(vz_hxiz), minval(vz_ixhz)), '(es)') &
                        //' ~ '//num2str(max(maxval(vz_hxiz), maxval(vz_ixhz)), '(es)'))
                end if
            end if

        end do

        ! Close boundary saving
        if (yn_reconstruct) then
            call close_boundary_saving
        end if

        ! Rotate receivers if necessary
        if (this%receiver_vertical_to_surface) then

            rec_slopex = ginterp(topox, slopex_i, sgmtr%recr(:)%x, 'cubic')
            rec_slopex = atan(rec_slopex)

            datax = this%seis_vx%to_array()
            dataz = this%seis_vz%to_array()

            !$omp parallel do private(i)
            do i = 1, sgmtr%nr

                this%seis_vx%trace(i)%data = cos(rec_slopex(i))*datax(:, i) - sin(rec_slopex(i))*dataz(:, i)
                this%seis_vz%trace(i)%data = sin(rec_slopex(i))*datax(:, i) + cos(rec_slopex(i))*dataz(:, i)

            end do
            !$omp end parallel do

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
