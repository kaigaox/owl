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


program forward

    use mod_parameters
    use mod_data_processing
    use mod_utility
    use mod_model
    use mod_forward_modeling

    implicit none

    call mpistart

    if (command_argument_count() == 0) then
        if (rankid == 0) then
            call warn('')
            call warn(date_time_compact()//' Error: Parameter file not found. Exiting. ')
            call warn('')
        end if
        call mpibarrier
        call mpistop
    end if

    if (rankid == 0) then
        call warn('')
        call warn(tile('=', 80))
        call warn(center_substring('Forward modeling begins', 80))
        call warn('')
        call print_date_time
        call warn('')
    end if

    call read_parameters

    call set_regular_space

    call load_geometry

    !    ! Initialize source encoding
    !    if (yn_source_encoding) then
    !        call init_source_encoding
    !    end if

    call mpistart_group

    ! divide shots
    call divide_shots

    ! Prepare models
    call prepare_model

    call mpibarrier

    call forward_modeling

    call mpibarrier

    if (rankid == 0) then
        call warn('')
        call print_date_time
        call warn('')
        call warn(center_substring('Forward modeling completed', 80))
        call warn(tile('=', 80))
        call warn('')
    end if

    call mpiend_group
    call mpiend

end program
