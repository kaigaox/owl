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


module mod_grid

    use libflit

    implicit none

    ! 1D grid
    type, public :: grid1

        real    :: o1 = 0.0
        real    :: d1 = 1.0
        integer :: n1
        real    :: l1
        character(len=64) :: label1 = 'Axis 1'
        character(len=64) :: unit1 = 'unitless'
        character(len=80), dimension(1:50) :: header
        logical :: initialized = .false.
        real, allocatable, dimension(:) :: array

    contains
        procedure :: init => init_grid_1d
        procedure :: from_array => init_grid_from_1d
        procedure :: output => save_grid_1d
        procedure :: input => read_grid_1d
        procedure :: stdin => stdin_grid_1d
        procedure :: stdout => stdout_grid_1d
        procedure :: to_array => extract_grid_1d
        procedure :: free => free_grid_1d

    end type grid1

    ! 2D grid
    type, public ::  grid2

        real    :: o1 = 0.0, o2 = 0.0
        real    :: d1 = 1.0, d2 = 1.0
        integer :: n1, n2
        real    :: l1, l2
        character(len=64) :: label1 = 'Axis 1', label2 = 'Axis 2'
        character(len=64) :: unit1 = 'unitless', unit2 = 'unitless'
        character(len=80), dimension(1:50) :: header
        logical :: initialized = .false.
        real, allocatable, dimension(:, :) :: array

    contains
        procedure :: init => init_grid_2d
        procedure :: from_array => init_grid_from_2d
        procedure :: output => save_grid_2d
        procedure :: input => read_grid_2d
        procedure :: stdin => stdin_grid_2d
        procedure :: stdout => stdout_grid_2d
        procedure :: to_array => extract_grid_2d
        procedure :: free => free_grid_2d

    end type grid2

    ! 3D grid
    type, public ::  grid3

        real    :: o1 = 0.0, o2 = 0.0, o3 = 0.0
        real    :: d1 = 1.0, d2 = 1.0, d3 = 1.0
        integer :: n1, n2, n3
        real    :: l1, l2, l3
        character(len=64) :: label1 = 'Axis 1', label2 = 'Axis 2', label3 = 'Axis 3'
        character(len=64) :: unit1 = 'unitless', unit2 = 'unitless', unit3 = 'unitless'
        character(len=80), dimension(1:50) :: header
        logical :: initialized = .false.
        real, allocatable, dimension(:, :, :) :: array

    contains
        procedure :: init => init_grid_3d
        procedure :: from_array => init_grid_from_3d
        procedure :: output => save_grid_3d
        procedure :: input => read_grid_3d
        procedure :: stdin => stdin_grid_3d
        procedure :: stdout => stdout_grid_3d
        procedure :: to_array => extract_grid_3d
        procedure :: free => free_grid_3d

    end type grid3

    !    ! 4D grid
    !    type, public ::  grid4
    !
    !        real    :: o1 = 0.0, o2 = 0.0, o3 = 0.0, o4 = 0.0
    !        real    :: d1 = 1.0, d2 = 1.0, d3 = 1.0, d4 = 1.0
    !        integer :: n1, n2, n3, n4
    !        real    :: l1, l2, l3, l4
    !        character(len=64) :: label1 = 'Axis 1', label2 = 'Axis 2', label3 = 'Axis 3', label4 = 'Axis 4'
    !        character(len=64) :: unit1 = 'unitless', unit2 = 'unitless', unit3 = 'unitless', unit4 = 'unitless'
    !        real :: polar, azimuth
    !
    !
    !    end type grid4

contains

    subroutine init_grid_1d(this, n1, d1, o1, label1, unit1, const)

        class(grid1), intent(inout) :: this
        integer, intent(in) :: n1
        real, intent(in), optional :: d1, o1, const
        character(len=*), intent(in), optional :: label1, unit1

        character(len=80) :: h
        integer :: i

        this%n1 = n1

        if (present(d1)) then
            this%d1 = d1
        end if
        if (present(o1)) then
            this%o1 = o1
        end if
        if (present(label1)) then
            this%label1 = label1
        end if
        if (present(unit1)) then
            this%unit1 = unit1
        end if
        this%l1 = (this%n1 - 1)*this%d1 + this%o1

        call alloc_array(this%array, [1, this%n1])
        if (present(const)) then
            this%array = const
        end if

        this%header = ''
        this%header(1) = '# Start 4000-byte ASCII header '
        this%header(2) = 'n1 = '//num2str(this%n1)
        this%header(3) = 'd1 = '//num2str(this%d1, '(es)')
        this%header(4) = 'o1 = '//num2str(this%o1, '(es)')
        this%header(5) = 'l1 = '//num2str(this%o1 + (this%n1 - 1)*this%d1, '(es)')
        this%header(6) = 'label1 = '//tidy(this%label1)
        this%header(7) = 'unit1 = '//tidy(this%unit1)
        this%header(size(this%header)) = '# End 4000-byte ASCII header'
        do i = 1, size(this%header)
            h = this%header(i)
            h(79:79) = achar(13)
            h(80:80) = achar(10)
            this%header(i) = h
        end do

        this%initialized = .true.

    end subroutine init_grid_1d

    subroutine init_grid_from_1d(this, w)

        class(grid1), intent(inout) :: this
        real, dimension(:), intent(in) :: w

        integer :: i
        character(len=80) :: h

        this%n1 = size(w)

        this%header = ''
        this%header(1) = '# Start 4000-byte ASCII header '
        this%header(2) = 'n1 = '//num2str(this%n1)
        this%header(3) = 'd1 = '//num2str(this%d1, '(es)')
        this%header(4) = 'o1 = '//num2str(this%o1, '(es)')
        this%header(5) = 'l1 = '//num2str(this%o1 + (this%n1 - 1)*this%d1, '(es)')
        this%header(6) = 'label1 = '//tidy(this%label1)
        this%header(7) = 'unit1 = '//tidy(this%unit1)
        this%header(size(this%header)) = '# End 4000-byte ASCII header'
        do i = 1, size(this%header)
            h = this%header(i)
            h(79:79) = achar(13)
            h(80:80) = achar(10)
            this%header(i) = h
        end do

        this%initialized = .true.

        this%array = w

    end subroutine init_grid_from_1d

    subroutine read_grid_1d(this, filename)

        class(grid1), intent(inout) :: this
        character(len=*), intent(in) :: filename

        integer :: funit, i

        funit = nint(rand()*1000)
        open (funit, file=tidy(filename), access='stream', form='unformatted', status='old')
        ! Read header part
        do i = 1, size(this%header)
            read(funit) this%header(i)
        end do

        call parsepar_int(this%header, 'n1', this%n1, int((get_file_size(filename) - 4000.0)/4))
        call parsepar_float(this%header, 'd1', this%d1, 1.0)
        call parsepar_float(this%header, 'o1', this%o1, 0.0)
        this%l1 = (this%n1 - 1)*this%d1 + this%o1
        call parsepar_string(this%header, 'label1', this%label1, 'Axis 1')
        call parsepar_string(this%header, 'unit1', this%unit1, 'unitless')

        this%array = zeros(this%n1)

        ! Read array part
        read (funit, pos=4001) this%array
        close (funit)

    end subroutine read_grid_1d

    subroutine save_grid_1d(this, filename)

        class(grid1), intent(inout) :: this
        character(len=*), intent(in) :: filename

        integer :: funit, i

        funit = nint(rand()*1000)
        open (funit, file=tidy(filename), access='stream', form='unformatted', status='replace')
        ! Write header part as ASCII
        do i = 1, size(this%header)
            write (funit) this%header(i)
        end do
        ! Write array part as raw binary
        write (funit, pos=4001) this%array
        close (funit)

    end subroutine save_grid_1d

    subroutine stdin_grid_1d(this)

        class(grid1), intent(inout) :: this

        integer :: i

        close (input_unit)
        open (input_unit, access='stream', form='unformatted')

        ! Read header part
        do i = 1, size(this%header)
            read(input_unit) this%header(i)
        end do

        ! Parse header
        call parsepar_int(this%header, 'n1', this%n1, int((get_stdin_size() - 4000.0)/4))
        call parsepar_float(this%header, 'd1', this%d1, 1.0)
        call parsepar_float(this%header, 'o1', this%o1, 0.0)
        this%l1 = (this%n1 - 1)*this%d1 + this%o1
        call parsepar_string(this%header, 'label1', this%label1, 'Axis 1')
        call parsepar_string(this%header, 'unit1', this%unit1, 'unitless')
        call alloc_array(this%array, [1, this%n1])

        ! Read array part
        read (input_unit, pos=4001) this%array
        close (input_unit)

    end subroutine stdin_grid_1d

    subroutine stdout_grid_1d(this)

        class(grid1), intent(inout) :: this

        integer :: i

        close (output_unit)
        open (output_unit, access='stream', form='unformatted', status='replace')
        ! Write header part as ASCII
        do i = 1, size(this%header)
            write (output_unit) this%header(i)
        end do
        ! Write array part as raw binary
        write (output_unit, pos=4001) this%array
        close (output_unit)

    end subroutine stdout_grid_1d

    function extract_grid_1d(this) result(w)

        class(grid1), intent(in) :: this

        real, allocatable, dimension(:) :: w

        w = this%array

    end function extract_grid_1d

    subroutine free_grid_1d(this)

        class(grid1), intent(inout) :: this

        if (allocated(this%array)) then
            deallocate(this%array)
        end if

    end subroutine free_grid_1d

    !===========================================================================

    subroutine init_grid_2d(this, n, d, o, label, unit, const)

        class(grid2), intent(inout) :: this
        integer, dimension(:), intent(in) :: n
        real, dimension(:), intent(in), optional :: d, o
        real, intent(in), optional :: const
        character(len=*), dimension(:), intent(in), optional :: label, unit

        character(len=80) :: h
        integer :: i

        call assert(size(n) == 2, 'size(n) /= 2')
        this%n1 = n(1)
        this%n2 = n(2)

        if (present(d)) then
            call assert(size(d) == 2, 'size(d) /= 2')
            this%d1 = d(1)
            this%d2 = d(2)
        end if
        if (present(o)) then
            call assert(size(o) == 2, 'size(d) /= 2')
            this%o1 = o(1)
            this%o2 = o(2)
        end if
        if (present(label)) then
            call assert(size(label) == 2, 'size(label) /= 2')
            this%label1 = label(1)
            this%label2 = label(2)
        end if
        if (present(unit)) then
            call assert(size(unit) == 2, 'size(unit) /= 2')
            this%unit1 = unit(1)
            this%unit2 = unit(2)
        end if
        this%l1 = (this%n1 - 1)*this%d1 + this%o1
        this%l2 = (this%n2 - 1)*this%d2 + this%o2

        call alloc_array(this%array, [1, this%n1, 1, this%n2])
        if (present(const)) then
            this%array = const
        end if

        this%header = ''
        this%header(1) = '# Start 4000-byte ASCII header '
        this%header(2) = 'n1 = '//num2str(this%n1)
        this%header(3) = 'd1 = '//num2str(this%d1, '(es)')
        this%header(4) = 'o1 = '//num2str(this%o1, '(es)')
        this%header(5) = 'l1 = '//num2str(this%o1 + (this%n1 - 1)*this%d1, '(es)')
        this%header(6) = 'label1 = '//tidy(this%label1)
        this%header(7) = 'unit1 = '//tidy(this%unit1)
        this%header(8) = 'n2 = '//num2str(this%n2)
        this%header(9) = 'd2 = '//num2str(this%d2, '(es)')
        this%header(10) = 'o2 = '//num2str(this%o2, '(es)')
        this%header(11) = 'l2 = '//num2str(this%o2 + (this%n2 - 1)*this%d2, '(es)')
        this%header(12) = 'label2 = '//tidy(this%label2)
        this%header(13) = 'unit2 = '//tidy(this%unit2)
        this%header(size(this%header)) = '# End 4000-byte ASCII header'
        do i = 1, size(this%header)
            h = this%header(i)
            h(79:79) = achar(13)
            h(80:80) = achar(10)
            this%header(i) = h
        end do

        this%initialized = .true.

    end subroutine init_grid_2d

    subroutine init_grid_from_2d(this, w)

        class(grid2), intent(inout) :: this
        real, dimension(:, :), intent(in) :: w

        integer :: i
        character(len=80) :: h

        this%n1 = size(w, 1)
        this%n2 = size(w, 2)

        this%header = ''
        this%header(1) = '# Start 4000-byte ASCII header '
        this%header(2) = 'n1 = '//num2str(this%n1)
        this%header(3) = 'd1 = '//num2str(this%d1, '(es)')
        this%header(4) = 'o1 = '//num2str(this%o1, '(es)')
        this%header(5) = 'l1 = '//num2str(this%o1 + (this%n1 - 1)*this%d1, '(es)')
        this%header(6) = 'label1 = '//tidy(this%label1)
        this%header(7) = 'unit1 = '//tidy(this%unit1)
        this%header(8) = 'n2 = '//num2str(this%n2)
        this%header(9) = 'd2 = '//num2str(this%d2, '(es)')
        this%header(10) = 'o2 = '//num2str(this%o2, '(es)')
        this%header(11) = 'l2 = '//num2str(this%o2 + (this%n2 - 1)*this%d2, '(es)')
        this%header(12) = 'label2 = '//tidy(this%label2)
        this%header(13) = 'unit2 = '//tidy(this%unit2)
        this%header(size(this%header)) = '# End 4000-byte ASCII header'
        do i = 1, size(this%header)
            h = this%header(i)
            h(79:79) = achar(13)
            h(80:80) = achar(10)
            this%header(i) = h
        end do

        this%initialized = .true.

        this%array = w

    end subroutine init_grid_from_2d

    subroutine read_grid_2d(this, filename)

        class(grid2), intent(inout) :: this
        character(len=*), intent(in) :: filename

        integer :: funit, i

        funit = nint(rand()*1000)
        open (funit, file=tidy(filename), access='stream', form='unformatted', status='old')
        ! Read header part
        do i = 1, size(this%header)
            read(funit) this%header(i)
        end do

        call parsepar_int(this%header, 'n1', this%n1, int((get_file_size(filename) - 4000.0)/4))
        call parsepar_float(this%header, 'd1', this%d1, 1.0)
        call parsepar_float(this%header, 'o1', this%o1, 0.0)
        this%l1 = (this%n1 - 1)*this%d1 + this%o1
        call parsepar_string(this%header, 'label1', this%label1, 'Axis 1')
        call parsepar_string(this%header, 'unit1', this%unit1, 'unitless')
        call parsepar_int(this%header, 'n2', this%n2, int((get_file_size(filename) - 4000.0)/this%n1/4))
        call parsepar_float(this%header, 'd2', this%d2, 1.0)
        call parsepar_float(this%header, 'o2', this%o2, 0.0)
        this%l2 = (this%n2 - 1)*this%d2 + this%o2
        call parsepar_string(this%header, 'label2', this%label2, 'Axis 2')
        call parsepar_string(this%header, 'unit2', this%unit2, 'unitless')

        this%array = zeros(this%n1, this%n2)

        ! Read array part
        read (funit, pos=4001) this%array
        close (funit)

    end subroutine read_grid_2d

    subroutine save_grid_2d(this, filename)

        class(grid2), intent(inout) :: this
        character(len=*), intent(in) :: filename

        integer :: funit, i

        funit = nint(rand()*1000)
        open (funit, file=tidy(filename), access='stream', form='unformatted', status='replace')
        ! Write header part as ASCII
        do i = 1, size(this%header)
            write (funit) this%header(i)
        end do
        ! Write array part as raw binary
        write (funit, pos=4001) this%array
        close (funit)

    end subroutine save_grid_2d

    subroutine stdin_grid_2d(this)

        class(grid2), intent(inout) :: this

        integer :: i

        close (input_unit)
        open (input_unit, access='stream', form='unformatted')

        ! Read header part
        do i = 1, size(this%header)
            read(input_unit) this%header(i)
        end do

        ! Parse header
        call parsepar_int(this%header, 'n1', this%n1, int((get_stdin_size() - 4000.0)/4))
        call parsepar_float(this%header, 'd1', this%d1, 1.0)
        call parsepar_float(this%header, 'o1', this%o1, 0.0)
        this%l1 = (this%n1 - 1)*this%d1 + this%o1
        call parsepar_string(this%header, 'label1', this%label1, 'Axis 1')
        call parsepar_string(this%header, 'unit1', this%unit1, 'unitless')
        call parsepar_int(this%header, 'n2', this%n2, int((get_stdin_size() - 4000.0)/this%n1/4))
        call parsepar_float(this%header, 'd2', this%d2, 1.0)
        call parsepar_float(this%header, 'o2', this%o2, 0.0)
        this%l2 = (this%n2 - 1)*this%d2 + this%o2
        call parsepar_string(this%header, 'label2', this%label2, 'Axis 2')
        call parsepar_string(this%header, 'unit2', this%unit2, 'unitless')
        call alloc_array(this%array, [1, this%n1, 1, this%n2])

        ! Read array part
        read (input_unit, pos=4001) this%array
        close (input_unit)

    end subroutine stdin_grid_2d

    subroutine stdout_grid_2d(this)

        class(grid2), intent(inout) :: this

        integer :: i

        close (output_unit)
        open (output_unit, access='stream', form='unformatted', status='replace')
        ! Write header part as ASCII
        do i = 1, size(this%header)
            write (output_unit) this%header(i)
        end do
        ! Write array part as raw binary
        write (output_unit, pos=4001) this%array
        close (output_unit)

    end subroutine stdout_grid_2d

    function extract_grid_2d(this) result(w)

        class(grid2), intent(in) :: this

        real, allocatable, dimension(:, :) :: w

        w = this%array

    end function extract_grid_2d

    subroutine free_grid_2d(this)

        class(grid2), intent(inout) :: this

        if (allocated(this%array)) then
            deallocate(this%array)
        end if

    end subroutine free_grid_2d

    !===========================================================================

    subroutine init_grid_3d(this, n, d, o, label, unit, const)

        class(grid3), intent(inout) :: this
        integer, dimension(:), intent(in) :: n
        real, dimension(:), intent(in), optional :: d, o
        real, intent(in), optional :: const
        character(len=*), dimension(:), intent(in), optional :: label, unit

        character(len=80) :: h
        integer :: i

        call assert(size(n) == 3, 'size(n) /= 3')
        this%n1 = n(1)
        this%n2 = n(2)
        this%n3 = n(3)

        if (present(d)) then
            call assert(size(d) == 3, 'size(d) /= 3')
            this%d1 = d(1)
            this%d2 = d(2)
            this%d3 = d(3)
        end if
        if (present(o)) then
            call assert(size(o) == 3, 'size(d) /= 3')
            this%o1 = o(1)
            this%o2 = o(2)
            this%o3 = o(3)
        end if
        if (present(label)) then
            call assert(size(label) == 3, 'size(label) /= 3')
            this%label1 = label(1)
            this%label2 = label(2)
            this%label3 = label(3)
        end if
        if (present(unit)) then
            call assert(size(unit) == 3, 'size(unit) /= 3')
            this%unit1 = unit(1)
            this%unit2 = unit(2)
            this%unit3 = unit(3)
        end if
        this%l1 = (this%n1 - 1)*this%d1 + this%o1
        this%l2 = (this%n2 - 1)*this%d2 + this%o2
        this%l3 = (this%n3 - 1)*this%d3 + this%o3

        call alloc_array(this%array, [1, this%n1, 1, this%n2, 1, this%n3])
        if (present(const)) then
            this%array = const
        end if

        this%header = ''
        this%header(1) = '# Start 4000-byte ASCII header '
        this%header(2) = 'n1 = '//num2str(this%n1)
        this%header(3) = 'd1 = '//num2str(this%d1, '(es)')
        this%header(4) = 'o1 = '//num2str(this%o1, '(es)')
        this%header(5) = 'l1 = '//num2str(this%o1 + (this%n1 - 1)*this%d1, '(es)')
        this%header(6) = 'label1 = '//tidy(this%label1)
        this%header(7) = 'unit1 = '//tidy(this%unit1)
        this%header(8) = 'n2 = '//num2str(this%n2)
        this%header(9) = 'd2 = '//num2str(this%d2, '(es)')
        this%header(10) = 'o2 = '//num2str(this%o2, '(es)')
        this%header(11) = 'l2 = '//num2str(this%o2 + (this%n2 - 1)*this%d2, '(es)')
        this%header(12) = 'label2 = '//tidy(this%label2)
        this%header(13) = 'unit2 = '//tidy(this%unit2)
        this%header(14) = 'n3 = '//num2str(this%n3)
        this%header(15) = 'd3 = '//num2str(this%d3, '(es)')
        this%header(16) = 'o3 = '//num2str(this%o3, '(es)')
        this%header(17) = 'l3 = '//num2str(this%o3 + (this%n3 - 1)*this%d3, '(es)')
        this%header(18) = 'label3 = '//tidy(this%label3)
        this%header(19) = 'unit3 = '//tidy(this%unit3)
        this%header(size(this%header)) = '# End 4000-byte ASCII header'
        do i = 1, size(this%header)
            h = this%header(i)
            h(79:79) = achar(13)
            h(80:80) = achar(10)
            this%header(i) = h
        end do

        this%initialized = .true.

    end subroutine init_grid_3d

    subroutine init_grid_from_3d(this, w)

        class(grid3), intent(inout) :: this
        real, dimension(:, :, :), intent(in) :: w

        integer :: i
        character(len=80) :: h

        this%n1 = size(w, 1)
        this%n2 = size(w, 2)
        this%n3 = size(w, 3)

        this%header = ''
        this%header(1) = '# Start 4000-byte ASCII header '
        this%header(2) = 'n1 = '//num2str(this%n1)
        this%header(3) = 'd1 = '//num2str(this%d1, '(es)')
        this%header(4) = 'o1 = '//num2str(this%o1, '(es)')
        this%header(5) = 'l1 = '//num2str(this%o1 + (this%n1 - 1)*this%d1, '(es)')
        this%header(6) = 'label1 = '//tidy(this%label1)
        this%header(7) = 'unit1 = '//tidy(this%unit1)
        this%header(8) = 'n2 = '//num2str(this%n2)
        this%header(9) = 'd2 = '//num2str(this%d2, '(es)')
        this%header(10) = 'o2 = '//num2str(this%o2, '(es)')
        this%header(11) = 'l2 = '//num2str(this%o2 + (this%n2 - 1)*this%d2, '(es)')
        this%header(12) = 'label2 = '//tidy(this%label2)
        this%header(13) = 'unit2 = '//tidy(this%unit2)
        this%header(14) = 'n3 = '//num2str(this%n3)
        this%header(15) = 'd3 = '//num2str(this%d3, '(es)')
        this%header(16) = 'o3 = '//num2str(this%o3, '(es)')
        this%header(17) = 'l3 = '//num2str(this%o3 + (this%n3 - 1)*this%d3, '(es)')
        this%header(18) = 'label3 = '//tidy(this%label3)
        this%header(19) = 'unit3 = '//tidy(this%unit3)
        this%header(size(this%header)) = '# End 4000-byte ASCII header'
        do i = 1, size(this%header)
            h = this%header(i)
            h(79:79) = achar(13)
            h(80:80) = achar(10)
            this%header(i) = h
        end do

        this%initialized = .true.

        this%array = w

    end subroutine init_grid_from_3d

    subroutine read_grid_3d(this, filename)

        class(grid3), intent(inout) :: this
        character(len=*), intent(in) :: filename

        integer :: funit, i

        funit = nint(rand()*1000)
        open (funit, file=tidy(filename), access='stream', form='unformatted', status='old')
        ! Read header part
        do i = 1, size(this%header)
            read(funit) this%header(i)
        end do

        call parsepar_int(this%header, 'n1', this%n1, int((get_file_size(filename) - 4000.0)/4))
        call parsepar_float(this%header, 'd1', this%d1, 1.0)
        call parsepar_float(this%header, 'o1', this%o1, 0.0)
        this%l1 = (this%n1 - 1)*this%d1 + this%o1
        call parsepar_string(this%header, 'label1', this%label1, 'Axis 1')
        call parsepar_string(this%header, 'unit1', this%unit1, 'unitless')
        call parsepar_int(this%header, 'n2', this%n2, int((get_file_size(filename) - 4000.0)/this%n1/4))
        call parsepar_float(this%header, 'd2', this%d2, 1.0)
        call parsepar_float(this%header, 'o2', this%o2, 0.0)
        this%l2 = (this%n2 - 1)*this%d2 + this%o2
        call parsepar_string(this%header, 'label2', this%label2, 'Axis 2')
        call parsepar_string(this%header, 'unit2', this%unit2, 'unitless')
        call parsepar_int(this%header, 'n3', this%n3, int((get_file_size(filename) - 4000.0)/this%n1/this%n2/4))
        call parsepar_float(this%header, 'd3', this%d3, 1.0)
        call parsepar_float(this%header, 'o3', this%o3, 0.0)
        this%l3 = (this%n3 - 1)*this%d3 + this%o3
        call parsepar_string(this%header, 'label3', this%label3, 'Axis 2')
        call parsepar_string(this%header, 'unit3', this%unit3, 'unitless')

        this%array = zeros(this%n1, this%n2, this%n3)

        ! Read array part
        read (funit, pos=4001) this%array
        close (funit)

    end subroutine read_grid_3d

    subroutine save_grid_3d(this, filename)

        class(grid3), intent(inout) :: this
        character(len=*), intent(in) :: filename

        integer :: funit, i

        funit = nint(rand()*1000)
        open (funit, file=tidy(filename), access='stream', form='unformatted', status='replace')
        ! Write header part as ASCII
        do i = 1, size(this%header)
            write (funit) this%header(i)
        end do
        ! Write array part as raw binary
        write (funit, pos=4001) this%array
        close (funit)

    end subroutine save_grid_3d

    subroutine stdin_grid_3d(this)

        class(grid3), intent(inout) :: this

        integer :: i

        close (input_unit)
        open (input_unit, access='stream', form='unformatted')

        ! Read header part
        do i = 1, size(this%header)
            read(input_unit) this%header(i)
        end do

        ! Parse header
        call parsepar_int(this%header, 'n1', this%n1, int((get_stdin_size() - 4000.0)/4))
        call parsepar_float(this%header, 'd1', this%d1, 1.0)
        call parsepar_float(this%header, 'o1', this%o1, 0.0)
        this%l1 = (this%n1 - 1)*this%d1 + this%o1
        call parsepar_string(this%header, 'label1', this%label1, 'Axis 1')
        call parsepar_string(this%header, 'unit1', this%unit1, 'unitless')
        call parsepar_int(this%header, 'n2', this%n2, int((get_stdin_size() - 4000.0)/this%n1/4))
        call parsepar_float(this%header, 'd2', this%d2, 1.0)
        call parsepar_float(this%header, 'o2', this%o2, 0.0)
        this%l2 = (this%n2 - 1)*this%d2 + this%o2
        call parsepar_string(this%header, 'label2', this%label2, 'Axis 2')
        call parsepar_string(this%header, 'unit2', this%unit2, 'unitless')
        call parsepar_int(this%header, 'n3', this%n3, int((get_stdin_size() - 4000.0)/this%n1/this%n2/4))
        call parsepar_float(this%header, 'd3', this%d3, 1.0)
        call parsepar_float(this%header, 'o3', this%o3, 0.0)
        this%l3 = (this%n3 - 1)*this%d3 + this%o3
        call parsepar_string(this%header, 'label3', this%label3, 'Axis 2')
        call parsepar_string(this%header, 'unit3', this%unit3, 'unitless')
        call alloc_array(this%array, [1, this%n1, 1, this%n2, 1, this%n3])

        ! Read array part
        read (input_unit, pos=4001) this%array
        close (input_unit)

    end subroutine stdin_grid_3d

    subroutine stdout_grid_3d(this)

        class(grid3), intent(inout) :: this

        integer :: i

        close (output_unit)
        open (output_unit, access='stream', form='unformatted', status='replace')
        ! Write header part as ASCII
        do i = 1, size(this%header)
            write (output_unit) this%header(i)
        end do
        ! Write array part as raw binary
        write (output_unit, pos=4001) this%array
        close (output_unit)

    end subroutine stdout_grid_3d

    function extract_grid_3d(this) result(w)

        class(grid3), intent(in) :: this

        real, allocatable, dimension(:, :, :) :: w

        w = this%array

    end function extract_grid_3d

    subroutine free_grid_3d(this)

        class(grid3), intent(inout) :: this

        if (allocated(this%array)) then
            deallocate(this%array)
        end if

    end subroutine free_grid_3d

end module
