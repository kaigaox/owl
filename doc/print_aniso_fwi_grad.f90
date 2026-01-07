
program test

    use libflit

    character(len=8), allocatable, dimension(:) :: param
    character(len=32) :: label

    integer :: i, ci, cj

    !    ! Thomsen epsilon-delta parameterization
    !    label = 'thomsen'
    !    param = ['vp', 'vs', 'epsilon', 'delta', 'gamma']

    ! Alkhalifah-Tsvankin epsilon-eta parameterization
    label = 'alkhalifah_tsvankin'
    param = ['vp', 'vs', 'epsilon', 'eta', 'gamma']

    open (3, file='./grad_'//tidy(label)//'.txt')

    do i = 1, size(param)

        write (3, *) 'grad = &'
        do ci = 1, 6
            do cj = 1, 6

                !                ! For general anisotropy case
                !                if (cj >= ci) then

                ! For isotropic-VTI-HTI
                if ((ci <= 3 .and. cj >= ci .and. cj <= 3) .or. (ci >=4 .and. ci == cj)) then

                    write(3, *) ' + '//tidy(label)//'_dc'//num2str(ci)//num2str(cj)//'_d' &
                        //tidy(param(i)(1:3))//'*grad_c'//num2str(ci)//num2str(cj)//' &'

                end if

            end do
        end do

        write (3, *) 'grd%array = permute(grad, 321)'
        write (3, '(a)') "call grd%output(tidy(dir_working)// &
            '/shot_'//num2str(sgmtr%id)//'_grad_"//tidy(param(i))//".grd')"
        write (3, *) ''

    end do

end program

