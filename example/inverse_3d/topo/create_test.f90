
!
! Create the models and the geometry for 3-D elastic FWI of Mount Rainier
!
! The model simplifies the geology of Mount Rainier. The volcanic edifice, a stack of
! andesite lava flows and breccias parallel to the flanks of the cone, overlies an
! unconformity at an elevation of about 2 km. Below the unconformity lie layered and
! faulted Tertiary volcanic rocks, intruded by the granodiorite Tatoosh pluton south of
! the summit. Hydrothermal alteration weakens the upper edifice beneath the summit, and
! weathering lowers the velocities near the surface. dem_to_topo.py creates the
! topography from a DEM, and the model depth z is measured downward from the summit.
!
program test

    use libflit
    use librgm

    implicit none

    ! Model grid
    integer, parameter :: nx = 201, ny = 201, nz = 101
    real, parameter :: d = 60.0
    ! Grid points of the edifice layers, which follow the smoothed surface
    integer, parameter :: nze = 45

    type(rgm3_curved) :: pb, pe
    real, allocatable, dimension(:, :) :: topo, zs, zsm, eu, thick
    real, allocatable, dimension(:, :, :) :: vp, vs, rho, ratio, m, sp, ss, vp0, vs0
    real :: hmax, xc, yc, x, y, z, e, r, rn, theta, dep, w
    integer :: i, j, k, l, ic(2), ishot, is1, is2, ir1, ir2

    real :: f0 = 1.5
    real :: sz = 25.0
    real :: rz = 0.0

    call make_directory('./model')

    ! Topography (elevation above sea level), the depth of the surface below the summit,
    ! and the location of the summit
    topo = load('./model/topo.bin', ny, nx)
    hmax = maxval(topo)
    zs = hmax - topo
    ic = maxloc(topo)
    yc = (ic(1) - 1)*d
    xc = (ic(2) - 1)*d

    ! Tertiary basement: layered and faulted volcanic rocks
    pb%n1 = nz
    pb%n2 = ny
    pb%n3 = nx
    pb%seed = 1357
    pb%nl = 12
    pb%refl_shape = 'perlin'
    pb%refl_shape_top = 'perlin'
    pb%refl_smooth = 3
    pb%refl_smooth_top = 0
    pb%lwv = 0.3
    pb%lwh = 0.3
    pb%refl_height = [0, 15]
    pb%refl_height_top = [0, 3]
    pb%nf = 4
    pb%disp = [3.0, 8.0]
    pb%yn_vary_disp = .true.
    pb%yn_elastic = .true.
    pb%vpvsratio = [1.70, 1.80]
    call pb%generate
    pb%vs = pb%vp/pb%vs
    pb%vp = rescale(pb%vp, [4200.0, 5600.0])

    ! Edifice: lava flows and breccias, whose velocities increase with depth
    pe%n1 = nze
    pe%n2 = ny
    pe%n3 = nx
    pe%seed = 2468
    pe%nl = 20
    pe%refl_shape = 'perlin'
    pe%refl_shape_top = 'perlin'
    pe%refl_smooth = 3
    pe%refl_smooth_top = 0
    pe%lwv = 0.5
    pe%lwh = 0.5
    pe%refl_height = [0, 4]
    pe%refl_height_top = [0, 2]
    pe%nf = 0
    pe%yn_fault = .false.
    pe%yn_elastic = .true.
    pe%vpvsratio = [1.72, 1.85]
    call pe%generate
    pe%vs = pe%vp/pe%vs
    pe%vp = rescale(pe%vp, [3000.0, 4200.0])

    ! The edifice layers follow the surface smoothed over 600 m. The edifice lies above an
    ! unconformity at an elevation of 1900 m plus 15% of the difference between the
    ! topography smoothed over 2 km and 1900 m; in the deep valleys, the basement crops out
    zsm = hmax - gauss_filt(topo, [10.0, 10.0])
    eu = 1900.0 + 0.15*(gauss_filt(topo, [33.0, 33.0]) - 1900.0)
    vp = zeros(nz, ny, nx)
    ratio = zeros(nz, ny, nx)
    !$omp parallel do private(i, j, k, l, z, e, w) collapse(2)
    do j = 1, nx
        do i = 1, ny
            do k = 1, nz
                e = hmax - (k - 1)*d
                if (e > eu(i, j)) then
                    z = min(max((k - 1) - zsm(i, j)/d + 1.0, 1.0), nze - 1.0e-3)
                    l = floor(z)
                    w = z - l
                    vp(k, i, j) = (1.0 - w)*pe%vp(l, i, j) + w*pe%vp(l + 1, i, j)
                    ratio(k, i, j) = (1.0 - w)*pe%vs(l, i, j) + w*pe%vs(l + 1, i, j)
                else
                    vp(k, i, j) = pb%vp(k, i, j)
                    ratio(k, i, j) = pb%vs(k, i, j)
                end if
            end do
        end do
    end do
    !$omp end parallel do

    vs = vp/ratio
    rho = 310.0*vp**0.25

    ! Weathered zone: the velocities increase from 80% (P) and 76% (S) of those of the
    ! fresh rock at the surface to 100% at the base of the zone; the zone is 60 m thick
    ! on the ridges and up to 180 m thick in the valleys, where the topography is concave
    thick = 60.0 + 120.0*clip((gauss_filt(topo, [2.5, 2.5]) - topo)/40.0, 0.0, 1.0)
    !$omp parallel do private(i, j, k, dep, w) collapse(2)
    do j = 1, nx
        do i = 1, ny
            do k = 1, nz
                dep = (k - 1)*d - zs(i, j)
                if (dep >= 0 .and. dep < thick(i, j)) then
                    w = 0.8 + 0.2*dep/thick(i, j)
                    vp(k, i, j) = vp(k, i, j)*w
                    vs(k, i, j) = vs(k, i, j)*w**1.2
                    rho(k, i, j) = 310.0*vp(k, i, j)**0.25
                end if
            end do
        end do
    end do
    !$omp end parallel do

    ! Above the surface, fill each column with the values just below the surface,
    ! and mask the gradients there
    m = ones(nz, ny, nx)
    do j = 1, nx
        do i = 1, ny
            l = ceiling(zs(i, j)/d) + 1
            vp(1:l - 1, i, j) = vp(l, i, j)
            vs(1:l - 1, i, j) = vs(l, i, j)
            rho(1:l - 1, i, j) = rho(l, i, j)
            m(1:l - 1, i, j) = 0.0
        end do
    end do

    ! Initial models: smooth the slownesses of the model without the pluton and the
    ! alteration, in coordinates that follow the surface, with Gaussian widths of 180 m
    ! in depth and 1.2 km horizontally
    sp = zeros(nz, ny, nx)
    ss = zeros(nz, ny, nx)
    !$omp parallel do private(i, j, k, l, z, w) collapse(2)
    do j = 1, nx
        do i = 1, ny
            do k = 1, nz
                z = min(k + zs(i, j)/d, nz - 1.0e-3)
                l = floor(z)
                w = z - l
                sp(k, i, j) = 1.0/((1.0 - w)*vp(l, i, j) + w*vp(l + 1, i, j))
                ss(k, i, j) = 1.0/((1.0 - w)*vs(l, i, j) + w*vs(l + 1, i, j))
            end do
        end do
    end do
    !$omp end parallel do
    sp = gauss_filt(sp, [3.0, 20.0, 20.0])
    ss = gauss_filt(ss, [3.0, 20.0, 20.0])
    vp0 = zeros(nz, ny, nx)
    vs0 = zeros(nz, ny, nx)
    !$omp parallel do private(i, j, k, l, z, w) collapse(2)
    do j = 1, nx
        do i = 1, ny
            do k = 1, nz
                z = min(max(k - zs(i, j)/d, 1.0), nz - 1.0e-3)
                l = floor(z)
                w = z - l
                vp0(k, i, j) = 1.0/((1.0 - w)*sp(l, i, j) + w*sp(l + 1, i, j))
                vs0(k, i, j) = 1.0/((1.0 - w)*ss(l, i, j) + w*ss(l + 1, i, j))
            end do
        end do
    end do
    !$omp end parallel do

    call output_array(vp0, './model/vp_init.bin')
    call output_array(vs0, './model/vs_init.bin')

    ! Tatoosh pluton: granodiorite centered 3 km south of the summit, with semi-axes of
    ! 4.5 km (east-west) and 3.5 km (north-south) and irregular flanks; its domed top
    ! lies at an elevation of 1300 m and descends below the model toward the flanks
    !$omp parallel do private(i, j, k, x, y, e, rn, theta) collapse(2)
    do j = 1, nx
        do i = 1, ny
            x = (j - 1)*d
            y = (i - 1)*d
            theta = atan2((y - yc + 3000.0)/3500.0, (x - xc)/4500.0)
            rn = sqrt(((x - xc)/4500.0)**2 + ((y - yc + 3000.0)/3500.0)**2) &
                *(1.0 + 0.1*cos(3*theta + 0.4) + 0.05*cos(5*theta + 1.9))
            do k = 1, nz
                e = hmax - (k - 1)*d
                if (e <= 1300.0 - 3200.0*rn**2) then
                    vp(k, i, j) = 5900.0
                    vs(k, i, j) = 3400.0
                    rho(k, i, j) = 2700.0
                end if
            end do
        end do
    end do
    !$omp end parallel do

    ! Hydrothermal alteration beneath the summit: within 900 m of the axis, from 150 m
    ! below the surface down to an elevation of 2800 m, the velocities decrease by up to
    ! 15% (P) and 25% (S), with smooth edges
    !$omp parallel do private(i, j, k, x, y, e, r, dep, w) collapse(2)
    do j = 1, nx
        do i = 1, ny
            x = (j - 1)*d
            y = (i - 1)*d
            r = sqrt((x - xc)**2 + (y - yc)**2)
            if (r < 900.0) then
                do k = 1, nz
                    e = hmax - (k - 1)*d
                    dep = (k - 1)*d - zs(i, j)
                    w = (1.0 - (r/900.0)**2)**2*clip((e - 2800.0)/300.0, 0.0, 1.0)*clip((dep - 150.0)/150.0, 0.0, 1.0)
                    if (w > 0) then
                        vp(k, i, j) = vp(k, i, j)*(1.0 - 0.15*w)
                        vs(k, i, j) = vs(k, i, j)*(1.0 - 0.25*w)
                        rho(k, i, j) = 310.0*vp(k, i, j)**0.25
                    end if
                end do
            end if
        end do
    end do
    !$omp end parallel do

    call output_array(vp, './model/vp.bin')
    call output_array(vs, './model/vs.bin')
    call output_array(rho, './model/rho.bin')
    call output_array(m, './model/mask.bin')

    print *, ' Summit at x, y = ', xc, yc, ', elevation ', hmax
    print *, ' vp range: ', minval(vp), maxval(vp)
    print *, ' vs range: ', minval(vs), maxval(vs)
    print *, ' rho range: ', minval(rho), maxval(rho)

    ! Geometry: 36 explosive sources 25 m below the surface on a 2 km grid, and
    ! 1521 three-component receivers on the surface on a 300 m grid
    call make_directory('./geometry')
    open (3, file='./geometry/geometry.txt')
    do ishot = 1, 36
        write (3, '(a)') 'shot_'//num2str(ishot)//'_geometry.txt'
    end do
    close (3)

    ishot = 1
    do is2 = 1, 6
        do is1 = 1, 6
            open (3, file='./geometry/shot_'//num2str(ishot)//'_geometry.txt')
            write (3, *) ishot
            write (3, *)
            write (3, *) 1
            write (3, *) 1000.0 + (is2 - 1)*2000.0, 1000.0 + (is1 - 1)*2000.0, sz
            write (3, *) 'explosion'
            write (3, *) 'ricker', f0, 1e9, 0.0
            write (3, *) 0, 0
            write (3, *)
            write (3, *) 39*39
            do ir2 = 1, 39
                do ir1 = 1, 39
                    write (3, '(4es)') ir2*300.0, ir1*300.0, rz, 1.0
                end do
            end do
            close (3)
            ishot = ishot + 1
        end do
    end do

end program
