program epicentral
    implicit none
    character(len = 100) :: line

    integer :: n, m, i
    integer :: io_status, unit_num

    real, allocatable :: CONTOUR_X(:), CONTOUR_Y(:)
    real :: x, y, z, w
    real, allocatable :: EQ_X(:), EQ_Y(:), EQ_Z(:), EQ_M(:)
    real, allocatable :: EQ_X_proj(:), EQ_Y_proj(:)
    real, allocatable :: EQ_Z_inrange(:), EQ_M_inrange(:)
    real, allocatable :: EQ_X_proj_lon(:), EQ_Y_proj_lat(:)
    real, allocatable :: EQ_T(:), EQ_Dist_Along_Profile(:)

    real :: ax_x_max, ax_x_min, ax_y_max, ax_y_min, margen

    real :: EQ_M_min, EQ_M_max
    real :: EQ_Z_min, EQ_Z_max
    real :: symbol_size, size_min, size_max

    real :: R_C0, G_C0, B_C0
    real :: R_C1, G_C1, B_C1
    real :: symbol_color_R, symbol_color_G, symbol_color_B

    real :: pi, cos_lat_ref, lon_ref, lat_ref
    real :: lon1, lat1, lon2, lat2
    real :: x1, y1, x2, y2, x0, y0
    real :: dx_AB, dy_AB, dx_AP, dy_AP
    real :: distance, t, numerator, denominator
    real :: x_proj, y_proj
    real :: length_profile
    integer :: m_inrange

    real :: profile_lon(2), profile_lat(2)

    ! Get contour data -------------------------------------------------------------------
    open(newunit = unit_num, file = "../data/Taiwan.txt", status = "old", action = "read")
    n = 0
    do
        read(unit_num, '(A)', iostat=io_status) line
        if (io_status /= 0) exit
        n = n + 1
    end do
    close(unit_num)
    allocate(CONTOUR_X(n), CONTOUR_Y(n))
    open(newunit=unit_num, file = "../data/Taiwan.txt", status = "old", action = "read")
    do i = 1, n
        read(unit_num, *, iostat=io_status) CONTOUR_X(i), CONTOUR_Y(i)
        if (io_status /= 0) exit
    end do
    close(unit_num)

    ! Initialize constants ----------------------------------------------------------------
    pi = 4.0 * atan(1.0)
    lon_ref = 120.0
    lat_ref = 23.0
    cos_lat_ref = cos(lat_ref * pi / 180.0)

    lon1 = 120.0
    lat1 = 23.0
    lon2 = 122.0
    lat2 = 25.0

    ! Compute x1, y1, x2, y2 in km
    x1 = (lon1 - lon_ref) * 111.0 * cos_lat_ref
    y1 = (lat1 - lat_ref) * 111.0
    x2 = (lon2 - lon_ref) * 111.0 * cos_lat_ref
    y2 = (lat2 - lat_ref) * 111.0

    dx_AB = x2 - x1
    dy_AB = y2 - y1

    length_profile = sqrt(dx_AB**2 + dy_AB**2)

    ! Prepare profile line for plotting
    profile_lon(1) = lon1
    profile_lat(1) = lat1
    profile_lon(2) = lon2
    profile_lat(2) = lat2

    ! Get earthquake data ----------------------------------------------------------------
    open(newunit = unit_num, file = "../data/1999.lis", status = "old", action = "read")
    m = 0
    do
        read(unit_num, '(A)', iostat=io_status) line
        if (io_status /= 0) exit
        m = m + 1
    end do
    close(unit_num)
    allocate(EQ_X(m), EQ_Y(m), EQ_Z(m), EQ_M(m))
    allocate(EQ_X_proj(m), EQ_Y_proj(m))
    allocate(EQ_Z_inrange(m), EQ_M_inrange(m))
    allocate(EQ_X_proj_lon(m), EQ_Y_proj_lat(m))
    allocate(EQ_T(m), EQ_Dist_Along_Profile(m))

    m_inrange = 0

    open(newunit=unit_num, file = "../data/1999.lis", status = "old", action = "read")
    do i = 1, m
        read(unit_num, '(18X, F2.0, F5.2, F3.0, F5.2, F6.2, F5.2)', iostat=io_status) &
            x, y, z, w, EQ_Z(i), EQ_M(i)
        EQ_X(i) = (z + w / 60.0)
        EQ_Y(i) = (x + y / 60.0)
        if (io_status /= 0) exit

        ! Compute x0, y0 in km
        x0 = (EQ_X(i) - lon_ref) * 111.0 * cos_lat_ref
        y0 = (EQ_Y(i) - lat_ref) * 111.0

        dx_AP = x0 - x1
        dy_AP = y0 - y1

        numerator = abs(dx_AB * dy_AP - dx_AP * dy_AB)
        denominator = sqrt(dx_AB**2 + dy_AB**2)
        distance = numerator / denominator

        t = (dx_AB * dx_AP + dy_AB * dy_AP) / (dx_AB**2 + dy_AB**2)

        if (distance <= 50.0 .and. t >= 0.0 .and. t <= 1.0) then
            m_inrange = m_inrange + 1

            x_proj = x1 + t * dx_AB
            y_proj = y1 + t * dy_AB

            EQ_X_proj(m_inrange) = x_proj
            EQ_Y_proj(m_inrange) = y_proj
            EQ_Z_inrange(m_inrange) = EQ_Z(i)
            EQ_M_inrange(m_inrange) = EQ_M(i)

            EQ_X_proj_lon(m_inrange) = x_proj / (111.0 * cos_lat_ref) + lon_ref
            EQ_Y_proj_lat(m_inrange) = y_proj / 111.0 + lat_ref

            EQ_T(m_inrange) = t
            EQ_Dist_Along_Profile(m_inrange) = t * length_profile
        end if
    end do
    close(unit_num)

    print *, 'Number of events within 50 km of profile:', m_inrange

    if (m_inrange == 0) then
        print *, 'No events within 50 km of the profile.'
        stop
    end if

    ! Set x, y lims ----------------------------------------------------------------------
    margen = 0.25  ! 0.25 degree larger
    ax_x_min = minval(EQ_X_proj_lon(1:m_inrange)) - margen
    ax_x_max = maxval(EQ_X_proj_lon(1:m_inrange)) + margen
    ax_y_min = minval(EQ_Y_proj_lat(1:m_inrange)) - margen
    ax_y_max = maxval(EQ_Y_proj_lat(1:m_inrange)) + margen

    ! Set scatter size lims --------------------------------------------------------------
    EQ_M_min = minval(EQ_M_inrange(1:m_inrange))
    EQ_M_max = maxval(EQ_M_inrange(1:m_inrange))
    size_min = 0.5
    size_max = 3.0

    ! Set scatter color
    R_C0 = 0.121
    G_C0 = 0.466
    B_C0 = 0.706

    R_C1 = 1.00
    G_C1 = 0.498
    B_C1 = 0.055

    EQ_Z_min = minval(EQ_Z_inrange(1:m_inrange))
    EQ_Z_max = maxval(EQ_Z_inrange(1:m_inrange))

    ! Plot Map View -----------------------------------------------------------------------
    call pgopen('../plot/1999_event_distribution.ps/VCPS')
    call pgsci(1)
    call pgenv(ax_x_min, ax_x_max, ax_y_min, ax_y_max, 0, 1)
    call pgscf(1)
    call pglab('Longitude (E)', 'Latitude (N)', '1999 Projected Events within 50 km of Profile')

    ! Plot profile line
    call pgsci(2)
    call pgline(2, profile_lon, profile_lat)

    ! Plot symbols -----------------------------------------------------------------------
    do i = 1, m_inrange
        symbol_size = size_min + (EQ_M_inrange(i) - EQ_M_min) * (size_max - size_min) / (EQ_M_max - EQ_M_min)

        symbol_color_R = R_C0 + (EQ_Z_inrange(i) - EQ_Z_min) * (R_C1 - R_C0) / (EQ_Z_max - EQ_Z_min)
        symbol_color_G = G_C0 + (EQ_Z_inrange(i) - EQ_Z_min) * (G_C1 - G_C0) / (EQ_Z_max - EQ_Z_min)
        symbol_color_B = B_C0 + (EQ_Z_inrange(i) - EQ_Z_min) * (B_C1 - B_C0) / (EQ_Z_max - EQ_Z_min)

        call pgscr(42, symbol_color_R, symbol_color_G, symbol_color_B)
        call pgsci(42)
        call pgsch(symbol_size)
        call pgpt1(EQ_X_proj_lon(i), EQ_Y_proj_lat(i), 22)
    end do

    ! Plot contour
    call pgsci(1)
    call pgline(n, CONTOUR_X, CONTOUR_Y)

    call pgclos()

    ! Add cross-sectional plot -------------------------------------------------------------

    call pgopen('../plot/1999_event_cut.ps/VCPS')

    ! Set up axes for cross-sectional plot
    ax_x_min = 0.0
    ax_x_max = length_profile
    ax_y_min = EQ_Z_max + 10.0
    ax_y_max = EQ_Z_min - 10.0

    call pgsci(1)
    call pgenv(ax_x_min, ax_x_max, ax_y_min, ax_y_max, 0, -1)
    call pgscf(1)
    call pglab('Distance along profile (km)', 'Depth (km)', 'Cross-sectional View of Projected Events')

    ! Plot events in cross-sectional view
    do i = 1, m_inrange
        symbol_size = size_min + (EQ_M_inrange(i) - EQ_M_min) * (size_max - size_min) / (EQ_M_max - EQ_M_min)

        symbol_color_R = R_C0 + (EQ_Z_inrange(i) - EQ_Z_min) * (R_C1 - R_C0) / (EQ_Z_max - EQ_Z_min)
        symbol_color_G = G_C0 + (EQ_Z_inrange(i) - EQ_Z_min) * (G_C1 - G_C0) / (EQ_Z_max - EQ_Z_min)
        symbol_color_B = B_C0 + (EQ_Z_inrange(i) - EQ_Z_min) * (B_C1 - B_C0) / (EQ_Z_max - EQ_Z_min)

        call pgscr(42, symbol_color_R, symbol_color_G, symbol_color_B)
        call pgsci(42)
        call pgsch(symbol_size)
        call pgpt1(EQ_Dist_Along_Profile(i), EQ_Z_inrange(i), 22)
    end do

    ! Draw profile line
    call pgsci(2)
    call pgline(2, (/0.0, length_profile/), (/0.0, 0.0/))

    call pgclos()

    ! Deallocate arrays
    deallocate(CONTOUR_X, CONTOUR_Y)
    deallocate(EQ_X, EQ_Y, EQ_Z, EQ_M)
    deallocate(EQ_X_proj, EQ_Y_proj)
    deallocate(EQ_Z_inrange, EQ_M_inrange)
    deallocate(EQ_X_proj_lon, EQ_Y_proj_lat)
    deallocate(EQ_T, EQ_Dist_Along_Profile)
end program epicentral
