 program epicentral
    implicit none
    character(len = 100) :: line
    
    integer :: n, m, i
    integer :: io_status, unit_num

    real, allocatable :: CONTOUR_X(:), CONTOUR_Y(:)
    real :: x, y, z, w
    real, allocatable :: EQ_X(:), EQ_Y(:), EQ_Z(:), EQ_M(:)
    
    real :: ax_x_max, ax_x_min, ax_y_max, ax_y_min, margen

    real :: EQ_M_min, EQ_M_max
    real :: EQ_Z_min, EQ_Z_max
    real :: symbol_size, size_min, size_max

    real :: R_C0, G_C0, B_C0
    real :: R_C1, G_C1, B_C1
    real :: symbol_color_R, symbol_color_G, symbol_color_B

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

    open(newunit=unit_num, file = "../data/1999.lis", status = "old", action = "read")
    do i = 1, m
        read(unit_num, '(18X, F2.0, F5.2, F3.0, F5.2, F6.2, F5.2)', iostat=io_status) x, y, z, w, EQ_Z(i), EQ_M(i)
        EQ_X(i) = (z + w / 60)
        EQ_Y(i) = (x + y / 60)


        if (io_status /= 0) exit
    end do
    close(unit_num)


    ! Set x, y lims ----------------------------------------------------------------------
    margen = 0.25  ! 0.25 degree larger
    ax_x_min = minval(EQ_X) - margen
    ax_x_max = maxval(EQ_X) + margen
    ax_y_min = minval(EQ_Y) - margen * 2
    ax_y_max = maxval(EQ_Y) + margen * 2

    ! Set scatter size lims --------------------------------------------------------------
    EQ_M_min = minval(EQ_M)
    EQ_M_max = maxval(EQ_M)
    size_min = 0.5
    size_max = 3.0

    ! Set scatter color
    R_C0 = 0.121
    G_C0 = 0.466
    B_C0 = 0.706

    R_C1 = 1.00
    G_C1 = 0.498
    B_C1 = 0.055

    EQ_Z_min = minval(EQ_Z)
    EQ_Z_max = maxval(EQ_Z)


    ! Plot Contour -----------------------------------------------------------------------
    call pgopen('../plot/1999_event_distribution_all.ps/VCPS')
    call pgsci(1)
    call pgenv(ax_x_min, ax_x_max, ax_y_min, ax_y_max, 0, 1)
    call pgscf(1)
    call pglab('Longitude (E)', 'Latitude (N)', '1999 Event Distribution')
    ! call pgline(n, CONTOUR_X, CONTOUR_Y)
    

    ! Plot symbol size with "EQ_Z", chahge to EQ_M later ---------------------------------
    ! call pgscr(42, 31/255, 119/255, 180/255)
    ! call pgscr(42, R_C0, G_C0, B_C0)
    ! call pgsci(13)
    ! call pgsci(42)

    do i = 1, m
        symbol_size = size_min + (EQ_M(i) - EQ_M_min) * (size_max - size_min) / (EQ_M_max - EQ_M_min)
        
        symbol_color_R = R_C0 + (EQ_Z(i) - EQ_Z_min) * (R_C1 - R_C0) / (EQ_Z_max - EQ_Z_min)
        symbol_color_G = G_C0 + (EQ_Z(i) - EQ_Z_min) * (G_C1 - G_C0) / (EQ_Z_max - EQ_Z_min)
        symbol_color_B = B_C0 + (EQ_Z(i) - EQ_Z_min) * (B_C1 - B_C0) / (EQ_Z_max - EQ_Z_min)


        call pgscr(42, symbol_color_R, symbol_color_G, symbol_color_B)
        call pgsci(42)
        call pgsch(symbol_size)
        

        call pgpt1(EQ_X(i), EQ_Y(i), 22)
    end do
    
    call pgsci(1)
    call pgline(n, CONTOUR_X, CONTOUR_Y)
    
    call pgclos()


    deallocate(CONTOUR_X, CONTOUR_Y)
    deallocate(EQ_X, EQ_Y, EQ_Z, EQ_M)
end program epicentral
