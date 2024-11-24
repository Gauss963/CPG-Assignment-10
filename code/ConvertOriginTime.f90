program ConvertOriginTime
    use iso_c_binding
    implicit none

    real(c_double) :: origin_time_seconds
    integer :: year, month, day, hour, minute, second

    origin_time_seconds = 757382400.0

    call unix_time_to_date(origin_time_seconds, year, month, day, hour, minute, second)

    print *, 'Origin Date and Time:', year, '-', month, '-', day, ' ', hour, ':', minute, ':', second

contains

    subroutine unix_time_to_date(seconds_since_epoch, year, month, day, hour, minute, second)
        use iso_fortran_env, only: int64
        implicit none
        real(c_double), intent(in) :: seconds_since_epoch
        integer, intent(out) :: year, month, day, hour, minute, second
        integer(int64) :: total_seconds, days, leap_days
        integer :: y, m, d
        integer, parameter :: seconds_in_day = 86400

        total_seconds = int(seconds_since_epoch, int64)

        ! Calculate days and remaining seconds
        days = total_seconds / seconds_in_day
        second = mod(total_seconds, seconds_in_day)

        ! Now convert days to date
        call days_to_date(days, year, month, day)

        ! Convert remaining seconds to hours, minutes, seconds
        hour = second / 3600
        minute = mod(second, 3600) / 60
        second = mod(second, 60)
    end subroutine unix_time_to_date

    subroutine days_to_date(days_since_epoch, year, month, day)
        implicit none
        integer(int64), intent(in) :: days_since_epoch
        integer, intent(out) :: year, month, day
        integer :: n, i
        integer, parameter :: month_days(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
        integer :: y
        integer :: days_in_year

        y = 1970
        n = days_since_epoch

        do
            if (is_leap_year(y)) then
                days_in_year = 366
            else
                days_in_year = 365
            end if
            if (n >= days_in_year) then
                n = n - days_in_year
                y = y + 1
            else
                exit
            end if
        end do

        year = y

        if (is_leap_year(year)) then
            month_days(2) = 29
        else
            month_days(2) = 28
        end if

        do i = 1, 12
            if (n >= month_days(i)) then
                n = n - month_days(i)
            else
                month = i
                day = n + 1
                exit
            end if
        end do
    end subroutine days_to_date

    logical function is_leap_year(year)
        implicit none
        integer, intent(in) :: year
        if (mod(year, 400) == 0) then
            is_leap_year = .true.
        else if (mod(year, 100) == 0) then
            is_leap_year = .false.
        else if (mod(year, 4) == 0) then
            is_leap_year = .true.
        else
            is_leap_year = .false.
        end if
    end function is_leap_year

end program ConvertOriginTime
