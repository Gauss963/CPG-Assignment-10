program PlotIFFTBP
    use iso_c_binding
    implicit none

    type, bind(c) :: Header
        character(kind=c_char), dimension(4) :: Code
        real(c_double)                       :: orgintime
        integer(c_int16_t)                   :: ncom
        integer(c_int32_t)                   :: ndata
        real(c_float)                        :: dt
    end type Header

    ! Declare variables
    type(Header) :: wh
    real(c_float), allocatable :: wd(:,:)
    integer :: i, j, k
    real(c_float), allocatable :: T(:)
    character(len=4) :: CodeString

    ! Variables for FFT
    integer :: n_fft, n_freq
    complex, allocatable :: fft_data(:)
    real, allocatable :: freq(:)
    real, allocatable :: amplitude(:,:)
    real, parameter :: pi = 3.14159265359

    ! Variables for PGPLOT
    integer :: pgid
    real :: min_T, max_T, min_E, max_E
    real :: min_F, max_F, min_A, max_A

    ! Variables for filtered data
    real(c_float), allocatable :: filtered_data(:)
    real, allocatable :: amplitude_filtered(:,:)

    ! Filter parameters
    real :: lowcorner = 0.5
    real :: highcorner = 3.0
    character(len=6) :: filter_type = 'BUTTER'
    integer :: order = 4
    character(len=3) :: pass_type = 'BP '  ! Band-pass
    integer :: passes = 2
    integer :: np

    ! Open and read binary file
    open(8, file='../data/seismicdata.bin', status='old', &
        form='unformatted', access='stream', convert='little_endian')
    read(8) wh

    ! Convert Code array to string
    CodeString = ''
    do i = 1, 4
        CodeString(i:i) = wh%Code(i)
    end do

    ! Print header info
    print *, 'Code:', trim(CodeString)
    print *, 'Origin Time:', wh%orgintime
    print *, 'Number of Components:', wh%ncom
    print *, 'Number of Data Points:', wh%ndata
    print *, 'Sampling Interval:', wh%dt

    ! Allocate and read waveform data
    allocate(wd(wh%ncom, wh%ndata))
    read(8) wd
    close(8)

    ! Create time array
    allocate(T(wh%ndata))
    do i = 1, wh%ndata
        T(i) = (i - 1) * wh%dt
    end do

    np = wh%ndata  ! Number of points

    ! Remove mean from each component
    do i = 1, wh%ncom
        wd(i,:) = wd(i,:) - (sum(wd(i,:)) / real(wh%ndata))
    end do

    ! Prepare for FFT
    ! Find next power of 2 for FFT
    n_fft = 2 ** ceiling(log(real(wh%ndata)) / log(2.0))
    n_freq = n_fft / 2 + 1

    allocate(fft_data(n_fft))
    allocate(freq(n_freq))
    allocate(amplitude(wh%ncom, n_freq))

    ! Calculate frequency array
    do i = 1, n_freq
        freq(i) = (i-1) * (1.0/(wh%dt * n_fft))
    end do

    ! Initialize PGPLOT
    call pgopen('../plot/waveform_fft.ps/VCPS')
    call pgsubp(2, 3)  ! 2 rows, 3 columns

    ! Process and plot each component
    do i = 1, wh%ncom
        ! Copy data to complex array and pad with zeros
        fft_data = cmplx(0.0, 0.0)
        do j = 1, wh%ndata
            fft_data(j) = cmplx(wd(i,j), 0.0)
        end do

        ! Perform FFT using built-in Fortran function
        fft_data = fft_data / sqrt(real(n_fft))  ! Normalize
        call fft(fft_data)

        ! Calculate amplitude spectrum
        do j = 1, n_freq
            amplitude(i,j) = sqrt(real(fft_data(j))**2 + aimag(fft_data(j))**2)
            ! Convert to dB
            if (amplitude(i,j) > 0.0) then
                amplitude(i,j) = 20.0 * log10(amplitude(i,j))
            else
                amplitude(i,j) = -200.0  ! Set minimum dB value
            end if
        end do

        ! Plot time domain
        call pgsci(1)
        min_E = minval(wd(i,:))
        max_E = maxval(wd(i,:))
        call pgenv(minval(T), maxval(T), min_E, max_E, 0, 0)
        call pglab('Time (sec)', 'Amplitude', 'Waveform E'//char(i+48))
        call pgsci(i+1)
        call pgline(wh%ndata, T, wd(i,:))

        ! Plot frequency domain
        call pgsci(1)
        min_A = minval(amplitude(i,:))
        max_A = maxval(amplitude(i,:))
        call pgenv(0.0, min(50.0, maxval(freq)), min_A, max_A, 0, 0)  ! Limit to 50 Hz
        call pglab('Frequency (Hz)', 'Amplitude (dB)', 'Spectrum E'//char(i+48))
        call pgsci(i+1)
        call pgline(n_freq, freq, amplitude(i,:))
    end do

    call pgclos()

    ! Apply IIR filter to the data in time domain
    ! Allocate filtered data array
    allocate(filtered_data(wh%ndata))
    allocate(amplitude_filtered(wh%ncom, n_freq))

    ! Open new PGPLOT device
    call pgopen('../plot/waveform_fft_BP.ps/VCPS')
    call pgsubp(2, 3)  ! 2 rows, 3 columns

    ! Process and plot each component
    do i = 1, wh%ncom
        ! Copy the original data to filtered_data array
        filtered_data = wd(i,:)

        ! Apply IIR filter
        call IIRFILT(filtered_data, np, filter_type, order, pass_type, &
            lowcorner, highcorner, wh%dt, passes)

        ! Copy filtered data to complex array for FFT
        fft_data = cmplx(0.0, 0.0)
        do j = 1, wh%ndata
            fft_data(j) = cmplx(filtered_data(j), 0.0)
        end do

        ! Perform FFT on filtered data
        fft_data = fft_data / sqrt(real(n_fft))  ! Normalize
        call fft(fft_data)

        ! Calculate amplitude spectrum of filtered data
        do j = 1, n_freq
            amplitude_filtered(i,j) = sqrt(real(fft_data(j))**2 + aimag(fft_data(j))**2)
            ! Convert to dB
            if (amplitude_filtered(i,j) > 0.0) then
                amplitude_filtered(i,j) = 20.0 * log10(amplitude_filtered(i,j))
            else
                amplitude_filtered(i,j) = -200.0  ! Set minimum dB value
            end if
        end do

        ! Plot filtered time domain
        call pgsci(1)
        min_E = minval(filtered_data)
        max_E = maxval(filtered_data)
        call pgenv(minval(T), maxval(T), min_E, max_E, 0, 0)
        call pglab('Time (sec)', 'Amplitude', 'Filtered Waveform E'//char(i+48))
        call pgsci(i+1)
        call pgline(wh%ndata, T, filtered_data)

        ! Plot filtered frequency domain
        call pgsci(1)
        min_A = minval(amplitude_filtered(i,:))
        max_A = maxval(amplitude_filtered(i,:))
        call pgenv(0.0, min(50.0, maxval(freq)), min_A, max_A, 0, 0)  ! Limit to 50 Hz
        call pglab('Frequency (Hz)', 'Amplitude (dB)', 'Filtered Spectrum E'//char(i+48))
        call pgsci(i+1)
        call pgline(n_freq, freq, amplitude_filtered(i,:))
    end do

    call pgclos()

    ! Cleanup
    deallocate(wd, T, fft_data, freq, amplitude, filtered_data, amplitude_filtered)

contains
    ! FFT subroutine using recursive algorithm
    recursive subroutine fft(x)
        complex, intent(inout) :: x(:)
        complex, allocatable :: even(:), odd(:)
        complex :: w
        integer :: n, k

        n = size(x)
        if (n <= 1) return

        allocate(even(n/2), odd(n/2))

        ! Split into even and odd
        even = x(1:n:2)
        odd = x(2:n:2)

        ! Recursive FFT
        call fft(even)
        call fft(odd)

        ! Combine
        do k = 1, n/2
            w = cmplx(cos(-2.0*pi*(k-1)/n), sin(-2.0*pi*(k-1)/n))
            x(k) = even(k) + w * odd(k)
            x(k + n/2) = even(k) - w * odd(k)
        end do

        deallocate(even, odd)
    end subroutine fft

end program PlotIFFTBP