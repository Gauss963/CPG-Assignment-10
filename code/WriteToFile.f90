subroutine write_to_file(filename, data)
    
    integer, parameter :: sp = kind(1.0)
    type :: SeismicData
        real(sp), allocatable :: waveform(:)
        real(sp), allocatable :: time(:)
        integer :: n_points
    end type SeismicData

    character(len=*), intent(in) :: filename
    type(SeismicData), intent(in) :: data
    integer :: i

    
    open(10, file=filename, status='replace')
    do i = 1, data%n_points
        write(10, *) data%waveform(i)
    end do
    close(10)
    print *, "Waveform data written to:", filename
end subroutine write_to_file