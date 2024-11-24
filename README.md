# Assignment 9

## Task: Read, Process, and Plot Seismic Waveform Data

This project involves reading seismic waveform data, processing it, and visualizing the results using the PGPLOT library. The tasks include waveform analysis, time conversion, Fast Fourier Transform (FFT) for frequency analysis, and generating an epicentral cross-section plot. Additionally, a Python script is provided for converting plots to PDF format.

## Files

### 1. `PlotWaveform.f90`
This Fortran program reads a binary seismic data file (`seismicdata.bin`), extracts the header information, and plots the waveform data using PGPLOT.

**Key Features:**
- Reads a binary file in little-endian format.
- Extracts header information, including:
  - Code: Identification code (4 characters).
  - Origin Time: The origin time of the seismic event.
  - Number of Components (`ncom`): Typically 3 (East, North, Vertical).
  - Number of Data Points (`ndata`).
  - Sampling Interval (`dt`).
- Uses PGPLOT for visualization and saves the waveform plot as a PostScript file (`../plot/Code.ps`).

### 2. `ConvertOriginTime.f90`
This Fortran program converts the origin time from Julian Day format to a human-readable date format (YYYY-MM-DD).

**Key Features:**
- Accepts origin time as a floating-point value (Julian date).
- Converts it to a standard date format considering leap years.

### 3. `WriteToFile.f90`
This Fortran subroutine writes seismic waveform data to a text file.

**Key Features:**
- Accepts a filename and seismic data structure as input.
- Writes waveform data points to the specified file.

### 4. `PlotWaveformFFT.f90`
This Fortran program reads binary seismic waveform data, performs a Fast Fourier Transform (FFT), and plots both time-domain and frequency-domain representations using PGPLOT.

**Key Features:**
- Reads binary file in little-endian format and extracts header and waveform data.
- Removes the mean value from waveform data to eliminate DC offset.
- Performs FFT and calculates the amplitude spectrum.
- Plots the time-domain waveform and the frequency spectrum, saving them as a PostScript file (`../plot/Code_fft.ps`).

**Explanation: Fast Fourier Transform (FFT)**
- The FFT is an efficient algorithm used to compute the Discrete Fourier Transform (DFT) of a signal, converting time-domain data into the frequency domain.
- In seismic analysis, the FFT helps identify dominant frequencies and analyze the spectral content of the waveform, which can provide insights into the seismic source characteristics and subsurface structures.
- The program plots both the original waveform (time-domain) and its corresponding frequency spectrum (frequency-domain), allowing for a comprehensive analysis of the seismic signal.

### 5. `epicentral_cut.f90`
This Fortran program generates an epicentral cross-sectional plot, showing the depth distribution of seismic events along a specified profile.

**Key Features:**
- Reads earthquake data, including the coordinates and depths of seismic events.
- Projects the earthquake locations along a chosen profile line, aligning them based on the distance from the profile start point.
- Uses PGPLOT to create a contour plot for the surface distribution and a cross-sectional view of earthquake depths.
- Adjusts the symbol size and color based on the depth and magnitude of the events, making it easier to visualize the earthquake characteristics.

**Explanation: Epicentral Cross-Section**
- An epicentral cross-section is a two-dimensional representation of seismic events projected onto a vertical plane along a chosen profile. It helps in visualizing the spatial distribution of earthquake hypocenters (locations beneath the surface).
- By analyzing the cross-section, seismologists can identify patterns such as fault planes, subduction zones, and other structural features related to seismic activity.
- This program plots both the surface distribution (map view) and the depth distribution (cross-sectional view), providing a clear visualization of seismic event locations relative to the Earth's surface.

### 6. `ps2pdf.py`
This Python script automates the conversion of PostScript (`.ps`) files generated by PGPLOT to PDF format.

**Key Features:**
- Searches for `.ps` files in the `../plot/` directory.
- Converts each `.ps` file to a `.pdf` file using the `ps2pdf` utility.
- Outputs the conversion status for each file.

## How to Run

### Prerequisites
- **Fortran Compiler:** Ensure a Fortran compiler (e.g., `gfortran`) is installed.
- **PGPLOT Library:** The Fortran programs rely on the PGPLOT library for plotting. Ensure it is correctly installed and linked.
- **Python 3:** Required for running the `ps2pdf.py` script.
- **Ghostscript (`ps2pdf`):** Required for converting `.ps` files to `.pdf`.

### Compilation
To compile the Fortran programs, use the following commands:

```bash
gfortran -o PlotWaveform PlotWaveform.f90 -lpgplot
gfortran -o ConvertOriginTime ConvertOriginTime.f90
gfortran -c WriteToFile.f90
gfortran -o PlotWaveformFFT PlotWaveformFFT.f90 -lpgplot
gfortran -o epicentral_cut epicentral_cut.f90 -lpgplot