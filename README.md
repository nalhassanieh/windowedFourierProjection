# windowedFourierProjection
1D Scattering from Springs on a String using the Windowed Fourier Projection Method

## Dependencies

This project requires  [Matlab](https://www.mathworks.com/help/install/ug/install-products-with-internet-connection.html)
and [finufft](https://github.com/flatironinstitute/finufft.git). 

## Example runs
Enter experiment number (expNum) in the driver file main.m to run built-in examples. 
Parameters for example runs are found in get_experimentParameters(expNum):
expNum = 1, interpolation order 3 manufactured solution convergence study with 10 randomly located sources in [-2,2] up to time 3*pi
expNum = 2, interpolation order 8 scattering from 150 springs located randomly on [-2,2] up to time 10*pi
expNum = 3, interpolation order 8 scattering from 1000 springs located randomly on [-2,2] up to time 3*pi with narrow incident wave
expNum = 4, interpolation order 8 scattering from 10000 springs located randomly on [-2,2] up to time 3*pi with narrow incident wave
expNum = 5, interpolation order 8 scattering from 2 stiff springs located on -0.5 and 0.5 up to time 30*pi with wide incident wave to simulate Fabry-Perot interferometer
expNum = 6, interpolation order 8 scattering from 10 stiff springs located uniformly between -2 to 2 up to time 40*pi final time with wide incident wave
expNum = 7, interpolation order 8 scattering from 200 springs located uniformly between -2 to 2 up to time 40*pi final time with wide incident wave
expNum = 999, experiment number reserved for running tests and debugging 