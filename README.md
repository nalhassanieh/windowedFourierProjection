# windowedFourierProjection
1D Scattering from Springs on a String using the Windowed Fourier Projection Method

## Dependencies

This project requires  [Matlab](https://www.mathworks.com/help/install/ug/install-products-with-internet-connection.html)
and [finufft](https://github.com/flatironinstitute/finufft.git). 

## Running Built-in Examples

Enter the experiment number (`expNum`) in the driver file `main.m` to run built-in examples.  
Parameters for example runs are found in `get_experimentParameters(expNum)`:

- `expNum = 1`: Interpolation order 3 — manufactured solution convergence study with 10 randomly located sources in `[-2, 2]`, up to time `3π`.
- `expNum = 2`: Interpolation order 8 — scattering from 150 springs located randomly on `[-2, 2]`, up to time `10π`.
- `expNum = 3`: Interpolation order 8 — scattering from 1000 springs located randomly on `[-2, 2]`, up to time `3π` with narrow incident wave.
- `expNum = 4`: Interpolation order 8 — scattering from 10000 springs located randomly on `[-2, 2]`, up to time `3π` with narrow incident wave.
- `expNum = 5`: Interpolation order 8 — scattering from 2 stiff springs located at `-0.5` and `0.5`, up to time `30π` with wide incident wave (simulates Fabry-Perot interferometer).
- `expNum = 6`: Interpolation order 8 — scattering from 10 stiff springs located uniformly between `-2` to `2`, up to time `40π` with wide incident wave.
- `expNum = 7`: Interpolation order 8 — scattering from 200 springs located uniformly between `-2` to `2`, up to time `40π` with wide incident wave.
- `expNum = 999`: Reserved for running tests and debugging.
