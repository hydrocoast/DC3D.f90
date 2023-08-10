<p align="center">
<img src="https://github.com/hydrocoast/DC3D.f90/blob/master/test/uxuyuz_strike-slip.png", width="250">
<img src="https://github.com/hydrocoast/DC3D.f90/blob/master/test/uxuyuz_dip-slip.png", width="250">
<img src="https://github.com/hydrocoast/DC3D.f90/blob/master/test/uxuyuz_tensile.png", width="250">
</p>  

# DC3D0/DC3D Fortran

## Overview
DC3D.f90 is a subroutine package in Fortran free-form for calculating deformation due to a fault model.  
This is just a file converted from the original code [DC3D0/DC3D](http://www.bosai.go.jp/study/application/dc3d/DC3Dhtml_E.html) in FORTRAN77 fixed-form
and no additional function is given.
Any bug report would be appreciated.

## Usage
To compile [DC3D.f90](https://github.com/hydrocoast/DC3D.f90/blob/master/DC3D.f90) and generate an object file, clone this repository and run the following command on the top directory:
```bash
make obj
```
The shared object can also be generated using `make so`.
Note that a Fortran compiler (e.g. `ifort` and `gfortran`) and an environment variable `FC` are required to run these commands.  

## Tests
You can quickly test the codes as follows:
```bash
make test
```
Files in the directory `test` can help you for further testing.
Some of them need to be run with the Julia language.
```bash
make so
cd test
julia -q
julia> include("./plot_strike-dip-tensile.jl")
```


## Validation
`test/validation_DC3D.jl` compares return values of the original code (FORTRAN77) and the converted one (Fortran90/95).  
<p align="center">
<img src="https://github.com/hydrocoast/DC3D.f90/blob/master/test/diffu.svg", width="900">
</p>  


## License
MIT  

The author got permission in an e-mail massage to upload this repository from National Research Institute for Earth Science and Disaster Resilience (NIED), the copyright holder of the original file.
NIED also takes no responsibility for any damage and loss.

## Author
Takuya Miyashita  
Doctoral student, Kyoto University  
miyashita@hydrocoast.jp  
http://hydrocoast.jp   
