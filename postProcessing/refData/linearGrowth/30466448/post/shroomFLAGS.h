![ This is a header file for shroom.F90.

![ Serial (=0) or MPI (>0) run.
#define useMPI 1

![ Output mode. Note that ADIOS requires MPI.
![   =0 binary.
![   >0 ADIOS.
#define ioMode 1

![ Indicate whether to use adiabatic electrons
![ or adiabatic ions (can't use both), or neither.
![   =0 use fluid species.
![   >0 use adiabatic species.
#define adiabaticElc 1
#define adiabaticIon 0

![ Hyperdiffusion model (not all may be coded).
![ options 1-6 act on the vorticity, and options
![ 7-12 are equivalent but act on phi.
![   =0 no hyperdiffusion.
![   =1,7 constant amplitude & order (Laplacian-based).
![   =2,8 constant amplitude & order (no mixed derivatives).
![   =3,9 constant order, adjusted amplitude (Laplacian-based). 
![   =4,10 constant order, adjusted amplitude (no mixed derivatives). 
![   =5,11 adjusted amplitude & order (Laplacian-based).
![   =6,12 adjusted amplitude & order (no mixed derivatives).
#define HDmodel 9
![ Option for time discretization of hyperdiffusion.
![   =1 Implicit, applied every time step and in every stage.
![   =2 Implicit, applied every time step, end of time step only.
#define HDtimeOP 2 

![ Time stepper.
![   =1 forward Euler.
![   =2 trapezoidal leap-frog.
![   =3 Runge-Kutta 3rd order.
![   =4 Runge-Kutta 4th order.
![   =5 SSP-RK3.
#define STEPPER 4

![ Adjust time step (CURRENTLY NOT WORKING):
![   =0 No.
![   >0 Yes.
#define ADJUSTdt 0

![ FFT mode
![   =0 Gather data in one dimension + use parallel 2D r2c FFT.
![   >0 Use parallel c2c 1D FFTs (data remains decomposed in 2D).
#define FFTmode 1
