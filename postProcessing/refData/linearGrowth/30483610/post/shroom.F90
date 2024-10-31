![ ********************************************************************!
![  MuSHrooM
![  Multiscale reduced moment 2D solver.
![  Darin Ernst, Manaure Francisquez.
![
![ Guiding principles for efficiency, reliability and bug minimization:
![ 1) Use an editor that allows you to edit by search & replace, and
![    copy and paste. And become fluent in these.
![ 2) Comments:
![    2a. Start with ![ followed by a space. We can change this simply
![        with a search & replace.
![    2b. Capitalize first letter after ![.
![    2c. End the comment with a period. That way we know the comment
![        finishes there.
![ 3) Commented code is commented just with ! (not ![).
![ 4) Please pull & commit often. Don't wait to commit a whole lot of
![    big changes because it might interfere with someone else editing
![    at the same time.
![ 5) When hardcoding a constant float, make sure to use double
![    precision notation, e.g. 2.0d0, to avoid extra precision
![    conversions.
![
![
![ Additional notes for reader/developer:
![ - At the moment the source code is in a single file (shroom.F90),
![   and there are additional inputs in shroom.in and shroomFLAGS.h.
![ - The source in shroom.F90 is organized as follows.
![     module parameters
![     module fields
![     module MPItools
![       subroutine init_MPI
![       subroutine init_COMMs
![       subroutine terminate_MPI
![     module utilities
![       subroutine abortSimulation
![       subroutine distributeDOFs2D
![       subroutine alloc_1DFourierX 
![       subroutine alloc_1DFourierY
![       subroutine alloc_1DFourierAliasedX
![       subroutine alloc_1DFourierAliasedY
![       subroutine alloc_1DrealX
![       subroutine alloc_1DrealY 
![       subroutine alloc_2DFourier 
![       subroutine alloc_2DFourierAliased 
![       subroutine alloc_2Dreal 
![       subroutine alloc_2DrealAliased 
![     module FLRfunctions
![       subroutine calci0
![       subroutine calci1
![       function besselI0Exp
![       function besselI1Exp
![     module IOtools
![       subroutine init_IO
![       subroutine setupIOfiles
![       subroutine readRestartInputs
![       subroutine setup2DFieldComplex
![       subroutine setup2DFieldReal
![       subroutine setup2DFieldFourierReal
![       subroutine setup2DFieldRealAliased
![       subroutine writeAttributes
![       subroutine out2DFieldComplex
![       subroutine out2DFieldReal
![       subroutine out2DFieldFourierReal
![       subroutine out2DFieldRealAliased
![       subroutine outRestart
![       subroutine readRestartFields
![       subroutine writeFields
![       subroutine terminate_IO
![     module FFTmodule
![       subroutine init_FFTs
![       subroutine init_appendRemoveConjugates
![       subroutine FFT2Da_c2r 
![       subroutine FFT2Da_r2c 
![       subroutine FFT2D_c2r 
![       subroutine FFT2D_r2c 
![       subroutine terminate_FFTs
![       subroutine testFFTs
![     module timeUpdates
![       subroutine compute_steppingFactors 
![       function dtUpdate 
![       subroutine adjustHyperDiff 
![     module initialize
![       subroutine read_inputs 
![       subroutine init_grids
![       subroutine allocate_fields 
![       function closest_power_of_two
![       subroutine init_preFactors
![       subroutine init_random_seed 
![       subroutine setInitialCondition
![       subroutine init_all 
![     module terminate
![       subroutine isNaNorInf
![       subroutine deallocate_parameters 
![       subroutine deallocate_fields 
![       subroutine terminate_all
![     module nonlinearities
![       subroutine placeInAliased
![       subroutine takeFromAliased
![       subroutine poisson_bracket 
![     module timeSteppers
![       subroutine timeRateOfChange 
![       subroutine eulerStep
![       subroutine forwardEuler
![       subroutine trapezoidalLeapFrog 
![       subroutine rk3
![       subroutine rk4
![       subroutine sspRK3
![     program shroom
![
![ ********************************************************************!

module parameters
  ![ This module contains scalar and array variables used throughout the program
  ![ that remain constant in time. Field arrays (e.g. phi) are in the fields module.

![Include user-defined preprocessor variables and directives.
#include "shroomFLAGS.h"

  ![ Universal constant parameters.
  double complex              :: Im1 = dcmplx(0.0d0, 1.0d0)   ![ Imaginary number.
  double precision, parameter :: pi  = 4.0d0*atan(1.0_8)

  ![ Input variable declaration.
  integer           :: NkxG, NkyG             ![ Number of distinct absolute amplitude wavenumbers (de-alised).
  double precision  :: kxMin, kyMin           ![ Minimum finite absolute amplitude wavenumbers.
  double precision  :: dt                     ![ Time step.
  double precision  :: endTime                ![ Absolute end time.
  integer           :: nFrames                ![ Absolute frames to output.
  integer           :: xyProcs(2)             ![ Number of MPI processes along x and y.
  double precision  :: Lnorm                  ![ Reference L divided by reference rho, typically L_n/rho_s.
  double precision  :: omSte                  ![ Coefficient in electron diamagnetic frequency (rho_ref*T_{e0}/(L_n*T_ref)).
  double precision  :: omde                   ![ Coefficient in electron drift frequency (rho_ref*T_{e0}/(L_B*T_ref)).
  double precision  :: tau                    ![ Ion to electron temperature ratio.
  double precision  :: mu                     ![ Square root of ion mass to electron mass ratio.
  double precision  :: deltae                 ![ Electron (T_parallel + T_perp)/T.
  double precision  :: deltaPerpe             ![ Electron T_perp/T.
  double precision  :: eta_e                  ![ Ratio of density gradient scale length to electron temperature gradient scale length.
  double precision  :: deltai                 ![ Ion (T_parallel + T_perp)/T.
  double precision  :: deltaPerpi             ![ Ion T_perp/T.
  double precision  :: eta_i                  ![ Ratio of density gradient scale length to ion temperature gradient scale length.
  double precision  :: lambdaD                ![ Normalized Debye length (for Debye shielding term in Poisson).
  integer           :: initialOp              ![ Initial condition option. =0 noise, =1 power law in k-space.
  double precision  :: initAuxX, initAuxY     ![ Auxiliary variable for initial condition. If initialOp=1, it is the power law exponent. 
  double precision  :: initA                  ![ Amplitude of initial condition.

  ![ Amplitude and order of artificial diffusion.
  double precision :: hDiffne,hDiffTe,hDiffni,hDiffTi
  ![ Hyperdiffusion order (may not be an integer).
  double precision :: hDiffOrder

  namelist/spaceIn/NkxG,NkyG,kxMin,kyMin,dt,endTime,nFrames,xyProcs ![ Namelist with discrete space parameters.
  namelist/sourceIn/Lnorm,omSte,omde,tau,mu,deltae,deltaPerpe,eta_e,deltai,deltaPerpi, &
                    eta_i,lambdaD,hDiffOrder,hDiffne,hDiffTe,hDiffni,hDiffTi  ![ Namelist with source/drive parameters.
  namelist/icsIn/initialOp,initAuxX,initAuxY,initA       ![ Namelist with initial condition parameters.
   
  character(len=100) :: inputFile       ![ Name of input file.
  character(len=100) :: outputDir       ![ Address of output directory.
  character(len=100) :: restartDir      ![ Address of restart directory.

  ![ Boolean indicating whether to use restart0 or restart1 restart file.
  logical            :: useRestartFile0
  ![ When reading inputs, this boolean indicates whether to take
  ![ those from the input file or from the restart file.
  logical, parameter :: favorRestart=.true.

  ![ Other variables derived from inputs.
  double precision  :: dtMax             ![ Maximum allowed time step.
  integer           :: NxaG, NyaG        ![ Number of cells in aliased configuration space.
  integer           :: NkxaG, NkyaG      ![ Number of distinct (absolute) aliased wavenumbers.
  integer           :: NekxG, NekyG      ![ Number of elements in a de-aliased k-space array.
  integer           :: NekxaG, NekyaG    ![ Number of elements in an aliased k-space array.
  integer           :: NxG, NyG          ![ Number of cells in de-aliased configuration space.
  double precision  :: Lx, Ly            ![ Length of real space domain.
  logical           :: isRestart         ![ Is this simulation a restart of a previous one?
  double precision  :: nFramesD          ![ Double precision version of nFrames.

  integer :: myRank      ![ Rank when using MPI, =0 otherwise.
  integer :: cartRank    ![ Rank in cartCOMM when using MPI, =0 otherwise.
  logical :: fftMember   ![ Boolean indicating if this rank participates in FFT.
#if ((useMPI > 0) && (FFTmode == 0))
  integer :: fftRank     ![ Rank within the fftCOMM.
#endif
  integer :: xyIDs(2)    ![ ID (rank) in the 1D yx communicators.
  integer :: Nxa(2)      ![ Number of local cells in aliased configuration space.
  integer :: Nkxa(2)     ![ Number of local distinct (absolute) aliased wavenumbers.
  integer :: Nekx(2)     ![ Number of local elements in a de-aliased k-space array.
  integer :: Nekxa(2)    ![ Number of local elements in an aliased k-space array.
  integer :: Nx(2)       ![ Number of local cells in de-aliased configuration space.

  integer :: firstkG(2)   ![ First element of global de-aliased k-space array stored in this rank.
  integer :: firstkaG(2)  ![ First element of global aliased k-space array stored in this rank.
  integer :: firstxG(2)   ![ First element of global real space array stored in this rank.
  integer :: firstxaG(2)  ![ First element of global aliased real space array stored in this rank.

  ![ Useful prefactors for normalizing FFTs.
  double precision :: rSqrtNxyG, rSqrtNxyaG

  ![ Spacing in each grid
  double precision :: dx, dy
  double precision :: dxa, dya
  double precision :: dkx, dky

  ![ Arrays with global grids.
  double precision, allocatable :: kxaG(:),kyaG(:)                  ![ Arrays for aliased kx/ky values.
  double precision, allocatable :: kxG(:),kyG(:),xG(:),yG(:)        ![ Arrays for de-aliased kx/ky/x/y values.

  ![ Local grids.
  double precision, allocatable :: kxa(:),kya(:)                  ![ Arrays for local aliased kx/ky values.
  double precision, allocatable :: kx(:),ky(:),x(:),y(:)          ![ Arrays for local de-aliased kx/ky/x/y values.

  double complex, allocatable   :: ikxa(:),ikya(:)              ![ Multiply kxa and kya by 1i.

  ![ Magnitude of k^2, i*ky.
  double precision, allocatable :: kSq(:,:)
  double complex, allocatable   :: iky(:)

  ![ Diamagnetic and curvature drift frequency factors.
  double complex, allocatable   :: iOmdElc(:),iOmStElc(:)
  double complex, allocatable   :: iOmdIon(:),iOmStIon(:)

  ![ Quantities dependent on k, used for FLR effects.
  double precision, allocatable :: Gamma0Ion(:,:),Gamma0Elc(:,:)
  double precision, allocatable :: avgJ0Ion(:,:),avgJ0Elc(:,:)
  double precision, allocatable :: hatLapAvgJ0D2Ion(:,:),hatLapAvgJ0D2Elc(:,:)
  double precision, allocatable :: poiSbIon(:,:),poiSbElc(:,:)
  double precision, allocatable :: DbIon(:,:),DbElc(:,:)

  ![ Function mediating the transition between adiabatic and fluid electrons.
  double precision, allocatable :: chiElc(:,:),oMchiElc(:,:)
  ![ Parameter that allows us to use adiabatic ions. 
  double precision :: chiIon

  ![ Reciprocal of the function multiplying the potential in the field equation.
  double precision, allocatable :: rPoiPhikFac(:,:)

  integer :: lc, lckx, lcky, lcx, lcy, lcF,lcTs    ![ Loop controllers.

  ![ Quantities used for Smith+Hammett 1996 hyperviscosity.
  double precision, allocatable :: kxdcSq(:),kydcSq(:),kDkcSq(:,:)    ![ Wavenumbers divided by cutoff wavenumbers along x and y, squared.
  double precision :: kc, kc1p7, kc0p1    ![ Cutoff wavenumber, multiplied by 1.7 and 0.1.

  logical          :: timeStepMore  ![ Boolean indicating time-loop completion.
  double precision :: simTime       ![ Current simulation time.
  integer          :: timeSteps     ![ Current number of completed time steps.

  ![ Time rate at which to print message to screen/log
  ![ given as a (decimal) percentage of endTime.
  double precision :: tRateLogEntry=0.01d0
  integer          :: logEntries    ![ Number of messages printed to log/screen.

  ![ Time rate at which to output, adjust hyperdiffusion and adjust time step.
  double precision :: tRateOutput, tRateAdjustHD, tRateAdjustDt
  double precision :: tTol  ![ Time-tolerance used in time-loop.
  integer          :: framesOut     ![ Frames written so far (not counting ICs).
  integer          :: hdAdjusts     ![ Adjustments to hyperdiffusion so far.
  integer          :: dtAdjusts     ![ Adjustments to time step so far.

  integer, parameter :: intBytes = 4         ![ Bytes in an integer.
  integer, parameter :: dbleBytes = 8        ![ Bytes in a double.
  integer, parameter :: dbleCmplxBytes = 16  ![ Bytes in a double complex.

end module parameters

![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

module fields
  ![ Module declaring the fields variables used through out the code.

  ![ Electrostatic potential.
  double complex, allocatable   :: phik(:,:)
  double precision, allocatable :: phi(:,:)
  ![ Gyroaveraged potential.
  double complex, allocatable   :: elcGyroPhik(:,:),ionGyroPhik(:,:)
  ![ Hat Laplacian divided by 2, operating on the gyroaveraged potential.
  double complex, allocatable   :: hatLapD2elcGyroPhik(:,:),hatLapD2ionGyroPhik(:,:)

  ![ Densities and (perpendicular) temperatures.
  double complex, allocatable   :: denElck(:,:),tempElck(:,:)
  double complex, allocatable   :: denIonk(:,:),tempIonk(:,:)
  double precision, allocatable :: denElc(:,:),tempElc(:,:)
  double precision, allocatable :: denIon(:,:),tempIon(:,:)
  
  ![ Vorticity.
  double precision, allocatable :: vort(:,:)
  double complex, allocatable   :: vortk(:,:)

  ![ Time rates of change.
  double complex, allocatable   :: denElcDot(:,:),tempElcDot(:,:)
  double complex, allocatable   :: denIonDot(:,:),tempIonDot(:,:)
#if ((STEPPER == 3) || (STEPPER == 4) || (STEPPER == 5))
  ![ For RK3 time stepping need df/dt at two more times.
  double complex, allocatable   :: denElcDot1(:,:),tempElcDot1(:,:)
  double complex, allocatable   :: denIonDot1(:,:),tempIonDot1(:,:)
  double complex, allocatable   :: denElcDot2(:,:),tempElcDot2(:,:)
  double complex, allocatable   :: denIonDot2(:,:),tempIonDot2(:,:)
#if (STEPPER == 4)
  ![ For RK4 time stepping need df/dt at three more times.
  double complex, allocatable   :: denElcDot3(:,:),tempElcDot3(:,:)
  double complex, allocatable   :: denIonDot3(:,:),tempIonDot3(:,:)
#endif
#endif

  ![ Poisson Brackets.
  double complex, allocatable   :: PBk(:,:),flrPBk(:,:)
  ![ These fields are used in the Poisson bracket calculation.
  double complex, allocatable   :: fka(:,:), gka(:,:)
  double complex, allocatable   :: fka_x(:,:), gka_x(:,:), &
                                   fka_y(:,:), gka_y(:,:)
  double precision, allocatable :: fRa_x(:,:), gRa_x(:,:), &
                                   fRa_y(:,:), gRa_y(:,:)
  double precision, allocatable :: PBfgRa(:,:)
  double complex, allocatable   :: PBfgka(:,:)

#if (STEPPER == 2)
  ![ Quantities at previous time step (needed by leapfrog-trapezoidal).
  double complex, allocatable   :: phiks(:,:)
  double complex, allocatable   :: denElcks(:,:),tempElcks(:,:)
  double complex, allocatable   :: denIonks(:,:),tempIonks(:,:)
#endif

  ![ Temporary field used by time steppers.
  double complex, allocatable   :: denElckTmp(:,:),tempElckTmp(:,:)
  double complex, allocatable   :: denIonkTmp(:,:),tempIonkTmp(:,:)

  ![ Implicit hyperdiffusion function (depends on hyperdiffusion model).
  double precision, allocatable :: iHDdenElc(:,:),iHDtempElc(:,:)
  double precision, allocatable :: iHDdenIon(:,:),iHDtempIon(:,:)
  ![ Factor applying implicit hyperdiffusion (e.g. 1/(1 - dt*hDiff*kSq)).
  double precision, allocatable :: denElcDiffFac(:,:),tempElcDiffFac(:,:)
  double precision, allocatable :: denIonDiffFac(:,:),tempIonDiffFac(:,:)
  double precision, allocatable :: denElcDotDiffFac(:,:),tempElcDotDiffFac(:,:)
  double precision, allocatable :: denIonDotDiffFac(:,:),tempIonDotDiffFac(:,:)

end module fields

![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

#if (useMPI > 0)
module MPItools
  use parameters
  include 'mpif.h'
  ![ Variables, arrays and routines used to manage MPI communications.
  integer :: errFlagMPI
  integer :: nProcesses, nProcs
  integer :: cartCOMM    ![ 2D Cartesian communicator.
  integer :: xyCOMMs(2)  ![ 1D communicators along each direction.
  ![ Rank of next process and previous process along each direction.
  integer :: prevNeighbor(2), nextNeighbor(2)
#if (FFTmode == 0)
  integer            :: fftCOMM     ![ Communicator of ranks participating in FFT.
  integer, parameter :: fftKyID=0   ![ ky rank ID (within xyCOMMs(1) or along second dimension
                                    ![ of cartCOMM) of processes participating in FFT.
#endif

  contains

  subroutine init_MPI
  ![ Initialize MPI, get rank of this process and total number of processes.
  implicit none
  integer :: threadsProvided
  ![ Initialize MPI
  call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,threadsProvided,errFlagMPI)
  ![ Rank of this MPI process.
  call MPI_COMM_RANK(MPI_COMM_WORLD,myRank,errFlagMPI)
  ![ Number of MPI processes.
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nProcesses,errFlagMPI)
  end subroutine init_MPI
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine init_COMMs
  ![ Initialize sub-communicators used for FFTs and IO.
  implicit none
  logical              :: cartBCs(2),remain(2),reorder
  integer              :: blockID(2),coordTMP
#if (FFTmode == 0)
  integer              :: cartGroup,fftGroup
  integer, allocatable :: fftRanks(:)
#endif

  reorder    = .true.
  ![ MPI_CART boundary conditions along x and y. True=periodic
  cartBCs(1) = .true.    ![ Periodic along x.
  cartBCs(2) = .true.    ![ Periodic along y.

  ![ Create a Cartesian topology.
  call MPI_Cart_create(MPI_COMM_WORLD,2,xyProcs,cartBCs, &
                       reorder,cartCOMM,errFlagMPI)
  ![ Find my coordinate parameters in the Cartesial topology
  ![ and use them to compute the arrays x,y.
  call MPI_Comm_rank(cartCOMM, cartRank, errFlagMPI)
  call MPI_Cart_coords(cartCOMM, myRank, 2, blockID, errFlagMPI)

  ![ X and Y SUBCOMMUNICATORS.
  ![ The first entry in xyCOMMs, xyIDs and prev/next Neighbor
  ![ will correspond to ky, the second to kx.
  do lc = 1,2
    remain       = .false.
    remain(3-lc) = .true.     !.Don't keep lc direction.
    call MPI_Cart_sub(cartCOMM, remain, xyCOMMs(lc), errFlagMPI)
    call MPI_Comm_rank(xyCOMMs(lc), xyIDs(lc), errFlagMPI)
    call MPI_Cart_coords(xyCOMMs(lc), xyIDs(lc), 1, coordTMP, errFlagMPI)
    call MPI_Cart_Shift(xyCOMMs(lc),0,1,prevNeighbor(lc),nextNeighbor(lc),errFlagMPI)
  enddo

#if (FFTmode == 0)
  ![ We need a separate communicator for the processes participating in
  ![ the 2D FFTs decomposed along the second dimension (for Poisson brackets).
  call MPI_Comm_group(cartCOMM, cartGroup, errFlagMPI)

  allocate(fftRanks(xyProcs(1)))
  do lckx = 1,xyProcs(1)
    blockID = (/lckx-1, fftKyID/)
    call MPI_Cart_rank(cartCOMM, blockID, fftRanks(lckx), errFlagMPI)
  enddo
  call MPI_Group_incl(cartGroup, xyProcs(1), fftRanks, fftGroup, errFlagMPI)
  call MPI_Comm_create_group(cartCOMM, fftGroup, 0, fftCOMM, errFlagMPI)
  fftMember = .false.
  if (fftCOMM .ne. MPI_COMM_NULL) fftMember = .true.
  ![ Will need the rank of this process in fftCOMM.
  call MPI_Comm_rank(fftCOMM, fftRank, errFlagMPI)


  deallocate(fftRanks)
  call MPI_Group_free(cartGroup, errFlagMPI)
  call MPI_Group_free(fftGroup, errFlagMPI)
#endif

  end subroutine init_COMMs
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine terminate_MPI
  ![ Close MPI up, thoroughly.
  implicit none
  call MPI_BARRIER(MPI_COMM_WORLD, errFlagMPI) ![ To avoid premature deallocations.

#if (FFTmode == 0)
  if (fftMember) call MPI_Comm_free(fftCOMM, errFlagMPI)
#endif
  call MPI_Comm_free(xyCOMMs(1), errFlagMPI)
  call MPI_Comm_free(xyCOMMs(2), errFlagMPI)
  call MPI_Comm_free(cartCOMM, errFlagMPI)

  call MPI_BARRIER(MPI_COMM_WORLD, errFlagMPI) ![ To avoid premature deallocations.
  call MPI_FINALIZE(errFlagMPI)
  end subroutine terminate_MPI

end module MPItools
#endif

![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

module utilities
  use parameters
  use fields

  ![ The following interfaces allow calling by a single function name
  ![ in order to allocate both real and complex arrays.

  ![ Routines for allocating 1D arrays.
  interface alloc_1DFourierX
    module procedure alloc_1DFourierDbleX
    module procedure alloc_1DFourierCmplxX
  end interface alloc_1DFourierX
  interface alloc_1DFourierY
    module procedure alloc_1DFourierDbleY
    module procedure alloc_1DFourierCmplxY
  end interface alloc_1DFourierY
  interface alloc_1DFourierAliasedX
    module procedure alloc_1DFourierAliasedDbleX
    module procedure alloc_1DFourierAliasedCmplxX
  end interface alloc_1DFourierAliasedX
  interface alloc_1DFourierAliasedY
    module procedure alloc_1DFourierAliasedDbleY
    module procedure alloc_1DFourierAliasedCmplxY
  end interface alloc_1DFourierAliasedY
  interface alloc_1DrealX
    module procedure alloc_1DrealDbleX
    module procedure alloc_1DrealCmplxX
  end interface alloc_1DrealX
  interface alloc_1DrealY
    module procedure alloc_1DrealDbleY
    module procedure alloc_1DrealCmplxY
  end interface alloc_1DrealY
  ![ Routines for allocating 2D arrays.
  interface alloc_2DFourier
    module procedure alloc_2DFourierDble
    module procedure alloc_2DFourierCmplx
  end interface alloc_2DFourier
  interface alloc_2DFourierAliased
    module procedure alloc_2DFourierAliasedDble
    module procedure alloc_2DFourierAliasedCmplx
  end interface alloc_2DFourierAliased
  interface alloc_2Dreal
    module procedure alloc_2DrealDble
    module procedure alloc_2DrealCmplx
  end interface alloc_2Dreal
  interface alloc_2DrealAliased
    module procedure alloc_2DrealAliasedDble
    module procedure alloc_2DrealAliasedCmplx
  end interface alloc_2DrealAliased

  contains

  subroutine abortSimulation(errorString)
  ![ Abrutly terminate the simulation (presumably due to an error).
#if (useMPI > 0)
  use MPItools
#endif
  implicit none
  character(len=*), intent(in) :: errorString
  print*,'  ',errorString
#if (useMPI > 0)
  call MPI_ABORT(MPI_COMM_WORLD, errFlagMPI)
#else
  error stop
#endif
  end subroutine abortSimulation
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine distributeDOFs2D(IDs,Nprocs,Nglobal,Nlocal,firstElement)
  ![ Distribute Nglobal degrees of freedom in each direction among
  ![ Nprocs processes with ranks given by IDs. The local number of
  ![ degrees of freedom is given in Nlocal, and firstElement is the index
  ![ of the first element in the global array stored in this process.
  implicit none
  integer, intent(in)    :: IDs(2), Nprocs(2), Nglobal(2)
  integer, intent(inout) :: Nlocal(2),firstElement(2)
  integer                :: remainingDOFs(2)

  Nlocal        = [ floor(dble(Nglobal(1))/dble(Nprocs(1))), &
                    floor(dble(Nglobal(2))/dble(Nprocs(2))) ]
  remainingDOFs = [ Nglobal(1) - Nlocal(1)*Nprocs(1), &
                    Nglobal(2) - Nlocal(2)*Nprocs(2) ]

  ![ Store the index of the first element allocated to this process.
  firstElement = [IDs(1)*Nlocal(1), IDs(2)*Nlocal(2)]

  do lcF = 1,2
    do while (remainingDOFs(lcF) > 0)
      do lc = 0,Nprocs(lcF)-1
        remainingDOFs(lcF) = remainingDOFs(lcF)-1
        if (lc == IDs(lcF)) then
          Nlocal(lcF) = Nlocal(lcF)+1
        elseif (lc < IDs(lcF)) then
          firstElement(lcF) = firstElement(lcF)+1
        endif
        if (remainingDOFs(lcF)==0) exit
      enddo
    enddo
  enddo

  firstElement = firstElement + 1

  end subroutine distributeDOFs2D
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_1DFourierDbleX(aIn)
  ![ Allocate a 1D double precision array with the dimension of de-aliased Fourier space along kx.
  implicit none
  double precision, allocatable, intent(inout) :: aIn(:)
  allocate(aIn(Nekx(2)))
  end subroutine alloc_1DFourierDbleX
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_1DFourierDbleY(aIn)
  ![ Allocate a 1D double precision array with the dimension of de-aliased Fourier space along ky.
  implicit none
  double precision, allocatable, intent(inout) :: aIn(:)
  allocate(aIn(Nekx(1)))
  end subroutine alloc_1DFourierDbleY
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_1DFourierAliasedDbleX(aIn)
  ![ Allocate a 1D double precision array with the dimension of aliased Fourier space along kx.
  implicit none
  double precision, allocatable, intent(inout) :: aIn(:)
  allocate(aIn(Nekxa(2)))
  end subroutine alloc_1DFourierAliasedDbleX
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_1DFourierAliasedDbleY(aIn)
  ![ Allocate a 1D double precision array with the dimension of aliased Fourier space along ky.
  implicit none
  double precision, allocatable, intent(inout) :: aIn(:)
#if ((useMPI > 0) && (FFTmode==0))
  allocate(aIn(NekyaG))
#else
  allocate(aIn(Nekxa(1)))
#endif
  end subroutine alloc_1DFourierAliasedDbleY
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_1DrealDbleX(aIn)
  ![ Allocate a 1D double precision array with the dimension of de-aliased space along x.
  implicit none
  double precision, allocatable, intent(inout) :: aIn(:)
  allocate(aIn(Nx(2)))
  end subroutine alloc_1DrealDbleX
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_1DrealDbleY(aIn)
  ![ Allocate a 1D double precision array with the dimension of de-aliased space along y.
  implicit none
  double precision, allocatable, intent(inout) :: aIn(:)
  allocate(aIn(Nx(1)))
  end subroutine alloc_1DrealDbleY
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_1DFourierCmplxX(aIn)
  ![ Allocate a 1D complex array with the dimension of de-aliased Fourier space along kx.
  implicit none
  double complex, allocatable, intent(inout) :: aIn(:)
  allocate(aIn(Nekx(2)))
  end subroutine alloc_1DFourierCmplxX
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_1DFourierCmplxY(aIn)
  ![ Allocate a 1D complex array with the dimension of de-aliased Fourier space along ky.
  implicit none
  double complex, allocatable, intent(inout) :: aIn(:)
  allocate(aIn(Nekx(1)))
  end subroutine alloc_1DFourierCmplxY
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_1DFourierAliasedCmplxX(aIn)
  ![ Allocate a 1D complex array with the dimension of aliased Fourier space along kx.
  implicit none
  double complex, allocatable, intent(inout) :: aIn(:)
  allocate(aIn(Nekxa(2)))
  end subroutine alloc_1DFourierAliasedCmplxX
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_1DFourierAliasedCmplxY(aIn)
  ![ Allocate a 1D complex array with the dimension of aliased Fourier space along ky.
  implicit none
  double complex, allocatable, intent(inout) :: aIn(:)
#if ((useMPI > 0) && (FFTmode == 0))
  allocate(aIn(NekyaG))
#else
  allocate(aIn(Nekxa(1)))
#endif
  end subroutine alloc_1DFourierAliasedCmplxY
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_1DrealCmplxX(aIn)
  ![ Allocate a 1D complex array with the dimension of de-aliased space along x.
  implicit none
  double complex, allocatable, intent(inout) :: aIn(:)
  allocate(aIn(Nx(2)))
  end subroutine alloc_1DrealCmplxX
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_1DrealCmplxY(aIn)
  ![ Allocate a 1D complex array with the dimension of de-aliased space along y.
  implicit none
  double complex, allocatable, intent(inout) :: aIn(:)
  allocate(aIn(Nx(1)))
  end subroutine alloc_1DrealCmplxY
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_2DFourierDble(aIn)
  ![ Allocate a 2D double precision array with the dimensions of de-aliased Fourier space.
  implicit none
  double precision, allocatable, intent(inout) :: aIn(:,:)
  allocate(aIn(Nekx(1),Nekx(2)))
  end subroutine alloc_2DFourierDble
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_2DFourierAliasedDble(aIn)
  ![ Allocate a 2D double precision array with the dimensions of aliased Fourier space.
  implicit none
  double precision, allocatable, intent(inout) :: aIn(:,:)
#if ((useMPI > 0) && (FFTmode == 0))
  allocate(aIn(NekyaG,Nekxa(2)))
#else
  allocate(aIn(Nekxa(1),Nekxa(2)))
#endif
  end subroutine alloc_2DFourierAliasedDble
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_2DrealDble(aIn)
  ![ Allocate a 2D double precision array with the dimensions of de-aliased real space.
  implicit none
  double precision, allocatable, intent(inout) :: aIn(:,:)
  allocate(aIn(Nx(2),Nx(1)))
  end subroutine alloc_2DrealDble
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_2DrealAliasedDble(aIn)
  ![ Allocate a 2D double precision array with the dimensions of aliased real space.
  implicit none
  double precision, allocatable, intent(inout) :: aIn(:,:)
#if ((useMPI > 0) && (FFTmode == 0))
  allocate(aIn(NyaG,Nxa(2)))
#else
  allocate(aIn(Nxa(2),Nxa(1)))
#endif
  end subroutine alloc_2DrealAliasedDble
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_2DFourierCmplx(aIn)
  ![ Allocate a 2D complex array with the dimensions of de-aliased Fourier space.
  implicit none
  double complex, allocatable, intent(inout) :: aIn(:,:)
  allocate(aIn(Nekx(1),Nekx(2)))
  end subroutine alloc_2DFourierCmplx
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_2DFourierAliasedCmplx(aIn)
  ![ Allocate a 2D complex array with the dimensions of aliased Fourier space.
  implicit none
  double complex, allocatable, intent(inout) :: aIn(:,:)
#if ((useMPI > 0) && (FFTmode == 0))
  allocate(aIn(NekyaG,Nekxa(2)))
#else
  allocate(aIn(Nekxa(1),Nekxa(2)))
#endif
  end subroutine alloc_2DFourierAliasedCmplx
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_2DrealCmplx(aIn)
  ![ Allocate a 2D complex array with the dimensions of de-aliased real space.
  implicit none
  double complex, allocatable, intent(inout) :: aIn(:,:)
  allocate(aIn(Nx(2),Nx(1)))
  end subroutine alloc_2DrealCmplx
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine alloc_2DrealAliasedCmplx(aIn)
  ![ Allocate a 2D complex array with the dimensions of aliased real space.
  implicit none
  double complex, allocatable, intent(inout) :: aIn(:,:)
#if ((useMPI > 0) && (FFTmode == 0))
  allocate(aIn(NyaG,Nxa(2)))
#else
  allocate(aIn(Nxa(2),Nxa(1)))
#endif
  end subroutine alloc_2DrealAliasedCmplx

end module utilities

![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

module FLRfunctions
  ![ Calculate things arising from finite Larmor radius effects (FLR),
  ![ typically depending on Bessel functions.

  contains

  subroutine calci0( arg, result, jint )
  ![ CALCI0 computes various I0 Bessel functions.
  ![ Origin:
  ![   SPECFUN, see https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html.
  ![ Discussion:
  ![   This routine computes modified Bessel functions of the first kind
  ![   and order zero, I0(X) and EXP(-ABS(X))*I0(X), for real
  ![   arguments X.
  ![   The main computation evaluates slightly modified forms of
  ![   minimax approximations generated by Blair and Edwards, Chalk
  ![   River (Atomic Energy of Canada Limited) Report AECL-4928,
  ![   October, 1974.
  ![ Licensing:
  ![   This code is distributed under the GNU LGPL license.
  ![ Modified:
  ![   03 April 2007
  ![ Author:
  ![   Original FORTRAN77 version by William Cody, Laura Stoltz.
  ![   FORTRAN90 version by John Burkardt.
  ![ Parameters:
  ![   - Input, real ( kind = 8 ) ARG, the argument.  If JINT = 1, then
  ![     the argument must be less than XMAX.
  ![   - Output, real ( kind = 8 ) RESULT, the value of the function,
  ![     which depends on the input value of JINT:
  ![     1, RESULT = I0(x);
  ![     2, RESULT = exp(-x) * I0(x);
  ![   - Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
  ![     1, I0(x);
  ![     2, exp(-x) * I0(x);
  implicit none
  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) exp40
  real ( kind = 8 ) forty
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jint
  real ( kind = 8 ) one5
  real ( kind = 8 ) p(15)
  real ( kind = 8 ) pp(8)
  real ( kind = 8 ) q(5)
  real ( kind = 8 ) qq(7)
  real ( kind = 8 ) result
  real ( kind = 8 ) rec15
  real ( kind = 8 ) sump
  real ( kind = 8 ) sumq
  real ( kind = 8 ) two25
  real ( kind = 8 ) x
  real ( kind = 8 ) xinf
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xsmall
  real ( kind = 8 ) xx
  ![ Mathematical constants.
  data one5 /15.0d0/
  data exp40 /2.353852668370199854d17/
  data forty /40.0d0/
  data rec15 /6.6666666666666666666d-2/
  data two25 /225.0d0/
  ![ Machine-dependent constants.
  data xsmall /5.55d-17/
  data xinf /1.79d308/
  data xmax /713.986d0/
  ![ Coefficients for XSMALL <= ABS(ARG) < 15.0.
  data  p/-5.2487866627945699800d-18,-1.5982226675653184646d-14, &
          -2.6843448573468483278d-11,-3.0517226450451067446d-08, &
          -2.5172644670688975051d-05,-1.5453977791786851041d-02, &
          -7.0935347449210549190d+00,-2.4125195876041896775d+03, &
          -5.9545626019847898221d+05,-1.0313066708737980747d+08, &
          -1.1912746104985237192d+10,-8.4925101247114157499d+11, &
          -3.2940087627407749166d+13,-5.5050369673018427753d+14, &
          -2.2335582639474375249d+15/
  data  q/-3.7277560179962773046d+03, 6.5158506418655165707d+06, &
          -6.5626560740833869295d+09, 3.7604188704092954661d+12, &
          -9.7087946179594019126d+14/
  ![ Coefficients for 15.0 <= ABS(ARG).
  data pp/-3.9843750000000000000d-01, 2.9205384596336793945d+00, &
          -2.4708469169133954315d+00, 4.7914889422856814203d-01, &
          -3.7384991926068969150d-03,-2.6801520353328635310d-03, &
           9.9168777670983678974d-05,-2.1877128189032726730d-06/
  data qq/-3.1446690275135491500d+01, 8.5539563258012929600d+01, &
          -6.0228002066743340583d+01, 1.3982595353892851542d+01, &
          -1.1151759188741312645d+00, 3.2547697594819615062d-02, &
          -5.5194330231005480228d-04/

  x = abs ( arg )

  if ( x < xsmall ) then
    result = 1.0D+00
  else if ( x < one5 ) then
    ![ XSMALL <= ABS(ARG) < 15.0.

    xx = x * x
    sump = p(1)
    do i = 2, 15
      sump = sump * xx + p(i)
    end do
    xx = xx - two25

    sumq = (((( &
        xx + q(1) ) &
      * xx + q(2) ) &
      * xx + q(3) ) &
      * xx + q(4) ) &
      * xx + q(5)

    result = sump / sumq

    if ( jint == 2 ) then
      result = result * exp ( - x )
    end if

  else if ( one5 <= x ) then

    if ( jint == 1 .and. xmax < x ) then
      result = xinf
    else
      ![ 15.0 <= ABS(ARG).
      xx = 1.0D+00 / x - rec15

      sump = (((((( &
               pp(1) &
        * xx + pp(2) ) &
        * xx + pp(3) ) &
        * xx + pp(4) ) &
        * xx + pp(5) ) &
        * xx + pp(6) ) &
        * xx + pp(7) ) &
        * xx + pp(8)

      sumq = (((((( &
          xx + qq(1) ) &
        * xx + qq(2) ) &
        * xx + qq(3) ) &
        * xx + qq(4) ) &
        * xx + qq(5) ) &
        * xx + qq(6) ) &
        * xx + qq(7)

      result = sump / sumq

      if ( jint == 2 ) then
        result = ( result - pp(1) ) / sqrt ( x )
      else
        ![ Calculation reformulated to avoid premature overflow.
        if ( x .le.( xmax - one5 ) ) then
          a = exp ( x )
          b = 1.0D+00
        else
          a = exp ( x - forty )
          b = exp40
        end if

        result = ( ( result * a - pp(1) * a ) / sqrt ( x ) ) * b

      end if

    end if

  end if

  return
  end subroutine calci0
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine calci1( arg, result, jint )
  ![ CALCI1 computes various I1 Bessel functions.
  ![ Origin:
  ![   SPECFUN, see https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html.
  ![ Discussion:
  ![   This routine computes modified Bessel functioons of the first kind
  ![   and order one, I1(X) and EXP(-ABS(X))*I1(X), for real
  ![   arguments X.
  ![   The main computation evaluates slightly modified forms of
  ![   minimax approximations generated by Blair and Edwards, Chalk
  ![   River (Atomic Energy of Canada Limited) Report AECL-4928,
  ![   October, 1974.
  ![ Licensing:
  ![   This code is distributed under the GNU LGPL license.
  ![ Modified:
  ![   03 April 2007
  ![ Author:
  ![   Original FORTRAN77 version by William Cody, Laura Stoltz.
  ![   FORTRAN90 version by John Burkardt.
  ![ Parameters:
  ![   - Input, real ( kind = 8 ) ARG, the argument.  If JINT = 1, then
  ![     the argument must be less than XMAX.
  ![   - Output, real ( kind = 8 ) RESULT, the value of the function,
  ![     which depends on the input value of JINT:
  ![     1, RESULT = I1(x);
  ![     2, RESULT = exp(-x) * I1(x);
  ![   - Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
  ![     1, I1(x);
  ![     2, exp(-x) * I1(x);
  implicit none
  real ( kind = 8 ) a
  real ( kind = 8 ) arg
  real ( kind = 8 ) b
  real ( kind = 8 ) exp40
  real ( kind = 8 ) forty
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jint
  real ( kind = 8 ) one5
  real ( kind = 8 ) p(15)
  real ( kind = 8 ) pbar
  real ( kind = 8 ) pp(8)
  real ( kind = 8 ) q(5)
  real ( kind = 8 ) qq(6)
  real ( kind = 8 ) rec15
  real ( kind = 8 ) result
  real ( kind = 8 ) sump
  real ( kind = 8 ) sumq
  real ( kind = 8 ) two25
  real ( kind = 8 ) x
  real ( kind = 8 ) xinf
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xsmall
  real ( kind = 8 ) xx
  ![ Mathematical constants.
  data one5 / 15.0d0 /
  data exp40 / 2.353852668370199854d17 /
  data forty / 40.0d0 /
  data rec15 / 6.6666666666666666666d-2 /
  data two25 / 225.0d0 /
  ![ Machine-dependent constants.
  data xsmall /5.55d-17/
  data xinf /1.79d308/
  data xmax /713.987d0/
  ![ Coefficients for XSMALL <= ABS(ARG) < 15.0.
  data p/-1.9705291802535139930d-19,-6.5245515583151902910d-16, &
         -1.1928788903603238754d-12,-1.4831904935994647675d-09, &
         -1.3466829827635152875d-06,-9.1746443287817501309d-04, &
         -4.7207090827310162436d-01,-1.8225946631657315931d+02, &
         -5.1894091982308017540d+04,-1.0588550724769347106d+07, &
         -1.4828267606612366099d+09,-1.3357437682275493024d+11, &
         -6.9876779648010090070d+12,-1.7732037840791591320d+14, &
         -1.4577180278143463643d+15/
  data q/-4.0076864679904189921d+03, 7.4810580356655069138d+06, &
         -8.0059518998619764991d+09, 4.8544714258273622913d+12, &
         -1.3218168307321442305d+15/
  ![ Coefficients for 15.0 <= ABS(ARG).
  data pp/-6.0437159056137600000d-02, 4.5748122901933459000d-01, &
          -4.2843766903304806403d-01, 9.7356000150886612134d-02, &
          -3.2457723974465568321d-03,-3.6395264712121795296d-04, &
           1.6258661867440836395d-05,-3.6347578404608223492d-07/
  data qq/-3.8806586721556593450d+00, 3.2593714889036996297d+00, &
          -8.5017476463217924408d-01, 7.4212010813186530069d-02, &
          -2.2835624489492512649d-03, 3.7510433111922824643d-05/
  data pbar/3.98437500d-01/

  x = abs ( arg )
  if ( x < xsmall ) then
    ![ Return for ABS(ARG) < XSMALL.
    result = 0.5D+00 * x
  else if ( x < one5 ) then
    ![ XSMALL <= ABS(ARG) < 15.0.

    xx = x * x
    sump = p(1)
    do j = 2, 15
      sump = sump * xx + p(j)
    end do
    xx = xx - two25

    sumq = (((( &
        xx + q(1) ) &
      * xx + q(2) ) &
      * xx + q(3) ) &
      * xx + q(4) ) &
      * xx + q(5)

    result = ( sump / sumq ) * x

    if ( jint == 2 ) then
      result = result * exp ( -x )
    end if

  else if ( jint == 1 .and. xmax < x ) then

    result = xinf

  else
    ![ 15.0 <= ABS(ARG).
    xx = 1.0D+00 / x - rec15

    sump = (((((( &
             pp(1) &
      * xx + pp(2) ) &
      * xx + pp(3) ) &
      * xx + pp(4) ) &
      * xx + pp(5) ) &
      * xx + pp(6) ) &
      * xx + pp(7) ) &
      * xx + pp(8)

    sumq = ((((( &
        xx + qq(1) ) &
      * xx + qq(2) ) &
      * xx + qq(3) ) &
      * xx + qq(4) ) &
      * xx + qq(5) ) &
      * xx + qq(6)

    result = sump / sumq

    if ( jint /= 1 ) then
      result = ( result + pbar ) / sqrt ( x )
    else
      ![ Calculation reformulated to avoid premature overflow.
      if ( xmax - one5 < x ) then
        a = exp ( x - forty )
        b = exp40
      else
        a = exp ( x )
        b = 1.0D+00
      end if

      result = ( ( result * a + pbar * a ) / sqrt ( x ) ) * b

    end if
  end if

  if ( arg < 0.0D+00 ) then
    result = -result
  end if

  return
  end subroutine calci1
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function besselI0Exp(x)
  ![ besselI0Exp (originally named BESEI0) evaluates the exponentially scaled Bessel I0(X) function.
  ![ Origin:
  ![   SPECFUN, see https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html.
  ![ Discussion:
  ![   This routine computes approximate values for the modified Bessel
  ![   function of the first kind of order zero multiplied by EXP(-ABS(X)).
  ![ Licensing:
  ![   This code is distributed under the GNU LGPL license.
  ![ Modified:
  ![   03 April 2007
  ![ Author:
  ![   Original FORTRAN77 version by William Cody.
  ![   FORTRAN90 version by John Burkardt.
  ![ Parameters:
  ![   - Input, real ( kind = 8 ) X, the argument of the function.
  ![   - Output, real ( kind = 8 ) BESEI0, the value of the function.
  implicit none
  real ( kind = 8 ) besselI0Exp
  integer ( kind = 4 ) jint
  real ( kind = 8 ) result
  real ( kind = 8 ) x
  
  jint = 2
  call calci0 ( x, result, jint )
  besselI0Exp = result
  
  return
  end function besselI0Exp
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function besselI1Exp(x)
  ![ besselI1Exp (originally named BESEI1) evaluates the exponentially scaled Bessel I1(X) function.
  ![ Origin:
  ![   SPECFUN, see https://people.sc.fsu.edu/~jburkardt/f_src/special_functions/special_functions.html.
  ![ Discussion:
  ![   This routine computes approximate values for the
  ![   modified Bessel function of the first kind of order one
  ![   multiplied by EXP(-ABS(X)).
  ![ Licensing:
  ![   This code is distributed under the GNU LGPL license.
  ![ Modified:
  ![   03 April 2007
  ![ Author:
  ![   Original FORTRAN77 version by William Cody.
  ![   FORTRAN90 version by John Burkardt.
  ![ Parameters:
  ![   - Input, real ( kind = 8 ) X, the argument of the function.
  ![   - Output, real ( kind = 8 ) BESEI1, the value of the function.
  implicit none
  real ( kind = 8 ) besselI1Exp
  integer ( kind = 4 ) jint
  real ( kind = 8 ) result
  real ( kind = 8 ) x

  jint = 2
  call calci1 ( x, result, jint )
  besselI1Exp = result

  return
  end function besselI1Exp
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

end module FLRfunctions

![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

module IOtools
  ![ Initialize, read, write and close I/O files.
  ![ Both ADIOS and binary IO are supported.
  use parameters
  use fields
#if (ioMode > 0)
  use MPItools
  use adios_write_mod
  use adios_read_mod
  character(len=3), parameter  :: sFileExt='.bp'
  character(len=20), parameter :: bpSchemaVersion='1.1'
  integer, parameter           :: bpMaxBufferSize = 50    ![ Maximum buffer size in MB.
  integer                      :: bpReadMethod=ADIOS_READ_METHOD_BP
  integer                      :: errFlagADIOS
  ![ Below bpG stands for bp Group.
  integer*8 :: bpHandle
  integer*8 :: bpGRestart,bpGsizeRestart,bpGtotalSizeRestart
  ![ The following has two entries, one for 2DComplex and one for 2DReal,
  ![ but users should refer to them with complex2Didx and real2Didx for safety.
  integer*8 :: bpG2D(4),bpGsize2D(4),bpGtotalSize2D(4)
  integer, parameter          :: complex2Didx=1, real2Didx=2, fourierReal2Didx=3, realAliased2Didx=4
#else
  character(len=4), parameter :: sFileExt='.bin'
  integer, parameter          :: restartFileIdx=10
#endif
  character(len=7), parameter :: restartFileName='restart'
  character(len=150)          :: restartFile

  ![ Package write routines so no distinction between types has to be
  ![ made in the code that calls the routine.
  interface out2DField
    module procedure out2DFieldComplex
    module procedure out2DFieldReal
  end interface out2DField

  contains

  subroutine init_IO
  ![ Initialize IO interface.
  implicit none
#if (ioMode == 0)
  if (myRank == 0) then
    write(*,*) ' '
    write(*,'(" Performing IO with binary files. Yay. ")')
    write(*,'(" Output directory: ",A50)') trim(outputDir)
    write(*,*) ' '
  endif
#else
  if (myRank == 0) then
    write(*,*) ' '
    write(*,'(" Performing IO with ADIOS. ")')
    write(*,'(" Output directory: ",A50)') trim(outputDir)
    write(*,*) ' '
  endif

  ![ Initialize ADIOS output interface without XML config file.
  call adios_init_noxml(MPI_COMM_WORLD,errFlagADIOS)
  ![ Maximum buffer size used for read/write.
  call adios_set_max_buffer_size(bpMaxBufferSize)

  ![ Initialize ADIOS read interface.
  call adios_read_init_method(bpReadMethod,MPI_COMM_WORLD,"verbose=3",errFlagADIOS)
#endif
  end subroutine init_IO
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine setupIOfiles
  ![ Create output files (isRestart=T or isRestart=F but file is missing)
  ![ or connect with existing output file (isRestart=T).
  implicit none
#if (ioMode > 0)
  character(len=5)  :: cN1L,cN2L,cN1G,cN2G,cN1off,cN2off
  character(len=13) :: cDimsL,cDimsG,cOffsets 
  integer*8         :: bpVarID
  integer           :: i
#endif

  if (myRank == 0) then
    write(*,'(" Creating output data files... ")')
    write(*,*) ' '
  endif

#if (ioMode != 0)
  ![ .......... System for ADIOS files containing 2D complex and real fields ......... ]!
  do i = 1,size(bpG2D) 
    ![ Group used to write 2D field.
    if (i==complex2Didx) then
      call adios_declare_group(bpG2D(i),"complex2d","",ADIOS_STAT_DEFAULT,errFlagADIOS)
      ![ Convert local and global dimensions of field data to string.
      write(cN1L,'(i0)') Nekx(1)
      write(cN2L,'(i0)') Nekx(2)
      write(cN1G,'(i0)') NekyG
      write(cN2G,'(i0)') NekxG
      write(cN1off,'(i0)') firstkG(1)-1
      write(cN2off,'(i0)') firstkG(2)-1
    elseif (i==real2Didx) then
      call adios_declare_group(bpG2D(i),"real2d","",ADIOS_STAT_DEFAULT,errFlagADIOS)
      ![ Convert local and global dimensions of field data to string.
      write(cN1L,'(i0)') Nx(2)
      write(cN2L,'(i0)') Nx(1)
      write(cN1G,'(i0)') NxG
      write(cN2G,'(i0)') NyG
      write(cN1off,'(i0)') firstxG(2)-1
      write(cN2off,'(i0)') firstxG(1)-1
    elseif (i==fourierReal2Didx) then
      call adios_declare_group(bpG2D(i),"fourierReal2d","",ADIOS_STAT_DEFAULT,errFlagADIOS)
      ![ Convert local and global dimensions of field data to string.
      write(cN1L,'(i0)') Nekx(1)
      write(cN2L,'(i0)') Nekx(2)
      write(cN1G,'(i0)') NekyG
      write(cN2G,'(i0)') NekxG
      write(cN1off,'(i0)') firstkG(1)-1
      write(cN2off,'(i0)') firstkG(2)-1
    elseif (i==realAliased2Didx) then
      call adios_declare_group(bpG2D(i),"realAliased2d","",ADIOS_STAT_DEFAULT,errFlagADIOS)
      ![ Convert local and global dimensions of field data to string.
#if ((useMPI > 0) && (FFTmode == 0))
      write(cN1L,'(i0)') NyaG
      write(cN2L,'(i0)') Nxa(2)
      write(cN1G,'(i0)') NyaG
      write(cN2G,'(i0)') NxaG
      write(cN1off,'(i0)') 0
      write(cN2off,'(i0)') firstxaG(2)-1
#else
      write(cN1L,'(i0)') Nxa(2)
      write(cN2L,'(i0)') Nxa(1)
      write(cN1G,'(i0)') NxaG
      write(cN2G,'(i0)') NyaG
      write(cN1off,'(i0)') firstxaG(2)-1
      write(cN2off,'(i0)') firstxaG(1)-1
#endif
    endif
    cDimsL   = cN1L//","//cN2L
    cDimsG   = cN1G//","//cN2G
    cOffsets = cN1off//","//cN2off
    call adios_select_method(bpG2D(i),"MPI","","",errFlagADIOS)
    ![ Define the variables that we wish to include in the output file.
    ![ For each new variable defined (via adios_define_var) one must
    ![ include the appropriate contribution to bpGsize below.
    ![ Simulation time for this frame.
    call adios_define_var(bpG2D(i),"simulationTime","",adios_double,"","","",bpVarID)
    ![ Time steps completed by this frame.
    call adios_define_var(bpG2D(i),"timeSteps","",adios_integer,"","","",bpVarID)
    ![ (Absolute) number of frames outputted so far.
    call adios_define_var(bpG2D(i),"framesOut","",adios_integer,"","","",bpVarID)
    ![ (Absolute) number of times hyperdiffusion has been adjusted.
    call adios_define_var(bpG2D(i),"hdAdjusts","",adios_integer,"","","",bpVarID)
    ![ (Absolute) number of times dt has been adjusted.
    call adios_define_var(bpG2D(i),"dtAdjusts","",adios_integer,"","","",bpVarID)
    ![ Time steps size.
    call adios_define_var(bpG2D(i),"dt","",adios_double,"","","",bpVarID)
    ![ Hyper-diffusion order.
    call adios_define_var(bpG2D(i),"hDiffOrder","",adios_double,"","","",bpVarID)
    ![ Hyper-diffusion coefficient.
    call adios_define_var(bpG2D(i),"hDiffne","",adios_double,"","","",bpVarID)
    call adios_define_var(bpG2D(i),"hDiffTe","",adios_double,"","","",bpVarID)
    call adios_define_var(bpG2D(i),"hDiffni","",adios_double,"","","",bpVarID)
    call adios_define_var(bpG2D(i),"hDiffTi","",adios_double,"","","",bpVarID)
    ![ Two-dimensional field frame.
    if (i==complex2Didx) then
      call adios_define_var(bpG2D(i),"complexField","",adios_double_complex, &
                            trim(cDimsL),trim(cDimsG),trim(cOffsets),bpVarID)
    elseif (i==real2Didx) then
      call adios_define_var(bpG2D(i),"realField","",adios_double, &
                            trim(cDimsL),trim(cDimsG),trim(cOffsets),bpVarID)
    elseif (i==fourierReal2Didx) then
      call adios_define_var(bpG2D(i),"fourierRealField","",adios_double, &
                            trim(cDimsL),trim(cDimsG),trim(cOffsets),bpVarID)
    elseif (i==realAliased2Didx) then
      call adios_define_var(bpG2D(i),"realAliasedField","",adios_double, &
                            trim(cDimsL),trim(cDimsG),trim(cOffsets),bpVarID)
    endif
    ![ Add other quantities defining this simulation as attributes in ADIOS file.
    call writeAttributes(bpG2D(i),i)

    call adios_define_schema_version(bpG2D(i), bpSchemaVersion)
  enddo

  bpGsize2D(complex2Didx) = dbleBytes*7 &                   ![ Doubles for simulationTime, dt, hDiffOrder, hDiffs.
                           +intBytes*4 &                    ![ Integer for timeSteps, framesOut, hdAdjusts, dtAdjusts.
                           +dbleCmplxBytes*Nekx(1)*Nekx(2)  ![ Complex field data.

  bpGsize2D(real2Didx) = dbleBytes*7 &          ![ Doubles for simulationTime, dt, hDiffOrder, hDiffs.
                        +intBytes*4 &           ![ Integer for timeSteps, framesOut, hdAdjusts, dtAdjusts.
                        +dbleBytes*Nx(2)*Nx(1)  ![ Real field data.

  bpGsize2D(fourierReal2Didx) = dbleBytes*7 &              ![ Doubles for simulationTime, dt, hDiffOrder, hDiffs.
                               +intBytes*4 &               ![ Integer for timeSteps, framesOut, hdAdjusts, dtAdjusts.
                               +dbleBytes*Nekx(1)*Nekx(2)  ![ Real field data with Fourier space dimensions.

#if ((useMPI > 0) && (FFTmode == 0))
  bpGsize2D(realAliased2Didx) = dbleBytes*7 &          ![ Doubles for simulationTime, dt, hDiffOrder, hDiffs.
                               +intBytes*4 &           ![ Integer for timeSteps, framesOut, hdAdjusts, dtAdjusts.
                               +dbleBytes*NyaG*Nxa(2)  ![ Real field data.
#else
  bpGsize2D(realAliased2Didx) = dbleBytes*7 &            ![ Doubles for simulationTime, dt, hDiffOrder, hDiffs.
                               +intBytes*4 &             ![ Integer for timeSteps, framesOut, hdAdjusts, dtAdjusts.
                               +dbleBytes*Nxa(2)*Nxa(1)  ![ Real field data.
#endif

  ![ .......... System for ADIOS file used for restarts ......... ]!
  useRestartFile0 = .true.
  ![ Group used to write 2D complex field.
  call adios_declare_group(bpGRestart,"restart","",ADIOS_STAT_DEFAULT,errFlagADIOS)
  call adios_select_method(bpGRestart,"MPI","","",errFlagADIOS)
  ![ Convert local and global dimensions of field data to string.
  write(cN1L,'(i0)') Nekx(1)
  write(cN2L,'(i0)') Nekx(2)
  write(cN1G,'(i0)') NekyG
  write(cN2G,'(i0)') NekxG
  write(cN1off,'(i0)') firstkG(1)-1
  write(cN2off,'(i0)') firstkG(2)-1
  cDimsL   = cN1L//","//cN2L
  cDimsG   = cN1G//","//cN2G
  cOffsets = cN1off//","//cN2off
  ![ Define the variables that we wish to include in the output file.
  ![ For each new variable defined (via adios_define_var) one must
  ![ include the appropriate contribution to bpGsize below.
  ![ Simulation time for this frame.
  call adios_define_var(bpGRestart,"simulationTime","",adios_double,"","","",bpVarID)
  ![ Time steps completed by this frame.
  call adios_define_var(bpGRestart,"timeSteps","",adios_integer,"","","",bpVarID)
  ![ (Absolute) number of frames outputted so far.
  call adios_define_var(bpGRestart,"framesOut","",adios_integer,"","","",bpVarID)
  ![ (Absolute) number of times hyperdiffusion has been adjusted.
  call adios_define_var(bpGRestart,"hdAdjusts","",adios_integer,"","","",bpVarID)
  ![ (Absolute) number of times dt has been adjusted.
  call adios_define_var(bpGRestart,"dtAdjusts","",adios_integer,"","","",bpVarID)
  ![ Time steps size.
  call adios_define_var(bpGRestart,"dt","",adios_double,"","","",bpVarID)
  ![ Hyper-diffusion order.
  call adios_define_var(bpGRestart,"hDiffOrder","",adios_double,"","","",bpVarID)
  ![ Hyper-diffusion coefficient.
  call adios_define_var(bpGRestart,"hDiffne","",adios_double,"","","",bpVarID)
  call adios_define_var(bpGRestart,"hDiffTe","",adios_double,"","","",bpVarID)
  call adios_define_var(bpGRestart,"hDiffni","",adios_double,"","","",bpVarID)
  call adios_define_var(bpGRestart,"hDiffTi","",adios_double,"","","",bpVarID)
  ![ Save the densities, temperatures and potential at the last time step.
  call adios_define_var(bpGRestart,"denElck","",adios_double_complex, &
                        trim(cDimsL),trim(cDimsG),trim(cOffsets),bpVarID)
  call adios_define_var(bpGRestart,"tempElck","",adios_double_complex, &
                        trim(cDimsL),trim(cDimsG),trim(cOffsets),bpVarID)
  call adios_define_var(bpGRestart,"denIonk","",adios_double_complex, &
                        trim(cDimsL),trim(cDimsG),trim(cOffsets),bpVarID)
  call adios_define_var(bpGRestart,"tempIonk","",adios_double_complex, &
                        trim(cDimsL),trim(cDimsG),trim(cOffsets),bpVarID)
  call adios_define_var(bpGRestart,"phik","",adios_double_complex, &
                        trim(cDimsL),trim(cDimsG),trim(cOffsets),bpVarID)
#if (STEPPER == 2)
  ![ For trapezoidal leapfrogp save the fields at the second to last time step.
  call adios_define_var(bpGRestart,"denElcks","",adios_double_complex, &
                        trim(cDimsL),trim(cDimsG),trim(cOffsets),bpVarID)
  call adios_define_var(bpGRestart,"tempElcks","",adios_double_complex, &
                        trim(cDimsL),trim(cDimsG),trim(cOffsets),bpVarID)
  call adios_define_var(bpGRestart,"denIonks","",adios_double_complex, &
                        trim(cDimsL),trim(cDimsG),trim(cOffsets),bpVarID)
  call adios_define_var(bpGRestart,"tempIonks","",adios_double_complex, &
                        trim(cDimsL),trim(cDimsG),trim(cOffsets),bpVarID)
  call adios_define_var(bpGRestart,"phiks","",adios_double_complex, &
                        trim(cDimsL),trim(cDimsG),trim(cOffsets),bpVarID)
#endif
  ![ Add other quantities defining this simulation as attributes in ADIOS file.
  call writeAttributes(bpGRestart,complex2Didx)

  call adios_define_schema_version(bpGRestart, bpSchemaVersion)

#if (STEPPER == 2)
  ![ Trapezoidal leapfrog has one extra complex field.
  bpGsizeRestart = dbleBytes*7 &                       ![ Doubles for simulationTime, dt, hDiffOrder, hDiffs.
                  +intBytes*4 &                        ![ Integer for timeSteps, framesOut, hdAdjusts, dtAdjusts.
                  +2*5*dbleCmplxBytes*Nekx(1)*Nekx(2)  ![ Complex field data.
#else
  bpGsizeRestart = dbleBytes*7 &                     ![ Doubles for simulationTime, dt, hDiffOrder, hDiffs.
                  +intBytes*4 &                      ![ Integer for timeSteps, framesOut, hdAdjusts, dtAdjusts.
                  +5*dbleCmplxBytes*Nekx(1)*Nekx(2)  ![ Complex field data.
#endif

  ![ ............. USER CAN CREATE MORE OUTPUT FILES HERE ........... ]!
  ![ Create a file for each field one wishes to save.

  call setup2DFieldComplex('phik')    
  call setup2DFieldComplex('denElck') 
  call setup2DFieldComplex('tempElck')
  call setup2DFieldComplex('denIonk') 
  call setup2DFieldComplex('tempIonk')

  call setup2DFieldFourierReal('Gamma0Ion')
  call setup2DFieldFourierReal('Gamma0Elc')
  call setup2DFieldFourierReal('avgJ0Ion')
  call setup2DFieldFourierReal('avgJ0Elc')
  call setup2DFieldFourierReal('poiSbIon')
  call setup2DFieldFourierReal('poiSbElc')

!  call setup2DFieldReal('denElc')
!  call setup2DFieldReal('tempElc')
!  call setup2DFieldReal('denIon')
!  call setup2DFieldReal('tempIon')

!  ![ These are used to test the FFTs (see testFFTs).
!  if (fftMember) call setup2DFieldRealAliased('phi')
!  if (fftMember) call setup2DFieldRealAliased('vort')

#else
  ![ Initialize binary output files.
  ![ To output other fields simply add other calls to:
  ![   a) open (and binAttributes) with the corresponding variable name here.
  ![   b) out2DField(Complex/Real) with the correct unit number and corresponding
  ![      variable to output in writeFields. 
  ![   c) close the file for each output variable in terminate_IO.
  ![ Use file unit numbers > 10 and < 99.

  call setup2DFieldComplex(11,'phik')    
  call setup2DFieldComplex(12,'denElck') 
  call setup2DFieldComplex(13,'tempElck')
  call setup2DFieldComplex(14,'denIonk') 
  call setup2DFieldComplex(15,'tempIonk')

!  call setup2DFieldReal(16,'phi')
!  call setup2DFieldReal(17,'denElc')
!  call setup2DFieldReal(18,'tempElc')
!  call setup2DFieldReal(19,'denIon')
!  call setup2DFieldReal(20,'tempIon')

#endif

  end subroutine setupIOfiles
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
#if (ioMode == 0)    
  subroutine readRestartInputs
  ![ When restarting from a previous simulation, we read the
  ![ inputs from the restart file.
  implicit none
  integer, parameter  :: rfIdx=99
  integer             :: intBuf0(2),intBuf1(6)
  double precision    :: doubleBuf0(17)
  logical             :: fileExists
  double precision    :: restartSimTime(2)

  do lc = 0,1
    ![ Full path to restart file.
    if (lc==0) restartFile = trim(restartDir)//restartFileName//'0'//sFileExt
    if (lc==1) restartFile = trim(restartDir)//restartFileName//'1'//sFileExt

    inquire(file=trim(restartFile), exist=fileExists)

    if (fileExists) then
      ![ Read simulation time of this restart file.
      open(unit=rfIdx,file=trim(restartFile),status='unknown',form='unformatted')
      read(unit=rfIdx) intBuf0, doubleBuf0, intBuf1   ![ Skip the attributes. 
      read(unit=rfIdx) restartSimTime(lc+1)
      close(unit=rfIdx)    ![ Close restart file.
    else
      restartSimTime(lc+1) = -9.99e11  ![ Arbitrary small number.
    endif
  enddo

  ![ Select full path to restart file with latest simulation time.
  if (restartSimTime(1) > restartSimTime(2)) then
    restartFile = trim(restartDir)//restartFileName//'0'//sFileExt
  else
    restartFile = trim(restartDir)//restartFileName//'1'//sFileExt
  endif

  ![ Open restart file and read attributes.
  ![ IMPORTANT: Restart file also has pre-processor variables as attributes
  ![ but currently we do not read those in and assume the user is careful.
  open(unit=rfIdx,file=trim(restartFile),status='unknown',form='unformatted')
  read(unit=rfIdx) intBuf0, doubleBuf0, intBuf1
  close(unit=rfIdx)    ![ Close restart file.
  ![ Inputs defining the grid have to be the same as those in the
  ![ restart file (not currently supporting grid changes upon restart).
  NkxG = intBuf0(1)
  NkyG = intBuf0(2)

  kxMin = doubleBuf0(1)
  kyMin = doubleBuf0(2)

  !( Other inputs are chosen from the restart file if favorRestart=T.
  if (favorRestart) then
    Lnorm         = doubleBuf0(3)
    omSte         = doubleBuf0(4)
    omde          = doubleBuf0(5)
    tau           = doubleBuf0(6)
    mu            = doubleBuf0(7)
    deltae        = doubleBuf0(8)
    deltaPerpe    = doubleBuf0(9)
    eta_e         = doubleBuf0(10)
    deltai        = doubleBuf0(11)
    deltaPerpi    = doubleBuf0(12)
    eta_i         = doubleBuf0(13)
    lambdaD       = doubleBuf0(14)
    tRateOutput   = doubleBuf0(15)
    tRateAdjustHD = doubleBuf0(16)
    tRateAdjustDt = doubleBuf0(17)
  endif

  end subroutine readRestartInputs
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine setup2DFieldComplex(fIdxIn,cVarName)
  ![ Create binary file used to store 2D complex field.
  implicit none
  integer, intent(in)          :: fIdxIn
  character(len=*), intent(in) :: cVarName

  if (isRestart) then
    open(unit=fIdxIn,file=trim(outputDir)//cVarName//sFileExt, &
         status='old',form='unformatted',position='append')
  else
    open(unit=fIdxIn,file=trim(outputDir)//cVarName//sFileExt, &
         status='unknown',form='unformatted')
    call writeAttributes(fIdxIn,complex2Didx)    ![ Write some key variables at the top.
  endif

  end subroutine setup2DFieldComplex
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine setup2DFieldReal(fIdxIn,cVarName)
  ![ Create binary file used to store 2D real field.
  implicit none
  integer, intent(in)          :: fIdxIn
  character(len=*), intent(in) :: cVarName

  if (isRestart) then
    open(unit=fIdxIn,file=trim(outputDir)//cVarName//sFileExt, &
         status='old',form='unformatted',position='append')
  else
    open(unit=fIdxIn,file=trim(outputDir)//cVarName//sFileExt, &
         status='unknown',form='unformatted')
    call writeAttributes(fIdxIn,real2Didx)    ![ Write some key variables at the top.
  endif

  end subroutine setup2DFieldReal
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine setup2DFieldFourierReal(fIdxIn,cVarName)
  ![ Create binary file used to store 2D real field.
  implicit none
  integer, intent(in)          :: fIdxIn
  character(len=*), intent(in) :: cVarName

  if (isRestart) then
    open(unit=fIdxIn,file=trim(outputDir)//cVarName//sFileExt, &
         status='old',form='unformatted',position='append')
  else
    open(unit=fIdxIn,file=trim(outputDir)//cVarName//sFileExt, &
         status='unknown',form='unformatted')
    call writeAttributes(fIdxIn,fourierReal2Didx)    ![ Write some key variables at the top.
  endif

  end subroutine setup2DFieldFourierReal
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine setup2DFieldRealAliased(fIdxIn,cVarName)
  ![ Create binary file used to store 2D real field defined in aliased real space.
  implicit none
  integer, intent(in)          :: fIdxIn
  character(len=*), intent(in) :: cVarName

  if (isRestart) then
    open(unit=fIdxIn,file=trim(outputDir)//cVarName//sFileExt, &
         status='old',form='unformatted',position='append')
  else
    open(unit=fIdxIn,file=trim(outputDir)//cVarName//sFileExt, &
         status='unknown',form='unformatted')
    call writeAttributes(fIdxIn,realAliased2Didx)    ![ Write some key variables at the top.
  endif

  end subroutine setup2DFieldRealAliased
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine writeAttributes(fileUnit,fileDataType)
  ![ Add quantities (attributes) to the top of a binary file.
  ![ These are used by post-processing to reconstruct details of the simulation.
  implicit none
  integer, intent(in) :: fileUnit,fileDataType
  integer             :: intBuf0(2),intBuf1(7)
  double precision    :: doubleBuf0(17)


  if ((fileDataType == complex2Didx) .or. (fileDataType == fourierReal2Didx)) then
    intBuf0    = [NkxG,NkyG]
    doubleBuf0 = [kxMin,kyMin,Lnorm,omSte,omde,tau,mu,deltae, &
                  deltaPerpe,eta_e,deltai,deltaPerpi,eta_i,lambdaD, &
                  tRateOutput,tRateAdjustHD,tRateAdjustDt]
  elseif (fileDataType == real2Didx) then
    intBuf0    = [NxG,NyG]
    doubleBuf0 = [dx,dy,Lnorm,omSte,omde,tau,mu,deltae, &
                  deltaPerpe,eta_e,deltai,deltaPerpi,eta_i,lambdaD, &
                  tRateOutput,tRateAdjustHD,tRateAdjustDt]
  elseif (fileDataType == realAliased2Didx) then
    intBuf0    = [NxaG,NyaG]
    doubleBuf0 = [dxa,dya,Lnorm,omSte,omde,tau,mu,deltae, &
                  deltaPerpe,eta_e,deltai,deltaPerpi,eta_i,lambdaD, &
                  tRateOutput,tRateAdjustHD,tRateAdjustDt]
  endif
  intBuf1    = [adiabaticElc,adiabaticIon,useMPI,ioMode,HDmodel,STEPPER,ADJUSTdt]
  write(unit=fileUnit) intBuf0, doubleBuf0, intBuf1
  end subroutine writeAttributes
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine out2DFieldComplex(fileUnit,kData)
  ![ Write complex field data to already opened binary output file.
  ![ We will also write the simulation time and time steps completed. 
  implicit none
  integer, intent(in)        :: fileUnit
  double complex, intent(in) :: kData(:,:)

  ![ Also save other variables changing in time useful for
  ![ restarts and analysis.
  write(unit=fileUnit) simTime, timeSteps, &
                       framesOut, hdAdjusts, dtAdjusts, &
                       dt, hDiffOrder, hDiffne, hDiffTe, hDiffni, hDiffTi, &
                       kData

  end subroutine out2DFieldComplex
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine out2DFieldReal(fileUnit,rData)
  ![ Write real field data to already opened binary output file.
  ![ We will also write the simulation time and time steps completed. 
  implicit none
  integer, intent(in)          :: fileUnit
  double precision, intent(in) :: rData(:,:)

  ![ Also save other variables changing in time useful for
  ![ restarts and analysis.
  write(unit=fileUnit) simTime, timeSteps, &
                       framesOut, hdAdjusts, dtAdjusts, &
                       dt, hDiffOrder, hDiffne, hDiffTe, hDiffni, hDiffTi, &
                       rData

  end subroutine out2DFieldReal
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine out2DFieldFourierReal(fileUnit,rData)
  ![ Write real field (defined inn Fourier space) data to already opened binary
  ![ output file. We will also write the simulation time and time steps completed. 
  implicit none
  integer, intent(in)          :: fileUnit
  double precision, intent(in) :: rData(:,:)

  ![ Also save other variables changing in time useful for
  ![ restarts and analysis.
  write(unit=fileUnit) simTime, timeSteps, &
                       framesOut, hdAdjusts, dtAdjusts, &
                       dt, hDiffOrder, hDiffne, hDiffTe, hDiffni, hDiffTi, &
                       rData

  end subroutine out2DFieldFourierReal
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine out2DFieldRealAliased(fileUnit,rData)
  ![ Write real field (in aliased real sapce) data to already opened binary output
  ![ file. We will also write the simulation time and time steps completed. 
  implicit none
  integer, intent(in)          :: fileUnit
  double precision, intent(in) :: rData(:,:)

  ![ Also save other variables changing in time useful for
  ![ restarts and analysis.
  write(unit=fileUnit) simTime, timeSteps, &
                       framesOut, hdAdjusts, dtAdjusts, &
                       dt, hDiffOrder, hDiffne, hDiffTe, hDiffni, hDiffTi, &
                       rData

  end subroutine out2DFieldRealAliased
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine outRestart
  ![ Write the restart file with latest data.
  implicit none
  character(len=100)           :: filePathName

  if (useRestartFile0) then
    filePathName = trim(outputDir)//restartFileName//'0'//sFileExt
  else
    filePathName = trim(outputDir)//restartFileName//'1'//sFileExt
  endif
  useRestartFile0 = .not. useRestartFile0

  ![ Open restart file.
  open(unit=restartFileIdx,file=trim(filePathName), &
       status='unknown',form='unformatted')
  
  call writeAttributes(restartFileIdx,complex2Didx)    ![ Write some key variables at the top.
  
  ![ Also save other variables changing in time useful for
  ![ restarts and analysis.
#if (STEPPER == 2)
  ![ Trapezoidal leapfrog requires one extra field.
  write(unit=restartFileIdx) simTime, timeSteps, 
                             framesOut, hdAdjusts, dtAdjusts, &
                             dt, hDiffOrder, hDiffne, hDiffTe, hDiffni, hDiffTi, & 
                             phik, denElck, tempElck, denIonk, tempIonk, &
                             phiks, denElcks, tempElcks, denIonks, tempIonks
#else
  write(unit=restartFileIdx) simTime, timeSteps, &
                             framesOut, hdAdjusts, dtAdjusts, &
                             dt, hDiffOrder, hDiffne, hDiffTe, hDiffni, hDiffTi, &
                             phik, denElck, tempElck, denIonk, tempIonk
#endif

  close(unit=restartFileIdx)    ![ Close restart file.

  end subroutine outRestart
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine readRestartFields
  ![ When restarting from a previous simulation, we read the
  ![ dynamic complex fields from the restart file.
  implicit none
  integer, parameter  :: rfIdx=99
  integer             :: intBuf0(2),intBuf1(7)
  double precision    :: doubleBuf0(17)

  ![ Open restart file.
  open(unit=rfIdx,file=trim(restartFile),status='unknown',form='unformatted')

  read(unit=rfIdx) intBuf0, doubleBuf0, intBuf1   ![ Skip the attributes. 

#if (STEPPER == 2)
  ![ Trapezoidal leapfrog requires one extra field.
  read(unit=rfIdx) simTime, timeSteps, &
                   framesOut, hdAdjusts, dtAdjusts, &
                   dt, hDiffOrder, hDiffne, hDiffTe, hDiffni, hDiffTi, &
                   phik, denElck, tempElck, denIonk, tempIonk, &
                   phiks, denElcks, tempElcks, denIonks, tempIonks
#else
  read(unit=rfIdx) simTime, timeSteps, &
                   framesOut, hdAdjusts, dtAdjusts, &
                   dt, hDiffOrder, hDiffne, hDiffTe, hDiffni, hDiffTi, &
                   phik, denElck, tempElck, denIonk, tempIonk
#endif

  close(unit=rfIdx)    ![ Close restart file.

  end subroutine readRestartFields

#else

  subroutine readRestartInputs
  ![ When restarting from a previous simulation, we read the
  ![ inputs from the restart file.
  implicit none
  logical          :: fileExists
  double precision :: restartSimTime(2)
  integer*8        :: bpBlockSel

  do lc = 0,1
    ![ Full path to restart file.
    if (lc==0) restartFile = trim(restartDir)//restartFileName//'0'//sFileExt
    if (lc==1) restartFile = trim(restartDir)//restartFileName//'1'//sFileExt

    inquire(file=trim(restartFile), exist=fileExists)

    if (fileExists) then
      ![ Read simulation time of this restart file.
      call adios_read_open_file(bpHandle,trim(restartFile),bpReadMethod, &
                                MPI_COMM_WORLD,errFlagADIOS)
      bpBlockSel = 0
      call adios_schedule_read(bpHandle,bpBlockSel,"simulationTime",0,1,restartSimTime(lc+1),errFlagADIOS)
      call adios_perform_reads(bpHandle,errFlagADIOS)  ![ Perform scheduled reads.
      call adios_read_close(bpHandle,errFlagADIOS)     ![ Close restart file.
    else
      restartSimTime(lc+1) = -9.99e11  ![ Arbitrary small number.
    endif
  enddo

  ![ Select full path to restart file with latest simulation time.
  if (restartSimTime(1) > restartSimTime(2)) then
    restartFile = trim(restartDir)//restartFileName//'0'//sFileExt
  else
    restartFile = trim(restartDir)//restartFileName//'1'//sFileExt
  endif

  if (myRank == 0) then
    print*,'Reading inputs from ',trim(restartFile)
    write(*,*) ' '
  endif

  ![ Open restart file.
  call adios_read_open_file(bpHandle,trim(restartFile),bpReadMethod, &
                            MPI_COMM_WORLD,errFlagADIOS)
  
  ![ Inputs defining the grid have to be the same as those in the
  ![ restart file (not currently supporting grid changes upon restart).
  call adios_get_attr(bpHandle,  "Nkx", NkxG,errFlagADIOS)
  call adios_get_attr(bpHandle,  "Nky", NkyG,errFlagADIOS)
  call adios_get_attr(bpHandle,"kxMin",kxMin,errFlagADIOS)
  call adios_get_attr(bpHandle,"kyMin",kyMin,errFlagADIOS)

  ![ Other inputs are chosen from the restart file if favorRestart=T.
  ![ IMPORTANT: Restart file also has pre-processor variables as attributes
  ![ but currently we do not read those in and assume the user is careful.
  if (favorRestart) then
    call adios_get_attr(bpHandle,        "Lnorm",        Lnorm,errFlagADIOS)
    call adios_get_attr(bpHandle,        "omSte",        omSte,errFlagADIOS)
    call adios_get_attr(bpHandle,         "omde",         omde,errFlagADIOS)
    call adios_get_attr(bpHandle,          "tau",          tau,errFlagADIOS)
    call adios_get_attr(bpHandle,           "mu",           mu,errFlagADIOS)
    call adios_get_attr(bpHandle,       "deltae",       deltae,errFlagADIOS)
    call adios_get_attr(bpHandle,   "deltaPerpe",   deltaPerpe,errFlagADIOS)
    call adios_get_attr(bpHandle,        "eta_e",        eta_e,errFlagADIOS)
    call adios_get_attr(bpHandle,       "deltai",       deltai,errFlagADIOS)
    call adios_get_attr(bpHandle,   "deltaPerpi",   deltaPerpi,errFlagADIOS)
    call adios_get_attr(bpHandle,        "eta_i",        eta_i,errFlagADIOS)
    call adios_get_attr(bpHandle,      "lambdaD",      lambdaD,errFlagADIOS)
    call adios_get_attr(bpHandle,  "tRateOutput",  tRateOutput,errFlagADIOS)
    call adios_get_attr(bpHandle,"tRateAdjustHD",tRateAdjustHD,errFlagADIOS)
    call adios_get_attr(bpHandle,"tRateAdjustDt",tRateAdjustDt,errFlagADIOS)
  endif

  ![ Close restart file.
  call adios_read_close(bpHandle,errFlagADIOS)

  end subroutine readRestartInputs
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine setup2DFieldComplex(cVarName)
  ![ Create ADIOS file used to store 2D complex field.
  implicit none
  character(len=*), intent(in) :: cVarName
  character(len=100)           :: filePathName

  ![ Open output file.
  filePathName = trim(outputDir)//cVarName//sFileExt
  if (isRestart) then
    call adios_open(bpHandle,"complex2d",trim(filePathName),'a',MPI_COMM_WORLD,errFlagADIOS)
  else
    call adios_open(bpHandle,"complex2d",trim(filePathName),'w',MPI_COMM_WORLD,errFlagADIOS)
  endif
  ![ Pass group size to the internal ADIOS transport structure.
  call adios_group_size(bpHandle,bpGsize2D(complex2Didx),bpGtotalSize2D(complex2Didx),errFlagADIOS)
  ![ Close ADIOS file.
  call adios_close(bpHandle,errFlagADIOS)

  end subroutine setup2DFieldComplex
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine setup2DFieldReal(cVarName)
  ![ Create ADIOS file used to store 2D real field.
  implicit none
  character(len=*), intent(in) :: cVarName
  character(len=100)           :: filePathName

  ![ Open output file.
  filePathName = trim(outputDir)//cVarName//sFileExt
  if (isRestart) then
    call adios_open(bpHandle,"real2d",trim(filePathName),'a',MPI_COMM_WORLD,errFlagADIOS)
  else
    call adios_open(bpHandle,"real2d",trim(filePathName),'w',MPI_COMM_WORLD,errFlagADIOS)
  endif
  ![ Pass group size to the internal ADIOS transport structure.
  call adios_group_size(bpHandle,bpGsize2D(real2Didx),bpGtotalSize2D(real2Didx),errFlagADIOS)
  ![ Close ADIOS file.
  call adios_close(bpHandle,errFlagADIOS)

  end subroutine setup2DFieldReal
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine setup2DFieldFourierReal(cVarName)
  ![ Create ADIOS file used to store 2D real field with the dimensions
  ![ of a Fourier space field.
  implicit none
  character(len=*), intent(in) :: cVarName
  character(len=100)           :: filePathName

  ![ Open output file.
  filePathName = trim(outputDir)//cVarName//sFileExt
  if (isRestart) then
    call adios_open(bpHandle,"fourierReal2d",trim(filePathName),'a',MPI_COMM_WORLD,errFlagADIOS)
  else
    call adios_open(bpHandle,"fourierReal2d",trim(filePathName),'w',MPI_COMM_WORLD,errFlagADIOS)
  endif
  ![ Pass group size to the internal ADIOS transport structure.
  call adios_group_size(bpHandle,bpGsize2D(fourierReal2Didx),bpGtotalSize2D(fourierReal2Didx),errFlagADIOS)
  ![ Close ADIOS file.
  call adios_close(bpHandle,errFlagADIOS)

  end subroutine setup2DFieldFourierReal
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine setup2DFieldRealAliased(cVarName)
  ![ Create ADIOS file used to store 2D real field (defined in aliased real space).
  implicit none
  character(len=*), intent(in) :: cVarName
  character(len=100)           :: filePathName

  ![ Open output file.
  filePathName = trim(outputDir)//cVarName//sFileExt
#if ((useMPI > 0) && (FFTmode == 0))
  if (isRestart) then
    call adios_open(bpHandle,"realAliased2d",trim(filePathName),'a',fftCOMM,errFlagADIOS)
  else
    call adios_open(bpHandle,"realAliased2d",trim(filePathName),'w',fftCOMM,errFlagADIOS)
  endif
#else
  if (isRestart) then
    call adios_open(bpHandle,"realAliased2d",trim(filePathName),'a',MPI_COMM_WORLD,errFlagADIOS)
  else
    call adios_open(bpHandle,"realAliased2d",trim(filePathName),'w',MPI_COMM_WORLD,errFlagADIOS)
  endif
#endif
  ![ Pass group size to the internal ADIOS transport structure.
  call adios_group_size(bpHandle,bpGsize2D(realAliased2Didx),bpGtotalSize2D(realAliased2Didx),errFlagADIOS)
  ![ Close ADIOS file.
  call adios_close(bpHandle,errFlagADIOS)

  end subroutine setup2DFieldRealAliased
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine writeAttributes(bpG,fileDataType)
  ![ Add quantities (attributes) to the ADIOS file.
  ![ These are used by post-processing to reconstruct details of the simulation.
  implicit none
  integer*8, intent(in) :: bpG
  integer, intent(in)   :: fileDataType

  if ((fileDataType == complex2Didx) .or. (fileDataType == fourierReal2Didx)) then
    call adios_define_attribute_byvalue(bpG,        "Nkx","",1,         NkxG,errFlagADIOS)
    call adios_define_attribute_byvalue(bpG,        "Nky","",1,         NkyG,errFlagADIOS)
    call adios_define_attribute_byvalue(bpG,      "kxMin","",1,        kxMin,errFlagADIOS)
    call adios_define_attribute_byvalue(bpG,      "kyMin","",1,        kyMin,errFlagADIOS)
  elseif (fileDataType == real2Didx) then
    call adios_define_attribute_byvalue(bpG,         "Nx","",1,          NxG,errFlagADIOS)
    call adios_define_attribute_byvalue(bpG,         "Ny","",1,          NyG,errFlagADIOS)
    call adios_define_attribute_byvalue(bpG,         "dx","",1,           dx,errFlagADIOS)
    call adios_define_attribute_byvalue(bpG,         "dy","",1,           dy,errFlagADIOS)
  elseif (fileDataType == realAliased2Didx) then
    call adios_define_attribute_byvalue(bpG,        "Nxa","",1,         NxaG,errFlagADIOS)
    call adios_define_attribute_byvalue(bpG,        "Nya","",1,         NyaG,errFlagADIOS)
    call adios_define_attribute_byvalue(bpG,        "dxa","",1,          dxa,errFlagADIOS)
    call adios_define_attribute_byvalue(bpG,        "dya","",1,          dya,errFlagADIOS)
  endif
  call adios_define_attribute_byvalue(bpG,        "Lnorm","",1,        Lnorm,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,        "omSte","",1,        omSte,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,         "omde","",1,         omde,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,          "tau","",1,          tau,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,           "mu","",1,           mu,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,       "deltae","",1,       deltae,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,   "deltaPerpe","",1,   deltaPerpe,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,        "eta_e","",1,        eta_e,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,       "deltai","",1,       deltai,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,   "deltaPerpi","",1,   deltaPerpi,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,        "eta_i","",1,        eta_i,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,      "lambdaD","",1,      lambdaD,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,  "tRateOutput","",1,  tRateOutput,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,"tRateAdjustHD","",1,tRateAdjustHD,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,"tRateAdjustDt","",1,tRateAdjustDt,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG, "adiabaticElc","",1, adiabaticElc,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG, "adiabaticIon","",1, adiabaticIon,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,       "useMPI","",1,       useMPI,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,       "ioMode","",1,       ioMode,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,      "HDmodel","",1,      HDmodel,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,      "STEPPER","",1,      STEPPER,errFlagADIOS)
  call adios_define_attribute_byvalue(bpG,     "ADJUSTdt","",1,     ADJUSTdt,errFlagADIOS)
  end subroutine writeAttributes
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine out2DFieldComplex(cVarName,kData)
  ![ Open ADIOS file and write complex field data to it. 
  implicit none
  character(len=*), intent(in) :: cVarName
  double complex, intent(in)   :: kData(:,:)

  ![ Open output file.
  call adios_open(bpHandle,"complex2d",trim(outputDir)//cVarName//sFileExt, &
                  'a',MPI_COMM_WORLD,errFlagADIOS)
  ![ Pass group size to the internal ADIOS transport structure.
  call adios_group_size(bpHandle,bpGsize2D(complex2Didx),bpGtotalSize2D(complex2Didx),errFlagADIOS)
  ![ Write other variables changing in time we wish to save.
  call adios_write(bpHandle,"simulationTime",   simTime,errFlagADIOS)
  call adios_write(bpHandle,     "timeSteps", timeSteps,errFlagADIOS)
  call adios_write(bpHandle,     "framesOut", framesOut,errFlagADIOS)
  call adios_write(bpHandle,     "hdAdjusts", hdAdjusts,errFlagADIOS)
  call adios_write(bpHandle,     "dtAdjusts", dtAdjusts,errFlagADIOS)
  call adios_write(bpHandle,            "dt",        dt,errFlagADIOS)
  call adios_write(bpHandle,    "hDiffOrder",hDiffOrder,errFlagADIOS)
  call adios_write(bpHandle,       "hDiffne",   hDiffne,errFlagADIOS)
  call adios_write(bpHandle,       "hDiffTe",   hDiffTe,errFlagADIOS)
  call adios_write(bpHandle,       "hDiffni",   hDiffni,errFlagADIOS)
  call adios_write(bpHandle,       "hDiffTi",   hDiffTi,errFlagADIOS)
  ![ Write 2D complex field variable.
  call adios_write(bpHandle,  "complexField",     kData,errFlagADIOS)
  ![ Close ADIOS file.
  call adios_close(bpHandle,errFlagADIOS)

  end subroutine out2DFieldComplex
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine out2DFieldReal(cVarName,rData)
  ![ Open ADIOS file and write real field data to it. 
  implicit none
  character(len=*), intent(in) :: cVarName
  double precision, intent(in) :: rData(:,:)

  ![ Open output file.
  call adios_open(bpHandle,"real2d",trim(outputDir)//cVarName//sFileExt, &
                  'a',MPI_COMM_WORLD,errFlagADIOS)
  ![ Pass group size to the internal ADIOS transport structure.
  call adios_group_size(bpHandle,bpGsize2D(real2Didx),bpGtotalSize2D(real2Didx),errFlagADIOS)
  ![ Write other variables changing in time we wish to save.
  call adios_write(bpHandle,"simulationTime",   simTime,errFlagADIOS)
  call adios_write(bpHandle,     "timeSteps", timeSteps,errFlagADIOS)
  call adios_write(bpHandle,     "framesOut", framesOut,errFlagADIOS)
  call adios_write(bpHandle,     "hdAdjusts", hdAdjusts,errFlagADIOS)
  call adios_write(bpHandle,     "dtAdjusts", dtAdjusts,errFlagADIOS)
  call adios_write(bpHandle,            "dt",        dt,errFlagADIOS)
  call adios_write(bpHandle,    "hDiffOrder",hDiffOrder,errFlagADIOS)
  call adios_write(bpHandle,       "hDiffne",   hDiffne,errFlagADIOS)
  call adios_write(bpHandle,       "hDiffTe",   hDiffTe,errFlagADIOS)
  call adios_write(bpHandle,       "hDiffni",   hDiffni,errFlagADIOS)
  call adios_write(bpHandle,       "hDiffTi",   hDiffTi,errFlagADIOS)
  ![ Write 2D real field variable.
  call adios_write(bpHandle,     "realField",     rData,errFlagADIOS)
  ![ Close ADIOS file.
  call adios_close(bpHandle,errFlagADIOS)

  end subroutine out2DFieldReal
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine out2DFieldFourierReal(cVarName,rData)
  ![ Open ADIOS file and write real field (defined in Fourier space) data to it. 
  ![ This means that it is a real data array, but has the dimensions of a
  ![ complex field in Fourier space.
  implicit none
  character(len=*), intent(in) :: cVarName
  double precision, intent(in) :: rData(:,:)

  ![ Open output file.
  call adios_open(bpHandle,"fourierReal2d",trim(outputDir)//cVarName//sFileExt, &
                  'a',MPI_COMM_WORLD,errFlagADIOS)
  ![ Pass group size to the internal ADIOS transport structure.
  call adios_group_size(bpHandle,bpGsize2D(fourierReal2Didx),bpGtotalSize2D(fourierReal2Didx),errFlagADIOS)
  ![ Write other variables changing in time we wish to save.
  call adios_write(bpHandle,  "simulationTime",   simTime,errFlagADIOS)
  call adios_write(bpHandle,       "timeSteps", timeSteps,errFlagADIOS)
  call adios_write(bpHandle,       "framesOut", framesOut,errFlagADIOS)
  call adios_write(bpHandle,       "hdAdjusts", hdAdjusts,errFlagADIOS)
  call adios_write(bpHandle,       "dtAdjusts", dtAdjusts,errFlagADIOS)
  call adios_write(bpHandle,              "dt",        dt,errFlagADIOS)
  call adios_write(bpHandle,      "hDiffOrder",hDiffOrder,errFlagADIOS)
  call adios_write(bpHandle,         "hDiffne",   hDiffne,errFlagADIOS)
  call adios_write(bpHandle,         "hDiffTe",   hDiffTe,errFlagADIOS)
  call adios_write(bpHandle,         "hDiffni",   hDiffni,errFlagADIOS)
  call adios_write(bpHandle,         "hDiffTi",   hDiffTi,errFlagADIOS)
  ![ Write 2D real field variable defined in Fourier space.
  call adios_write(bpHandle,"fourierRealField",     rData,errFlagADIOS)
  ![ Close ADIOS file.
  call adios_close(bpHandle,errFlagADIOS)

  end subroutine out2DFieldFourierReal
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine out2DFieldRealAliased(cVarName,rData)
  ![ Open ADIOS file and write real field (defined in aliased real space) data to it. 
  implicit none
  character(len=*), intent(in) :: cVarName
  double precision, intent(in) :: rData(:,:)

  ![ Open output file.
#if ((useMPI > 0) && (FFTmode == 0))
  call adios_open(bpHandle,"realAliased2d",trim(outputDir)//cVarName//sFileExt, &
                  'a',fftCOMM,errFlagADIOS)
#else
  call adios_open(bpHandle,"realAliased2d",trim(outputDir)//cVarName//sFileExt, &
                  'a',MPI_COMM_WORLD,errFlagADIOS)
#endif
  ![ Pass group size to the internal ADIOS transport structure.
  call adios_group_size(bpHandle,bpGsize2D(realAliased2Didx),bpGtotalSize2D(realAliased2Didx),errFlagADIOS)
  ![ Write other variables changing in time we wish to save.
  call adios_write(bpHandle,  "simulationTime",   simTime,errFlagADIOS)
  call adios_write(bpHandle,       "timeSteps", timeSteps,errFlagADIOS)
  call adios_write(bpHandle,       "framesOut", framesOut,errFlagADIOS)
  call adios_write(bpHandle,       "hdAdjusts", hdAdjusts,errFlagADIOS)
  call adios_write(bpHandle,       "dtAdjusts", dtAdjusts,errFlagADIOS)
  call adios_write(bpHandle,              "dt",        dt,errFlagADIOS)
  call adios_write(bpHandle,      "hDiffOrder",hDiffOrder,errFlagADIOS)
  call adios_write(bpHandle,         "hDiffne",   hDiffne,errFlagADIOS)
  call adios_write(bpHandle,         "hDiffTe",   hDiffTe,errFlagADIOS)
  call adios_write(bpHandle,         "hDiffni",   hDiffni,errFlagADIOS)
  call adios_write(bpHandle,         "hDiffTi",   hDiffTi,errFlagADIOS)
  ![ Write 2D real field variable.
  call adios_write(bpHandle,"realAliasedField",     rData,errFlagADIOS)
  ![ Close ADIOS file.
  call adios_close(bpHandle,errFlagADIOS)

  end subroutine out2DFieldRealAliased
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine outRestart
  ![ Open and write ADIOS file used for restarting the simulation.
  implicit none
  character(len=100)           :: filePathName

  if (useRestartFile0) then
    filePathName = trim(outputDir)//restartFileName//'0'//sFileExt
  else
    filePathName = trim(outputDir)//restartFileName//'1'//sFileExt
  endif
  useRestartFile0 = .not. useRestartFile0

  ![ Open output file.
  call adios_open(bpHandle,"restart",trim(filePathName), &
                  'w',MPI_COMM_WORLD,errFlagADIOS)
  ![ Pass group size to the internal ADIOS transport structure.
  call adios_group_size(bpHandle,bpGsizeRestart,bpGtotalSizeRestart,errFlagADIOS)
  call writeAttributes(bpGRestart,complex2Didx)    ![ Write attributes into restart file.
  ![ Write other variables changing in time we wish to save.
  call adios_write(bpHandle,"simulationTime",   simTime,errFlagADIOS)
  call adios_write(bpHandle,     "timeSteps", timeSteps,errFlagADIOS)
  call adios_write(bpHandle,     "framesOut", framesOut,errFlagADIOS)
  call adios_write(bpHandle,     "hdAdjusts", hdAdjusts,errFlagADIOS)
  call adios_write(bpHandle,     "dtAdjusts", dtAdjusts,errFlagADIOS)
  call adios_write(bpHandle,            "dt",        dt,errFlagADIOS)
  call adios_write(bpHandle,    "hDiffOrder",hDiffOrder,errFlagADIOS)
  call adios_write(bpHandle,       "hDiffne",   hDiffne,errFlagADIOS)
  call adios_write(bpHandle,       "hDiffTe",   hDiffTe,errFlagADIOS)
  call adios_write(bpHandle,       "hDiffni",   hDiffni,errFlagADIOS)
  call adios_write(bpHandle,       "hDiffTi",   hDiffTi,errFlagADIOS)
  ![ Write the potential at the last time step.
  call adios_write(bpHandle,    "phik",    phik,errFlagADIOS)
  call adios_write(bpHandle, "denElck", denElck,errFlagADIOS)
  call adios_write(bpHandle,"tempElck",tempElck,errFlagADIOS)
  call adios_write(bpHandle, "denIonk", denIonk,errFlagADIOS)
  call adios_write(bpHandle,"tempIonk",tempIonk,errFlagADIOS)
#if (STEPPER == 2)
  call adios_write(bpHandle,    "phiks",    phiks,errFlagADIOS)
  call adios_write(bpHandle, "denElcks", denElcks,errFlagADIOS)
  call adios_write(bpHandle,"tempElcks",tempElcks,errFlagADIOS)
  call adios_write(bpHandle, "denIonks", denIonks,errFlagADIOS)
  call adios_write(bpHandle,"tempIonks",tempIonks,errFlagADIOS)
#endif
  ![ Close ADIOS file.
  call adios_close(bpHandle,errFlagADIOS)

  end subroutine outRestart
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine readRestartFields
  ![ When restarting from a previous simulation, we read the
  ![ inputs from the restart file.
  implicit none
  integer*8 :: bpBlockSel,selStart(2),selCount(2)

  ![ Open restart file.
  call adios_read_open_file(bpHandle,trim(restartFile),bpReadMethod, &
                            MPI_COMM_WORLD,errFlagADIOS)

  ![ Schedule reading of restart data variables.
  bpBlockSel = 0
  call adios_schedule_read(bpHandle,bpBlockSel,"simulationTime",0,1,   simTime,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,     "timeSteps",0,1, timeSteps,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,     "framesOut",0,1, framesOut,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,     "hdAdjusts",0,1, hdAdjusts,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,     "dtAdjusts",0,1, dtAdjusts,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,            "dt",0,1,        dt,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,    "hDiffOrder",0,1,hDiffOrder,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,       "hDiffne",0,1,   hDiffne,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,       "hDiffTe",0,1,   hDiffTe,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,       "hDiffni",0,1,   hDiffni,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,       "hDiffTi",0,1,   hDiffTi,errFlagADIOS)

  ![ Select this process' data to read.
  selStart(1) = firstkG(1)-1
  selStart(2) = firstkG(2)-1
  selCount(1) = Nekx(1)
  selCount(2) = Nekx(2)
  call adios_selection_boundingbox(bpBlockSel,2,selStart,selCount)
  
  call adios_schedule_read(bpHandle,bpBlockSel,          "phik",0,1,      phik,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,       "denElck",0,1,   denElck,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,      "tempElck",0,1,  tempElck,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,       "denIonk",0,1,   denIonk,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,      "tempIonk",0,1,  tempIonk,errFlagADIOS)
#if (STEPPER == 2)
  call adios_schedule_read(bpHandle,bpBlockSel,         "phiks",0,1,     phiks,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,      "denElcks",0,1,  denElcks,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,     "tempElcks",0,1, tempElcks,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,      "denIonks",0,1,  denIonks,errFlagADIOS)
  call adios_schedule_read(bpHandle,bpBlockSel,     "tempIonks",0,1, tempIonks,errFlagADIOS)
#endif

  call adios_perform_reads(bpHandle,errFlagADIOS)  ![ Perform scheduled reads.

  call adios_read_close(bpHandle,errFlagADIOS)     ![ Close restart file.

  end subroutine readRestartFields
#endif
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine writeFields
  ![ Write field variables to files.
  implicit none

  ![ To output other fields simply add more calls to
  ![ out2DField with the corresponding variable you wish to write.
  ![ Each field must have its own file created with calls to
  ![ setup2D in init_IO.
#if (ioMode == 0)
  ![ Write to already opened binary output files.
  ![ File unit number has to match that in setupIOfiles and terminate_IO.

  call out2DField(11,phik)
  call out2DField(12,denElck)
  call out2DField(13,tempElck)
  call out2DField(14,denIonk)
  call out2DField(15,tempIonk)

!  call out2DField(16,phi)
!  call out2DField(17,denElc)
!  call out2DField(18,tempElc)
!  call out2DField(19,denIon)
!  call out2DField(20,tempIon)

#else
  ![ Open ADIOS output files, write data to them and close them.

  call out2DField('phik',phik)
  call out2DField('denElck',denElck)
  call out2DField('tempElck',tempElck)
  call out2DField('denIonk',denIonk)
  call out2DField('tempIonk',tempIonk)

  if (timeSteps<1) then
    call out2DFieldFourierReal('Gamma0Ion',Gamma0Ion)
    call out2DFieldFourierReal('avgJ0Ion',avgJ0Ion)
    call out2DFieldFourierReal('poiSbIon',poiSbIon)
    call out2DFieldFourierReal('Gamma0Elc',Gamma0Elc)
    call out2DFieldFourierReal('avgJ0Elc',avgJ0Elc)
    call out2DFieldFourierReal('poiSbElc',poiSbElc)
  endif

!  call out2DField('phi',phi)
!  call out2DField('denElc',denElc)
!  call out2DField('tempElc',tempElc)
!  call out2DField('denIon',denIon)
!  call out2DField('tempIon',tempIon)

#endif

  end subroutine writeFields
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine terminate_IO
  ![ Close files (if needed) and output interface.
  ![ Free memory allocated for I/O.
  implicit none

#if (ioMode == 0)
  ![ Close binary output files.
  ![ File unit number has to match that in setupIOfiles and writeFields.

  close(unit=11)
  close(unit=12)
  close(unit=13)
  close(unit=14)
  close(unit=15)

!  close(unit=16)
!  close(unit=17)
!  close(unit=18)
!  close(unit=19)
!  close(unit=20)

#else
  ![ Wait for all processes to finish IO operations.
  call MPI_BARRIER(MPI_COMM_WORLD, errFlagMPI)
  ![ Terminate ADIOS read interface.
  call adios_read_finalize_method(bpReadMethod,errFlagADIOS)
  ![ Terminate ADIOS interface.
  call adios_finalize(myRank,errFlagADIOS)
#endif

  end subroutine terminate_IO
end module IOtools

![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

module FFTmodule
  use parameters
#if (useMPI > 0)
  use MPItools
#endif
  use utilities
  ![ Must include this module in any subroutine that uses FFTW routines.
  use, intrinsic :: iso_c_binding
  ![ The following include statement points to the directory where FFTW was installed.
#if (useMPI > 0)
  include 'fftw3-mpi.f03'
#else
  include 'fftw3.f03'
#endif

  ![ Buffer for storing aliased complex data (before backward xy FFT).
  type(C_PTR)                        :: funka_ptr
  complex(C_DOUBLE_COMPLEX), pointer :: funka(:,:)
  ![ Buffer for storing de-aliased complex data (after forward xy FFT).
  type(C_PTR)                        :: funk_ptr
  complex(C_DOUBLE_COMPLEX), pointer :: funk(:,:)
  ![ Pointers to arrays holding the real data (after backward xy FFT in the
  ![ case of aliased, and before forward xy FFT in the case of de-aliased).
  type(C_PTR)                        :: funRa_ptr
  type(C_PTR)                        :: funR_ptr

  logical, parameter :: outputWisdom = .true.
  ![ Threads used in each FFT.
  integer, parameter :: FFTWthreads=1

#if (useMPI > 0)

#if (FFTmode == 0)
![ else for if (FFTmode==0) then else.

  ![ Variables containing FFTW plans for 2D FFT 1D domain decomposed.
  type(C_PTR) :: plan2Da_c2r, plan2Da_r2c    ![ Plans for 2D transform of aliased data.
  type(C_PTR) :: plan2D_c2r, plan2D_r2c      ![ Plans for 2D transform of de-aliased data.

  ![ Buffer which stores the real aliased data (after backward xy FFT).
  real(C_DOUBLE), pointer :: funRa(:,:)
  ![ Buffer which stores the real de-aliased data (after before forward xy FFT).
  real(C_DOUBLE), pointer :: funR(:,:)

  ![ Number of aliased kx's in this rank.
  integer(C_INTPTR_T) :: lockxa
  ![ Index of the last x-element contained in the previous process along x.
  integer(C_INTPTR_T) :: lockxa0 

  ![ Array to contain all the de-aliased ky's, but distributed along kx.
  double complex, allocatable ::  funkyG(:,:)

  ![ User-defined type used to receive data from other processes along ky.
  type kyGatherInfoType
    integer, allocatable :: rankID(:), kyIdxLower(:), kyIdxUpper(:)
    integer, allocatable :: mpiDDtype(:)
    integer, allocatable :: req(:),stat(:,:)
  end type kyGatherInfoType

  ![ Each entry in kyGatherType corresponds to a process this rank receives
  ![ data from, when gathering along ky. This is done prior to the 2D FFT.
  type(kyGatherInfoType) :: kyGatherType

#else
![ else for if (FFTmode==0) then else.

  ![ Variables containing FFTW plans. For 2D domain decomposition
  ![ we have to use 1D complex to complex FFTs.
  type(C_PTR) :: plan1Da_c2r_x, plan1Da_r2c_y    ![ Plans for c2r transform of aliased data.
  type(C_PTR) :: plan1Da_r2c_x, plan1Da_c2r_y    ![ Plans for r2c transform of aliased data.
  type(C_PTR) :: plan1D_c2r_x, plan1D_r2c_y      ![ Plans for c2r transform of de-aliased data.
  type(C_PTR) :: plan1D_r2c_x, plan1D_c2r_y      ![ Plans for r2c transform of de-aliased data.

  ![ Buffer which stores the real aliased data (after backward xy FFT).
  complex(C_DOUBLE_COMPLEX), pointer :: funRa(:,:)
  ![ Buffer which stores the real de-aliased data (after before forward xy FFT).
  complex(C_DOUBLE_COMPLEX), pointer :: funR(:,:)
  ![ Intermediate buffers holding the aliased data that is complex along y
  ![ and its transpose (the latter with the conjugate pairs along ky).
  type(C_PTR)                        :: funkya_ptr, funkyaTc_ptr
  complex(C_DOUBLE_COMPLEX), pointer :: funkya(:,:), funkyaTc(:,:)
  double complex, allocatable        :: funkyaT(:,:)  ![ Local transpose of funkya.
  ![ Intermediate buffers holding the de-aliased data that is complex along y
  ![ and its transpose (the latter with the conjugate pairs along ky).
  type(C_PTR)                        :: funky_ptr, funkyTc_ptr
  complex(C_DOUBLE_COMPLEX), pointer :: funky(:,:), funkyTc(:,:)
  double complex, allocatable        :: funkyT(:,:)  ![ Local transpose of funky.

  ![ Number of aliased real space elements in this rank.
  integer(C_INTPTR_T) :: locxa
  ![ Number of aliased ky's and kx's in this rank along y and x, respectively.
  integer(C_INTPTR_T) :: lockya,lockxa
  ![ Index of the last y/x-element contained in the previous process along y/x.
  integer(C_INTPTR_T) :: lockya0,lockxa0 
  ![ Analogs for de-aliased FFTs.
  integer(C_INTPTR_T) :: lockx0,locky0

  ![ User-defined type used to store information used when appending
  ![ complex conjugates along ky.
  type cpInfoType
    integer, allocatable :: rankID(:), kyIdxLower(:), kyIdxUpper(:)
    logical, allocatable :: isConjugate(:)
    integer, allocatable :: mpiDDtype(:)
    integer, allocatable :: req(:),stat(:,:)
  end type cpInfoType

  ![ In cpRankS each entry corresponds to a process this rank has to send data
  ![ to. In cpRankR store the equivalent information for receiving.
  type(cpInfoType) :: cpRanksS, cpRanksR
  integer          :: cpRanksSlen, cpRanksRlen
  integer          :: cpRanksSlenP, cpRanksRlenP

  ![ These are used to conjugate the necessary negative-k amplitudes along ky.
  logical :: doConjugation
  integer :: cpIdxLower, cpIdxUpper

  ![ The portion in cpRankS/cpRankR corresponding to the positive-k magnitudes
  ![ will also be used when discarding negative-k amplitudes in r2c aliased
  ![ FFTs. This information is stored in noNkyaRanksS/noNkRanksR
  type(cpInfoType) :: noNkyaRanksS, noNkyaRanksR

![ endif for if (FFTmode==0) then else.
#endif

  ![ Derived type storing an allocatable list of integers.
  type intAlloc1D
    integer, allocatable :: is(:)
  end type intAlloc1D

  ![ User-defined time storing the information on how to transfer
  ![ de-alised data to aliased arrays.
  type d2aInfoType
    integer, allocatable :: rankID(:), idxLower(:,:), idxUpper(:,:)
    integer, allocatable :: mpiDDtype(:)
    integer, allocatable :: req(:),stat(:,:)
    integer              :: selfIdx 
    type(intAlloc1D)     :: remap(3)
  end type d2aInfoType

  ![ In d2aRanksS each entry corresponds to a process this rank sends data
  ![ to. In d2aRanksR the entries store equivalent info for receiving data.
  type(d2aInfoType) :: d2aRanksS, d2aRanksR
  integer           :: d2aSubarraysS, d2aSubarraysR

#else

  type(C_PTR) :: plan2Da_c2r, plan2Da_r2c    ![ Plans for aliased data.
  type(C_PTR) :: plan2D_c2r, plan2D_r2c      ![ Plans for de-aliased data.

  ![ Buffer which stores the real aliased data (after backward xy FFT).
  real(C_DOUBLE), pointer :: funRa(:,:)
  ![ Buffer which stores the real de-aliased data (after before forward xy FFT).
  real(C_DOUBLE), pointer :: funR(:,:)

#endif

  contains

#if (useMPI > 0)
#if (FFTmode == 0)
  subroutine init_FFTs
  ![ Initialize plans and dynamically allocate FFTW-related arrays.
  ![ Use the 2D r2c and c2r parallel transforms, decomposed along the
  ![ second dimension. For this to work we must gather all the data
  ![ along the first dimension.
  ![ A further optimization will be to use the TRANSPOSED_OUT feature
  ![ in FFTW, but that requires switching to kx along the first
  ![ dimension and ky along the second throughout the code.
  implicit none
  integer                          :: fftErrFlag
  character(C_CHAR, len=54)        :: wisdomFile
  character(len=54)                :: wisdomFileF
  character(len=5)                 :: cN1,cN2,cProcX,cProcY,cThreads,cMode
  integer(C_INTPTR_T)              :: allocLocal,fftM,fftN
  integer(C_INTPTR_T)              :: d1k,d1R
  integer                          :: intBuf,kyOffset
  integer, allocatable             :: NekyAll(:)
  integer                          :: dims,sizes(2),subsizes(2),starts(2)

  ![ Allow support for multithreaded FFTs.
  fftErrFlag = fftw_init_threads()
  if (fftErrFlag == 0) call abortSimulation('Error initializing FFTW threads')
  call fftw_mpi_init()
  if (FFTWthreads>1) call fftw_plan_with_nthreads(FFTWthreads)

  if (myRank == 0) write(*,'(" Creating/Loading FFTW wisdom...")')

  ![ Allocata array that holds the data distributed along kx
  ![ only (after a gathering of data long ky).
  allocate(funkyG(NekyG,Nekx(2)))

  ![ In order to gather the data from other processes along ky
  ![ we need to know how many ky's each process has.
  allocate(NekyAll(xyProcs(2)))
  intBuf = Nekx(1) 
  call MPI_Allgather(intBuf,1,MPI_INTEGER, &
                     NekyAll,1,MPI_INTEGER,xyCOMMs(1),errFlagMPI)

  if (fftMember) then

    ![ Import wisdom from file if there is any.
    if (myRank == 0) then
      fftErrFlag = 0
      write(cN1,'(i0)') NekyG
      write(cN2,'(i0)') NekxG
      write(cProcX,'(i0)') xyProcs(1)
      write(cProcY,'(i0)') xyProcs(2)
      write(cThreads,'(i0)') FFTWthreads
      write(cMode,'(i0)') FFTmode
      wisdomFileF = 'shroomWisdom_'//trim(cN1)//'x'//trim(cN2)//'_procs'//trim(cProcX)//'x' &
                    //trim(cProcY)//'threads'//trim(cThreads)//'FFTmode'//trim(cMode)//'.dat'
      write(wisdomFile,'(A54)') wisdomFileF
      fftErrFlag = fftw_import_wisdom_from_filename(trim(wisdomFile)//C_NULL_CHAR)
      if (fftErrFlag==0) then
        write(*,*) ' '
        write(*,'(" --> *Could not import FFTW wisdom from file. Will create. <--")')
        write(*,*) ' '
      endif
    else
      fftErrFlag = 1
    endif
    call MPI_BCAST(fftErrFlag,1,MPI_INTEGER,0,fftCOMM,errFlagMPI)
    if (fftErrFlag /= 0) then
      call fftw_mpi_broadcast_wisdom(fftCOMM)
    endif

    ![ ~~~~~~~~~~ Arrays and plans for aliased FFT along kx and ky ~~~~~~~~~~ ]!
    fftM = NekyaG
    fftN = NekxaG
    allocLocal = fftw_mpi_local_size_2d(fftN, fftM, fftCOMM, lockxa, lockxa0)

    if (lockxa /= Nekxa(2)) then
      write(*,'(" Shroom and FFTW use different MPI decomposition. myRank=",i4," | lockxa=",i4," | Nekxa=",i4)') &
        myRank,lockxa,Nekxa(2)
      call abortSimulation(' ')
    endif

    ![ Dynamically allocate pointer to complex and real arrays.
    funka_ptr = fftw_alloc_complex(int(allocLocal,C_SIZE_T))
    funRa_ptr = fftw_alloc_real(int(2*allocLocal,C_SIZE_T))

    ![ Convert C pointers to FORTRAN pointers.
    d1k = NekyaG
    d1R = 2*d1k
    call c_f_pointer( funka_ptr,  funka, [d1k,lockxa])
    call c_f_pointer( funRa_ptr,  funRa, [d1R,lockxa])

    ![ Plans for backward and forward FFTs.
    fftM = NyaG
    fftN = NxaG
    plan2Da_c2r = fftw_mpi_plan_dft_c2r_2d(fftN, fftM, funka, funRa, fftCOMM, FFTW_EXHAUSTIVE)
    plan2Da_r2c = fftw_mpi_plan_dft_r2c_2d(fftN, fftM, funRa, funka, fftCOMM, FFTW_EXHAUSTIVE)

    ![ Save wisdom if there is none available
    if (outputWisdom) then
      call fftw_mpi_gather_wisdom(fftCOMM)
      if (fftRank == 0) then
        fftErrFlag = fftw_export_wisdom_to_filename(trim(wisdomFile)//C_NULL_CHAR)
        if (fftErrFlag == 0) call abortSimulation('Error exporting wisdom to file.')
      endif
    endif

    ![ Create the subarrays used to gather the data from other processes along ky.
    allocate(kyGatherType%rankID(xyProcs(2)), kyGatherType%mpiDDtype(xyProcs(2)), &
             kyGatherType%kyIdxLower(xyProcs(2)), kyGatherType%kyIdxUpper(xyProcs(2)), &
             kyGatherType%req(xyProcs(2)), kyGatherType%stat(xyProcs(2),MPI_STATUS_SIZE))
    kyOffset = 0
    do lc = 1,xyProcs(2)
      kyGatherType%rankID(lc) = lc-1

      kyGatherType%kyIdxLower(lc) = kyOffset+1
      kyGatherType%kyIdxUpper(lc) = kyGatherType%kyIdxLower(lc)-1+NekyAll(lc)

      dims     = 2
      sizes    = [NekyG, Nekx(2)]
      subsizes = [NekyAll(lc), Nekx(2)]
      starts   = [kyOffset, 0]    ![ starts are zero-indexed.
      call MPI_TYPE_CREATE_SUBARRAY(dims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, kyGatherType%mpiDDtype(lc), errFlagMPI)
      call MPI_TYPE_COMMIT(kyGatherType%mpiDDtype(lc), errFlagMPI)
      kyOffset = kyOffset+NekyAll(lc)
    enddo

    ![ Initalize MPI derived data types used for transferring de-aliased
    ![ data to larger arrays that are aliased along kx. This involves placing
    ![ a band of zeros in the middle of the data, therefore compressing the
    ![ de-aliased data towards the 0 and xyprocs(1)-1 ranks.
    call init_toAliased

  endif

  deallocate(NekyAll)

  end subroutine init_FFTs
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine init_toAliased
  ![ First gather the data along ky using an MPI subarray type.
  ![ Then assuming de-aliased data is distributed across kx only, we will
  ![ transfer it to a larger aliased array that is also distributed across
  ![ kx. This involves communication in 1D, via fftCOMM.
  implicit none
  integer              :: intBuf
  integer, allocatable :: lockxaAll(:),NekxAll(:)
  integer, allocatable :: lockxa0All(:)
  integer              :: rc
  logical              :: rankCounted
  integer, allocatable :: d2aRankIDsX(:),d2aBoundsX(:,:)
  integer              :: d2aSubarraysX
  integer              :: dims,sizes(2),subsizes(2),starts(2)
  integer              :: ksChecked

  ![ Collect the number of elements in each process.
  allocate(lockxaAll(xyProcs(1)),NekxAll(xyProcs(1)))
  intBuf = lockxa
  call MPI_Allgather(intBuf,1,MPI_INTEGER, &
                     lockxaAll,1,MPI_INTEGER,fftCOMM,errFlagMPI)
  intBuf = Nekx(2) 
  call MPI_Allgather(intBuf,1,MPI_INTEGER, &
                     NekxAll,1,MPI_INTEGER,fftCOMM,errFlagMPI)
  ![ Also collect the index of the last element in the previous process along kx
  ![ and the index of the first element (amongst the global ones) in each process along ky.
  allocate(lockxa0All(xyProcs(1)))
  intBuf = lockxa0 
  call MPI_Allgather(intBuf,1,MPI_INTEGER, &
                     lockxa0All,1,MPI_INTEGER,fftCOMM,errFlagMPI)

  ![ Temporary arrays with arbitrarily large number of elements.
  allocate(d2aRankIDsX(xyProcs(1)*3),d2aBoundsX(xyProcs(1)*3,2))

  ![ Establish sending pattern along kx.
  ![ First do the positive-k amplitudes.
  d2aSubarraysX = 0
  do rc = 1,xyProcs(1)
    rankCounted = .false.
    do lckx = 1,Nekx(2)
      if ( ((firstkG(2)+lckx-1) <= NkxG) .and. & 
           ((firstkG(2)+lckx-1) >= lockxa0All(rc)+1) .and. &
           ((firstkG(2)+lckx-1) <= (lockxa0All(rc)+lockxaAll(rc))) ) then
        if (.not. rankCounted) then
          d2aSubarraysX               = d2aSubarraysX + 1
          rankCounted                 = .true.
          d2aRankIDsX(d2aSubarraysX)  = rc-1
          d2aBoundsX(d2aSubarraysX,1) = lckx
        endif
        d2aBoundsX(d2aSubarraysX,2) = lckx
      endif
    enddo
  enddo
  ![ Now do the negative-k amplitudes.
  do rc = 1,xyProcs(1)
    rankCounted = .false.
    do lckx = 1,Nekx(2)
      if ( ((firstkG(2)+lckx-1) > NkxG) .and. & 
           ((NekxG-(firstkG(2)-1+lckx-1)) <= (NekxaG-lockxa0All(rc))) .and. & 
           ((NekxG-(firstkG(2)-1+lckx-1)) > (NekxaG-(lockxa0All(rc)+lockxaAll(rc)))) ) then
        if (.not. rankCounted) then
          d2aSubarraysX               = d2aSubarraysX + 1
          rankCounted                 = .true.
          d2aRankIDsX(d2aSubarraysX)  = rc-1
          d2aBoundsX(d2aSubarraysX,1) = lckx
        endif
        d2aBoundsX(d2aSubarraysX,2) = lckx
      endif
    enddo
  enddo

  ![ Store the number of messages to send along kx in d2aSubbarraysS.
  d2aSubarraysS = d2aSubarraysX

  allocate(d2aRanksS%rankID(d2aSubarraysS),d2aRanksS%idxLower(d2aSubarraysS,2), &
           d2aRanksS%idxUpper(d2aSubarraysS,2),d2aRanksS%mpiDDtype(d2aSubarraysS), &
           d2aRanksS%req(d2aSubarraysS),d2aRanksS%stat(d2aSubarraysS,MPI_STATUS_SIZE))
  do lckx = 1,d2aSubarraysX
    rc  = lckx
    d2aRanksS%rankID(rc) = d2aRankIDsX(lckx)
    if (d2aRanksS%rankID(rc) == fftRank) d2aRanksS%selfIdx = rc

    d2aRanksS%idxLower(rc,1) = 1
    d2aRanksS%idxLower(rc,2) = d2aBoundsX(lckx,1)
    d2aRanksS%idxUpper(rc,1) = NekyG
    d2aRanksS%idxUpper(rc,2) = d2aBoundsX(lckx,2)

    dims     = 2
    sizes    = [NekyG, Nekx(2)]
    subsizes = [d2aRanksS%idxUpper(rc,1)-d2aRanksS%idxLower(rc,1)+1, &
                d2aRanksS%idxUpper(rc,2)-d2aRanksS%idxLower(rc,2)+1] 
    starts   = [d2aRanksS%idxLower(rc,1)-1, d2aRanksS%idxLower(rc,2)-1]    ![ starts are zero-indexed.
    call MPI_TYPE_CREATE_SUBARRAY(dims, sizes, subsizes, starts, &
      MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, d2aRanksS%mpiDDtype(rc), errFlagMPI)
    call MPI_TYPE_COMMIT(d2aRanksS%mpiDDtype(rc),errFlagMPI)
  enddo

  deallocate(d2aRankIDsX,d2aBoundsX)

  ![ Temporary arrays with arbitrarily large number of elements.
  allocate(d2aRankIDsX(xyProcs(1)*3),d2aBoundsX(xyProcs(1)*3,2))

  ![ Determine the receiving patterns, in the same order as the sends.
  ![ Positive-kx amplitudes.
  d2aSubarraysX = 0
  ksChecked     = 0
  do rc = 1,xyProcs(1)
    rankCounted = .false.
    do lckx = 1,Nekxa(2)
      if ( ((firstkaG(2)+lckx-1) <= NkxG) .and. &
           ((lockxa0+lckx) > ksChecked) .and. (lockxa0+lckx <= ksChecked+NekxAll(rc))) then
        if (.not. rankCounted) then
          d2aSubarraysX               = d2aSubarraysX+1
          rankCounted                 = .true.
          d2aRankIDsX(d2aSubarraysX)  = rc-1
          d2aBoundsX(d2aSubarraysX,1) = lckx
        endif
        d2aBoundsX(d2aSubarraysX,2)   = lckx
      endif
    enddo
    ksChecked = ksChecked+NekxAll(rc)
  enddo
  ![ Now do the negative-kx amplitudes.
  ksChecked     = 0
  do rc = 1,xyProcs(1)
    rankCounted = .false.
    do lckx = 1,Nekxa(2)
      if ( ((NekxaG-(firstkaG(2)-1+lckx-1)) < NkxG) .and. &
           ((NekxaG-(lockxa0+lckx-1)) <= (NekxG-ksChecked)) .and. &
           ((NekxaG-(lockxa0+lckx-1)) > (NekxG-(ksChecked+NekxAll(rc)))) ) then
        if (.not. rankCounted) then
          d2aSubarraysX               = d2aSubarraysX+1
          rankCounted                 = .true.
          d2aRankIDsX(d2aSubarraysX)  = rc-1
          d2aBoundsX(d2aSubarraysX,1) = lckx
        endif
        d2aBoundsX(d2aSubarraysX,2) = lckx
      endif
    enddo
    ksChecked = ksChecked+NekxAll(rc)
  enddo

  ![ Store the number of messages to receive along kx in d2aSubbarraysR.
  d2aSubarraysR = d2aSubarraysX

  allocate(d2aRanksR%rankID(d2aSubarraysR),d2aRanksR%idxLower(d2aSubarraysR,2), &
           d2aRanksR%idxUpper(d2aSubarraysR,2),d2aRanksR%mpiDDtype(d2aSubarraysR), &
           d2aRanksR%req(d2aSubarraysR),d2aRanksR%stat(d2aSubarraysR,MPI_STATUS_SIZE))
  do lckx = 1,d2aSubarraysX
    rc  = lckx
    d2aRanksR%rankID(rc) = d2aRankIDsX(lckx)
    if (d2aRanksR%rankID(rc) == fftRank) d2aRanksR%selfIdx = rc

    d2aRanksR%idxLower(rc,1) = 1
    d2aRanksR%idxLower(rc,2) = d2aBoundsX(lckx,1)
    d2aRanksR%idxUpper(rc,1) = NekyG
    d2aRanksR%idxUpper(rc,2) = d2aBoundsX(lckx,2)

    dims     = 2
    sizes    = [NekyaG, Nekxa(2)]
    subsizes = [d2aRanksR%idxUpper(rc,1)-d2aRanksR%idxLower(rc,1)+1, &
                d2aRanksR%idxUpper(rc,2)-d2aRanksR%idxLower(rc,2)+1]
    starts   = [d2aRanksR%idxLower(rc,1)-1, d2aRanksR%idxLower(rc,2)-1]    ![ starts are zero-indexed.
    call MPI_TYPE_CREATE_SUBARRAY(dims, sizes, subsizes, starts, &
      MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, d2aRanksR%mpiDDtype(rc), errFlagMPI)
    call MPI_TYPE_COMMIT(d2aRanksR%mpiDDtype(rc),errFlagMPI)
  enddo


  deallocate(d2aRankIDsX,d2aBoundsX)
  deallocate(lockxa0All)
  deallocate(lockxaAll,NekxAll)

  end subroutine init_toAliased
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine FFT2Da_c2r(fkaIn, fRaOut)
  ![ Complex-to-real 2D FFT of aliased data.
  ![ This routine assumes that (for fftMember==.true.) fkaIn has
  ![ all the kys and the kx's are distributed across fftCOMM.
  implicit none
  double complex, intent(in)    :: fkaIn(:,:)
  double precision, intent(out) :: fRaOut(:,:)

  if (fftMember) then
    ![ Relocate data from (NekyaG x Nekxa) array field array
    ![ to (NekyaG x Nekxa)=(NekyaG x Nxa) buffer array.
    funka = fkaIn
  
    call fftw_mpi_execute_dft_c2r(plan2Da_c2r, funka, funRa)    ![ Inverse transform along kx and ky.
  
    forall(lcx=1:Nxa(2),lcy=1:NyaG) fRaOut(lcy,lcx) = funRa(lcy,lcx)*rSqrtNxyaG    ![ Unitary normalization.
  else
    fRaOut = 0.0d0
  endif

  end subroutine FFT2Da_c2r
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine FFT2Da_r2c(fRaIn, fkaOut)
  ![ Real-to-complex 2D FFT of aliased data.
  ![ This routine assumes that (for fftMember==.true.) fRaIn has all the
  ![ y-data and the x-data distributed across fftCOMM.
  implicit none
  double precision, intent(in) :: fRaIn(:,:)
  double complex, intent(out)  :: fkaOut(:,:)
  
  if (fftMember) then
    ![ Relocate data into FFT buffer.
    forall(lcx=1:Nxa(2),lcy=1:NyaG) funRa(lcy,lcx) = fRaIn(lcy,lcx)

    call fftw_mpi_execute_dft_r2c(plan2Da_r2c, funRa, funka)    ![ Forward transform along x and y.

    fkaOut = funka*rSqrtNxyaG    ![ Unitary normalization.
  else
    fkaOut = dcmplx(0.0d0,0.0d0)
  endif

  end subroutine FFT2Da_r2c

#else
![ else for if (FFTmode==0) then else.

  subroutine init_FFTs
  ![ Initialize plans and dynamically allocate FFTW-related arrays.
  ![ Unfortunately FFTW does not support 2D-domain-decomposed 2D FFTs,
  ![ so we have to use 1D FFTs. Furthermore, it does not have a parallel
  ![ real-to-complex 1D FFT, so we have to use complex-to-complex 1D
  ![ parallel FFTs.
  implicit none
  integer                          :: fftErrFlag
  character(kind=C_CHAR,len=54)    :: wisdomFile
  character(len=54)                :: wisdomFileF
  character(len=5)                 :: cN1,cN2,cProcX,cProcY,cThreads,cMode
  integer(C_INTPTR_T)              :: allocLocal,fftN,howMany,m1D(1)
  integer(C_INTPTR_T)              :: locxa0
  integer(C_INTPTR_T)              :: locya,locya0
  integer(C_INTPTR_T)              :: lockx,locx,locx0
  integer(C_INTPTR_T)              :: locky,locy,locy0
  integer(C_INTPTR_T)              :: d1

  ![ Allow support for multithreaded FFTs.
  fftErrFlag = fftw_init_threads()
  if (fftErrFlag == 0) call abortSimulation('Error initializing FFTW threads')
  call fftw_mpi_init()
  if (FFTWthreads>1) call fftw_plan_with_nthreads(FFTWthreads)

  if (myRank == 0) write(*,'(" Creating/Loading FFTW wisdom...")')

  ![ Import wisdom from file if there is any.
  if (myRank == 0) then
    fftErrFlag = 0
    write(cN1,'(i0)') NekyG
    write(cN2,'(i0)') NekxG
    write(cProcX,'(i0)') xyProcs(1)
    write(cProcY,'(i0)') xyProcs(2)
    write(cThreads,'(i0)') FFTWthreads
    write(cMode,'(i0)') FFTmode
    wisdomFileF = 'shroomWisdom_'//trim(cN1)//'x'//trim(cN2)//'_procs'//trim(cProcX)//'x' &
                  //trim(cProcY)//'threads'//trim(cThreads)//'FFTmode'//trim(cMode)//'.dat'
    write(wisdomFile,'(A54)') wisdomFileF
    fftErrFlag = fftw_import_wisdom_from_filename(trim(wisdomFile)//C_NULL_CHAR)
    if ((fftErrFlag==0) .and. (myRank==0)) then
      write(*,*) ' '
      write(*,'(" --> *Could not import FFTW wisdom from file. Will create. <--")')
      write(*,*) ' '
    endif
  else
    fftErrFlag = 1
  endif
  call MPI_BCAST(fftErrFlag,1,MPI_INTEGER,0,MPI_COMM_WORLD,errFlagMPI)
  if (fftErrFlag /= 0) then
    call fftw_mpi_broadcast_wisdom(MPI_COMM_WORLD)
  endif

  ![ ~~~~~~~~~~ Arrays and plans for aliased FFT along kx ~~~~~~~~~~ ]!
  fftN    = NekxaG
  howMany = Nekxa(1)
  d1      = Nekxa(1)
  allocLocal = fftw_mpi_local_size_many_1d(fftN,howMany,xyCOMMs(2), &
               FFTW_BACKWARD,FFTW_EXHAUSTIVE,lockxa,lockxa0,locxa,locxa0)

  if (lockxa /= Nekxa(2)) then
    write(*,'(" Shroom and FFTW use different MPI decomposition. myRank=",i4," | lockxa=",i4," | Nekxa=",i4)') &
      myRank,lockxa,Nekxa(2)
    call abortSimulation(' ')
  endif

  ![ Dynamically allocate pointer to complex arrays that hold the complex
  ![ data and the data that is real along x (after FFT along x).
  funka_ptr   = fftw_alloc_complex(int(allocLocal,C_SIZE_T))
  funkya_ptr  = fftw_alloc_complex(int(allocLocal,C_SIZE_T))
  
  ![ Convert C pointers to FORTRAN pointers.
  call c_f_pointer(  funka_ptr,   funka, [d1,lockxa])
  call c_f_pointer( funkya_ptr,  funkya, [d1,locxa])
  allocate(funkyaT(lockxa,d1))  ![ Buffer with local transpose of funkya.

  ![ Plans for backward and forward FFTs along kx.
  m1D = [ NekxaG ]
  plan1Da_c2r_x = fftw_mpi_plan_many_dft(1,m1D,howMany, &
                  FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
                  funka,funkya,xyCOMMs(2),FFTW_BACKWARD,FFTW_EXHAUSTIVE)
  plan1Da_r2c_x = fftw_mpi_plan_many_dft(1,m1D,howMany, &
                  FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
                  funkya,funka,xyCOMMs(2),FFTW_FORWARD,FFTW_EXHAUSTIVE)

  ![ ~~~~~~~~~~ Arrays and plans for aliased FFT along ky ~~~~~~~~~~ ]!

  fftN    = NyaG
  howMany = Nxa(2)
  d1      = Nxa(2)
  allocLocal = fftw_mpi_local_size_many_1d(fftN,howMany,xyCOMMs(1), &
               FFTW_BACKWARD,FFTW_EXHAUSTIVE,lockya,lockya0,locya,locya0)

  if (lockya /= Nxa(1)) then
    write(*,'(" Shroom and FFTW use different MPI decomposition. myRank=",i4," | lockya=",i4," | Nya=",i4)') &
      myRank,lockya,Nxa(1)
    call abortSimulation(' ')
  endif

  ![ Dynamically allocate pointer to complex arrays that hold the data
  ![ that is complex along ky and real along x, and the array for
  ![ the real data (after performing both x and y FFTs).
  funkyaTc_ptr = fftw_alloc_complex(int(allocLocal,C_SIZE_T))
  funRa_ptr    = fftw_alloc_complex(int(allocLocal,C_SIZE_T))

  ![ Convert C pointers to FORTRAN pointers.
  call c_f_pointer(funkyaTc_ptr, funkyaTc, [d1,lockya])
  call c_f_pointer(   funRa_ptr,    funRa,  [d1,locya])

  ![ Plans for backward and forward FFTs along ky.
  m1D = [ NyaG ]
  plan1Da_c2r_y = fftw_mpi_plan_many_dft(1,m1D,howMany, &
                  FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
                  funkyaTc,funRa,xyCOMMs(1),FFTW_BACKWARD,FFTW_EXHAUSTIVE)
  plan1Da_r2c_y = fftw_mpi_plan_many_dft(1,m1D,howMany, &
                  FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
                  funRa,funkyaTc,xyCOMMs(1),FFTW_FORWARD,FFTW_EXHAUSTIVE)

  ![ Initalize MPI derived data types used for transferring de-aliased
  ![ data to larger aliased arrays.
  call init_toAliased

  ![ Create MPI derived data types used for enlarging an array
  ![ by appending the complex conjugates along ky.
  call init_appendRemoveConjugates


  if (initialOp == 0) then
    ![ ~~~~~~~~~~ Arrays and plans for de-aliased FFT along y ~~~~~~~~~~ ]!
    fftN    = NyG
    howMany = Nx(2)
    d1      = Nx(2)
    allocLocal = fftw_mpi_local_size_many_1d(fftN,howMany,xyCOMMs(1), &
                 FFTW_FORWARD,FFTW_EXHAUSTIVE,locy,locy0,locky,locky0)

    ![ Dynamically allocate pointer to complex arrays that hold the real
    ![ data and the array with the data that is complex along ky and real
    ![ along x (after having performed the y FFT).
    funkyTc_ptr = fftw_alloc_complex(int(allocLocal,C_SIZE_T))
    funR_ptr    = fftw_alloc_complex(int(allocLocal,C_SIZE_T))

    ![ Convert C pointers to FORTRAN pointers.
    call c_f_pointer(   funR_ptr,    funR,  [d1,locy])
    call c_f_pointer(funkyTc_ptr, funkyTc, [d1,locky])

    ![ Plans for backward and forward FFTs along y.
    m1D = [ NyG ]
    plan1D_r2c_y = fftw_mpi_plan_many_dft(1,m1D,howMany, &
                   FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
                   funR,funkyTc,xyCOMMs(1),FFTW_FORWARD,FFTW_EXHAUSTIVE)
    plan1D_c2r_y = fftw_mpi_plan_many_dft(1,m1D,howMany, &
                   FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
                   funkyTc,funR,xyCOMMs(1),FFTW_BACKWARD,FFTW_EXHAUSTIVE)

    ![ ~~~~~~~~~~ Arrays and plans for de-aliased FFT along x ~~~~~~~~~~ ]!
    fftN    = NxG
    howMany = Nekx(1)
    d1      = Nekx(1)
    allocLocal = fftw_mpi_local_size_many_1d(fftN,howMany,xyCOMMs(2), &
                 FFTW_FORWARD,FFTW_EXHAUSTIVE,locx,locx0,lockx,lockx0)

    ![ Dynamically allocate pointer to complex arrays that hold the data
    ![ that is real along x and the complex data (after FFT along x).
    funky_ptr  = fftw_alloc_complex(int(allocLocal,C_SIZE_T))
    funk_ptr   = fftw_alloc_complex(int(allocLocal,C_SIZE_T))
    
    ![ Convert C pointers to FORTRAN pointers.
    call c_f_pointer(funky_ptr,funky,[d1,locx])
    call c_f_pointer(funk_ptr, funk, [d1,lockx])
    allocate(funkyT(lockx,d1))  ![ Buffer with local transpose of funky.

    ![ Plans for backward and forward FFTs along x.
    m1D = [ NxG ]
    plan1D_r2c_x = fftw_mpi_plan_many_dft(1,m1D,howMany, &
                   FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
                   funky,funk,xyCOMMs(2),FFTW_FORWARD,FFTW_EXHAUSTIVE)
    plan1D_c2r_x = fftw_mpi_plan_many_dft(1,m1D,howMany, &
                   FFTW_MPI_DEFAULT_BLOCK,FFTW_MPI_DEFAULT_BLOCK, &
                   funk,funky,xyCOMMs(2),FFTW_BACKWARD,FFTW_EXHAUSTIVE)
  endif

  ![ Save wisdom if there is none available
  if (outputWisdom) then
    call fftw_mpi_gather_wisdom(MPI_COMM_WORLD)
    if (myRank == 0) then
      fftErrFlag = fftw_export_wisdom_to_filename(trim(wisdomFile)//C_NULL_CHAR)
      if (fftErrFlag == 0) call abortSimulation('Error exporting wisdom to file.')
    endif
  endif

  end subroutine init_FFTs
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine init_toAliased
  ![ The de-aliased data is distributed across all MPI proceses. We will transfer it
  ![ to a larger alased array that is also distributed across all MPI processes.
  ![ This involve communication in 2D. We will first identify the communication
  ![ pattern in each direction, and that in 2D is a union of the two.
  implicit none
  integer              :: intBuf
  integer, allocatable :: lockxaAll(:),NekyaAll(:),NekxAll(:),NekyAll(:)
  integer, allocatable :: lockxa0All(:),firstkyaGAll(:)
  integer              :: rc,IDs(2),blockID(2),mssgsInDir(3)
  logical              :: rankCounted
  integer, allocatable :: d2aRankIDsX(:),d2aBoundsX(:,:)
  integer, allocatable :: d2aRankIDsY(:),d2aBoundsY(:,:)
  integer              :: d2aSubarraysX,d2aSubarraysY
  integer              :: dims,sizes(2),subsizes(2),starts(2)
  integer              :: ksChecked

  ![ Collect the number of elements in each process.
  allocate(lockxaAll(xyProcs(1)),NekyaAll(xyProcs(2)),NekxAll(xyProcs(1)),NekyAll(xyProcs(2)))
  intBuf = lockxa
  call MPI_Allgather(intBuf,1,MPI_INTEGER, &
                     lockxaAll,1,MPI_INTEGER,xyCOMMs(2),errFlagMPI)
  intBuf = Nekxa(1) 
  call MPI_Allgather(intBuf,1,MPI_INTEGER, &
                     NekyaAll,1,MPI_INTEGER,xyCOMMs(1),errFlagMPI)
  intBuf = Nekx(2) 
  call MPI_Allgather(intBuf,1,MPI_INTEGER, &
                     NekxAll,1,MPI_INTEGER,xyCOMMs(2),errFlagMPI)
  intBuf = Nekx(1) 
  call MPI_Allgather(intBuf,1,MPI_INTEGER, &
                     NekyAll,1,MPI_INTEGER,xyCOMMs(1),errFlagMPI)
  ![ Also collect the index of the last element in the previous process along kx
  ![ and the index of the first element (amongst the global ones) in each process along ky.
  allocate(lockxa0All(xyProcs(1)),firstkyaGAll(xyProcs(2)))
  intBuf = lockxa0 
  call MPI_Allgather(intBuf,1,MPI_INTEGER, &
                     lockxa0All,1,MPI_INTEGER,xyCOMMs(2),errFlagMPI)
  intBuf = firstkaG(1) 
  call MPI_Allgather(intBuf,1,MPI_INTEGER, &
                     firstkyaGAll,1,MPI_INTEGER,xyCOMMs(1),errFlagMPI)

  ![ Temporary arrays with arbitrarily large number of elements.
  allocate(d2aRankIDsX(xyProcs(1)*3),d2aBoundsX(xyProcs(1)*3,2))

  ![ Establish sending pattern along kx. This is easier than ky
  ![ because no conjugate-pair appending needs to happen along kx
  ![ (so lockxa and lockxa0 can be used here).
  ![ First do the positive-k amplitudes.
  d2aSubarraysX = 0
  do rc = 1,xyProcs(1)
    rankCounted = .false.
    do lckx = 1,Nekx(2)
      if ( ((firstkG(2)+lckx-1) <= NkxG) .and. & 
           ((firstkG(2)+lckx-1) >= lockxa0All(rc)+1) .and. &
           ((firstkG(2)+lckx-1) <= (lockxa0All(rc)+lockxaAll(rc))) ) then
        if (.not. rankCounted) then
          d2aSubarraysX               = d2aSubarraysX + 1
          rankCounted                 = .true.
          d2aRankIDsX(d2aSubarraysX)  = rc-1
          d2aBoundsX(d2aSubarraysX,1) = lckx
        endif
        d2aBoundsX(d2aSubarraysX,2) = lckx
      endif
    enddo
  enddo
  ![ Now do the negative-k amplitudes.
  do rc = 1,xyProcs(1)
    rankCounted = .false.
    do lckx = 1,Nekx(2)
      if ( ((firstkG(2)+lckx-1) > NkxG) .and. & 
           ((NekxG-(firstkG(2)-1+lckx-1)) <= (NekxaG-lockxa0All(rc))) .and. & 
           ((NekxG-(firstkG(2)-1+lckx-1)) > (NekxaG-(lockxa0All(rc)+lockxaAll(rc)))) ) then
        if (.not. rankCounted) then
          d2aSubarraysX               = d2aSubarraysX + 1
          rankCounted                 = .true.
          d2aRankIDsX(d2aSubarraysX)  = rc-1
          d2aBoundsX(d2aSubarraysX,1) = lckx
        endif
        d2aBoundsX(d2aSubarraysX,2) = lckx
      endif
    enddo
  enddo

  ![ Temporary arrays with arbitrarily large number of elements.
  allocate(d2aRankIDsY(xyProcs(2)*3),d2aBoundsY(xyProcs(2)*3,2))

  ![ Establish the communication pattern along ky.
  d2aSubarraysY = 0
  do rc = 1,xyProcs(2)
    rankCounted = .false.
    do lcky = 1,Nekx(1)
      if ( ((firstkG(1)+lcky-1) >= firstkyaGAll(rc)) .and. &
           ((firstkG(1)+lcky-1) < (firstkyaGAll(rc)+NekyaAll(rc))) ) then
        if (.not. rankCounted) then
          d2aSubarraysY               = d2aSubarraysY + 1
          rankCounted                 = .true.
          d2aRankIDsY(d2aSubarraysY)  = rc-1
          d2aBoundsY(d2aSubarraysY,1) = lcky
        endif
        d2aBoundsY(d2aSubarraysY,2) = lcky
      endif
    enddo
  enddo

  ![ Form the union of the communication pattern along kx and ky.
  d2aSubarraysS = d2aSubarraysX*d2aSubarraysY

  mssgsInDir = 0
  call MPI_CART_COORDS(cartCOMM, cartRank, 2, blockID, errFlagMPI)  ![ Get the 2D ID of this rank.

  allocate(d2aRanksS%rankID(d2aSubarraysS),d2aRanksS%idxLower(d2aSubarraysS,2), &
           d2aRanksS%idxUpper(d2aSubarraysS,2),d2aRanksS%mpiDDtype(d2aSubarraysS), &
           d2aRanksS%req(d2aSubarraysS),d2aRanksS%stat(d2aSubarraysS,MPI_STATUS_SIZE))
  do lckx = 1,d2aSubarraysX
    do lcky = 1,d2aSubarraysY
      rc  = (lckx-1)*d2aSubarraysY+lcky
      IDs = [d2aRankIDsX(lckx),d2aRankIDsY(lcky)]
      call MPI_CART_RANK(cartCOMM,IDs,d2aRanksS%rankID(rc),errFlagMPI)  ![ Recall reverse order of communicator.
      if (d2aRanksS%rankID(rc) == cartRank) d2aRanksS%selfIdx = rc

      if ((blockID(1) .ne. IDs(1)) .and. (blockID(2) .ne. IDs(2))) then
        mssgsInDir(3) = mssgsInDir(3)+1  ![ Messages sent along ky and kx.
      elseif ((blockID(1) .ne. IDs(1)) .and. (blockID(2) == IDs(2))) then
        mssgsInDir(2) = mssgsInDir(2)+1  ![ Messages sent along kx. 
      elseif ((blockID(1) == IDs(1)) .and. (blockID(2) .ne. IDs(2))) then
        mssgsInDir(1) = mssgsInDir(1)+1  ![ Messages sent along ky. 
      endif

      d2aRanksS%idxLower(rc,1) = d2aBoundsY(lcky,1)
      d2aRanksS%idxLower(rc,2) = d2aBoundsX(lckx,1)
      d2aRanksS%idxUpper(rc,1) = d2aBoundsY(lcky,2)
      d2aRanksS%idxUpper(rc,2) = d2aBoundsX(lckx,2)

      dims     = 2
      sizes    = Nekx
      subsizes = [d2aRanksS%idxUpper(rc,1)-d2aRanksS%idxLower(rc,1)+1, &
                  d2aRanksS%idxUpper(rc,2)-d2aRanksS%idxLower(rc,2)+1] 
      starts   = [d2aRanksS%idxLower(rc,1)-1, d2aRanksS%idxLower(rc,2)-1]    ![ starts are zero-indexed.
      call MPI_TYPE_CREATE_SUBARRAY(dims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, d2aRanksS%mpiDDtype(rc), errFlagMPI)
      call MPI_TYPE_COMMIT(d2aRanksS%mpiDDtype(rc),errFlagMPI)
    enddo
  enddo

  ![ In MIT's Engaging cluster posting all the MPI_IRECVs followed by all the
  ![ MPI_ISENDs can lead to a hang. Perhaps this is hardware dependent. We will
  ![ try to fix this by first posting the communications along one dimesion,
  ![ followed by the communications along the other dimension, and finishing
  ![ with the communication that happens along both directions.
  do lc = 1, 3
    allocate(d2aRanksS%remap(lc)%is(mssgsInDir(lc)))
  enddo
  mssgsInDir = 0
  do rc = 1, d2aSubarraysS
    call MPI_CART_COORDS(cartCOMM, d2aRanksS%rankID(rc), 2, IDs, errFlagMPI)
    if ((blockID(1) .ne. IDs(1)) .and. (blockID(2) .ne. IDs(2))) then
      mssgsInDir(3) = mssgsInDir(3)+1
      d2aRanksS%remap(3)%is(mssgsInDir(3)) = rc
    elseif ((blockID(1) .ne. IDs(1)) .and. (blockID(2) == IDs(2))) then
      mssgsInDir(2) = mssgsInDir(2)+1
      d2aRanksS%remap(2)%is(mssgsInDir(2)) = rc
    elseif ((blockID(1) == IDs(1)) .and. (blockID(2) .ne. IDs(2))) then
      mssgsInDir(1) = mssgsInDir(1)+1
      d2aRanksS%remap(1)%is(mssgsInDir(1)) = rc
    endif
  enddo

  deallocate(d2aRankIDsY,d2aBoundsY)
  deallocate(d2aRankIDsX,d2aBoundsX)

  ![ Temporary arrays with arbitrarily large number of elements.
  allocate(d2aRankIDsX(xyProcs(1)*3),d2aBoundsX(xyProcs(1)*3,2))

  ![ Determine the receiving patterns, in the same order as the sends.
  ![ Positive-kx amplitudes.
  d2aSubarraysX = 0
  ksChecked     = 0
  do rc = 1,xyProcs(1)
    rankCounted = .false.
    do lckx = 1,Nekxa(2)
      if ( ((firstkaG(2)+lckx-1) <= NkxG) .and. & 
           ((lockxa0+lckx) > ksChecked) .and. (lockxa0+lckx <= ksChecked+NekxAll(rc))) then
        if (.not. rankCounted) then
          d2aSubarraysX               = d2aSubarraysX+1
          rankCounted                 = .true.
          d2aRankIDsX(d2aSubarraysX)  = rc-1
          d2aBoundsX(d2aSubarraysX,1) = lckx
        endif
        d2aBoundsX(d2aSubarraysX,2)   = lckx
      endif
    enddo
    ksChecked = ksChecked+NekxAll(rc)
  enddo
  ![ Now do the negative-kx amplitudes.
  ksChecked     = 0
  do rc = 1,xyProcs(1)
    rankCounted = .false.
    do lckx = 1,Nekxa(2)
      if ( ((NekxaG-(firstkaG(2)-1+lckx-1)) < NkxG) .and. & 
           ((NekxaG-(lockxa0+lckx-1)) <= (NekxG-ksChecked)) .and. &
           ((NekxaG-(lockxa0+lckx-1)) > (NekxG-(ksChecked+NekxAll(rc)))) ) then
        if (.not. rankCounted) then
          d2aSubarraysX               = d2aSubarraysX+1
          rankCounted                 = .true.
          d2aRankIDsX(d2aSubarraysX)  = rc-1
          d2aBoundsX(d2aSubarraysX,1) = lckx
        endif
        d2aBoundsX(d2aSubarraysX,2) = lckx
      endif
    enddo
    ksChecked = ksChecked+NekxAll(rc)
  enddo

  ![ Temporary arrays with arbitrarily large number of elements.
  allocate(d2aRankIDsY(xyProcs(2)*3),d2aBoundsY(xyProcs(2)*3,2))
  ![ Receiving patterns along y.
  d2aSubarraysY = 0
  ksChecked     = 0
  do rc = 1,xyProcs(2)
    rankCounted = .false.
    do lcky = 1,Nekxa(1)
      if ( ((firstkaG(1)+lcky-1) <= NkyG) .and. & 
           ((firstkaG(1)-1+lcky) > ksChecked) .and. (firstkaG(1)-1+lcky <= ksChecked+NekyAll(rc))) then
        if (.not. rankCounted) then
          d2aSubarraysY               = d2aSubarraysY+1
          rankCounted                 = .true.
          d2aRankIDsY(d2aSubarraysY)  = rc-1
          d2aBoundsY(d2aSubarraysY,1) = lcky
        endif
        d2aBoundsY(d2aSubarraysY,2)   = lcky
      endif
    enddo
    ksChecked = ksChecked+NekyAll(rc)
  enddo

  ![ The subarrays sent/receive are given by the union of the
  ![ kx and ky communication patterns above.
  d2aSubarraysR = d2aSubarraysX*d2aSubarraysY

  mssgsInDir = 0

  allocate(d2aRanksR%rankID(d2aSubarraysR),d2aRanksR%idxLower(d2aSubarraysR,2), &
           d2aRanksR%idxUpper(d2aSubarraysR,2),d2aRanksR%mpiDDtype(d2aSubarraysR), &
           d2aRanksR%req(d2aSubarraysR),d2aRanksR%stat(d2aSubarraysR,MPI_STATUS_SIZE))
  do lckx = 1,d2aSubarraysX
    do lcky = 1,d2aSubarraysY
      rc  = (lckx-1)*d2aSubarraysY+lcky
      IDs = [d2aRankIDsX(lckx), d2aRankIDsY(lcky)]  ![ Recall reverse order of communicator.
      call MPI_CART_RANK(cartCOMM,IDs,d2aRanksR%rankID(rc),errFlagMPI)
      if (d2aRanksR%rankID(rc) == cartRank) d2aRanksR%selfIdx = rc

      if ((blockID(1) .ne. IDs(1)) .and. (blockID(2) .ne. IDs(2))) then
        mssgsInDir(3) = mssgsInDir(3)+1  ![ Messages received along ky and kx.
      elseif ((blockID(1) .ne. IDs(1)) .and. (blockID(2) == IDs(2))) then
        mssgsInDir(2) = mssgsInDir(2)+1  ![ Messages received along kx. 
      elseif ((blockID(1) == IDs(1)) .and. (blockID(2) .ne. IDs(2))) then
        mssgsInDir(1) = mssgsInDir(1)+1  ![ Messages received along ky. 
      endif

      d2aRanksR%idxLower(rc,1) = d2aBoundsY(lcky,1)
      d2aRanksR%idxLower(rc,2) = d2aBoundsX(lckx,1)
      d2aRanksR%idxUpper(rc,1) = d2aBoundsY(lcky,2)
      d2aRanksR%idxUpper(rc,2) = d2aBoundsX(lckx,2)

      dims     = 2
      sizes    = Nekxa
      subsizes = [d2aRanksR%idxUpper(rc,1)-d2aRanksR%idxLower(rc,1)+1, &
                  d2aRanksR%idxUpper(rc,2)-d2aRanksR%idxLower(rc,2)+1] 
      starts   = [d2aRanksR%idxLower(rc,1)-1, d2aRanksR%idxLower(rc,2)-1]    ![ starts are zero-indexed.
      call MPI_TYPE_CREATE_SUBARRAY(dims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, d2aRanksR%mpiDDtype(rc), errFlagMPI)
      call MPI_TYPE_COMMIT(d2aRanksR%mpiDDtype(rc),errFlagMPI)
    enddo
  enddo

  ![ In MIT's Engaging cluster posting all the MPI_IRECVs followed by all the
  ![ MPI_ISENDs can lead to a hang. Perhaps this is hardware dependent. We will
  ![ try to fix this by first posting the communications along one dimesion,
  ![ followed by the communications along the other dimension, and finishing
  ![ with the communication that happens along both directions.
  do lc = 1, 3
    allocate(d2aRanksR%remap(lc)%is(mssgsInDir(lc)))
  enddo
  mssgsInDir = 0
  do rc = 1, d2aSubarraysR
    call MPI_CART_COORDS(cartCOMM, d2aRanksR%rankID(rc), 2, IDs, errFlagMPI)
    if ((blockID(1) .ne. IDs(1)) .and. (blockID(2) .ne. IDs(2))) then
      mssgsInDir(3) = mssgsInDir(3)+1
      d2aRanksR%remap(3)%is(mssgsInDir(3)) = rc
    elseif ((blockID(1) .ne. IDs(1)) .and. (blockID(2) == IDs(2))) then
      mssgsInDir(2) = mssgsInDir(2)+1
      d2aRanksR%remap(2)%is(mssgsInDir(2)) = rc
    elseif ((blockID(1) == IDs(1)) .and. (blockID(2) .ne. IDs(2))) then
      mssgsInDir(1) = mssgsInDir(1)+1
      d2aRanksR%remap(1)%is(mssgsInDir(1)) = rc
    endif
  enddo


  deallocate(d2aRankIDsY,d2aBoundsY)
  deallocate(d2aRankIDsX,d2aBoundsX)
  deallocate(lockxa0All,firstkyaGAll)
  deallocate(lockxaAll,NekyaAll,NekxAll,NekyAll)

  end subroutine init_toAliased
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine init_appendRemoveConjugates
  ![ Initialize the infrastructure for appending the complex conjugates
  ![ along ky. In general case this requires a tricky communication pattern.
  ![ We will first try the following:
  ![   1. Give every rank the knowledge of what the first element of every
  ![      rank corresponds to in the global conjugate array (broadcast lockya).
  ![   2. Create an array of user-defined types (cpInfoS, for conjugate pair
  ![      send info). Each entry of this array corresponds to a process that
  ![      this rank needs to communicate with. Each cpInfoS type stores
  ![        a) The rank ID of the process it communicates with.
  ![        b) The first and last index (data bounds) of the Nekya-sized
  ![           array (funkyaT) that will be communicated to this rank.
  ![        c) A flag indicating whether it is sending it to be received
  ![           as a positive-k amplitude, or a negative-k amplitude (conjugate pair).
  ![ Then in the FFT2Da routines one would loop over the array of cpInfoS,
  ![ sending the right data to each process. At the same time one would be
  ![ receiving corresponding data and organizing it into funkyaTc.
  ![
  ![ This routine also creates the MPI derived data types used to remove
  ![ the complex conjugates (for the r2c transform).
  implicit none
  integer, allocatable :: lockyaAll(:),lockya0All(:), NekyaAll(:)
  integer              :: intBuf
  integer              :: rc, lcb
  logical              :: rankCounted
  integer, allocatable :: cpRankIDsTMP(:),cpBoundsTMP(:,:)
  logical, allocatable :: isConjugateTMP(:)
  integer              :: dims,sizes(2),subsizes(2),starts(2),cpBlocks
  integer, allocatable :: cpDisplS(:),cpDisplR(:)
  integer              :: kysChecked

  ![ Collect lockya and lockya0 from all processes.
  allocate(lockyaAll(xyProcs(2)),lockya0All(xyProcs(2)),NekyaAll(xyProcs(2)))
  intBuf = lockya
  call MPI_Allgather(intBuf,1,MPI_INTEGER, &
                     lockyaAll,1,MPI_INTEGER,xyCOMMs(1),errFlagMPI)
  intBuf = lockya0 
  call MPI_Allgather(intBuf,1,MPI_INTEGER, &
                     lockya0All,1,MPI_INTEGER,xyCOMMs(1),errFlagMPI)
  ![ Also collect the number of elements in each process.
  intBuf = Nekxa(1) 
  call MPI_Allgather(intBuf,1,MPI_INTEGER, &
                     NekyaAll,1,MPI_INTEGER,xyCOMMs(1),errFlagMPI)

  ![ Temporary arrays with arbitrarily large number of elements.
  allocate(cpRankIDsTMP(xyProcs(2)*3),cpBoundsTMP(xyProcs(2)*3,2),isConjugateTMP(xyProcs(2)*3))

  ![ Determine the number of processes this rank sends data to.
  cpRanksSlen = 0
  ![ Positive-k processes this rank sends data to.
  do rc = 1, xyProcs(2)
    rankCounted = .false.
    do lcky = 1,Nekxa(1)
      if ( ((firstkaG(1)+lcky-1) >= lockya0All(rc)+1) .and. &
           ((firstkaG(1)+lcky-1) <= (lockya0All(rc)+lockyaAll(rc))) ) then
        if (.not. rankCounted) then
          cpRanksSlen                 = cpRanksSlen + 1
          rankCounted                 = .true.
          cpRankIDsTMP(cpRanksSlen)   = rc-1
          cpBoundsTMP(cpRanksSlen,1)  = lcky
          isConjugateTMP(cpRanksSlen) = .false.
        endif
        cpBoundsTMP(cpRanksSlen,2) = lcky
      endif
    enddo
  enddo
  ![ Number of ranks communicating in positive-k amplitudes only.
  cpRanksSlenP = cpRanksSlen
  ![ Negative-k processes this rank sends data to.
  do rc = 1, xyProcs(2)
    rankCounted = .false.
    do lcky = 1,Nekxa(1)
      if ( (kya(lcky) > kyaG(1)) .and. (kya(lcky) < kyaG(NekyaG)) .and. &  ![ Exclude ky=0 and highest ky.
           (2*NekyaG-1-(firstkaG(1)+lcky-1) >= lockya0All(rc)) .and. &
           (2*NekyaG-1-(firstkaG(1)+lcky-1) < (lockya0All(rc)+lockyaAll(rc))) ) then
        if (.not. rankCounted) then
          cpRanksSlen                 = cpRanksSlen + 1
          rankCounted                 = .true.
          cpRankIDsTMP(cpRanksSlen)   = rc-1
          cpBoundsTMP(cpRanksSlen,1)  = lcky
          isConjugateTMP(cpRanksSlen) = .true.
        endif
        cpBoundsTMP(cpRanksSlen,2) = lcky
      endif
    enddo
  enddo

  allocate(cpRanksS%rankID(cpRanksSlen),cpRanksS%kyIdxLower(cpRanksSlen), &
           cpRanksS%kyIdxUpper(cpRanksSlen),cpRanksS%isConjugate(cpRanksSlen), &
           cpRanksS%mpiDDtype(cpRanksSlen),cpRanksS%req(cpRanksSlen), &
           cpRanksS%stat(cpRanksSlen,MPI_STATUS_SIZE))
  do rc = 1,cpRanksSlen
    cpRanksS%rankID(rc)      = cpRankIDsTMP(rc)
    cpRanksS%kyIdxLower(rc)  = cpBoundsTMP(rc,1)
    cpRanksS%kyIdxUpper(rc)  = cpBoundsTMP(rc,2)
    cpRanksS%isConjugate(rc) = isConjugateTMP(rc)

    cpBlocks = cpRanksS%kyIdxUpper(rc)-cpRanksS%kyIdxLower(rc)+1
    if (.not. cpRanksS%isConjugate(rc)) then
      ![ The sends corresponding to positive-k amplitudes can be done
      ![ wth the MPI derived datatype 'subarray'.
      dims     = 2
      sizes    = [Nxa(2), Nekxa(1)]
      subsizes = [Nxa(2), cpBlocks] 
      starts   = [0, cpRanksS%kyIdxLower(rc)-1]    ![ starts are zero-indexed.
      call MPI_TYPE_CREATE_SUBARRAY(dims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, cpRanksS%mpiDDtype(rc), errFlagMPI)
      call MPI_TYPE_COMMIT(cpRanksS%mpiDDtype(rc),errFlagMPI)
    else
      ![ The sends corresponding to complex conjugation consists of
      ![ contiguous blocks of length Nxa, placed in decreasing order on
      ![ the receiving end.
      ![ For this we will use the MPI derived datatype 'indexed_block'.
      ![ MF: is it possible for the sender to use a subarray type, but
      ![     the receiver uses an indexed_block type?
      allocate(cpDisplS(cpBlocks))
      do lcb = 1, cpBlocks
        cpDisplS(lcb) = (cpRanksS%kyIdxLower(rc)-1+lcb-1)*Nxa(2)
      enddo 
      call MPI_TYPE_CREATE_INDEXED_BLOCK(cpBlocks,Nxa(2),cpDisplS, &
        MPI_DOUBLE_COMPLEX,cpRanksS%mpiDDtype(rc),errFlagMPI)
      call MPI_TYPE_COMMIT(cpRanksS%mpiDDtype(rc),errFlagMPI)
      deallocate(cpDisplS)
    endif
  enddo

  ![ These processes doing the sending are also the ones that receive
  ![ data in the same pattern when discarding negative-k amplitudes
  ![ as part of the r2c FFT.
  allocate(noNkyaRanksR%rankID(cpRanksSlenP),noNkyaRanksR%kyIdxLower(cpRanksSlenP), &
           noNkyaRanksR%kyIdxUpper(cpRanksSlenP),noNkyaRanksR%isConjugate(cpRanksSlenP), &
           noNkyaRanksR%mpiDDtype(cpRanksSlenP),noNkyaRanksR%req(cpRanksSlenP), &
           noNkyaRanksR%stat(cpRanksSlenP,MPI_STATUS_SIZE))
  do rc = 1,cpRanksSlenP
    noNkyaRanksR%rankID(rc)      = cpRankIDsTMP(rc)
    noNkyaRanksR%kyIdxLower(rc)  = cpBoundsTMP(rc,1)
    noNkyaRanksR%kyIdxUpper(rc)  = cpBoundsTMP(rc,2)
    noNkyaRanksR%isConjugate(rc) = isConjugateTMP(rc)

    cpBlocks = noNkyaRanksR%kyIdxUpper(rc)-noNkyaRanksR%kyIdxLower(rc)+1
    if (.not. noNkyaRanksR%isConjugate(rc)) then
      ![ These receives corresponding to positive-k amplitudes can be done
      ![ wth the MPI derived datatype 'subarray'.
      dims     = 2
      sizes    = [Nxa(2), Nekxa(1)]
      subsizes = [Nxa(2), cpBlocks] 
      starts   = [0, noNkyaRanksR%kyIdxLower(rc)-1]    ![ starts are zero-indexed.
      call MPI_TYPE_CREATE_SUBARRAY(dims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, noNkyaRanksR%mpiDDtype(rc), errFlagMPI)
      call MPI_TYPE_COMMIT(noNkyaRanksR%mpiDDtype(rc),errFlagMPI)
    endif
  enddo

  deallocate(cpRankIDsTMP,cpBoundsTMP,isConjugateTMP)

  ![ Now establish the processes this rank receives data from and
  ![ store the required information to place the data in the right place.

  ![ Temporary arrays with arbitrarily large number of elements.
  allocate(cpRankIDsTMP(xyProcs(2)*3),cpBoundsTMP(xyProcs(2)*3,2),isConjugateTMP(xyProcs(2)*3))

  ![ Given lockya0All and the number of elements in each process, NekyaAll,
  ![ determine the processes that this rank receives data from corresponding
  ![ to positive-k amplitudes.
  cpRanksRlen = 0
  kysChecked  = 0
  do rc = 1,xyProcs(2)
    rankCounted = .false.
    do lcky = 1,Nxa(1) 
      if ((lockya0+lcky > kysChecked) .and. (lockya0+lcky <= kysChecked+NekyaAll(rc))) then
        if (.not. rankCounted) then
          cpRanksRlen                 = cpRanksRlen+1
          rankCounted                 = .true.
          cpRankIDsTMP(cpRanksRlen)   = rc-1
          cpBoundsTMP(cpRanksRlen,1)  = lcky
          isConjugateTMP(cpRanksRlen) = .false.
        endif
        cpBoundsTMP(cpRanksRlen,2)  = lcky
      endif
    enddo
    kysChecked = kysChecked+NekyaAll(rc)
  enddo
  ![ Number of ranks communicating in positive-k amplitudes only.
  cpRanksRlenP = cpRanksRlen
  ![ Processes this rank receives data from in order to fill the negative-k amplitudes.
  kysChecked = 0
  do rc = 1,xyProcs(2)
    rankCounted = .false.
    do lcky = Nxa(1),1,-1
      if ( (lockya0+lcky > NkyaG) .and. &    ![ Exclude kya<=kyaMax.
           (NyaG-(lockya0+lcky)+1 > kysChecked-1) .and. &
           (NyaG-(lockya0+lcky)+1 <= kysChecked+NekyaAll(rc)-1) ) then ![ Minus one to exclude kya=0.
        if (.not. rankCounted) then
          cpRanksRlen                 = cpRanksRlen+1
          rankCounted                 = .true.
          cpRankIDsTMP(cpRanksRlen)   = rc-1
          cpBoundsTMP(cpRanksRlen,2)  = lcky  ![ Upper bound here. We will place blocks in decending order.
          isConjugateTMP(cpRanksRlen) = .true.
        endif
        cpBoundsTMP(cpRanksRlen,1)  = lcky
      endif
    enddo
    kysChecked = kysChecked+NekyaAll(rc)
  enddo

  allocate(cpRanksR%rankID(cpRanksRlen),cpRanksR%kyIdxLower(cpRanksRlen), &
           cpRanksR%kyIdxUpper(cpRanksRlen),cpRanksR%isConjugate(cpRanksRlen), &
           cpRanksR%mpiDDtype(cpRanksRlen),cpRanksR%req(cpRanksRlen), &
           cpRanksR%stat(cpRanksRlen,MPI_STATUS_SIZE))
  do rc = 1,cpRanksRlen
    cpRanksR%rankID(rc)      = cpRankIDsTMP(rc)
    cpRanksR%kyIdxLower(rc)  = cpBoundsTMP(rc,1)
    cpRanksR%kyIdxUpper(rc)  = cpBoundsTMP(rc,2)
    cpRanksR%isConjugate(rc) = isConjugateTMP(rc)

    cpBlocks = cpRanksR%kyIdxUpper(rc)-cpRanksR%kyIdxLower(rc)+1
    if (.not. cpRanksR%isConjugate(rc)) then
      ![ The receives corresponding to positive-k amplitudes can be done
      ![ wth the MPI derived datatype 'subarray'.
      dims     = 2
      sizes    = [Nxa(2), Nxa(1)]
      subsizes = [Nxa(2), cpBlocks] 
      starts   = [0, cpRanksR%kyIdxLower(rc)-1]    ![ starts are zero-indexed.
      call MPI_TYPE_CREATE_SUBARRAY(dims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, cpRanksR%mpiDDtype(rc), errFlagMPI)
      call MPI_TYPE_COMMIT(cpRanksR%mpiDDtype(rc),errFlagMPI)
    else
      ![ The receives corresponding to complex conjugation consists of
      ![ contiguous blocks of length Nxa, placed in decreasing order on
      ![ the receiving end.
      ![ For this we will use the MPI derived datatype 'indexed_block'.
      ![ MF: is it possible for the sender to use a subarray type, but
      ![     the receiver uses an indexed_block type?
      allocate(cpDisplR(cpBlocks))
      do lcb = 1, cpBlocks
        cpDisplR(lcb) = (cpRanksR%kyIdxUpper(rc)-1-(lcb-1))*Nxa(2)
      enddo 
      call MPI_TYPE_CREATE_INDEXED_BLOCK(cpBlocks,Nxa(2),cpDisplR, &
        MPI_DOUBLE_COMPLEX,cpRanksR%mpiDDtype(rc),errFlagMPI)
      call MPI_TYPE_COMMIT(cpRanksR%mpiDDtype(rc),errFlagMPI)
      deallocate(cpDisplR)
    endif
  enddo

  ![ The receive pattern for the positive-k amplitudes is the same as that
  ![ for sends upon discaring negative-k amplitudes as part of the r2c FFT.
  allocate(noNkyaRanksS%rankID(cpRanksRlenP),noNkyaRanksS%kyIdxLower(cpRanksRlenP), &
           noNkyaRanksS%kyIdxUpper(cpRanksRlenP),noNkyaRanksS%isConjugate(cpRanksRlenP), &
           noNkyaRanksS%mpiDDtype(cpRanksRlenP),noNkyaRanksS%req(cpRanksRlenP), &
           noNkyaRanksS%stat(cpRanksRlenP,MPI_STATUS_SIZE))
  do rc = 1,cpRanksRlenP
    noNkyaRanksS%rankID(rc)      = cpRankIDsTMP(rc)
    noNkyaRanksS%kyIdxLower(rc)  = cpBoundsTMP(rc,1)
    noNkyaRanksS%kyIdxUpper(rc)  = cpBoundsTMP(rc,2)
    noNkyaRanksS%isConjugate(rc) = isConjugateTMP(rc)

    cpBlocks = noNkyaRanksS%kyIdxUpper(rc)-noNkyaRanksS%kyIdxLower(rc)+1
    if (.not. noNkyaRanksS%isConjugate(rc)) then
      ![ The sends corresponding to positive-k amplitudes can be done
      ![ wth the MPI derived datatype 'subarray'.
      dims     = 2
      sizes    = [Nxa(2), Nxa(1)]
      subsizes = [Nxa(2), cpBlocks] 
      starts   = [0, noNkyaRanksS%kyIdxLower(rc)-1]    ![ starts are zero-indexed.
      call MPI_TYPE_CREATE_SUBARRAY(dims, sizes, subsizes, starts, &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, noNkyaRanksS%mpiDDtype(rc), errFlagMPI)
      call MPI_TYPE_COMMIT(noNkyaRanksS%mpiDDtype(rc),errFlagMPI)
    endif
  enddo

  deallocate(cpRankIDsTMP,cpBoundsTMP,isConjugateTMP)

  deallocate(lockyaAll,lockya0All,NekyaAll)

  ![ After the data is sent and received, some processed need to perform
  ![ a complex conjugation. Identify them before hand and save the bounds
  ![ of the data that needs to be conjugated.
  doConjugation = .false.
  if ((NkyaG <= lockya0+1) .and. (NkyaG < lockya0+Nxa(1))) then
    doConjugation = .true.
    cpIdxLower    =  9999
    cpIdxUpper    = -9999
    do lc = 1,Nxa(1)
      if ((lockya0+lc > NkyaG) .and. (cpIdxLower > lockya0+lc)) cpIdxLower = lc
      if ((lockya0+lc > NkyaG) .and. (cpIdxUpper < lockya0+lc)) cpIdxUpper = lc
    enddo
  endif

  end subroutine init_appendRemoveConjugates
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine FFT2Da_c2r(fkaIn, fRaOut)
  ![ Complex-to-real 2D FFT of aliased data.
  implicit none
  double complex, intent(in)    :: fkaIn(:,:)
  double precision, intent(out) :: fRaOut(:,:)

  ![ Relocate data from (Nekya x Nekxa) array field array to
  ![ (Nekya x Nekxa)=(Nekya x Nxa) buffer array.
  funka = fkaIn

  call fftw_mpi_execute_dft(plan1Da_c2r_x, funka, funkya)    ![ Inverse transform along kx.

  ![ Transpose the data locally. No global transpose is needed because
  ![ the second c2r transform is performed with a different communicator.
  forall(lcx=1:locxa,lcky=1:Nekxa(1)) funkyaT(lcx,lcky) = funkya(lcky,lcx)

  ![ Enlarge the data (funkyaT -> funkyaTc) by appending the complex conjugates along ky.
  do lc = 1,cpRanksRlen
    ![ Post a receive for each process that this rank expects to receive data from.
    call MPI_IRECV(funkyaTc,1,cpRanksR%mpiDDtype(lc),cpRanksR%rankID(lc),1,xyCOMMs(1),cpRanksR%req(lc),errFlagMPI)
  enddo
  do lc = 1,cpRanksSlen
    call MPI_ISEND(funkyaT,1,cpRanksS%mpiDDtype(lc),cpRanksS%rankID(lc),1,xyCOMMs(1),cpRanksS%req(lc),errFlagMPI)
  enddo
  call MPI_WAITALL(cpRanksRlen,cpRanksR%req,cpRanksR%stat,errFlagMPI)
  ![ Data was received as it was for positive-k modes. Complex conjugate it
  ![ for negative-k modes.
  if (doConjugation) funkyaTc(:,cpIdxLower:cpIdxUpper) = dconjg(funkyaTc(:,cpIdxLower:cpIdxUpper))
  call MPI_WAITALL(cpRanksSlen,cpRanksS%req,cpRanksS%stat,errFlagMPI)
 
  call fftw_mpi_execute_dft(plan1Da_c2r_y, funkyaTc, funRa)    ![ Inverse transform along ky.

  fRaOut = dble(funRa)*rSqrtNxyaG    ![ Unitary normalization.

  end subroutine FFT2Da_c2r
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine FFT2Da_r2c(fRaIn, fkaOut)
  ![ Real-to-complex 2D FFT of aliased data.
  implicit none
  double precision, intent(in) :: fRaIn(:,:)
  double complex, intent(out)  :: fkaOut(:,:)
  
  funRa = dcmplx(fRaIn,0.0d0)

  call fftw_mpi_execute_dft(plan1Da_r2c_y, funRa, funkyaTc)    ![ Forward transform along y.

  ![ Shrink the data (funkyaTc -> funkyaT) by disregarding the negative ky's (assume Hermitian symmetry).
  do lc = 1,cpRanksSlenP
    ![ Post a receive for each process that this rank expects to receive data from.
    call MPI_IRECV(funkyaT,1,noNkyaRanksR%mpiDDtype(lc),noNkyaRanksR%rankID(lc), &
                   1,xyCOMMs(1),noNkyaRanksR%req(lc),errFlagMPI)
  enddo
  do lc = 1,cpRanksRlenP
    call MPI_ISEND(funkyaTc,1,noNkyaRanksS%mpiDDtype(lc),noNkyaRanksS%rankID(lc), &
                   1,xyCOMMs(1),noNkyaRanksS%req(lc),errFlagMPI)
  enddo

  call MPI_WAITALL(cpRanksSlenP,noNkyaRanksR%req,noNkyaRanksR%stat,errFlagMPI)
  ![ Transpose the data locally. No global transpose is needed because
  ![ the second c2r transform is performed with a different communicator.
  forall(lcky=1:Nekxa(1),lcx=1:locxa) funkya(lcky,lcx) = funkyaT(lcx,lcky)
  call MPI_WAITALL(cpRanksRlenP,noNkyaRanksS%req,noNkyaRanksS%stat,errFlagMPI)

  call fftw_mpi_execute_dft(plan1Da_r2c_x, funkya, funka)    ![ Forward transform along x.

  fkaOut = funka*rSqrtNxyaG    ![ Unitary normalization.

  end subroutine FFT2Da_r2c
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  subroutine FFT2D_c2r(fkIn, fROut)
!  ![ Complex-to-real 2D FFT of de-aliased data.
!  implicit none
!  double complex, intent(in)    :: fkIn(:,:)
!  double precision, intent(out) :: fROut(:,:)
!
!  ![ Relocate data from (Neky x Nekx) array field array to
!  ![ (Neky x Nekx)=(Neky x Nx) buffer array.
!  funk = fkIn
!
!  call fftw_mpi_execute_dft(plan1D_c2r_x, funk, funky)    ![ Inverse transform along kx.
!
!  ![ Transpose the data locally. No global transpose is needed because
!  ![ the second c2r transform is performed with a different communicator.
!  forall(lcx=1:locx,lcky=1:Nekx(1)) funkyT(lcx,lcky) = funky(lcky,lcx)
!
!  ![ Enlarge the data (funkyT -> funkyTc) by appending the complex conjugates along ky.
!
!
!!  ![ Transpose the data, and use Hermitian symmetry to fill the negative ky's.
!!  forall(lcx=1:NxaF,lcky=1:Nekya)        funkyaT(lcx,lcky) = funkya(lcky,lcx)
!!  forall(lcx=1:NxaF,lcky=Nekya+1:NekyaF) funkyaT(lcx,lcky) = dconjg(funkya(Nekya-(lcky-Nekya),lcx))
!
!  call fftw_mpi_execute_dft(plan1Da_c2r_y, funkyTc, funR)    ![ Inverse transform along ky.
!
!  fROut = dble(funR)*rSqrtNxyG    ![ Unitary normalization.
!
!  end subroutine FFT2D_c2r
!![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  subroutine FFT2D_r2c(fRIn, fkOut)
!  ![ Real-to-complex 2D FFT of de-aliased data.
!  implicit none
!  double precision, intent(in) :: fRIn(:,:)
!  double complex, intent(out)  :: fkOut(:,:)
!  
!  funR = dcmplx(fRIn,0.0d0)
!
!  call fftw_mpi_execute_dft(plan1D_r2c_y, funR, funkyTc)    ![ Forward transform along y.
!
!!  ![ Transpose the data, and disregard the negative ky's (assume Hermitian symmetry).
!!  forall(lcx=1:NxF,lcky=1:Neky) funky(lcky,lcx) = funkyT(lcx,lcky)
!
!  ![ Shrink the data (funkyTc -> funkyT) by disregarding the negative ky's (assume Hermitian symmetry).
!
!
!  ![ Transpose the data locally. No global transpose is needed because
!  ![ the second c2r transform is performed with a different communicator.
!  forall(lcky=1:Nekx(1),lcx=1:locx) funky(lcky,lcx) = funkyT(lcx,lcky)
!
!  call fftw_mpi_execute_dft(plan1D_r2c_x, funky, funk)    ![ Forward transform along x.
!
!  fkOut = funk*rSqrtNxyG    ![ Unitary normalization.
!
!  end subroutine FFT2D_r2c

![ endif for if (FFTmode==0) then else.
#endif 

#else
![ else for if (useMPI>0) then else. 

  subroutine init_FFTs
  ![ Initialize plans and dynamically allocate FFTW-related arrays.
  ![ Unfortunately FFTW does not support 2D-domain-decomposed 2D FFTs,
  ![ so we have to use 1D FFTs. Furthermore, it does not have a parallel
  ![ real-to-complex 1D FFT, so we have to use complex-to-complex 1D
  ![ parallel FFTs.
  implicit none

  ![ ~~~~ Arrays and plans for aliased FFTs ~~~~ ]!
  ![ Dynamically allocate pointer to real and complex arrays
  funka_ptr = fftw_alloc_complex(int(NekyaG*NekxaG,C_SIZE_T))
  funRa_ptr = fftw_alloc_real(int(NyaG*NxaG,C_SIZE_T))
  ![ Convert C pointers to FORTRAN pointers
  call c_f_pointer(funRa_ptr,funRa,[NyaG,NxaG])
  call c_f_pointer(funka_ptr,funka,[NekyaG,NekxaG])
  ![ The order of the sizes has to be reversed here.
  plan2Da_c2r = fftw_plan_dft_c2r_2d(NxaG, NyaG, funka, funRa, FFTW_EXHAUSTIVE)
  plan2Da_r2c = fftw_plan_dft_r2c_2d(NxaG, NyaG, funRa, funka, FFTW_EXHAUSTIVE)


  if (initialOp == 0) then
    ![ ~~~~ Arrays and plans for de-aliased FFTs ~~~~ ]!
    ![ Dynamically allocate pointer to real and complex arrays
    funR_ptr = fftw_alloc_real(int(NyG*NxG,C_SIZE_T))
    funk_ptr = fftw_alloc_complex(int(NekyG*NekxG,C_SIZE_T))
    ![ Convert C pointers to FORTRAN pointers
    call c_f_pointer(funR_ptr,funR,[NyG,NxG])
    call c_f_pointer(funk_ptr,funk,[NekyG,NekxG])
    ![ The order of the sizes has to be reversed here.
    plan2D_c2r = fftw_plan_dft_c2r_2d(NxG, NyG, funk, funR, FFTW_EXHAUSTIVE)
    plan2D_r2c = fftw_plan_dft_r2c_2d(NxG, NyG, funR, funk, FFTW_EXHAUSTIVE)
  endif

  end subroutine init_FFTs
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine FFT2Da_c2r(fkaIn, fRaOut)
  ![ Complex-to-real 2D FFT of aliased data.
  implicit none
  double complex, intent(in)    :: fkaIn(:,:)
  double precision, intent(out) :: fRaOut(:,:)

  ![ Relocate data from (Nekya x Nekxa) array field array to
  ![ (Nekya x Nekxa)=(Nekya x Nxa) buffer array.
  funka = fkaIn

  call fftw_execute_dft_c2r(plan2Da_c2r, funka, funRa)

  fRaOut = funRa*rSqrtNxyaG    ![ Unitary normalization.

  end subroutine FFT2Da_c2r
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine FFT2Da_r2c(fRaIn, fkaOut)
  ![ Real-to-complex 2D FFT of aliased data.
  implicit none
  double precision, intent(in) :: fRaIn(:,:)
  double complex, intent(out)  :: fkaOut(:,:)
  
  funRa = fRaIn

  call fftw_execute_dft_r2c(plan2Da_r2c, funRa, funka)    ![ To be normalized.

  fkaOut = funka*rSqrtNxyaG    ![ Unitary normalization.

  end subroutine FFT2Da_r2c
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine FFT2D_c2r(fkIn, fROut)
  ![ Complex-to-real 2D FFT of de-aliased data.
  implicit none
  double complex, intent(in)    :: fkIn(:,:)
  double precision, intent(out) :: fROut(:,:)

  ![ Need to shift arrays so they are ordered as FFTW expects them.
  forall(lckx=1:NkxG,lcky=1:NekyG)      funk(lcky,lckx) = fkIn(lcky,lckx+NkxG-1)
  forall(lckx=NkxG+1:NekxG,lcky=1:NekyG) funk(lcky,lckx) = fkIn(lcky,lckx-NkxG)

  call fftw_execute_dft_c2r(plan2D_c2r, funk, funR)

  forall(lcx=1:NxG,lcy=1:NyG) fROut(lcy,lcx) = funR(lcy,lcx)*rSqrtNxyG    ![ Unitary normalization.

  end subroutine FFT2D_c2r
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine FFT2D_r2c(fRIn, fkOut)
  ![ Real-to-complex 2D FFT of de-aliased data.
  implicit none
  double precision, intent(in) :: fRIn(:,:)
  double complex, intent(out)  :: fkOut(:,:)
  
  forall(lcx=1:NxG,lcy=1:NyG) funR(lcy,lcx) = fRIn(lcy,lcx)

  call fftw_execute_dft_r2c(plan2D_r2c, funR, funk)

  forall(lckx=1:NekxG,lcky=1:NekyG)  fkOut(lcky,lckx) = funk(lcky,lckx)*rSqrtNxyG    ![ Unitary normalization.

  end subroutine FFT2D_r2c
#endif
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine terminate_FFTs
  ![ Destroy plans and free memory allocated for FFTW.
  implicit none

  ![ Destroy FFTW plans.
#if (useMPI > 0)
#if (FFTmode == 0)
  if (fftMember) then
    call fftw_destroy_plan(plan2Da_c2r)
    call fftw_destroy_plan(plan2Da_r2c)
    if (initialOp == 0) then
      call fftw_destroy_plan(plan2D_c2r)
      call fftw_destroy_plan(plan2D_r2c)
    endif
  endif
#else
  call fftw_destroy_plan(plan1Da_c2r_x)
  call fftw_destroy_plan(plan1Da_c2r_y)
  call fftw_destroy_plan(plan1Da_r2c_x)
  call fftw_destroy_plan(plan1Da_r2c_y)
  if (initialOp == 0) then
    call fftw_destroy_plan(plan1D_c2r_x)
    call fftw_destroy_plan(plan1D_c2r_y)
    call fftw_destroy_plan(plan1D_r2c_x)
    call fftw_destroy_plan(plan1D_r2c_y)
  endif
#endif
#else
  call fftw_destroy_plan(plan2Da_c2r)
  call fftw_destroy_plan(plan2Da_r2c)
  call fftw_destroy_plan(plan2D_c2r)
  call fftw_destroy_plan(plan2D_r2c)
#endif

  ![ Free dynamically allocated memory for FFTs and nullify other pointers.
  if (fftMember) then
    call fftw_free(funka_ptr)
    call fftw_free(funRa_ptr)
    nullify(funka,funRa)
    if (initialOp == 0) then
      call fftw_free(funk_ptr)
      call fftw_free(funR_ptr)
      nullify(funk,funR)
    endif
  endif

#if (useMPI > 0)
  ![ Arrays used for transfering de-aliased data to aliased arrays.
  if (fftMember) then
    do lc = 1,d2aSubarraysS
      call MPI_TYPE_free(d2aRanksS%mpiDDtype(lc),errFlagMPI)
    enddo
    deallocate(d2aRanksS%rankID,d2aRanksS%idxLower,d2aRanksS%idxUpper, &
               d2aRanksS%mpiDDtype,d2aRanksS%req,d2aRanksS%stat)
    do lc = 1,d2aSubarraysR
      call MPI_TYPE_free(d2aRanksR%mpiDDtype(lc),errFlagMPI)
    enddo
    deallocate(d2aRanksR%rankID,d2aRanksR%idxLower,d2aRanksR%idxUpper, &
               d2aRanksR%mpiDDtype,d2aRanksR%req,d2aRanksR%stat)
  endif
#if (FFTmode == 0)
  deallocate(funkyG)
  if (fftMember) then
    deallocate(kyGatherType%rankID, kyGatherType%mpiDDtype, &
               kyGatherType%kyIdxLower, kyGatherType%kyIdxUpper, &
               kyGatherType%req, kyGatherType%stat)
  endif
#else
  call fftw_free(funkya_ptr)
  call fftw_free(funkyaTc_ptr)
  nullify(funkya,funkyaTc)
  deallocate(funkyaT)
  ![ Arrays used to organize the communications in aliasing and de-aliasing.
  do lc = 1,3
    deallocate(d2aRanksS%remap(lc)%is)
    deallocate(d2aRanksR%remap(lc)%is)
  enddo
  ![ Arrays used for appending conjugate pairs (along ky).
  do lc = 1,cpRanksSlen
    call MPI_TYPE_free(cpRanksS%mpiDDtype(lc),errFlagMPI)
  enddo
  deallocate(cpRanksS%rankID,cpRanksS%kyIdxLower,cpRanksS%kyIdxUpper, &
             cpRanksS%isConjugate,cpRanksS%mpiDDtype,cpRanksS%req,cpRanksS%stat)
  do lc = 1,cpRanksRlen
    call MPI_TYPE_free(cpRanksR%mpiDDtype(lc),errFlagMPI)
  enddo
  deallocate(cpRanksR%rankID,cpRanksR%kyIdxLower,cpRanksR%kyIdxUpper, &
             cpRanksR%isConjugate,cpRanksR%mpiDDtype,cpRanksR%req,cpRanksR%stat)
  do lc = 1,cpRanksSlenP
    call MPI_TYPE_free(noNkyaRanksR%mpiDDtype(lc),errFlagMPI)
  enddo
  deallocate(noNkyaRanksR%rankID,noNkyaRanksR%kyIdxLower,noNkyaRanksR%kyIdxUpper, &
             noNkyaRanksR%isConjugate,noNkyaRanksR%mpiDDtype,noNkyaRanksR%req,noNkyaRanksR%stat)
  do lc = 1,cpRanksRlenP
    call MPI_TYPE_free(noNkyaRanksS%mpiDDtype(lc),errFlagMPI)
  enddo
  deallocate(noNkyaRanksS%rankID,noNkyaRanksS%kyIdxLower,noNkyaRanksS%kyIdxUpper, &
             noNkyaRanksS%isConjugate,noNkyaRanksS%mpiDDtype,noNkyaRanksS%req,noNkyaRanksS%stat)
  if (initialOp == 0) then
    call fftw_free(funky_ptr)
    call fftw_free(funkyTc_ptr)
    nullify(funky,funkyTc)
    deallocate(funkyT)
  endif
#endif

  call fftw_mpi_cleanup()
#endif

  end subroutine terminate_FFTs
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine testFFTs
  ![ Test FFTs by transforming an aliased real space field into aliased Fourier
  ![ space, and inverse transforming it back. The result should be zero change,
  ![ to machine precision.
  use IOtools
  implicit none
  double precision, allocatable :: xaG(:), yaG(:), xa(:), ya(:)
  double precision, allocatable :: flda(:,:)
  double complex, allocatable   :: fldka(:,:)

  allocate( xaG(NxaG), yaG(NyaG) )
  forall(lc=1:NxaG) xaG(lc) = dble(lc-1)*dxa+dble(1-mod(NxaG,2))*0.50d0*dxa-0.50d0*Lx
  forall(lc=1:NyaG) yaG(lc) = dble(lc-1)*dya+dble(1-mod(NyaG,2))*0.50d0*dya-0.50d0*Ly

  allocate(xa(Nxa(2)))
  xa  = xaG(firstxaG(2):firstxaG(2)+Nxa(2)-1)

  call alloc_2DFourierAliased(fldka)

#if ((useMPI > 0) && (FFTmode == 0))
  allocate(ya(NyaG))
  ya  = yaG

  allocate(flda(NyaG,Nxa(2)))
  forall(lcx=1:Nxa(2),lcy=1:NyaG) flda(lcy,lcx) = cos(2.0d0*(2.0d0*pi/Lx)*xa(lcx))*sin(3.0d0*(2.0d0*pi/Ly)*ya(lcy))

  if (fftMember) call out2DFieldRealAliased('phi',flda)

  call FFT2Da_r2c(flda, fldka)
  call FFT2Da_c2r(fldka, flda) 

  if (fftMember) call out2DFieldRealAliased('vort',flda)
#else
  allocate(ya(Nxa(1)))
  ya  = yaG(firstxaG(1):firstxaG(1)+Nxa(1)-1)

  allocate(flda(Nxa(2),Nxa(1)))
  forall(lcx=1:Nxa(2),lcy=1:Nxa(1)) flda(lcx,lcy) = cos(2.0d0*(2.0d0*pi/Lx)*xa(lcx))*sin(3.0d0*(2.0d0*pi/Ly)*ya(lcy))

  call out2DFieldRealAliased('phi',flda)

  call FFT2Da_r2c(flda, fldka)
  call FFT2Da_c2r(fldka, flda) 

  call out2DFieldRealAliased('vort',flda)
#endif

  deallocate(flda,fldka)
  deallocate(xa, ya, xaG, yaG)
  end subroutine testFFTs

end module FFTmodule

![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

module timeUpdates
  ![ Routines use for adaptive methods called on the fly.
  use parameters
  use fields

  contains

  subroutine compute_steppingFactors
  ![ This subroutine calculates the prefactors in front of f and df/dt in
  ![ the time-steppers. If hyperdiffusion is adjusted in time this gets
  ![ called every time step
  implicit none

  ![ Factors needed for implicit hyperdiffusion.
  do lckx = 1, Nekx(2); do lcky = 1, Nekx(1)
    denElcDiffFac(lcky,lckx)     = 1.0d0/(1.0d0-iHDdenElc(lcky,lckx))
    tempElcDiffFac(lcky,lckx)    = 1.0d0/(1.0d0-iHDtempElc(lcky,lckx))
    denIonDiffFac(lcky,lckx)     = 1.0d0/(1.0d0-iHDdenIon(lcky,lckx))
    tempIonDiffFac(lcky,lckx)    = 1.0d0/(1.0d0-iHDtempIon(lcky,lckx))
#if (HDtimeOP == 1)
    denElcDotDiffFac(lcky,lckx)  = denElcDiffFac(lcky,lckx)
    tempElcDotDiffFac(lcky,lckx) = tempElcDiffFac(lcky,lckx)
    denIonDotDiffFac(lcky,lckx)  = denIonDiffFac(lcky,lckx)
    tempIonDotDiffFac(lcky,lckx) = tempIonDiffFac(lcky,lckx)
#elif (HDtimeOP == 2)
    denElcDotDiffFac(lcky,lckx)  = 1.0d0
    tempElcDotDiffFac(lcky,lckx) = 1.0d0
    denIonDotDiffFac(lcky,lckx)  = 1.0d0
    tempIonDotDiffFac(lcky,lckx) = 1.0d0
#endif
  enddo; enddo

  end subroutine compute_steppingFactors
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
#if (ADJUSTdt > 0)
  function dtUpdate(phikIn, dtIn) result(newDt)
  ![ Check CFL condition (in aliased real space for now) and adjust the time step.
  ![ Only considers advection CFL constraint at the moment. There can also
  ![ be a diffusive limit if explicit diffusion is used.
  implicit none
  double precision             :: newDt
  double precision, intent(in) :: dtIn
  double complex, intent(in)   :: phikIn
  double precision, parameter  :: tol=0.50d0, chng=0.50d0
  double precision, parameter  :: CFL=0.20d0
  double precision             :: omegaMax, omegaMaxG, minDt, propDt
  integer                      :: specI
#if (adiabaticElc != 0)
  integer, parameter           :: specInit=1, specEnd=1
#elif (adiabaticIon != 0)
  integer, parameter           :: specInit=2, specEnd=2
#else
  integer, parameter           :: specInit=1, specEnd=2
#endif

  omegaMax = 0.0d0
  do specInd = specInit,specEnd 
    if (specInd == 1) then
      ![ Compute the Courant constraint from the ion ExB advection.
      ionGyroPhik = avgJ0Ion*phikIn
      call placeInAliased(ionGyroPhik,gka)
    endif
    
    if (specInd == 2) then
      ![ Compute the Courant constraint from the electron ExB advection.
      elcGyroPhik = avgJ0Elc*phikIn
      call placeInAliased(elcGyroPhik,gka)
    endif
  
    do lckx = 1,Nekxa(2)
      do lcky = 1,Nekxa(1)
        gka_x(lcky, lckx) = ikxa(lckx)*gka(lcky,lckx)
        gka_y(lcky, lckx) = ikya(lcky)*gka(lcky,lckx)
      enddo
    enddo
  
    call FFT2Da_c2r(gka_x, gRa_x) 
    call FFT2Da_c2r(gka_y, gRa_y) 
  
    do lcx=1,Nxa(2); do lcy=1,Nxa(1)
      omegaMax = max(omegaMax, abs(gRa_y(lcx,lcy))/dxa, abs(gRa_x(lcx,lcy))/dya)
    enddo; enddo
  enddo

  ![ Compute the maximum frequency over the whole domain.
  call MPI_Allreduce(omegaMax,omegaMaxG,1,MPI_DOUBLE_PRECISION,MPI_MAX,cartCOMM,errFlagMPI)

  omegaMax = max(omegaMaxG, CFL/dtMax)

  minDt = min(dtMax, CFL/omegaMax)

!  if (dtIn > minDt) then
    newDt = minDt
!    write(*,'("    > Time step lowered to dt= ",E14.4E3," <    ")') newDt
!  elseif (dtIn < minDt*(1.0d0-tol)) then
!    propDt = dtIn*(1.0d0+chng)
!    if (propDt < minDt) newDt = propDt
!    write(*,'("    > Time step raised to dt = ",E14.4E3," <    ")') newDt
!  endif

  end function dtUpdate
#endif
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
#if (((HDmodel > 2) && (HDmodel < 7)) || (HDmodel > 8))
  subroutine adjustHyperDiff(phikIn,vortkIn)
  ![ Adjust the order and/or the amplitude of the hyperdiffusion.
#if (useMPI > 0)
  use MPItools
#endif
  implicit none
  double complex, intent(in) :: phikIn(:,:), vortkIn(:,:)
#if (((HDmodel > 2) && (HDmodel < 7)) || (HDmodel > 8))
  double precision :: hDiffS, hDiffrkav     ![ Volume average shearing rate and reciprocal of average wavenumber.
  double precision :: hDiffSL, hDiffrkavL   ![ Local sum.
#endif
  double precision :: hDm,iHDfunc

  ![ Compute Smith+Hammett 1996 hyperviscosity parameters.
  hDiffSL    = 0.0d0
  hDiffrkavL = 0.0d0
  do lckx = 1,Nekx(2)
    do lcky = 1,Nekx(1)
      hDiffSL    = hDiffSL+abs(vortkIn(lcky,lckx))**2
      hDiffrkavL = hDiffrkavL+kSq(lcky,lckx)*(abs(phikIn(lcky,lckx))**2)
    enddo
  enddo
#if (useMPI > 0)
  call MPI_Allreduce(hDiffSL,hDiffS,1,MPI_DOUBLE_PRECISION,MPI_SUM,cartCOMM,errFlagMPI)
  call MPI_Allreduce(hDiffrkavL,hDiffrkav,1,MPI_DOUBLE_PRECISION,MPI_SUM,cartCOMM,errFlagMPI)
#else
  hDiffS    = hDiffSL
  hDiffrkav = hDiffrkavL
#endif
  hDiffS    = sqrt(0.50d0*hDiffS)
  hDiffrkav = sqrt(hDiffrkav)/hDiffS

#if (((HDmodel > 4) && (HDmodel < 7)) || (HDmodel > 10))
  ![ Change coefficient and order.
  hDiffOrder = kc1p7*hDiffrkav+2.40d0
  hDiffne    = kc0p1*hDiffS*hDiffrkav
#else
  ![ Change coefficient.
  hDiffne    = hDiffS !/(kc**hDiffOrder)
#endif
  hDiffTe    = hDiffne 
  hDiffni    = hDiffne 
  hDiffTi    = hDiffne 

  hDm = hDiffOrder/2.0d0
  do lckx = 1, Nekx(2)
    do lcky = 1, Nekx(1)

      ![ HDmodel=1-6 assume diffusion on vorticity (nabla^2 phi) not to phi.
      ![ HDmodel=7-12 assume diffusion on phi.
#if ((HDmodel == 3) || (HDmodel == 5))
      iHDfunc = -kDkcSq(lcky,lckx)**(hDm+1.0d0)
#elif ((HDmodel == 4) || (HDmodel == 6))
      iHDfunc = -(kxdcSq(lckx)**hDm+kydcSq(lcky)**hDm)*kSq(lcky,lckx)
#elif ((HDmodel == 9) || (HDmodel == 11))
      iHDfunc = -kDkcSq(lcky,lckx)**hDm
#elif ((HDmodel == 10) || (HDmodel == 12))
      iHDfunc = -(kxdcSq(lckx)**hDm+kydcSq(lcky)**hDm)
#endif

      iHDdenElc(lcky,lckx)  = hDiffne*iHDfunc
      iHDtempElc(lcky,lckx) = hDiffTe*iHDfunc
      iHDdenIon(lcky,lckx)  = hDiffni*iHDfunc
      iHDtempIon(lcky,lckx) = hDiffTi*iHDfunc

    enddo
  enddo

  ![ Update fnFac and fDotFac, factors multiplying f and df/dt in time steppers.
  call compute_steppingFactors

  end subroutine adjustHyperDiff
#endif
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine getPhik(denElckIn,tempElckIn,denIonkIn,tempIonkIn)
  ![ Solve the field equation (e.g. Poisson, quasineutrality) to obtain
  ![ the potential in Fourier space (phik).
  implicit none
  double complex, intent(in) :: denElckIn(:,:), tempElckIn(:,:)
  double complex, intent(in) :: denIonkIn(:,:), tempIonkIn(:,:)

  ![ Compute the potential from Poisson's equation.
  do lckx = 1, Nekx(2)
    do lcky = 1, Nekx(1)
      phik(lcky,lckx) = (chiIon*(poiSbIon(lcky,lckx)*denIonkIn(lcky,lckx) &
        +(1.0d0/(tau*deltaPerpi*DbIon(lcky,lckx)))*hatLapAvgJ0D2Ion(lcky,lckx)*tempIonkIn(lcky,lckx)) & 
        -oMchiElc(lcky,lckx)*(poiSbElc(lcky,lckx)*denElckIn(lcky,lckx) &
        +(1.0d0/(deltaPerpe*DbElc(lcky,lckx)))*hatLapAvgJ0D2Elc(lcky,lckx)*tempElckIn(lcky,lckx)))*rPoiPhikFac(lcky,lckx)
    enddo
  enddo

  ![ Also need some gyroaverage potentials.
#if (adiabaticElc==0)
  elcGyroPhik         = avgJ0Elc*phik
  hatLapD2elcGyroPhik = hatLapAvgJ0D2Elc*phik
#endif
#if (adiabaticIon==0)
  ionGyroPhik         = avgJ0Ion*phik
  hatLapD2ionGyroPhik = hatLapAvgJ0D2Ion*phik
#endif
  
  end subroutine getPhik

end module timeUpdates
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

module initialize
  use parameters
  use fields
  use IOtools
#if (useMPI > 0)
  use MPItools
#endif
  use timeUpdates
  use FFTmodule

  contains 

  subroutine read_inputs
  ![ Read the input variables, either from the input file (isRestart=F)
  ![ or from the restart file (isRestart=T).
  use utilities
  implicit none
  integer, parameter :: ifIdx=1    ![ File indexing.

  ![ Zeroth rank reads from input file.
  if (myRank==0) then
    open(ifIdx, file=trim(inputFile), status='old')    ![ Open file.
    
    ![ Read in namelists.
    read(ifIdx, NML=spaceIn)
    read(ifIdx, NML=sourceIn)
    read(ifIdx, NML=icsIn)
  
    close(ifIdx)    ![ Close file.
  endif

#if (useMPI > 0)
  ![ Broadcast from zeroth rank to all others.
  call MPI_BCAST(NkxG,       1, MPI_INTEGER,          0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(NkyG,       1, MPI_INTEGER,          0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(kxMin,      1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(kyMin,      1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(dt,         1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(endTime,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(nFrames,    1, MPI_INTEGER,          0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(xyProcs,    2, MPI_INTEGER,          0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(Lnorm,      1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(omSte,      1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(omde,       1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(tau,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(mu,         1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(deltae,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(deltaPerpe, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(eta_e,      1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(deltai,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(deltaPerpi, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(eta_i,      1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(lambdaD,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(hDiffOrder, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(hDiffne,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(hDiffTe,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(hDiffni,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(hDiffTi,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(initialOp,  1, MPI_INTEGER,          0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(initAuxX,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(initAuxY,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)
  call MPI_BCAST(initA,      1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, errFlagMPI)

  nProcs = xyProcs(1)*xyProcs(2)

  if (nProcs /= nProcesses) &
    call abortSimulation('nProcs parameter in code differs from that in mpirun -np xxx')
#else

  xyProcs = [1, 1] ![ Make sure the number of processes along kx and ky is 1.

#endif

  nFramesD = dble(nFrames)    ![ Used in time loop.

  ![ Maximum allowed time step.
  dtMax = dt*10.0d0

  if (isRestart) then
    ![ Re-write some inputs using those from restart file.
    call readRestartInputs
  endif
  
  end subroutine read_inputs
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine init_grids
  implicit none
  integer :: uniformNt

  ![ Relationship between de-aliased grid, aliased grid
  ![ and real space grid.

  ![ Given user-input number of distinct dealised wavenumbers, Nkx,
  ![ the number of real-space cells is Nx = 3*(Nkx-1). We prefer this
  ![ to be a power of 2, so we may need to adjust Nkx.
  if (myrank==0) &
    write(*,'(" User requested  NkxG=",i5, " and NkyG=",i5, " distinct wavenumbers (absolute magnitude).")') NkxG, NkyG
  NxaG = closest_power_of_two(3*(NkxG-1))
  NyaG = closest_power_of_two(3*(NkyG-1))

  ![ Number of distinct aliased (absolute) wavenumbers.
  NkxaG = NxaG/2+1
  NkyaG = NyaG/2+1
  ![ Length of aliased arrays along kx and ky.
  NekxaG = NxaG
  NekyaG = NkyaG

  ![ Recompute the number of distinct de-aliased (absolute) wavenumbers.
  NkxG = 2*(NkxaG-1)/3+1
  NkyG = 2*(NkyaG-1)/3+1
  ![ Length of de-aliased arrays along kx and ky.
  NekxG = 2*(NkxG-1)+1
  NekyG = NkyG

  ![ Number of cells in de-aliased real-space.
  NxG = 2*(NkxG-1)+1
  NyG = 2*(NkyG-1)

  ![ Expected number of time steps if taking uniform steps of size dt.
  uniformNt = floor(endTime/dt)
  
  if (myrank==0) then
    write(*,*) ' '
    write(*,'(" Proceeding with : ")')
    write(*,'(" Number of distinct de-aliased absolute wavenumbers: NkxG   =",i6, "  |  NkyG   =",i6)') NkxG, NkyG
    write(*,'(" Length of de-aliased k-space arrays:                NekxG  =",i6, "  |  NekyG  =",i6)') NekxG, NekyG
    write(*,'(" Number of distinct aliased absolute wavenumbers:    NkxaG  =",i6, "  |  NkyaG  =",i6)') NkxaG, NkyaG
    write(*,'(" Length of aliased k-space arrays:                   NekxaG =",i6, "  |  NekyaG =",i6)') NekxaG, NekyaG
    write(*,'(" Number of aliased real space cells:                 NxaG   =",i6, "  |  NyaG   =",i6)') NxaG, NyaG
    write(*,'(" Number of de-aliased real space cells:              NxG    =",i6, "  |  NyG    =",i6)') NxG, NyG
  endif

  ![ The size of the k-space cells is just kmin.
  dkx = kxMin
  dky = kyMin

  ![ Size of real-space domains (aliased and de-aliased).
  Lx = 2.0d0*pi/kxMin
  Ly = 2.0d0*pi/kyMin

  ![ Length of aliased real-space cell.
  dxa = Lx/dble(NxaG)
  dya = Ly/dble(NyaG)
  if (myrank==0) then
    write(*,*) ' '
    write(*,'(" Minimum absolute magnitude of wavenumbers: kxMin =",F9.4,"        |    kyMin     =", F9.4)') kxMin, kyMin
    write(*,'(" Number of time steps and frames:           Nt    =",i8, "         |  nFrames     =",i6)') uniformNt, nFrames
    write(*,'(" Time step:                                 dt    =",E13.4E2)') dt

    write(*,*) ' '
    write(*,'(" Normalized scale length        :           Lnorm =",E13.4E2)') Lnorm
    write(*,'(" Drift frequency factor         :           omSte =",E13.4E2)') omSte
    write(*,'(" Diamagnetic frequency factor   :            omde =",E13.4E2)') omde
    write(*,'(" Ti/Te ratio                    :             tau =",E13.4E2)') tau
    write(*,'(" Sqrt of m_i/m_e                :              mu =",E13.4E2)') mu
    write(*,'(" Elc & ion (Tpar+Tperp)/T       :          deltae =",E13.4E2,"    |  deltai     =", E13.4E2)') deltae, deltai
    write(*,'(" Elc & ion Tperp/T              :      deltaPerpe =",E13.4E2,"    |  deltaPerpi =", E13.4E2)') deltaPerpe, deltaPerpi
    write(*,'(" Elc & ion L_n/L_T              :           eta_e =",E13.4E2,"    |  eta_i      =", E13.4E2)') eta_e, eta_i
    write(*,'(" Debye shielding factor         :        lambdaD  =",E13.4E2)') lambdaD
    write(*,'(" Diffusion order                :      hDiffOrder =",E13.4E2)') hDiffOrder
    write(*,'(" Diffusion coefficients         :         hDiffne =",E13.4E2,"    |   hDiffni   =", E13.4E2)') hDiffne, hDiffni
    write(*,'("                                          hDiffTe =",E13.4E2,"    |   hDiffTi   =", E13.4E2)') hDiffTe, hDiffTi
  endif

  ![ Length of aliased real-space cell.
  dx = Lx/dble(NxG-mod(NxG,2))
  dy = Ly/dble(NyG-mod(NyG,2))

  ![ Global aliased and de-aliased real-space grids. 
  allocate( xG(NxG), yG(NyG) )
  forall(lc=1:NxG)   xG(lc) = dble(lc-1)*dx+dble(1-mod(NxG,2))*0.50d0*dx-0.50d0*Lx
  forall(lc=1:NyG)   yG(lc) = dble(lc-1)*dy+dble(1-mod(NyG,2))*0.50d0*dy-0.50d0*Ly

  ![ Global aliased and de-aliased k-space grids
  allocate( kxG(NekxG), kyG(NekyG), kxaG(NekxaG), kyaG(NekyaG) )
  forall(lc=1:NkxG)         kxG(lc)  = (lc-1)*dkx
  forall(lc=NkxG+1:NekxG)   kxG(lc)  = -( NkxG-1-(lc-(NkxG+1)) )*dkx
  forall(lc=1:NekyG)        kyG(lc)  = (lc-1)*dky
  forall(lc=1:NkxaG)        kxaG(lc) = (lc-1)*dkx
  forall(lc=NkxaG+1:NekxaG) kxaG(lc) = -( NkxaG-2-(lc-(NkxaG+1)) )*dkx
  forall(lc=1:NekyaG)       kyaG(lc) = (lc-1)*dky

  ![ Useful prefactors for normalizing FFTs.
  rSqrtNxyG  = 1.0/sqrt(dble(NxG*NyG))
  rSqrtNxyaG = 1.0/sqrt(dble(NxaG*NyaG))

  timeStepMore = .true.    ![ Take more time steps?
  simTime      = 0.0d0     ![ Current simulation time.
  timeSteps    = 0         ![ Time steps completed.

  ![ Time rate at which to output, adjust hyperdiffusion and adjust time step.
  tRateOutput   = endTime/dble(nFrames)
  tRateAdjustHD = dt
  tRateAdjustDt = 1.0d3*dt

  framesOut = 0    ![ Frames written so far (not counting ICs).
  hdAdjusts = 0    ![ Adjustments to hyperdiffusion so far.
  dtAdjusts = 0    ![ Adjustments to time step so far.

  tTol = dt*1.0d-1    ![ Time tolerance used in time-loop.

  ![ Establish the MPI decomposition and the local grids.
  ![ Give every process at least the floor of the division of the
  ![ degrees of freedom by the number of processors. The distribute
  ![ the remaining elements one at a time.

  ![ Also store the index of the first element in kx and ky
  ![ that is allocated to this process.
  call distributeDOFs2D(xyIDs,[xyProcs(2),xyProcs(1)],[NekyG,NekxG],Nekx,firstkG)

  ![ Distribute real-space de-aliased degrees of freedom (recall transpose in FFT2D is local).
  call distributeDOFs2D(xyIDs,[xyProcs(2),xyProcs(1)],[NyG,NxG],Nx,firstxG)

  ![ Distribute aliased degrees of freedom.
  call distributeDOFs2D(xyIDs,[xyProcs(2),xyProcs(1)],[NekyaG,NekxaG],Nekxa,firstkaG)

  ![ Distribute real-space aliased degrees of freedom (recall transpose in FFT2Da).
  call distributeDOFs2D(xyIDs,[xyProcs(2),xyProcs(1)],[NyaG,NxaG],Nxa,firstxaG)

  ![ Allocate arrays for local wavenumbers. 
  call alloc_1DrealX( x )
  call alloc_1DrealY( y )
  call alloc_1DFourierX( kx )
  call alloc_1DFourierY( ky )
  call alloc_1DFourierAliasedX( kxa )
  call alloc_1DFourierAliasedY( kya )

  y   = yG(firstxG(1):firstxG(1)+Nx(1)-1)
  x   = xG(firstxG(2):firstxG(2)+Nx(2)-1)

  ky  = kyG(firstkG(1):firstkG(1)+Nekx(1)-1)
  kx  = kxG(firstkG(2):firstkG(2)+Nekx(2)-1)
#if ((useMPI > 0) && (FFTmode == 0))
  kya = kyaG
#else
  kya = kyaG(firstkaG(1):firstkaG(1)+Nekxa(1)-1)
#endif
  kxa = kxaG(firstkaG(2):firstkaG(2)+Nekxa(2)-1)

  call alloc_1DFourierAliasedX( ikxa )
  call alloc_1DFourierAliasedY( ikya )
  ikxa = Im1*kxa
  ikya = Im1*kya

  end subroutine init_grids
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine allocate_fields
  use utilities
  implicit none

  ![ Allocate field arrays.
  call alloc_2DFourier(phik)
  call alloc_2Dreal(phi)
  call alloc_2DFourier(elcGyroPhik)
  call alloc_2DFourier(ionGyroPhik)
  call alloc_2DFourier(hatLapD2elcGyroPhik)
  call alloc_2DFourier(hatLapD2ionGyroPhik)

  call alloc_2DFourier(denElck)
  call alloc_2DFourier(tempElck)
  call alloc_2DFourier(denIonk)
  call alloc_2DFourier(tempIonk)
  call alloc_2Dreal(denElc)
  call alloc_2Dreal(tempElc)
  call alloc_2Dreal(denIon)
  call alloc_2Dreal(tempIon)

  call alloc_2Dreal(vort)
  call alloc_2DFourier(vortk)

  call alloc_2DFourier(denElcDot)
  call alloc_2DFourier(tempElcDot)
  call alloc_2DFourier(denIonDot)
  call alloc_2DFourier(tempIonDot)
#if ((STEPPER == 3) || (STEPPER == 4) || (STEPPER == 5))
  call alloc_2DFourier(denElcDot1)
  call alloc_2DFourier(tempElcDot1)
  call alloc_2DFourier(denIonDot1)
  call alloc_2DFourier(tempIonDot1)
  call alloc_2DFourier(denElcDot2)
  call alloc_2DFourier(tempElcDot2)
  call alloc_2DFourier(denIonDot2)
  call alloc_2DFourier(tempIonDot2)
#if (STEPPER == 4)
  call alloc_2DFourier(denElcDot3)
  call alloc_2DFourier(tempElcDot3)
  call alloc_2DFourier(denIonDot3)
  call alloc_2DFourier(tempIonDot3)
#endif
#endif

  call alloc_2DFourier(PBk)
  call alloc_2DFourier(flrPBk)
#if (NONLINEAR == 0)
  PBk    = dcmplx(0.0d0, 0.0d0)
  flrPBk = dcmplx(0.0d0, 0.0d0)
#endif
  ![ Arrays defined on aliased Fourier grid.
  call alloc_2DFourierAliased(fka)
  call alloc_2DFourierAliased(gka)
  call alloc_2DFourierAliased(fka_x)
  call alloc_2DFourierAliased(gka_x)
  call alloc_2DFourierAliased(fka_y)
  call alloc_2DFourierAliased(gka_y)
  call alloc_2DFourierAliased(PBfgka)
  ![ Arrays defined on aliased real-space grid.
  call alloc_2DrealAliased(fRa_x)
  call alloc_2DrealAliased(gRa_x)
  call alloc_2DrealAliased(fRa_y)
  call alloc_2DrealAliased(gRa_y)
  call alloc_2DrealAliased(PBfgRa)

  ![ For multi-step integrators need solution at previous step.
#if (STEPPER == 2)
  call alloc_2DFourier(phiks)
  call alloc_2DFourier(denElcks)
  call alloc_2DFourier(tempElcks)
  call alloc_2DFourier(denIonks)
  call alloc_2DFourier(tempIonks)
#endif

  ![ Temporary field used by time steppers.
  call alloc_2DFourier(denElckTmp)
  call alloc_2DFourier(tempElckTmp)
  call alloc_2DFourier(denIonkTmp)
  call alloc_2DFourier(tempIonkTmp)

  ![ Implicit hyperdiffusion function (depends on hyperdiffusion model).
  call alloc_2DFourier( iHDdenElc )
  call alloc_2DFourier( iHDtempElc )
  call alloc_2DFourier( iHDdenIon )
  call alloc_2DFourier( iHDtempIon )
  ![ Factor applying implicit hyperdiffusion (e.g. 1/(1 - dt*hDiff*kSq)).
  call alloc_2DFourier( denElcDiffFac )
  call alloc_2DFourier( tempElcDiffFac )
  call alloc_2DFourier( denIonDiffFac )
  call alloc_2DFourier( tempIonDiffFac )
  call alloc_2DFourier( denElcDotDiffFac )
  call alloc_2DFourier( tempElcDotDiffFac )
  call alloc_2DFourier( denIonDotDiffFac )
  call alloc_2DFourier( tempIonDotDiffFac )

  end subroutine allocate_fields
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  function closest_power_of_two(intIn)    ![ Find the closest power of 2.
  implicit none
  integer, intent(in) :: intIn    ![ Input integer.
  integer             :: prev_power_of_two    ![ Power of 2 before intIn.
  integer             :: next_power_of_two    ![ Power of 2 after intIn.
  integer             :: closest_power_of_two 
  integer             :: lc                   ![ Dummy counter.
  integer             :: prevDist, nextDist   ![ Distance from intIn.

  lc = 0      
  do while (((2**lc) < intIn) .and. (lc < 1000000))    ![ Safety check included for big lc.
    lc = lc + 1
  end do
  prev_power_of_two = 2**(lc-1)
  next_power_of_two = 2**lc

  prevDist = abs(prev_power_of_two-intIn)
  nextDist = abs(next_power_of_two-intIn)
  if (prevDist < nextDist) then
    closest_power_of_two = prev_power_of_two
  else
    closest_power_of_two = next_power_of_two
  endif

  end function closest_power_of_two
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine init_preFactors
  ![ Pre-compute some arrays used in time evolution.
  use FLRfunctions
  implicit none
  double precision, allocatable :: bIon(:,:), bElc(:,:)
  double precision              :: kperp,krhoIon,krhoElc,GammaRat, &
                                   hatLap,hatLapSq,hathatLap,Nb
  double precision              :: hDm,iHDfunc

  ![ Simple functions of kx and ky.
  call alloc_1DFourierY( iky )
  call alloc_2DFourier( kSq )
  call alloc_2DFourier( bIon )
  call alloc_2DFourier( bElc )
  forall(lcky=1:Nekx(1)) iky(lcky) = Im1*ky(lcky)
  do lckx = 1, Nekx(2)
    do lcky = 1, Nekx(1)
      kSq(lcky,lckx) = kx(lckx)*kx(lckx) + ky(lcky)*ky(lcky)

      kperp   = sqrt(kSq(lcky,lckx))
      krhoIon = sqrt(tau)*kperp
      krhoElc = kperp/mu

      bIon(lcky,lckx) = krhoIon**2 
      bElc(lcky,lckx) = krhoElc**2 
    enddo
  enddo

  ![ Diamagnetic and curvature drift frequencies.
  call alloc_1DFourierY( iOmdElc )
  call alloc_1DFourierY( iOmStElc )
  call alloc_1DFourierY( iOmdIon )
  call alloc_1DFourierY( iOmStIon )
  forall(lcky=1:Nekx(1))
    iOmdElc(lcky)  = Im1*omde*ky(lcky)
    iOmStElc(lcky) = -Im1*omSte*ky(lcky)
    iOmdIon(lcky)  = -tau*iOmdElc(lcky)
    iOmStIon(lcky) = tau*iOmStElc(lcky)
  end forall

  ![ Quantities used for FLR effects.
  call alloc_2DFourier( Gamma0Ion )
  call alloc_2DFourier( Gamma0Elc )
  call alloc_2DFourier( avgJ0Elc )
  call alloc_2DFourier( avgJ0Ion )
  call alloc_2DFourier( hatLapAvgJ0D2Elc )
  call alloc_2DFourier( hatLapAvgJ0D2Ion )
  call alloc_2DFourier( poiSbElc )
  call alloc_2DFourier( poiSbIon )
  call alloc_2DFourier( DbElc )
  call alloc_2DFourier( DbIon )
  do lckx = 1, Nekx(2)
    do lcky = 1, Nekx(1)
      ![ Ion FLR functions.
      Gamma0Ion(lcky,lckx) = besselI0Exp(bIon(lcky,lckx))
      avgJ0Ion(lcky,lckx)  = sqrt(Gamma0Ion(lcky,lckx))

      GammaRat  = besselI1Exp(bIon(lcky,lckx))/Gamma0Ion(lcky,lckx)
      hatLap    = bIon(lcky,lckx)*(GammaRat-1.0d0)
      hatLapSq  = hatLap**2
      hathatLap = bIon(lcky,lckx)*((0.50d0*GammaRat-1.0d0)-0.25d0*bIon(lcky,lckx)*(3.0d0+GammaRat)*(GammaRat-1.0d0))
      Nb        = 1.0d0+hathatLap-0.50d0*hatLapSq
    
      hatLapAvgJ0D2Ion(lcky,lckx) = 0.50d0*hatLap*avgJ0Ion(lcky,lckx)
      DbIon(lcky,lckx)            = 1.0d0+hathatLap-0.25d0*hatLapSq
      poiSbIon(lcky,lckx)         = (avgJ0Ion(lcky,lckx)*Nb)/DbIon(lcky,lckx)

      ![ Electron FLR functions.
      Gamma0Elc(lcky,lckx) = besselI0Exp(bElc(lcky,lckx))
      avgJ0Elc(lcky,lckx)  = sqrt(Gamma0Elc(lcky,lckx))

      GammaRat  = besselI1Exp(bElc(lcky,lckx))/Gamma0Elc(lcky,lckx)
      hatLap    = bElc(lcky,lckx)*(GammaRat-1.0d0)
      hatLapSq  = hatLap**2
      hathatLap = bElc(lcky,lckx)*((0.50d0*GammaRat-1.0d0)-0.25d0*bElc(lcky,lckx)*(3.0d0+GammaRat)*(GammaRat-1.0d0))
      Nb        = 1.0d0+hathatLap-0.50d0*hatLapSq
    
      hatLapAvgJ0D2Elc(lcky,lckx) = 0.50d0*hatLap*avgJ0Elc(lcky,lckx)
      DbElc(lcky,lckx)            = 1.0d0+hathatLap-0.25d0*hatLapSq
      poiSbElc(lcky,lckx)         = (avgJ0Elc(lcky,lckx)*Nb)/DbElc(lcky,lckx)
    enddo
  enddo
  call alloc_2DFourier( chiElc )
  call alloc_2DFourier( oMchiElc )
#if (adiabaticElc > 0)
  chiElc = 1.0d0
  chiIon = 1.0d0
#elif (adiabaticIon > 0)
  chiElc = 0.0d0
  chiIon = 0.0d0
#else
  chiElc = Gamma0Ion
  chiIon = 1.0
#endif
  oMchiElc = 1.0d0-chiElc

  ![ Reciprocal of the factor multiplying the potential in the field equation.
  call alloc_2DFourier( rPoiPhikFac )
  do lckx = 1, Nekx(2)
    do lcky = 1, Nekx(1)
      rPoiPhikFac(lcky,lckx) = (lambdaD**2)*kSq(lcky,lckx) &
                              +(1.0d0/(tau*Lnorm))*(1.0d0-chiIon*Gamma0Ion(lcky,lckx)) &
                              -oMchiElc(lcky,lckx)*(1.0d0/Lnorm)*(1.0d0-Gamma0Elc(lcky,lckx)) &
                              +(1.0d0/Lnorm)*chiElc(lcky,lckx)
      if ((xyIDs(1)==0) .and. (lcky == 1)) then
        ![ Subtract the zonal flow term.
        rPoiPhikFac(lcky,lckx) = rPoiPhikFac(lcky,lckx)-(1.0d0/Lnorm)*chiElc(lcky,lckx)
      endif
      if (abs(rPoiPhikFac(lcky,lckx)) .gt. epsilon(kx(lckx))) then 
        rPoiPhikFac(lcky,lckx) = 1.0d0/rPoiPhikFac(lcky,lckx)
      else
        rPoiPhikFac(lcky,lckx) = 0.0d0
      endif
    enddo
  enddo

  ![ For Smith+Hammett 1996 hyperviscosity need the following.
  call alloc_1DFourierY( kydcSq )
  call alloc_1DFourierX( kxdcSq )
  forall(lc=1:Nekx(1)) kydcSq(lc) = (abs(ky(lc))/kyG(NkyG))**2
  forall(lc=1:Nekx(2)) kxdcSq(lc) = (abs(kx(lc))/kxG(NkxG))**2
  kc    = sqrt(kxG(NkxG)**2+kyG(NkyG)**2)
  kc0p1 = 0.10d0*kc
  kc1p7 = 1.70d0*kc

  call alloc_2DFourier( kDkcSq )
  kDkcSq = kSq/(kc**2)

  hDm = hDiffOrder/2.0d0
  do lckx = 1, Nekx(2)
    do lcky = 1, Nekx(1)

      ![ HDmodel=1-6 assume diffusion on the Laplacian of the dynamic quantity (e.g. nabla^2 n_e).
      ![ HDmodel=7-12 assume diffusion on the dynamic variable itself.
#if (HDmodel == 0)
      iHDfunc = 0.0d0
#elif (HDmodel == 1)
      iHDfunc = -kDkcSq(lcky,lckx)**(hDm+1.0d0)
#elif (HDmodel == 2)
      iHDfunc = -(kxdcSq(lckx)**hDm+kydcSq(lcky)**hDm)*kSq(lcky,lckx)
#elif ((HDmodel == 3) || (HDmodel == 5))
      iHDfunc = -kDkcSq(lcky,lckx)**(hDm+1.0d0)
#elif ((HDmodel == 4) || (HDmodel == 6))
      iHDfunc = -(kxdcSq(lckx)**hDm+kydcSq(lcky)**hDm)*kSq(lcky,lckx)
#elif (HDmodel == 7)
      iHDfunc = -kDkcSq(lcky,lckx)**hDm
#elif (HDmodel == 8)
      iHDfunc = -(kxdcSq(lckx)**hDm+kydcSq(lcky)**hDm)
#elif ((HDmodel == 9) || (HDmodel == 11))
      iHDfunc = -kDkcSq(lcky,lckx)**hDm
#elif ((HDmodel == 10) || (HDmodel == 12))
      iHDfunc = -(kxdcSq(lckx)**hDm+kydcSq(lcky)**hDm)
#endif

#if ((HDmodel < 3) || (HDmodel == 7) || (HDmodel == 8))
      iHDdenElc(lcky,lckx)  = dt*hDiffne*iHDfunc
      iHDtempElc(lcky,lckx) = dt*hDiffTe*iHDfunc
      iHDdenIon(lcky,lckx)  = dt*hDiffni*iHDfunc
      iHDtempIon(lcky,lckx) = dt*hDiffTi*iHDfunc
#elif (((HDmodel > 2) && (HDmodel < 7)) || (HDmodel > 8))
      iHDdenElc(lcky,lckx)  = iHDfunc
      iHDtempElc(lcky,lckx) = iHDfunc
      iHDdenIon(lcky,lckx)  = iHDfunc
      iHDtempIon(lcky,lckx) = iHDfunc
#endif

    enddo
  enddo

  call compute_steppingFactors

  deallocate(bIon,bElc)

  end subroutine init_preFactors
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine init_random_seed()
  ![ Initialize a random seed based on system clock for random
  ![ initial fluctuations.
  implicit none
  integer              :: i, n, clock
  integer, allocatable :: seed(:)
  
  call random_seed(size = n)
  allocate(seed(n))
  
  call system_clock(count=clock)
  
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)
  
  deallocate(seed)
  end subroutine init_random_seed
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine setInitialCondition
  ![ Initialize the Fourier transform of the electrostatic potential. 
  use utilities
  implicit none

  if (isRestart) then
    ![ This simulation is a restart of a previous one.
    ![ Read fields from restart file.
    call readRestartFields

    if (abs(simTime-endTime) <= tTol) &
      call abortSimulation('simTime >= endTime; do nothing. Terminating...')

    if (myRank == 0) then
      write(*,'(" Restarting a simulation from")')
      write(*,*) "     ",adjustl(trim(restartFile))
      write(*,'(" at time t = ",E15.4E2)') simTime
      write(*,'( )')
    endif

  else
    ![ This simulation is not a restart. Apparate fields from nothing!
    call init_random_seed()

    
    if (initialOp == 0) then

      if (myRank == 0) then
        write(*,'(" Starting a simulation from time t = 0.0")')
        write(*,'(" Initial condition: random noise with amplitude", E14.4E3)') initA
        write(*,*) ' '
      endif

      ![ Just initial noise in real (dealiased) space.
      call random_number(phi)
      denElc = initA*phi
      denIon = initA*phi
      call random_number(phi)
      tempElc = initA*phi
      tempIon = initA*phi

      phi = 0.0d0  ![ Potential gets computed by the field solver.

      ![ Convert real-space to Fourier space.
!      call FFT2D_r2c(denElc, denElck) 
!      call FFT2D_r2c(denIon, denIonk) 
!      call FFT2D_r2c(tempElc, tempElck) 
!      call FFT2D_r2c(tempIon, tempIonk) 

    elseif (initialOp == 1) then

      if (myRank == 0) then
        write(*,'(" Starting a simulation from time t = 0.0")')
        write(*,'(" Initial condition: k-space power law with amplitude", E14.4E3)') initA
        write(*,*) ' '
      endif

      ![ Initialize densities with power-law in k-space.
      do lckx = 1,Nekx(2)
        do lcky = 1,Nekx(1)
          if ((abs(kx(lckx)) .lt. epsilon(kx(lckx))) .and. (abs(ky(lcky)) .lt. epsilon(ky(lcky)))) then
            denElck(lcky,lckx) = dcmplx(0.0d0,0.0d0)
          else
            denElck(lcky,lckx) = initA*(((kxMin+abs(kx(lckx)))/kxMin)**initAuxX)*(((kyMin+abs(ky(lcky)))/kyMin)**initAuxY)
          endif
        enddo
      enddo
      denIonk = denElck
      ![ Initialize temperatures to zero.
      tempElck = dcmplx(0.0d0,0.0d0)
      tempIonk = dcmplx(0.0d0,0.0d0)

    endif

#if (STEPPER == 2)
    ![ Densities and temperatures at previous time step.
    denElcks = denElck
    denIonks = denIonk
    tempElcks = tempElck
    tempIonks = tempIonk
#endif

  endif    ![ End isRestart if statement.

  ![ Compute the electrostatic potential.
  call getPhik(denElck,tempElck,denIonk,tempIonk)

  vortk = -kSq*phik    ![ Compute the vorticity in Fourier space.

#if (((HDmodel > 2) && (HDmodel < 7)) || (HDmodel > 8))
  call adjustHyperDiff(phik,vortk)
#endif

  end subroutine setInitialCondition
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine init_all
  ![ Initialize all variables, arrays and functions.
  use utilities
  implicit none
  integer :: args

  ![ Check for commandline arguments.
  ![ Currently we expect two arguments in this order:
  ![   1) name of the input file.
  ![   2) absolute address of the output directory.
  ![ For restarts add the restart directory as a 3rd argument.
  args = 0
  args = command_argument_count()
  if (args < 2) then 
    if (myrank==0) then
      write(*,'( )')
      write(*,'(" --> Not enough inputs. Need input file and output directory as command line arguments." )')
      write(*,'(" TERMINATING ")')
    endif
    call abortSimulation(' ')
  else
    call get_command_argument(1, inputFile)
    call get_command_argument(2, outputDir)
    isRestart = .false.
    if (args > 2) then    ![ Simulation is a restart of a previous one.
      call get_command_argument(3, restartDir)
      isRestart = .true.
    endif
  endif

  ![ Check that electrons and ions are not both adiabatic
  if ((adiabaticElc > 0) .and. (adiabaticIon > 0)) then
    write(*,'( )')
    write(*,'(" --> Cannot have adiabaticElc and adiabaticIon both >0." )')
    write(*,'(" TERMINATING ")')
    call abortSimulation(' ')
  endif

  call init_IO        ![ Initialize I/O interface.

  ![ Read inputs and initialize variables and arrays.
  call read_inputs

  fftMember = .true.  ![ Flag indicating if this MPI process participates in FFT.
#if (useMPI > 0)
  call init_COMMs  ![ Create subcommunicators and MPI types.
#else
  xyIDs    = [0, 0]    ![ No MPI. Only one rank/process.
  cartRank = 0
#endif

  call init_grids

  call init_FFTs

  call allocate_fields

  call init_preFactors

  call setupIOfiles    ![ Create or locate files for IO.

  ![ Initialize field arrays.
  call setInitialCondition

  ![ Set the time rate at which to enter messages in the log and
  ![ determine the number of entries in log files/screen so far.
  tRateLogEntry = tRateLogEntry*endTime
  if (isRestart) then
    logEntries  = nint(simTime/tRateLogEntry)
  else
    logEntries  = 0
  endif

  end subroutine init_all
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

end module initialize

![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

module terminate
  use parameters
  use fields
  use FFTmodule
  use IOtools

  contains

  subroutine isNaNorInf
  ![ Check if there is a NaN or an extremely large number (that may
  ![ become an Inf) in the main quantities.
  use utilities
  implicit none
  logical          :: infF, nanF
  double precision :: amp

  infF = .FALSE.
  nanF = .FALSE.
  do lckx=1,Nekx(2); do lcky=1,Nekx(1)
    amp = abs(phik(lcky,lckx))
    if (amp>1.0d16) infF = .TRUE.
    if (isnan(amp)) nanF = .TRUE.
  enddo; enddo
  if (infF) call abortSimulation('Value >1e16 found. Terminating simulation.')
  if (nanF) call abortSimulation('NaN value found. Terminating simulation.')

  end subroutine isNaNorInf
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine deallocate_parameters
  ![ Subroutine that deallocates all dynamically allocated arrays except for field arrays.
  implicit none
  deallocate(  xG,  yG )
  deallocate( kxG, kyG, kxaG, kyaG )
  deallocate(  x,  y )
  deallocate( kx, ky, kxa, kya )
  deallocate( ikxa, ikya )

  deallocate( kSq )
  deallocate( iky )

  deallocate( iOmdElc, iOmStElc )
  deallocate( iOmdIon, iOmStIon )

  deallocate( Gamma0Ion, Gamma0Elc )
  deallocate( avgJ0Elc, avgJ0Ion )
  deallocate( hatLapAvgJ0D2Elc, hatLapAvgJ0D2Ion )
  deallocate( poiSbElc, poiSbIon )
  deallocate( DbElc, DbIon )
  deallocate( chiElc, oMchiElc )

  deallocate( rPoiPhikFac )

  deallocate( kxdcSq, kydcSq )
  deallocate( kDkcSq )
  end subroutine deallocate_parameters
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine deallocate_fields
  ![ Subroutine that deallocates all dynamically allocated field arrays.
  implicit none
  deallocate( phik, phi )
  deallocate( elcGyroPhik, ionGyroPhik )
  deallocate( hatLapD2elcGyroPhik, hatLapD2ionGyroPhik )

  deallocate( denElck, tempElck )
  deallocate( denIonk, tempIonk )
  deallocate( denElc, tempElc )
  deallocate( denIon, tempIon )
  deallocate( vort, vortk )

  deallocate( denElcDot, tempElcDot )
  deallocate( denIonDot, tempIonDot )
#if ((STEPPER == 3) || (STEPPER == 4) || (STEPPER == 5))
  deallocate( denElcDot1, tempElcDot1 )
  deallocate( denIonDot1, tempIonDot1 )
  deallocate( denElcDot2, tempElcDot2 )
  deallocate( denIonDot2, tempIonDot2 )
#if (STEPPER == 4)
  deallocate( denElcDot3, tempElcDot3 )
  deallocate( denIonDot3, tempIonDot3 )
#endif
#endif

  deallocate( PBk, flrPBk )
  deallocate( fka, gka )
  deallocate( fka_x, gka_x, &
              fka_y, gka_y )
  deallocate( fRa_x, gRa_x, &
              fRa_y, gRa_y )
  deallocate( PBfgRa )
  deallocate( PBfgka )

#if (STEPPER == 2)
  deallocate( phiks )
  deallocate( denElcks, tempElcks )
  deallocate( denIonks, tempIonks )
#endif

  deallocate(denElckTmp, tempElckTmp)
  deallocate(denIonkTmp, tempIonkTmp)

  deallocate( iHDdenElc, iHDtempElc )
  deallocate( iHDdenIon, iHDtempIon )
  deallocate( denElcDiffFac, tempElcDiffFac )
  deallocate( denIonDiffFac, tempIonDiffFac )
  deallocate( denElcDotDiffFac, tempElcDotDiffFac )
  deallocate( denIonDotDiffFac, tempIonDotDiffFac )
  end subroutine deallocate_fields
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine terminate_all
  ![ Call routines to deallocate arrays, pointers and free memory.
  implicit none

  call terminate_IO    ![ Close files and output interface.
 
  call terminate_FFTs           ![ Free FFTW memory and deallocate pointers.
  call deallocate_parameters    ![ Deallocate other time-constant arrays.
  call deallocate_fields        ![ Deallocate field arrays.

  end subroutine terminate_all

end module terminate

![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

module nonlinearities
  use FFTmodule
  use fields

  contains

![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine placeInAliased(fkIn,fkaOut)
  ![ Take de-aliased data and place it in an aliased array.
  implicit none
  double complex, intent(in)    :: fkIn(:,:)
  double complex, intent(inout) :: fkaOut(:,:)
#if (useMPI > 0)
#if (FFTmode == 0)
  integer                       :: idx,sendCount
#else
  integer                       :: idx
#endif
#endif

#if (useMPI > 0)
#if (FFTmode == 0)
  ![ First perform a gather of the data along ky.
  if (xyIDs(1) == fftKyID) then
    do lc = 1, xyProcs(2)
      if (kyGatherType%rankID(lc) == xyIDs(1)) then
        ![ Copy instead of self-messaging.
        funkyG(kyGatherType%kyIdxLower(lc):kyGatherType%kyIdxUpper(lc),:) = fkIn
        kyGatherType%req(lc) = MPI_REQUEST_NULL
      else
        call MPI_Irecv(funkyG,1,kyGatherType%mpiDDtype(lc),kyGatherType%rankID(lc), &
                       1,xyCOMMs(1),kyGatherType%req(lc),errFlagMPI)
      endif
    enddo
    call MPI_Waitall(xyProcs(2),kyGatherType%req,kyGatherType%stat,errFlagMPI)
  else
    sendCount = Nekx(1)*Nekx(2)
    call MPI_Send(fkIn,sendCount,MPI_DOUBLE_COMPLEX,fftKyID,1,xyCOMMs(1),errFlagMPI)
  endif

  if (fftMember) then
    ![ Then place the data in a larger, de-aliased array which
    ![ is now distributed accross fftCOMM only.
    do lc = 1, d2aSubarraysR
      ![ Post a receive for each process that this rank expects to receive data from.
      if (d2aRanksR%rankID(lc) == fftRank) then
        ![ Copy instead of self-messaging.
        idx = d2aRanksS%selfIdx
        fkaOut(d2aRanksR%idxLower(lc,1):d2aRanksR%idxUpper(lc,1),d2aRanksR%idxLower(lc,2):d2aRanksR%idxUpper(lc,2)) &
          = funkyG(d2aRanksS%idxLower(idx,1):d2aRanksS%idxUpper(idx,1),d2aRanksS%idxLower(idx,2):d2aRanksS%idxUpper(idx,2))
        d2aRanksR%req(lc) = MPI_REQUEST_NULL
      else
        call MPI_Irecv(fkaOut,1,d2aRanksR%mpiDDtype(lc),d2aRanksR%rankID(lc), &
                       1,fftCOMM,d2aRanksR%req(lc),errFlagMPI)
      endif
    enddo
    do lc = 1, d2aSubarraysS
      if (d2aRanksS%rankID(lc) == fftRank) then
        ![ Do not send to self.
        d2aRanksS%req(lc) = MPI_REQUEST_NULL
      else
        call MPI_Isend(funkyG,1,d2aRanksS%mpiDDtype(lc),d2aRanksS%rankID(lc), &
                       1,fftCOMM,d2aRanksS%req(lc),errFlagMPI)
      endif
    enddo
  
    call MPI_Waitall(d2aSubarraysR,d2aRanksR%req,d2aRanksR%stat,errFlagMPI)
    call MPI_Waitall(d2aSubarraysS,d2aRanksS%req,d2aRanksS%stat,errFlagMPI)
  endif

#else
  ![ Placing the data in the larger, aliased arrays requires re-organizing
  ![ the data across the 2D communicator.
  ![ NOTE: instead of posting all the IRECVs and all the ISEND (see commits
  ![       prior to August 10 2020), we instead use the 'remap' to post
  ![       directionwise IRECV-SEND pairs. 
  do lc = 1, d2aSubarraysR
    ![ Post a receive for each process that this rank expects to receive data from.
    if (d2aRanksR%rankID(lc) == cartRank) then
      ![ Copy instead of self-messaging.
      idx = d2aRanksS%selfIdx
      fkaOut(d2aRanksR%idxLower(lc,1):d2aRanksR%idxUpper(lc,1),d2aRanksR%idxLower(lc,2):d2aRanksR%idxUpper(lc,2)) &
        = fkIn(d2aRanksS%idxLower(idx,1):d2aRanksS%idxUpper(idx,1),d2aRanksS%idxLower(idx,2):d2aRanksS%idxUpper(idx,2))
      d2aRanksR%req(lc) = MPI_REQUEST_NULL
    else
      call MPI_IRECV(fkaOut,1,d2aRanksR%mpiDDtype(lc),d2aRanksR%rankID(lc), &
                     1,cartCOMM,d2aRanksR%req(lc),errFlagMPI)
    endif
  enddo
  do lc = 1, d2aSubarraysS
    if (d2aRanksS%rankID(lc) == cartRank) then
      ![ Do not send to self.
      d2aRanksS%req(lc) = MPI_REQUEST_NULL
    else
      call MPI_ISEND(fkIn,1,d2aRanksS%mpiDDtype(lc),d2aRanksS%rankID(lc), &
                     1,cartCOMM,d2aRanksS%req(lc),errFlagMPI)
    endif
  enddo
  call MPI_WAITALL(d2aSubarraysS,d2aRanksS%req,d2aRanksS%stat,errFlagMPI)
  call MPI_WAITALL(d2aSubarraysR,d2aRanksR%req,d2aRanksR%stat,errFlagMPI)
#endif
#else
  fkaOut(1:NekyG, 1:NkxG)                         = fkIn(:,1:NkxG)
  fkaOut(1:NekyG, (NekxaG-(NekxG-NkxG)+1):NekxaG) = fkIn(:,NkxG+1:NekxG)
#endif

  end subroutine placeInAliased
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine takeFromAliased(fkaIn,fkOut) 
  ![ Take de-alised data from an aliased array and place it in a de-aliased array.
  implicit none
  double complex, intent(in)    :: fkaIn(:,:)
  double complex, intent(inout) :: fkOut(:,:)
#if (useMPI > 0)
#if (FFTmode == 0)
  integer                       :: idx,recvCount,recvStat(MPI_STATUS_SIZE)
#else
  integer                       :: idx
#endif
#endif

#if (useMPI > 0)
#if (FFTmode == 0)
  if (fftMember) then
    ![ First place the data into the smaller, de-aliased arrays.
    do lc = 1, d2aSubarraysS
      ![ Post a receive for each process that this rank expects to receive data from.
      if (d2aRanksS%rankID(lc) == fftRank) then
        ![ Copy instead of self-messaging.
        idx = d2aRanksR%selfIdx
        funkyG(d2aRanksS%idxLower(lc,1):d2aRanksS%idxUpper(lc,1),d2aRanksS%idxLower(lc,2):d2aRanksS%idxUpper(lc,2)) &
          = fkaIn(d2aRanksR%idxLower(idx,1):d2aRanksR%idxUpper(idx,1),d2aRanksR%idxLower(idx,2):d2aRanksR%idxUpper(idx,2))
        d2aRanksS%req(lc) = MPI_REQUEST_NULL
      else
        call MPI_IRECV(funkyG,1,d2aRanksS%mpiDDtype(lc),d2aRanksS%rankID(lc), &
                       1,fftCOMM,d2aRanksS%req(lc),errFlagMPI)
      endif
    enddo
    do lc = 1, d2aSubarraysR
      if (d2aRanksR%rankID(lc) == fftRank) then
        ![ Do not send to self.
        d2aRanksR%req(lc) = MPI_REQUEST_NULL
      else
        call MPI_ISEND(fkaIn,1,d2aRanksR%mpiDDtype(lc),d2aRanksR%rankID(lc), &
                       1,fftCOMM,d2aRanksR%req(lc),errFlagMPI)
      endif
    enddo
  
    call MPI_WAITALL(d2aSubarraysR,d2aRanksR%req,d2aRanksR%stat,errFlagMPI)
    call MPI_WAITALL(d2aSubarraysS,d2aRanksS%req,d2aRanksS%stat,errFlagMPI)
  endif

  ![ Then scatter the data along ky.
  if (xyIDs(1) == fftKyID) then
    do lc = 1, xyProcs(2)
      if (kyGatherType%rankID(lc) == xyIDs(1)) then
        ![ Copy instead of self-messaging.
        fkOut = funkyG(kyGatherType%kyIdxLower(lc):kyGatherType%kyIdxUpper(lc),:)
        kyGatherType%req(lc) = MPI_REQUEST_NULL
      else
        call MPI_Isend(funkyG,1,kyGatherType%mpiDDtype(lc),kyGatherType%rankID(lc), &
                       1,xyCOMMs(1),kyGatherType%req(lc),errFlagMPI)
      endif
    enddo
    call MPI_Waitall(xyProcs(2),kyGatherType%req,kyGatherType%stat,errFlagMPI)
  else
    recvCount = Nekx(1)*Nekx(2)
    call MPI_Recv(fkOut,recvCount,MPI_DOUBLE_COMPLEX,fftKyID,1,xyCOMMs(1),recvStat,errFlagMPI)
  endif

#else
  ![ Placing the data in the smaller, d2-aliased arrays requires re-organizing
  ![ the data across the 2D communicator.
  ![ NOTE: instead of posting all the IRECVs and all the ISEND (see commits
  ![       prior to August 10 2020), we instead use the 'remap' to post
  ![       directionwise IRECV-SEND pairs. 
  do lc = 1, d2aSubarraysS
    ![ Post a receive for each process that this rank expects to receive data from.
    if (d2aRanksS%rankID(lc) == cartRank) then
      ![ Copy instead of self-messaging.
      idx = d2aRanksR%selfIdx
      fkOut(d2aRanksS%idxLower(lc,1):d2aRanksS%idxUpper(lc,1),d2aRanksS%idxLower(lc,2):d2aRanksS%idxUpper(lc,2)) &
        = fkaIn(d2aRanksR%idxLower(idx,1):d2aRanksR%idxUpper(idx,1),d2aRanksR%idxLower(idx,2):d2aRanksR%idxUpper(idx,2))
      d2aRanksS%req(lc) = MPI_REQUEST_NULL
    else
      call MPI_IRECV(fkOut,1,d2aRanksS%mpiDDtype(lc),d2aRanksS%rankID(lc), &
                     1,cartCOMM,d2aRanksS%req(lc),errFlagMPI)
    endif
  enddo
  do lc = 1, d2aSubarraysR
    if (d2aRanksR%rankID(lc) == cartRank) then
      ![ Do not send to self.
      d2aRanksR%req(lc) = MPI_REQUEST_NULL
    else
      call MPI_ISEND(fkaIn,1,d2aRanksR%mpiDDtype(lc),d2aRanksR%rankID(lc), &
                     1,cartCOMM,d2aRanksR%req(lc),errFlagMPI)
    endif
  enddo
  call MPI_WAITALL(d2aSubarraysR,d2aRanksR%req,d2aRanksR%stat,errFlagMPI)
  call MPI_WAITALL(d2aSubarraysS,d2aRanksS%req,d2aRanksS%stat,errFlagMPI)
#endif
#else
  fkOut(:,1:NkxG)       = fkaIn(1:NekyG,1:NkxG)
  fkOut(:,NkxG+1:NekxG) = fkaIn(1:NekyG,(NekxaG-(NekxG-NkxG)+1):NekxaG)
#endif

  end subroutine takeFromAliased
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine poisson_bracket(fk, gk, PBfgk)
  ![ Compute the Poisson Bracket [f, g] = d(f)/dx*d(g)/dy - d(f)/dy*d(g)/dx.
  ![ Do this by computing the derivatives spectrally, then transforming each of them
  ![ to real space, calculating [f,g] in real space and transforming the result
  ![ back to Fourier space.
  ![ It is unclear whether one should de-alias the derivatives or de-alias the end
  ![ product.
  implicit none
  double complex, intent(in)  :: fk(:,:), gk(:,:)
  double complex, intent(out) :: PBfgk(:,:)

  fka = dcmplx(0.0d0,0.0d0)
  gka = dcmplx(0.0d0,0.0d0)

  ![ Place de-aliased arrays into larger, aliased arrays.
  call placeInAliased(fk,fka)
  call placeInAliased(gk,gka)

  do lckx = 1,Nekxa(2)
    do lcky = 1,Nekxa(1)
      fka_x(lcky, lckx) = ikxa(lckx)*fka(lcky,lckx)
      fka_y(lcky, lckx) = ikya(lcky)*fka(lcky,lckx)
      gka_x(lcky, lckx) = ikxa(lckx)*gka(lcky,lckx)
      gka_y(lcky, lckx) = ikya(lcky)*gka(lcky,lckx)
    enddo
  enddo

  ![ Perform complex to real transforms of the derivatives in Poisson bracket. 
  call FFT2Da_c2r(fka_x, fRa_x) 
  call FFT2Da_c2r(fka_y, fRa_y) 
  call FFT2Da_c2r(gka_x, gRa_x) 
  call FFT2Da_c2r(gka_y, gRa_y)

#if ((useMPI > 0) && (FFTmode == 0))
  forall(lcy=1:NyaG,lcx=1:Nxa(2)) PBfgRa(lcy,lcx) = fRa_x(lcy,lcx)*gRa_y(lcy,lcx) &
                                                   -fRa_y(lcy,lcx)*gRa_x(lcy,lcx) 
#else
  forall(lcx=1:Nxa(2),lcy=1:Nxa(1)) PBfgRa(lcx,lcy) = fRa_x(lcx,lcy)*gRa_y(lcx,lcy) &
                                                     -fRa_y(lcx,lcy)*gRa_x(lcx,lcy) 
#endif

  call FFT2Da_r2c(PBfgRa, PBfgka)

  ![ Take data from aliased-Fourier array into de-aliased Fourier array:
  call takeFromAliased(PBfgka,PBfgk)

  end subroutine poisson_bracket

end module nonlinearities

![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

module timeSteppers
  ![ Module with different time integrators to be tried.
  use parameters
  use fields
  use nonlinearities
  use timeUpdates, only: getPhik

  contains

#if (adiabaticElc == 0)
  subroutine denElcTimeRateOfChange(denIn,tempIn,denDotOut)
  ![ Compute d(n_e)/dt.
  implicit none
  double complex, intent(in)  :: denIn(:,:),tempIn(:,:)
  double complex, intent(out) :: denDotOut(:,:)

#if (NONLINEAR != 0)
  ![ Compute the Poisson bracket in real-space.
  call poisson_bracket(denIn, elcGyroPhik, PBk) 
#endif

  forall(lcky=1:Nekx(1),lckx=1:Nekx(2)) &
    denDotOut(lcky,lckx) = (PBk(lcky,lckx)+(iOmStElc(lcky)-2.0d0*iOmdElc(lcky))*elcGyroPhik(lcky,lckx) &
      -Lnorm*iOmdElc(lcky)*(deltae*denIn(lcky,lckx)+tempIn(lcky,lckx)))*denElcDotDiffFac(lcky,lckx)

  end subroutine denElcTimeRateOfChange
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine tempElcTimeRateOfChange(denIn,tempIn,tempDotOut)
  ![ Compute d(T_{perp e})/dt.
  implicit none
  double complex, intent(in)  :: denIn(:,:),tempIn(:,:)
  double complex, intent(out) :: tempDotOut(:,:)

#if (NONLINEAR != 0)
  ![ Compute the Poisson bracket in real-space.
  call poisson_bracket(tempIn, elcGyroPhik, PBk) 
  call poisson_bracket(denIn, hatLapD2elcGyroPhik, flrPBk) 
#endif

  forall(lcky=1:Nekx(1),lckx=1:Nekx(2)) &
    tempDotOut(lcky,lckx) = (PBk(lcky,lckx)+deltaPerpe*(flrPBk(lcky,lckx) &
      +((1.0d0+eta_e)*iOmStElc(lcky)-3.0d0*iOmdElc(lcky))*elcGyroPhik(lcky,lckx)))*tempElcDotDiffFac(lcky,lckx)

  end subroutine tempElcTimeRateOfChange
#endif
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
#if (adiabaticIon == 0)
  subroutine denIonTimeRateOfChange(denIn,tempIn,denDotOut)
  ![ Compute d(n_i)/dt.
  implicit none
  double complex, intent(in)  :: denIn(:,:),tempIn(:,:)
  double complex, intent(out) :: denDotOut(:,:)

#if (NONLINEAR != 0)
  ![ Compute the Poisson bracket in real-space.
  call poisson_bracket(denIn, ionGyroPhik, PBk) 
#endif

  forall(lcky=1:Nekx(1),lckx=1:Nekx(2)) &
    denDotOut(lcky,lckx) = (PBk(lcky,lckx)+(1.0d0/tau)*(iOmStIon(lcky)-2.0d0*iOmdIon(lcky))*ionGyroPhik(lcky,lckx) &
      -Lnorm*(1.0d0/tau)*iOmdIon(lcky)*(tau*deltai*denIn(lcky,lckx)+tempIn(lcky,lckx)))*denIonDotDiffFac(lcky,lckx)

  end subroutine denIonTimeRateOfChange
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine tempIonTimeRateOfChange(denIn,tempIn,tempDotOut)
  ![ Compute d(T_{perp i})/dt.
  implicit none
  double complex, intent(in)  :: denIn(:,:),tempIn(:,:)
  double complex, intent(out) :: tempDotOut(:,:)

#if (NONLINEAR != 0)
  ![ Compute the Poisson bracket in real-space.
  call poisson_bracket(tempIn, ionGyroPhik, PBk) 
  call poisson_bracket(denIn, hatLapD2ionGyroPhik, flrPBk) 
#endif

  forall(lcky=1:Nekx(1),lckx=1:Nekx(2)) &
    tempDotOut(lcky,lckx) = (PBk(lcky,lckx)+deltaPerpi*(tau*flrPBk(lcky,lckx) &
      +((1.0d0+eta_i)*iOmStIon(lcky)-3.0d0*iOmdIon(lcky))*ionGyroPhik(lcky,lckx)))*tempIonDotDiffFac(lcky,lckx)

  end subroutine tempIonTimeRateOfChange
#endif
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine eulerStep(fldIn, fldDot, dtIn, fldOut)
  ![ Take an Euler step of size dtIn.
  implicit none
  double complex, intent(in)   :: fldIn(:,:), fldDot(:,:)
  double precision, intent(in) :: dtIn
  double complex, intent(out)  :: fldOut(:,:)

  fldOut = fldIn+dtIn*fldDot

  end subroutine eulerStep
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
#if (STEPPER == 1)
  subroutine forwardEuler(dtIn)
  ![ Simple forward Euler stepping.
  implicit none
  double precision, intent(in)  :: dtIn

#if (adiabaticElc==0)
  call denElcTimeRateOfChange( denElck,tempElck, denElcDot)
  call tempElcTimeRateOfChange(denElck,tempElck,tempElcDot)
#endif
#if (adiabaticIon==0)
  call denIonTimeRateOfChange( denIonk,tempIonk, denIonDot)
  call tempIonTimeRateOfChange(denIonk,tempIonk,tempIonDot)
#endif

#if (HDtimeOP == 1)
#if (adiabaticElc==0)
  call eulerStep( denElcDiffFac* denElck, denElcDot,dtIn, denElck)
  call eulerStep(tempElcDiffFac*tempElck,tempElcDot,dtIn,tempElck)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonDiffFac* denIonk, denIonDot,dtIn, denIonk)
  call eulerStep(tempIonDiffFac*tempIonk,tempIonDot,dtIn,tempIonk)
#endif

  call getPhik(denElck,tempElck,denIonk,tempIonk)

#elif (HDtimeOP == 2)
#if (adiabaticElc==0)
  call eulerStep( denElck, denElcDot,dtIn, denElck)
  call eulerStep(tempElck,tempElcDot,dtIn,tempElck)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonk, denIonDot,dtIn, denIonk)
  call eulerStep(tempIonk,tempIonDot,dtIn,tempIonk)
#endif
#endif

  end subroutine forwardEuler
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
#elif (STEPPER == 2)
  subroutine trapezoidalLeapFrog(dtIn)
  ![ Second-order accurate trapezoidal leapfrog.
  implicit none
  double precision, intent(in)  :: dtIn

  ![ Predictor step:
#if (adiabaticElc==0)
  call denElcTimeRateOfChange( denElck,tempElck, denElcDot)
  call tempElcTimeRateOfChange(denElck,tempElck,tempElcDot)
#endif
#if (adiabaticIon==0)
  call denIonTimeRateOfChange( denIonk,tempIonk, denIonDot)
  call tempIonTimeRateOfChange(denIonk,tempIonk,tempIonDot)
#endif

#if (HDtimeOP == 1)
#if (adiabaticElc==0)
  call eulerStep(0.50d0* denElcDiffFac*( denElck+ denElcks), denElcDot,dtIn, denElckTmp)
  call eulerStep(0.50d0*tempElcDiffFac*(tempElck+tempElcks),tempElcDot,dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep(0.50d0* denIonDiffFac*( denIonk+ denIonks), denIonDot,dtIn, denIonkTmp)
  call eulerStep(0.50d0*tempIonDiffFac*(tempIonk+tempIonks),tempIonDot,dtIn,tempIonkTmp)
#endif
#elif (HDtimeOP == 2)
#if (adiabaticElc==0)
  call eulerStep(0.50d0*( denElck+ denElcks), denElcDot,dtIn, denElckTmp)
  call eulerStep(0.50d0*(tempElck+tempElcks),tempElcDot,dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep(0.50d0*( denIonk+ denIonks), denIonDot,dtIn, denIonkTmp)
  call eulerStep(0.50d0*(tempIonk+tempIonks),tempIonDot,dtIn,tempIonkTmp)
#endif
#endif

  call getPhik(denElckTmp,tempElckTmp,denIonkTmp,tempIonkTmp)

  denElcks  =  denElck
  tempElcks = tempElck
  denIonks  =  denIonk
  tempIonks = tempIonk

  ![ Corrector step:
#if (adiabaticElc==0)
  call denElcTimeRateOfChange( denElckTmp,tempElckTmp, denElcDot)
  call tempElcTimeRateOfChange(denElckTmp,tempElckTmp,tempElcDot)
#endif
#if (adiabaticIon==0)
  call denIonTimeRateOfChange( denIonkTmp,tempIonkTmp, denIonDot)
  call tempIonTimeRateOfChange(denIonkTmp,tempIonkTmp,tempIonDot)
#endif

#if (HDtimeOP == 1)
#if (adiabaticElc==0)
  call eulerStep( denElcDiffFac* denElcks, denElcDot,dtIn, denElck)
  call eulerStep(tempElcDiffFac*tempElcks,tempElcDot,dtIn,tempElck)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonDiffFac* denIonks, denIonDot,dtIn, denIonk)
  call eulerStep(tempIonDiffFac*tempIonks,tempIonDot,dtIn,tempIonk)
#endif

  call getPhik(denElck,tempElck,denIonk,tempIonk)
#elif (HDtimeOP == 2)
#if (adiabaticElc==0)
  call eulerStep( denElcks, denElcDot,dtIn, denElck)
  call eulerStep(tempElcks,tempElcDot,dtIn,tempElck)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonks, denIonDot,dtIn, denIonk)
  call eulerStep(tempIonks,tempIonDot,dtIn,tempIonk)
#endif
#endif

  end subroutine trapezoidalLeapFrog
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
#elif (STEPPER == 3)
  subroutine rk3(dtIn)
  ![ Third-order accurate Runge-Kutta.
  implicit none
  double precision, intent(in)  :: dtIn

  ![ First stage.
#if (adiabaticElc==0)
  call denElcTimeRateOfChange( denElck,tempElck, denElcDot)
  call tempElcTimeRateOfChange(denElck,tempElck,tempElcDot)
#endif
#if (adiabaticIon==0)
  call denIonTimeRateOfChange( denIonk,tempIonk, denIonDot)
  call tempIonTimeRateOfChange(denIonk,tempIonk,tempIonDot)
#endif

#if (HDtimeOP == 1)
#if (adiabaticElc==0)
  call eulerStep( denElcDiffFac* denElck, denElcDot,0.50d0*dtIn, denElckTmp)
  call eulerStep(tempElcDiffFac*tempElck,tempElcDot,0.50d0*dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonDiffFac* denIonk, denIonDot,0.50d0*dtIn, denIonkTmp)
  call eulerStep(tempIonDiffFac*tempIonk,tempIonDot,0.50d0*dtIn,tempIonkTmp)
#endif
#elif (HDtimeOP == 2)
#if (adiabaticElc==0)
  call eulerStep( denElck, denElcDot,0.50d0*dtIn, denElckTmp)
  call eulerStep(tempElck,tempElcDot,0.50d0*dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonk, denIonDot,0.50d0*dtIn, denIonkTmp)
  call eulerStep(tempIonk,tempIonDot,0.50d0*dtIn,tempIonkTmp)
#endif
#endif

  call getPhik(denElckTmp,tempElckTmp,denIonkTmp,tempIonkTmp)

  ![ Second stage.
#if (adiabaticElc==0)
  call denElcTimeRateOfChange( denElckTmp,tempElckTmp, denElcDot1)
  call tempElcTimeRateOfChange(denElckTmp,tempElckTmp,tempElcDot1)
#endif
#if (adiabaticIon==0)
  call denIonTimeRateOfChange( denIonkTmp,tempIonkTmp, denIonDot1)
  call tempIonTimeRateOfChange(denIonkTmp,tempIonkTmp,tempIonDot1)
#endif

#if (HDtimeOP == 1)
#if (adiabaticElc==0)
  call eulerStep( denElcDiffFac* denElck,2.0d0* denElcDot1- denElcDot,dtIn, denElckTmp)
  call eulerStep(tempElcDiffFac*tempElck,2.0d0*tempElcDot1-tempElcDot,dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonDiffFac* denIonk,2.0d0* denIonDot1- denIonDot,dtIn, denIonkTmp)
  call eulerStep(tempIonDiffFac*tempIonk,2.0d0*tempIonDot1-tempIonDot,dtIn,tempIonkTmp)
#endif
#elif (HDtimeOP == 2)
#if (adiabaticElc==0)
  call eulerStep( denElck,2.0d0* denElcDot1- denElcDot,dtIn, denElckTmp)
  call eulerStep(tempElck,2.0d0*tempElcDot1-tempElcDot,dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonk,2.0d0* denIonDot1- denIonDot,dtIn, denIonkTmp)
  call eulerStep(tempIonk,2.0d0*tempIonDot1-tempIonDot,dtIn,tempIonkTmp)
#endif
#endif

  call getPhik(denElckTmp,tempElckTmp,denIonkTmp,tempIonkTmp)

  ![ Third stage.
#if (adiabaticElc==0)
  call denElcTimeRateOfChange( denElckTmp,tempElckTmp, denElcDot2)
  call tempElcTimeRateOfChange(denElckTmp,tempElckTmp,tempElcDot2)
#endif
#if (adiabaticIon==0)
  call denIonTimeRateOfChange( denIonkTmp,tempIonkTmp, denIonDot2)
  call tempIonTimeRateOfChange(denIonkTmp,tempIonkTmp,tempIonDot2)
#endif

#if (HDtimeOP == 1)
#if (adiabaticElc==0)
  call eulerStep( denElcDiffFac* denElck, denElcDot+4.0d0* denElcDot1+ denElcDot2,dtIn/6.0d0, denElck)
  call eulerStep(tempElcDiffFac*tempElck,tempElcDot+4.0d0*tempElcDot1+tempElcDot2,dtIn/6.0d0,tempElck)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonDiffFac* denIonk, denIonDot+4.0d0* denIonDot1+ denIonDot2,dtIn/6.0d0, denIonk)
  call eulerStep(tempIonDiffFac*tempIonk,tempIonDot+4.0d0*tempIonDot1+tempIonDot2,dtIn/6.0d0,tempIonk)
#endif

  call getPhik(denElck,tempElck,denIonk,tempIonk)
#elif (HDtimeOP == 2)
#if (adiabaticElc==0)
  call eulerStep( denElck, denElcDot+4.0d0* denElcDot1+ denElcDot2,dtIn/6.0d0, denElck)
  call eulerStep(tempElck,tempElcDot+4.0d0*tempElcDot1+tempElcDot2,dtIn/6.0d0,tempElck)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonk, denIonDot+4.0d0* denIonDot1+ denIonDot2,dtIn/6.0d0, denIonk)
  call eulerStep(tempIonk,tempIonDot+4.0d0*tempIonDot1+tempIonDot2,dtIn/6.0d0,tempIonk)
#endif
#endif
  end subroutine rk3
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
#elif (STEPPER == 4)
  subroutine rk4(dtIn)
  ![ Fourth-order accurate Runge-Kutta.
  implicit none
  double precision, intent(in)  :: dtIn

  ![ First stage.
#if (adiabaticElc==0)
  call denElcTimeRateOfChange( denElck,tempElck, denElcDot)
  call tempElcTimeRateOfChange(denElck,tempElck,tempElcDot)
#endif
#if (adiabaticIon==0)
  call denIonTimeRateOfChange( denIonk,tempIonk, denIonDot)
  call tempIonTimeRateOfChange(denIonk,tempIonk,tempIonDot)
#endif

#if (HDtimeOP == 1)
#if (adiabaticElc==0)
  call eulerStep( denElcDiffFac* denElck, denElcDot,0.50d0*dtIn, denElckTmp)
  call eulerStep(tempElcDiffFac*tempElck,tempElcDot,0.50d0*dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonDiffFac* denIonk, denIonDot,0.50d0*dtIn, denIonkTmp)
  call eulerStep(tempIonDiffFac*tempIonk,tempIonDot,0.50d0*dtIn,tempIonkTmp)
#endif
#elif (HDtimeOP == 2)
#if (adiabaticElc==0)
  call eulerStep( denElck, denElcDot,0.50d0*dtIn, denElckTmp)
  call eulerStep(tempElck,tempElcDot,0.50d0*dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonk, denIonDot,0.50d0*dtIn, denIonkTmp)
  call eulerStep(tempIonk,tempIonDot,0.50d0*dtIn,tempIonkTmp)
#endif
#endif

  call getPhik(denElckTmp,tempElckTmp,denIonkTmp,tempIonkTmp)

  ![ Second stage.
#if (adiabaticElc==0)
  call denElcTimeRateOfChange( denElckTmp,tempElckTmp, denElcDot1)
  call tempElcTimeRateOfChange(denElckTmp,tempElckTmp,tempElcDot1)
#endif
#if (adiabaticIon==0)
  call denIonTimeRateOfChange( denIonkTmp,tempIonkTmp, denIonDot1)
  call tempIonTimeRateOfChange(denIonkTmp,tempIonkTmp,tempIonDot1)
#endif

#if (HDtimeOP == 1)
#if (adiabaticElc==0)
  call eulerStep( denElcDiffFac* denElck, denElcDot1,0.50d0*dtIn, denElckTmp)
  call eulerStep(tempElcDiffFac*tempElck,tempElcDot1,0.50d0*dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonDiffFac* denIonk, denIonDot1,0.50d0*dtIn, denIonkTmp)
  call eulerStep(tempIonDiffFac*tempIonk,tempIonDot1,0.50d0*dtIn,tempIonkTmp)
#endif
#elif (HDtimeOP == 2)
#if (adiabaticElc==0)
  call eulerStep( denElck, denElcDot1,0.50d0*dtIn, denElckTmp)
  call eulerStep(tempElck,tempElcDot1,0.50d0*dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonk, denIonDot1,0.50d0*dtIn, denIonkTmp)
  call eulerStep(tempIonk,tempIonDot1,0.50d0*dtIn,tempIonkTmp)
#endif
#endif

  call getPhik(denElckTmp,tempElckTmp,denIonkTmp,tempIonkTmp)

  ![ Third stage.
#if (adiabaticElc==0)
  call denElcTimeRateOfChange( denElckTmp,tempElckTmp, denElcDot2)
  call tempElcTimeRateOfChange(denElckTmp,tempElckTmp,tempElcDot2)
#endif
#if (adiabaticIon==0)
  call denIonTimeRateOfChange( denIonkTmp,tempIonkTmp, denIonDot2)
  call tempIonTimeRateOfChange(denIonkTmp,tempIonkTmp,tempIonDot2)
#endif

#if (HDtimeOP == 1)
#if (adiabaticElc==0)
  call eulerStep( denElcDiffFac* denElck, denElcDot2,dtIn, denElckTmp)
  call eulerStep(tempElcDiffFac*tempElck,tempElcDot2,dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonDiffFac* denIonk, denIonDot2,dtIn, denIonkTmp)
  call eulerStep(tempIonDiffFac*tempIonk,tempIonDot2,dtIn,tempIonkTmp)
#endif
#elif (HDtimeOP == 2)
#if (adiabaticElc==0)
  call eulerStep( denElck, denElcDot2,dtIn, denElckTmp)
  call eulerStep(tempElck,tempElcDot2,dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonk, denIonDot2,dtIn, denIonkTmp)
  call eulerStep(tempIonk,tempIonDot2,dtIn,tempIonkTmp)
#endif
#endif

  call getPhik(denElckTmp,tempElckTmp,denIonkTmp,tempIonkTmp)

  ![ Fourth stage.
#if (adiabaticElc==0)
  call denElcTimeRateOfChange( denElckTmp,tempElckTmp, denElcDot3)
  call tempElcTimeRateOfChange(denElckTmp,tempElckTmp,tempElcDot3)
#endif
#if (adiabaticIon==0)
  call denIonTimeRateOfChange( denIonkTmp,tempIonkTmp, denIonDot3)
  call tempIonTimeRateOfChange(denIonkTmp,tempIonkTmp,tempIonDot3)
#endif

#if (HDtimeOP == 1)
#if (adiabaticElc==0)
  call eulerStep( denElcDiffFac* denElck, denElcDot+2.0d0*( denElcDot1+ denElcDot2)+ denElcDot3,dtIn/6.0d0, denElck)
  call eulerStep(tempElcDiffFac*tempElck,tempElcDot+2.0d0*(tempElcDot1+tempElcDot2)+tempElcDot3,dtIn/6.0d0,tempElck)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonDiffFac* denIonk, denIonDot+2.0d0*( denIonDot1+ denIonDot2)+ denIonDot3,dtIn/6.0d0, denIonk)
  call eulerStep(tempIonDiffFac*tempIonk,tempIonDot+2.0d0*(tempIonDot1+tempIonDot2)+tempIonDot3,dtIn/6.0d0,tempIonk)
#endif

  call getPhik(denElck,tempElck,denIonk,tempIonk)
#elif (HDtimeOP == 2)
#if (adiabaticElc==0)
  call eulerStep( denElck, denElcDot+2.0d0*( denElcDot1+ denElcDot2)+ denElcDot3,dtIn/6.0d0, denElck)
  call eulerStep(tempElck,tempElcDot+2.0d0*(tempElcDot1+tempElcDot2)+tempElcDot3,dtIn/6.0d0,tempElck)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonk, denIonDot+2.0d0*( denIonDot1+ denIonDot2)+ denIonDot3,dtIn/6.0d0, denIonk)
  call eulerStep(tempIonk,tempIonDot+2.0d0*(tempIonDot1+tempIonDot2)+tempIonDot3,dtIn/6.0d0,tempIonk)
#endif
#endif

  end subroutine rk4
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
#elif (STEPPER == 5)
  subroutine sspRK3(dtIn)
  ![ Strong stability preserving RK3.
  implicit none
  double precision, intent(in)  :: dtIn

  ![ First stage.
#if (adiabaticElc==0)
  call denElcTimeRateOfChange( denElck,tempElck, denElcDot)
  call tempElcTimeRateOfChange(denElck,tempElck,tempElcDot)
#endif
#if (adiabaticIon==0)
  call denIonTimeRateOfChange( denIonk,tempIonk, denIonDot)
  call tempIonTimeRateOfChange(denIonk,tempIonk,tempIonDot)
#endif

#if (HDtimeOP == 1)
#if (adiabaticElc==0)
  call eulerStep( denElcDiffFac* denElck, denElcDot,dtIn, denElckTmp)
  call eulerStep(tempElcDiffFac*tempElck,tempElcDot,dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonDiffFac* denIonk, denIonDot,dtIn, denIonkTmp)
  call eulerStep(tempIonDiffFac*tempIonk,tempIonDot,dtIn,tempIonkTmp)
#endif
#elif (HDtimeOP == 2)
#if (adiabaticElc==0)
  call eulerStep( denElck, denElcDot,dtIn, denElckTmp)
  call eulerStep(tempElck,tempElcDot,dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonk, denIonDot,dtIn, denIonkTmp)
  call eulerStep(tempIonk,tempIonDot,dtIn,tempIonkTmp)
#endif
#endif

  call getPhik(denElckTmp,tempElckTmp,denIonkTmp,tempIonkTmp)

  ![ Second stage.
#if (adiabaticElc==0)
  call denElcTimeRateOfChange( denElckTmp,tempElckTmp, denElcDot1)
  call tempElcTimeRateOfChange(denElckTmp,tempElckTmp,tempElcDot1)
#endif
#if (adiabaticIon==0)
  call denIonTimeRateOfChange( denIonkTmp,tempIonkTmp, denIonDot1)
  call tempIonTimeRateOfChange(denIonkTmp,tempIonkTmp,tempIonDot1)
#endif

#if (HDtimeOP == 1)
#if (adiabaticElc==0)
  call eulerStep( denElcDiffFac* denElck, denElcDot+ denElcDot1,0.250d0*dtIn, denElckTmp)
  call eulerStep(tempElcDiffFac*tempElck,tempElcDot+tempElcDot1,0.250d0*dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonDiffFac* denIonk, denIonDot+ denIonDot1,0.250d0*dtIn, denIonkTmp)
  call eulerStep(tempIonDiffFac*tempIonk,tempIonDot+tempIonDot1,0.250d0*dtIn,tempIonkTmp)
#endif
#elif (HDtimeOP == 2)
#if (adiabaticElc==0)
  call eulerStep( denElck, denElcDot+ denElcDot1,0.250d0*dtIn, denElckTmp)
  call eulerStep(tempElck,tempElcDot+tempElcDot1,0.250d0*dtIn,tempElckTmp)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonk, denIonDot+ denIonDot1,0.250d0*dtIn, denIonkTmp)
  call eulerStep(tempIonk,tempIonDot+tempIonDot1,0.250d0*dtIn,tempIonkTmp)
#endif
#endif

  call getPhik(denElckTmp,tempElckTmp,denIonkTmp,tempIonkTmp)

  ![ Third stage.
#if (adiabaticElc==0)
  call denElcTimeRateOfChange( denElckTmp,tempElckTmp, denElcDot2)
  call tempElcTimeRateOfChange(denElckTmp,tempElckTmp,tempElcDot2)
#endif
#if (adiabaticIon==0)
  call denIonTimeRateOfChange( denIonkTmp,tempIonkTmp, denIonDot2)
  call tempIonTimeRateOfChange(denIonkTmp,tempIonkTmp,tempIonDot2)
#endif

#if (HDtimeOP == 1)
#if (adiabaticElc==0)
  call eulerStep( denElcDiffFac* denElck, denElcDot+ denElcDot1+4.0d0* denElcDot2,dtIn/6.0d0, denElck)
  call eulerStep(tempElcDiffFac*tempElck,tempElcDot+tempElcDot1+4.0d0*tempElcDot2,dtIn/6.0d0,tempElck)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonDiffFac* denIonk, denIonDot+ denIonDot1+4.0d0* denIonDot2,dtIn/6.0d0, denIonk)
  call eulerStep(tempIonDiffFac*tempIonk,tempIonDot+tempIonDot1+4.0d0*tempIonDot2,dtIn/6.0d0,tempIonk)
#endif

  call getPhik(denElck,tempElck,denIonk,tempIonk)
#elif (HDtimeOP == 2)
#if (adiabaticElc==0)
  call eulerStep( denElck, denElcDot+ denElcDot1+4.0d0* denElcDot2,dtIn/6.0d0, denElck)
  call eulerStep(tempElck,tempElcDot+tempElcDot1+4.0d0*tempElcDot2,dtIn/6.0d0,tempElck)
#endif
#if (adiabaticIon==0)
  call eulerStep( denIonk, denIonDot+ denIonDot1+4.0d0* denIonDot2,dtIn/6.0d0, denIonk)
  call eulerStep(tempIonk,tempIonDot+tempIonDot1+4.0d0*tempIonDot2,dtIn/6.0d0,tempIonk)
#endif
#endif

  end subroutine sspRK3
#endif
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  subroutine stepSolution
  ![ Advance the solution one time step.
  implicit none

#if (STEPPER == 1)
  call forwardEuler(dt)

#elif (STEPPER == 2)
  call trapezoidalLeapFrog(dt)

#elif (STEPPER == 3)
  call rk3(dt)

#elif (STEPPER == 4)
  call rk4(dt)

#elif (STEPPER == 5)
  call sspRK3(dt)

#endif


#if (HDtimeOP == 2)
  ![ Apply hyperdiffusion.
#if (adiabaticElc==0)
  denElck  =  denElcDiffFac* denElck
  tempElck = tempElcDiffFac*tempElck
#endif
#if (adiabaticIon==0)
  denIonk  =  denIonDiffFac* denIonk
  tempIonk = tempIonDiffFac*tempIonk
#endif

  ![ Compute the electrostatic potential.
  call getPhik(denElck,tempElck,denIonk,tempIonk)
#endif

  ![ Compute the vorticity in Fourier space. 
  vortk = -kSq*phik

  end subroutine stepSolution

end module timeSteppers

![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
![ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

program shroom
  use initialize
  use terminate
  use timeSteppers
  use IOtools
  use timeUpdates
#if (useMPI > 0)
  use MPItools
#endif
  implicit none
  character(8)  :: rdate
  character(10) :: rtime

#if (useMPI == 0)

#if (ioMode==1)
  ![ Check that ADIOS is requested along with MPI.
  if (myrank==0) then
    write(*,'( )')
    write(*,'(" --> Need MPI for ADIOS IO." )')
    write(*,'(" TERMINATING ")')
  endif
  error stop
#endif

  myRank = 0

#else

  call init_MPI    ![ Initialize MPI interface.

#endif

  if (myrank==0) then
    write(*,'( )')
    write(*,'("     --> Welcome to MuSHrooM <--    " )')
    write(*,'( )')
  endif

  call init_all    ![ Initialize arrays, functions and ICs.

  if (.not. isRestart) call writeFields       ![ Write initial conditions.

  call date_and_time(rdate,rtime)
  if (myRank == 0) then
    write(*,'( )')
    write(*,'(" Entering time loop on ",A8," at ",A10)') rdate, rtime
    write(*,'("            FASTEN YOUR SEATBELT!         ")')
    write(*,'( )')
  endif
  ![ ............... START OF TIME LOOP ............... ]!
  do while (timeStepMore)

    call stepSolution

    timeSteps = timeSteps + 1    ![ Time steps taken.
    simTime   = simTime + dt     ![ Current simulation time.

    if (abs(simTime - dble(framesOut+1)*tRateOutput) <= tTol) then
      call isNaNorInf    ![ Check if solutions diverged.

      framesOut = framesOut+1

      call writeFields    ![ Append fields out output files.

      call outRestart     ![ Write file used for restarts.
    endif

    if (abs(simTime - dble(logEntries+1)*tRateLogEntry) <= tTol) then
      call date_and_time(rdate,rtime)
      if (myRank == 0) write(*,'(" Completed ",i9," steps on ",A8," at ",A10," | dt=",E14.4E3," | frames out=",i5)') & 
        timeSteps, rdate, rtime, dt, framesOut
      logEntries = logEntries + 1    ![ Number of messages to log/screen.
    endif

    if (abs(simTime-endTime) <= tTol) then    ![ Time loop completed. Punch out.
      timeStepMore = .false.
      exit
    endif

#if (((HDmodel > 2) && (HDmodel < 7)) || (HDmodel > 8))
    if (abs(simTime - dble(hdAdjusts+1)*tRateAdjustHD) <= tTol) then
      ![ Dynamically adjust hyperdiffusion.
      call adjustHyperDiff(phik,vortk)
      hdAdjusts = hdAdjusts + 1
    endif
#endif

#if (ADJUSTdt > 0)
    if (abs(simTime - dble(dtAdjusts+1)*tRateAdjustDt) <= tTol) then
      ![ Dynamically adjust the time step.
      dt = dtUpdate(phik,dt)
      dtAdjusts = dtAdjusts + 1
    endif
#endif

  end do
  ![ ............... END OF TIME LOOP ............... ]!

  call date_and_time(rdate,rtime)
  if (myRank == 0) then
    write(*,'( )')
    write(*,'(" Time loop finished on ",A8," at ",A10,". Did ",i9," steps | frames out=",i5)') & 
      rdate, rtime, timeSteps, framesOut
    write(*,'( )')
    write(*,'(" Freeing MuSHrooM memory ...   " )')
  endif
  call terminate_all    ![ Free memory, deallocate arrays and pointers.
 
  if (myRank == 0) then
    write(*,'( )')
    write(*,'("     <-- MuSHrooM Finalized -->    " )')
    write(*,'( )')
  endif

#if (useMPI > 0)
  call terminate_MPI    ![ Finalize MPI.
#endif

end program shroom


