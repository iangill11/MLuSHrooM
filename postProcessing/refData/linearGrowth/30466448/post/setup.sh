
##[ Jay Lang's Kali computer:
#export FCOMPILER=/home/linuxbrew/.linuxbrew/Cellar/gcc/5.5.0_4/bin/gfortran-5
#export FFTW_DIR=/home/linuxbrew/.linuxbrew/Cellar/fftw/3.3.8/

##[ Darin Ernst's computer:
#export FCOMPILER=gfortran-mp-6
#export FFTW_DIR=/opt/local

##[ Mana's TheFarm:
#export FCOMPILER=/usr/local/Cellar/gcc/8.2.0/bin/gfortran-8
#export FFTW_DIR=/usr/local/Cellar/fftw/3.3.8/

#[ Plankton:
#export FCOMPILER=/opt/intel/compilers_and_libraries_2017.7.259/linux/bin/intel64/ifort
#export FFTW_DIR=/usr/local/

##[ Discovery:
#module unload intel-compilers
#module unload openmpi
##module load intel-compilers/13.0
##module load mpich3/3.0.4-intel13.0
##module load fftw/3.3.3-mpich3.0-intel13.0
##export FFTW_DIR=/opt/fftw/3.3.3-mpich3-intel13.0/
##export ADIOS_DIR=/ihome/mana/documents/multiscale/code/adios-1.13.1_intel13.0-mpich3.0.4/
#module load intel-compilers/17.0
#module load openmpi/1.10.1-intel17.0
#module load fftw/3.3.4-openmpi1.10.1-intel17.0
#export FFTW_DIR=/opt/fftw/3.3.4-openmpi1.10.1-intel17.0/
#export ADIOS_DIR=/ihome/mana/documents/multiscale/code/adios-1.13.1_intel17.0-openmpi.1.10.1/
#export FCOMPILER=mpif90
#export ADIOS_INC=`$ADIOS_DIR/bin/adios_config -fc`
#export ADIOS_FLIB=`$ADIOS_DIR/bin/adios_config -fl`
#export ADIOS_FREADLIB=`$ADIOS_DIR/bin/adios_config -flr`

#[ MIT's Engaging:
#module use /home/software/psfc/modulefiles
#module add psfc/config
#module load intel/2017-01
#module load impi/2017-01
#module load psfc/fftw/intel17/3.3.5
#module load psfc/adios/intel-2017/1.13.1
#export FFTW_DIR=$FFTW3DIR
#export ADIOS_DIR=$ADIOS_ROOT
#export MPI_DIR=
#export FCOMPILER=mpiifort
#module load gcc/8.3.0
##export FFTW_DIR=$HOME/multiscale/code/fftw3.3.8/
##export ADIOS_DIR=$HOME/multiscale/code/adios1.13.1/
##export MPI_DIR=$HOME/multiscale/code/openmpi4.0.4/
#export FFTW_DIR=$HOME/multiscale/code/fftw3.3.8-openmpi3.1.6/
#export ADIOS_DIR=$HOME/multiscale/code/adios1.13.1-openmpi3.1.6/
#export MPI_DIR=$HOME/multiscale/code/openmpi3.1.6/
#export FCOMPILER=$MPI_DIR/bin/mpif90
#export FFTW3DIR=$FFTW_DIR
#export FFTW3LIB=$FFTW_DIR/lib/
#export FFTW3INCLUDE=$FFTW_DIR/include/
module load intel/2018-01
module load impi/2018-01
export FFTW_DIR=$HOME/multiscale/code/fftw3.3.8-intel18/
export ADIOS_DIR=$HOME/multiscale/code/adios1.13.1-intel18/
export MPI_DIR=
export FCOMPILER=mpiifort

export ADIOS_INC=`$ADIOS_DIR/bin/adios_config -fc`
export ADIOS_FLIB=`$ADIOS_DIR/bin/adios_config -fl`
export ADIOS_FREADLIB=`$ADIOS_DIR/bin/adios_config -flr`
