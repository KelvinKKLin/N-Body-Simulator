# N-Body-Simulator
This is a parallel C++ based N-body Simulator created with OpenMP and MPI for Sfwr Eng 4F03 (Parallel Computing). The code implements the
basic N-body algorithm presented in Pacheco's Introduction to Parallel Programming. 

## Running the Code
Running the code requires OpenMP and MPI (OpenRTE 1.10.6) to be installed on your computer. To run the code, run the makefile. The
makefile produces the Linux executable x.project.

x.project takes the following parameters:
- **numberOfLightParticles**: the number of light particles to be used in the simulation.
- **numberOfMediumParticles**: the number of medium particles to be used in the simulation.
- **numberOfHeavyParticles**: the number of heavy particles to be used in the simulation.
- **numberOfSteps**: the number of major steps in the simulation. Images will be produced at every major step.
- **numberOfSubsteps**: the number of substeps between each major step.
- **imageWidth**: the width of the output image
- **imageHeight**: the height of the output image
- **imageFilenamePrefix**: the prefix of the image file
