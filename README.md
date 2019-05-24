## Least Squares Meshes
This is an libigl implementation of [Least-squares Meshes paper](https://igl.ethz.ch/projects/Laplacian-mesh-processing/ls-meshes/ls-meshes.pdf)

# Building
On Unix systems, g++ and CMake is required. The steps are as follows
cd $PROJECT_FOLDER
mkdir Build && cd build
cmake ..
make


# Usage
cd $PROJECT_FOLDER/build
./my_project_bin [F] [M]

F: desired input mesh's path
M: number of sample points to be used for mesh reconsturction
