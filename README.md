## Least Squares Meshes
This is an libigl implementation of [Least-squares Meshes paper](https://igl.ethz.ch/projects/Laplacian-mesh-processing/ls-meshes/ls-meshes.pdf)

# Building
On Unix systems, g++ and CMake is required. The steps are as follows <br/>
cd $PROJECT_FOLDER <br/>
mkdir Build && cd build <br/>
cmake .. <br/>
make 


# Usage
`cd $PROJECT_FOLDER/build` <br/>
If you compile without any change main.cpp file which includes basic implementation with random and FPS point selection methods will be built. If you make changes in `CMakeLists.txt` file and compile with `main_2.cpp` the brute force algorithm will run. If `main_3.cpp` file is build point selection based on maximum number of edges and maximum curvature can be used. <br/>
`./my_project_bin [F] [N] [M]` <br/>

F: desired input mesh's path <br/>
N: number of sample points to be used for mesh reconstruction <br/>
M: denotes the method for choosing sample points. 0 for FPS and
1 for random point sampling, If you compile main_3.cpp you can also choose 2 for maximum edge point selection and 3 for maximum curvature selection
