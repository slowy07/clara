![clara-image](.github/clara.png)

![Clara ubuntu test](https://img.shields.io/github/actions/workflow/status/slowy07/clara/cpp-testing.yml?style=flat-square&logo=github&label=Clara%20Ubuntu%20Test) 
![GitHub Workflow Status (with event)](https://img.shields.io/github/actions/workflow/status/slowy07/clara/docker-testing.yml?style=flat-square&logo=docker&label=Docker%20build)


Quantum computing is a type of computing that harnesses the power of quantum mechanics to solve problems that are too difficult for classical computers. Quantum computers use qubits, which are quantum bits, to store information. Qubits can be in a superposition of states, which means that they can be both 0 and 1 at the same time. This allows quantum computers to perform calculations that are exponentially faster than classical computers.

Clara are open sourec library for quantum computing processing that is designed to be faster on classical computer. the library provides a high-level interface for writing quantum algorithm and it uses a variety of techniques to optimize the performance of these algorithm on classical computer. the lib has the potential to make quantum algorithm more accessible to a wider range of researcher and developers and could help to accelerate the develompent of quantum computing applications

## Package used on clara

### Eigen

Eigen is a high-performance C++ template library for linear algebra. It provides a wide range of features, including matrices, vectors, numerical solvers, and related algorithms. Eigen is designed to be fast, efficient, and flexible. It is used by a wide range of developers, including those working in the fields of computer graphics, machine learning, and scientific computing.

website : [eigen.tuxfamily.org](https://eigen.tuxfamily.org/index.php?title=Main_Page)

repository: [gitlab](https://gitlab.com/libeigen/eigen)

*installation for clara*

- debian package (apt)
    for installing lib eigen on ubuntu
    ```
    sudo apt-get update
    sudo apt-get install libeigen3-dev
    cp -r /usr/local/include/eigen3/Eigen /usr/local/include
    ```
- arch package
    for installing on arch based
    ```
    yay -S eigen
    ```
    information about package on arch : [here](https://archlinux.org/packages/extra/any/eigen/)
- MacOS
    for installing on MacOS based
    ```
    brew install eigen
    ```
- windows
    recommended for install [Cygwin](https://www.cygwin.com/) with [CMake](http://www.cmake.org/) and [g++](https://gcc.gnu.org/) 

## OpenMp

OpenMP is a widely adopted, portable, and scalable application programming interface (API) for shared-memory parallel programming in C, C++, and Fortran. It provides a high-level, portable abstraction for expressing parallel computation on shared-memory systems.

website [openmp.org](https://www.openmp.org/)

*installation for clara*

- debian package (apt)
    for installing libeigen on ubuntu
    ```
    sudo apt-get install
    ```
- arch package
    for installing libeigen on arch based
    ```
    yay -S openmp
    ```
- MacOS
    ```
    brew install llvm
    brew install libomp
    ```

## Test clara

Clara using [GoogleTest](https://github.com/google/googletest) framework and mocking tools, for installation you can check on [here](https://github.com/google/googletest), or you can just find on the testing script on [``clara_test/run_test``](clara_test/run_test)

## running test clara

```
g++ -pedantic -std=c++11 -Wall -Wextra -Weffc++ -fopenmp -g3 -DDEBUG -isystem $HOME/eigen -I $HOME/clara/include testing/channels.cpp -o channels
```

## running test clara with docker
```
docker build -t clara-test -f ./docker/Dockerfile.clara_test . && docker run -it --name=clara-test clara-test
```

information:
- ``-pedantic`` ensure that the code conforms to the strictest standard of the C++ language
- ``-std=c++11`` specifies that the code should be compiled using c++11
- ``-Wall`` enable all warnings
- ``-Wextra`` enables additional warnings that are not enabled by default
- ``-Weffc++`` enables warnings about potential erroes in the code
- ``-fopenmp`` enables ``OpenMP`` support
- ``-g3`` enable debugging symbol
- ``docker`` best pratices for testing application if you are using windows

## Reference

reference (journal/article/paper) used in clara can check in [``REFERENCE.md``](REFERENCE.md)
