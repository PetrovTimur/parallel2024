# Conjugate Gradient Solver

A C++ project for solving large sparse systems using parallel algorithms (OpenMP, MPI, CUDA). Includes utilities for logging, argument parsing, and matrix operations.

## Requirements

- CMake >= 3.12
- C++14 compiler (GCC, Clang, etc.)
- OpenMP (for parallel CPU execution)
- MPI (optional, for distributed execution)
- CUDA Toolkit (optional, for GPU support)

## Building the Project

1. **Clone the repository:**
    ```sh
    git clone https://github.com/PetrovTimur/parallel2024.git
    cd parallel2024
    ```

2. **Create a build directory:**
   ```sh
   mkdir build
   cd build
   ```

3. **Configure with CMake:**
    ```sh
    cmake ..
    ```

4. **Build:**
    ```sh
    make -j$(nproc)
    ```

## Build Options

| Option       | Description                          | Default |
|--------------|--------------------------------------|---------|
| DEBUG_MODE   | Enable debug logging and checks      | OFF     |
| USE_MPI      | Enable MPI parallelization           | OFF     |
| USE_CUDA     | Enable CUDA GPU support              | OFF     |
| ENABLE_TESTS | Build and enable tests               | ON      |
| MEASURES     | Build performance measuring binaries | OFF     |

Set options via `cmake -DOPTION=ON ..`.

## Running

- Executables are placed in `build/bin/`.
- Logs are written to `build/log/` by default.
- Example run (OpenMP):
    ```sh
    ./bin/CGSolver --Nx 100 --Ny 100 --K1 1 --K2 1
    ```
- Example run (MPI):
    ```sh
    mpirun -np 4 ./bin/CGSolver --Nx 100 --Ny 100 --K1 1 --K2 1 --Px 2 --Py 2
    ```

## Directory Structure

- `Solver/` — Main solvers and parallel implementations
- `Utilities/` — Logging, argument parsing, input utilities
- `Tests/` — Unit and integration tests