# OdeVec

OdeVec is a special-purpose ODE solver designed for solving stiff equation systems with small- to medium-sized Jacobians, optimized for high throughput. This is achieved through vectorization and a Python pre-processor that symbolically prepares the stiff and sparse equation system, then writes the optimized code in Fortran syntax using source-file templates.

    Version : 0.1.0
    Author  : Jannes Klee
    Contact : jklee@astrophysik.uni-kiel.de
    License : GNU GPL v3

OdeVec was originally developed to solve chemical networks within the vectorized hydro-code [Fosite](https://github.com/tillense/fosite).

The implementation is based on the [BDF](https://en.wikipedia.org/wiki/Backward_differentiation_formula) (Backward Differentiation Formula) method described in [LSODE](https://computing.llnl.gov/sites/default/files/ODEPACK_pub2_u113855.pdf) from [ODEPACK](https://computing.llnl.gov/projects/odepack/publications) by A. C. Hindmarsh.

## Key Features

- **Vectorized ODE solving**: Processes many ODE systems of the same form simultaneously for improved throughput
- **Symbolic preprocessing**: Python-based symbolic computation of Jacobians and LU decompositions
- **Fortran backend**: Generates optimized Fortran code for high performance
- **Stiff equation support**: Specialized for stiff differential equations
- **Sparse matrix support**: Efficient handling of sparse Jacobians with CSC (Compressed Sparse Column) format
- **Multiple reordering algorithms**: CMK (Cuthill-McKee), INVERT, and FF (Fewest First) for optimal matrix structure
- **Built-in examples**: Robertson's problem, Primordial chemical network, and Oregonator system

## Requirements

### Python Dependencies
- Python 3.x
- See `requirements.txt` for specific package versions:
  - SymPy (for symbolic mathematics)
  - SciPy (for sparse matrix operations)
  - NumPy (for numerical computing)

### Fortran Compiler
- GNU Fortran (gfortran) or compatible Fortran compiler
- Tested with gfortran

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/jannesklee/odevec.git
   cd odevec
   ```

2. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```
   
   This will install all required packages:
   - `sympy` for symbolic mathematics
   - `scipy` for sparse matrix operations  
   - `numpy` for numerical computing

## Quick Start

### Running the Preprocessor

The preprocessor generates Fortran source files tailored to your specific ODE system:

```bash
./pre_odevec.py --nvector=32 --example="ROBER" && cd build
```

This command:
- Sets vector length to 32 (process 32 ODE systems simultaneously)
- Uses Robertson's example problem
- Generates preprocessed Fortran files in the `build/` directory

### Available Examples

- **ROBER**: Robertson's chemical reaction problem (3 species, very stiff)
- **PRIMORDIAL**: Primordial chemical network (16 species)
- **OREGO**: Oregonator chemical system (3 species)
- **KROME**: Custom networks from KROME package (requires setup file)

Example usage with different parameters:
```bash
# Primordial network with vector length 64
./pre_odevec.py --nvector=64 --example="PRIMORDIAL"

# Oregonator with CMK reordering
./pre_odevec.py --nvector=16 --example="OREGO" --ordering="CMK"

# Custom KROME network
./pre_odevec.py --nvector=32 --example="KROME" --krome_setupfile="path/to/setup.py"
```

### Compilation

Compile the generated Fortran files:

```bash
cd build
gfortran -O2 -c -g -cpp -fcheck=all -fno-range-check -Wall -fbacktrace odevec.f90 odevec_commons.f90 test.f90
gfortran -O2 -o rober *.o
```

### Execution

Run the compiled executable:

```bash
./rober
```

## Advanced Usage

### Command Line Options

The preprocessor supports various options for customization:

```bash
./pre_odevec.py --help
```

Key options:

- `--nvector`: Vector length (required, e.g., 32, 64, 128)
- `--example`: Predefined network (ROBER, PRIMORDIAL, OREGO, KROME)
- `--ordering`: Reordering algorithm (CMK, INVERT, FF, None)
- `--packaging`: Matrix packaging (DENSE, CSC)
- `--LUmethod`: LU decomposition method (1, 2, 3)
- `--dt_min`: Minimum timestep
- `--heatcool`: Enable cooling/heating for KROME networks
- `--sparsity_structure`: Print sparsity information

### Reordering Algorithms

OdeVec supports several matrix reordering algorithms to optimize the LU decomposition:

- **CMK**: Cuthill-McKee algorithm (reduces bandwidth)
- **INVERT**: Simple index inversion
- **FF**: Fewest First (orders by combined non-zeros)
- **None**: No reordering (default)

Example with CMK reordering:
```bash
./pre_odevec.py --nvector=32 --example="PRIMORDIAL" --ordering="CMK"
```

### Sparse Matrix Packaging

For efficient memory usage and computation:

- **DENSE**: Standard dense matrix format (default)
- **CSC**: Compressed Sparse Column format (for sparse systems)

Example with CSC packaging:
```bash
./pre_odevec.py --nvector=32 --example="ROBER" --packaging="CSC"
```

## Integration with KROME

OdeVec can be integrated with the [KROME](https://bitbucket.org/kkrome/krome/src/master/) package for astrochemical networks:

1. Generate a KROME setup file that dumps the network in Python syntax
2. Use the `--krome_setupfile` option:
   ```bash
   ./pre_odevec.py --nvector=32 --example="KROME" --krome_setupfile="krome_network.py"
   ```

3. For networks with cooling/heating, use the `--heatcool` flag:
   ```bash
   ./pre_odevec.py --nvector=32 --example="KROME" --krome_setupfile="krome_network.py" --heatcool=1
   ```

## Technical Details

### Preprocessing Workflow

1. **System Setup**: Load the ODE system (RHS and variables)
2. **Jacobian Computation**: Symbolic computation of the Jacobian matrix
3. **Matrix Reordering**: Apply selected reordering algorithm
4. **LU Decomposition**: Symbolic LU decomposition for preconditioning
5. **Code Generation**: Write optimized Fortran code with pragmas replaced

### Fortran Code Structure

The preprocessor generates two main files:
- `odevec.f90`: Main solver implementation
- `odevec_commons.f90`: Common data structures and utilities

### Performance Optimization

- **Vectorization**: Processes multiple ODE systems in parallel
- **Symbolic Precomputation**: Reduces runtime computational overhead
- **Sparse Matrix Handling**: Optimized memory usage for sparse systems
- **Preconditioning**: LU decomposition prepared symbolically

## Troubleshooting

### Common Issues

**Compilation Errors**:
- Ensure you have a compatible Fortran compiler installed
- Check that all generated files are in the build directory
- Verify that the vector length matches your system requirements

**Runtime Errors**:
- Check that the ODE system is properly defined
- Verify that initial conditions are physically reasonable
- Adjust tolerance parameters if convergence issues occur

**Performance Issues**:
- Experiment with different vector lengths
- Try different reordering algorithms
- Consider switching between DENSE and CSC packaging

## Contributing

Contributions are welcome! Please follow these guidelines:

1. Fork the repository
2. Create a feature branch
3. Implement your changes
4. Add appropriate tests
5. Submit a pull request

### Development Setup

```bash
git clone https://github.com/jannesklee/odevec.git
cd odevec
pip install -r requirements.txt
```

## License

OdeVec is licensed under the GNU General Public License v3. See [LICENSE.md](LICENSE.md) for details.

## Citation

If you use OdeVec in your research, please cite:

```bibtex
@misc{odevec,
  author = {Jannes Klee},
  title = {OdeVec: A vectorized ODE solver for stiff chemical networks},
  year = {2019},
  howpublished = {\url{https://github.com/jannesklee/odevec}}
}
```

## Contact

For questions, issues, or collaboration opportunities:
- **Author**: Jannes Klee
- **Email**: jklee@astrophysik.uni-kiel.de
- **GitHub**: [jannesklee/odevec](https://github.com/jannesklee/odevec)

## Acknowledgements

OdeVec builds upon the LSODE algorithm from ODEPACK and is inspired by the needs of the Fosite hydro-code for efficient chemical network integration.
