# density: Calculating Density in Gromacs Simulations

Calculates density of selected particles in a simulation box.

## Dependencies

`density` requires you to have groan library installed. You can get groan from [here](https://github.com/Ladme/groan). See also the [installation instructions](https://github.com/Ladme/groan#installing) for groan.

## Installation

1) Run `make groan=PATH_TO_GROAN` to create a binary file `density` that you can place wherever you want. `PATH_TO_GROAN` is a path to the directory containing groan library (containing `groan.h` and `libgroan.a`).
2) (Optional) Run `make install` to copy the the binary file `density` into `${HOME}/.local/bin`.

## Options

```
Usage: density -c GRO_FILE -s SELECTION [OPTION]...

OPTIONS
-h               print this message and exit
-c STRING        gro file to read
-f STRING        xtc file to read (optional)
-n STRING        ndx file to read (optional, default: index.ndx)
-o STRING        output file (default: density.dx)
-s STRING        selection of atoms which density shall be calculated
-x FLOAT-FLOAT   grid dimension in x axis (default: box size from gro file)
-y FLOAT-FLOAT   grid dimension in y axis (default: box size from gro file)
-z FLOAT-FLOAT   grid dimension in z axis (default: box size from gro file)
-d INTEGER       density of the grid used for calculation in points per nm (default: 10)
```

Use [groan selection language](https://github.com/Ladme/groan#groan-selection-language) to select the atoms for density calculation.

## Example usage

```
density -c md.gro -f md.xtc -s "resname POPC" -d 5
```

The program will analyze trajectory `md.xtc` and calculate the average density of atoms with residue name POPC (flag `-s`). Information about the atoms will be read from `md.gro` but density will not be calculated for the snapshot in the `gro` file. Density will be calculated using a grid with grid points separated by 0.2 nm (flag `-d`). The resulting density will be written into an OpenDX file `density.dx` (default option of the flag `-o`) and can be visualized using VMD or other software.

## Limitations

`density` calculates density in a box which means that it can not treat periodic boundary conditions.

Only tested on Linux. Probably will not work on anything that is not UNIX-like.
