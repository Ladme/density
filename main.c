
#include <unistd.h>
#include <groan.h>

const char VERSION[] = "v2022/09/19";

// frequency of printing during the calculation
static const int PROGRESS_FREQ = 10000;

/*
 * Gets user-defined grid dimension.
 */
int get_range(const char *optarg, float array[2])
{
    if (sscanf(optarg, "%f-%f", &array[0], &array[1]) != 2 && 
        sscanf(optarg, "%f - %f", &array[0], &array[1]) != 2 &&
        sscanf(optarg, "%f %f", &array[0], &array[1]) != 2) {
        fprintf(stderr, "Could not understand grid dimension specifier.\n");
        return 1;
    }

    return 0;
}

/*
 * Parses command line arguments.
 * Returns zero, if parsing has been successful. Else returns non-zero.
 */
int get_arguments(
        int argc, 
        char **argv,
        char **gro_file,
        char **xtc_file,
        char **ndx_file,
        char **output_file,
        char **atoms,
        float *array_dimx,
        float *array_dimy,
        float *array_dimz,
        int *grid_density)
{
    int gro_specified = 0, atoms_specified = 0;

    int opt = 0;
    while((opt = getopt(argc, argv, "c:f:n:o:s:x:y:z:d:h")) != -1) {
        switch (opt) {
        // help
        case 'h':
            return 1;
        // gro file to read
        case 'c':
            *gro_file = optarg;
            gro_specified = 1;
            break;
        // xtc file to read
        case 'f':
            *xtc_file = optarg;
            break;
        // ndx file to read
        case 'n':
            *ndx_file = optarg;
            break;
        // output file
        case 'o':
            *output_file = optarg;
            break;
        // selection to work with
        case 's':
            *atoms = optarg;
            atoms_specified = 1;
            break;
        // grid dimensions
        case 'x':
            if (get_range(optarg, array_dimx) != 0) return 1;
            break;
        case 'y':
            if (get_range(optarg, array_dimy) != 0) return 1;
            break;
        case 'z':
            if (get_range(optarg, array_dimz) != 0) return 1;
            break;
        // grid density
        case 'd':
            sscanf(optarg, "%d", grid_density);
            if (*grid_density <= 0) {
                fprintf(stderr, "Grid density must be > 0.\n");
                return 1;
            }
            break;
        default:
            //fprintf(stderr, "Unknown command line option: %c.\n", opt);
            return 1;
        }
    }

    if (!gro_specified || !atoms_specified) {
        fprintf(stderr, "Gro file and selection must always be supplied.\n");
        return 1;
    }
    return 0;
}

void print_usage(const char *program_name)
{
    printf("Usage: %s -c GRO_FILE -s SELECTION [OPTION]...\n", program_name);
    printf("\nOPTIONS\n");
    printf("-h               print this message and exit\n");
    printf("-c STRING        gro file to read\n");
    printf("-f STRING        xtc file to read (optional)\n");
    printf("-n STRING        ndx file to read (optional, default: index.ndx)\n");
    printf("-o STRING        output file (default: density.dx)\n");
    printf("-s STRING        selection of atoms which density shall be calculated\n");
    printf("-x FLOAT-FLOAT   grid dimension in x axis (default: box size from gro file)\n");
    printf("-y FLOAT-FLOAT   grid dimension in y axis (default: box size from gro file)\n");
    printf("-z FLOAT-FLOAT   grid dimension in z axis (default: box size from gro file)\n");
    printf("-d INTEGER       density of the grid used for calculation in points per nm (default: 10)\n");
    printf("\n");
}

/*
 * Prints parameters that the program will use for the calculation.
 */
void print_arguments(
        FILE *stream,
        const char *gro_file, 
        const char *xtc_file,
        const char *ndx_file,
        const char *output_file,
        const char *atoms,
        const float array_dimx[2],
        const float array_dimy[2],
        const float array_dimz[2],
        const int grid_density)
{
    fprintf(stream, "Parameters for Density calculation:\n");
    fprintf(stream, ">>> gro file:         %s\n", gro_file);
    if (xtc_file == NULL) fprintf(stream, ">>> xtc file:         ----\n");
    else fprintf(stream, ">>> xtc file:         %s\n", xtc_file);
    fprintf(stream, ">>> ndx file:         %s\n", ndx_file);
    fprintf(stream, ">>> output file:      %s\n", output_file);
    fprintf(stream, ">>> selection:        %s\n", atoms);
    fprintf(stream, ">>> grid dimensions:  x: %.1f - %.1f nm, y: %.1f - %.1f nm, z: %.1f - %.1f nm\n", 
            array_dimx[0], array_dimx[1], 
            array_dimy[0], array_dimy[1], 
            array_dimz[0], array_dimz[1]);
    fprintf(stream, ">>> grid density:     %d points per nm\n\n", grid_density);
}

/* 
 * Converts index of an array to coordinate.
 */
static inline float index2coor(int x, float minx, int grid_density)
{
    return (float) x / grid_density + minx;
}

/*
 * Converts coordinate to an index array.
 */
static inline size_t coor2index(float x, float minx, int grid_density)
{
    return (size_t) roundf((x - minx) * grid_density);
}

void calc_density_frame(
        const atom_selection_t *selection,
        size_t *grid,
        const size_t xsize,
        const size_t ysize,
        const size_t zsize,
        const float array_dimx[2],
        const float array_dimy[2],
        const float array_dimz[2],
        const int grid_density)
{
    for (size_t i = 0; i < selection->n_atoms; ++i) {
        atom_t *atom = selection->atoms[i];
        size_t xindex = coor2index(atom->position[0], array_dimx[0], grid_density);
        size_t yindex = coor2index(atom->position[1], array_dimy[0], grid_density);
        size_t zindex = coor2index(atom->position[2], array_dimz[0], grid_density);

        if (xindex >= xsize || yindex >= ysize || zindex >= zsize) continue;

        ++grid[zindex * xsize * ysize + yindex * xsize + xindex];
    }
}

void write_density(
        FILE *output,
        const size_t *grid,
        const size_t xsize,
        const size_t ysize,
        const size_t zsize,
        const float array_dimx[2],
        const float array_dimy[2],
        const float array_dimz[2],
        const int grid_density,
        const size_t n_samples)
{
    float rec_density = 10.0 * 1.0 / grid_density;

    fprintf(output, "object 1 class gridpositions counts %lu %lu %lu\n", xsize, ysize, zsize);
    fprintf(output, "origin %.5f %.5f %.5f\n", array_dimx[0] * 10.0, array_dimy[0] * 10.0, array_dimz[0] * 10.0);
    fprintf(output, "delta %.5f 0.00000 0.00000\n", rec_density);
    fprintf(output, "delta 0.00000 %.5f 0.00000\n", rec_density);
    fprintf(output, "delta 0.00000 0.00000 %.5f\n", rec_density);
    fprintf(output, "object 2 class gridconnections counts %lu %lu %lu\n", xsize, ysize, zsize);
    fprintf(output, "object 3 class array type \"double\" rank 0 items %lu data follows\n", xsize * ysize * zsize);

    size_t iterator = 0;
    for (size_t x = 0; x < xsize; ++x) {
        for (size_t y = 0; y < ysize; ++y) {
            for (size_t z = 0; z < zsize; ++z) {
                fprintf(output, "%.5f ", 10.0 * (float) grid[z * xsize * ysize + y * xsize + x] / n_samples);
                ++iterator;
                if (iterator % 3 == 0) fprintf(output, "\n");
            }
        }
    }

    if (iterator % 3 != 0) fprintf(output, "\n");

    fprintf(output, "attribute \"dep\" string \"positions\"\n");
    fprintf(output, "object \"density\" class field\n");
    fprintf(output, "component \"positions\" value 1\n");
    fprintf(output, "component \"connections\" value 2\n");
    fprintf(output, "component \"data\" value 3\n");
}

int main(int argc, char **argv)
{
    // get arguments
    char *gro_file = NULL;
    char *xtc_file = NULL;
    char *ndx_file = "index.ndx";
    char *output_file = "density.dx";
    char *atoms = NULL;
    float array_dimx[2] = {0.};
    float array_dimy[2] = {0.};
    float array_dimz[2] = {0.};
    int grid_density = 10;

    if (get_arguments(argc, argv, &gro_file, &xtc_file, &ndx_file, &output_file, &atoms, array_dimx, array_dimy, array_dimz, &grid_density) != 0) {
        print_usage(argv[0]);
        return 1;
    }

    printf("\n");

    // read gro file
    system_t *system = load_gro(gro_file);
    if (system == NULL) {
        return 1;
    }

    // if array dimensions were not set, get them from gro file
    if (array_dimx[0] == 0 && array_dimx[1] == 0) {
        array_dimx[1] = system->box[0];  
    }
    if (array_dimy[0] == 0 && array_dimy[1] == 0) {
        array_dimy[1] = system->box[1];
    }
    if (array_dimz[0] == 0 && array_dimz[1] == 0) {
        array_dimz[1] = system->box[2];
    }

    // check that the array dimensions don't have nonsensical values
    if (array_dimx[0] >= array_dimx[1] || array_dimy[0] >= array_dimy[1] || array_dimz[0] >= array_dimz[1]) {
        fprintf(stderr, "Nonsensical array dimensions.\n");
        free(system);
        return 1;
    }

    print_arguments(stdout, gro_file, xtc_file, ndx_file, output_file, atoms, array_dimx, array_dimy, array_dimz, grid_density);

    // try opening output file
    FILE *output = fopen(output_file, "w");
    if (output == NULL) {
        fprintf(stderr, "Output file could not be opened.\n");
        free(system);
        return 1;
    }

    // writing header for the output file
    fprintf(output, "# OpenDX density file written by density (C Density Calculator) %s\n", VERSION);
    fprintf(output, "# Command line: ");
    for (int i = 0; i < argc; ++i) {
        fprintf(output, "%s ", argv[i]);
    }
    fprintf(output, "\n");
    fprintf(output, "# Resulting data are in Ångström and can be loaded into VMD.\n");
    fprintf(output, "# Data is written in C array order: In grid[x,y,z] the axis z is fastest\n");
    fprintf(output, "# varying, then y, then finally x, i.e. z is the innermost loop.\n");

    // try reading ndx file (ignore if this fails)
    dict_t *ndx_groups = read_ndx(ndx_file, system);

    // select all atoms
    atom_selection_t *all = select_system(system);

    // select specified atoms
    atom_selection_t *selected = smart_select(all, atoms, ndx_groups);
    free(all);
    all = NULL;

    // check that the selection has been successful
    if (selected == NULL || selected->n_atoms == 0) {
        fprintf(stderr, "No atoms ('%s') found.\n", atoms);

        dict_destroy(ndx_groups);
        free(system);
        free(selected);
        fclose(output);
        return 1;
    }

    // prepare grid
    size_t xsize = (size_t) roundf( (array_dimx[1] - array_dimx[0]) * grid_density ) + 1;
    size_t ysize = (size_t) roundf( (array_dimy[1] - array_dimy[0]) * grid_density ) + 1;
    size_t zsize = (size_t) roundf( (array_dimz[1] - array_dimz[0]) * grid_density ) + 1;
    size_t n_points = xsize * ysize * zsize;

    size_t *grid = calloc(n_points, sizeof(size_t));

    if (grid == NULL) {
        fprintf(stderr, "Could not allocate memory (grid too large?)\n");
        dict_destroy(ndx_groups);
        free(system);
        free(selected);
        fclose(output);
        return 1;
    }

    // if xtc file is not provided, analyze only the gro file
    int return_code = 0;
    size_t n_frames = 0;
    if (xtc_file == NULL) {
        calc_density_frame(selected, grid, xsize, ysize, zsize, array_dimx, array_dimy, array_dimz, grid_density);

        write_density(output, grid, xsize, ysize, zsize, array_dimx, array_dimy, array_dimz, grid_density, 1);
        printf("File %s has been written.\n", output_file);

    } else {
        XDRFILE *xtc = xdrfile_open(xtc_file, "r");
        if (xtc == NULL) {
            fprintf(stderr, "File %s could not be read as an xtc file.\n", xtc_file);
            return_code = 1;
            goto program_end;
        }

        if (!validate_xtc(xtc_file, (int) system->n_atoms)) {
            fprintf(stderr, "Number of atoms in %s does not match %s.\n", xtc_file, gro_file);
            xdrfile_close(xtc);
            return_code = 1;
            goto program_end;
        }

        while (read_xtc_step(xtc, system) == 0) {
            // print info about the progress of reading
            if ((int) system->time % PROGRESS_FREQ == 0) {
                printf("Step: %d. Time: %.0f ps\r", system->step, system->time);
                fflush(stdout);
            }

            calc_density_frame(selected, grid, xsize, ysize, zsize, array_dimx, array_dimy, array_dimz, grid_density);
            ++n_frames;
        }

        
        write_density(output, grid, xsize, ysize, zsize, array_dimx, array_dimy, array_dimz, grid_density, n_frames);
        printf("File %s has been written.\n", output_file);
        xdrfile_close(xtc);
    }

    program_end:
    free(grid);
    dict_destroy(ndx_groups);
    free(system);
    free(selected);
    fclose(output);

    return return_code;

}