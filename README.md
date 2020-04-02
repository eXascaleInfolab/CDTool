# CD_tool

## Prerequisites

- Ubuntu 16 or higher (including Ubuntu derivatives, e.g., Xubuntu).
- `sudo apt install build-essential`
- `sudo apt install clang`

## Building:

- `make all`

## Running:

- `cd cmake-build-debug`
- `./cdec -task [arg] -act [arg] -input [filename]`
- Mandatory command-line arguments for the program:

 | -task | -act | -input |
 | -------- | -------- | -------- | 
 | dec    | res        | [filename] |
 | rec    | prec       | |
 | norm   | runtime    | |

Actions:
- dec - performs decomposition
- rec - performs recovery
- norm - performs z-score normalization

Tests:
- res - provides output of the action
- prec - measures precision (only valid for decomposition)
- runtime - measures runtime (in microseconds)

Result is written into the file, decomposition output is written into three files with resp. suffixes (Loading matrix, Relevance matrix, Centroid Values)

Optional commands:

- `--help` displays build-in help file with all commands and their usage.
- `-n [size]`, `-m [size]` makes program read only specifies number of rows/columns from the input matrix.
- `-k [size]` makes decomposition or recovery use this truncation parameter. default (decomposition): run decomposition to maximum. default (recovery): use built-in auto-detection of truncation.
- `-cdvar [cd-variant]` sets the execution of CD to use a non-default algorithm variant. See built-in help file for a list of variants.
- `-output [filename]` redirects the output from a file decided by a program into a different one. Keep in mind, in case of decomposition output your file name will be appended with `.Load`, `.Rel` and `.Centroid` suffixes for three outputs.
- other parameters can be found in built-in help file.

# Examples

- Run CD to decompose a matrix `example.txt` and store the result
```bash
    $ ./cdec -task dec -act res -input example.txt
```
You will find `example.txt.Load` with loading matrix,  `example.txt.Rel` with relevance matrix and `example.txt.Centroid` with centroid values that contain the result of the decomposition.

- Run CD to recover missing values in `example_mis.txt` and store the runtime
```bash
    $ ./cdec -task rec -act runtime -input example_mis.txt
```
You will find `example.txt.runtime` with the running time of the algorithm (in microseconds)

- Normalize the matrix `example.txt` and store the result in a custom file
```bash
    $ ./cdec -task norm -act res -input example.txt -output example_z-score_normalized.txt
```

- Run truncated decomposition on the matrix `example.txt` with truncation factor of 3
```bash
    $ ./cdec -task dec -act res -input example.txt -k 3
```

- Test the runtime of decomposition of the first 100, then 500 rows of the matrix `example.txt`.
```bash
    $ ./cdec -task dec -act runtime -input example.txt -n 100
    $ ./cdec -task dec -act runtime -input example.txt -n 500
```

- Recover missing values in the first 500 rows, 10 columns using truncation of 3 and store the reult in a custom file
```bash
    $ ./cdec -task rec -act res -input example_mis.txt -n 500 -m 10 -k 3 -output example_500_10_3_recovered.txt
```
