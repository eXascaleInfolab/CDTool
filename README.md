# CD_tool

## Prerequisites

- Ubuntu 16 or higher (including Ubuntu derivatives, e.g., Xubuntu).
- `sudo apt install build-essential`
- `sudo apt install clang`

## Build:

- `make all`

## Execution:

```bash
    $ cd cmake-build-debug
    $ ./cdec -task [arg] -act [arg] -input [filename]
```
- Command-line arguments for the program:

 | -task | -act | -input |
 | -------- | -------- | -------- | 
 | dec    | res        | [filename] |
 | rec    | prec       | |
 | norm   | runtime    | |

Tasks:
- dec - performs decomposition
- rec - performs recovery
- norm - performs z-score normalization

Actions:
- res - provides output of the action
- prec - measures precision (only valid for decomposition)
- runtime - measures runtime (in microseconds)


The result is written into the file, decomposition output is written into three files with resp. suffixes (Loading matrix, Relevance matrix, Centroid Values)

Optional commands:

- `-n [size]`, `-m [size]` makes the program read only specific number of rows/columns from the input matrix.
- `-k [size]` makes decomposition or recovery use this truncation parameter. 
    - Default truncation for decomposition: full matrix. 
    - Default truncation for recovery: use built-in auto-detection of truncation.
- `-output [filename]` redirects the output from a file decided by a program into a different one. 
- `-cdvar [cd-variant]` sets the execution of CD to use a differnt variant of the algorithm. See built-in help file for a list of variants.
- other parameters can be found in built-in help file: `./cdec --help`

# Examples

- Decompose the matrix `example.txt` and store the result of the decomposition
```bash
    $ ./cdec -task dec -act res -input example.txt
```
You will find `example.txt.Load` with loading matrix,  `example.txt.Rel` with relevance matrix and `example.txt.Centroid` with centroid values that contain the result of the decomposition.

- Recover the missing values in `example_mis.txt` and store the runtime
```bash
    $ ./cdec -task rec -act runtime -input example_mis.txt
```
You will find `example.txt.runtime` with the running time of the algorithm (in microseconds)

- Normalize the columns of matrix `example.txt` and store the result in a custom file
```bash
    $ ./cdec -task norm -act res -input example.txt -output example_norm.txt
```

- Decompose the matrix `example.txt` with a truncation factor 3 (keep only 3 dimensions)
```bash
    $ ./cdec -task dec -act res -input example.txt -k 3
```

- Decompose the first 100 of matrix `example.txt` and store the runtime.
```bash
    $ ./cdec -task dec -act runtime -input example.txt -n 100
```

- Recover the missing values of the first 500 rows, 10 columns of matrix `example_mis.txt` using a truncation of 3
```bash
    $ ./cdec -task rec -act res -input example_mis.txt -n 500 -m 10 -k 3
```
