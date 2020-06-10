# CDTool
CDTool is a set of CMake primitives to process time series data and peform the following tasks:
- Normalize a time series matrix using z-score normalization technique.
- Decompose a time series matrix to expose its rank.
- Recover the missing blocks inside time series data. 

The decomposition and the recovery tasks rely on an efficient matrix decomposition technique called the Centroid Decomposition (CD). Different variants of the CD technique have been integrated in the tool.

## Prerequisites

- Ubuntu 16 or 18 (including Ubuntu derivatives, e.g., Xubuntu).
```bash
    $ sudo apt install build-essential
    $ sudo apt install clang
```

## Build:


```bash
    $ make all
```

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

- **Tasks**: list of tasks
    - dec: performs decomposition
    - rec: performs recovery
    - norm: performs z-score normalization
- **Act**: list of actions
    - res: returns the result of the action
    - prec: returns the precision of the action (valid only for decomposition)
    - runtime: returns the runtime of the action (in microseconds)

<!---
The result is written into the file, decomposition output is written into three files with resp. suffixes (Loading matrix, Relevance matrix, Centroid Values)
-->

- **Input**: Input matrix of time series



## Optional commands:

- `-n [size]`, `-m [size]` makes the program read only specific number of rows/columns from the input matrix.
- `-k [size]` makes decomposition or recovery use this truncation parameter. 
    - Default truncation for decomposition: full matrix. 
    - Default truncation for recovery: use built-in auto-detection of truncation.
- `-output [filename]` redirects the output from a file decided by a program into a different one. 
- `-cdvar [cd-variant]` sets the execution of CD to use a differnt variant of the algorithm. See built-in help file for a list of variants.
- other parameters can be found in the built-in help file: `./cdec --help`

# Examples

- Decompose the matrix `example.txt` and store the result of its decomposition
```bash
    $ ./cdec -task dec -act res -input example.txt
```
You will find `example.txt.Load` with loading matrix,  `example.txt.Rel` with relevance matrix and `example.txt.Centroid` with centroid values that contain the result of the decomposition.

- Recover the missing values in `example_mis.txt` and store the runtime (in microseconds)
```bash
    $ ./cdec -task rec -act runtime -input example_mis.txt
```

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
