# CD_tool

## Prerequisites

- Ubuntu 16 or higher (including Ubuntu derivatives, e.g., Xubuntu).
- `sudo apt install build-essential`
- `sudo apt install clang`

## Building:

- `make all`

## Running:

- `cd cmake-build-debug`
- `./incCD [arguments]`
- Mandatory command-line arguments for the program:

 | -act | -test | -input |
 | -------- | -------- | -------- | 
 | dec    | out        | [filename] |
 | rec    | prec       | |
 | norm   | runtime    | |

Actions:
- dec - performs decomposition
- rec - performs recovery
- norm - performs normalization

Tests:
- out - provides output of the action
- prec - measures precision (only valid for decomposition)
- runtime - measures runtime

Result is written into the file, decomposition output is written into three files (Loading matrix, Relevance matrix, Centroid Values)

Optional commands:

- `-output [filename]` redirects the output from a file decided by a program into a different one. Keep in mind, in case of decomposition output your file name will be appended with `.Load`, `.Rel` and `.Centroid` suffixes for three outputs
- `-n [size]`, `-m [size]` makes program read only specifies number of rows/columns from the input matrix
- `-k [size]` makes decomposition or recovery use this truncation parameter. default (decomposition): run decomposition to maximum. default (recovery): use built-in auto-detection of truncation.
- other parameters can be found in built-in help file, which can be accessed by running `./incCD --help`
