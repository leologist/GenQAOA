# maxsat_solvers

This directory contains implementation of state-of-the-art classical algorithms for solving (weighted) (partial) MaxSAT.

## Solvers
1. `akmaxsat` - a complete solver for MaxSAT, guaranteed to find the optimal assignment. Winner of 2010
2. `CCLS2014` - an incomplete solver for MaxSAT, use heuristics such as stochastic local search (SLS) to find a good assignment.
Can repeat for multiple tries and keep the best one.
3. `CCLS_to_akmaxsat` - first calls `CCLS_2014` to get a good lower bound, and then pass the result to `akmaxsat` to find the best assignment.

## Compiling
To compile the solvers for use, simply run

`sh ./build_all_solvers.sh`

This will create a directory `bin` where all the binary executables are located

## Usage
Navigate to the `bin` directory and call the appropriate executable corresponding to your desired solver


`cd bin`

`./akmaxsat path_to_instance_file.txt`

or 

`./CCLS2014 path_to_instance_file.txt <seed> <max_tries>`

or

`./CCLS_to_akmaxsat path_to_instance_file.txt <lower_bound>`

Note
* Some good default choices: \<seed> = 1; <max_tries> = 10; <lower_bound> = 0
* It is only necessary to navigate to the `bin` directory if you want to use `CCLS_to_akmaxsat`,
which assumes that both `CCLS2014` and `akmaxsat` are in the current directory.
