# maxsat_solvers

This directory contains implementations of state-of-the-art classical algorithms for solving (weighted) (partial) MaxSAT.

## MaxSAT Problem
MaxSAT is the following problem: Given a CNF formula, which is a conjunction of many clauses (i.e., an AND of ORs), find an assignment that satisfies the maximum number of clauses possible.
An example CNF: (x_1 OR ~x_2) AND (x_2 OR x_3 OR ~x4) AND (~x_1 OR  ~x_3) AND ...

Weighted MaxSAT assigns a weight to each clause, and find an assignment that maximizes the total weight of clauses satisfied.

Weighted Partial MaxSAT is like Weighted MaxSAT, except some clauses are hard (must be satisfied in the output assignment) and others are soft (may be violated in the assignment).

Each instance of these problems may be written in the DIMACS format, according to specification at [http://www.maxsat.udl.cat/14/requirements/](http://www.maxsat.udl.cat/14/requirements/).


## Solvers
1. `akmaxsat` - a complete solver for MaxSAT, guaranteed to find the optimal assignment. Winner of [2010 MaxSAT Evaluation competition](http://www.maxsat.udl.cat/10/results/#wms-random).
2. `CCLS2014` - an incomplete solver for MaxSAT, use heuristics such as stochastic local search (SLS) to find a good assignment.
Can repeat for multiple tries and keep the best one.
3. `CCLS_to_akmaxsat` - first calls `CCLS_2014` to get a good lower bound, and then pass the result to `akmaxsat` to find the best assignment. Winner of various categories in weighted MaxSAT in [2014 MaxSAT Evaluation Competition](http://www.maxsat.udl.cat/14/results/index.html#wpms-random).

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
