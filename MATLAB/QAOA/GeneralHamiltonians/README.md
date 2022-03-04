## Explanation of the different versions of QAOA

1. `GenQAOA`, `GenQAOAGrad`, and `GenQAOAGradSmall` are for the QAOA with two general drivers HamC and HamB, where the objective Hamiltonian is HamC
2. `ExtQAOAGrad` is for a slight extension of the above QAOA with two general drivers HamC and HamB, where the objective Hamiltonian is separately specified by HamObj
3. `MultiQAOA` and `MultiQAOAGrad` are the most general form of the QAOA with arbitrary number of drivers, where the objective Hamiltonian is HamObj

