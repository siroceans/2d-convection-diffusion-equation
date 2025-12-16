# Solving the steady-state 2D convection-diffusion equation. 
***For more details on this project and the solution methodology please read the pdf of the full project report***

This project report presents the development of a solution to the 2D Steady-State Convection-Diffusion equation using Finite Difference techniques in MATLAB. The equation is solved in a square domain of side length 𝒍 = 𝟏, where the boundary conditions at each of the edges of the domain are given. The diffusion term is discretized using a central difference scheme for the second derivative. The discretization of the convection term is done using three separate techniques: first order UDS, second order UDS, and a central difference scheme; allowing us to do a qualitative comparison on the effects that each of these discretization techniques has on the approximated solution. The developed MATLAB program allows for a change in input values of the flow direction angle, diffusivity, and domain dimensions, 𝑵𝒙 and 𝑵𝒚.

### Methods Used
- The difussion term is discretized using a Central Difference Scheme (CDS)
- The Convective term is discretized using three different approaches: 
    1. 1st order Upwind Scheme (UD1)
    2. 2nd order Upwind Scheme (UD2)
    3. Central Difference Scheme (CDS)

### Results
The results show that discretizing the convection term using CDS shows unstable oscillatory results. UD1 and UD2 show stable results and a more detailed discussion, as well as a grid refinement study is presented in the full project report. 