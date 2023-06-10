# DynamicNetworks

DESCRIPTION
-----------
This is a Julia program (version 1.5.3) designed for the production of paper titled Strategy evolution on dynamic networks.

FILE
-----
<pre>
DynamicNetworks.jl      obtain the critical benefit-to-cost ratio (b/c)^* for dynamic networks that are comprised of any number of networks and involve any transition between them
                        provide detailed instruction of each function and illustrate the application with a simple example from Figure 2a in the main text
</pre>
                        
INSTALLATION
------------
<pre>
*) Download the open-source software Julia 1.5.3 or a newer version from https://julialang.org/downloads/.
*) Install Julia on your local PC. The expected installation time is approximately five minutes.
*) Download the open-source software VS Code 1.77.3 or a newer version from https://code.visualstudio.com/download.
*) Install VS Code on your local PC. The expected installation time is approximately ten minutes.
*) Install the Julia VS Code extension from https://code.visualstudio.com/docs/languages/julia. The expected installation time is approximately five minutes.
*) Run VS Code, add DynamicNetworks.jl, and install a solver package: 
   1) Type the following command in TEMINAL, and then press ENTER:
      julia
   2) Type the following command, and then press ENTER:  
      using Pkg; Pkg.add("IterativeSolvers");
*) Click "Run and Debug" to execute the program.
*) The expected outputs are
   1) Two critical benefit-to-cost ratios of examples from Figure 2a in the main text (by analytical computations), one with \alpha=0.5, the other with \alpha=0.8, corresponding to the two vertical lines in Figure 2b. The expected run time is approximately three minutes;
   2) A set of critical benefit-to-cost ratios presented in the inset of Figure 3a in the main text (by analytical computations). The expected run time is approximately five minutes.
</pre>

QUESTIONS
---------
For any question about this program, please contact
Dr. Qi Su, Email: qisu@sjtu.edu.cn, qisu1991@sas.upenn.edu
