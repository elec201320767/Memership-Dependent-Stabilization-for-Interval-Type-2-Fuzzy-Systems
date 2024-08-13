This code is for the paper "H∞ control for interval type-2 Takagi–Sugeno fuzzy systems via the membership-quadratic framework" in Information Sciences.

The author of the code is KyungSoo Kim _kyungsoo@postech.ac.kr_.

The code considers Example 1 of the paper
Here are the brief descriptions:

  Cor1.m provides 
  
            a graphical illustration of distributions of lower, embedded, and upper membership functions in 3-dimensional space - see file "Membership distributions in 3-dimensional spaces.png"
  
            LMI conditions and their solution
           
            the Simulink 'Cor1_Ex1' result           
           
  Cor1_deflmivar_cvx.m provides
  
             LMI declaration
  
  Cor1_Ex1.slx provides
  
             overall system configuration in MATLAB Simulink
  
  functions provides 
  
             tools for experiment

[Requirement]
The convex programming CVX should be set in your MATLAB
CVX can be installed from https://cvxr.com/cvx/doc/install.html

[Environment]
MATLAB R2020b 
CVX MATLAB
