% Create matrices M and N of Example 9 in M. HUHTANEN, V. KOTILA
% AND P. UUSITALO: FASTEST QUOTIENT ITERATION FOR GENERALIZED
% SELF-ADJOINT EIGENVALUE PROBLEMS
% with Matlab's PDE Toolbox.

% The dimensions of the Z-shape:
R=15;
L=3;

% Geometry description matrix:
% Row 1:	2 (indicates a polygon)
% Row 2:	Number of line segments n
% Rows 3 through 3+n-1:	x-coordinate of edge starting points
% Rows 3+n through 2*n+2:	y-coordinate of edge starting points

gd = [2
8
-R
-1/2
-1/2
R
R
1/2
1/2
-R
L/2-1 
L/2-1
-L/2
-L/2
1-L/2
1-L/2
L/2
L/2
];

dl = decsg(gd);

% Create a pde-model and import the geometry into it

model = createpde;
geometryFromEdges(model,dl); % geometryFromEdges for 2-D
applyBoundaryCondition(model,"dirichlet","Edge",1:8,"u",0);
specifyCoefficients(model,"d",1,"c",1,"a",0,"f",0,"m",0);
generateMesh(model,"GeometricOrder","quadratic","Hmax",0.05,"Hmin",0.0025,"Hgrad",1.3);


% For comparison, the results-object contains the least eigenvalues with  
% the corresponding eigenvectors, solved with pde-toolbox.
results = solvepdeeig(model,[0,10]); 

% Uncomment these two lines to view the mesh for the model:
% figure
% pdemesh(model)

% Use assembleFEMatrices to export the global finite element mass and stiffness
% matrices with boundary conditions imposed using nullspace approach.

FEMn = assembleFEMatrices(model,"nullspace")
M = FEMn.Kc;
N = FEMn.M;



