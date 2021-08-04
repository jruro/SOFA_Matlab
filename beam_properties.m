%
%   properties of the beam material
%


L = 0.732;
E = 500e+9;
height = 0.015; 
width = 0.00302;
Area = width*height;
Inertia = (height*(width^3))/12;
rho = 2700;
mass = rho*Area*L;
M = 0.0464;
mass_eq = 0.23*mass + M;


k = (3*E*Inertia)/(L^3)