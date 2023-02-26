function [c,s] = stumps_CS(z)
%{
This function evaluates the Stumpff functions C(z)& S(z) pg. 169
z = 1/a * chi^2
c is C(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% functions required: none
%}

if z > 0     % Ellipse
c = (1 - cos(sqrt(z)))/z;      
s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;

elseif z < 0 % Hyperbola
c = (cosh(sqrt(-z)) - 1)/(-z); 
s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;

elseif z == 0% Parabola
c = 1/2;                      
s = 1/6;
end
end

