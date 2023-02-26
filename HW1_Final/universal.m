function chi_x = universal(dt, r0, vr0, a)
%{
This function uses Newton’s method to solve the universal
Kepler equation for the universal anomaly.
mu - gravitational parameter (km^3/s^2)
x - the universal anomaly (km^0.5)
dt - time since x = 0 (s)
ro - radial position (km) when x = 0
vro - radial velocity (km/s) when x = 0
a - reciprocal of the semimajor axis (1/km)
z - auxiliary variable (z = a*x^2)
c - value of Stumpff function C(z)
s - value of Stumpff function S(z)
itr - number of iterations for convergence
itrMax - maximum allowable number of iterations
%%%%%%%%%%%%%%%% Functions required: stumps_CS
%}
% ––––––––––––––––––––––––––––––––––––––––––––––
global mu
%error tolerance     limit of the number of iterations:
 error = 1.e-9;      itrMax = 1000;
%Initial guess for the universal anomaly 𝝌
chi_x = sqrt(mu)*abs(a)*dt;

%...Iterate on Equation 3.65 until until convergence occurs within
%...the error tolerance:
itr = 0;
F = 1;
while abs(F) > error && itr <= itrMax
itr = itr + 1;

[c,~] = stumps_CS(a*chi_x^2);
[~,s] = stumps_CS(a*chi_x^2);

F =  r0*vr0/sqrt(mu)*chi_x^2*c +...
     (1 - a*r0)*chi_x^3*s +...
     r0*chi_x - sqrt(mu)*dt;
 
dFdx =  r0*vr0/sqrt(mu)*chi_x*...
       (1 - a*chi_x^2*s) +...
       (1 - a*r0)*chi_x^2*c + r0;
   
chi_x_new = chi_x - F/dFdx;

%Update
chi_x = chi_x_new;

end
%...Deliver a value for chi_x, but report that itrMax was reached:
if itr > itrMax
fprintf('\n **No. iterations of Kepler’s equation = %g', itr);
fprintf('\n F/dFdx = %g\n', F/dFdx);
end
%