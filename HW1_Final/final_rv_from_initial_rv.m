function [Rf,Vf] = final_rv_from_initial_rv(R0, V0, t)
%{
This is the Universal Variable two-body orbit propagator function

mu - the gravitational parameter (km^3/s^2)
a - reciprocal of the semimajor axis (1/km)
t - the time elapsed since ro (s)
x - the universal anomaly after time t (km^0.5)
f - the Lagrange f coefficient (dimensionless)
g - the Lagrange g coefficient (s)
%%%%%%%%%%%%%%% Functions required: stumps_CS and universal
%}
% ––––––––––––––––––––––––––––––––––––––––––––––
global mu

% Given initial Position and velocity vectors
r0 = norm(R0);
v0 = norm(V0);

% Radial velocity component:
vr0 = dot(R0, V0)/r0;

% Inverse of semimajor axis pg 175:
alpha = 2/r0 - v0^2/mu;

% Universal anomaly 𝝌:
x = universal(t, r0, vr0, alpha);

% f and g functions from Universal anomaly and Stumpff function:
z = alpha*x^2;
[c,~] = stumps_CS(z);
[~,s] = stumps_CS(z);

f = 1 - x^2/r0*c;         %...Equation 3.69a Pg. 174
g = t - 1/sqrt(mu)*x^3*s; %...Equation 3.69b Pg. 174

 % Final position vector:
  Rf = f*R0 + g*V0;

Rf_mag = norm(Rf); % magnitude of R


%%%%%%%%%%%%%%%%%%%%%%%%%% F & G Dots 

fdot = sqrt(mu)/(Rf_mag*r0) * (alpha*(x^3)*s - x); %...Equation 3.69c Pg. 174
gdot = 1 - (x^2/Rf_mag)*c;                         %...Equation 3.69d Pg. 174

 % Final velocity vector:
   Vf = fdot*R0 + gdot*V0;




end