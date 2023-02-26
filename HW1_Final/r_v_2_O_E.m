function Orbital_Elements = r_v_2_O_E(r0,v0,mu)
%{
% Function converts position and state vectors to classical orbital elements

mu - gravitational parameter (km^3/s^2)
r0 - position vector(km)
v0 - velocity vector (km)
r, v - the magnitudes of R and V
vr - 
H - the angular momentum vector (km^2/s)
h - the magnitude of H (km^2/s)
i - inclination of the orbit (deg)
N - the node line vector (km^2/s)
Omega - right ascension of the ascending node (deg)
e - eccentricity 
w - argument of perigee (deg)
theta - true anomaly (deg)
a - semimajor axis (km)

Orbital_Elements - vector of orbital elements [h e Omega i w theta a]
%%%%%%%%%%%%%%%%%%%%%%%% Functions required: None
%}
% –––––––––––––––––––––––––––––––––––––––––––––

r0_magnitude = norm(r0); %initial position magnitude
v0_magnitude = norm(v0); %initial velocity magnitude

vr = dot(r0,v0)/r0_magnitude; %Radial velocity component (km/s)

h_vector = cross(r0,v0); %angular momentum
h_magnitude = norm(h_vector);

I_unit = [1;0;0]; %I unit vector
J_unit = [0;1;0]; %J unit vector
K_unit = [0;0;1]; %K unit vector

node_vector = cross(K_unit , h_vector/h_magnitude); % node vector
node_magnitude = norm(node_vector);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% eccentricity (e)
e_vector = ((v0_magnitude)^2 / mu - 1/r0_magnitude)*r0 - (1/mu)*(dot(r0,v0))*v0;
e_magnitude = norm(e_vector);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inclination
i = acos(dot(h_vector/h_magnitude , K_unit)); %inclination
i_degree = 180/pi * i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Longitude of ascending node Ω (omega)
if dot( node_vector , J_unit) < 0 
    omega = 2*pi - acos((dot(node_vector , I_unit))/node_magnitude); 
else 
    omega = acos((dot(node_vector , I_unit))/node_magnitude); 
end 
omega_degree = 180/pi * omega;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Argument of periapse (w)
if dot( e_vector , K_unit) < 0
    w = 2*pi - acos(dot(node_vector , e_vector) / (node_magnitude * e_magnitude) );
else
    w = acos(dot(node_vector , e_vector) / (node_magnitude * e_magnitude) );
end
w_degree = 180/pi * w;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% True anomaly (f)
if vr < 0
    f = 2*pi - acos(dot(r0 , e_vector) / (r0_magnitude * e_magnitude) );
elseif vr >= 0
    f = acos(dot(r0 , e_vector) / (r0_magnitude * e_magnitude) );
end
f_degree = 180/pi * f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% semi-major axis (a)
if e_magnitude > 1 % Hyperbola
    a = h_magnitude^2/mu/(e_magnitude^2 - 1);
else
    a = 1/(2/r0_magnitude - (v0_magnitude)^2 / mu);
end
Orbital_Elements = [h_magnitude e_magnitude omega_degree i_degree w_degree f_degree a];
end

