clear all; clc
%%%initial state vectors Distance: AU | Velocity: AU/day
%Time 
current_jdate = 2457754.5; % Julian date
current_Caldate = datetime(current_jdate,'convertfrom','juliandate','Format','dd-MMM-yyyy HH:mm:ss');
%Ouamuamua
r_oum1 = [3.515868886595499*10^-2,-3.162046390773074, 4.493983111703389];
v_oum1 = [-2.317577766980901*10^-3, 9.843360903693031*10^-3, -1.541856855538041*10^-2];
%Borisov
r_bor1 = [7.249472033259724, 14.61063037906177, 14.24274452216359];
v_bor1 = [-8.241709369476881*10^-3,-1.156219024581502*10^-2,-1.317135977481448*10^-2];
%Earth
r_e1 = [-1.796136509111975*10^-1, 9.667949206859814*10^-1,-3.668681017942158*10^-5];
v_e1 = [-1.720038360888334*10^-2,-3.211186197806460*10^-3, 7.927736735960840*10^-7];

global mu;
mu = 1.32712440018 * 10^11; % km^3 / s^2 (mu sun)

% Conversion to metric units to make use of r_v_2_0_E Function
day2sec = 24*60*60;   % 1 Day = 86,400 seconds
AU2km = 149597870; % 1 AU = 149597870 km (NASA JPL)

%Converted Initial states:
%Ouamuamua
r_oum = r_oum1*(AU2km);
v_oum = v_oum1*(AU2km/day2sec);
%Borisov
r_bor = r_bor1*(AU2km);
v_bor = v_bor1*(AU2km/day2sec);
%Earth
r_e = r_e1*(AU2km);
v_e = v_e1*(AU2km/day2sec);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Earth Initial state
Orbital_Elements1 = r_v_2_O_E(r_e,v_e,mu);
fprintf('\n\nEarth Orbital elements with respect to:')
fprintf('\n-Center body name:\t Sun')
fprintf('\n-Ref.Epoch:\t\t\t J2000')
fprintf('\n-Reference frame:\t Ecliptic of J2000')
fprintf('\n-At Calendar Date:\t %s ',current_Caldate)
fprintf('\n #Angular momentum\t\t = %g(km^2/s)', Orbital_Elements1(1))
fprintf('\n #Eccentricity\t\t\t = %g', Orbital_Elements1(2))
fprintf('\n #Inclination\t\t\t = %g(deg)', Orbital_Elements1(4))
fprintf('\n #RA of ascending node\t = %g(deg)', Orbital_Elements1(3))
fprintf('\n #Argument of perigee\t = %g(deg)', Orbital_Elements1(5))
fprintf('\n #True anomaly\t\t\t = %g(deg)', Orbital_Elements1(6))
fprintf('\n #Semimajor axis\t\t = %g(km)', Orbital_Elements1(7))
fprintf('\n #Periapse radius\t\t = %g(km)', Orbital_Elements1(1)^2/mu/(1 + Orbital_Elements1(2)))
%...If the orbit is an ellipse, output its period:
if Orbital_Elements1(2)<1
    T = 2*pi/sqrt(mu)*Orbital_Elements1(7)^1.5;
fprintf('\n Period:')
fprintf('\n Seconds = %g', T)
fprintf('\n Minutes = %g', T/60)
fprintf('\n Hours = %g', T/3600)
fprintf('\n Days = %g', T/24/3600)
else
fprintf('\n Object is interstellar!')  
end
fprintf('\n–––––––––––––––––––––––––––––––––––––––––––––––––––––\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ouamuamua Initial state
Orbital_Elements2 = r_v_2_O_E(r_oum,v_oum,mu);
fprintf('\n\n Ouamuamua Orbital elements with respect to:')
fprintf('\n-Center body name:\t Sun')
fprintf('\n-Reference frame:\t Ecliptic of J2000')
fprintf('\n-Ref.Epoch:  2458080.5 | 23-Nov-2017')
fprintf('\n-At Calendar Date:\t\t %s ',current_Caldate)
fprintf('\n #Angular momentum\t\t = %g(km^2/s)', Orbital_Elements2(1))
fprintf('\n #Eccentricity\t\t\t = %g', Orbital_Elements2(2))
fprintf('\n #Inclination\t\t\t = %g(deg)', Orbital_Elements2(4))
fprintf('\n #RA of ascending node\t = %g(deg)', Orbital_Elements2(3))
fprintf('\n #Argument of perigee\t = %g(deg)', Orbital_Elements2(5))
fprintf('\n #True anomaly\t\t\t = %g(deg)', Orbital_Elements2(6))
fprintf('\n #Semimajor axis\t\t = %g(km)', Orbital_Elements2(7))
fprintf('\n #Periapse radius\t\t = %g(km)', Orbital_Elements2(1)^2/mu/(1 + Orbital_Elements2(2)))
%...If the orbit is an ellipse, output its period:
if Orbital_Elements2(2)<1
    T = 2*pi/sqrt(mu)*Orbital_Elements2(7)^1.5;
fprintf('\n Period:')
fprintf('\n Seconds = %g', T)
fprintf('\n Minutes = %g', T/60)
fprintf('\n Hours = %g', T/3600)
fprintf('\n Days = %g', T/24/3600)
else
fprintf('\n Object is interstellar!')  
end
fprintf('\n–––––––––––––––––––––––––––––––––––––––––––––––––––––\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Borisov Initial state
Orbital_Elements3 = r_v_2_O_E(r_bor,v_bor,mu);
fprintf('\n\n Borisov Orbital elements with respect to:')
fprintf('\n-Center body name:\t Sun')
fprintf('\n-Reference frame:\t Ecliptic of J2000')
fprintf('\n-Ref.Epoch:  2459062.5 | 01-Aug-2020')
fprintf('\n-At Calendar Date:\t\t %s ',current_Caldate)
fprintf('\n #Angular momentum\t\t = %g(km^2/s)', Orbital_Elements3(1))
fprintf('\n #Eccentricity\t\t\t = %g', Orbital_Elements3(2))
fprintf('\n #Inclination\t\t\t = %g(deg)', Orbital_Elements3(4))
fprintf('\n #RA of ascending node\t = %g(deg)', Orbital_Elements3(3))
fprintf('\n #Argument of perigee\t = %g(deg)', Orbital_Elements3(5))
fprintf('\n #True anomaly\t\t\t = %g(deg)', Orbital_Elements3(6))
fprintf('\n #Semimajor axis\t\t = %g(km)', Orbital_Elements3(7))
fprintf('\n #Periapse radius\t\t = %g(km)', Orbital_Elements3(1)^2/mu/(1 + Orbital_Elements3(2)))
%...If the orbit is an ellipse, output its period:
if Orbital_Elements3(2)<1
    T = 2*pi/sqrt(mu)*Orbital_Elements3(7)^1.5;
fprintf('\n Period:')
fprintf('\n Seconds = %g', T)
fprintf('\n Minutes = %g', T/60)
fprintf('\n Hours = %g', T/3600)
fprintf('\n Days = %g', T/24/3600)
else
fprintf('\n Object is interstellar!')  
end
fprintf('\n–––––––––––––––––––––––––––––––––––––––––––––––––––––\n')
