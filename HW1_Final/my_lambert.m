function [V1, V2] = my_lambert(R1, R2, t)

global mu
%...Magnitudes of R1 and R2:
r1 = norm(R1);
r2 = norm(R2);


theta = acos(dot( R1 , R2 ) / (norm(R1) * norm(R2)));

 if t < 0     %then assume longer route. 
     t = t*-1;
theta = 2*pi - theta; 
 end 
 

A = sin(theta) * sqrt((r1*r2)/(1 - cos(theta))); %eqn 5.35 pg.241


for z1 = -5:0.2:1*10^8
[c,~] = stumps_CS(z1);
[~,s] = stumps_CS(z1);
y1 = r1 + r2 + A*(z1*s - 1)/sqrt(c);
F1 = (y1/c)^1.5*s + A*sqrt(y1) - sqrt(mu)*t;
if F1 >= 0
break;
end
end
% z = 0.3;
% F = 1;

% z-correction factor;
  z = z1 - 0.05;
  
itr = 0;
ratio = 1;
while (abs(ratio) > 1e-8) && itr <= 1000
    
itr = itr + 1;
[c,~] = stumps_CS(z);
[~,s] = stumps_CS(z);

y = r1 + r2 + A*(z*s - 1)/sqrt(c);

F = (y/c)^1.5*s + A*sqrt(y) - sqrt(mu)*t;

    if z == 0
    y_0 = r1 + r2 + A*(z*(1/6 - z/120) - 1)/(1/2 - z/24);    
    dFdz = sqrt(2)/40*y_0^1.5 + A/8*(sqrt(y_0) + A*sqrt(1/2/y_0));
    else
    dFdz = (y/c)^1.5*(1/2/z*(c - 3*s/2/c) ...
         + 3*s^2/4/c) + A/8*(3*s/c*sqrt(y) ...
         + A*sqrt(c/y));
    end
     ratio = F/dFdz;
z_new = z - ratio;
% Update
z = z_new;

end

%...Equation 5.46a:
f = 1 - y/r1;
%...Equation 5.46b:
g = A*sqrt(y/mu);
%...Equation 5.46d:
gdot = 1 - y/r2;
%...Equation 5.28:
V1 = 1/g*(R2 - f*R1);
%...Equation 5.29:
V2 = 1/g*(gdot*R2 - R1);

end
