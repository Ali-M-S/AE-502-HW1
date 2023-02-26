clear;clc
global mu
mu = 1.32712440018 * 10^11; % km^3 / s^2 (mu sun)

max_v1 = 50; % (km/s) maximum delta v for rendez-vous
max_v2 = 20; % (km/s) maximum delta v for fly-by

[err,data1] = LoadData('C:\Users\Asus\Desktop\AE 502 Advanced Orbital Mechanics\HW\HW1_Final\Earth2Oumuamua.txt');
  je =data1.jDate(:,:); % Julian dates earth
  jcalenum1=data1.calDate(:,:);
  T1 = datetime(jcalenum1,'ConvertFrom','datenum');
  re =data1.r(:,1:3);   % Position earth
  ve =data1.v(:,1:3);   % Velocity earth
[err,data2] = LoadData('C:\Users\Asus\Desktop\AE 502 Advanced Orbital Mechanics\HW\HW1_Final\Oumuamua_arrival.txt');
  jcom =data2.jDate(:,:);% Julian dates comet
  jcalenum2=data2.calDate(:,:);
  T2 = datetime(jcalenum2,'ConvertFrom','datenum');
  rcom =data2.r(:,1:3);  % Position comet
  vcom =data2.v(:,1:3);  % Velocity comet
  
  % total possible combinations of arrival and departures
  total = length(je)*length(jcom);
  
  fprintf('Departure Dates: %g', length(je))
  fprintf('\n Departure Dates: %g', length(jcom))
  fprintf('\n Total combinations: %g', total)
  
  v1_mat_short = zeros(length(jcom),length(je));
  v1_mat_long  = zeros(length(jcom),length(je));
  v2_mat_short = zeros(length(jcom),length(je));
  v2_mat_long  = zeros(length(jcom),length(je));
  tof = zeros(length(jcom),length(je));
  
  x = length(je);
  y = length(jcom);
  
for na = 1:y % length of arrival
      for nd = 1:x % length of departure         
          
%%%%%%%%%%%%%%%%%%%%%% Print current iteration
fprintf('\n departure = %d/365 | arrival= %d/549', nd, na)
        
      tof_seconds = (jcom(na,1) - je(nd,1))*24*60*60;  % Time of flight
   
      
       Orbital_Elements = r_v_2_O_E(re(nd,:), ve(nd,:), mu);
       if Orbital_Elements(4) <= 90 %prograde
           
         try % short
      [v_depart_short, v_arrive_short] = my_lambert(re(nd,:), rcom(na,:), tof_seconds); %Solve lambert short transfer problem
         catch
           v_depart_short = [1000 1000 1000];
           v_arrive_short = [1000 1000 1000];
           fprintf('\n no success short');
         end
         try % long
         [v_depart_long, v_arrive_long]   = my_lambert(re(nd,:), rcom(na,:), -tof_seconds); %Solve lambert long transfer problem
         catch
           v_depart_long = [1000 1000 1000];
           v_arrive_long = [1000 1000 1000];
           fprintf('\n no success long');
         end
         %%%%%%%%%%%%%%%%%%%%
       elseif Orbital_Elements(4) > 90
         try % short
      [v_depart_short, v_arrive_short] = my_lambert(re(nd,:), rcom(na,:), tof_seconds); %Solve lambert short transfer problem
         catch
           v_depart_short = [1000 1000 1000];
           v_arrive_short = [1000 1000 1000];
           fprintf('\n no success short');
         end
         try % long
         [v_depart_long, v_arrive_long]   = my_lambert(re(nd,:), rcom(na,:), -tof_seconds); %Solve lambert long transfer problem
         catch
           v_depart_long = [1000 1000 1000];
           v_arrive_long = [1000 1000 1000];
           fprintf('\n no success long');
         end
        end
       
        % departing velocity - earth velocity (Leaving earth)
         v1_short = norm(v_depart_short - ve(nd,:));
         v1_long = norm(v_depart_long - ve(nd,:));
         % arrival velocity - comet velocity (arriving to comet)
         v2_short = norm(v_arrive_short - vcom(na,:));
         v2_long = norm(v_arrive_long - vcom(na,:));
        
         tof(na,nd) = tof_seconds/24/60/60; % indexing current time of flight
              
        v1_mat_short(na,nd) =v1_short; %put in matrix the departing velocity
        v2_mat_short(na,nd) =v2_short; %put in matrix the arrivig velocity
        v1_mat_long(na,nd)  =v1_long; %put in matrix the departing velocity
        v2_mat_long(na,nd)  =v2_long; %put in matrix the arrivig velocity
      
      end
end

%%%%%%% Delta V for Rendevouz
delta_v_short = v1_mat_short + v2_mat_short; %(Total delta v)
delta_v_short(delta_v_short > max_v1) = NaN; % impose velocity limit
delta_v_long = v1_mat_long + v2_mat_long; %(Total delta v)
delta_v_long(delta_v_long > max_v1) = NaN; % impose velocity limit

%%%%%%% Delta V for fly-by
delta_v_short_f = v1_mat_short;
delta_v_long_f  = v1_mat_long;
delta_v_short_f(delta_v_short_f > max_v2) = NaN; % impose velocity limit
delta_v_long_f(delta_v_long_f > max_v2) = NaN; % impose velocity limit

%%%%%Plotting
figure;
subplot(1,2,1)
grid on;
hold on;
c = contour(delta_v_short,'ShowText','on','lineWidth',2);
h = contour(delta_v_long,'ShowText','on','lineWidth',2);
title('Rendevouz Extremal Field Map (Earth to Ouamuamua)');
xlabel('Departure Days after January 1st 2017');
ylabel('Arrival Days after August 1st 2017');
legend('Total Δv (km/s)', 'Location', 'best');
hold off

subplot(1,2,2)
delta_v_small = min(delta_v_short,delta_v_long);
f = surf(T1,T2,delta_v_small);
set(f,'LineStyle','none')
view (0,90)
colormap turbo
a1 = colorbar;
title(a1,'Total Δv (km/s)','bold')
title('Rendevouz Porkchop plot (Earth to Ouamuamua)');
xlabel('Earth departure date');
ylabel('Oumuamua arrival date');
                      
figure
subplot(1,2,1)
grid on;
hold on;
c2 = contour(delta_v_short_f,'ShowText','on','lineWidth',2);
h2 = contour(delta_v_long_f,'ShowText','on','lineWidth',2);
title('Fly-by Extremal Field Map (Earth to Ouamuamua)');
xlabel('Departure Days after January 1st 2017');
ylabel('Arrival Days after August 1st 2017');
legend('Total Δv (km/s)', 'Location', 'best');
hold off
subplot(1,2,2)
delta_v_small_f = min(delta_v_short_f,delta_v_long_f);
f2 = surf(T1,T2,delta_v_small_f);
set(f2,'LineStyle','none')
view (0,90)
colormap turbo
a2 = colorbar;
title(a2,'Total Δv (km/s)','bold')
title('Fly-by Porkchop Plot (Earth to Ouamuamua)');
xlabel('Earth departure date');
ylabel('Oumuamua arrival date');

    