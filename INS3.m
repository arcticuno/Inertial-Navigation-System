%% EHG446 INS 3  
clear all; close all; clc; 

%% Task One

% open all data
load('sensors_noisy.mat');
load('sensors_clean.mat');
load('att_truth.mat');
load('pos_llh.mat');

%extract clean & noisy body rate data
ax = sensors_clean(5,:);            %gyroscope x
ay = sensors_clean(6,:);            %gyroscope y 
az = sensors_clean(7,:);            %gyroscope z 

axn = sensors_noisy(5,:);            %gyroscope x
ayn = sensors_noisy(6,:);            %gyroscope y 
azn = sensors_noisy(7,:);            %gyroscope z 

%create plot vector
t0 = length(sensors_clean);

% Initialise all given values from INS 2 Task Sheet
initPosition = deg2rad([-27.41, 153.120, 300]);         % u latitude & u Lon in rad
h = 300;                                                    % initial height of aircraft 
u_v = [30.6 -12.7 -10.54];                                  % Initial Vleocity
initial_attitude = [0, 0, 5.8407];                          % Initial Attitude
g = [0;0;9.79];                                             % gravitational acceleration
omega = 7.292115e-5;                                          % Omega
r0 = 6378137;                                               % earth radius in meters
nr0 = [0; 0; r0];                                           %earth radius in matrix form


%% Task Two
%Using Initial Values from Task Sheet 

phi = 0; 
theta = 0; 
psi = 5.8407;

%As per Chapter 3 pg(45-46)
a = (cos(psi/2)*cos(theta/2)*cos(phi/2)) + (sin(psi/2)*sin(theta/2)*sin(phi/2)); 
b = (cos(psi/2)*cos(theta/2)*sin(phi/2)) - (sin(psi/2)*sin(theta/2)*cos(phi/2)); 
c = (cos(psi/2)*sin(theta/2)*cos(phi/2)) + (sin(psi/2)*cos(theta/2)*sin(phi/2)); 
d = (sin(psi/2)*cos(theta/2)*cos(phi/2)) - (cos(psi/2)*sin(theta/2)*sin(phi/2)); 

%Converting now to single Quaternion 
q1= [a b c d];
q1 = transpose(q1);

%Initialising Cn
DCM = [a^2+b^2-c^2-d^2 2*(b*c-a*d) 2*(b*d+a*c);
       2*(b*c+a*d) a^2-b^2+c^2-d^2 2*(c*d-a*b);
       2*(b*d-a*c) 2*(c*d+a*b) a^2-b^2-c^2+d^2];

   
%initialize wn_b from task sheet 
wn_b = zeros(3,t0);
wn_b(:,1) = [30.6 -12.7 -10.54]';

%initialise qdot matrix and new quaternion
q_dot = zeros(4,t0);
quaternion2 = zeros(4,t0);
quaternion2(:,1)= q1;

%Clean WIB from slides for clean 
wib_b= sensors_clean(5:7,:);

%initalize pnb as per DemoV3 with qdot
Pnbi = zeros(1,t0);
Pnb= [Pnbi; wn_b];                      %put extra rows into wnb matrix

%initalize final update 1
update_Exp1 = zeros(3,t0);
update_Exp1(:,1) = [initial_attitude]';

%initalize final update 2
T= length(att_truth);                   %Use new time vector 

%Original WEN & WIE without correction 
wen_n = [0 0 0]';                 

wie_n =[omega*cos(initPosition(1)) 0 -omega*sin(initPosition(1))]';

update_Exp2=zeros(3,T); 
update_Exp2(:,1) = [-27.41*pi/180 153.120*pi/180 300];

%Original GLN as per Slides
gl_n = zeros(3,(t0));
gl_n(:,1) = 9-wie_n.^2*rad2deg(initPosition(:,1));

fb_c = sensors_clean(2:4,:);                 %Clean  Body Axis Accel                  
fb_n = sensors_noisy(2:4,:);                 %Noisy Body Axis Accel

%initialize v_n from task sheet 
v_n  = zeros(3,T);
v_n(:,1) = u_v';

%Initiate sample time 
sampt = 0.01;               

%Initialise zeros for LAT,LON, FN & VDOT
latdot=zeros(3,T);             
londot= zeros(3,T);  
f_n = zeros(3,T);
v_dot= zeros(3,T);

%Updating to new values 
for k =2:12001
    %Task Three & Four --- Incoperated in loop
    
    %New Equations with corrections as per slides
    wen_n = [ v_n(2,k-1)/(r0 + update_Exp2(3,k-1)); -v_n(1,k-1)/(r0 + update_Exp2(3,k-1)); (-v_n(2,k-1)*tan(update_Exp2(1,k-1)))/(r0 + update_Exp2(3,k-1)) ];
    
    wie_n =[ omega*cos(update_Exp2(1,k-1)); 0; -omega*sin(update_Exp2(1,k-1)) ];

    gl_n = g - ((omega)^2*(r0 + update_Exp2(3,k-1))/2)*[sin(2*update_Exp2(1,k-1)); 0; 1+cos(2*update_Exp2(1,k-1))];
    
    %Obtain Cn_b (From slide transpose DCM) 
    
    DCM = [a^2+b^2-c^2-d^2 2*(b*c-a*d) 2*(b*d+a*c);
           2*(b*c+a*d) a^2-b^2+c^2-d^2 2*(c*d-a*b);
           2*(b*d-a*c) 2*(c*d+a*b) a^2-b^2-c^2+d^2];
       
    Cn_b = transpose(DCM);
    
    %As per demov3 sheet update wn_b
    wn_b(:,k)=wib_b(:,k)-Cn_b*(wie_n+wen_n);
    
    %As per slides find Q & qdot 
    Q_q = [a -b -c -d;
          b a -d c;
          c d a -b;
          d -c b a];
    
    %As per demov3 sheet update Pn_b & qdot 
    Pnb= [0, wn_b(:,k)']';
          
    q_dot(:,k-1) = 0.5.* Q_q* Pnb;
    
    %Normalising Quaternion
    quaternion2(:,k) = quaternion2(:,k-1) + sampt * q_dot(:,k-1);
    quaternion2(:,k) = quaternion2(:,k)/norm(quaternion2(:,k));
   
   %Updated values 
    a = quaternion2(1,k);
    b = quaternion2(2,k);
    c = quaternion2(3,k);
    d = quaternion2(4,k);
    
    %Computing Euler Expressions after update 
    % As per INS LEC 1 Slide 19
    tx = atan2(2*(c*d+a*b),(a^2-b^2-c^2+d^2));
    sy = asin(-2*(b*d-a*c));
    tz = atan2(2*(b*c+a*d),(a^2+b^2-c^2-d^2));
    
    %Update Complete Final Expression1 
    update_Exp1(:,k) = [tx; sy; tz]; 
         
    %New Vlaues to initialise Lab Nine
    %Initialising DCM using updated a,b,c&d
    
    DCM_new = [a^2+b^2-c^2-d^2 2*(b*c-a*d) 2*(b*d+a*c);
              2*(b*c+a*d) a^2-b^2+c^2-d^2 2*(c*d-a*b);
              2*(b*d-a*c) 2*(c*d+a*b) a^2-b^2-c^2+d^2];
   
    %Clean FN
    f_n(:,k-1) = DCM_new * fb_c(:,k-1);
   
    %As per powerpoint
    v_dot(:,k-1) = f_n(:,k-1) - cross((2*wie_n+wen_n),v_n(:,k-1)) + gl_n;
    
    %Final Updadte 
    v_n(:,k) = v_n(:,k-1)+sampt*v_dot(:,k-1);

    % updating the change in lattitude propogation
    latdot= v_n(1,k)/(r0+update_Exp2(3,k-1)); 
    
    % updating the change in longditude propogation
    londot= v_n(2,k)*sec(update_Exp2(1,k-1))/(r0+update_Exp2(3,k-1));

    % updating the change in  height propogation
    hdot= -v_n(3,k); 
    
    % Determining the new lattitude and longditude positions
    lat=update_Exp2(1,k-1)+sampt*latdot;            % updated lattitude
    lon= update_Exp2(2,k-1)+sampt*londot;           % updated the longditude
    h= update_Exp2(3,k-1)+ sampt*hdot;              % updating the height
    
    %Update Complete expression 2 
    update_Exp2(:,k)= [lat; lon; h];                % list of updated positions 
  
    
end 

%Task 5,6&7
%Plotting all Relvant Information 

% Obtaining now the Truth data, not used in task 3 
new_update_Exp2= zeros(1,length(pos_llh)); 
new_update_Exp1= zeros(1,length(pos_llh)); 

%Creating a loop to comute solutions within
for k= 1:3
    
    new_update_Exp2(k,:)= update_Exp2(k,1:100:T);      
    new_update_Exp1(k,:)= update_Exp1(k,1:100:T);      
    
end

t= 1:length(new_update_Exp2);               %New time vector 

%Loading in Truth 
roll_truth= att_truth(2,:);
pitch_truth= att_truth(3,:); 
yaw_truth= att_truth(4,:);

%Plots

figure(1)
hold on; 
plot(new_update_Exp2(1,:))
plot(pos_llh(2,:));
hold off 
title('Clean Data')
xlabel('Step')
ylabel('Latitude(Deg)')
legend('Solution','Truth')

figure(2)
hold on; 
plot(new_update_Exp2(2,:))
plot(pos_llh(3,:));
hold off 
title('Clean Data')
xlabel('Step (Sec)')
ylabel('Longitude (Deg)')
legend('Solution','Truth')

figure(3)
hold on; 
plot(new_update_Exp2(3,:))
plot(pos_llh(4,:));
hold off 
title('Clean Data')
xlabel('Step (Sec)')
ylabel('Height (m)')
legend('Solution','Truth')




%% Task Eight 

%Nosiy WIB from slides for clean 
wib_b= sensors_noisy(5:7,:);

%Updating to new values 
for k =2:12001
    %Task Three & Four --- Incoperated in loop
    
    %New Equations with corrections as per slides
    wen_n = [ v_n(2,k-1)/(r0 + update_Exp2(3,k-1)); -v_n(1,k-1)/(r0 + update_Exp2(3,k-1)); (-v_n(2,k-1)*tan(update_Exp2(1,k-1)))/(r0 + update_Exp2(3,k-1)) ];
    
    wie_n =[ omega*cos(update_Exp2(1,k-1)); 0; -omega*sin(update_Exp2(1,k-1)) ];

    gl_n = g - ((omega)^2*(r0 + update_Exp2(3,k-1))/2)*[sin(2*update_Exp2(1,k-1)); 0; 1+cos(2*update_Exp2(1,k-1))];
    
    %Obtain Cn_b (From slide transpose DCM) 
    
    DCM = [a^2+b^2-c^2-d^2 2*(b*c-a*d) 2*(b*d+a*c);
           2*(b*c+a*d) a^2-b^2+c^2-d^2 2*(c*d-a*b);
           2*(b*d-a*c) 2*(c*d+a*b) a^2-b^2-c^2+d^2];
       
    Cn_b = transpose(DCM);
    
    %As per demov3 sheet update wn_b
    wn_b(:,k)=wib_b(:,k)-Cn_b*(wie_n+wen_n);
    
    %As per slides find Q & qdot 
    Q_q = [a -b -c -d;
          b a -d c;
          c d a -b;
          d -c b a];
    
    %As per demov3 sheet update Pn_b & qdot 
    Pnb= [0, wn_b(:,k)']';
          
    q_dot(:,k-1) = 0.5.* Q_q* Pnb;
    
    %Normalising Quaternion
    quaternion2(:,k) = quaternion2(:,k-1) + sampt * q_dot(:,k-1);
    quaternion2(:,k) = quaternion2(:,k)/norm(quaternion2(:,k));
   
   %Updated values 
    a = quaternion2(1,k);
    b = quaternion2(2,k);
    c = quaternion2(3,k);
    d = quaternion2(4,k);
    
    %Computing Euler Expressions after update 
    % As per INS LEC 1 Slide 19
    tx = atan2(2*(c*d+a*b),(a^2-b^2-c^2+d^2));
    sy = asin(-2*(b*d-a*c));
    tz = atan2(2*(b*c+a*d),(a^2+b^2-c^2-d^2));
    
    %Update Complete Final Expression1 
    update_Exp1(:,k) = [tx; sy; tz]; 
         
    %New Vlaues to initialise Lab Nine
    %Initialising DCM using updated a,b,c&d
    
    DCM_new = [a^2+b^2-c^2-d^2 2*(b*c-a*d) 2*(b*d+a*c);
              2*(b*c+a*d) a^2-b^2+c^2-d^2 2*(c*d-a*b);
              2*(b*d-a*c) 2*(c*d+a*b) a^2-b^2-c^2+d^2];
   
    %Clean FN
    f_n(:,k-1) = DCM_new * fb_n(:,k-1);
   
    %As per powerpoint
    v_dot(:,k-1) = f_n(:,k-1) - cross((2*wie_n+wen_n),v_n(:,k-1)) + gl_n;
    
    %Final Updadte 
    v_n(:,k) = v_n(:,k-1)+sampt*v_dot(:,k-1);

    % updating the change in lattitude propogation
    latdot= v_n(1,k)/(r0+update_Exp2(3,k-1)); 
    
    % updating the change in longditude propogation
    londot= v_n(2,k)*sec(update_Exp2(1,k-1))/(r0+update_Exp2(3,k-1));

    % updating the change in  height propogation
    hdot= -v_n(3,k); 
    
    % Determining the new lattitude and longditude positions
    lat=update_Exp2(1,k-1)+sampt*latdot;            % updated lattitude
    lon= update_Exp2(2,k-1)+sampt*londot;           % updated the longditude
    h= update_Exp2(3,k-1)+ sampt*hdot;              % updating the height
    
    %Update Complete expression 2 
    update_Exp2(:,k)= [lat; lon; h];                % list of updated positions 
  
    
end 

%Task 5,6&7
%Plotting all Relvant Information 

% Obtaining now the Truth data, not used in task 3 
new_update_Exp2= zeros(1,length(pos_llh)); 
new_update_Exp1= zeros(1,length(pos_llh)); 

%Creating a loop to comute solutions within
for k= 1:3
    
    new_update_Exp2(k,:)= update_Exp2(k,1:100:T);      
    new_update_Exp1(k,:)= update_Exp1(k,1:100:T);      
    
end

t= 1:length(new_update_Exp2);               %New time vector 

%Loading in Truth 
roll_truth= att_truth(2,:);
pitch_truth= att_truth(3,:); 
yaw_truth= att_truth(4,:);

%Plots for Noisy 

figure(8)
hold on; 
plot(new_update_Exp2(1,:))
plot(pos_llh(2,:));
hold off 
title('Noisy Data')
xlabel('Step')
ylabel('Latitude(Deg)')
legend('Solution','Truth')

figure(9)
hold on; 
plot(new_update_Exp2(2,:))
plot(pos_llh(3,:));
hold off 
title('Noisy Data')
xlabel('Step (Sec)')
ylabel('Longitude (Deg)')
legend('Solution','Truth')

figure(10)
hold on; 
plot(new_update_Exp2(3,:))
plot(pos_llh(4,:));
hold off 
title('Noisy Data')
xlabel('Step (Sec)')
ylabel('Height (m)')
legend('Solution','Truth')
%
figure(11)
subplot(3,1,1);
hold on
plot(update_Exp1(1,:));
plot(roll_truth);
hold off
title('Noisy Roll')
xlabel('Time (Sec) ')
ylabel('Roll (Phi)')
legend('Solution','Truth')

subplot(3,1,2);
hold on
plot(update_Exp1(2,:));
plot(pitch_truth);
hold off
title('Noisy Pitch')
xlabel('Time (Sec)')
ylabel('Roll (Theta)')
legend('Solution','Truth')


subplot(3,1,3);
hold on
plot(update_Exp1(3,:));
plot(yaw_truth);
hold off
title('Noisy Yaw')
xlabel('Time (Sec)')
ylabel('Roll (Psi)')
legend('Solution','Truth')

%3D Plot
figure(14)
plot3(update_Exp2(1,1:100:T),update_Exp2(2,1:100:T),update_Exp2(3,1:100:T))
hold on
plot3(pos_llh(2,:),pos_llh(3,:),pos_llh(4,:))
title('Noisy V True 3D Rep');
xlabel('Latitude(Deg)');
ylabel('Longitude(Deg)');
zlabel('Altitude(m)');
legend('Solution', 'Truth');


