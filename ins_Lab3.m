clear all; close all; clc;

% load sensor data
load ('att_truth'); load ('sensors_clean'); load ('sensors_noisy'); 
load('vel_e');load pos_llh.mat; load pos_ecef.mat;

% store data for att_truth
% T = length(att_truth);
% datav = att_truth(2:end,:);
% roll = datav(1,:); pitch = datav(2,:); yaw = datav(3,:);
% dataTruth = [roll; pitch; yaw];

% store data for sensors_clean
dataClean = sensors_clean(2:end,:);
ax_clean = dataClean(2,:); ay_clean = dataClean(3,:); az_clean = dataClean(4,:);
wx_clean = dataClean(4,:); wy_clean = dataClean(5,:); wz_clean = dataClean(6,:);
fb = [ax_clean; ay_clean; az_clean]; 
wibb = [wx_clean; wy_clean; wz_clean];

% store data for sensors_noisy
dataNoisy = sensors_noisy(2:end,:);
ax_noisy = dataNoisy(2,:); ay_noisy = dataNoisy(3,:); az_noisy = dataNoisy(4,:);
wx_noisy = dataNoisy(4,:); wy_noisy = dataNoisy(5,:); wz_noisy = dataNoisy(6,:);
fbNoisy = [ax_noisy; ay_noisy; az_noisy]; 
wibbNoisy = [wx_noisy; wy_noisy; wz_noisy];

%% initializing variables
omega = 7.292115e-5; r0 = 6378137;
latitude = zeros(1,12001); longditude = zeros(1,12001); height = zeros(1,12001);
latitude(1) = -0.4784; longditude(1) = 2.6724; height(1) = 300;
VN(1) = 30.6; VE(1) = -12.7; VD(1) = -10.54; 
vn(:,1) = [30.6; -12.7; -10.54];%wn_b
g = [0; 0; 9.79];

% convert euler to quarternion
phi(1) = 0; theta(1) = 0; psi(1) = 5.8407;

a = (cos(psi/2)*cos(theta/2)*cos(phi/2)) + (sin(psi/2)*sin(theta/2)*sin(phi/2));
b = (cos(psi/2)*cos(theta/2)*sin(phi/2)) - (sin(psi/2)*sin(theta/2)*cos(phi/2));
c = (cos(psi/2)*sin(theta/2)*cos(phi/2)) + (sin(psi/2)*cos(theta/2)*sin(phi/2));
d = (sin(psi/2)*cos(theta/2)*cos(phi/2)) - (cos(psi/2)*sin(theta/2)*sin(phi/2));

% store quaternion value in column vector
q1 = [a; b; c; d];%quaternion1
% qk_clean = zeros(4,T);
% qk_clean(:,1) = q1;

dcm1 = abcd2dcm(a,b,c,d);
q2(:,1) = q1;



%% for clean sensors
eulerClean = zeros(3,12001);
for k = 2:12001
    % new equations updated each loop
    gln = g - (((omega^2*(r0+height(k-1)))/2).*[(sin(2*latitude(k-1))); 0; (1+cos(2*latitude(k-1)))]);
    wien = [omega * cos(latitude(k-1)); 0; -omega*sin(latitude(k-1))];
    wenn = [vn(2,k-1)/(r0+height(k-1)); -vn(2,k-1)/(r0+height(k-1)); -vn(2,k-1)*tan(latitude(k-1))/(r0+height(k-1))];
    
    % create DCM from new quaternion vector
    Cbn = abcd2dcm(a, b, c, d);
    Cnb = Cbn';
    
    % Compute body rates w.r.t the navigation frame
    wnbb(:,k) = wibb(:,k) - (Cnb * ((wien + wenn)));
%     wnbb = wnbb';
 
    % Propagate attitudes in quaternion forms
    Pnbb= [0, wnbb(:,k)']';
    
    Qq = [a, -b, -c, -d;
          b,  a, -d,  c;
          c,  d,  a, -b;
          d, -c,  b,  a]; 
        
            
    % Update quaternion using Euler Rule
    q_dot(:,k-1) = (0.5).* Qq * Pnbb;
    q2(:,k) = [a; b; c; d] + (0.01 * q_dot(:,k-1));
    q2(:,k)= q2(:,k)/norm(q2(:,k));
    
    % store quaternion value in column vector
    a = q2(1,k);
    b = q2(2,k);
    c = q2(3,k);
    d = q2(4,k);
    
    % Convert quaternion to Euler     
    psi(k) = atan2(2*(c*d+a*b),(a^2-b^2-c^2+d^2));
    theta(k) = asin(-2*(b*d-a*c));
    phi(k) = atan2(2*(b*c+a*d),(a^2+b^2-c^2-d^2));
    
    % Store converted euler
    eulerClean(:,k) = [psi(k);theta(k);phi(k)];
    
    dcm2 = abcd2dcm(a,b,c,d);
    
    %convert body frame to navigation frame
    fn(:,k-1) = dcm2 * fb(:,k-1);

    % true acceleration
    vdot(:,k-1) = fn(:,k-1) - cross((2*wien+wenn),vn(:,k-1))+gln;

    %euler rule
    vn(:,k) = vn(:,k-1) + 0.01*vdot(:,k-1);
    
    
    %euler rule second time to update the lat & lon
    latdot = vn(1,k)/(r0+height(k-1));
    londot = vn(2,k)*sec(latitude(k-1))/(r0+height(k-1));
    hdot = -vn(3,k);

    %euler rule second update on latitude and longditude
    latitude(k) = latitude(k-1) + 0.01*latdot; 
    longditude(k) = longditude(k-1) + 0.01*londot; 
    height(k) = height(k-1) + 0.01*hdot; 
    
end

VN(1) = vn(1,1);
VE(1) = vn(2,1);
VD(1) = vn(3,1);
latplot(1) = latitude(1);
lonplot(1) = longditude(1);
hplot(1) = height(1);

for i = 1:120;
    latplot(i+1) = latitude(i*100+1);
    lonplot(i+1) = longditude(i*100+1);
    hplot (i+1) = height(i*100+1);
    VN(i+1) = vn(1,i*100+1);
    VE(i+1) = vn(2,i*100+1);
    VD(i+1) = vn(3,i*100+1);  
end 

%%
figure
subplot(3,1,1)
plot(phi); hold on 
plot(att_truth(4,:));
title ('yaw attitude - true vs clean')
legend ('clean ','truth ')


subplot(3,1,2)
plot(theta); hold on
plot(att_truth(3,:));
title ('pitch attitude  - true vs clean')
legend ('clean ','truth ')


subplot(3,1,3)
plot(psi); hold on
plot(att_truth(2,:));
title ('roll attitude  - true vs clean')
legend ('clean ','truth ')

figure
plot (VN)
hold on
plot (vel_e(2,:))
title ('Velocity North')
legend ('Clean Data', 'True Data')
ylabel('m/s')
xlabel('step')

figure
plot (VE)
hold on
plot (vel_e(3,:))
title ('Velocity East')
legend ('Clean Data', 'True Data')
ylabel('m/s')
xlabel('step')

figure
plot (VD)
hold on
plot (vel_e(4,:))
title ('Velocity Down')
legend ('Clean Data', 'True Data')
ylabel('m/s')
xlabel('step')

figure
plot (latplot*180/pi)
hold on
plot (pos_llh(2,:)*180/pi)
title ('Latitude')
legend ('Clean Data', 'True Data')
ylabel('degrees')
xlabel('step')

figure
plot (lonplot*180/pi)
hold on
plot (pos_llh(3,:)*180/pi)
title ('Longditude')
legend ('Clean Data', 'True Data')
ylabel('degrees')
xlabel('step')

figure
plot (hplot*180/pi)
hold on
plot (pos_llh(4,:)*180/pi)
title ('Height')
legend ('Clean Data', 'True Data')
ylabel('degrees')
xlabel('step')



figure
plot (latplot*180/pi,lonplot*180/pi)
hold on
plot (pos_llh(2,:)*180/pi,pos_llh(3,:)*180/pi)
title ('Latitude vs Longditude')
legend ('Clean Data', 'True Data')
ylabel('degrees')
xlabel('degrees')

figure
plot3(latplot*180/pi,lonplot*180/pi,hplot)
hold on
plot3 (pos_llh(2,:)*180/pi,pos_llh(3,:)*180/pi,pos_llh(4,:))
title ('Latitude vs Longditude vs  Height')
legend ('Clean Data', 'True Data')
ylabel('degrees')
xlabel('degrees')
zlabel('meters')
