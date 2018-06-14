clc; 
global e a u R startIndex estimationInterval motionStartIndex includeRealConfiguration; 

e = 6.69437999013e-3; % squared eccentricity
a = 6378137; % big axle 
u = 7292115.8553e-11; % angular velocity
R = 6370000; % earth radius
startIndex = 1;   % initial index 
estimationInterval = 15270;  % interval average
motionStartIndex = estimationInterval + 1;   % first index pf data, when motion starts
includeRealConfiguration = true; % include optional plots and calculations

% deserialize input data (local measures) in start position t = 0
data = load('start.dat'); % initial conditions

phi = rad2deg(data(1)); % latitude
lambda = rad2deg(data(2)); % longtitude

h0 = data(3); % initial height
course = data(4); % course angle

% deserialize input data (accelerometers and angular velocity measures)
data = load('imu.dat'); 
t = data(:, 1);

omega1 = deg2rad(data(:, 2));  
omega2 = deg2rad(data(:, 3));  
omega3 = deg2rad(data(:, 4));  
    
z1 = data(:, 5); 
z2 = data(:, 6);
z3 = data(:, 7);

% get mean on estimation interval
fz = mean(data(1 : estimationInterval, 5 : 7));
disp(fz);
uz = mean(data(1 : estimationInterval, 2 : 4));  
disp(uz);

% deserialize correct values of coordinates and velocities
data = load('trj.dat');
t2 = data(:,1); % time

lat = deg2rad(data(:,2));  
long = deg2rad(data(:,3)); 
h = data(:,4); 

Ve = data(:,5); % velocity to East
Vn = data(:,6); % velocity to Nord
Vh = data(:,7); % velocity with respect to surface normal

% Krylov's angles
head = deg2rad(data(:,8));  
pitch = deg2rad(data(:,9));   
roll = deg2rad(data(:,10)); 

% Init values of velocities and coordinates
phi0 = lat(startIndex); 
lambda0 = long(startIndex);
Ve0 = Ve(startIndex); 
Vn0 = Vn(startIndex);
Vh0 = Vh(startIndex);

% calculating real free fall acceleration
g = 9.78030*(1 + 0.0053020*(sin(phi0))^2 - 0.0000070*(sin(2*phi0))^2) - 0.00014 - 0.000003086*h0;
 
% calculating orientation matrix L
L = zeros(3, 3);
% properTeta  = 0;
% properGamma = 0;
% properPsi = -1.483529864195180;
% L(:, 3) = fz / g;
L(:, 3) = [0 0 1];
sinTeta = L(2, 3);
properTeta = asin(sinTeta);
cosGamma = L(3, 3) / cos(properTeta);
properGamma = -abs(acos(cosGamma));
%re-calculating 3rd column of L by true values of angles
L( :, 3) = [
            -cos(properTeta) * sin(properGamma); 
            sin(properTeta); 
            cos(properTeta) * cos(properTeta)
            ];
        
L(:, 2) = [
            (uz(1)-L(1,3)*norm(uz)*sin(phi0))/(cos(phi0)*norm(uz));
            (uz(2)-L(2,3)*norm(uz)*sin(phi0))/(cos(phi0)*norm(uz));
            (uz(3)-L(3,3)*norm(uz)*sin(phi0))/(cos(phi0)*norm(uz))
            ];

cosPsi = L(2, 2) / cos(properTeta);
properPsi = -acos(cosPsi);

L(:, 2) = [
           sin(properPsi)*cos(properTeta) + cos(properPsi)*sin(properTeta)*sin(properGamma); 
           cos(properPsi)*cos(properTeta); 
           sin(properPsi)*sin(properGamma) - cos(properPsi)*sin(properTeta)*cos(properGamma)
           ];

L(:, 1) = [
           cos(properPsi)*cos(properGamma) - sin(properPsi)*sin(properTeta)*sin(properGamma); 
           -sin(properPsi)*cos(properTeta); 
           cos(properPsi)*sin(properGamma) + sin(properPsi)*sin(properTeta)*cos(properGamma)
           ]; 

% start 2nd stage of solution
% initializing matrices 
Az = L; 
Ax = eye(3);

% 
h_s = zeros(size(z1));
Ve_s = zeros(size(z1));
Vn_s = zeros(size(z1));
Vh_s = zeros(size(z1));
psi = zeros(size(z1));
gam = zeros(size(z1));
tet = zeros(size(z1));
for i = 1 : motionStartIndex
   lambda(i) = lambda0;
   phi(i) = phi0;
   h_s(i) = h0;
   Ve_s(i) = Ve0;
   Vn_s(i) = Vn0;
   Vh_s(i) = Vh0;
   psi(i) = properPsi;
   gam(i) = properGamma;
   tet(i) = properTeta;
end

if includeRealConfiguration == true
    deltaTDotV1 = zeros(size(Ve) - 3);
    deltaTDotV2 = zeros(size(Ve) - 3);
    deltaTDotV3 = zeros(size(Ve) - 3);
end

%dec size(t) to avoid exceptions in getting k+1 element
for k = motionStartIndex : (size(t)-1)
    % deltaT time
    deltaT = t(k+1) - t(k);

    omega = sqrt(omega1(k)^2 + omega2(k)^2 + omega3(k)^2);
    % coeffiecients from (3.12) for Az calculating
    c1 = sin(omega*deltaT) / omega;
    c2 = (1-cos(omega*deltaT)) / (omega^2) ;

    % (3.4) 
    Omega = [0 ,            omega3(k),      - omega2(k);...  
            -omega3(k),     0 ,               omega1(k);...
             omega2(k),      - omega1(k) ,           0];
    
    % A_i+1  recurrent method(3.12)
    % using previous Az calculated 1 iteration before
    Az = (eye(3) + c1*Omega+ c2* Omega^2) * Az;  

    % (3.7)
    % Curve radiuses (east and nord)
    Re = a / sqrt(1 - e*(sin(phi(k)))^2);
    Rn = a * (1 - e) / (1 - e * (sin(phi(k)))^2)^(3/2);
    
    % calculating current free fall acceleration
    g = 9.78030*(1 + 0.0053020 * (sin(phi(k)))^2 - 0.0000070*(sin(2*phi(k)))^2) - 0.00014 - 0.000003086*h_s(k);
    
    % angular velocity with respect to geographical axles
    uE = 0;  
    uN = u * cos(phi(k));
    uH = u * sin(phi(k));

    % (3.6)
    % components of geographical repper angular velocity with respect to Earth
    omegaE = -Vn_s(k) / (Rn + h_s(k));  
    omegaN = Ve_s(k) / (Re + h_s(k));   
    omegaH = omegaN * tan(phi(k));  

    % Filling matrix 
    Omega = [0 ,            omegaH+uH,      - omegaN-uN;...
            -omegaH-uH,     0 ,               omegaE+uE;...
             omegaN+uN,     -omegaE-uE ,             0];

    % calculating omega + ux (defined above)
    omega = norm([omegaE+uE, omegaN+uN, omegaH + uH]);
 
    c1 = sin(omega*deltaT) / omega;
    c2 = (1 - cos(omega*deltaT)) / (omega^2);
    F = eye(3) + c1*Omega + c2* Omega^2;
    Ax = F * Ax; % A_i+1

    % Re-calculating L in new iteration
    % Has to be always ortogonal
    L = Az * Ax'; 

    % (3.7)
    fz = [z1(k); z2(k); z3(k)];
    fy = L' * fz;  

    % (3.14 (2)) equations
    dotVe =  (omegaH + 2*uH)*Vn_s(k) - (omegaN + 2*uN)*Vh_s(k) + fy(1,1);
    dotVn = -(omegaH + 2*uH)*Ve_s(k) + (omegaE + 2*uE)*Vh_s(k) + fy(2,1);
    dotVh =  (omegaN + 2*uN)*Ve_s(k) - (omegaE + 2*uE)*Vn_s(k) + fy(3,1) - g;

    for indexVe = 1 : size(Ve)
        % Optional calculations for plots
        % k - index of outer loop
         if includeRealConfiguration == true
            if t(k) == t2(indexVe)
                teta2 = pitch(indexVe);
                gamma2 = roll(indexVe);
                psi2 = head(indexVe);
                L1( :,3) = [-cos(teta2) * sin(gamma2); sin(teta2); cos(teta2)*cos(gamma2)];
                L1(:, 2) = [sin(psi2)*cos(gamma2) + cos(psi2)*sin(teta2)*sin(gamma2); cos(psi2)*cos(teta2); sin(psi2)*sin(gamma2) - cos(psi2)*sin(teta2)*cos(gamma2)];
                L1(:, 1) = [cos(psi2)*cos(gamma2) - sin(psi2)*sin(teta2)*sin(gamma2); -sin(psi2)*cos(teta2); cos(psi2)*sin(gamma2) + sin(psi2)*sin(teta2)*cos(gamma2)]; 
             
                fy1 = L1 * fz;
                
                dotVe1 = (omegaH + 2*uH)*Vn(indexVe) - (omegaN + 2*uN)*Vh(indexVe) + fy1(1,1);
                dotVn1 = -(omegaH + 2*uH)*Ve(indexVe) + (omegaE + 2*uE)*Vh(indexVe) + fy1(2,1);
                dotVh1 = (omegaN + 2*uN)*Ve(indexVe) - (omegaE + 2*uE)*Vn(indexVe) + fy1(3,1) - g;

                deltaTDotV1(indexVe) = (dotVe - dotVe1);
                deltaTDotV2(indexVe) = (dotVn - dotVn1);
                deltaTDotV3(indexVe) = (dotVh - dotVh1);
            end
         end
    end

    % filling last delta's of velocities
    Ve_s(k+1)=  Ve_s(k) + deltaT*dotVe;  
    Vn_s(k+1)=  Vn_s(k) + deltaT*dotVn;  
    Vh_s(k+1)=  Vh_s(k) + deltaT*dotVh;
    
    % calculating new values of geographical coordinates
    lambda(k+1) = lambda(k) + deltaT*Ve_s(k)/((Re+h_s(k))*cos(phi(k)));
    phi(k+1)= phi(k) + deltaT*Vn_s(k)/(Rn+h_s(k));
    h_s(k+1) = h_s(k) + deltaT*Vh_s(k);

end

% transform to Cartesian coordinates and move to (0,0)
% h - already good
x_0 = (R + h(1)) * long(1);
y_0 = (R + h(1)) * log(tan(pi/4 + 0.5*lat(1)));
x = zeros(size(long));
y = zeros(size(long));
for i1 = 1 : size(long) 
    x(i1) = (R + h(i1)) * long(i1) - x_0;
    y(i1) = (R + h(i1)) * log(tan(pi/4 + 0.5*lat(i1))) - y_0;
end

% transform to Cartesian coordinates and mobe to (0,0)
x_1_0 = (R + h_s(1)) * lambda(1);
y_1_0 = (R + h_s(1)) * log(tan(pi/4 + 0.5*phi(1)));
x1 = zeros(size(lambda'));
y1 = zeros(size(lambda'));
h1 = zeros(size(lambda'));
for i1 = 1 : size(lambda')
        x1(i1) = (R + h_s(i1)) * lambda(i1) - x_1_0;
        y1(i1) = (R + h_s(i1)) * log(tan(pi/4 + 0.5*phi(i1))) - y_1_0;
        h1(i1) = h_s(i1);
end    


% Printig trajectories - true and calculated
figure('Name', 'Trajectories');
clf;
plot3(x,y,h,'k','LineWidth' ,3);
hold on;
plot3(x1,y1,h1,'g','LineWidth' ,3);
grid on;

% Calculating error beetween positions in the end of motion
r1 = [  (R + h_s(k)) * sin(phi(k)) * cos(lambda(k))
        (R + h_s(k)) * sin(phi(k)) * sin(lambda(k))
        (R + h_s(k)) * cos(phi(k))];
    
r2 = [  (R + h(901)) * sin(lat(901)) * cos(long(901))
        (R + h(901)) * sin(lat(901)) * sin(long(901))
        (R + h(901)) * cos(lat(901))];
    
% Distance beetween endpoints
r = norm(r1 - r2);
disp(r);

if includeRealConfiguration == true
    deltaTV1 = zeros(901, 1);
    deltaTV2 = zeros(901, 1);
    deltaTV3 = zeros(901, 1);
    
    % calculating velocity differences in last time (optional)
    dotVe1 =  (omegaH + 2*uH)*Vn(end) - (omegaN + 2*uN)*Vh(end) + fy(1,1);
    dotVn1 = -(omegaH + 2*uH)*Ve(end) + (omegaE + 2*uE)*Vh(end) + fy(2,1);
    dotVh1 =  (omegaN + 2*uN)*Ve(end) - (omegaE + 2*uE)*Vn(end) + fy(3,1) - g;
    deltaTDotV1(end) = (dotVe - dotVe1);
    deltaTDotV2(end) = (dotVn - dotVn1);
    deltaTDotV3(end) = (dotVh - dotVh1);

    % Differnces beetween velocities
    for i = 1 : size(Ve)
        for indexVe = 1 : size(Ve_s)
             if t(indexVe) == t2(i)
                 deltaTV1(i) = (Ve(i) - Ve_s(indexVe));
                 deltaTV2(i) = (Vh(i) - Vh_s(indexVe));
                 deltaTV3(i) = (Vn(i) - Vn_s(indexVe));
             end
        end
    end

    deltaTV1(end) = (Ve(end) - Ve_s(end));
    deltaTV2(end) = (Vh(end) - Vh_s(end));
    deltaTV3(end) = (Vn(end) - Vn_s(end));

    % Printing plots of dotV differnces
    figure('Name', 'Delta dotV east');clf;
    plot(t2(3 : 901), deltaTDotV1);
    hold on;
    grid on;

    figure('Name', 'Delta dotV height');clf;
    plot(t2(3 : 901), deltaTDotV2);
    hold on;
    grid on;

    figure('Name', 'Delta dotV north');clf;
    plot(t2(3 : 901), deltaTDotV3);
    hold on;
    grid on;

    % Printing plots of V differences
    figure('Name', 'Delta V east');clf;
    plot(t2, deltaTV1);
    hold on;
    grid on;

    figure('Name', 'Delta V height');clf;
    plot(t2, deltaTV2);
    hold on;
    grid on;

    figure('Name', 'Delta V north');clf;
    plot(t2, deltaTV3);
    hold on;
    grid on;

    % Printing velocities 
    figure('Name', 'East Velocity (true and calculated)');clf;
    plot(t2, Ve)
    hold on
    plot(t,Ve_s)
    hold on

    figure('Name', 'Height Velocity (true and calculated)');clf;
    plot(t2, Vh)
    hold on
    plot(t, Vh_s)
    hold on

    figure('Name', 'North Velocity (true and calculated)');clf;
    plot(t2, Vn)
    hold on
    plot(t, Vn_s)
    hold on
    grid on;

end
  
   