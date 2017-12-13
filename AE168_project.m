%% AE 168 Project
 % Team HALO minus Alexis
 % 12/13/2017

clear all
close all
clc
%% Defining Constants

% Geometric data
W = 607.53;      % Empty weight lbs
S = 260.94;      % Wing surface area ft^2
b = 75.77;       % Wing span ft
sweep_LE = 10;   % Sweep angle at the LE deg
taper = 0.5;     % Taper ratio
I_x = 74.02;     % slugs-ft^2
I_y = 2844.27;   % slugs-ft^2
I_z = 2893.58;   % slugs-ft^2
m = W/32.2;      % Mass slugs
x_cm = -.3880;   % LE to CM in ft
x_Wac = 0.8929;  % LE to .25c wing ft

% Basic Aircraft Data
rho = 0.0001705; % Density slugs/ft^3
U_1 = 170;       % True airspeed ft/sec
e = 0.80;        % Oswald's efficiency 
C_D0 = 0.003392; % Zero lift drag
c_la = 5.73;     % 2D wing lift curve slope rad^-1
eta_h = 1;       % Dynamic pressure ratio - No tail
tau = 0.125;     % Elevator effectiveness 
g = 32.2;        % Standard gravity ft/sec^2
M = 0.176;       % Mach number
alpha0 = -13.2;  % Zero lift AOA deg 

%% Calculating constants

% Basic wing dimmensions
c_root = 2*S/(b*(1+taper));     % Chord at the root in ft
AR_W = (b^2)/S;                 % Wing aspect ratio
c_barW = (2/3)*c_root*((1+taper+taper^2)/(1+taper));  % Wing mean geometric chord ft

% Dynamic Pressure
q_bar = 0.5*rho*(U_1)^2;

% Wing lift coefficient variation with AOA
beta = (1-M^2)^.5;
kappa = c_la/(2*pi/beta);

%% Calculating Longitudinal Stabilty & Control Coefficients

tan_sweep_50 = tand(sweep_LE)-((4/AR_W)*((1-taper)/(2*(1+taper))));
C_LAW = 2*pi*AR_W/(2+(((AR_W^2*beta^2)/kappa^2)*(1+(((tan_sweep_50)^2)/beta^2))+4)^.5);

% Pitching Moment Coefficient variation AOA
C_Ma = C_LAW*((x_cm/c_barW)-(x_Wac/c_barW));

% Pitch Damping due to elevator deflection
C_Mde = -(S/S)*eta_h*((x_Wac/c_barW)-(x_cm/c_barW))*C_LAW*tau;

% Lift Coefficient of entire A/C
C_L1 = W/(q_bar*S);

% Drag Coefficient of entire A/C 
C_D1 = C_D0 + ((C_L1^2)/(pi*AR_W*e));

% Drag Coefficient variation due to forward velocity perturbation
C_DU = 0;

% Lift coefficient variation with AOA
C_LA = C_LAW; % Since there's no tail

% Drag coeffcient variation with AOA
C_DA = ((2*C_L1*C_LA)/(pi*AR_W*e));

% Lift coefficient variation with forward velocity perturbation
C_L0 = C_LA*(-alpha0)*pi/180;
C_LU = M^2*C_L0/(1-M^2); 

% Lift coefficient variation with elevator deflection
C_Lde = C_LAW*eta_h*(S/S)*tau;

% Pitch Damping due to pitch rate
l_W = x_Wac-x_cm;
tan_sweep_25 = tand(sweep_LE)-((4/AR_W)*((1-taper)/(4*(1+taper))));
sweep_25 = atand(tan_sweep_25);
C_Mq = -1*c_la*cosd(sweep_25)*((AR_W*(.5*l_W/c_barW+2*(l_W/c_barW)^2) / ...
    (AR_W+2*cosd(sweep_25)))+(1/24)*(AR_W^3*(tand(sweep_25))^2/(AR_W+6*cosd(sweep_25)))+1/8);

% Pitching moment coefficient due to AOA rate
C_Madot1 = -81/32*(x_Wac/c_root)^2*C_LAW;
C_Ladot = 2*C_LAW*(l_W/c_barW)*eta_h*(S/S);
C_Madot = C_Madot1+(x_cm/c_barW)*C_Ladot;
%% Calculating Longitudinal Stability & Control Forces and Moments

% Pitching Moment resulting from AOA
M_a = (q_bar*S*c_barW*C_Ma)/I_y;

% Pitching Moment resulting from change in elevator deflection
M_de = (q_bar*S*c_barW*C_Mde)/I_y;

% Change in X-force resulting from change in forward velocity
X_U = (1/(m*U_1))*((-C_DU-2*C_D1)*q_bar*S);

% Change in X-force resulting from change in AOA
X_alpha = (1/m)*(-C_DA+C_L1)*q_bar*S;

% Change in Z-force resulting from change in forward velocity
Z_U = (1/(m*U_1)*(-C_LU-(2*C_L1))*q_bar*S);

% Change in Z-force resulting from change in AOA
Z_alpha = (1/m)*(-C_LA-C_D1)*q_bar*S;

% Change in Z-force resulting from elevator deflection
Z_de = -(1/m)*C_Lde*q_bar*S;

% Change in Pitching Moment resulting from pitch rate
M_q = (q_bar*S*(c_barW^2)*C_Mq)/(2*I_y*U_1);

% Pitching moment due to AOA rate
M_adot = (q_bar*S*(c_barW^2)*C_Madot)/(2*I_y*U_1);


%% Longitudinal Transfer Functions 

s = tf('s');

% Phugoid Approximation
airspeed = ((g*(Z_de))/U_1)/(s^2 + (-X_U)*s + (-g*Z_U)/U_1); % Airspeed TF
pitch_att = ((Z_de/U_1)*s + (-Z_de*X_U)/U_1) / ...
    (s^2 + (-X_U)*s + (-g*Z_U)/U_1); % Theta TF

%Short Period Approximation
pitch_rate = (( (U_1*M_de + Z_de*M_adot)*s + (M_a*Z_de - Z_alpha*M_de)) / ...
    (U_1*(s^2 - (M_q+M_adot+(Z_alpha/U_1))*s + (M_q*(Z_alpha/U_1)-M_a))));

AOA = ((Z_de*s + (U_1*M_de-(M_q*Z_de))) / (U_1*(s^2 - ...
    (M_q+M_adot + (Z_alpha/U_1))*s + (M_q*(Z_alpha/U_1)-M_a))));

%% Lateral Stability & Control

% elevon position
y1 = 23.414;
y2 = 32.885;

%                           Lateral Coefficients
% Roll Moment Coefficient due to Aileron Deflection
C_roll_delA = (2*C_LAW*tau/(S*b))*c_root*((y2^2-y1^2)/2 - (2/3)*((1-taper)/b)*(y2^3 - y1^3));
% Roll Damping Moment Coefficient
C_roll_p = -(4*C_LAW/(S*b^2))*(1/8)*c_root*(b^3)*((1/3) - (1/4)*(1-taper));

%                           Lateral Moments
% Roll Moment due to Aileron Deflection
L_delA = q_bar*S*b*C_roll_delA/I_x;
% Roll Moment due to Roll Rate
L_p = q_bar*S*b^2*C_roll_p/(2*I_y*U_1);

% Transfer Function for Roll Rate due to Aileron Deflection
roll_rate = L_delA/(s - L_p);

% Transfer Function for Bank Angle due to Aileron Defleciton
bank_angle = roll_rate/s;

%% Servo Actuator Transfer Function





%% Simulink

SimTime = 200; % seconds
stepSize = 0.1; % step size
t = 0:stepSize:SimTime; % time vector

[pitch_num,pitch_den] = tfdata(pitch_att,'v');

sim('pitch_attitude')

figure
plot(pitch_att_sim(:,1),pitch_att_sim(:,2))
legend('show')
xlabel('Time (seconds)')
ylabel('Amplitude')
title('Step Response')

%%


