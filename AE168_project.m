%% AE 168 Project
 % Team HALO minus Alexis

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
rho = 0.0001705;  % Density slugs/ft^3
U_1 = 170;        % True airspeed ft/sec
e = 0.80;         % Oswald's efficiency 
C_D0 = 0.003392;  % Zero lift drag
c_la = 5.73;      % 2D wing lift curve slope rad^-1
eta_h = 1;        % Dynamic pressure ratio - No tail
tau_elev = 0.125; % Elevator effectiveness 
g = 32.2;         % Standard gravity ft/sec^2
M = 0.176;        % Mach number
alpha0 = -13.2;   % Zero lift AOA deg 

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
C_Mde = -(S/S)*eta_h*((x_Wac/c_barW)-(x_cm/c_barW))*C_LAW*tau_elev;

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
C_Lde = C_LAW*eta_h*(S/S)*tau_elev;

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
C_roll_delA = (2*C_LAW*tau_elev/(S*b))*c_root*((y2^2-y1^2)/2 - (2/3)*((1-taper)/b)*(y2^3 - y1^3));
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

%% System Variables

SimTime = 15; % seconds
stepSize = 0.001; % step size
t = 0:stepSize:SimTime; % time vector
theta_cmd = 10;

%% Servo Actuator Transfer Function

tau = 0.1; % time constant
servo = 1/(tau*s + 1); % elevon servo gain

%% Pitch Attitude ATTEMPT 1

%% Pitch Attitude Open Loop Root Locus

% Using Phugoid Approximation
figure
% step(servo*pitch_att,t)
rlocus(-servo*pitch_att)

% Using Short Period Approximation
figure
pitch_att_OL = servo*(-pitch_rate)*(1/s);
rlocus(pitch_att_OL)

% Pitch Attitude Root Locus Analysis
 % Analyze Short Period Pitch Attitude Root Locus HERE

%% Pitch Attitude Closed Loop Feedback (w/o Closed Loop Pitch Damping)
K_theta = 17;
pitch_att_CL = feedback(servo*(-pitch_rate)*(1/s)*K_theta,1);

figure
step(pitch_att_CL,t)

% Pitch Attitude Closed Loop Analysis (w/o pitch rate)
 % Analyze Short Period Pitch Attitude Closed Loop HERE
 
%% Pitch Attitude ATTEMPT 2a

%% Pitch Rate Open Loop Root Locus

figure
rlocus(-pitch_rate*servo)

% Pitch Rate Root Locus Analysis
 % Analyze Short Period Pitch Rate Root Locus HERE

 % K_thetadot gains can range from 6.47 to 11.2
%% Pitch Damping Closed Loop Feedback

K_thetadot = 6.5;
pitch_rate_CL = feedback(servo*(-pitch_rate),K_thetadot);

figure
rlocus(pitch_rate_CL*(1/s))

% Pitch Damping Closed Loop Analysis
 % Write Analysis

 % K_thetadot = 6.47, root locus crosses imaginary axis at K_theta =~123
 % K_thetadot = 11.2, root locus crosses imaginary axis at K_theta =~ 149 
    % also no damping from 0 to 8.14
%% Pitch Attitude Closed Loop (w/ Pitch Damping) in Simulink
K_theta_damping = 20;
[q_num,q_den] = tfdata(-pitch_rate,'v');
sim('pitch_attitude_SP')

figure
plot(pitch_att_damp_sim(:,1),pitch_att_damp_sim(:,2))
legend('Simulink Pitch Attitude')
xlabel('Time (seconds)')
ylabel('Pitch Attitude')
title('Pitch Attitude Closed Loop w/ Pitch Damping Step Response')

%% Pitch Attitude Attempt 2b

%%
clc
close all
figure
plot([0 SimTime],[theta_cmd theta_cmd],'r')
legend_text = {'ref'};
%% PID Tuning

% Ziegler-Nichols Method Attempt
% Ku = 130;
% Tu = 1.25;
% Kp = 0.7*Ku;
% Ki = Tu/2.5;
% Kd = 3*Tu/20;

% Manual Tuning
Kp = 100;
Ki = 25;
Kd = 20;

%% Pitch Attitude with PID Tuner

sim('pitch_attitude_SP_with_PID')
hold on
plot(pitch_PID_sim(:,1),pitch_PID_sim(:,2))
legend_text = [legend_text {['Kp ' num2str(Kp) '  Ki ' num2str(Ki) '  Kd ' num2str(Kd)]}];
legend(legend_text)

%% ROLL

%% Roll Rate Open Loop Root Locus

figure
rlocus(roll_rate*servo)

% Roll Rate Root Locus Analysis
 % Analyze Short Period Pitch Rate Root Locus HERE

 % K_thetadot gains can range from 0 to 0.0393
%% Roll Damping Closed Loop Feedback

K_phidot = 0.0393;
roll_rate_CL = feedback(servo*(roll_rate),K_phidot);

figure
rlocus(roll_rate_CL*(1/s))

% Roll Damping Closed Loop Analysis
 % Write Analysis

 % K_phidot = 0.0393, root locus crosses imaginary axis at K_phi =~1.45
 % No damping from 0 to 0.116

%% Bank Angle with PID Tuner

K_phi = 0.116;
phi_disturb = 10;
[r_num,r_den] = tfdata(roll_rate,'v');
sim('bank_angle_w_damping')

figure
plot(bank_angle_sim(:,1),bank_angle_sim(:,2))
legend('Bank Angle from Simulink')
xlabel('Time (seconds)')
ylabel('Bank Angle')
title('Bank Angle Closed Loop w/ Roll Damping Step Response')
% axis([0 SimTime -2 phi_disturb*2])








