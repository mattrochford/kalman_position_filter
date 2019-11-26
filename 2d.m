%% EE 525 Final Project
% Kalman Filter for 2D Position Tracking

clc
clear all 
close all

%% Import Dataset
% Important position data
fileID = fopen('2Ddata_r.txt','r');
formatSpec = '%*s %d %d';
sizeA = [2 Inf];
A_pos = fscanf(fileID,formatSpec,sizeA); % Position values in meters

% Import file as single string
fileID = fopen('2Ddata_r.txt','r');
formatSpec = '%c';
sizeA = [1 Inf];
A_time = fscanf(fileID,formatSpec,sizeA);

% Initialize values
colon = 0;
i = 1;
j = 1;

% Extract time data from string
while i < length(A_time)
    if (A_time(i) == ':')
        colon = colon + 1; % Increment colon flag
    end
    if colon == 2 
        colon = 0; % Reset colon flag
        time(j) = str2num(A_time(i+1:i+6)); % Store relevant data
        j = j + 1; 
    end
    i = i + 1;
end

% Convert time data to dt values
for i = 1:length(time)-1
    dt(i) = time(i+1)-time(i);
end
% Correct negative values
dt(dt<0) = dt(dt<0)+60;

% Calculate velocities
%for i = 1:length(dt)
%    A_vel(1,i) = (A_pos(1,i+1)-A_pos(1,i))/dt(i); % Y velocity
%    A_vel(2,i) = (A_pos(2,i+1)-A_pos(2,i))/dt(i); % X velocity
%end

% Assign final measurement values
len = length(A_pos);
x0 = A_pos(2,1);
y0 = A_pos(1,1);
y_pos = A_pos(1,2:len);
x_pos = A_pos(2,2:len);
%y_vel = A_vel(1,:);
%x_vel = A_vel(2,:);

% Clear unneeded variables
clearvars -except y_pos x_pos x0 y0 dt

%% Design state model
W1 = 1; % PSD for X acceleration
W2 = 1; % PSD for Y acceleration

F = zeros(4);
F(1,2) = 1;
F(3,4) = 1;

G = zeros(4,2);
G(2,1) = sqrt(W1);
G(4,2) = sqrt(W2);

H = [1 0 0 0; 0 0 1 0];

%% Initialize Phi and Q
sigma = var(dt);
dt = mean(dt); % Sampling interval
n = 4; % Dimension of x(t)

% Compute A matrix
A = [-F*dt G*W1*G'*dt; zeros(n) F'*dt];

% Compute B matrix
B = expm(A);

% Compute Phi and Q
phi = B(n+1:2*n,n+1:2*n)';
Q = phi*B(1:n,n+1:2*n);

% Initialize Rk (measurement variance)
R = 10000*eye(n/2);

%% Implement Kalman Filter

P_ = 0.1*eye(n); % Initialize error
X_ = [x0; 0; y0; 0]; % Initial guess
Z = [x_pos; y_pos];

for k = 1:length(x_pos)
    % Measurement Update
    K = P_*H'*inv(H*P_*H'+R); % Compute Kalman Gain
    X(:,k) = X_+K*(Z(:,k)-H*X_); % Update estimate
    P = (eye(n)-K*H)*P_;
    
    % Time Update
    X_ = phi*X(:,k); % Project state ahead
    P_ = phi*P*phi'+Q; % Project error ahead
     
end

figure
plot(x_pos,y_pos,X(1,:),X(3,:))
legend('Measured Position','Filtered Position')
title('Kalman Filtered 2D Position')
xlabel('X Position (meters)')
ylabel('Y Position (meters)')
grid on
