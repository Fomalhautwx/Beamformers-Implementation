% Generalized Sidelobe Canceller Implementation
% Griffith & Jim
% 2024.02.21

close all; 
clear all; 

M = 4;  % num of mics

Wc = ones(M, 1) / M;                % Fixed Beamformer
Ws = [1,1,-1,-1;1,-1,-1,1;1,-1,1,1]; % Blocking Matrix
Al = randn(3, 1);                   % Random initialization of MC

SIG_LENGTH = 10000;
tau = 1;    % define frame gap
% F_curly = [1, -2, 1.5, 2]';  % define the temporal response
F_curly = [1, 0, 0, 0]';
% Pseudo-Gaussian signal 
s1 = 0.5 * randn(SIG_LENGTH, 1); % Look direction signal
bpFilt = designfilt('bandpassfir', 'FilterOrder', 100, ...
             'CutoffFrequency1', 0.22, 'CutoffFrequency2', 0.28,...
             'SampleRate', 1);
s1 = filter(bpFilt,s1);
s2 = 0.5 * randn(SIG_LENGTH, 1); % Interferer #1
bpFilt = designfilt('bandpassfir', 'FilterOrder', 100, ...
             'CutoffFrequency1', 0.1, 'CutoffFrequency2', 0.12,...
             'SampleRate', 1);
s2 = filter(bpFilt,s2);
s3 = 0.5 * randn(SIG_LENGTH, 1); % Interferer #2
bpFilt = designfilt('bandpassfir', 'FilterOrder', 100, ...
             'CutoffFrequency1', 0.35, 'CutoffFrequency2', 0.4,...
             'SampleRate', 1);
s3 = filter(bpFilt,s3);

% Steering: t_k = klsin(theta)/c
l = 0.1; c = 340; sample_rate = 13600; noise_weight = 0.01;
theta0 = 0*pi/180; theta1 = 60*pi/180; theta2 = -30*pi/180; % define the impinging directions
Delay0 = round([0:3]' * (l * sin(theta0) * sample_rate / c)); 
Delay1 = round([0:3]' * (l * sin(theta1) * sample_rate / c)); 
Delay2 = round([0:3]' * (l * sin(theta2) * sample_rate / c)); 
% Delay and comb
x1 = s1 + s2 + s3 + noise_weight * randn(length(s1), 1); 
x2 = circshift(s1, Delay0(2))+circshift(s2, Delay1(2))+circshift(s3, Delay2(2))+noise_weight * randn(length(s1), 1); 
x3 = circshift(s1, Delay0(3))+circshift(s2, Delay1(3))+circshift(s3, Delay2(3))+noise_weight * randn(length(s1), 1);
x4 = circshift(s1, Delay0(4))+circshift(s2, Delay1(4))+circshift(s3, Delay2(4))+noise_weight * randn(length(s1), 1);
Received_signal = [x1, x2, x3, x4];
clear x1 x2 x3 x4;

mu = 0.01;  % Learning rate
Power = zeros(SIG_LENGTH, 1);
Target_power = s1'*s1 / SIG_LENGTH;
Yo = zeros(SIG_LENGTH, 1);

for index=1:SIG_LENGTH
    input = Received_signal(index, :);
    % FBF
    yc = input*Wc;
    % BM
    X_prime = Ws*input';
    % MC
    ya = Al'*X_prime;
    yo = yc - ya;   % Yield the output
    Al = Al+mu*yo*X_prime / (X_prime'*X_prime); % NLMS
    
    Power(index) = yo'*yo / Target_power;   % Output power normalized by target signal power.
    Yo(index) = yo;
end

figure;
plot(Power, '.')

