% FROST linearly constrained adaptive filter LMS
% Implementation 
% 2024.2.14

close all; 
clear all; 

K = 4;  % num of mics
J = 4;  % num of taps. Symbols corresponding the original paper.

SIG_LENGTH = 20000;
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

C = zeros(16, 4);
c = ones(4, 1);
C(1:4, 1) = c; C(5:8, 2) = c; C(9:12, 3) = c; C(13:16, 4) = c;
% define F and P
F = C * inv(C'*C) * F_curly;
I_16 = eye(K*J);
P = I_16 - C * inv(C'*C) * C';

mu = 0.02;% Learning rate
W = F;  % Initialization
% W = randn(16, 1);
Y = zeros(SIG_LENGTH, 1);
k = 1;  % Output index
Power = zeros(SIG_LENGTH, 1);
Mismatch = zeros(SIG_LENGTH, 1);

lambda = 0.9985;    % coefficient of forgetting
X0 = [Received_signal(4, :), Received_signal(3, :), Received_signal(2, :), Received_signal(1, :)]';
Rxx = X0*X0';

%index = 1;  % Iteration index

for index = 4:SIG_LENGTH
    Xk = [Received_signal(index, :), Received_signal(index-1, :), Received_signal(index-2, :), Received_signal(index-3, :)]';
    Rxx = lambda*Rxx + (1-lambda)*(Xk*Xk'); % Adaptive learning of Rxx
    Y(index) = W'*Xk;
    if mod(index, 70)==0
        %mu = 0.99 * mu;  % Learning rate decline😅
        W = P * (W - mu * Rxx * Xk) + F;
    end
    Power(index) = W'*Rxx*W;
    Mismatch(index) = norm(C'*W-F_curly);
end


figure;
plot(Power, '.')
