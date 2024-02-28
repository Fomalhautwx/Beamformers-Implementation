M = 4;
f = 1000;   % Frequecy of impinging signal
d = 0.1;    % Gap between sensors;
c = 340;    % Speed of Sound
rate = 8000;% define sample rate
SIG_LENGTH = 2e4;

Wc = ones(M, 1) / M;                % Fixed Beamformer
Ws = [1,1,-1,-1;1,-1,-1,1;1,-1,1,1]; % Blocking Matrix
Al = randn(3, 1);                   % Random initialization of MC

theta0 = 0*pi/180;  % Target direction
theta1 = 30*pi/180; theta2 = -60*pi/180;
s0=exp(2j*pi*[0:M-1]'*f*d*sin(theta0)/c); % Steering vector of Target direction. Column vector get.
s1 = exp(2j*pi*[0:M-1]'*f*d*sin(theta1)/c); s2 = exp(2j*pi*[0:M-1]'*f*d*sin(theta2)/c);

theta = pi*[-0.5:0.001:0.5];
s_theta = exp(-2j*pi*[0:M-1]'*f*d*sin(theta)/c);   % Test directions

sig1 = 0.5 * randn(SIG_LENGTH, 1); % Look direction signal
bpFilt = designfilt('bandpassfir', 'FilterOrder', 100, ...
             'CutoffFrequency1', 0.12, 'CutoffFrequency2', 0.13,...
             'SampleRate', 1);
sig1 = filter(bpFilt,sig1);
input_signal = sig1*s0'+sig1*s1'+sig1*s2'+randn(SIG_LENGTH, 1);% 4-channel delayed input signal

PreSteering = exp(2j*pi*[0:M-1]'*f*d*sin(theta0)/c);
Steered = input_signal*diag(PreSteering);

mu = 0.01;  % Learning rate
Power = zeros(SIG_LENGTH, 1);
Target_power = sig1'*sig1 / SIG_LENGTH;
Yo = zeros(SIG_LENGTH, 1);

input = Steered(1, :);
% FBF
yc = input*Wc;
% BM
X_prime = Ws*input';
% MC
ya = Al'*X_prime;
yo = real(yc - ya);   % Yield the output
Power(1) = yo'*yo / Target_power;
lambda = 0.0005; % Coefficient of forget

for index=2:SIG_LENGTH
    input = Steered(index, :);
    % FBF
    yc = input*Wc;
    % BM
    X_prime = Ws*input';
    % MC
    ya = Al'*X_prime;
    yo = real(yc - ya);   % Yield the output
    if mod(index, 20) == 0
        Al = Al+mu*real(yo)*real(X_prime) / (X_prime'*X_prime); % NLMS
    end    
    Power(index) = lambda*(yo'*yo) / Target_power+(1-lambda)*Power(index-1);   % Output power normalized by target signal power.
    Yo(index) = yo;
end

%%Test
Xm = s_theta'*diag(PreSteering);
Output = real(Xm*Wc-(Al'*(Ws*Xm'))');

plot(theta, Output)
xticks(-pi/2:pi/4:pi/2); % 设置刻度为0、π/2、π、3π/2、2π
xticklabels({'-π/2', '-π/4', '0', 'π/4', 'π/2'}); % 设置标签

% figure;
% plot(Power, '.')
