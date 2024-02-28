% BeamPattern Formation attempt
% Implemantation of Conventional BF and MVDR from 'Adaptive Filters',
% 18.1.4.
% 2024.02.07

% Conventional BF
M=10;
theta0=45*pi/180;
theta1=20*pi/180;
theta2=330*pi/180;
s0=exp(1j*pi*[0:M-1]'*sin(theta0)); % steering vectors. Column vector get.
s1=exp(1j*pi*[0:M-1]'*sin(theta1));
s2=exp(1j*pi*[0:M-1]'*sin(theta2));

theta = 2*pi*[0:0.001:0.999];
s_theta = exp(1j*pi*[0:M-1]'*sin(theta));

% w = s0 / M;
% Gain = w'*s_theta;
% figure(1)
% f=polarplot(theta,Gain);
% set(f,'LineWidth',1)


% MVDR implementation
P0=1;
P1=1;
P2=1;   % standard variance
sigma_nu=0.1;
X=[];
for n=1:100
    X=[X sqrt(P0)*randn*s0+sqrt(P1)*randn*s1+sqrt(P2)*randn*s2...
        +sigma_nu*randn(M,1)];
end         %simulated input signal
R=X*X'/100; % Covariance Matrix 
w_mvdr = (R \ s0)/(s0'*inv(R)*s0);
figure(1)
Gain_mvdr = w_mvdr'*s_theta;
f=polarplot(theta,Gain_mvdr);
set(f,'LineWidth',1)
figure(2)
plot(theta, Gain_mvdr)
xticks(0:pi/2:2*pi); % 设置刻度为0、π/2、π、3π/2、2π
xticklabels({'0', 'π/2', 'π', '3π/2', '2π'}); % 设置标签
