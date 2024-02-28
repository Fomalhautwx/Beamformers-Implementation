%
% Direction of arrival based on the MVDR method
% 
M=10;
theta0=45*pi/180;
theta1=0*pi/180;
theta2=0*pi/180;
s0=exp(1j*pi*[0:M-1]'*sin(theta0)); % steering vectors. Column vector get.
s1=exp(1j*pi*[0:M-1]'*sin(theta1));
s2=exp(1j*pi*[0:M-1]'*sin(theta2));

P0=1;
P1=1;
P2=1;
sigma_nu=0.1;
X=[];
for n=1:100
    X=[X sqrt(P0)*randn*s0+sqrt(P1)*randn*s1+sqrt(P2)*randn*s2...
        +sigma_nu*randn(M,1)];
end         %simulated input signal
R=X*X'/100; % Covariance Matrix 
for n=1:100
    theta=2 * (n-1) * pi / 100;
    s=exp(1j*pi*[0:M-1]'*sin(theta));
    S(n)=1/(s'*inv(R)*s);
end
% figure,axes('position',[0.25 0.25 0.5 0.5])
% plot([-50:50],S,'k'),hold on
% plot(20*[1 1],[0 1],'--')
% text(20.5,0.9,'\theta_0')
% plot(25*[1 1],[0 1],'--')
% text(25.5,0.9,'\theta_1')
% plot(-30*[1 1],[0 1],'--')
% text(-29.5,0.9,'\theta_2')
% xlabel('\theta')5
% hold off

S_norm = normalize(real(S), 'range') + 1i * normalize(imag(S), 'range');
theta=2*pi*[0:0.01:0.99];
gain=abs(S_norm);
f=polarplot(theta,gain);
set(f,'LineWidth',1)
