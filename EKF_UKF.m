clear

c_v=10;
c_c=7.475;
m=3;
k=0.7;

T=0.01;
EndTime = 19.99;
time=0:T:EndTime;
% N = EndTime/T+1;
N=2000;
n=3;
R=0.1;
Q=diag([1e-5,1e-5,1e-5]);
B=[0; 1/m; 0];

t=(T*(0:(N-1)))';
u=@(t) 4*sawtooth(t*sqrt(2))+10*sin(t);

dxdt = @(t,x)[x(2); -x(3)*x(2)/m-k*x(1)/m; 0]+[0; 1/m; 0] * u(t);
f=c2d_rk4(dxdt,T); 
% f = @(t,x) [x(1)+T*x(2);x(2)+(-x(3)*x(2)/m-k*x(1)/m);x(3)]+[0; 1/m; 0]*u(t);
h    = @(x) x(1);

A = @(x) [1 T 0;
    -(k/m)*T 1-c_v/m*T -x(2)/m*T;
    0 0 1];

C = @(x) [1; 0; 0];

x=zeros(n,3);
y0=zeros(n,1);

x(1,:)=[0;0;10];
y0(1)=h(x(1,:));

for i=2:N
    x(i,:)=f((i-1)*T,x(i-1,:)');
    y0(i,:)=h(x(i,:));
end
w=randn(N,1)*sqrtm(R);

y=y0+w;

% xhat_ekf(1,:) = x(1,:);
% xhat_ekf(1,:) = x(1,:);
xhat_ekf(1,:)=[0; 0; 0.1*c_v];
yhat_ekf(1,:)=h(xhat_ekf(1,:));

P_ekf=diag([10,10,10]);

for i=2:N
    [xhat_ekf(i,:),P_ekf] = ...
        ekf(@(x) f((i-1)*T,x),h,A,B,C,Q,R,y(i,:),xhat_ekf(i-1,:),P_ekf);
end

xhat_ukf=zeros(N,3);
yhat_ukf=zeros(N,1);

xhat_ukf(1,:)=[0; 0; 0.1*c_v];
yhat_ukf(1,:)=h(xhat_ukf(1,:));
P_ukf=diag([10,10,10]);

for i=2:N
    [xhat_ukf(i,:),P_ukf] = ...
        ukf(@(x) f((i-1)*T,x),h,1,Q,R,y(i,:),xhat_ukf(i-1,:),P_ukf);
    yhat_ukf(i,:)=h(xhat_ukf(i,:));
end

figure(1),clf
for p=1:3
    subplot(3,1,p)
    plot(time,x(:,p));
    xlabel('Time[s]'),ylabel(sprintf('x%d',p))
end

figure(2),clf
for p=1:3
    subplot(3,1,p)
    plot(time,x(:,p),'k',...
        time,xhat_ekf(:,p),'b-.',...
        time,xhat_ukf(:,p),'r:');
    xlabel('Time[s]'),ylabel(sprintf('x%d',p))
    legend('true','ekf','ukf')
end
% % % ylabels = {'Position','Velocity','Parameter c_v'};
% % % for p=1:3
% % %     subplot(3,1,p);
% % %     plot(t,x(:,p),'r',t,xhat(:,p),'b');
% % %     xlim([min(t) max(t)]);
% % %     ylabel(ylabels{p});
% % %     xlabel('Time[s]');
% % %     legend('true','estimated','Location','SouthEast');
% % % end
