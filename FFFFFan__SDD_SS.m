
clc;
clear all;
close all;
tic
Qx = diag([200, 1000]);
Rx=diag([1, 10^-4]);
% F = diag([20, 20]);

mg=0.46;
mx=4.9;
my=8.5;
r=0.12;
J=0.05;
d=1.2;
g=9.81;

t0=0;
tf=20;
DeltaT=0.01;
t=0:0.01:9.999;

it= round((tf-t0)/DeltaT); %Number of steps during each simulation

% Memory allocation
x = zeros(6,it); % History of states
u_SD=zeros(2,it);
u_dis=zeros(2,it);

x(:,1) =  [-0.2 0.25 -pi/10 0.1 -pi/6 0]';% Initial state vector/ Initial condition;
n1=0;
n2=0;
N1=0;
N2=0;
M1=0;
u = zeros(2,it); % History of control
Area1=0;
Area2=0;
d1=0;
dxt=[0 0]';
for i=2:it
    
    t(i)=(i-1)*DeltaT;
    if t(i)>8 && t(i)<9
        d1=3*sin(t(i))^2+x(1,i-1)^2+x(1,i-1)^2;
        M1=M1+1;
        d1=sin(t(i));
        dxt=[d1 10*d1]';
    end
    
    Ax=[0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1
        0 0 -mg*g*sin(x(3,i-1))/(mx*x(3,i-1)) -d/mx 0 0
        0 0 mg*g*(cos(x(3,i-1))-1)/(my*x(3,i-1)) 0 -d/my 0
        0 0 0 0 0 0];
    
    Bx=[0 0;0 0;0 0;cos(x(3,i-1))/mx -sin(x(3,i-1))/mx;
        sin(x(3,i-1))/my cos(x(3,i-1))/my;r/J 0];
    Qc=ctrb(Ax,Bx);
    
    if rank(Qc)==6
        n1=n1+1;
    else
        n2=n2+1;
    end
    
    C=[1 0 0 0 0 0;0 1 0 0 0 0];
    P = care(Ax, Bx, C'*Qx*C, Rx);
    x_n=x(4:5,i-1);
    
    A=[-mg*g*sin(x(3,i-1))/(mx*x(4,i-1))-d/mx 0
        0 mg*g*(cos(x(3,i-1))-1)/(my*x(5,i-1))-d/my];
    
    B=[cos(x(3,i-1))/mx -sin(x(3,i-1))/mx
        sin(x(3,i-1))/my cos(x(3,i-1))/my];
    
    u_SD(:,i-1)=-Rx^-1*Bx'*P*x(:,i-1);
    
    A_n=A+B*diag([u_SD(1,i-1)/x(4,i-1) u_SD(2,i-1)/x(5,i-1)]);
    
    Qc=ctrb(A,B);
    if rank(Qc)==2
        n1=n1+1;
    else
        n2=n2+1;
    end
    
    P_n = care(A_n, B, Qx, Rx);
    
    %     P_dot=-A'*P-P*A-Qx+P*B*Rx^-1*B'*P;
    
    %% defintion sliding Surface
    
    Delta_XT=Qx-Rx;
    S_xt=B'*P_n;
    Sigma_xt=S_xt*x_n;
    delta_xt=sqrt(x_n'*Delta_XT*x_n);
    S_t= abs(Sigma_xt) <= delta_xt;
    %     k1=30;
    %     k2=10;
    k=0.3;
    
    %% Design of control law SSC
    if ~ismember(x_n,S_t)
        u(:,i-1)=u_SD(:,i-1)-k*tanh(Sigma_xt);
        %       u(:,i-1)=-(k1+k2*norm(Sigma_xt)^2*a-1)*tanh(Sigma_xt);
        N1=N1+1;
        %                 else
        %             Beta=0.1;
        %             k3=15;
        %             k4=0;
        %             h=B'*P*B;
        %             G=B'*(P_dot+P*A)*x_n;
        %             u(:,i-1)=-h^-1*(G+(k3*h+k4*norm(Sigma_xt)^2*Beta*-1)*tanh(Sigma_xt));
        %             N2=N2+1;
    end
    u=min(u,+50);
    u=max(u,-50);
    
    x(:,i) = x(:,i-1) + DeltaT*(Ax*x(:,i-1)+(Bx)*(u(:,i-1)+dxt));
    x(3,i-1)=max(x(3,i-1),-pi/20);
    x(3,i-1)=min(x(3,i-1),+pi/20);
%     if t(i)>18
%         x(1,i)=x(1,i)+0.005;
%     end
    Area1=Area1+abs(x(:,i-1))*DeltaT;
    Area2=Area2+abs(u(:,i-1))*DeltaT;
end
Area_X=Area1
Area_U=Area2
fprintf('Done \n')

figure;
subplot(2,2,1);
plot(t,x(1,1:it),'r','LineWidth',1);
xlabel('time(s)');
ylabel('The First State (x)');
% hleg1 = legend('x_1 : SDDRE+SMC');
grid on;

subplot(2,2,2);
plot(t,x(2,1:it),'b','LineWidth',1);
xlabel('time(s)');
ylabel('The second state (y)');
% hleg2 = legend('x_2 : SDDRE+SMC');
grid on;

subplot(2,2,3);
plot(t,x(4,1:it),'b','LineWidth',1);
xlabel('time(s)');
ylabel('pitch Angle(Teta)');
hleg3 = legend('x_3 : SDDRE+SMC');
grid on;

subplot(2,2,4);
plot(t,u(1,1:it),'r','LineWidth',1);
hold on;
plot(t,u(2,1:it),'LineWidth',1);
legend('u_1','u_2')
xlabel('time(s)');
ylabel('u');
hleg5 = legend('u_1','u_2');
grid on;
toc
% 
% subplot(2,2,1);
% plot(t,u(1,1:it),'r','LineWidth',1);
% grid on;
% 
% subplot(2,2,2);
% plot(t,u(2,1:it),'LineWidth',1);
% legend('u_1','u_2')
% xlabel('time(s)');
% ylabel('u');
% hleg5 = legend('u_1','u_2');
% grid on;
% toc
