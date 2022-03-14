function [X0,T,R,rLum,force,stiffness,phosphorylation,m10,m20,nmp,nm]=RK4ZahM_ASM()

%N number of time points
%trange specifies the interval of time steps
%T time grid
%X0 space grid
%gamma=25 rho=1 s=airway generation lambda=force activation

trange = [0 200];
N = 2000; M = 1000;
xrange = [-30 30]; 
P0 = 25; Pmin = 2;
 lambda = 11;
%  lambda = 11.75;


T=linspace(trange(1),trange(2),N);
X0=linspace(xrange(1),xrange(2),M);
dt=(trange(2)-trange(1))/(N); %dt <= 0.1

A=xlsread('airways1.xlsx');
Ri_sq=A(:,14);
rmax_sq=A(:,16);
N1=A(:,11);
N2=A(:,12);
P1=A(:,13);
P2=A(:,18);
rmax=A(:,8);
s = 7;
% Run fsolve to find equilibria

%Initial conditions

M10=0.309120293063990;
M11=0.118325676501984;
M12=0.066855519431672;
M20=0.134379171029593;
M21=-0.045132660622671;
M22=0.101673063419951;
C=0.6;
r55=0.442685665740035;
R=[M10;M11;M12;M20;M21;M22;C;r55];

force = zeros(1,N);
stiffness = zeros(1,N);
phosphorylation = zeros(1,N);
m10 = zeros(1,N);
m11 = zeros(1,N);
m21 = zeros(1,N);
m20 = zeros(1,N);
nmp = zeros(1,N);
nm = zeros(1,N);

for i=1:N %loop over number of timesteps to update then run RK4
    
    
    t=trange(1)+(i)*dt;
 
    [F1,~,~,~,~]=DM_funcs(t,R,lambda,N1(s),N2(s),Ri_sq(s),rmax_sq(s),rmax(s),P1(s),P2(s),P0,Pmin);%,t1,t2);
    [F2,~,~,~,~]=DM_funcs(t+(dt/2),R+((dt/2)*F1),lambda,N1(s),N2(s),Ri_sq(s),rmax_sq(s),rmax(s),P1(s),P2(s),P0,Pmin);%,t1,t2);
    [F3,~,~,~,~]=DM_funcs(t+(dt/2),R+((dt/2)*F2),lambda,N1(s),N2(s),Ri_sq(s),rmax_sq(s),rmax(s),P1(s),P2(s),P0,Pmin);%,t1,t2);
    [F4,q1,q2,p1,p2]=DM_funcs(t+(dt),R+((dt)*F3),lambda,N1(s),N2(s),Ri_sq(s),rmax_sq(s),rmax(s),P1(s),P2(s),P0,Pmin);%,t1,t2);
  
    Rnew = R + (dt/6)*(F1+(2*F2)+(2*F3)+F4);
    R=Rnew;
    
% A rebuild of the distribution approximation for Namp (M10:M12) and
    
    
    n1=(R(1)/((sqrt(2*pi))*q1)).*exp(-((X0-p1).^(2))/(2*(q1.^(2))));
    n2=(R(4)/((sqrt(2*pi))*q2)).*exp(-((X0-p2).^(2))/(2*(q2.^(2))));
    
    stiffness(i)=R(1)+R(4);
    
    rLum(i)=R(8) ;
    nmp(i)=1-R(7)-R(1);
    nm(i)=1-R(4)-R(1)-nmp(i);
    m10(i)=R(1);
    m11(i)=R(2);
    m21(i)=R(5);
    m20(i)=R(4);
    phosphorylation(i)=(nmp(i)+m10(i));
    force(i)=lambda*(R(2)+R(5));

end

plot(T,rLum)
xlabel('Time (s)')
ylabel('Radius of airway lumen (cm)')
hold on
