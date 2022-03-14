function [X0,T,R,force,stiffness,phosphorylation,m10,m20,nmp,nm]=ZahM()
%N number of time points
%trange specifies the interval of time steps
%T time grid
%X0 space grid
trange = [0 300];
N = 3000; M = 1000;
xrange = [-30 30]; f = 0.33;

T=linspace(trange(1),trange(2),N);

X0=linspace(xrange(1),xrange(2),M);

%dt=timestep
dt=(trange(2)-trange(1))/(N);

force = zeros(1,N);
stiffness = zeros(1,N);
phosphorylation = zeros(1,N);
m10 = zeros(1,N);
m11 = zeros(1,N);
m21 = zeros(1,N);
m20 = zeros(1,N);
nmp = zeros(1,N);
nm = zeros(1,N);


% Run fsolve to find equilibria

M10=0.05;
M11=0.001;
M12=0.1;
M20=0.05;
M21=0.001;
M22=0.1;
C=1;

R=[M10;M11;M12;M20;M21;M22;C];


for i=1:N %loop over number of timesteps to update then run euler
    
    
    t=trange(1)+i*dt;
    X=X0;
    
    [F1,~,~,~,~]=DM_funcs_1(t,R,f);
    [F2,~,~,~,~]=DM_funcs_1(t+(dt/2),R+((dt/2)*F1),f);
    [F3,~,~,~,~]=DM_funcs_1(t+(dt/2),R+((dt/2)*F2),f);
    [F4,q1,q2,p1,p2]=DM_funcs_1(t+dt,R+(dt*F3),f);
    
    
    
    Rnew = R + (dt/6)*(F1+(2*F2)+(2*F3)+F4);
    R=Rnew;
    
    % A rebuild of the distribution approximation for
    n1=(R(1)/((sqrt(2*pi))*q1)).*exp(-((X-p1).^(2))/(2*(q1.^(2))));
    n2=(R(4)/((sqrt(2*pi))*q2)).*exp(-((X-p2).^(2))/(2*(q2.^(2))));
    nmp(i)=1-R(7)-R(1);
    nm(i)=1-R(4)-R(1)-nmp(i);
    m10(i)=R(1);
    m11(i)=R(2);
    m21(i)=R(5);
    m20(i)=R(4);
    force(i)=R(2)+R(5);
    stiffness(i)=R(1)+R(4);
    phosphorylation(i)=(nmp(i)+m10(i));
    
    
end
figure(1)
plot(T,force,T,stiffness,T,phosphorylation)
legend('Force','Stiffness','Phosphorylation')

figure(2)
plot(T,m10,T,m20,T,nm,T,nmp)
legend('AMp','AM','M','Mp')

