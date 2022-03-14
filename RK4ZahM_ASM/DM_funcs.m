function [F,q1,q2,p1,p2,P0,V]=DM_funcs(t,R,lambda,N1,N2,Ri_sq,rmax_sq,rmax,P1,P2,P01,Pmin)

k2=0.1;
k1=0.35.*(t<5)+0.06.*(t>=5);
% k1=0.06;
k7=0.005;
g1=2.*k7;
g2=20.*g1;
g3=3.*g1;
fp1=0.88;
gp1=0.22;
gp2=4.*(fp1+gp1);
gp3=3.*gp1;

f = 0.045;
%Using R as initial conditions
r(1)=R(1);
r(2)=R(2);
r(3)=R(3);
r(4)=R(4);
r(5)=R(5);
r(6)=R(6);
r(7)=R(7);
r(8)=R(8);
force=lambda*(r(2)+r(5));
rho =1;
gamma = 25;

P0=(((P01))*triangle(2*pi*f*t)+Pmin); %to be used for DM ASM coupled system


%Transmural pressure as a function of force.

Ptm=((P0)-(force./r(8)));

if Ptm<=0
    
    rad=sqrt(Ri_sq*(1-(Ptm/P1)).^-N1);
    
elseif Ptm>=0
    % a=(1-((1-a0).*(1-(Ptm./P2)).^(-N2)));
    rad=sqrt(rmax_sq-(rmax_sq-Ri_sq)*(1-(Ptm/P2)).^-N2);
    
end



drdt=rho*(rad-(r(8)));  % this is to be used for the orginal DM coupled system.
V=-((gamma)/(2*pi*rmax))*drdt;



%Mean for cdf
p1=r(2)/r(1);
p2=r(5)/r(4);

%Standard deviation for cdf
q1=(sqrt((r(3)/r(1))-((r(2)/r(1)).^(2))));
q2=(sqrt((r(6)/r(4))-((r(5)/r(4)).^(2))));



%r,phi and I values for 1st PDE M1_lambda
r0=-p1/q1;
r1=(1-p1)/q1;
%     rinf=(Inf-p1)/q1;


phi0=0.5*(1+erf((r0-p1)/(q1*sqrt(2))));


phi1=0.5*(1+erf((r1-p1)/(q1*sqrt(2))));

I0=-(exp(-((-p1./q1).^(2))/2))/(sqrt(2*pi));
I1=-(exp(-((((1-p1)./q1)).^(2))/2))/(sqrt(2*pi));
%     I2=-(exp(-(((Inf-p1)./q1)).^(2))/2)/(sqrt(2*pi));

%r,phi and I values for 2nd PDE M2_lambda
r20=-p2/q2;
r21=(1-p2)/q2;

phi20=0.5*(1+erf((r20-p2)/(q2*sqrt(2))));
phi21=0.5*(1+erf((r21-p2)/(q2*sqrt(2))));

I20=-(exp(-((-p2./q2).^(2))/2))/(sqrt(2*pi));
I21=-(exp(-((((1-p2)./q2)).^(2))/2))/(sqrt(2*pi));

%Functions for the rhs of the first PDE M1_lambda
      J0=phi0;
%     J01=phi1;
%     J0inf=phinf;

       J10=((p1.*phi0)+(q1.*I0));
       J11=(p1.*phi1)+(q1.*I1);
%     J12=(p1.*phinf)+(q1.*I2);

      J20=((p1.^(2)).*phi0)+((2*p1.*q1).*I0)+((q1.^(2)).*(phi0+(r0*I0)));
      J21=((p1.^(2)).*phi1)+((2*p1.*q1).*I1)+((q1.^(2)).*(phi1+(r1*I1)));
%     J22=((p1.^(2)))+((2*p1.*q1).*I2)+((q1.^(2)));

      J30=(p1.^(3).*phi0)+((3.*p1.^(2).*q1).*I0)+((3.*p1.*q1.^(2)).*((phi0)+(r0*I0)))+((q1.^(3).*(2+r0.^(2)).*I0));
      J31=(p1.^(3).*phi1)+((3.*p1.^(2).*q1).*I1)+((3.*p1.*q1.^(2)).*(phi1+(r1*I1)))+(q1.^(3).*(2+(r1.^(2))).*I1);
%     J32=(p1.^(3))+((3.*p1.^(2).*q1).*I2)+((3.*p1.*q1.^(2)));

%Functions defined for the RHS of the second PDE M2_lambda
      K0=phi20;
%     K01=phi21;
%     K0inf=phi2inf;

      K10=((p2.*phi20)+(q2.*I20));
      K11=(p2.*phi21)+(q2.*I21);
%     K12=(p2.*phi2inf)+(q2.*I2in);

      K20=((p2.^(2)).*phi20)+((2*p2.*q2).*I20)+((q2.^(2)).*(phi20+(r20*I20)));
      K21=((p2.^(2)).*phi21)+((2*p2.*q2).*I21)+((q2.^(2)).*(phi21+(r21*I21)));
%     K22=((p2.^(2)))+((2*p2.*q2).*I2in)+((q2.^(2)));

      K30=(p2.^(3).*phi20)+((3.*p2.^(2).*q2).*I20)+((3.*p2.*q2.^(2)).*((phi20)+(r20*I20)))+((q2.^(3).*(2+r20.^(2)).*I20));
      K31=(p2.^(3).*phi21)+((3.*p2.^(2).*q2).*I21)+((3.*p2.*q2.^(2)).*(phi21+(r21*I21)))+(q2.^(3).*(2+(r21.^(2))).*I21);
%     K32=(p2.^(3))+((3.*p2.^(2).*q2).*I2in)+((3.*p2.*q2.^(2)));


%Components for the matrix F that will represent each moment,
%M1_lambda and M2_lambda

      A0= ((fp1*(1-r(7)))/2)-(fp1*(J11-J10)*r(1));%-2*(fp1*(K11-K10)*r(4));
      A1=((fp1*(1-r(7)))/3)-(fp1*(J21-J20)*r(1));%-2*(fp1*(K21-K20)*r(4));
      A2=((fp1*(1-r(7)))/4)-(fp1*(J31-J30)*r(1));%-2*(fp1*(K31-K30)*r(4));

      B0=(gp2*J0)+gp1*(J11-J10)+(gp1+gp3)*(p1-J11);
      B1=(gp2*J10)+gp1*(J21-J20)+(gp1+gp3)*((p1.^(2)+q1.^(2))-J21);
      B2=(gp2*J20)+gp1*(J31-J30)+(gp1+gp3)*((p1.^(3)+3*p1*q1.^(2))-J31);

      C0=k1;
      C1=k1*p2;
      C2=k1*(p2.^(2)+q2.^(2));

      D0=k2;
      D1=k2*p1;
      D2=k2*(p1.^(2)+q1.^(2));

      E0=(g2*K0)+g1*(K11-K10)+(g1+g3)*(p2-K11);
      E1=(g2*K10)+g1*(K21-K20)+(g1+g3)*((p2.^(2)+q2.^(2))-K21);
      E2=(g2*K20)+g1*(K31-K30)+(g1+g3)*(p2.^(3)+(3*p2*q2.^(2))-K31);


% Putting together the matrix F for the RHS of the system of
% equations

F=[A0-B0*r(1)+C0*r(4)-k2*r(1);A1-B1*r(1)+C1*r(4)-k2*r(2)-V*r(1);A2-B2*r(1)+C2*r(4)-k2*r(3)-2*V*r(2);D0*r(1)-E0*r(4)-k1*r(4);D1*r(1)-E1*r(4)-k1*r(5)-V*r(4);D2*r(1)-E2*r(4)-k1*r(6)-2*V*r(5);-k1*r(7)+(1-r(7))*k2; rho*(rad-r(8))];
