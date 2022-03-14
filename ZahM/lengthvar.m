function [alpha,V] =lengthvar(t,f)

L=1*(t<180)+(1+(0.01*sin(2*pi*f*t)))*(t>=180 && t<=245)+(1+(0.04*sin(2*pi*f*t)))*(t>245);
gamma=25;
alpha=-gamma*(L-1);

%This
%is the V to use for DM

V=0*(t<180)+(-2*pi*gamma*0.01*f*cos(2*pi*f*t))*(t>=180 && t<=245)+(-2*pi*gamma*0.04*f*cos(2*pi*f*t))*(t>245);