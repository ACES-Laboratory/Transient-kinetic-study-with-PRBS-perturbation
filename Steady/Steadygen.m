clear all
clc
%Operation Condition
mc=66.8/10^6;
void=0.5716;
dp=0.1/1000;
Di=3.9/1000;
Hb=10/1000;
P=1.06*101325;
Rid=8.314;
T=150+273.15;
QS=0/10^6/60;
QN=110/10^6/60;
QT=QS+QN;
tao=(Di^2*pi/4*Hb)/QT;
CT=P/Rid/T;
%fraction
fCOprofile=linspace(0.01,0.1,10);
fO2profile=linspace(0.01,0.1,10);
COexp=zeros(length(fCOprofile),length(fO2profile));
O2exp=zeros(length(fCOprofile),length(fO2profile));
CO2exp=zeros(length(fCOprofile),length(fO2profile));
fCO2=0;
%Inlet
te=500;
np=8;
%deltat0=1+(te/np)*rand(1,np*2);
%deltat0(end)=te-sum(deltat0)+deltat0(end);
%while deltat0(end)<=0
%deltat0=1+(te/np)*rand(1,np*2);
%deltat0(end)=te-sum(deltat0)+deltat0(end);
%end
deltat0=[600];
deltat=deltat0;
tRs=zeros(length(deltat)*2+2,1);
for i=1:length(deltat)
tRs(2*i+1)=tRs(2*i)+deltat(i);
tRs(2*i+2)=tRs(2*i+1)+0.01;
end
tRs(2)=tRs(1)+0.01;
COins=ones(1,length(tRs));
tRs=tRs(1:end-1);
COins=COins(1:end-1);
ns=1000;
tR=linspace(tRs(1),tRs(end),ns);
COin=pchip(tRs,COins,tR);
CO2in=COin;
%O2in=COin;
O2in=COin;
%kinetics
k5step=[3.5*10^-5*5 0.1*5 3.5*10^-5*5 0.1*5 0.8*5 0.08*5 0.1 NaN*5 0.1 3.5*10^-5];
Ct=0.4;
%NIST: O2:100-700K CO:298-1300K, CO2:298-1200
delHCO=(25.56759*T/1000+6.096130*(T/1000)^2/2+4.054656*(T/1000)^3/3-2.671301*(T/1000)^4/4-0.131021/(T/1000)-118.0089+110.5271)*1000;	
SCOs1=(25.56759*log(T/1000)+6.096130*(T/1000)+4.054656*(T/1000)^2/2-2.671301*(T/1000)^3/3-0.131021/2/(T/1000)^2+227.3665);
SCOs2=(25.56759*log(298.15/1000)+6.096130*(298.15/1000)+4.054656*(298.15/1000)^2/2-2.671301*(298.15/1000)^3/3-0.131021/2/(298.15/1000)^2+227.3665);
delGCO=delHCO-(T*SCOs1-298.15*SCOs2);
delHO2=(31.32234*T/1000-20.23531*(T/1000)^2/2+57.86644*(T/1000)^3/3-36.50624*(T/1000)^4/4+0.007374/(T/1000)-8.903471)*1000;	
SO2s1=(31.32234*log(T/1000)-20.23531*(T/1000)+57.86644*(T/1000)^2/2-36.50624*(T/1000)^3/3+0.007374/2/(T/1000)^2+246.7945);
SO2s2=(31.32234*log(298.15/1000)-20.23531*(298.15/1000)+57.86644*(298.15/1000)^2/2-36.50624*(298.15/1000)^3/3+0.007374/2/(298.15/1000)^2+246.7945);
delGO2=delHO2-(T*SO2s1-298.15*SO2s2);
delHCO2=(24.99735*T/1000+55.18696*(T/1000)^2/2-33.69137*(T/1000)^3/3+7.948387*(T/1000)^4/4+0.136638/(T/1000)-403.6075+393.5224)*1000;	
SCO2s1=(24.99735*log(T/1000)+55.18696*(T/1000)-33.69137*(T/1000)^2/2+7.948387*(T/1000)^3/3+0.136638/2/(T/1000)^2+228.2431);
SCO2s2=(24.99735*log(298.15/1000)+55.18696*(298.15/1000)-33.69137*(298.15/1000)^2/2+7.948387*(298.15/1000)^3/3+0.136638/2/(298.15/1000)^2+228.2431);
delGCO2=delHCO2-(T*SCO2s1-298.15*SCO2s2);
GR=(-2*394.373+2*137.168)*1000-(2*delGCO+delGO2)+2*delGCO2
Kr=exp(-GR/Rid/T)/10^5
k5step(8)=k5step(7)*k5step(1)/k5step(2)*k5step(9)/k5step(10)*sqrt(k5step(3)/k5step(4)*k5step(5)/k5step(6)/Kr);
tspan=tR';
n=101;
diffz=Hb/(n-1);
xs=linspace(0,Hb,n);
x05step=zeros(1,7*n-3);
%fraction
for k=1:length(fCOprofile)
    for l=1:length(fO2profile)
       yCO2=fCO2*QN/QT;
       yO2=fO2profile(l)*QN/QT;
       yCO=fCOprofile(k)*QN/QT;
       COinput=yCO*CT*COin;
       CO2input=yCO2*CT*CO2in;
       O2input=yO2*CT*O2in;
       options=odeset('RelTol',1e-5);
       [ts,Cpro]=ode23tb(@(t,C)odefun5step(t,C,diffz,Di,QT,mc,Hb,CT,T,P,QS,void,dp,n,tspan,Rid,COinput,O2input,CO2input,k5step,Ct,yCO2,yCO,yO2),tspan,x05step,options);
       Cprofile=[COinput',Cpro(:,1:n-1),O2input',Cpro(:,n:2*n-2),CO2input',Cpro(:,2*n-1:end)];
       i=1:n;
       CCO=Cprofile(:,i);
       CO2=Cprofile(:,i+1*n);
       CCO2=Cprofile(:,i+2*n);
       theCOa=Cprofile(:,i+3*n);
       theO2a=Cprofile(:,i+4*n);
       theOa=Cprofile(:,i+5*n);
       theCO2a=Cprofile(:,i+6*n);
       A=pi*Di^2/4;
       us0=QT/A;
       densb=mc/(A*Hb);
       Cexp=[awgn(CCO(:,n),35);awgn(CO2(:,n),35);awgn(CCO2(:,n),35)];
       COexp(k,l)=mean(Cexp(length(ts)-50:length(ts)));
       O2exp(k,l)=mean(Cexp(2*length(ts)-50:2*length(ts)));
       CO2exp(k,l)=mean(Cexp(3*length(ts)-50:3*length(ts)));
    end
end
save('Steadyfit.mat')