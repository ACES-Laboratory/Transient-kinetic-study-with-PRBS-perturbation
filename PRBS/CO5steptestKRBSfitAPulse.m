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
QS=10/10^6/60;
QN=100/10^6/60;
QT=QS+QN;
tao=(Di^2*pi/4*Hb)/QT;
CT=P/Rid/T;
%fraction
fCO=1;
fO2=0.05;
fCO2=0;
yCO2=fCO2*QS/QT;
yO2=fO2*QN/QT;
yCO=fCO*QS/QT;
%Inlet
tpr=0.01;
np=8;
te=110-(2*np-1)*tpr;
deltat0=1+(te/np)*rand(1,np*2);
deltat0(end)=te-sum(deltat0)+deltat0(end);
while deltat0(end)<=0
deltat0=1+(te/np)*rand(1,np*2);
deltat0(end)=te-sum(deltat0)+deltat0(end);
end
%deltat0=[6.51732656345635	2.04311763527987	4.29434868348367	2.69332312395537	3.52528381785293	4.29484811803620	6.72959844236475	1.68181739641035	13.3954205842742	13.9731090988654	7.74017756970232	7.71805029103026	5.63730964610978	13.3588643786225	6.07021986325695	10.1771847872991];
steady=30;
deltat=deltat0;
tRs=zeros(length(deltat)*2+2,1);
COins=zeros(length(deltat)*2+2,1);
for i=1:length(deltat)
tRs(2*i+1)=tRs(2*i)+deltat(i);
tRs(2*i+2)=tRs(2*i+1)+tpr;
end
tRs(2)=tRs(1)+tpr;
i=1:length(tRs);
COins(mod(i,4)==2)=1;
COins(mod(i,4)==3)=1;
tRs=tRs(1:end-1);
COins=COins(1:end-1);
nsample=1001;
tR=linspace(tRs(1),tRs(end),nsample);
COin=pchip(tRs,COins,tR);
tR=[0,tR+steady];
COin=[0,COin];
CO2in=zeros(1,nsample+1);
%O2in=COin;
O2in=ones(1,nsample+1);
%realinputwith dilution
COinput=yCO*CT*COin*QT./(QN+COin*QS);
CO2input=yCO2*CT*CO2in*QT./(QN+COin*QS);
O2input=yO2*CT*O2in*QT./(QN+COin*QS);
%space & time
n=101;
diffz=Hb/(n-1);
xs=linspace(0,Hb,n);
tspan=tR';
%% 
%kinetics
%reversible original
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
%%Equation
x05step=zeros(1,7*n-3);
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
figure()
plot(ts(2:end)-ts(2),CCO(2:end,1),'LineWidth',1.5,'Color',[0.9290 0.6940 0.1250])
hold on
plot(ts(2:end)-ts(2),CO2(2:end,1),'LineWidth',1.5,'Color',[0.4940 0.1840 0.5560])
%plot(ts,CCO2(:,1),'LineWidth',1.5)
%ylim([-0.1 round(1.1*max([CCO;CO2;CCO2],[],'all'),1)])
ylim([0,4.5])
%legend('CO inlet','O_2 inlet','CO_2 inlet')
lgd=legend('A inlet','B_2 inlet');
xlabel('t/s','FontSize',15)
ylabel('C/(mol/m^3)','FontSize',15)
lgd.FontSize=15;
ax = gca;
ax.FontSize=15;
hold off
figure()
plot(ts(2:end)-ts(2),CCO(2:end,n),'LineWidth',1.5)
hold on
plot(ts(2:end)-ts(2),CO2(2:end,n),'LineWidth',1.5)
plot(ts(2:end)-ts(2),CCO2(2:end,n),'LineWidth',1.5)
lgd=legend('A outlet','B_2 outlet','AB_2 outlet');
ylim([0 4])
xlabel('t/s','FontSize',15)
ylabel('C/(mol/m^3)','FontSize',15)
lgd.FontSize=15;
ax = gca;
ax.FontSize=15;
hold off
figure()
surf(ts(2:end)-ts(2),xs,CCO(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('C_{A}/(mol/m^3)')
figure()
surf(ts(2:end)-ts(2),xs,CO2(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('C_{B_2}/(mol/m^3)')
figure()
surf(ts(2:end)-ts(2),xs,CCO2(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('C_{AB}/(mol/m^3)')
figure()
surf(ts(2:end)-ts(2),xs,theCOa(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('\theta_{A}')
figure()
surf(ts(2:end)-ts(2),xs,theO2a(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('\theta_{B_2}')
figure()
surf(ts(2:end)-ts(2),xs,theOa(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('\theta_{B}')
figure()
surf(ts(2:end)-ts(2),xs,theCO2a(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('\theta_{AB}')
figure()
surf(ts(2:end)-ts(2),xs,theCO2a(2:end,:)'+theOa(2:end,:)'+theO2a(2:end,:)'+theCOa(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('\theta_{T}')
%% 
Cexp=[awgn(CCO(2:end,n),35),awgn(CO2(2:end,n),35),awgn(CCO2(2:end,n),35)];
figure()
plot(ts(2:end)-ts(2),Cexp(:,1),'LineWidth',1.5)
hold on
plot(ts(2:end)-ts(2),Cexp(:,2),'LineWidth',1.5)
plot(ts(2:end)-ts(2),Cexp(:,3),'LineWidth',1.5)
lgd=legend('Predicted A outlet','Predicted B_2 outlet','Predicted AB outlet');
%ylim([0 3.5])
xlabel('t/s','FontSize',15)
ylabel('C/(mol/m^3)','FontSize',15)
lgd.FontSize=15;
ax = gca;
ax.FontSize=15;
hold off
%
optionsop=optimset('Display','iter','Algorithm','trust-region-reflective','UseParallel',true,'TolX',1e-7,'MaxFunEvals',100000,'MaxIter',100000,'TolFun',1e-8,'FinDiffRelStep',0.0001);
ko=[k5step(1:7),k5step(9:10),Ct];
k0=0.1*ones(1,10);
lb=zeros(1,10);
ub=ones(1,9)*inf;
ub=[ub,10];
[kfit,resnorm,residual,exitflag,output,lambda,J] = lsqnonlin(@(k)objlqsp(k,diffz,Di,QT,mc,Hb,CT,T,P,QS,void,dp,n,tspan,Rid,COinput,O2input,CO2input,yCO2,yCO,yO2,x05step,Cexp,Kr),k0,lb,ub,optionsop);
%% 
k5stepf=[kfit(1:7),NaN,kfit(8:9)];
k5stepf(1)=k5stepf(1)/10^5;
k5stepf(3)=k5stepf(3)/10^5;
k5stepf(10)=k5stepf(10)/10^5;
k5stepf(8)=k5stepf(7)*k5stepf(1)/k5stepf(2)*k5stepf(9)/k5stepf(10)*sqrt(k5stepf(3)/k5stepf(4)*k5stepf(5)/k5stepf(6)/Kr);
Ctf=kfit(10);
cif=nlparci(kfit,residual(1:length(Cexp)*3),'jacobian',J(1:length(Cexp)*3,:));
cif(1,:)=cif(1,:)/10^5;
cif(3,:)=cif(3,:)/10^5;
cif(9,:)=cif(9,:)/10^5;
kof=[k5stepf(1:7),k5stepf(9:10),Ctf];
error=kof-ko
relerror=error./ko
[ts,Cprof]=ode23tb(@(t,C)odefun5step(t,C,diffz,Di,QT,mc,Hb,CT,T,P,QS,void,dp,n,tspan,Rid,COinput,O2input,CO2input,k5stepf,Ctf,yCO2,yCO,yO2),tspan,x05step,options);
Cprofilef=[COinput',Cprof(:,1:n-1),O2input',Cprof(:,n:2*n-2),CO2input',Cprof(:,2*n-1:end)];
i=1:n;
CCOf=Cprofilef(:,i);
CO2f=Cprofilef(:,i+1*n);
CCO2f=Cprofilef(:,i+2*n);
theCOaf=Cprofilef(:,i+3*n);
theO2af=Cprofilef(:,i+4*n);
theOaf=Cprofilef(:,i+5*n);
theCO2af=Cprofilef(:,i+6*n);
figure()
plot(ts(2:end)-ts(2),Cexp(:,1),'LineWidth',1.5)
hold on
plot(ts(2:end)-ts(2),CCOf(2:end,n),'LineWidth',1.5)
%ylim([0 round(1.1*max([CCO;CCOf],[],'all'),1)])
xlabel('t/s','FontSize',15)
ylabel('C/(mol/m^3)','FontSize',15)
lgd=legend('Synthetic A outlet','Fitted A outlet');
lgd.FontSize=15;
ax = gca;
ax.FontSize=15;
hold off
figure()
plot(ts(2:end)-ts(2),Cexp(:,2),'LineWidth',1.5)
hold on
plot(ts(2:end)-ts(2),CO2f(2:end,n),'LineWidth',1.5)
%ylim([0 round(1.1*max([CO2;CO2f],[],'all'),1)])
xlabel('t/s','FontSize',15)
ylabel('C/(mol/m^3)','FontSize',15)
lgd=legend('Synthetic B_2 outlet','Fitted B_2 outlet');
lgd.FontSize=15;
ax = gca;
ax.FontSize=15;
hold off
figure()
plot(ts(2:end)-ts(2),Cexp(:,3),'LineWidth',1.5)
hold on
plot(ts(2:end)-ts(2),CCO2f(2:end,n),'LineWidth',1.5)
lgd=legend('Synthetic AB outlet','Fitted AB outlet');
%ylim([0 round(1.1*max([CCO2 CCO2f],[],'all'),2)])
xlabel('t/s','FontSize',15)
ylabel('C/(mol/m^3)','FontSize',15)
lgd.FontSize=15;
ax = gca;
ax.FontSize=15;
hold off
figure()
surf(ts(2:end)-ts(2),xs,CCOf(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('C_{A}/(mol/m^3)')
figure()
surf(ts(2:end)-ts(2),xs,CO2f(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('C_{B_2}/(mol/m^3)')
figure()
surf(ts(2:end)-ts(2),xs,CCO2f(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('C_{AB}/(mol/m^3)')
figure()
surf(ts(2:end)-ts(2),xs,theCOaf(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('\theta_{A}')
figure()
surf(ts(2:end)-ts(2),xs,theO2af(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('\theta_{B_2}')
figure()
surf(ts(2:end)-ts(2),xs,theOaf(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('\theta_{B}')
figure()
surf(ts(2:end)-ts(2),xs,theCO2af(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('\theta_{AB}')
figure()
surf(ts(2:end)-ts(2),xs,theCO2af(2:end,:)'+theOaf(2:end,:)'+theO2af(2:end,:)'+theCOaf(2:end,:)')
shading interp
xlabel('t/s')
ylabel('x/m')
zlabel('\theta_{T}')
figure()
plot(ts(2:end)-ts(2),CCOf(2:end,1),'LineWidth',1.5)
hold on
plot(ts(2:end)-ts(2),CO2f(2:end,1),'LineWidth',1.5)
%plot(ts,CCO2(:,1),'LineWidth',1.5)
ylim([0 round(1.1*max([CCOf;CO2f;CCO2f],[],'all'),1)])
%legend('CO inlet','O_2 inlet','CO_2 inlet')
legend('Fitted A inlet','Fitted B_2 inlet')
xlabel('t/s')
ylabel('C/(mol/m^3)')
hold off
figure()
plot(ts(2:end)-ts(2),CCOf(2:end,n),'LineWidth',1.5)
hold on
plot(ts(2:end)-ts(2),CO2f(2:end,n),'LineWidth',1.5)
plot(ts(2:end)-ts(2),CCO2f(2:end,n),'LineWidth',1.5)
legend('Fitted A outlet','Fitted B_2 outlet','Fitted AB outlet')
ylim([0 round(1.1*max([CCOf;CO2f;CCO2f],[],'all'),1)])
xlabel('t/s')
ylabel('C/(mol/m^3)')
hold off
save('kfitPRBSfitAPnoise.mat')