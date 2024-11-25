function dCdt=odefun5step(t,C,diffz,Di,QT,mc,Hb,CT,T,P,QS,void,dp,n,tspan,Rid,COinput,O2input,CO2input,k5step,Ct,yCO2f,yCOf,yO2f)
CCOint=pchip(tspan,COinput,t);
CCO2int=pchip(tspan,CO2input,t);
CO2int=pchip(tspan,O2input,t);
Cvec=[CCOint;C(1:n-1);CO2int;C(n:2*n-2);CCO2int;C(2*n-1:end)];
i=1:n;
CCO=Cvec(i)';
CO2=Cvec(i+n)';
CCO2=Cvec(i+2*n)';
theCOa=Cvec(i+3*n)';
theO2a=Cvec(i+4*n)';
theOa=Cvec(i+5*n)';
theCO2a=Cvec(i+6*n)';
thev=1-theCOa-theO2a-theOa-theCO2a;
yCO2=CCO2/CT;
yCO=CCO/CT;
yO2=CO2/CT;
yHe=1-yCO2-yCO-yO2;
A=pi*Di^2/4;
densb=mc/(A*Hb);
us0=QT/A;
DCO2He=1*10^(-3)*T^1.75*sqrt(1/44+1/4)/(P/101325)/(2.88^(1/3)+26.9^(1/3))^2*10^(-4);
DCOHe=1*10^(-3)*T^1.75*sqrt(1/28+1/4)/(P/101325)/(2.88^(1/3)+18.9^(1/3))^2*10^(-4);
DO2He=1*10^(-3)*T^1.75*sqrt(1/32+1/4)/(P/101325)/(2.88^(1/3)+16.6^(1/3))^2*10^(-4);
DCO2CO=1*10^(-3)*T^1.75*sqrt(1/44+1/28)/(P/101325)/(26.9^(1/3)+18.9^(1/3))^2*10^(-4);
DCO2O2=1*10^(-3)*T^1.75*sqrt(1/44+1/32)/(P/101325)/(26.9^(1/3)+16.6^(1/3))^2*10^(-4);
DCOO2=1*10^(-3)*T^1.75*sqrt(1/28+1/32)/(P/101325)/(18.9^(1/3)+16.6^(1/3))^2*10^(-4);
DCO=(1-yCO)./(yCO2/DCO2CO+yO2/DCOO2+yHe/DCOHe);
DO2=(1-yO2)./(yCO/DCOO2+yCO2/DCO2O2+yHe/DO2He);
DCO2=(1-yCO2)./(yCO/DCO2CO+yO2/DCO2O2+yHe/DCO2He);
us=us0*CT*(1-yCO2f-yCOf-yO2f)./(CT-CCO-CCO2-CO2);
DLCO=(0.45+0.55*void)*DCO+2*0.5*dp/2*us/void;
DLO2=(0.45+0.55*void)*DO2+2*0.5*dp/2*us/void;
DLCO2=(0.45+0.55*void)*DCO2+2*0.5*dp/2*us/void;
%kenetic
k1=k5step(1);
kr1=k5step(2);
k2=k5step(3);
kr2=k5step(4);
k3=k5step(5);
kr3=k5step(6);
k4=k5step(7);
kr4=k5step(8);
k5=k5step(9);
kr5=k5step(10);
r1=k1*thev*Rid*T.*CCO-kr1*theCOa;
r2=k2*thev*Rid*T.*CO2-kr2*theO2a;
r3=k3*theO2a.*thev-kr3*theOa.^2;
r4=k4*theCOa.*theOa-kr4*thev.*theCO2a;
r5=k5*theCO2a-kr5*thev*Rid*T.*CCO2;
%
dCdt=zeros(length(C),1);
i=2;
dCdt(i-1)=DLCO(i)*1/diffz^2*(CCO(i+1)-2*CCO(i)+CCO(i-1))-us(i)/void*(CCO(i+1)-CCO(i-1))/2/diffz-CCO(i)/void*us(i)/(CT-CCO(i)-CCO2(i)-CO2(i))*(CCO(i+1)-CCO(i-1)+CCO2(i+1)-CCO2(i-1)+CO2(i+1)-CO2(i-1))/2/diffz-(1)/void*densb*r1(i);
dCdt(i-1+n-1)=DLO2(i)*1/diffz^2*(CO2(i+1)-2*CO2(i)+CO2(i-1))-us(i)/void*(CO2(i+1)-CO2(i-1))/2/diffz-CO2(i)/void*us(i)/(CT-CCO(i)-CCO2(i)-CO2(i))*(CCO(i+1)-CCO(i-1)+CCO2(i+1)-CCO2(i-1)+CO2(i+1)-CO2(i-1))/2/diffz-(1)/void*densb*r2(i);
dCdt(i-1+2*n-2)=DLCO2(i)*1/diffz^2*(CCO2(i+1)-2*CCO2(i)+CCO2(i-1))-us(i)/void*(CCO2(i+1)-CCO2(i-1))/2/diffz-CCO2(i)/void*us(i)/(CT-CCO(i)-CCO2(i)-CO2(i))*(CCO(i+1)-CCO(i-1)+CCO2(i+1)-CCO2(i-1)+CO2(i+1)-CO2(i-1))/2/diffz+(1)/void*densb*r5(i);
i=3:n-1;
dCdt(i-1)=DLCO(i)*1/diffz^2.*(CCO(i+1)-2*CCO(i)+CCO(i-1))-us(i)/void.*(3*CCO(i)-4*CCO(i-1)+CCO(i-2))/2/diffz-CCO(i)/void.*us(i)./(CT-CCO(i)-CCO2(i)-CO2(i)).*(3*CCO(i)-4*CCO(i-1)+CCO(i-2)+3*CO2(i)-4*CO2(i-1)+CO2(i-2)+3*CCO2(i)-4*CCO2(i-1)+CCO2(i-2))/2/diffz-(1)/void*densb*r1(i);
dCdt(i-1+n-1)=DLO2(i)*1/diffz^2.*(CO2(i+1)-2*CO2(i)+CO2(i-1))-us(i)/void.*(3*CO2(i)-4*CO2(i-1)+CO2(i-2))/2/diffz-CO2(i)/void.*us(i)./(CT-CCO(i)-CCO2(i)-CO2(i)).*(3*CCO(i)-4*CCO(i-1)+CCO(i-2)+3*CO2(i)-4*CO2(i-1)+CO2(i-2)+3*CCO2(i)-4*CCO2(i-1)+CCO2(i-2))/2/diffz-(1)/void*densb*r2(i);
dCdt(i-1+2*n-2)=DLCO2(i)*1/diffz^2.*(CCO2(i+1)-2*CCO2(i)+CCO2(i-1))-us(i)/void.*(3*CCO2(i)-4*CCO2(i-1)+CCO2(i-2))/2/diffz-CCO2(i)/void.*us(i)./(CT-CCO(i)-CCO2(i)-CO2(i)).*(3*CCO(i)-4*CCO(i-1)+CCO(i-2)+3*CO2(i)-4*CO2(i-1)+CO2(i-2)+3*CCO2(i)-4*CCO2(i-1)+CCO2(i-2))/2/diffz+(1)/void*densb*r5(i);
%end
dCdt(n-1)=DLCO(n)*1/diffz^2*(-2*CCO(n)+2*CCO(n-1))-(1)/void*densb*r1(n);
dCdt(2*n-2)=DLO2(n)*1/diffz^2*(-2*CO2(n)+2*CO2(n-1))-(1)/void*densb*r2(n);
dCdt(3*n-3)=DLCO2(n)*1/diffz^2*(-2*CCO2(n)+2*CCO2(n-1))+(1)/void*densb*r5(n);
i=1:n;
dCdt(i+3*n-3)=(r1-r4)'/Ct;
dCdt(i+4*n-3)=(r2-r3)'/Ct;
dCdt(i+5*n-3)=(2*r3-r4)'/Ct;
dCdt(i+6*n-3)=(r4-r5)'/Ct;
end