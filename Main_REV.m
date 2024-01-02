clear all
close all
warning('off');
REV_Functions;
theta_A=50/180*pi;
theta_R=40/180*pi;
% theta_A=50/180*pi;
% theta_R=40/180*pi;
theta=(theta_A+theta_R)/2;
% theta=theta_A
%% INPUT
h=5e-9; 
R=h*5;
Sw=0.15;
% Ccr=1000/5e9;
% E=1e9; v=0.3; m=0.001; Ccr=1/(E*m/2/(1-v^2));
Ccr=1e-6;
gamma=72e-3;
Kw   = 2.3e9;                % Bulk modulus of liquid, Pa
Kg   = 2e5; 
% Kg=Kw;
% Kw=Kg
pg   = Kg;
time=0:pi*1e-4:2*pi*3;
dStress_amp = 5e5
dStress      = dStress_amp*(sin(time)).^1; % sinusoidal wave
Solution=zeros(9,length(time));
[V,R1,R2,pcap]=Functions(theta,h,R,gamma);
Ns=h*Sw/V;
peff=pg-Ns*(2*pi*R2*gamma+pi*R2^2*pcap);
% [V,R1,R2,pcap]=Functions(theta,h,R,gamma)
Solution(1,1)=peff;
Solution(2,1)=h;
Solution(3,1)=V;
Solution(4,1)=R;
Solution(5,1)=theta;
Solution(6,1)=pg;
Solution(7,1)=pcap;
Solution(8,1)=R1;
Solution(9,1)=R2;
dS=1;
[IncrSol1]=Pinned(Ns,gamma,Kw,Kg,Ccr,dS,Solution(:,1));
Slope1=IncrSol1(2)/dS
dtheta=IncrSol1(5); dScr=(theta_R-theta_A)/2/dtheta*dS
[IncrSol2]=Slipping(Ns,gamma,Kw,Kg,Ccr,dS,Solution(:,1));
Slope2=IncrSol2(2)/dS
DSCR=(1-Slope1/Slope2)*dScr
for i=2:length(time)
    dS=dStress(i)-dStress(i-1);
    [IncrSol]=Pinned(Ns,gamma,Kw,Kg,Ccr,dS,Solution(:,i-1));
    theta_new=Solution(5,i-1)+IncrSol(5);
    if theta_new  <= theta_R || theta_new >= theta_A
         [IncrSol]=Slipping(Ns,gamma,Kw,Kg,Ccr,dS,Solution(:,i-1));
    end
    Solution(:,i)=Solution(:,i-1)+IncrSol;
    
end
% load
theta_ar=Solution(5,:);
h_ar=Solution(2,:);
R_ar=Solution(4,:);
pcap_ar=Solution(7,:);
p_ar=Solution(6,:);
p_eff=Solution(1,:);

Sign=1;


ind1=find(time<=2*pi);
ind2=find(time>2*pi);

figure(1)
subplot(611)
plot(time(ind1)/pi/2,Sign*dStress(ind1)/1e5,'--r','Linewidth',2), grid on
hold on
plot(time(ind2)/pi/2,Sign*dStress(ind2)/1e5,'-k','Linewidth',2), grid on
xlabel('time, cycle #','FontSize',12)
ylabel('\Delta\sigma, bars','FontSize',12)
subplot(612)
plot(time(ind1)/pi/2,R_ar(ind1)*1e9,'--r','Linewidth',2), grid on
hold on
plot(time(ind2)/pi/2,R_ar(ind2)*1e9,'-k','Linewidth',2), grid on
xlabel('time, cycle #','FontSize',12)
ylabel('R, nm','FontSize',12)
subplot(613)
plot(time(ind1)/pi/2,theta_ar(ind1)*180/pi,'--r','Linewidth',2)
hold on
plot(time(ind2)/pi/2,theta_ar(ind2)*180/pi,'-k','Linewidth',2), grid on
xlabel('time, cycle #','FontSize',12)
ylabel('\theta, degrees','FontSize',12)
subplot(614)
plot(time(ind1)/pi/2,h_ar(ind1)*1e9,'--r','Linewidth',2), grid on
hold on
plot(time(ind2)/pi/2,h_ar(ind2)*1e9,'-k','Linewidth',2), grid on
xlabel('time, cycle #','FontSize',12)
ylabel('h, nm','FontSize',12)
subplot(615)
plot(time(ind1)/pi/2,p_ar(ind1)/1e5,'--r','Linewidth',2), grid on
hold on
plot(time(ind2)/pi/2,p_ar(ind2)/1e5,'-k','Linewidth',2), grid on
xlabel('time, cycle #','FontSize',12)
ylabel('\Deltap_g_a_s, bars','FontSize',12)
subplot(616)
plot(time(ind1)/pi/2,(pcap_ar(ind1))/1e5,'--r','Linewidth',2), grid on
hold on
plot(time(ind2)/pi/2,(pcap_ar(ind2))/1e5,'-k','Linewidth',2), grid on
xlabel('time, cycle #','FontSize',12)
ylabel('p_c_a_p, bars','FontSize',12)


figure(2)
subplot(311)
plot(Sign*dStress(time>=2*pi)/1e5,Sign*(h_ar(time>=2*pi)-h)/h,'k','LineWidth',2)
ylabel('Crack Volume Strain','FontSize',13)
xlabel('\Delta\sigma, bars','FontSize',13)
hold on
plot(Sign*dStress(time<2*pi)/1e5,Sign*(h_ar(time<2*pi)-h)/h,'--r','LineWidth',2)
grid on
legend('Cycle #>1', 'Cycle #=1','FontSize',13)
xticks([-dStress_amp -dStress_amp/2 0 dStress_amp/2 dStress_amp]/1e5)

Cm=1/20e9;
%%
phi_cr=0.0005;
Strain_REV=Cm*dStress+phi_cr*(h_ar-h)/h;

subplot(312)
plot(dStress(ind2)/1e5,Strain_REV(ind2),'k','Linewidth',2), grid on
hold on, box on, grid on
plot(dStress(ind1)/1e5,Strain_REV(ind1),'--r','Linewidth',2), grid on
axis tight
% legend('Cycle #>1', 'Cycle #=1','FontSize',13)
ylabel('REV Volume Strain','FontSize',13)
xlabel('\Delta\sigma, bars','FontSize',13)
title('\phi_c_r=0.0005','FontSize',13)
xticks([-dStress_amp -dStress_amp/2 0 dStress_amp/2 dStress_amp]/1e5)

%%
phi_cr=0.00005;
Strain_REV=Cm*dStress+phi_cr*(h_ar-h)/h;

subplot(313)
plot(dStress(ind2)/1e5,Strain_REV(ind2),'k','Linewidth',2), grid on
hold on
plot(dStress(ind1)/1e5,Strain_REV(ind1),'--r','Linewidth',2), grid on
axis tight, box on, grid on
% legend('Cycle #>1', 'Cycle #=1','FontSize',13)
ylabel('REV Volume Strain','FontSize',13)
xlabel('\Delta\sigma, bars','FontSize',13)
title('\phi_c_r=0.00005','FontSize',13)
xticks([-dStress_amp -dStress_amp/2 0 dStress_amp/2 dStress_amp]/1e5)

X=1./(h/Ccr/gamma/cos((theta_A+theta_R)/2)/Sw)
Y=Kg*Ccr/(1-Sw)
