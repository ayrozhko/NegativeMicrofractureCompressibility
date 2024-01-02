clear all
close all
theta_A=72/180*pi;
theta_R=40/180*pi;
% theta=(theta_A+theta_R)/2;
theta=50/180*pi;
h=2.2e-5; 
R=h*2;
dh=1e-5;
Rmin1=h/2*(1/cos(theta)-tan(theta));
Rmin2=h/2*(2/cos(theta)-tan(theta));
[R Rmin1 Rmin2]
gamma=72e-3;
time=0:1e-4:2*pi*5;
h_ar=h-dh*sin(time).^1;
theta_ar=h_ar*0; theta_ar(1)=theta;
R_ar=h_ar*0; R_ar(1)=R;
  V_ar=h_ar*0;
for i=2:length(h_ar)
    
    [dV_dh,dV_dR,dV_dtheta,V_ar(i-1)]=ParallelPlates_Functions(theta_ar(i-1),h_ar(i-1),R_ar(i-1));
    dh=h_ar(i)-h_ar(i-1);
    dtheta =-dh*dV_dh/dV_dtheta;
    dR=0;
    theta_new=theta_ar(i-1)+dtheta;
    
    if theta_new  <= theta_R || theta_new >= theta_A
        dtheta=0;
        dR=-dh*dV_dh/dV_dR;
    end
    theta_ar(i)=theta_ar(i-1)+dtheta;
    R_ar(i)=R_ar(i-1)+dR;
  
    
end

R1=h_ar./(2*cos(theta_ar));
R2=R_ar+h/2*(tan(theta_ar)-  1./cos(theta_ar));
pcap=gamma*(1./R1-1./R2);
Fcap=2*pi*R_ar*gamma.*sin(theta_ar)+pi*R_ar.^2*gamma.*(1./R1 -1./R2 );
ind1=find(time<=2*pi);
ind2=find(time>2*pi);
ind3=find(time>=0);
figure(1)
subplot(311)

plot(h_ar(ind2),theta_ar(ind2)*180/pi,'-k','Linewidth',2), grid on
hold on
plot(h_ar(ind1),theta_ar(ind1)*180/pi,'--r','Linewidth',2)
% hold on
xlabel('h, m')
ylabel('\theta, degrees')
subplot(312)
plot(h_ar(ind2),R_ar(ind2),'-k','Linewidth',2), grid on
hold on
plot(h_ar(ind1),R_ar(ind1),'--r','Linewidth',2), grid on
xlabel('h, m')
ylabel('R, m')
subplot(313)
plot(h_ar(ind2),Fcap(ind2),'-k','Linewidth',2), grid on
hold on
plot(h_ar(ind1),Fcap(ind1),'--r','Linewidth',2), grid on
xlabel('h, m')
ylabel('Fcap, N')


figure(2)
subplot(411)
plot(time(ind1)/pi/2,theta_ar(ind1)*180/pi,'--r','Linewidth',2)
hold on
plot(time(ind2)/pi/2,theta_ar(ind2)*180/pi,'-k','Linewidth',2), grid on
xlabel('time, cycle #')
ylabel('\theta, degrees')
subplot(412)
plot(time(ind1)/pi/2,R_ar(ind1),'--r','Linewidth',2), grid on
hold on
plot(time(ind2)/pi/2,R_ar(ind2),'-k','Linewidth',2), grid on
xlabel('time, cycle #')
ylabel('R, m')
subplot(413)
plot(time(ind1)/pi/2,Fcap(ind1),'--r','Linewidth',2), grid on
hold on
plot(time(ind2)/pi/2,Fcap(ind2),'-k','Linewidth',2), grid on
xlabel('time, cycle #')
ylabel('Fcap, N')
subplot(414)
plot(time(ind1)/pi/2,h_ar(ind1),'--r','Linewidth',2), grid on
hold on
plot(time(ind2)/pi/2,h_ar(ind2),'-k','Linewidth',2), grid on
xlabel('time, cycle #')
ylabel('h, m')
% figure(4)
% plot(h_ar(ind3),Fcap(ind3),'-','Linewidth',2), grid on
% hold on

% figure(3)
% k=9.287*0.5;
% Fel=(h_ar-h_ar(1))*k-Fcap(1);
% subplot(121)
% plot(h_ar(ind3),Fcap(ind3)+Fel(ind3),'-b','Linewidth',2), grid on
% xlabel('h, m')
% ylabel('Fcap+Fel, N')
% subplot(122)
% plot(h_ar(ind3),Fcap(ind3),'-k','Linewidth',2), grid on
% hold on
% plot(h_ar(ind3),Fel(ind3),'-r','Linewidth',2), grid on
% xlabel('h, m')
% ylabel('Fcap and Fel, N')
% 
% % y=-9.287*h+0.0334;
