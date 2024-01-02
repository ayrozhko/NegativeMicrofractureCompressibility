function REV_Functions
assignin('base','Derivatives',@Derivatives);
assignin('base','Functions',@Functions);
assignin('base','Pinned',@Pinned);
assignin('base','Slipping',@Slipping);
assignin('base','AssignValues',@AssignValues);
end
function [peff,h,V,R,theta,pg,pcap,R1,R2]=AssignValues(Sol)
peff=Sol(1,1);
h=Sol(2,1);
V=Sol(3,1);
R=Sol(4,1);
theta=Sol(5,1);
pg=Sol(6,1);
pcap=Sol(7,1);
R1=Sol(8,1);
R2=Sol(9,1);
end

function [dR1_dh,dR1_dtheta,dR2_dh,dR2_dtheta,dV_dh,dR2_dR,dV_dR,dV_dtheta]=Derivatives(theta,h,R)
% V=2*pi*(h^3/8*(1/cos(theta)^2+(2*R/h+tan(theta))^2)-h^3/24-h^3/16*1/cos(theta)^2*(2*R/h+tan(theta))*(sin(2*theta)+(pi-2*theta)));
dV_dh=pi*(((5*h^2-4*h*pi*R+4*R^2+8*h*R*theta-(h^2-4*R^2 )*cos(2*theta)+2*h*R*sin(3*theta)/cos(theta)+h^2*(2*R/h+6*theta-3*pi)*tan(theta)))/(8*cos(theta)^2 ));
dV_dR=pi*((-h^2*pi+4*h*R+2*h^2*theta+4*h*R*cos(2*theta)+h^2*sin(2*theta))/(2*(1+cos(2*theta))));
dV_dtheta=pi*(h^2*((-2*h*pi+4*R+4*h*theta+(4*R+h*(pi-2*theta))*cos(2*theta)+(3*h-2*pi*R+4*R*theta)*sin(2*theta)))/(8*cos(theta)^4 ));
dR2_dR=1;
dR2_dtheta=h/2/cos(theta)*(1/cos(theta)-tan(theta));
dR2_dh=1/2*(tan(theta)-1/cos(theta));
dR1_dtheta=h/2/cos(theta)*tan(theta);
dR1_dh=1/2/cos(theta);
end
%%
function [V,R1,R2,pcap]=Functions(theta,h,R,gamma)
V=2*pi*(h^3/8*(1/cos(theta)^2+(2*R/h+tan(theta))^2)-h^3/24-h^3/16*1/cos(theta)^2*(2*R/h+tan(theta))*(sin(2*theta)+(pi-2*theta)));
R1=h/2/cos(theta);
R2=R+h/2*(tan(theta)-1/cos(theta));
pcap=gamma*(1/R1-1/R2);
% peff=pg-Ns*(2*pi*R*gamma*sin(theta)+pi*R^2*pcap);
end
%%
function [IncrSol]=Pinned(Ns,gamma,Kw,Kg,Ccr,ds,Sol)
[~,h,V,R,theta,~,pcap,R1,R2]=AssignValues(Sol);
[dR1_dh,dR1_dtheta,dR2_dh,dR2_dtheta,dV_dh,~,~,dV_dtheta]=Derivatives(theta,h,R);

Matr         = zeros(8,8);
Matr(1,:)    = [ 0,     Kg,    -Ns*Kg, 0,         (h-Ns*V), 0,          0,          0                       ];
Matr(2,:)    = [ 0,     0, Kw,  0,     V,         -V,       0,                      0                       ];
Matr(3,:)    = [ 0,     dV_dh, -1,     dV_dtheta,  0,       0,          0,          0                       ];
Matr(4,:)    = [ 1,     0,      0,     0,         -1,       pi*R2^2*Ns, 0,          2*pi*(gamma+pcap*R2)*Ns ];
Matr(5,:)    = [-Ccr*h, 1,      0,     0,          0,       0,          0,          0                       ];
Matr(6,:)    = [ 0,     dR1_dh, 0,     dR1_dtheta, 0,       0,         -1,          0                       ];
Matr(7,:)    = [ 0,     dR2_dh, 0,     dR2_dtheta, 0,       0,          0,         -1                       ];
Matr(8,:)    = [ 0,     0,      0,     0,          0,       1,          gamma/R1^2, gamma/R2^2              ];
% 
Vector       = [0, 0, 0, 0 , Ccr*h*ds, 0, 0, 0]';
dSol         = (Matr)\Vector;

dpeff  = dSol(1);    dh     = dSol(2);
dV     = dSol(3);    dtheta = dSol(4);
dpg    = dSol(5);    dpcap  = dSol(6);
dR1    = dSol(7);    dR2    = dSol(8);
dR     = 0;
IncrSol=zeros(9,1);
IncrSol(1) = dpeff;   IncrSol(2) = dh;
IncrSol(3) = dV;      IncrSol(4) = dR;
IncrSol(5) = dtheta;  IncrSol(6) = dpg;
IncrSol(7) = dpcap;   IncrSol(8) = dR1;
IncrSol(9) = dR2;
end
%%
function [IncrSol]=Slipping(Ns,gamma,Kw,Kg,Ccr,ds,Sol)
[~,h,V,R,theta,~,pcap,R1,R2]=AssignValues(Sol);
[dR1_dh,~,dR2_dh,~,dV_dh,dR2_dR,dV_dR,~]=Derivatives(theta,h,R);
Matr         = zeros(8,8);
Matr(1,:)    = [ 0,     Kg,    -Ns*Kg, 0,         (h-Ns*V), 0,          0,          0                       ];
Matr(2,:)    = [ 0,     0, Kw,  0,     V,         -V,       0,                      0                       ];
Matr(3,:)    = [ 0,     dV_dh, -1,     dV_dR,      0,       0,          0,          0                       ];
Matr(4,:)    = [ 1,     0,      0,     0,         -1,       pi*R2^2*Ns, 0,          2*pi*(gamma+pcap*R2)*Ns ];
Matr(5,:)    = [-Ccr*h, 1,      0,     0,          0,       0,          0,          0                       ];
Matr(6,:)    = [ 0,     dR1_dh, 0,     0,          0,       0,         -1,          0                       ];
Matr(7,:)    = [ 0,     dR2_dh, 0,     dR2_dR,     0,       0,          0,         -1                       ];
Matr(8,:)    = [ 0,     0,      0,     0,          0,       1,          gamma/R1^2, gamma/R2^2              ];
% 
Vector       = [0, 0, 0, 0 , Ccr*h*ds, 0, 0, 0]';
dSol         = (Matr)\Vector;

dpeff  = dSol(1);    dh     = dSol(2);
dV     = dSol(3);    dR     = dSol(4);
dpg    = dSol(5);    dpcap  = dSol(6);
dR1    = dSol(7);    dR2    = dSol(8);
dtheta = 0;
IncrSol=zeros(9,1);
IncrSol(1) = dpeff;   IncrSol(2) = dh;
IncrSol(3) = dV;      IncrSol(4) = dR;
IncrSol(5) = dtheta;  IncrSol(6) = dpg;
IncrSol(7) = dpcap;   IncrSol(8) = dR1;
IncrSol(9) = dR2;
end
