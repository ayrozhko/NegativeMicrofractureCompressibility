function [dV_dh,dV_dR,dV_dtheta,V]=ParallelPlates_Functions(theta,h,R)
V=2*pi*(h^3/8*(1/cos(theta)^2+(2*R/h+tan(theta))^2)-h^3/24-h^3/16*1/cos(theta)^2*(2*R/h+tan(theta))*(sin(2*theta)+(pi-2*theta)));
dV_dh=pi*(((5*h^2-4*h*pi*R+4*R^2+8*h*R*theta-(h^2-4*R^2 )*cos(2*theta)+2*h*R*sin(3*theta)/cos(theta)+h^2*(2*R/h+6*theta-3*pi)*tan(theta)))/(8*cos(theta)^2 ));
dV_dR=pi*((-h^2*pi+4*h*R+2*h^2*theta+4*h*R*cos(2*theta)+h^2*sin(2*theta))/(2*(1+cos(2*theta))));
dV_dtheta=pi*(h^2*((-2*h*pi+4*R+4*h*theta+(4*R+h*(pi-2*theta))*cos(2*theta)+(3*h-2*pi*R+4*R*theta)*sin(2*theta)))/(8*cos(theta)^4 ));
end