% semi-angle at half power
theta=35;
%Lambertian order of emission
m=-log10(2)/log10(cosd(theta));

%Center luminous intensity y total según número de LEDs
I0=0.73;
I0_total_1=20*I0;
I0_total_2=30*I0;
% room dimension in metre
lx=10; ly=10;
%the distance between source and receiver plane
h=10;
% position of LED1
XTrans1=4; YTrans1=4;
XTrans2 = -4; YTrans2 = -4;
% number of grid in the receiver plane (cada 2 centímetros)
Nx=lx*20; Ny=ly*20;
x=-lx:lx/Nx:lx; 
y=-ly:ly/Ny:ly;
[XRec,YRec]=meshgrid(x,y);
% distance vector from source 1
Vector_Distancia_1=sqrt((XRec-XTrans1(1,1)).^2+(YRec-YTrans1(1,1)).^2+h^2);
Vector_Distancia_2=sqrt((XRec-XTrans2(1,1)).^2+(YRec-YTrans2(1,1)).^2+h^2);
% vector que contiene el angulo de irradiancia (del LED) para cada posición
coseno_phi_1=h./Vector_Distancia_1;
coseno_phi_2=h./Vector_Distancia_2;
% FLUJO LUMINOSO o ILUMINANCIA (en lux)
E_lux_1=(I0_total_1*coseno_phi_1.^m)./(Vector_Distancia_1.^2);
E_lux_2=(I0_total_2*coseno_phi_2.^m)./(Vector_Distancia_2.^2);
figure (1)
meshc(x,y,E_lux_1);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Illuminance(lx)');

figure (2)
meshc(x,y,E_lux_2);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Illuminance(lx)');

figure (3)
meshc(x,y,E_lux_1+E_lux_2);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Illuminance(lx)');
%%

%% LOS Channel
theta= 70; % semi-angle at half power 
m=-log10(2)/log10(cosd(theta)); %Lambertian order of emission 
Adet=1e-4; %detector physical area of a PD
h = 0.2; %heigh
angle_vector = -pi/2:0.1:pi/2; %angle in rad from led to plain
distance = h./cos(angle_vector); % hipt=alt/cos(ang)

% Optics parameters
Ts=1; %gain of an optical filter; ignore if no filter is used 
index=1.5; %refractive index of a lens at a PD; ignore if no lens is used 
FOV=60*pi/180; %FOV of a receiver

%Calculate H and R
G_Con = (index^2)/sin(FOV); %gain of an optical concentrator
R = (m+1)/(2*pi).*cos(angle_vector).^m;
A_eff = Adet*cos(angle_vector);
H_DC = R.*A_eff.*Ts.*G_Con./(distance.^2); %H_DC for different angles


