%% LOS CHANNEL MODEL OWC
clear all;
theta= 70; % semi-angle at half power 
m=-log10(2)/log10(cosd(theta)); %Lambertian order of emission 
P_total=20; %tranmistted optical power by individeal LED 
Adet=1e-2; %detector physical area of a PD

% Optics parameters
Ts=1; %gain of an optical filter; ignore if no filter is used 
index=1.5; %refractive index of a lens at a PD; ignore if no lens is used 
FOV=60*pi/180; %FOV of a receiver
G_Con=(index^2)/sin(FOV); %gain of an optical concentrator

% Room dimension
lx=10; ly=10; % room dimension in meter
h=1; %the distance between source and receiver plane 
XT=4.5; YT=4.5;% position of LED (Mid Point); 
gridPoints = 20;
[XR, YR] = generateRoomGrid(gridPoints, lx, ly);

D1=sqrt((XR-XT(1,1)).^2+(YR-YT(1,1)).^2+h^2); % distance verctor from source 1
%XT=2.5; YT=2.5;% position of LED (Mid Point); 
%D1=sqrt((XR-XT(1,1)).^2+(YR-YT(1,1)).^2+h^2) + D1; % distance verctor from source 1
cosphi_A1=h./D1; % angle vector

%
H_A1=(m+1)*Adet.*cosphi_A1.^(m+1)./(2*pi.*D1.^2);
% channel DC gain for source 1
P_rec=P_total.*H_A1.*Ts.*G_Con; % received power from source 1; 
P_rec_dBm=10*log10(P_rec);

meshc(XR(1,:),YR(:,1),P_rec_dBm);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Received power (dBm)');
axis([0 lx 0 ly min(min(P_rec_dBm)) max(max(P_rec_dBm))]);


%%



