clear all;
close all;
%%

theta=70;
m=-log10(2)/log10(cosd(theta));
Pt=20;
Adet=1e-4;
Ts=1;
index=1.5;
FOV=60*pi/180;
G_Con=(index^2)/sin(FOV);
h=10;
responsivity = 0.59;
%%
angleVector = 0:0.01:FOV;

cosPhi = cos(angleVector);

distanceVector = h./cosPhi;

H0 = ((Adet*(m+1))./(2*pi*distanceVector.^2)).*cosPhi.^(m+1)*Ts*G_Con;

%%
L = 2.^(1:1:5); %Niveles de modulación

RL = 50; %resistencia
Rb = 100*10^6;
kB = 1.380649e-23; %Cte. Boltzmann
Tk = 295; %Kelvin, temperatura absoluta
Gol = 10; %Open-loop voltage gain
Cpd = 112; %pF/cm^2 %Capacitancia por unidad de area
A = Adet; %Área del captador
I2 = 0.562; %Factor de ruido
B = 1./log2(L) * Rb; %Ancho de banda
n = 1.5; %FET channel factor
gm = 30; %FET transconductancia
I3 = 0.0868; %Factor de ruido
E_c = 1.602e-19; %Electronic charge

sky_irra = 1e-3; % at 850nm wavelength, in W/cm^2-um-sr
sun_irra = 550e-4; %at 850nm wavelength, in W/cm^2-um
OBP = 1e-3; %Optical filter bandwidth in micrometre
Isky = sky_irra*OBP*(4/pi)*FOV^2; %Sky irradiance
Isun = sun_irra*OBP; %Sun irradiance

varShot = 2*E_c*Rb;
varBackgroundRad = 2*E_c*Rb*responsivity*(Isky + Isun);
%varThermal = ((8*pi*kB*Tk)/Gol)*Cpd*A*I2*Rb.^2;
varThermal = 4*kB*Tk*Rb/RL;
%%
noise = varShot + varThermal + varBackgroundRad;

signal = responsivity^2*Pt^2*H0.^2;
snr = signal/noise;

%%

Ber_pam = [];
for i=1:length(L)
    
    Ber = (2*(L(i)-1))/(L(i)*log2(L(i))) * erfc(sqrt(log2(L(i))/(L(i)-1)^2)*(responsivity*H0*Pt)/sqrt(noise));
    Ber_pam = [Ber_pam;Ber];
    M{i} = sprintf('%d-PAM', L(i));
    semilogy(flip(angleVector), flip(Ber_pam(i,:)))
    hold on
end

legend(M);
xlim([0,angleVector(end)])
set ( gca, 'xdir', 'reverse' )

%% BER para un ángulo de 0.4rad
angulo = 1;
cosPhi2 = cos(angulo);

d = h/cosPhi2;
H = ((Adet*(m+1))./(2*pi*d)).*cosPhi2^(m+1)*Ts*G_Con;

P_total = 10:0.5:100;
signal2 = responsivity^2*P_total.^2*H^2;
snr2 = signal2/noise;

Ber_pam2 = [];
figure(2)
for i=1:length(L)
    Ber = (2*(L(i)-1))/(L(i)*log2(L(i))) * erfc(sqrt(log2(L(i))/(L(i)-1)^2)*(responsivity*H*P_total)/sqrt(noise));
    Ber_pam2 = [Ber_pam2;Ber];
    M{i} = sprintf('%d-PAM', L(i));
    semilogy(10*log10(snr2), Ber_pam2(i,:))
    hold on
end

legend(M);
xlim([10*log10(snr2(1)),10*log10(snr2(end))])
ylim([10^-10, 1])


%% Simulación para comparar con el BER para OOK(2-PAM) Y 0.4rad
Rb = 100*10^6;
Tb = 1/Rb;

% Generación de bits
bits = round(rand(1,ceil(20e-8/Tb)));
bits_length = length(bits);
samplesPerBit = 10;
Tsamp = Tb/samplesPerBit;
% Alex alvarado
bits_padded = upsample(bits, samplesPerBit);
%Tx_filter=rcosdesign(0.2,6,samplesPerBit);
Tx_filter = ones(1,samplesPerBit);
%Tx_filter = 1;
h = [1,zeros(1,length(Tx_filter)) ];
Rx_filter=flipud(Tx_filter);

filter_sys = conv(Tx_filter, h);
filter_overall_sys = conv(filter_sys, Rx_filter);

[i, delay] = max(filter_overall_sys);

%%

received_symbols = conv(filter_sys, bits_padded);
ber = [];
 gauss_noise = sqrt(noise) * randn([1,length(received_symbols)]);
for i=1:1:length(P_total)
    Power_rec = P_total*H;
    Eb = snr2(i)*noise;
    received_symbols_noise = awgn(Eb*received_symbols, 10*log10(snr2(i)), 'measured');
    overall_symbols_noise = conv(Rx_filter, received_symbols_noise);
    real_bits_received = overall_symbols_noise(delay:samplesPerBit:delay+(bits_length*samplesPerBit)-1);
   estimated = (real_bits_received./max(real_bits_received)>=0.5);
  % real_bits_received = downsample(overall_symbols_noise,samplesPerBit);
  % estimated = zeros(1,length(bits));
  % estimated(find(real_bits_received>Eb/2))=1;
     
   [nerr, ber(i)] = biterr(bits, estimated);
end

figure(3)
semilogy(10*log10(snr2), Ber_pam2(1,:))
hold on
semilogy(10*log10(snr2), ber)
xlim([10*log10(snr2(1)),10*log10(snr2(end))])
ylim([10^-4, 1])
