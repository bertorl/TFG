%% SIMULACIÓN COMPLETA TFG
%

clear all;
close all;
%% PAR�?METROS GENERALES

modFactor = 32;
Rb = 2^(12)*log2(modFactor);
nbits = Rb;
bits = randi([0 1], 1, nbits);
samplesPerSymbol = 4;
Nportadoras = 256;

% Fases opcionales 1-> true, 0-> false
upsampling = 1;
convCodes = 1 ;
ecualizador = 0;
prefixCode = 0;
filtroTx = 0;

%% TRANSMISOR
%Estructura con la configuración del Transmisor
structTx = struct;
structTx.upSampling = 0;
structTx.samplesPerSymbol = 4;
structTx.convCodes = 0;
structTx.filtroTx = 0;

%Estructura con la configuración de la modulación
structMod = struct;
structMod.modFactor = 32;
structMod.Nportadoras = 256;
structMod.mu = 4;
structMod.mapType = 'QAM';
structMod.plotSignals = 0;
structMod.type = 'OFDM';
structMod.prefixCode = 0;

[dataClipped, ofdmSymb] = transmisorTX(structTx, bits, structMod);  

structTx = struct;
structTx.upSampling = 0;
structTx.samplesPerSymbol = 4;
structTx.convCodes = 0;
structTx.filtroTx = 0;

%Estructura con la configuración de la modulación
structMod = struct;
structMod.modFactor = 32;
structMod.Nportadoras = 256;
structMod.mu = 4;
structMod.mapType = 'QAM';
structMod.plotSignals = 0;
structMod.type = 'OFDM';
structMod.prefixCode = 0;

[dataClipped2, ofdmSymb] = transmisorTX(structTx, bits, structMod); 
%% CANAL
[psd,f] = periodogram(dataClipped(1:length(ofdmSymb)),... 
rectwin(length(ofdmSymb)), length(ofdmSymb)*2,1, 'centered');
%hFig1 = figure('Position', figposition([46 15 30 30]));
plot(f*length(ofdmSymb)*2,10*log10(psd));
hold on
plot(f*length(ofdmSymb)*2,10*log10(psd));
grid on
xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)');
ylim([-100,20])


dataClipped = awgn(dataClipped,20);

%% RECEPTOR
% bitsReceived = [];
% recoverSymb = [];
% count = 1;
% 
% for i=1:numeroStreams %Cada iteración es un Stream del OFDM
%     %Seleccionamos un stream y lo paralelizamos
%     %Cada stream tiene 2*N+2 subportadoras
%     dataStream = dataClipped(count:i*(2*Nportadoras+2));
%     count = i*(2*Nportadoras+2) + 1;
%     recoverOfdmFrame = fft(dataStream',2*Nportadoras+2);
%     recoverSymb = [recoverSymb; recoverOfdmFrame(2:Nportadoras+1)];
% end
% 
% recoverSymb = downsample(recoverSymb, samplesPerSymbol);
% bitsStream = qamdemod(recoverSymb,modFactor,'OutputType','bit','UnitAveragePower',false);
% bitsReceived = [bitsReceived bitsStream'];
% 
% % Decodificación Convolucional
% traceback = 5* length(g1);
% 
% if(convCodes == 1)
%     decodedData = vitdec(bitsReceived,trellis,traceback,'trunc','hard');
% else
%     decodedData = bitsReceived;
% end
% [nerr, ber] = biterr(bits, decodedData);
% 
% % DownSampling
% %bitsRecovered = downsample(decodedData,samplesPerBit);















