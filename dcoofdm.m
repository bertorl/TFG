function [dataClipped, ofdmSymb] = dcoofdm(codedBits,modFactor,Nportadoras,...
                    mu,mapType,plotSignals,upSampling,samplesPerSymbol)
    %2. MODULADOR DCO-OFDM
    %%% 2.1 MAPEO M-QAM O M-PSK
    if(mapType == 'QAM')
        symbols = qammod(codedBits',modFactor,'InputType','bit','UnitAveragePower',false);
    else
        symbols = pskmod(codedBits',modFactor,'InputType','bit','UnitAveragePower',false);
    end
    
    %%% 2.2 UPSAMPLING
    if(upSampling == 1)
        %symbols = upsample(symbols, samplesPerSymbol);
        symbols = rectpulse(symbols, samplesPerSymbol);
    end
    
    %%% 2.2 OPERACIONES OFDM
    data = [];
    ofdmSymb = [];
    count = 1;
    zeroPadding = 16;
    
    numeroStreams = length(symbols)/Nportadoras;
    
%     r=rem(s,N)
%     portadorasConInfo = Nportadoras - zeroPadding;
%     rows=floor(length(symbols)/N)+1;
    
    
    for i=1:numeroStreams %Cada iteración es un Stream del OFDM
        %Este bloque realiza la parte de paralelizar los datos
        symbolsStream = symbols(count:i*Nportadoras);
        
        count = i*Nportadoras + 1;
        %Transformación hermítica
        hermitSymbols = flipud(conj(symbolsStream));
        %Se añaden los 0s correspondientes a la hermítica
        ifftInputSymbols = [0;symbolsStream;zeros(zeroPadding,1);0;zeros(zeroPadding,1);hermitSymbols];
        dataStream = ifft(ifftInputSymbols, 2*Nportadoras+2+zeroPadding*2)';
        
        ofdmSymb = [ofdmSymb ifftInputSymbols];
        data = [data dataStream];
    end

    %%% 3.3 C�?LCULO DE COMPONENTE DC
    % mu --> Factor de atenuación de DC
    dc = mu*sqrt(var(data)); % ===== RMS
    dataBiased = data + dc;
    %Quitamos todos los valores menor que 0 -> clipping
    dataClipped = 1/2*(dataBiased + abs(dataBiased));
end

