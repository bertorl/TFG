function [dataToSend, ofdmSymb] = transmisorTX(structTx, bits, structMod)
    % 1. Convolutional Coding
    if(structTx.convCodes == 1)
        g1 = [1 1 1 0];
        g2 = [1 0 1 1];
        g1_oct = str2num(dec2base(bin2dec(int2str(g1)),8));
        g2_oct = str2num(dec2base(bin2dec(int2str(g2)),8));
        trellis = poly2trellis([length(g1)],[g1_oct g2_oct]);
        
        [codedBits,final_state] = convenc(bits,trellis);
    else
        codedBits = bits;
    end

    if(structMod.type == 'OFDM')
        [dataToSend, ofdmSymb] = dcoofdm(codedBits,structMod.modFactor,structMod.Nportadoras,...
        structMod.mu,structMod.mapType,structMod.plotSignals,structTx.upSampling, structTx.samplesPerSymbol);
    end
end