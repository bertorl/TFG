clear all;

N = 64;
  over_sample_factor = 2;
  M = N*over_sample_factor;
  Mod = 16;
  symbol = 1;
  bitlength = N*log2(Mod)*symbol;
  itr_num = 1000;
  fft_len = 2*M;
  signal_freq = zeros(itr_num,fft_len);
  for itr = 1:itr_num
      bit_data = randi([0,1],bitlength,1);
      zp_before = qammod(bit_data,Mod,'InputType','bit','UnitAveragePower',false);
      after_zp = zeros(1,M);
      after_zp(1:N/2) = zp_before(N/2+1:N);
      after_zp(M-N/2+1:M) = zp_before(1:N/2);
      ofdm_symbol = ifft(after_zp);
      signal_freq(itr,:) = abs(fft(ofdm_symbol,fft_len)).^2; 
  end
  PSD_mean = mean(signal_freq,1);
  plot(fftshift(10*log10(PSD_mean)));