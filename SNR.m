function y = SNR(x0,x)
y = 20*log10(norm(x0(:))/norm(x(:)-x0(:)));
% x0: data after filtering
% x noise free data
% snr = SNR(d2,d1)