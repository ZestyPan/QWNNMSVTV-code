function x = SNR(ref, sig)

% snr -- Compute Signal-to-Noise Ratio for images
%
% Usage:
%       x = SNR(ref, sig)
%
% Input:
%       ref         Reference image ��ʵͼ��
%       sig         Modified image ������
%  
% Output:
%       x           SNR value

mse = mean((ref(:)-sig(:)).^2);
dv = var(ref(:),1);
x = 10*log10(dv/mse);
