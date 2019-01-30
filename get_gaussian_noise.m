function [AF] = get_gaussian_noise(f,to,fo,sigma)
% Normal noise with complex mean and deviation
mu = fo;
% sigma = 2*pi*xi*fo;
%AF = exp(-(f-mu).^2/(2*sigma^2))/sigma.*exp(-1i*2*pi*f*to);
AF = exp(-(f-mu).^2/(2*sigma^2))/sigma.*exp(-1i*2*pi*f*to);
end