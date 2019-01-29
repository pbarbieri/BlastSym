function [F] = get_damped_armonic(f,to,fo,xi)
w = 2*pi*f;
wo = 2*pi*fo;
F = 1./(w.^2-2*1i*xi*w*wo-wo^2).*exp(-1i*w*to);
end