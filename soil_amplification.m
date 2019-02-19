function [AF,y] = soil_amplification(f,h,vs,shape,failure)

if size(f,2)>size(f,1), f =f.'; end
% bessel J0 zeos
a = [2.405,5.525,8.645,11.792];
% Natural period of the first 4 modes
switch lower(shape)
    case {'layer','l'}
        Tn = 4*[1;1/3;1/5;1/7];
    case {'dam','d'}
        Tn = 2.61*[1;a(1)/a(2);a(2)/a(3);a(3)/a(4)];
    case {'slope','s'}
        Tn = [0.462;0.173;0.107;0.077];
    otherwise
        error('Unknown shape. Try: layer, dam or slope.')
end
Tn = Tn*h/vs;
% Damping of the first four modes
xi = [0.2;0.5;0.5;0.5];
% SDOF amplification function
s(f,xi,T) = sqrt((1+4*(f*xi*T).^2)./((1-(f*T).^2)+4*(f*xi*T).^2));
Sa = [s(f,xi(1),Tn(1));s(f,xi(2),Tn(2));s(f,xi(3),Tn(3));s(f,xi(4),Tn(4))];

% Shape function
y = linspace(0,h,10^3+1).';
dy = y(2)-y(1);
switch lower(failure)
    case {'wedge','w'}
        I = @(a) 2./y.^2.*cumtrapz(2*besselj(0,y/h*a)/(a*besselj(1,a)).*y)*dy;
    case {'slip','s'}
        I = @(a) 1./y.*cumtrapz(2*besselj(0,y/h*a)/(a*besselj(1,a)))*dy;
    otherwise
        error('Unknown failure shape. Try: wedge or slip.')
end
Kn = [I(a(1));I(a(2));I(a(3));1-I(a(1))-I(a(2))-I(a(3))];

% Amplification factor
end