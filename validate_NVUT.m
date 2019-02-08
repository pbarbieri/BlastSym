function [] = validate_NVUT(ac)
t = linspace(0,1,10001).';
AT = sin(2*pi*t);
AT1 = AT-ac; AT1(AT<ac)=0;
[VT1,UT1] = Get_VUT(AT1,t);
[VT2,UT2,t2] = Get_NVUT(AT1,t,ac);

close all
figure(1)
plot(t,AT);
hold on
plot(t,AT1);
plot([0,1],[ac ac],'--k');
grid on
xlabel('t[s]');
ylabel('A[m/s2]');
legend({'AT','AT>ky','ky'},'location','SouthWest');
set(gcf,'Color',[1 1 1],'Position',[50,50,800,300]);
grid on

figure(2)
plot(t,VT1);
hold on
plot(t2,VT2);
grid on
set(gcf,'Color',[1 1 1],'Position',[50,50,800,300])
xlabel('t[s]');
ylabel('V[m/s]');
set(gcf,'Color',[1 1 1],'Position',[100,100,800,300])
legend({'Regular Newmark integration','Integration with friction braking'},'location','SouthEast');



figure(3)
plot(t,UT1);
hold on
plot(t2,UT2);
grid on
set(gcf,'Color',[1 1 1],'Position',[50,50,800,300])
xlabel('t[s]');
ylabel('U [m]');
set(gcf,'Color',[1 1 1],'Position',[150,150,800,300])
legend({'Regular Newmark integration','Integration with friction braking'},'location','SouthEast');
end

function [VT,UT,t] = Get_NVUT(AT,t,ac)
% NEWMARK INTEGRATION OF ACELERATION TIMESERIE(s) FOR SLINDING BLOCK WITH
% CRICITAL ACCELERATION
%   AT: Acceleration timeseries of record(s) by column.
%   t:  Time vector by column.
%   VT: Velocity timeseries of record(s) by column.
%   UT: Displacement timeseries of record(s) by column.
%   ac: Critical acceleration

% Default Newmark integration Parameters
gamma = 1/2;
beta = 1/4;

AT = abs(AT);

NP = size(AT,1);
dt = t(2) - t(1);
VT = zeros(NP,size(AT,2));
UT = VT;
for j = 2:NP
    if AT(j)>0
        VT(j,:) = VT(j-1,:)+(1-gamma)*dt*AT(j-1,:)+gamma*dt*AT(j,:);
        UT(j,:) = UT(j-1,:)+dt*VT(j-1,:)+(1/2-beta)*dt^2*AT(j-1,:)+beta*dt^2*AT(j,:);
    else
        VT(j,:) = VT(j-1,:)+(1-gamma)*dt*(-ac)+gamma*dt*(-ac);
        if VT(j,:)<0, VT(j,:) = 0;  end
        UT(j,:) = UT(j-1,:)+dt*VT(j-1,:)+(1/2-beta)*dt^2*(-ac)+beta*dt^2*(-ac);
    end
    
end
j = NP;
while VT(j,:)>0
    j = j+1;
    VT(j,:) = VT(j-1,:)+(1-gamma)*dt*(-ac)+gamma*dt*(-ac);
    if VT(j,:)<0, VT(j,:) = 0;  end
    UT(j,:) = UT(j-1,:)+dt*VT(j-1,:)+(1/2-beta)*dt^2*(-ac)+beta*dt^2*(-ac);
    t(j) = t(j-1) + dt;
end

end