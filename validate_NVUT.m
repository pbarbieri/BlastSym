function [] = validate_NVUT(t,AT,ac)
% t = linspace(0,1,10001).';
% AT = sin(2*pi*t);


[ATblock,VT,UT,t2] = Get_NVUT(AT,t,ac);

close all
figure(1)
plot(t,AT);
hold on
plot(t,ATblock);
plot([t(1),t(end)],[ac ac],'--k');
grid on
xlabel('t[s]');
ylabel('A[m/s2]');
legend({'AT','AT>ky','ky'},'location','SouthWest');
set(gcf,'Color',[1 1 1],'Position',[50,50,800,300]);
grid on

figure(2)
hold on
plot(t2,VT);
grid on
set(gcf,'Color',[1 1 1],'Position',[50,50,800,300])
xlabel('t[s]');
ylabel('V[m/s]');
set(gcf,'Color',[1 1 1],'Position',[100,100,800,300])



figure(3)
hold on
plot(t2,UT);
grid on
set(gcf,'Color',[1 1 1],'Position',[50,50,800,300])
xlabel('t[s]');
ylabel('U [m]');
set(gcf,'Color',[1 1 1],'Position',[150,150,800,300])
end

function [ATblock,VT,UT,t] = Get_NVUT(AT,t,ac)
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

NP = size(AT,1);
dt = t(2) - t(1);
VT = zeros(NP,1);
UT = zeros(NP,1);
ATblock = zeros(NP,1);

for j = 2:NP
    ATblock(j) = AT(j)-ac;
    VT(j) = VT(j-1)+(1-gamma)*dt*ATblock(j-1)+gamma*dt*ATblock(j);
    if VT(j)<0
        VT(j) = 0; 
        ATblock(j) = 0;
    end
    UT(j) = UT(j-1)+dt*VT(j-1)+(1/2-beta)*dt^2*ATblock(j-1)+beta*dt^2*ATblock(j);
    if UT(j)<UT(j-1), UT(j)=UT(j-1); end
end
j= NP;
while VT(j)>0
    j = j+1;
    VT(j) = VT(j-1)+(1-gamma)*dt*(-ac)+gamma*dt*(-ac);
    if VT(j)<0, VT(j) = 0;  end
    UT(j) = UT(j-1)+dt*VT(j-1)+(1/2-beta)*dt^2*(-ac)+beta*dt^2*(-ac);
    t(j) = t(j-1) + dt;
end


% for j = 2:NP
%     ATblock(j) = AT(j)-ac;
%     VT(j) = VT(j-1)+0.5*dt*(ATblock(j-1)+ATblock(j));
%     if VT(j)<0
%         VT(j) = 0;
%         ATblock(j) = 0;
%     end
%     UT(j) = UT(j-1)+dt*VT(j-1)+(2*ATblock(j-1)+ATblock(j))*dt^2/6;
%     if UT(j)<UT(j-1), UT(j)=UT(j-1); end
% end
% 
% j= NP;
% while VT(j)>0
%     j = j+1;
%     VT(j) = VT(j-1)-dt*ac;
%     if VT(j)<0, VT(j) = 0;  end
%     UT(j) = UT(j-1)+dt*VT(j-1)-(2*ac+ac)*dt^2/6;
%     t(j) = t(j-1) + dt;
% end

end