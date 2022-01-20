%% Parametros

clear
animation = 0;

he = 167;
hm = 100;

% Terra-Lua
Re = 6.3781366e3;
Rm = 1.7375e3;
Rs = 66182.9570725882;
D  = 3.84400e5;
wm = 2.649070880075291e-6*3600*24;

% Plutao-Caronte
% Re = 1.1883e3;
% Rm = 603.6;
% Rs = 8322.85883515962;
% D  = 19591;
% wm = 1.075432123802718e-5*3600*24;

%% Leitura do arquivo
data = readtable('trajetoria.dat');

x = data{:,1};
y = data{:,2};
t = data{:,3};

%% Uniformizacao do tempo para animacao
if animation
    npontos = 1000;
    x = interp1(t,x,linspace(0,t(end),npontos));
    y = interp1(t,y,linspace(0,t(end),npontos));
    t = linspace(0,t(end),npontos);
else
    npontos = length(x);
end

%% Animacao
figure
hold on

R_contact = Rs;
if ~animation
    Rs = 0;
end

rectangle('position',[-Re -Re 2*Re 2*Re],'curvature',[1 1],'facecolor','b','edgecolor','b')
rectangle('position',[-D -D 2*D 2*D],'curvature',[1 1],'edgecolor','b')
rectangle('position',[-(Re+he) -(Re+he) 2*(Re+he) 2*(Re+he)],'curvature',[1 1],'edgecolor','k')
soi = rectangle('position',[D*cos(wm*t(1))-Rs D*sin(wm*t(1))-Rs 2*Rs 2*Rs],'curvature',[1 1],'edgecolor','r');
lmo = rectangle('position',[D*cos(wm*t(1))-(Rm+hm) D*sin(wm*t(1))-(Rm+hm) 2*(Rm+hm) 2*(Rm+hm)],'curvature',[1 1],'edgecolor','k');
traj = plot(x(1),y(1),'k');
moon = rectangle('position',[D*cos(wm*t(1))-Rm D*sin(wm*t(1))-Rm 2*Rm 2*Rm],'curvature',[1 1],'facecolor','r','edgecolor','r');

% ax = gca;
% ax.XAxis.Exponent = 4;
% ax.YAxis.Exponent = 4;
% 
% xticks([0 0.25 0.5 0.75 1 1.25 1.5 1.75 2]*1e4)
% yticks([0 0.25 0.5 0.75 1 1.25 1.5 1.75 2]*1e4)

xlabel('km')
ylabel('km')
if animation
    title(['t = ' num2str(t(1)) ' dias'])
end
axis equal
grid on
box on
xlim([-0.1*D 1.05*(D+Rs)])
ylim([-0.1*D 1.05*(D+Rs)])

for i=2:npontos
    set(moon,'position',[D*cos(wm*t(i))-Rm D*sin(wm*t(i))-Rm 2*Rm 2*Rm])
    set(soi,'position',[D*cos(wm*t(i))-Rs D*sin(wm*t(i))-Rs 2*Rs 2*Rs])
    set(lmo,'position',[D*cos(wm*t(i))-(Rm+hm) D*sin(wm*t(i))-(Rm+hm) 2*(Rm+hm) 2*(Rm+hm)])
    set(traj,'xdata',x(1:i),'ydata',y(1:i))
    
    if abs((x(i)-D*cos(wm*t(i)))^2+(y(i)-D*sin(wm*t(i)))^2-R_contact^2) < R_contact^2/500
        plot(x(i),y(i),'ko','markerfacecolor','k','markersize',3)
    end
    
    if animation
        title(['t = ' num2str(t(i)) ' dias'])
        drawnow
    end
end

hold off