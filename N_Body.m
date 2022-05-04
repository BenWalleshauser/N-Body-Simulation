%B.W.
%3/3

%% N-Body Simulation
clear
clc
close all

NumBodies = 20;
xdomain = 100;
ydomain = 100;

t_run = 5000;
del_t = 0.01;
t = 0:del_t:t_run;

x = zeros(NumBodies,length(t));
y = zeros(NumBodies,length(t));

m_body = 20*rand(NumBodies,1);
m_body(1) = 20*1000;
G = 1e-1;
e_soft = 0.0001;

x(:,1) = xdomain*2*(rand(NumBodies,1)-0.5);
y(:,1) = ydomain*2*(rand(NumBodies,1)-0.5);

zx = zeros(NumBodies,length(t));
zy = zeros(NumBodies,length(t));

zx(2:NumBodies,1) = 10*(rand(NumBodies-1,1)-0.5);
zy(2:NumBodies,1) = 10*(rand(NumBodies-1,1)-0.5);

F_ij_x = 0;
F_ij_y = 0;

for z = 1:length(t)
    for i = 1:NumBodies
        for j = 1:NumBodies
            if j == i
                continue
            end
            r_ij = sqrt((x(i,z)-x(j,z))^2+(y(i,z)-y(j,z))^2);
            F_ij = G*m_body(i)*m_body(j)*r_ij/(r_ij^2+e_soft^2)^(3/2);
            F_ij_x = F_ij_x + F_ij*(x(j,z)-x(i,z))/r_ij;           
            F_ij_y = F_ij_y + F_ij*(y(j,z)-y(i,z))/r_ij;  
        end
        zx(i,z+1) = zx(i,z) + del_t*(1/m_body(i))*F_ij_x;
        zy(i,z+1) = zy(i,z) + del_t*(1/m_body(i))*F_ij_y;
        
        x(i,z+1) = x(i,z) + del_t*zx(i,z);
        y(i,z+1) = y(i,z) + del_t*zy(i,z);
        
        F_ij_x = 0;
        F_ij_y = 0;
    end
end


%%

figure(1)
figure('units','normalized','outerposition',[0 0 1 1])
grid off
m_size = normalize(m_body);
m_size = m_size+2*abs(min(m_size));
m_size = 10*(m_size);
set(gca,'Color','k')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
%%
for i=1:250:length(t)
    hold on
    set(gca,'Color','k')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    xlim([-xdomain*8 xdomain*8])
    ylim([-ydomain*8 ydomain*8])
    
    j=1;
    plot(x(j,i)-x(j,i),y(j,i)-y(j,i),'o','MarkerFaceColor', 'white','MarkerSize',4*max(m_size(2:NumBodies)))
    for j = 2:NumBodies
        plot(x(j,i)-x(1,i),y(j,i)-y(1,i),'o','MarkerFaceColor', 'white','MarkerSize',m_size(j))
    end
    pause(0.1)
    clf;
end

        













