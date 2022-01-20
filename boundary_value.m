clear

pdf = 1;

horario = readtable('horario.dat');
antihorario = readtable('antihorario.dat');

[~,horario_min_id] = min(horario{:,3});
[~,antihorario_min_id] = min(antihorario{:,3});

%% lambda1 x deltav
figure
hold on

plot(horario{:,1},horario{:,3},'b','linewidth',3)
plot(antihorario{:,1},antihorario{:,3},'r','linewidth',3)

xline(horario{horario_min_id,1},'b--','\Deltav mínimo','linewidth',1.5,'LabelHorizontalAlignment','left','LabelVerticalAlignment','top','FontSize',12)
xline(antihorario{antihorario_min_id,1},'r--','\Deltav mínimo','linewidth',1.5,'LabelHorizontalAlignment','right','LabelVerticalAlignment','top','FontSize',12)

xlabel('\lambda_1 [°]')
ylabel('\Deltav [km/s]')
hold off
grid on
grid minor
box on
legend('Chegada horária','Chegada anti-horária')
set(gca,'FontSize',12)
if pdf
    exportgraphics(gca,'deltav.pdf')
end

%% lambda1 x deltav1

[~,horario_min_deltav1_id] = min(horario{:,4});
[~,antihorario_min_deltav1_id] = min(antihorario{:,4});

figure
hold on

plot(horario{:,1},horario{:,4},'b','linewidth',3)
plot(antihorario{:,1},antihorario{:,4},'r','linewidth',3)

xline(horario{horario_min_id,1},'b--','\Deltav mínimo','linewidth',1.5,'LabelHorizontalAlignment','right','LabelVerticalAlignment','middle','FontSize',12)
xline(antihorario{antihorario_min_id,1},'r--','\Deltav mínimo','linewidth',1.5,'LabelHorizontalAlignment','right','LabelVerticalAlignment','middle','FontSize',12)

xline(horario{horario_min_deltav1_id,1},'b-.','\Deltav_1 mínimo','linewidth',1.5,'LabelHorizontalAlignment','left','LabelVerticalAlignment','middle','FontSize',12)
% xline(antihorario{antihorario_min_deltav1_id,1},'r-.','\Deltav_1 mínimo','linewidth',1.5,'LabelHorizontalAlignment','left','LabelVerticalAlignment','middle','FontSize',12)

xlabel('\lambda_1 [°]')
ylabel('\Deltav_1 [km/s]')
hold off
grid on
grid minor
box on
legend('Chegada horária','Chegada anti-horária')
set(gca,'FontSize',12)
if pdf
    exportgraphics(gca,'deltav1.pdf')
end

%% lambda1 x deltav2

[~,horario_min_deltav2_id] = min(horario{:,5});
[~,antihorario_min_deltav2_id] = min(antihorario{:,5});

figure
hold on

plot(horario{:,1},horario{:,5},'b','linewidth',3)
plot(antihorario{:,1},antihorario{:,5},'r','linewidth',3)

xline(horario{horario_min_id,1},'b--','\Deltav mínimo','linewidth',1.5,'LabelHorizontalAlignment','left','LabelVerticalAlignment','middle','FontSize',12)
xline(antihorario{antihorario_min_id,1},'r--','\Deltav mínimo','linewidth',1.5,'LabelHorizontalAlignment','right','LabelVerticalAlignment','middle','FontSize',12)

xline(horario{horario_min_deltav2_id,1},'b-.','\Deltav_2 mínimo','linewidth',1.5,'LabelHorizontalAlignment','right','LabelVerticalAlignment','middle','FontSize',12)
xline(antihorario{antihorario_min_deltav2_id,1},'r-.','\Deltav_2 mínimo','linewidth',1.5,'LabelHorizontalAlignment','right','LabelVerticalAlignment','middle','FontSize',12)

xlabel('\lambda_1 [°]')
ylabel('\Deltav_2 [km/s]')
hold off
grid on
grid minor
box on
legend('Chegada horária','Chegada anti-horária')
set(gca,'FontSize',12)
if pdf
    exportgraphics(gca,'deltav2.pdf')
end

%% lambda1 x deltat
figure
hold on

plot(horario{:,1},horario{:,6},'b','linewidth',3)
plot(antihorario{:,1},antihorario{:,6},'r','linewidth',3)

xline(horario{horario_min_id,1},'b--','\Deltav mínimo','linewidth',1.5,'LabelHorizontalAlignment','left','LabelVerticalAlignment','top','FontSize',12)
xline(antihorario{antihorario_min_id,1},'r--','\Deltav mínimo','linewidth',1.5,'LabelHorizontalAlignment','right','LabelVerticalAlignment','top','FontSize',12)

xlabel('\lambda_1 [°]')
ylabel('\Deltat [dias]')
hold off
grid on
grid minor
box on
legend('Chegada horária','Chegada anti-horária','Location','southeast')
set(gca,'FontSize',12)
if pdf
    exportgraphics(gca,'deltat.pdf')
end

%% lambda1 x thetaep0
figure
hold on

plot(horario{:,1},horario{:,7},'b','linewidth',3)
plot(antihorario{:,1},antihorario{:,7},'r','linewidth',3)

xline(horario{horario_min_id,1},'b--','\Deltav mínimo','linewidth',1.5,'LabelHorizontalAlignment','left','LabelVerticalAlignment','top','FontSize',12)
xline(antihorario{antihorario_min_id,1},'r--','\Deltav mínimo','linewidth',1.5,'LabelHorizontalAlignment','right','LabelVerticalAlignment','top','FontSize',12)

xlabel('\lambda_1 [°]')
ylabel('\theta_P_p_0 [°]')
hold off
grid on
grid minor
box on
legend('Chegada horária','Chegada anti-horária','Location','southeast')
set(gca,'FontSize',12)
if pdf
    exportgraphics(gca,'thetapp0.pdf')
end

%% lambda1 x energym

[~,horario_min_energy_id] = min(horario{:,8});
[~,antihorario_min_energy_id] = min(antihorario{:,8});

figure
hold on

plot(horario{:,1},horario{:,8},'b','linewidth',3)
plot(antihorario{:,1},antihorario{:,8},'r','linewidth',3)

xline(horario{horario_min_id,1},'b--','\Deltav mínimo','linewidth',1.5,'LabelHorizontalAlignment','left','LabelVerticalAlignment','middle','FontSize',12)
xline(antihorario{antihorario_min_id,1},'r--','\Deltav mínimo','linewidth',1.5,'LabelHorizontalAlignment','right','LabelVerticalAlignment','middle','FontSize',12)

xline(horario{horario_min_energy_id,1},'b-.','\epsilon_c mínimo','linewidth',1.5,'LabelHorizontalAlignment','right','LabelVerticalAlignment','middle','FontSize',12)
xline(antihorario{antihorario_min_energy_id,1},'r-.','\epsilon_c mínimo','linewidth',1.5,'LabelHorizontalAlignment','right','LabelVerticalAlignment','middle','FontSize',12)

xlabel('\lambda_1 [°]')
ylabel('\epsilon_C [km^2/s^2]')
hold off
grid on
grid minor
box on
legend('Chegada horária','Chegada anti-horária')
set(gca,'FontSize',12)
if pdf
    exportgraphics(gca,'energyc.pdf')
end

%% lambda1 x em
figure
hold on

plot(horario{:,1},horario{:,9},'b','linewidth',3)
plot(antihorario{:,1},antihorario{:,9},'r','linewidth',3)

xline(horario{horario_min_id,1},'b--','\Deltav mínimo','linewidth',1.5,'LabelHorizontalAlignment','left','LabelVerticalAlignment','top','FontSize',12)
xline(antihorario{antihorario_min_id,1},'r--','\Deltav mínimo','linewidth',1.5,'LabelHorizontalAlignment','right','LabelVerticalAlignment','top','FontSize',12)

xlabel('\lambda_1 [°]')
ylabel('e_C')
hold off
grid on
grid minor
box on
legend('Chegada horária','Chegada anti-horária')
set(gca,'FontSize',12)
if pdf
    exportgraphics(gca,'ec.pdf')
end