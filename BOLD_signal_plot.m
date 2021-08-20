%%%%%%%
%
% code to plot the variation of a subject's BOLD signal over time.
% 
% Written by João Sanguessuga, May 2021
%
%%%%%%%

N_areas=90; % Number of brain areas to consider from parcellation
load Aging_data.mat Subjects 
x=1; % Representative subject of the patient group (0)
y=98; % Representative subject of the control group (1)

disp(Subjects(x))            
% load BOLD signals in each scan
signal_patients_1=struct2array(load('Aging_data.mat', Subjects{x}));
signal_patients_2=signal_patients_1(:,1:N_areas)';

disp(Subjects(y))            
% load BOLD signals in each scan
signal_controls_1=struct2array(load('Aging_data.mat', Subjects{y}));
signal_controls_2=signal_controls_1(:,1:N_areas)';

figure
colormap(jet)
subplot(2,4,1:3)
plot(signal_patients_2')
xlabel('Time (TR=2s)')
ylabel('BOLD Signal')
ylim([-3 3])
xlim([0 175])
subplot(2,4,4)
Order=[1:2:N_areas N_areas:-2:2];
FC_patients=corrcoef(signal_patients_2');
imagesc(FC_patients(Order,Order),[-1 1])
title('Matrix FC Static')
axis square

figure
colormap(jet)
subplot(2,4,1:3)
plot(signal_controls_2')
xlabel('Time (TR=2s)')
ylabel('BOLD Signal')
ylim([-5 5])
xlim([0 175])
subplot(2,4,4)
FC_controls=corrcoef(signal_controls_2');
imagesc(FC_controls(Order,Order),[-1 1])
title('Matrix FC Static')
axis square



