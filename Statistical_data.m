%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to perform statistical data evaluation using ttest2
% The two-sample t-test (also known as the independent samples t-test) 
% is a method used to test whether the unknown population means of two groups are equal or not.
%
% Written by João Sanguessuga, April 2021
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
% Import data from excel to matlab
G_pacientes_age=xlsread('Dados_Limpos.xlsx',1,'AG2:AG44'); %Ages of Group 1 (poor cognitive performance)
G_controlo_age=xlsread('Dados_Limpos.xlsx',1,'AG45:AG99'); %Ages of Group 2 (good cognitive performance)

[hypothesis_age,p_valuet_age,ci_age,tstats_age] = ttest2(G_pacientes_age,G_controlo_age,'Vartype','unequal','Alpha',0.05);

disp('Age');
if p_valuet_age>0.05
    fprintf('p_valuet = %f | hypothesis = %d\n',p_valuet_age,hypothesis_age)
    fprintf('The data are normally distributed')
    fprintf('\n')

else
    fprintf('p_valuet = %f | hypothesis = %d\n',p_valuet_age,hypothesis_age)
    fprintf('The data are NOT normally distributed \n\n')
end

% Create a table of group statistics
Groups={'Grupo_pacientes (0)';'Grupo_controlo (1)'};
N_age=[length(G_pacientes_age);length(G_controlo_age)];
Mean_age=[mean(G_pacientes_age);mean(G_controlo_age)];
Median_age=[median(G_pacientes_age);median(G_controlo_age)];
Range_age=[range(G_pacientes_age);range(G_controlo_age)];
Std_age=[std(G_pacientes_age);std(G_controlo_age)];
Min_age=[min(G_pacientes_age);min(G_controlo_age)];
Max_age=[max(G_pacientes_age);max(G_controlo_age)];
Descriptive_Stats_age=table(N_age,Mean_age,Median_age,Range_age,Std_age,Min_age,Max_age,'RowNames',Groups);
fprintf('\n')

