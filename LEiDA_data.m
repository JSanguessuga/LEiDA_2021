function LEiDA_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS
%
%
% - Reads the BOLD data from the folders
% - Computes the BOLD phases
% - Calculates the instantaneous BOLD synchronization matrix
% - Saves the instantaneous Leading Eigenvectors
%
% Saves the Leading_Eigenvectors and FCD matrices into LEiDA_data.mat
%
% Written by Joana Cabral, May 2016 (joanacabral@med.uminho.pt)
% Modified by Joana Cabral, November 2020
% Modified by João Sanguessuga, April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USER: ADJUST THIS SECTION TO YOUR DATASET
N_areas=90; % Number of brain areas to consider from parcellation
Tmax=175; % Number of TRs per scan (if different in each scan, set maximum)

save_file='LEiDA_EigenVectors';

%%
% Create empty variables to save patterns 
% Define total number of scans that will be read
load Aging_data.mat Group
load Aging_data.mat Subjects 
total_scans=length(Subjects);

V1_all=zeros(total_scans*(Tmax-2),N_areas); % All leading eigenvectors
Time_sessions=zeros(1,total_scans*(Tmax-2)); % vector with subject nr and group at each t
t_all=0; % Index of time (starts at 0 and will be updated until total_scans*(Tmax-2))

for s=1:total_scans
       
            disp(Subjects(s))
            
            % load BOLD signals in each scan
            signal=struct2array(load('Aging_data.mat', Subjects{s}));
            signal=signal(:,1:N_areas)';
            %figure
            %plot(signal)

            % Get the BOLD phase using the Hilbert transform
            for seed=1:N_areas
                signal(seed,:)=angle(hilbert(signal(seed,:)));     
            end

            for t=2:size(signal,2)-1

            [v1,~]=eigs(cos(signal(:,t)-signal(:,t)'),1);

                if sum(v1)>0
                   v1=-v1;
                end

            t_all=t_all+1; %Update time
            V1_all(t_all,:)=v1; 
            Time_sessions(t_all)=s;
            end
end

% Reduce size in case some scans have less TRs than Tmax
V1_all(t_all+1:end,:)=[];
Time_sessions(:,t_all+1:end)=[];

save(save_file,'V1_all','Time_sessions','total_scans','Group')

disp(['BOLD Phase Eigenvectors saved successfully as ' save_file])

