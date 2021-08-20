function LEiDA_stats

%%%%%%%
%  
% Code to analyze the LEiDA clustering results 
% Comparing the occupancy of FC states between 2 Group
% Pre- vs Post-
% Can be adapted for more Group
%
% Modified by João Sanguessuga, April 2021
%
%%%%%%%

load_file = 'LEiDA_EigenVectors'; % generated with LEiDA_data_X
kmeans_file  = 'LEiDA_Kmeans_results'; % generated with LEiDA_cluster_X
save_file='LEiDA_Stats.mat';
number_permutations =10000; % to debug you can set a lower number, like 500
% but for final results choose around 10000 permutations to ensure stable
% results.
TR=2; % Repetition Time (seconds)

load(kmeans_file,'Kmeans_results','rangeK')
%load(kmeans_file,'Kmeans_placebo','rangeK')
%Kmeans_results=Kmeans_placebo;
%clear Kmeans_placebo

load(load_file,'Time_sessions','total_scans','Group')

%% Calculate the Occupancy of each state
       
% For every fMRI subject and group calculate probability of each state c.
P=zeros(total_scans,length(rangeK),rangeK(end));
LT=zeros(total_scans,length(rangeK),rangeK(end));

for k=1:length(rangeK) 
        for s=1:total_scans
            
            % Select the time points representing each scan
            T=Time_sessions==s;
            Ctime=Kmeans_results{k}.IDX(T);
            
            for c=1:rangeK(k)
                % Probability Time_sessions
                P(s,k,c)=mean(Ctime==c);
                
                % Mean Lifetime
                Ctime_bin=Ctime==c;
                    
                % Detect switches in and out of this state
                a=find(diff(Ctime_bin)==1);
                b=find(diff(Ctime_bin)==-1);
                    
                % We discard the cases where state sarts or ends ON
                if length(b)>length(a)
                    b(1)=[];
                elseif length(a)>length(b)
                    a(end)=[];
                elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
                    b(1)=[];
                    a(end)=[];
                end
                 
                if ~isempty(a) && ~isempty(b)
                    C_Durations=b-a;
                else
                    C_Durations=0;
                end
                LT(s,k,c)=mean(C_Durations)*TR;
            end
        end
end

%%  PERMUTATION STATISTICS
 
disp('Testing differences in Occupancy between Groups')
P_pval=zeros(length(rangeK),rangeK(end));
LT_pval=zeros(length(rangeK),rangeK(end));

Index_Poor=find(Group==0); % Number of subjects in the group 0 (poor performance)
Index_Good=find(Group==1); % Number of subjects in the group 1 (good performance)

    for k=1:length(rangeK)
        disp(['Now running statistics for ' num2str(rangeK(k)) ' FC states'])
        for c=1:rangeK(k)
            % Compare Probabilities
            a=squeeze(P(Index_Poor,k,c))';  % Vector containing Prob of c in g1 (good performance)
            b=squeeze(P(Index_Good,k,c))';  % Vector containing Prob of c in g2 (poor performance)
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],number_permutations,0.05,'ttest');
            P_pval(k,c)=min(stats.pvals);
            
            % Compare Lifetimes
            a=squeeze(LT(Index_Poor,k,c))'; % Vector containing Prob of c in g1 (good performance)
            b=squeeze(LT(Index_Good,k,c))'; % Vector containing Prob of c in g2 (poor performance)
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],number_permutations,0.05,'ttest');
            LT_pval(k,c)=min(stats.pvals);      
        end
    end
    
save(save_file,'a','b','P','P_pval','LT','LT_pval','rangeK','Kmeans_results','Index_Poor','Index_Good','kmeans_file','load_file')
   
disp('%%%%% STATISTICS SUCCESSFULLY COMPLETED %%%%%%%')


disp('Saving LEiDA results')

%% PLOT P-VALUES OVER THE RANGE OF K

figure
subplot(1,2,1)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK) 
    for c=1:rangeK(k)        
        if P_pval(k,c)>0.05
            semilogy(rangeK(k),P_pval(k,c),'*k');
            text(rangeK(k),P_pval(k,c),[' ' num2str(c)])
        end
        if P_pval(k,c)<0.05 && P_pval(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),P_pval(k,c),'*r');
            text(rangeK(k),P_pval(k,c),[' ' num2str(c)])
        end
        if P_pval(k,c)<(0.05/rangeK(k)) && P_pval(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),P_pval(k,c),'*g');
            text(rangeK(k),P_pval(k,c),[' ' num2str(c)])
        end
        if P_pval(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),P_pval(k,c),'*b');
            text(rangeK(k),P_pval(k,c),[' ' num2str(c)])
        end
    end
end

title('Probability Poor vs Good')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off

figure
subplot(1,2,1)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK) 
    for c=1:rangeK(k)        
        if LT_pval(k,c)>0.05
            semilogy(rangeK(k),LT_pval(k,c),'*k');
            text(rangeK(k),LT_pval(k,c),[' ' num2str(c)])
        end
        if LT_pval(k,c)<0.05 && LT_pval(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),LT_pval(k,c),'*r');
            text(rangeK(k),LT_pval(k,c),[' ' num2str(c)])
        end
        if LT_pval(k,c)<(0.05/rangeK(k)) && LT_pval(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),LT_pval(k,c),'*g');
            text(rangeK(k),LT_pval(k,c),[' ' num2str(c)])
        end
        if LT_pval(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),LT_pval(k,c),'*b');
            text(rangeK(k),LT_pval(k,c),[' ' num2str(c)])
        end
    end
end

title('Dweimes Poor vs Good')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off

%% Plot FC patterns and stastistics between groups

disp(' ')
disp('%%% PLOTS %%%%')
disp(['Choose number of clusters between ' num2str(rangeK(1)) ' and ' num2str(rangeK(end)) ])
Pmin_pval=min(P_pval(P_pval>0));
LTmin_pval=min(LT_pval(LT_pval>0));
if Pmin_pval<LTmin_pval
   [k,~]=ind2sub([length(rangeK),max(rangeK)],find(P_pval==Pmin_pval));
else
   [k,~]=ind2sub([length(rangeK),max(rangeK)],find(LT_pval==LTmin_pval));
end
disp(['Note: The most significant difference is detected with K=' num2str(rangeK(k)) ' (p=' num2str(min(Pmin_pval,LTmin_pval)) ')'])  

% To correct for multiple comparisons, you can divide p by the number of
% clusters considered

K = input('Number of clusters: ');
Best_Clusters=Kmeans_results{rangeK==K};
k=find(rangeK==K);

% Clusters are sorted according to their probability of occurrence
ProbC=zeros(1,K);
for c=1:K
    ProbC(c)=mean(Best_Clusters.IDX==c);
end
[~, ind_sort]=sort(ProbC,'descend'); 

% Get the K patterns
V=Best_Clusters.C(ind_sort,:);
[~, N]=size(Best_Clusters.C);
Order=[1:2:N N:-2:2];

figure
colormap(jet) 
% Pannel A - Plot the FC patterns over the cortex 
% Pannel B - Plot the FC patterns in matrix format
% Pannel C - Plot the probability of each state in each condition
% Pannel D - Plot the lifetimes of each state in each condition
   
for c=1:K
    subplot(4,K,c)
    plot_nodes_in_cortex(V(c,:))
    title({['PL State ' num2str(c), newline, 'V_C_' num2str(c) ]})
    subplot(4,K,K+c)
    FC_V=V(c,:)'*V(c,:);  
    li=max(abs(FC_V(:)));
    imagesc(FC_V(Order,Order),[-li li])   
    axis square
    title(['V_C_' num2str(c) '.V_C_' num2str(c) '^T'])
    ylabel('Brain area #')
    xlabel('Brain area #') 
    
    subplot(4,K,2*K+c)  
            Good=squeeze(P(Index_Good,k,ind_sort(c)))';
            Poor=squeeze(P(Index_Poor,k,ind_sort(c)))';
            bar([mean(Good) mean(Poor)],'EdgeColor','w','FaceColor',[.5 .5 .5])
            hold on
            % Error bar containing the standard error of the mean
            errorbar([mean(Good) mean(Poor)],[std(Good)/sqrt(numel(Good)) std(Poor)/sqrt(numel(Poor))],'LineStyle','none','Color','k')
            set(gca,'XTickLabel',{'Good', 'Poor'})
            if P_pval(k,ind_sort(c))<0.05 && P_pval(k,ind_sort(c))>(0.05/K)
               plot(1.5,max([mean(Good) mean(Poor)])+.01,'*r')
            end
            if P_pval(k,ind_sort(c))<(0.05/K)
               plot(1.5,max([mean(Good) mean(Poor)])+.01,'*g')
            end

            if c==1
                ylabel('Probability')
            end
            box off
            
     subplot(4,K,3*K+c)  
            Good=squeeze(LT(Index_Good,k,ind_sort(c)))';
            Poor=squeeze(LT(Index_Poor,k,ind_sort(c)))';
            bar([mean(Good) mean(Poor)],'EdgeColor','w','FaceColor',[.5 .5 .5])
            hold on
            errorbar([mean(Good) mean(Poor)],[std(Good)/sqrt(numel(Good)) std(Poor)/sqrt(numel(Poor))],'LineStyle','none','Color','k')
            set(gca,'XTickLabel',{'Good', 'Poor'})
            if LT_pval(k,ind_sort(c))<0.05 && LT_pval(k,ind_sort(c))>(0.05/K)
               plot(1.5,max([mean(Good) mean(Poor)])+.01,'*r')
            end
            if LT_pval(k,ind_sort(c))<(0.05/K)
               plot(1.5,max([mean(Good) mean(Poor)])+.01,'*g')
            end

            if c==1
                ylabel('Lifetime (seconds)')
            end
            box off           
end
    