%%%%%%%%
%
% Compare Psichological and Physiological markers with probability and 
% dwell times of LEiDA networks
%
%
% Modified by João Sanguessuga, May 2021
%
%%%%%%%%

load Scores.mat School_years
stats_file = 'LEiDA_stats';
load(stats_file,'rangeK','P','Index_Good','Index_Poor')
kmeans_file  = 'LEiDA_Kmeans_results';
load(kmeans_file,'Kmeans_results')
load AAL_labels.mat label90


%% Patient scores
Patient_Scores=School_years(1:43,:);

disp('School_years: Patients')
 
for k=1:length(rangeK)
        
    for c=1:rangeK(end)
            
            P_Patients=squeeze(P(Index_Poor,k,c));
            [cc_Patients,p_Patients]=corrcoef(P_Patients,Patient_Scores);
            
            Patients_P_pval(k,c)=p_Patients(2);
            
         if Patients_P_pval(k,c)<(0.05)
            disp(['signif correlation with prob K=' num2str(rangeK(k)) ' c=' num2str(c)])
            disp(['cc=' num2str(cc_Patients(2)) ' p=' num2str(p_Patients(2))])
         end   
            
    end
end

%% Plot Patient scores 

figure
subplot(1,2,1)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK) 
    for c=1:rangeK(k)        
        if Patients_P_pval(k,c)>0.05
            semilogy(rangeK(k),Patients_P_pval(k,c),'*k');
            text(rangeK(k),Patients_P_pval(k,c),[' ' num2str(c)])
        end
        if Patients_P_pval(k,c)<0.05 && Patients_P_pval(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),Patients_P_pval(k,c),'*r');
            text(rangeK(k),Patients_P_pval(k,c),[' ' num2str(c)])
        end
        if Patients_P_pval(k,c)<(0.05/rangeK(k)) && Patients_P_pval(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),Patients_P_pval(k,c),'*g');
            text(rangeK(k),Patients_P_pval(k,c),[' ' num2str(c)])
        end
        if Patients_P_pval(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),Patients_P_pval(k,c),'*b');
            text(rangeK(k),Patients_P_pval(k,c),[' ' num2str(c)])
        end
    end
end

title('Probability Networks vs School_Years in Patients Group')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off

%% Plots networks that relate significantly with each Patient score

disp('PLOT: School_years > Patients')

    % Bonferroni correction
    Patients_P_pcorr=squeeze(Patients_P_pval(:,:)).*rangeK';  % Correct p-values by the number of states;
    [k,c]=ind2sub([length(rangeK),max(rangeK)],find(Patients_P_pcorr<0.05 & Patients_P_pcorr>0));

    if k
        
        [~,b]=sort(k);
        
        for h=b'
            figure
            subplot(4,2,1)
            V=Kmeans_results{k(h)}.C(c(h),:);
            %V=round(V*1e2)/1e2;
            V=V/max(abs(V));
            
            Render_Net(V,'AAL116','g')
            title(['K=' num2str(rangeK(k(h))) ' C=' num2str(c(h))])
            
            subplot(4,2,3) 
            Render_Net(V,'AAL116','g')
            view(0,0)            
            
            subplot(4,2,5)
            plot(squeeze(P(Index_Poor,k(h),c(h))),'.')
            xlabel('Duration (s)')
            [cc, p]=corrcoef(squeeze(P(Index_Poor,k(h),c(h))));
            title({['r=' num2str(cc(2))],['p=' num2str(p(2),'%1.2e')]})
            ylabel('School_Years')
            
            subplot(4,2,[2 4 6 8])
            hold on
            barh(find(V<0),V(V<0),'FaceColor','k','EdgeColor','none','Barwidth',.5)
            barh(find(V>=0),V(V>=0),'FaceColor','green','EdgeColor','none','Barwidth',.5)
            set(gca,'YTick',1:90,'Fontsize',7)
            set(gca,'YTickLabel',label90(:,:))
            set(gca,'Ydir','reverse')
            ylim([0 117])
            xlim([-1 1])
            title({'BOLD phase','projection into V1'},'Fontsize',12)
            grid on   
            
        end
    end
    
%% Control scores
Control_Scores=School_years(44:98,:);

disp('School_years: Controls')
 
for k=1:length(rangeK)
        
    for c=1:rangeK(end)
            
            P_Controls=squeeze(P(Index_Good,k,c));
            [cc_Controls,p_Controls]=corrcoef(P_Controls,Control_Scores);
            
            Controls_P_pval(k,c)=p_Controls(2);
            
         if Controls_P_pval(k,c)<(0.05)
            disp(['signif correlation with prob K=' num2str(rangeK(k)) ' c=' num2str(c)])
            disp(['cc=' num2str(cc_Controls(2)) ' p=' num2str(p_Controls(2))])
         end 
         
    end
end

%% Plot Control scores

figure
subplot(1,2,1)

semilogy(rangeK(1)-1:rangeK(end)+1,0.05*ones(1,length(rangeK)+2),'r--','LineWidth',1)
hold on
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',1)
semilogy(rangeK(1)-1:rangeK(end)+1,0.05./sum(rangeK)*ones(1,length(rangeK)+2),'b--','LineWidth',1)

for k=1:length(rangeK) 
    for c=1:rangeK(k)        
        if Controls_P_pval(k,c)>0.05
            semilogy(rangeK(k),Controls_P_pval(k,c),'*k');
            text(rangeK(k),Controls_P_pval(k,c),[' ' num2str(c)])
        end
        if Controls_P_pval(k,c)<0.05 && Controls_P_pval(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),Controls_P_pval(k,c),'*r');
            text(rangeK(k),Controls_P_pval(k,c),[' ' num2str(c)])
        end
        if Controls_P_pval(k,c)<(0.05/rangeK(k)) && Controls_P_pval(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),Controls_P_pval(k,c),'*g');
            text(rangeK(k),Controls_P_pval(k,c),[' ' num2str(c)])
        end
        if Controls_P_pval(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),Controls_P_pval(k,c),'*b');
            text(rangeK(k),Controls_P_pval(k,c),[' ' num2str(c)])
        end
    end
end

title('Probability Networks vs School_Years in Controls Group')
ylabel('Difference (p-value)')
xlabel('Number of clusters K')
xlim([rangeK(1)-1 rangeK(end)+1])
box off

%% Plots networks that relate significantly with each Control score

disp('PLOT: School_years > Controls')

% Bonferroni correction
    Controls_P_pcorr=squeeze(Controls_P_pval(:,:)).*rangeK'*length(Control_Scores);  % Correct p-values by the number of states;
    [k,c]=ind2sub([length(rangeK),max(rangeK)],find(Controls_P_pcorr<0.05 & Controls_P_pcorr>0));

    if k

        [~,b]=sort(k);
        
        u=0;
        for h=b'
            u=u+1;
            subplot(numel(b),2,2*u-1)
            V=Kmeans_results{k(h)}.C(c(h),:);
            Render_Net(V,'AAL116','r')
            title(['K=' num2str(rangeK(k(h))) ' C=' num2str(c(h))])
            
            subplot(numel(b),2,2*u)
            plot(squeeze(P(:,k(h),c(h))),'.')
            xlabel('Probability')
            ylabel('School_Years')
            title(['p=' num2str(Controls_P_pcorr(k(h),c(h)),'%1.2e') ])
            
        end
    end
 