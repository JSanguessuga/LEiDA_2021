function Plot_Repertoire_DMTvsall

%%%%%%%%%
%
%  Function to plot the repertoire of LEiDA centroids for a given K
%  - Renders the areas in the smallest pole as patches in a 3D glass brain
% 
%
%  Joana Cabral Jan 2020 joanacabral@uminho.pt
%  Modified Dec 2020
%  Modified by João Sanguessuga, May 2021
%
%%%%%%%%%%%

% Define K
prompt='Qual é o K que pretende avaliar? K= ';
K=input(prompt);

stats_file = 'LEiDA_stats';
load(stats_file,'Kmeans_results','rangeK','P','P_pval','LT','LT_pval','Index_Good','Index_Poor')

Parcellation= 'AAL116';

% Same color code from the Yeo paper (CHANGE HERE THE COLORS)
% Modified Yeo Colors =[50 50 50; 240 90 0; 256 175 25; 125 50 125; 50 50 200; 150 200 100]./256;
YeoColor = [125 50 125; 50 50 200; 0 118 14; 196 58 250; 150 200 100; 256 175 25; 240 90 0]./256;
Yeo_names= {'Visual','Somatomotor','Dorsal Att.','Ventral Att.','Limbic','Frontoparietal','Default Mode'};

[cc_V_yeo7,p_V_yeo7] = Overlap_LEiDA_Yeo(Kmeans_results,rangeK);

% LEiDA networks colored according to closest RSN
Volume=struct2array(load('ParcelsMNI2mm',['V_' Parcellation]));
Brain_Mask=niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex=smooth3(Brain_Mask>0);

V=Kmeans_results{rangeK==K}.C;
V=V'./max(abs(V'));
clear Kmeans_results

N_areas=size(V,1);
Order=[1:2:N_areas N_areas:-2:2];

%%

figure
colormap(jet)
set(gcf, 'Position',  [100, 100, 1000, 800])

for c=1:K
    
    disp(['c=' num2str(c)])
    [~, net]= max(cc_V_yeo7(rangeK==K,c,:));
    
    % First Plot view from top
    subplot_tight(5,K,c)
    hold on
    cortexpatch=patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.01);
    isonormals(scortex,cortexpatch);
    
    % Color all areas with positive V with transparency proportional to
    % contribution
    n_pos=find(V(:,c)>0);
    if numel(n_pos)>0
        for region=1:length(n_pos)
            region_Vol=Volume==n_pos(region); % 3D volume with 1 only in the voxels of that area
            sregion=smooth3(region_Vol>0);
            
            if p_V_yeo7(rangeK==K,c,net)<0.05/K
                patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none')%,'FaceAlpha', V(n_pos(region),c))
            else
                patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none')%,'FaceAlpha', V(n_pos(region),c))
            end
        end
    else
        region_Vol=Volume>0 & Volume<=90; % 3D volume with 1 only in the voxels of that area
        sregion=smooth3(region_Vol>0);
        patch(isosurface(sregion,0.3), 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none')%,'FaceAlpha', V(n_pos(region),c))
    end
    if p_V_yeo7(rangeK==K,c,net)<0.05/K
        title({['V_C_' num2str(c), newline, Yeo_names{net}]} ) %,Yeo_names{c}})
    elseif c==1
        title({['V_C_' num2str(c), newline, 'Global Mode']} ) %,Yeo_names{c}})
    else
        title({['V_C_' num2str(c), newline, '']} ) %,Yeo_names{c}})
    end
    material dull; lighting gouraud
    daspect([1 1 1])
    view(-90,90)    % Top view    Side view:   view(0,0)
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
    %  Same but view from the side
    subplot_tight(5,K,c+K)
    hold on
    cortexpatch=patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
    reducepatch(cortexpatch,0.01);
    isonormals(scortex,cortexpatch);
    
    n_pos=find(V(:,c)>0);
    if n_pos
        for region=1:length(n_pos)
            region_Vol=Volume==n_pos(region); % 3D volume with 1 only in the voxels of that area
            sregion=smooth3(region_Vol>0);
            
            if p_V_yeo7(rangeK==K,c,net)<0.05/K
                patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none')%,'FaceAlpha', V(n_pos(region),c))
            else
                patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none')%,'FaceAlpha', V(n_pos(region),c))
            end
        end
    else
        region_Vol=Volume>0 & Volume<=90; % 3D volume with 1 only in the voxels of that area
        sregion=smooth3(region_Vol>0);
        patch(isosurface(sregion,0.3), 'FaceColor', [.6 .6 .6], 'EdgeColor', 'none')%,'FaceAlpha', V(n_pos(region),c))
        
        
    end
    material dull; lighting gouraud
    daspect([1 1 1])
    view(0,0)
    camlight;
    xlim([5 105]); ylim([5 85]); zlim([0 80])
    axis off
    
    % PLOT THE PROBABILITIES
    subplot_tight(5,K,c+2*K)
    
    Patients=squeeze(P(Index_Poor,rangeK==K,c))';  
    Controls=squeeze(P(Index_Good,rangeK==K,c))';
    
    if p_V_yeo7(rangeK==K,c,net)<0.05/K
        bar(1:2,[mean(Controls) mean(Patients)],'EdgeColor','w','FaceColor',YeoColor(net,:))
    else
        bar(1:2,[mean(Controls) mean(Patients)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    end
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Controls) mean(Patients)],[ std(Controls)/sqrt(numel(Controls)) std(Patients)/sqrt(numel(Patients))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'Good','Poor'})
    xtickangle(60)
    if P_pval(rangeK==K,c)<0.05 && P_pval(rangeK==K,c)>(0.05/K)
       plot(1.5,max([mean(Controls) mean(Patients)])+.01,'*r')
    end
    if P_pval(rangeK==K,c)<(0.05/K)
       plot(1.5,max([mean(Controls) mean(Patients)])+.01,'*g')
    end
    
    if c==1
        ylabel('Dominance probability')
    end
    hold off
    box off
    
    % PLOT THE LIFETIMES
    subplot_tight(5,K,c+3*K)
    
    Patients=squeeze(LT(Index_Poor,rangeK==K,c))';  
    Controls=squeeze(LT(Index_Good,rangeK==K,c))';
    
    if p_V_yeo7(rangeK==K,c,net)<0.05/K
        bar(1:2,[mean(Controls) mean(Patients)],'EdgeColor','w','FaceColor',YeoColor(net,:))
    else
        bar(1:2,[mean(Controls) mean(Patients)],'EdgeColor','w','FaceColor',[.5 .5 .5])
    end
    hold on
    % Error bar containing the standard error of the mean
    errorbar([mean(Controls) mean(Patients)],[ std(Controls)/sqrt(numel(Controls)) std(Patients)/sqrt(numel(Patients))],'LineStyle','none','Color','k')
    set(gca,'XTickLabel',{'Good','Poor'})
    xtickangle(60)
    if LT_pval(rangeK==K,c)<0.05 && LT_pval(rangeK==K,c)>(0.05/K)
       plot(1.5,max([mean(Controls) mean(Patients)])+.01,'*r')
    end
    if LT_pval(rangeK==K,c)<(0.05/K)
       plot(1.5,max([mean(Controls) mean(Patients)])+.01,'*g')
    end
    
    if c==1
        ylabel('Lifetime (seconds)')
    end
    hold off
    box off  
    
end

