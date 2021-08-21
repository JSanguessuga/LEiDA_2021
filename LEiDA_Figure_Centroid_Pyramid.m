% Script to see Pyramid of LEiDA centroids colored according to best fit to YEO RSNs
%
% Joana Cabral Nov 2019
% Modified by João Sanguessuga, April 2021
% 
%  Parcellation = 'AAL116' or 'glasser378' or 'dbs80'
%  N_areas      = 90 or 116 or 80 or 378
 
Parcellation='AAL116';
N_areas=90; 
mink=2;
maxk=20;
rangeK = mink:maxk;

load LEiDA_Kmeans_results  Kmeans_results 

[cc_V_yeo7,p_V_yeo7] = Overlap_LEiDA_Yeo(Kmeans_results,rangeK);

YeoColor = [120 18 134; 70 30 180; 0 118 14; 196 58 250; 220 248 164; 230 148 34; 205 62 78]./256;
%color_vecs=[0.8 0.8 0.8; 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0; 1 .5 0];

% LEiDA networks colored according to closest RSN
Volume=struct2array(load('ParcelsMNI2mm',['V_' Parcellation])); 
Brain_Mask=niftiread('MNI152_T1_2mm_brain_mask.nii');
scortex=smooth3(Brain_Mask>0);

close all
figure(1)
set(gcf, 'Position',  [100, 100, 1000, 1000])
% figure(2)

for k=1:length(rangeK)
    
    VLeida=Kmeans_results{k}.C;
    
    disp(['K= ' num2str(rangeK(k))])
    
    for Centroid=1:rangeK(k)
        
        V=VLeida(Centroid,:);
        V=round(V*1e2)/1e2;  

        [cc_net, net]= max(cc_V_yeo7(k,Centroid,:));    
        %[p_net, net]= min(p_V_yeo7(k,Centroid,:)); 
        
        Centroid_Vol=zeros(size(Volume));      
        % To color all areas above zero
        for n=find(V>=0)
            Centroid_Vol(Volume==n)=1;
        end 
        sregion=smooth3(Centroid_Vol>0);
        
        figure(1) % View networks from Top
        subplot_tight(length(rangeK),rangeK(end),Centroid+(k-1)*rangeK(end),0.01)

        hold on
        % First plot a transparent cortex
        cortexpatch=patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
        reducepatch(cortexpatch,0.01);
        isonormals(scortex,cortexpatch);      
           
        if p_V_yeo7(k,Centroid,net)<0.05/k
            patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %'FaceAlpha', V(n))
        else
            % Color in Black if no significant overlap with any Yeo RSN
            patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %'FaceAlpha', V(n))
        end
        
        material dull
        lighting gouraud
        view(-90,90) % Top view
        daspect([1 1 1])
        camlight;
        xlim([5 105])
        ylim([5 85])
        zlim([10 80])
        axis off
       
%         figure(2) % View networks  from the side
%         set(gcf, 'Position',  [100, 100, 1000, 1000])
%         subplot_tight(length(rangeK),rangeK(end),Centroid+(k-1)*rangeK(end),0.01)
% 
%         hold on
%         % First plot a transparent cortex
%         cortexpatch=patch(isosurface(scortex,0.1), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none','FaceAlpha', 0.2);
%         reducepatch(cortexpatch,0.01);
%         isonormals(scortex,cortexpatch);      
%            
%         if p_V_yeo7(k,Centroid,net)<0.05/k
%             patch(isosurface(sregion,0.3), 'FaceColor', YeoColor(net,:), 'EdgeColor', 'none') %'FaceAlpha', V(n))
%         else
%             % Color in Black if no significant overlap with any Yeo RSN
%             patch(isosurface(sregion,0.3), 'FaceColor', 'k', 'EdgeColor', 'none') %'FaceAlpha', V(n))
%         end
%         
%         material dull
%         lighting gouraud
%         view(0,0) % Side view
%         daspect([1 1 1])
%         camlight;
%         xlim([5 105])
%         ylim([5 85])
%         zlim([10 80])
%         axis off
%                      
    end
end


