function [cc_V_yeo7,p_V_yeo7] = Overlap_LEiDA_Yeo(Centroids, rangeK)

%%%%%%%%%
%
%  Function to calculate overlap between LEiDA K-means results obtained in Glasser, 
%  DBS80 or AAL parcellation with the 7 resting-state Networks from Yeo et al. 2011
%
%  Inputs are:
%  Parcellation = 'AAL116' or 'glasser378' or 'dbs80'
%  N_areas      = 90 or 116 or 80 or 378
%  Extensions   = [] or '_contCORRECTED' or '_contCORRECTEDfiltered'
%  save         = 0 or 1 
%
%  It returns the results over the range of k-means clustering solutions.
%  
%  Note that the Yeo parcellation does not cover some subcortical areas 
%  nor the cerebellum (so LEiDA results obtained including these strucutres 
%  it may make sence to compute the overlap removing these areas.
%  i.e. in aal116 compare only the first 90
%  i.e. in glasser378 compare only the first 360 
% 
%  Joana Cabral October 2019
%
%%%%%%%%%%%

Parcellation='AAL116';
N_areas=90;

% 1) Define the Yeo Networks in the new Parcellation

% First load the mask of the chosen parcellation in MNI space 2mm
V_Parcels=struct2array(load(['ParcelsMNI2mm'],['V_' Parcellation]));
V_Parcels((V_Parcels>N_areas))=0;
% For a new one, use V_Parcels = niftiread('dbs80symm_2mm.nii.gz');

% Then load the mask of the Yeo parcellation in MNI space 2mm
V_Yeo=struct2array(load('ParcelsMNI2mm','V_Yeo7'));

%N_areas=max(V_Parcels(:));
N_Yeo=max(V_Yeo(:));
Yeo_in_new_Parcel=zeros(N_Yeo,N_areas);

% Create 7 vectors representing the 7 Yeo RSNs in new Parcellation scheme
for n=1:N_areas
    indn=V_Parcels==n;
    for Net=1:7
        Yeo_in_new_Parcel(Net,n)=numel(find(V_Yeo(indn)==Net))/sum(indn(:));
    end
end
clear V_Yeo V_Parcels indn

% 2) Compare with the LEiDA results  

cc_V_yeo7=zeros(length(rangeK),max(rangeK),N_Yeo) ;
p_V_yeo7=ones(length(rangeK),max(rangeK),N_Yeo) ;

for k=1:length(rangeK)
    
    VLeida=Centroids{k}.C(:,1:N_areas);
    %disp(['K= ' num2str(rangeK(k))])
    
    for Centroid=1:rangeK(k)
        
        % Here define what part of the Centroids to be compared
        V=VLeida(Centroid,:);
        
%         % To rescale between -1 and 1;
%         V=V-min(V);
%         V=2*(V/max(V));
%         V=V-1;
%         V=round(V*1e3)/1e3;
% 
%         % To set negative elements to zero
         V=V.*(V>0);
    
        for NetYeo=1:N_Yeo
            
            [cc, p]  = corrcoef(V(:),Yeo_in_new_Parcel(NetYeo,:));
            
            cc_V_yeo7(k,Centroid,NetYeo) = cc(2);
            p_V_yeo7(k,Centroid,NetYeo)  = p(2);
            
        end
    end
end


