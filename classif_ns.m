% TP 1 TS 326
% Name 1
tic
clear all
close all
load('cloud_data_1.mat');

N = size(data,1);
classif_true = data(:,3);
K = length(unique(classif_true));
X = data(:,1:2)';

colors = [1 0 0; 0 1 0; 0 0 1; 0 0 0; 1 1 0; 1 0 1; 0 1 1];


figure(1)
subplot 121
hold on
for j=1:N
    color = zeros(1,3);
    color(data(j,3)) = 1;
    plot(data(j,2),data(j,1),'color',colors(classif_true(j),:),'Marker','*');
end


Gi = rand(K,2);

use_pp = 0;

%% K-means++
if use_pp
    % Selection d'un echantillon aleatoire
    Gi(:,1) = X(1:2,ceil(rand(1)*N));
    
    % Distance D^2(x) des données X au premier Gi
    D2 = norm(X-Gi(:,1));
    %D2 = sum((X - Gi(:,1)).^2);
    
    % Vecteur d'association
    L = ones(N,1) %Tout associé au premier moment
    

    
    for i=2:K
        % Somme cumulée H(L)
        H = cumsum(D2)
        
        % Normalisation
        H = H/H(end);
        
        %Tirage d'un aleatoire d'un nouveau Gi
        Gi(:,i) = X(:,find(rand(1)<H,1));
        
        %Distance des plus proche
        for j=1:i
            Dk(1:2,j) = sum((X-Gi(:,j)).^2);
        end
        
        [D2,~] = find(Dk, [], 2);
    end
    
    Gi = Gi';
    
else
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%% K-means %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    

    %%Init seeds Gi = zi
    %Gi = zeros(3,2); %To modify
    %Gi = rand(K,2);
    Gi(:,1) = Gi(:,1)*1.1-0.6;
    Gi(:,2) = Gi(:,2)*1.2-0.4;
    hold on, plot(Gi(:,1),Gi(:,2),'o','LineWidth',5)
end


iter = 10;

I_inter = [];
I_intra = [];
I_tot = [];
classif_exp = ones(N,1);


%%K-means loop
for i=1:iter

   %%Compute closest association
   
   %distance = zeros(N,K); % distance de chaque point au barycentre 
   
   for n=1:N
       distance = zeros(1,K);
       for k=1:K
          distance(1,k) = norm(Gi(k,:)-data(n,1:2)); 
       end
       [~,index] = min(distance);
       classif_exp(n,1) = index;
   end
   
   %dist_min = zeros(1,N); % distance min de chaque point au barycentre
   
   %for n=1:N
   %    classif_exp(n,1) = find(distance(:,n) == min(distance(:,n)));
   %end
   

   
   %%Update barycenters
   Gi = zeros(4,2);
   
   for n=1:N
       Gi(classif_exp(n),1) = Gi(classif_exp(n),1) + data(n,1);
       Gi(classif_exp(n),2) = Gi(classif_exp(n),2) + data(n,2);
   end
   
   
   
   for k=1:K
       Gi(k,1) = Gi(k,1)/ length(find(classif_exp(:,1) == k));
       Gi(k,2) = Gi(k,2)/ length(find(classif_exp(:,1) == k));
   end
   
   
   
  
   
   
   
   
   %%Score classif
   score = score_classif(K,classif_exp,classif_true);
   
   %%plot res
   figure(1)
   subplot 122
   hold on
   %To modify
   for j=1:N
        color = zeros(1,3);
        color(classif_exp(3,1)) = 1;
        plot(data(j,2),data(j,1),'color',colors(classif_exp(j),:),'Marker','*');
   end
   
   %%Inertia
   [I_inter_i, I_intra_i, I_tot_i] = compute_inertia(K, Gi, data, classif_exp, classif_true);
   I_inter(end+1) = I_inter_i; I_intra(end+1) = I_intra_i; I_tot(end+1) = I_tot_i;
   
   figure(2);
   plot(1:i,I_inter,1:i,I_intra,1:i,I_tot);
   legend('Inter classes', 'Intra classes', 'Totale');
   
end

toc
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Hierarchical Ascendant Classification %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure(1)
% subplot 121
% hold on
% for j=1:N
%     plot(data(j,2),data(j,1),'color',colors(classif_true(j),:),'Marker','*');
% end
% subplot 122
% hold on
% for j=1:N
%     color = rand(1,3);
%     plot(data(j,2),data(j,1),'color',color,'Marker','*');
% end
% 
% clustering_type = 'ward'; %'min', 'max', 'bary', 'ward'
% 
% classif_exp = zeros(N,N);
% classif_exp(1,:) = 1:N;
% 
% Gi = data(:,1:2);
% 
% I_inter = [];
% I_intra = [];
% I_tot = [];
% for i=2:N
%    
%    %Start from previous classification
%    
%    %Compute exhaustive distances between all clusters
%    
%    %find minimum distance
%    
%    %%update classif by merging the two labels
%    
%    %%compute barycenters
% 
%    %%plot res
%    figure(1)
%    s2 = subplot(122);
%    hold on
%    %To modify
%    
%    %%Inertia
%    [I_inter_i, I_intra_i, I_tot_i] = compute_inertia(K, seeds, data, classif(i,:));
%    I_inter(end+1) = I_inter_i; I_intra(end+1) = I_intra_i; I_tot(end+1) = I_tot_i;
%    
%    figure(2);
%    plot(2:i,I_inter,2:i,I_intra,2:i,I_tot);
%    legend('Inter classes', 'Intra classes', 'Totale');
%    
% end
% 
% %%Compute score
% score_classif(K,classif(end,:)',data(:,3));
% 

