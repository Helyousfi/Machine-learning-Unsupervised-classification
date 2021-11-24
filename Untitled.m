% TP 1 TS 326
% Name 1
tic
clear all
close all

%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%
%%% K-means %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
iter = 7;

%%Init seeds Gi = zi
Gi = zeros(4,2); % To modify
for i = 1:4
   Gi(i,1) = (max(data(:,1)) - min(data(:,1)))*rand(1,1) + min(data(:,1));
   Gi(i,2) = (max(data(:,2)) - min(data(:,2)))*rand(1,1) + min(data(:,2));
end

I_inter = [];
I_intra = [];
I_tot = [];
classif_exp = ones(N,1);

%%K-means loop
for i=1:iter
   %%Compute closest association
   for j = 1:N
        dist = zeros(1,4);
        for k = 1:4
            dist(1,k) = sqrt( (data(j,1) - Gi(k,1))^2 + (data(j,2) - Gi(k,2))^2 );
        end
        [~,idx] = min(dist);
        classif_exp(j, 1) = idx;
   end
   %%Update barycenters
   Gi = zeros(4,2);
   for i = 1:N
        Gi(classif_exp(i),1) = Gi(classif_exp(i),1)  + data(i,1);
        Gi(classif_exp(i),2) = Gi(classif_exp(i),2) + data(i,2);
   end
   for j = 1:4
        Gi(j,1) = Gi(j,1)/length(find(classif_exp(:,1) == j));
        Gi(j,2) = Gi(j,2)/length(find(classif_exp(:,1) == j));
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Hierarchical Ascendant Classification %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot 121
hold on
for j=1:N
    plot(data(j,2),data(j,1),'color',colors(classif_true(j),:),'Marker','*');
end
subplot 122
hold on
for j=1:N
    color = rand(1,3);
    plot(data(j,2),data(j,1),'color',color,'Marker','*');
end

clustering_type = 'ward'; %'min', 'max', 'bary', 'ward'

classif_exp = zeros(N,N);
classif_exp(1,:) = 1:N;

Gi = data(:,1:2);

I_inter = [];
I_intra = [];
I_tot = [];
for i=2:N
   
   %Start from previous classification
   
   %Compute exhaustive distances between all clusters
   
   %find minimum distance
   
   %%update classif by merging the two labels
   
   %%compute barycenters

   %%plot res
   figure(1)
   s2 = subplot(122);
   hold on
   %% To modify
   
   %% Inertia
   [I_inter_i, I_intra_i, I_tot_i] = compute_inertia(K, seeds, data, classif(i,:));
   I_inter(end+1) = I_inter_i; I_intra(end+1) = I_intra_i; I_tot(end+1) = I_tot_i;
   
   figure(2);
   plot(2:i,I_inter,2:i,I_intra,2:i,I_tot);
   legend('Inter classes', 'Intra classes', 'Totale');
   
end

%%Compute score
score_classif(K,classif(end,:)',data(:,3));
