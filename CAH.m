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
   classif_exp(i,:) = classif_exp(i-1,:);
   
   %Compute exhaustive distances between all clusters
   dist = zeros(N);
   for n=1:N
       for m=1:N
           if norm(Gi(classif_exp(i,n),:)-Gi(classif_exp(i,m),:)) == 0
               dist(n,m) = inf;
           else
               dist(n,m) = norm(Gi(classif_exp(i,n),:)-Gi(classif_exp(i,m),:));
           end
       end
   end
  
   
   %find minimum distance
   for j=1:N
       [~,index] = min(dist(j,:));
       
       %%update classif by merging the two labels
       classif_exp(i,j) = classif_exp(i,index);
       
       %%compute barycenters
       Gi(i,:) =  Gi(classif_exp(i,j),:);
   end
       

   
   
   %%plot res
   figure(1)
   s2 = subplot(122);
   hold on
   %To modify

   
   %%Inertia
%    [I_inter_i, I_intra_i, I_tot_i] = compute_inertia(K, seeds, data, classif(i,:));
%    I_inter(end+1) = I_inter_i; I_intra(end+1) = I_intra_i; I_tot(end+1) = I_tot_i;
%    
%    figure(2);
%    plot(2:i,I_inter,2:i,I_intra,2:i,I_tot);
%    legend('Inter classes', 'Intra classes', 'Totale');
   
end

%%Compute score
%score_classif(K,classif(end,:)',data(:,3));