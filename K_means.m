% Algorithme de Llyod des K-moyennes
clear all
close all
clc

% Load the data : 160 lines and 3 columns
load('cloud_data_1.mat');

N = size(data,1);
classif_true = data(:,3);
K = length(unique(classif_true));
X = data(:,1:2)';

colors = [1 0 0; 
          0 1 0; 
          0 0 1; 
          0 0 0; 
          1 1 0; 
          1 0 1; 
          0 1 1];


figure(1)
subplot 121
hold on
for j=1:N
    color = zeros(1,3);
    color(data(j,3)) = 1;
    plot(data(j,2),data(j,1),'color',colors(classif_true(j),:),'Marker','*');
end


Gi = zeros(2,K);

%%
init = 2; % 0 : K-means, 1 : ginput  2 : K-means++,

switch init
    case 0       
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% K-means %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%        
        Gi(1,:) = rand(1,K)*(max(X(1,:))-min(X(1,:)))+min(X(1,:))*1.4;
        Gi(2,:) = rand(1,K)*(max(X(2,:))-min(X(2,:)))+min(X(2,:));
        hold on, plot(Gi(:,2),Gi(:,1),'o','LineWidth',5)
    case 1
    case 2
        % Selection d'un echantillon aleatoire
        Gi(:,1) = X(1:2,ceil(rand(1)*N));
        % Distance D^2(x) des données X au premier Gi
        D2 = norm(X-Gi(:,1));
        %D2 = sum((X - Gi(:,1)).^2);
        % Vecteur d'association
        %L = ones(N,1) %Tout associé au premier moment
        for i=2:K
            % Somme cumulée H(L)
            H = cumsum(D2)
            % Normalisation
            H = H/H(end);
            %Tirage d'un aleatoire d'un nouveau Gi
            Gi(:,i) = X(:,find(rand(1)<H,1));
            %Distance des plus proche
            for j=1:i
                Dk(1:2,j) = norm(X-Gi(:,j));
            end
            [D2,~] = min(Dk, [], 2);
        end
end




iter = 10;

I_inter = [];
I_intra = [];
I_tot = [];
classif_exp = ones(N,1);


%% K-means loop
for i=1:iter

   %%Compute closest association (we take a point and find the nearest Gi)  
   for n=1:N
       distance = zeros(1,K);
       for k=1:K
          distance(1,k) = norm(Gi(:,k)-X(1:2,n)); 
       end
       [~,index] = min(distance);
       classif_exp(n,1) = index;
   end
  
   
   %%Update barycenters
   
   
   for n=1:N
       Gi(1,classif_exp(n)) = Gi(1,classif_exp(n)) + X(1,n);
       Gi(2,classif_exp(n)) = Gi(2,classif_exp(n)) + X(2,n);
   end
   
   for k=1:K
       Gi(1,k) = Gi(1,k)/ length(find(classif_exp(:,1) == k));
       Gi(2,k) = Gi(2,k)/ length(find(classif_exp(:,1) == k));
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
   
   for k=1:K
       plot(Gi(2,k), Gi(1,k), 'o', 'color', colors(k,:), 'MarkerSize', 10); 
   end
   %title(sprintf('Iteration %1.2d %1.3f',i, score(i)));
   
   %%Inertia
   [I_inter_i, I_intra_i, I_tot_i] = compute_inertia(K, Gi, data, classif_exp, classif_true);
   I_inter(end+1) = I_inter_i; I_intra(end+1) = I_intra_i; I_tot(end+1) = I_tot_i;
   
   figure(2);
   plot(1:i,I_inter,1:i,I_intra,1:i,I_tot);
   legend('Inter classes', 'Intra classes', 'Totale');
   
end





