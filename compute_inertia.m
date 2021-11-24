function [I_inter, I_intra, I_tot] = compute_inertia(K, seeds, data, classif_exp, classif_true)

%Variations
G = mean(seeds);
N = size(data,1);

%%% To modify
%Inter classes
I_inter = 0;

%Intra classes
I_intra = 0;

%Totale
I_tot = 0;

end
