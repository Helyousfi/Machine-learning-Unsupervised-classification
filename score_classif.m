

function [score_max] = score_classif(K,classif_exp,classif_real)


%Need to test all class label permutations to evaluate the highest possible score
%Use "perms" to generate label permutations

%reset between 1 and K
c = 1;
for i=unique(classif_exp)'
    classif_exp(classif_exp == i) = c;
    c = c + 1;
end

score_max = 0;
all_perms = perms(1:K);
for j=1:size(all_perms,1)
    classif_j = classif_exp;
    for p=1:K
        classif_j(classif_exp==p) = all_perms(j,p);
    end
    score = mean(classif_j == classif_real);
    if (score > score_max)
        score_max = score;
    end
end



end

