% AMATH 423
% Levi Davis and Oliver Speltz
close all; clc;

% Default parameter values
j1 = 1:1:150; %80; default
j2 = 0.05;
j3 = 0.5;
j4 = 0.16;
j5 = 0;
ta1 = 1500;
ta2 = 0.5;
ta3 = 60;
ta4 = 100;
% ta5 = ;

[pkm, actin, rna, hs, stability, indexes] = computeSS(j1,j2,j3,j4,j5);

stable = [];
stableIndex = [];
unstable = [];
unstableIndex = [];
for j = 1:length(pkm)
    for i = 1:length(pkm{j})
        if stability{j}(i) > 0
            stable = [stable, pkm{j}(i)];
            stableIndex = [stableIndex, indexes{j}(i)];
        else
            unstable = [unstable, pkm{j}(i)];
            unstableIndex = [unstableIndex, indexes{j}(i)];
        end
    end
end


figure()
scatter(stableIndex, stable, '.'); hold on;
scatter(unstableIndex, unstable, '.');

