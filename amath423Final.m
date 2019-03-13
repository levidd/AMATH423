% AMATH 423
% Levi Davis and Oliver Speltz
close all; clc;

% Default parameter values
j1 = 50; %80; default
j2 = 0.05;
j3 = 0.5;
j4 = 0.16;
j5 = 0;
ta1 = 1500;
ta2 = 0.5;
ta3 = 60;
ta4 = 100;
% ta5 = ;

[pkm, actin, rna, hs, stability] = computeSS(j1,j2,j3,j4,j5);

pkm{1}
stability{1}

