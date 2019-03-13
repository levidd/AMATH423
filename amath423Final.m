% AMATH 423
% Levi Davis and Oliver Speltz
close all; clc;

% Default parameter values
j1 = 57; %80; default
j2 = 0.05;
j3 = 0.5;
j4 = 0.16;
j5 = 0;
ta1 = 1500;
ta2 = 0.5;
ta3 = 60;
ta4 = 100;
% ta5 = ;

steadystates = @(pkm) j1*j4.*(j2+j3.*pkm).*pkm.*(1-pkm) - pkm.*(1+j2+j3.*pkm)...
    .*(1+j4.*(j2+j3.*pkm)./(1+j2+j3.*pkm).*pkm) - (j5*j2+j5.*j3.*pkm).*(1+j4.*...
    (j2+j3.*pkm)./(1+j2+j3.*pkm).*pkm);

ssActin = @(pkm) (j2+j3.*pkm)./(1+j2+j3.*pkm);
ssRNA = @(pkm) (j4.*ssActin(pkm).*pkm)./(1+j4.*ssActin(pkm).*pkm);
ssHS = @(pkm) j5.*ssActin(pkm);

% Jacobian entries
jacMat = @(pkm, actin, rna, hs, stim) ...
        [-j1.*rna-1, 0, j1.*(1-pkm), -1;
         1-actin, -(j2+j3.*pkm)-1, 0, 0;
         j4.*actin.*(1-rna), j4.*(pkm+stim).*(1-rna), -j4.*actin.*(pkm+stim)-1, 0;
         0, j5, 0, -1];

% symbolic solving for roots
syms x;
ssAns = double(vpasolve(steadystates(x)==0, x));

allSteadyStates = cell(length(ssAns) + 1, 5);
allSteadyStates{1,1} = "PKM";
allSteadyStates{1,2} = "Actin";
allSteadyStates{1,3} = "RNA";
allSteadyStates{1,4} = "HS";
allSteadyStates{1,5} = "Stability of SS";

for i = 1:length(ssAns)
    allSteadyStates{i+1,1} = ssAns(i);
    allSteadyStates{i+1,2} = ssActin(ssAns(i));
    allSteadyStates{i+1,3} = ssRNA(ssAns(i));
    allSteadyStates{i+1,4} = ssHS(ssAns(i));
    jac = jacMat(allSteadyStates{i+1,1}(1), allSteadyStates{i+1,2}(1),...
        allSteadyStates{i+1,3}(1), allSteadyStates{i+1,4}(1), 0);
    [~,~,isStable] = find(sign(eigs(jac)) == 1, 1); % returns -empty if no postive eigenvalues
    if isempty(isStable)
        allSteadyStates{i+1,5} = "Stable";
    else
        allSteadyStates{i+1,5} = "Unstable";
    end
end

disp(allSteadyStates);

