% takes in paramter values and computes the steady states and the stability
% of the found steady states. Returns everything via cells
function [pkm, actin, rna, hs, stability] = computeSS(j1,j2,j3,j4,j5)

% define returning vectors of length max of j1 and j2 (for looping through
% different parameter values. Only set up for changing j1 and j2
pkm = zeros(max(length(j1),length(j2)),1);
actin = zeros(max(length(j1),length(j2)),1);
rna = zeros(max(length(j1),length(j2)),1);
hs = zeros(max(length(j1),length(j2)),1);
stability = zeros(max(length(j1),length(j2)),1);

for j = 1:length(j1)
    for q = 1:length(j2)
        steadystates = @(pkm) j1(j)*j4.*(j2(q)+j3.*pkm).*pkm.*(1-pkm) - ...
            pkm.*(1+j2(q)+j3.*pkm).*(1+j4.*(j2(q)+j3.*pkm)./(1+j2(q)+...
            j3.*pkm).*pkm) - (j5*j2(q)+j5.*j3.*pkm).*(1+j4.*(j2(q)+j3.*pkm)...
            ./(1+j2(q)+j3.*pkm).*pkm);

        ssActin = @(pkm) (j2(q)+j3.*pkm)./(1+j2(q)+j3.*pkm);
        ssRNA = @(pkm) (j4.*ssActin(pkm).*pkm)./(1+j4.*ssActin(pkm).*pkm);
        ssHS = @(pkm) j5.*ssActin(pkm);

        % Jacobian entries
        jacMat = @(pkm, actin, rna, hs, stim) ...
                [-j1(j).*rna-1, 0, j1(j).*(1-pkm), -1;
                 1-actin, -(j2(q)+j3.*pkm)-1, 0, 0;
                 j4.*actin.*(1-rna), j4.*(pkm+stim).*(1-rna), -j4.*actin.*(pkm+stim)-1, 0;
                 0, j5, 0, -1];

        % symbolic solving for roots
        syms x;
        ssAns = double(vpasolve(steadystates(x)==0, x));
        
        
        for i = 1:length(ssAns)
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
    end
end
end
