% takes in paramter values and computes the steady states and the stability
% of the found steady states. Returns everything via cells
function [pkm, actin, rna, hs, stability, indexes] = computeSS(j1,j2,j3,j4,j5)

% define returning vectors of length max of j1 and j2 (for looping through
% different parameter values. Only set up for changing j1 and j2
total = max(max(length(j1),length(j2)), length(j5));
pkm = cell(total,1);
actin = cell(total,1);
rna = cell(total,1);
hs = cell(total,1);
stability = cell(total,1);
indexes = cell(total, 1);

for j = 1:length(j1)
    for q = 1:length(j2)
        for k = 1:length(j5)
            % find index for ease of indexing
            index = max(max(j,q),k);
            %preload empty arrays for concatenating error bounds
            pkm{index} = [];
            actin{index} = [];
            rna{index} = [];
            hs{index} = [];
            stability{index} = [];

            % start the computations!
            steadystates = @(pkm) j1(j)*j4.*(j2(q)+j3.*pkm).*pkm.*(1-pkm) - ...
                pkm.*(1+j2(q)+j3.*pkm).*(1+j4.*(j2(q)+j3.*pkm)./(1+j2(q)+...
                j3.*pkm).*pkm) - (j5(k)*j2(q)+j5(k).*j3.*pkm).*(1+j4.*(j2(q)+j3.*pkm)...
                ./(1+j2(q)+j3.*pkm).*pkm);

            ssActin = @(pkm) (j2(q)+j3.*pkm)./(1+j2(q)+j3.*pkm);
            ssRNA = @(pkm) (j4.*ssActin(pkm).*pkm)./(1+j4.*ssActin(pkm).*pkm);
            ssHS = @(pkm) j5(k).*ssActin(pkm);

            % Jacobian entries
            jacMat = @(pkm, actin, rna, hs, stim) ...
                    [-j1(j).*rna-1, 0, j1(j).*(1-pkm), -1;
                     j3*(1-actin), -(j2(q)+j3.*pkm)-1, 0, 0;
                     j4.*actin.*(1-rna), j4.*(pkm+stim).*(1-rna), -j4.*actin.*(pkm+stim)-1, 0;
                     0, j5(k), 0, -1];

            % symbolic solving for roots
            syms x;
            ssAns = double(vpasolve(steadystates(x)==0, x));

            for i = 1:length(ssAns)
                p = ssAns(i);
                if isreal(p) && p>=0
                    a = ssActin(p);
                    r = ssRNA(p);
                    h = ssHS(p);
                    pkm{index} = [pkm{index}, p];
                    actin{index} = [actin{index}, a];
                    rna{index} = [rna{index}, r];
                    hs{index}= [hs{index}, h];
                    jac = jacMat(p, a, r, h, 0);
                    [~,~,isStable] = find(sign(eigs(jac)) == 1, 1); % returns -empty if no postive eigenvalues
                    if isempty(isStable)
                        stability{index} = [stability{index}, 1];
                    else
                        stability{index} = [stability{index}, -1];
                    end
                    indexes{index} = [indexes{index}, index];
                end
            end
        end
    end
end
end
