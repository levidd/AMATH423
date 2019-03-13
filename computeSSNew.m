% takes in paramter values and computes the steady states and the stability
% of the found steady states. Returns everything via cells
function [pkm, actin, rna, hs, stability, indexes] = computeSSNew(j1,j2,j3,j4,j5)

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
            ssActin = @(pkm) (j2(q)+ j3.*pkm)./(1+j2(q)+j3.*pkm);
            ssHs = @(pkm) (j5(k).*ssActin(pkm))./(1+j5(k).*ssActin(pkm));
            ssRNA = @(pkm) (j4.*ssActin(pkm).*pkm)./(1+j4.*ssActin(pkm).*pkm);
            ssPkm =@(pkm) (j1(j).*ssRNA(pkm) - ssHs(pkm))./(1+j1(j).*ssRNA(pkm));
            

            % Jacobian entries
            jacMat = @(pkm, actin, rna, hs) ...
                    [-j1(j).*rna-1, 0, j1(j)-j1(j).*pkm, -1;
                    1-actin, -j2(q)-j3.*pkm-1, 0, 0;
                    j4.*actin-j4.*actin.*rna, j4.*pkm-j4.*pkm.*rna, -j4...
                        .*actin.*pkm - 1, 0;
                    0, j5(k)-j5(k).*hs, 0, -j5(k).*actin-1];
%             syms p a r h j1 j2 j3 j4 j5
%             odes = [j1*r.*(1-p)-p-h, (j2+j3*p)*(1-a)-a,...
%                 j4*a*p*(1-r)-r, j5*a*(1-h)-h];
%             J = jacobian(odes, [p, a, r, h]);

            % symbolic solving for roots
            syms x;
            ssAns = double(vpasolve(ssPkm(x)==0, x));

            for i = 1:length(ssAns)
                p = ssAns(i);
                if isreal(p) %&& p>=0
                    a = ssActin(p);
                    r = ssRNA(p);
                    h = ssHs(p);
                    pkm{index} = [pkm{index}, p];
                    actin{index} = [actin{index}, a];
                    rna{index} = [rna{index}, r];
                    hs{index}= [hs{index}, h];
                    jac = jacMat(p, a, r, h);
%                     jac = subs(J);
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
