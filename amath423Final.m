% AMATH 423
% Levi Davis and Oliver Speltz
close all; clc;

% Default parameter values
j1 = 2:1:150; %80; default
j2 = 0.05;
j3 = 0.5;
j4 = 0.16;
j5 = 0;
ta1 = 1500;
ta2 = 0.5;
ta3 = 60;
ta4 = 100;
ta5 = 1;
% 
% strength = 125;
% stimFun = @(t) 0.003 + (t>=0).*strength - strength.*(t>30); % basal value
%         % of 0.003. Strength of square wave between time 0 and 30
% 
% initialGuess = [0, 0, 0];
% 
% % ODE solution for simple model. 
% [T,Y] = ode45(@(t,y) neuronFireODE(t,y,stimFun,j1,j2,j3,j4,j5,ta1,ta2,...
%     ta3,ta4,ta5), [0,1000], initialGuess);
% disp('Trying to solve new ODE system');
% opts = odeset('RelTol', 1e-1, 'AbsTol', 1e-2);
% [T2,Y2] = ode45(@(t,y) neuronFireODENewTerm(t,y,stimFun,j1,j2,j3,j4,j5,...
%     ta1,ta2,ta3,ta4,ta5), [0,500], [initialGuess, 0], opts);
% 
% figure(1)
% plot(T, Y(:,1));hold on;
% plot(T2, Y2(:,1));
% legend('Paper Model', 'OurModel');
% xlabel('Time (s)'); ylabel('[PKM]');



% make sure to set this to the parameter being varied. For plotting
% purposes
varying = j1;

% compute values and determine stability
[pkm, actin, rna, hs, stability, indexes] = computeSSNew(j1,j2,j3,j4,j5);

stable = [];
stableIndex = [];
unstable = [];
unstableIndex = [];
for j = 1:length(pkm)
    for i = 1:length(pkm{j})
        if stability{j}(i) > 0
            stable = [stable, pkm{j}(i)];
            stableIndex = [stableIndex, varying(indexes{j}(i))];
        else
            unstable = [unstable, pkm{j}(i)];
            unstableIndex = [unstableIndex, varying(indexes{j}(i))];
        end
    end
end


figure(1)
scatter(stableIndex, stable, '.'); hold on;
scatter(unstableIndex, unstable, '.');
legend('Stable', 'Unstable');
ylabel('[PKM]'); xlabel('j1');

