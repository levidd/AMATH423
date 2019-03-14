function dydt = neuronFireODENewTerm(t,y, stimFun,j1,j2,j3,j4,j5,j6,ta1,...
    ta2,ta3,ta4,ta5)
stim = stimFun(t);
dydt = neuronFireODE(t,y,stimFun,j1,j2,j3,j4,j5,ta1,ta2,ta3,ta4,ta5);
hs = y(4);
actin = y(2);
dydt = [dydt; 0] + [0;%-(j6*hs)./ta1;
            0;
            -j6*hs/ta3;
            (j5.*actin.*stim*(1 - hs)-hs)./ta5];

end