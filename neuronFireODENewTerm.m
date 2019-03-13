function dydt = neuronFireODENewTerm(t,y, stimFun,j1,j2,j3,j4,j5,ta1,...
    ta2,ta3,ta4,ta5)
dydt = neuronFireODE(t,y,stimFun,j1,j2,j3,j4,j5,ta1,ta2,ta3,ta4,ta5);
hs = y(4);
actin = y(2);
dydt = [dydt; 0] - [hs./ta1;
            0;
            0;
            (j5.*actin - hs)./ta5];

end