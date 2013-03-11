function out = tThreeHalvs(T,dt,tau)
% function that returns f,f' where 
% f = T - dt*T.^(3/2) -dt*T -tau;
% fp = (1-dt) - (3/2)*dt*T^(1/2);

f = T - dt*T.^(3/2) -dt*T -tau;
fp = (1-dt) - (3/2)*dt*T.^(1/2);

out = {f,fp};