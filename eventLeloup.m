function [value,isterminal,direction]=eventLeloup(t,y,N,sc,K,L)
% EVENTLELOUP event function for GOLELOUP

xthr=50;
value=max(y)-xthr;
isterminal=1;
direction=1;
