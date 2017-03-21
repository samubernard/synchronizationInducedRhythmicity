function [t,y]=sdeEuler(f,tspan,ic,noise,varargin);
% SDEEULER SDE solver Euler method

h=.02;
sh=sqrt(h);
t=tspan(1):h:tspan(2);
y(:,1)=ic;
for i=1:length(t)-1,
    y(:,i+1)=y(:,i)+h*feval(f,t,y(:,i),varargin{:})+sh*noise.*randn(size(noise)); % single-cell noise
%     y(:,i+1)=y(:,i)+h*feval(f,t,y(:,i),varargin{:})+sh*noise*randn(1); % global noise
end
y=y(:,1:5:end)';
t=t(1:5:end);