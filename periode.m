function [per, perstd, pmaxf]=periode(t,f)

tol_per=.8;
sampling_rate=1/mean(diff(t));

% filter
% [b,a]=butter(4,1/12/sampling_rate);
% ff=filtfilt(b,a,f);
ff=f; % if BUTTER not available
% find maxima of f
fmoins=[Inf ff(1:end-1)];
fplus=[ff(2:end) Inf];
pmaxf=find(ff>=fmoins & ff>=fplus & abs((ff-max(ff))/(max(ff)-min(ff)))<tol_per);

per=mean(diff(t(pmaxf)));
perstd=std(diff(t(pmaxf)));