function h=LDplot(tspan,ton,toff,ylo,yhi,options);
% LDPLOT plots light-dark bars for LD cycle

if nargin==6,
    switch options,
        case {'scaled'},
            a=gca;
            ylim=get(a,'YLim');
            height=0.1*(ylim(2)-ylim(1));
            yhi=ylim(1);
            ylo=yhi-height;
            co='k';
            shad=0;
        case {'shadow'},
            a=gca;
            ylim=get(a,'YLim');
            yhi=ylim(2);
            ylo=ylim(1);
            co=[0.85 0.85 0.85];
            shad=1;
    end
else
    co='k';
    shad=0;
end
      

t0=tspan(1);
tf=tspan(2);
ttot=tf-t0;
c0=mod(t0,24);
tons=t0+mod(ton-c0,24):24:tf;
toffs=t0+mod(toff-c0,24):24:tf;
if tons(1)<toffs(1) && tons(1)>t0,
    toffs=[t0 toffs];
elseif tons(1)==t0,
    tons=tons(2:end);
end
if toffs(end)>tons(end) && toffs(end)<tf,
    tons=[tons tf];
elseif toffs(end)==tf,
    toffs=toffs(1:end-1);
end

if length(tons)~=length(toffs)
    error('length(tons)~=length(toffs)')
end

for i=1:length(tons)
    h=fill([toffs(i) tons(i) tons(i) toffs(i)],[ylo ylo yhi yhi],co);
    if shad,
%         set(h,'EdgeColor','none');
        c=get(gca,'Children');
        set(gca,'Children',[c(2:end); c(1)]);
    else
        hold on 
        plot([t0 tf],[yhi yhi],'k');
    end
end