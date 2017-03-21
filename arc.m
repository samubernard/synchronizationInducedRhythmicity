function arc(pt1,pt2,orientation,col);
% ARC draws an arc between two points

midpt=(pt2+pt1)/2;
v=pt2-pt1;
if v(2)~=0,
    nor=[1 -v(1)/v(2)];
else
    nor=[0 -1];
end
nor=nor/norm(nor);
tau=1:10;

switch orientation,
    case 'left'
        midcirc=midpt-0.1*(norm(v))*nor;
        linearc=[pt1; midcirc; pt2];
        linterp=interp1([0 5 10],linearc,0:10,'spline');
        plot(linterp(:,1),linterp(:,2),'Color',col)
    case 'right'
        midcirc=midpt+0.1*(norm(v))*nor;
        linearc=[pt1; midcirc; pt2];
        linterp=interp1([0 5 10],linearc,0:10,'spline');
        plot(linterp(:,1),linterp(:,2),'Color',col)
end