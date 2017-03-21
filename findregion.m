function [region, lreg] = findregion(neurnum,vlleft,vlright,dmleft,dmright);
% FINDREGION finds the SCN of given neurons
%   region = findregion(neurnum); finds the region
%   of SCN containing neurons NEURNUM, and returns a
%   vector REGION with entry 1 to 9:
%   1 : VL left
%   2 : VL right
%   3 : VL Left AND Right
%   4 : DM left
%   5 : DM right
%   9 : DM Left AND Right
%   The second output LREG is the number of neuron in each region.

[m,n]=size(neurnum);
if m>1 && n==1,
    neurnum = neurnum';
elseif m>1 && n>1,
    error('NEURNUM must be a vector')
end
ln=length(neurnum);

cmpvll=repmat(vlleft,1,ln)-repmat(neurnum,length(vlleft),1);
cmpvlr=repmat(vlright,1,ln)-repmat(neurnum,length(vlright),1);
cmpdml=repmat(dmleft,1,ln)-repmat(neurnum,length(dmleft),1);
cmpdmr=repmat(dmright,1,ln)-repmat(neurnum,length(dmright),1);

[bid, neurvll]=find(~cmpvll);
[bid, neurvlr]=find(~cmpvlr);
[bid, neurdml]=find(~cmpdml);
[bid, neurdmr]=find(~cmpdmr);
region=zeros(1,ln);
region(neurvll)=1;
region(neurvlr)=2+region(neurvlr);
region(neurdml)=4+region(neurdml);
region(neurdmr)=5+region(neurdmr);

lreg(1)=length(find(region==1));
lreg(2)=length(find(region==2));
lreg(3)=length(find(region==3));
lreg(4)=length(find(region==4));
lreg(5)=length(find(region==5));
lreg(6)=length(find(region==9));

