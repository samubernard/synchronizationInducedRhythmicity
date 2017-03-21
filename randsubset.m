function y=randsubset(v,c);
% RANDSUBSET extracts a random subset
%   Y=RANDSUBSET(V,C) extracts a random subset
%   from the vector V. If C is an integer larger than one and
%   smaller than length(V), the output Y is random subset of length C. If C
%   is a number between 0 and 1, Y is a random subset where each element of
%   V is in Y with probability C.

n=length(v);
if c==round(c) && c>1 && c<=n,
    s=randperm(n);
    s=sort(s(1:c));
    y=v(s);
elseif c>=0 && c<=1,
    p=rand(n,1)<=c;
    s=(1:n)';
    s=s.*p;
    s=s(find(s));
    y=v(s);
else
    error('Bad input: v or c don''match');
end