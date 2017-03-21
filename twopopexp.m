% 2 pop experiment

rat=[0 0.1 0.2 0.5 0.8 1];
nrun=5;
spread=0.05;
N=100;
ic=repmat([0.7098 0.8036 1.3887 1.0724 0.6178 1.0190 0.9811 0.7144 0.35 0.35]',N,1); % initial conditions
mfper=[];
trans=72;
tend=200;

for i=1:6,
    for j=1:nrun,
        K=rand(N)<=0.1;
        connectivity = full((K>0)*ones(N,1)/N);
        c=repmat(connectivity,1,N);
        K=sparse(K);
        KK=K/N./c;
        KK(find(isnan(KK)))=0;
        sc = 1./(1+spread*randn(N,1));    % period scaling
        if rat(i)==0,
            sc=sc*1.2;
        elseif rat(i)<1,
            nrat=floor(rat(i)*N);
            sc(nrat:end)=sc(nrat:end)*1.2;
        end
        %%%%%% Simulation
        figure(2); clf;
        outx=10*floor(N*rand(5,1))+1;
        options=odeset('OutputS',outx,'OutputF','odeplot');
        if (trans>0)
            [t,y] = ode45(@BWnt,[0 trans/2 trans],ic,[],sc,KK,0,12,12);
            ic = y(end,:);
        end
        sol = ode45(@BWnt,[0 tend],ic,options,sc,KK,0,12,12);        % ODE
        %%%%%% Results with high resolution
        resol=0.1;
        xint=[0:resol:tend];
        yint=deval(sol,xint);
        mf=mean(yint(10:10:end,:));
        mfper(i,j)=periode(xint,mf);
        i
        mfper(i,j)
    end
end

