function scnplot2(C,dim,cnt);
% SCNPLOT plots a 2D or a 3D graph of the SCN with data C
%   p=scnplot(C); plots a2D or a 3D graph of the SCN with data C.
%   USES UNSCALED DATA, ONLY FOR TIME SERIES



if nargin==2 & dim==3,
    C=flipdim(C,3);
    S1=size(C,1);
    S2=size(C,2);
    zmin=min(C(find(C)));
    zmax=max(C(find(C)));
%     [x y z C]=subvolume(C,[nan floor(2*S1/3) nan nan nan nan ]);
    [x y z C]=subvolume(C,[floor(S1/2)+1 nan floor(S2/3)+1 nan nan nan ]);
%     p=patch(isosurface(C~=0),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
    p=patch(isosurface(x,y,z,abs(C),0));
    p2=patch(isocaps(x,y,z,abs(C),0));
    isocolors(x,y,z,C,p2);
    set(p,'FaceColor','red','EdgeColor','none');
    set(p2,'FaceColor','flat','EdgeColor','none');
%     p2 = patch(isocaps(C~=0),'FaceColor','flat','EdgeColor','none');
%     hold on
%     isocolors(C,p2);
%     isocolors(C,p);
%     caxis([zmin zmax])
    view(-51,36); daspect([1 1 1]);
    axis(volumebounds(x,y,z,C))
%     camlight headlight; lighting gouraud;
    colorbar;
    
    %%%=============================================================
    
    % S1=size(C,1);
    % S2=size(C,2);
    % S3=size(C,3);
    % zmin=min(C(find(C~=0)));
    % zmax=max(C(find(C~=0)));
    % 
    % p3 = slice(C,[S1/4 S1/2 3*S1/4],[],[]);
    % % isocolors(C,p2);
    % caxis([zmin zmax])
    % % view(-17,44); daspect([1 1 1]);
    % axis tight
    % %     camlight; lighting p;
    % colorbar;
    
    
    %%% ABOVE 3D
else 
    %%%==============================================================
    %%% BELOW 2D
    
    S2=size(C,2);
    if S2>=3,
        C=squeeze(C(:,floor(S2/2),:));
    elseif S2==2,
        [m,layer]=max([length(find(C(:,1,:))) length(find(C(:,2,:)))]);
        C=squeeze(C(:,layer,:));
    else
        C=squeeze(C);
    end
    zmin=min(C(find(C~=0)));
    zmax=max(C(find(C~=0)));
    [p,q]=find(C);
    tb=[min(p)-1, max(p)+1, min(q)-1, max(q)+1];
    if numel(tb)~=0 && tb(1)<tb(2) && tb(3)<tb(4),
        X=C(tb(1):tb(2),tb(3):tb(4))';
        lc=length(colormap);
        %ceil(256*X)
%         X(find(X))=X(find(X))+0.01;
        image(ceil(cnt*lc*X))
%             cmap=contrast(X);
%             colormap(cmap);
%         caxis([zmin zmax])
    end
end