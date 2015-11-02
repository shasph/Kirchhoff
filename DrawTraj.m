function DrawTraj()
clear all;
close all;

load('Traj.mat');
a = 15;
b = 9;
t = 0.6;
N = size(DD,1);
% N = 108;%deleteCell(DD);
VertexData_0 = [-a/2 -b/2 0;
    -a/2 b/2 0;
    a/2 b/2 0;
    a/2 -b/2 0;
    -a/2 -b/2 t;
    -a/2 b/2 t;
    a/2 b/2 t;
    a/2 -b/2 t];
n_side = 4;
n_ver = 2*n_side;
for k=1:N
    G = DD{k,1};
%     P(k,:)=G(1:3,4)';
    %Take the frames
    if(mod(k,6)==1)
        po = G(1:3,4)';
        po = [po(1);po(2);po(3)];
        for i_ver=1:n_ver
            VertexData(i_ver,:) = po' + VertexData_0(i_ver,:)*G(1:3,1:3)';%*ori(0,0,pi/2);
        end
        
        % Side Patches
        for i_pat=1:n_side-1
            Index_Patch1(i_pat,:) = [i_pat,i_pat+1,i_pat+1+n_side,i_pat+n_side];
        end
        Index_Patch1(n_side,:) = [n_side,1,1+n_side,2*n_side];
        
        for i_pat=1:n_side
            
            % Side patches data
            PatchData1_X(:,i_pat) = VertexData(Index_Patch1(i_pat,:),1);
            PatchData1_Y(:,i_pat) = VertexData(Index_Patch1(i_pat,:),2);
            PatchData1_Z(:,i_pat) = VertexData(Index_Patch1(i_pat,:),3);
        end
        
        % Draw side patches
        %     figure(1);
%         figure;
        set(gca,'fontsize',14);
%         hold on;
        h1 = patch(PatchData1_X,PatchData1_Y,PatchData1_Z,'y');
        set(h1,'FaceLighting','phong','EdgeLighting','phong','EdgeColor','none','FaceColor',[0.91,0.62,0.06]);
        set(h1,'EraseMode','normal');
        
        % Bottom Patches
        Index_Patch2(1,:) = [1:n_side];
        Index_Patch2(2,:) = [n_side+1:2*n_side];
        
        for i_pat=1:2
            % Bottom patches data
            PatchData2_X(:,i_pat) = VertexData(Index_Patch2(i_pat,:),1);
            PatchData2_Y(:,i_pat) = VertexData(Index_Patch2(i_pat,:),2);
            PatchData2_Z(:,i_pat) = VertexData(Index_Patch2(i_pat,:),3);
        end
        
        % Draw bottom patches
        %     figure(1);
        h2 = patch(PatchData2_X,PatchData2_Y,PatchData2_Z,'y');
        set(h2,'FaceLighting','phong','EdgeLighting','phong','FaceColor',[0.91,0.62,0.06]);%[0,0.5,1]
        set(h2,'EraseMode','normal');
        view(3);
%         camlight;
        hold on;
    end
    if (po(3)<0)
        break;
    end
end

set(gca,'FontSize',18);
axis equal;
end