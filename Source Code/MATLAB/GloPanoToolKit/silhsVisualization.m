%% Silhouette Visualization %%

load('silhs1.mat')

r = size(silhs,1);
c = size(silhs,2);
z = size(silhs,3);

deg = 1:72/4:71;

for n = 1:length(deg)
    
    figure
    imagesc(silhs(:,:,deg(n)))
    colormap('gray')
    set(gca,'XTick',[],'YTick',[])
% % %     bw = silhs(:,:,n);
% % %     outline = bwperim(bw,8);
% % %     [or,oc] = find(outline);
% % %     oc = oc-c/2;
% % %     or = (or-r/2)*-1;
% % %     oz = zeros(length(oc),1);
% % %     
% % %     out = [oz oc or];
% % %     
% % %     Rz = [cosd(deg(n)) -sind(deg(n)) 0;
% % %         sind(deg(n)) cosd(deg(n)) 0;
% % %         0 0 1];
% % %     
% % %     out = out*Rz;
% % %     
% % %     
% % %     plot3(out(:,1),out(:,2),out(:,3),'c.')
    
end