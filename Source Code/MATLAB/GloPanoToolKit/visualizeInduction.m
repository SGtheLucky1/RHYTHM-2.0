%% Visualizing Arrhythmia Induction %%
% % % 
% Calculating phase
% % % cmosData = a.cmosData;
% % % phase = cell(1,4);
% % % for n = 1:4
% % % phase{n} = phaseMap(cmosData{n},0,(size(cmosData{n},3)-1)/a.Fs,a.Fs);
% % % end
% % % 
% % % % Calculating phase singularities
% % % nt = cell(1,4);
% % % for n = 1:4
% % %     nt{n} = zeros(size(cmosData{n},1),size(cmosData{n},2),size(cmosData{n},3));
% % %     for m = 1:size(cmosData{n},3)
% % %         % Singularities
% % %         nt{n}(:,:,m) = phaseSingularity(phase{n}(:,:,m));
% % %     end
% % % end


figure
% % % n = 7800;
for n = 10300:10:size(cmosData{1},3)
    % Membrane potential
    subplot(3,4,1)
    imagesc(cmosData{1}(:,:,n))
    title('LV')
    colormap('jet')
    subplot(3,4,4)
    imagesc(cmosData{2}(:,:,n))
    title('Posterior')
    colormap('jet')
    subplot(3,4,3)
    imagesc(cmosData{3}(:,:,n))
    title('RV')
    colormap('jet')
    subplot(3,4,2)
    imagesc(cmosData{4}(:,:,n))
    title('Anterior')
    colormap('jet')   
    
    % Phase
    subplot(3,4,5)
    imagesc(phase{1}(:,:,n))
    colormap('jet')
    caxis([-pi pi])
    subplot(3,4,8)
    imagesc(phase{2}(:,:,n))
    colormap('jet')
    caxis([-pi pi])
    subplot(3,4,7)
    imagesc(phase{3}(:,:,n))
    colormap('jet')
    caxis([-pi pi])
    subplot(3,4,6)
    imagesc(phase{4}(:,:,n))
    colormap('jet')
    caxis([-pi pi])

    % Singularities
    subplot(3,4,9)
    imagesc(nt{1}(:,:,n))
    colormap('jet')
    caxis([-1 1])
    subplot(3,4,12)
    imagesc(nt{2}(:,:,n))
    colormap('jet')
    caxis([-1 1])
    subplot(3,4,11)
    imagesc(nt{3}(:,:,n))
    colormap('jet')
    caxis([-1 1])
    subplot(3,4,10)
    imagesc(nt{4}(:,:,n))
    caxis([-1 1])
    colormap('jet')
    
    disp(n)
    pause
    
% % %     % Display index and pause
% % %     clc
% % %     disp('Left arrow (L) to decrement and right arrow (R) to increment.')
% % %     k=0;
% % %     while k == 0
% % %         k= waitforbuttonpress;
% % %         if ~strcmp(get(gcf,'currentcharacter'),'29');
% % %             k = 0;
% % %             n = n+1;
% % %         elseif ~strcmp(get(gcf,'currentcharacter'),'28');
% % %             k = 0;
% % %             n = n-1;
% % %         end
% % %     end
end