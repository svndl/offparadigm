function [rcaDataALL, W, A] = rcaRunProject(eegSrc, rcaDirs, timeCourseLen)

    %run RCA on eegSrc_p5
    nReg = 7;
    nComp = 3;
    
    dirResData = rcaDirs.results_Data;
    dirResFigures = rcaDirs.results_Figures;
    
    
    if (~exist(dirResData, 'dir'))
        mkdir(dirResData);
    end
    if (~exist(dirResFigures, 'dir'))
        mkdir(dirResFigures);
    end
   
    fileRCA = fullfile(dirResData, 'resultRCA.mat');
    
    % run RCA
    if(~exist(fileRCA, 'file'))
        [rcaDataALL, W, A, ~, ~, ~, ~] = rcaRun(eegSrc', nReg, nComp);
        save(fileRCA, 'rcaDataALL', 'W', 'A');
    else
        load(fileRCA);
    end
    % do the plots
    
    catDataAll = cat(3, rcaDataALL{:});
    muDataAll = nanmean(catDataAll, 3);
    muDataAll = muDataAll - repmat(muDataAll(1, :), [size(muDataAll, 1) 1]);
    semDataAll = nanstd(catDataAll, [], 3)/(sqrt(size(catDataAll, 3)));
    timeCourse = linspace(0, timeCourseLen, size(muDataAll, 1));
    close all;
    
    FontSize = 25;
    
    %% run the projections            
    for c = 1:nComp
        subplot(nComp, 2, 2*c - 1);
        xlabel('Time, ms') % x-axis label
        ylabel('Amplitude, V') % y-axis label
        h_xlabel = get(gca,'XLabel');
        set(h_xlabel,'FontSize', FontSize); 
        h_ylabel = get(gca,'YLabel');
        set(h_ylabel,'FontSize', FontSize);
        
        
        shadedErrorBar(timeCourse, muDataAll(:, c), semDataAll(:, c));hold on;
        title(['RC' num2str(c) ' time course']);
        subplot(nComp, 2, 2*c);
        plotOnEgi(A(:,c)); hold on;
    end
    %% Font Job
    set(findobj(gcf,'Type','text'), 'FontSize', FontSize)
    
    saveas(gcf, fullfile(dirResFigures, 'rcaComponentsAll'), 'fig');
    close(gcf);
    
    %% Project increment/decrement
    projectRC(eegSrc, W, A, nComp, {'Increment', 'Decrement'}, timeCourseLen, dirResFigures);
end