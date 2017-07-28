function whm_analysis

    [curr_path, ~, ~] = fileparts(mfilename('fullpath'));
    rca_path.rootFolder = curr_path;
    %src data
    rca_path.srcEEG = fullfile(curr_path, 'EEG_whmHexFF_3Hz');
    %rca-ready EEG folder
    rca_path.rcaEEG = fullfile(curr_path, 'whmHex3Hz', 'rcaEEG');
        
    %% copy
    rca_path_c12 = rca_path; 
    rca_path_c12.results_Figures = fullfile(curr_path, 'results', 'whmHex3Hz_c12', 'Figures');
    rca_path_c12.results_Data = fullfile(curr_path, 'results', 'whmHex3Hz_c12', 'data');

    rca_path_c34 = rca_path;     
    rca_path_c34.results_Figures = fullfile(curr_path, 'results', 'whmHex3Hz_c34', 'Figures');
    rca_path_c34.results_Data = fullfile(curr_path, 'results', 'whmHex3Hz_c34', 'data');
    
    rca_path_c56 = rca_path; 
    rca_path_c56.results_Figures = fullfile(curr_path, 'results', 'whmHex3Hz_c56', 'Figures');
    rca_path_c56.results_Data = fullfile(curr_path, 'results', 'whmHex3Hz_c56', 'data');
    
    rca_path_c78 = rca_path;     
    rca_path_c78.results_Figures = fullfile(curr_path, 'results', 'whmHex3Hz_c78', 'Figures');
    rca_path_c78.results_Data = fullfile(curr_path, 'results', 'whmHex3Hz_c78', 'data');
    
    rca_path_c910 = rca_path;         
    rca_path_c910.results_Figures = fullfile(curr_path, 'results', 'whmHex3Hz_c910', 'Figures');
    rca_path_c910.results_Data = fullfile(curr_path, 'results', 'whmHex3Hz_c910', 'data');

    NS_DAR = 420;
    freqHz = 3;
    cycleDurationSamples = NS_DAR/freqHz;
    rcaDataOut = rcaReadRawEEG(rca_path);
    eegSrc = resampleData(rcaDataOut, cycleDurationSamples);

    eegSrc12 = eegSrc(:, [1 2]);
    eegSrc56 = eegSrc(:, [5 6]);
    eegSrc78 = eegSrc(:, [7 8]);
    
    eegSrc34 = eegSrc(:, [3 4]);
    eegSrc910 = eegSrc(:, [9 10]);
    
    
    timeCourseLen = 333;
    [rcaDataALL34, W34, ~]  = rcaRunProject(eegSrc34, rca_path_c34, timeCourseLen);
%     [rcaDataALL910, W910, A2]  = rcaRunProject(eegSrc910, rca_path_c910, timeCourseLen);
    
    [rcaDataALL12, W12, ~]  = rcaRunProject(eegSrc12, rca_path_c12, timeCourseLen);
    [rcaDataALL56, W56, ~]  = rcaRunProject(eegSrc56, rca_path_c56, timeCourseLen);
    [rcaDataALL78, W78, ~]  = rcaRunProject(eegSrc78, rca_path_c78, timeCourseLen);
    
    [mu12, s12] = prepData(rcaDataALL12);
    [mu34, s34] = prepData(rcaDataALL34);
    [mu56, s56] = prepData(rcaDataALL56);
    [mu78, s78] = prepData(rcaDataALL78);
        
    
    tc = linspace(0, timeCourseLen, size(rcaDataALL12{1, 1}, 1));

    hc1 = shadedErrorBar(tc, mu12(:, 1), s12(:, 1), 'r');hold on;
    hc2 = shadedErrorBar(tc, mu34(:, 1), s34(:, 1), 'b');hold on;
    hc3 = shadedErrorBar(tc, mu56(:, 1), s56(:, 1), 'g');hold on;
    hc4 = shadedErrorBar(tc, mu78(:, 1), s78(:, 1), 'k');hold on;   

    legend([hc1.patch, hc2.patch, hc3.patch, hc4.patch]);
end

function [muData, semData] = prepData(rcaDataIn)
    catData = cat(3, rcaDataIn{:});
    muData = nanmean(catData, 3);
    muData = muData - repmat(muData(1, :), [size(muData, 1) 1]);
    semData = nanstd(catData, [], 3)/(sqrt(size(catData, 3)));
end
