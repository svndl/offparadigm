function whm_analysis_SuperSet

    [curr_path, ~, ~] = fileparts(mfilename('fullpath'));
    rca_path.rootFolder = curr_path;
    %src data`
    rca_path.srcEEG = fullfile(curr_path, 'EEG_whmHexSS');
    %rca-ready EEG folder
    rca_path.rcaEEG = fullfile(curr_path, 'whmHexSuperSet', 'rcaEEG');
    
    % condition list
    
  
    %% F = 2.727Hz, f1-sup, inc/dec
    rca_path_c12 = rca_path; 
    rca_path_c12.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c12', 'Figures');
    rca_path_c12.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c12', 'data');

    %% F = 3Hz, f1-sup, inc/dec    
    rca_path_c34 = rca_path;     
    rca_path_c34.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c34', 'Figures');
    rca_path_c34.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c34', 'data');
    
    %% F = 3.75 Hz, f1-sup, inc/dec    
    rca_path_c56 = rca_path; 
    rca_path_c56.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c56', 'Figures');
    rca_path_c56.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c56', 'data');
    
    %% F = 2.727Hz, p1-sup, inc/dec        
    rca_path_c78 = rca_path;     
    rca_path_c78.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c78', 'Figures');
    rca_path_c78.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c78', 'data');
    
    %% F = 3Hz, p1-sup, inc/dec    
    rca_path_c910 = rca_path;         
    rca_path_c910.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c910', 'Figures');
    rca_path_c910.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c910', 'data');

    %% F = 3.75Hz, p1-sup, inc/dec
    rca_path_c1112 = rca_path; 
    rca_path_c1112.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c1112', 'Figures');
    rca_path_c1112.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c1112', 'data');
    
    %% F = 2.272Hz, p2-sup, inc/dec
    rca_path_c1314 = rca_path;     
    rca_path_c1314.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c1314', 'Figures');
    rca_path_c1314.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c1314', 'data');

    %% F = 3Hz, p2-sup, inc/dec    
    rca_path_c1516 = rca_path; 
    rca_path_c1516.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c1516', 'Figures');
    rca_path_c1516.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c1516', 'data');

    %% F = 3.75 Hz, p2-sup, inc/dec    
    rca_path_c1718 = rca_path;     
    rca_path_c1718.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c1718', 'Figures');
    rca_path_c1718.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c1718', 'data');
    
    %% F = 2.727Hz, f1-inf, inc/dec
    rca_path_c1920 = rca_path;         
    rca_path_c1920.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c1920', 'Figures');
    rca_path_c1920.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c1920', 'data');

    %% F = 3Hz, f1-inf, inc/dec
    rca_path_c2122 = rca_path; 
    rca_path_c2122.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c2122', 'Figures');
    rca_path_c2122.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c2122', 'data');
    
    %% F = 3.75Hz, f1-inf, inc/dec
    rca_path_c2324 = rca_path;     
    rca_path_c2324.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c2324', 'Figures');
    rca_path_c2324.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c2324', 'data');

    %% F = 2.727Hz, p1-inf, inc/dec    
    rca_path_c2526 = rca_path; 
    rca_path_c2526.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c2526', 'Figures');
    rca_path_c2526.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c2526', 'data');
    
    %% F = 3Hz, p1-inf, inc/dec
    rca_path_c2728 = rca_path;     
    rca_path_c2728.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c2728', 'Figures');
    rca_path_c2728.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c2728', 'data');
    
    %% F = 3.75Hz, p1-inf, inc/dec
    rca_path_c2930 = rca_path;         
    rca_path_c2930.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c2930', 'Figures');
    rca_path_c2930.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c2930', 'data');
    
    %% F = 2.727Hz, p2-inf, inc/dec    
    rca_path_c3132 = rca_path; 
    rca_path_c3132.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c3132', 'Figures');
    rca_path_c3132.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c3132', 'data');
    
    %% F = 3Hz, p2-inf, inc/dec
    rca_path_c3334 = rca_path;     
    rca_path_c3334.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c3334', 'Figures');
    rca_path_c3334.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c3334', 'data');
    
    %% F = 3.75Hz, p2-inf, inc/dec
    rca_path_c3536 = rca_path;         
    rca_path_c3536.results_Figures = fullfile(curr_path, 'results', 'whmHexSS_c3536', 'Figures');
    rca_path_c3536.results_Data = fullfile(curr_path, 'results', 'whmHexSS_c3536', 'data');

    NS_DAR = 420;
    freqHz = [2.727, 3, 3.75];
    cycleDurationSamples = round(NS_DAR./freqHz);
    rcaDataOut = rcaReadRawEEG(rca_path);
    
    rcaDataOut2_7Hz = rcaDataOut(:, [1 2 7 8 13 14 19 20 25 26 31 32]);
    rcaDataOut3Hz =  rcaDataOut(:, [3 4 9 10 15 16 21 22 27 28 33 34]);
    rcaDataOut3_75Hz =  rcaDataOut(:, [5 6 11 12 17 18 23 24 29 30 35 36]);
    
    eegSrc_27Hz = resampleData(rcaDataOut2_7Hz, cycleDurationSamples(1));
    eegSrc_3Hz = resampleData(rcaDataOut3Hz, cycleDurationSamples(2));
    eegSrc_375Hz = resampleData(rcaDataOut3_75Hz, cycleDurationSamples(3));
    

    eegSrc12 = eegSrc_27Hz(:, [1 2]);
    eegSrc34 = eegSrc_3Hz(:, [1 2]);    
    eegSrc56 = eegSrc_375Hz(:, [1 2]);
    eegSrc78 = eegSrc_27Hz(:, [3 4]);
    eegSrc910 = eegSrc_3Hz(:, [3 4]);
    eegSrc1112 = eegSrc_375Hz(:, [3 4]);
    eegSrc1314 = eegSrc_27Hz(:, [5 6]);
    eegSrc1516 = eegSrc_3Hz(:, [5 6]);    
    eegSrc1718 = eegSrc_375Hz(:, [5 6]);
    eegSrc1920 = eegSrc_27Hz(:, [7 8]);
    eegSrc2122 = eegSrc_3Hz(:, [7 8]);
    eegSrc2324 = eegSrc_375Hz(:, [7 8]);
    eegSrc2526 = eegSrc_27Hz(:, [9 10]);   
    eegSrc2728 = eegSrc_3Hz(:, [9 10]);
    eegSrc2930 = eegSrc_375Hz(:, [9 10]);
    eegSrc3132 = eegSrc_27Hz(:, [11 12]);
    eegSrc3334 = eegSrc_3Hz(:, [11 12]);
    eegSrc3536 = eegSrc_375Hz(:, [11 12]);
    
    
    timeCourseLen = round(1000./freqHz);
    
    [rcaDataALL12, ~, ~]  = rcaRunProject(eegSrc12, rca_path_c12, timeCourseLen(1));
    [rcaDataALL34, ~, ~]  = rcaRunProject(eegSrc34, rca_path_c34, timeCourseLen(2));    
    [rcaDataALL56, ~, ~]  = rcaRunProject(eegSrc56, rca_path_c56, timeCourseLen(3));
    [rcaDataALL78, ~, ~]  = rcaRunProject(eegSrc78, rca_path_c78, timeCourseLen(1));
    [rcaDataALL910, ~, ~]  = rcaRunProject(eegSrc910, rca_path_c910, timeCourseLen(2));
    [rcaDataALL1112, ~, ~]  = rcaRunProject(eegSrc1112, rca_path_c1112, timeCourseLen(3));
    [rcaDataALL1314, ~, ~]  = rcaRunProject(eegSrc1314, rca_path_c1314, timeCourseLen(1));    
    [rcaDataALL1516, ~, ~]  = rcaRunProject(eegSrc1516, rca_path_c1516, timeCourseLen(2));
    [rcaDataALL1718, ~, ~]  = rcaRunProject(eegSrc1718, rca_path_c1718, timeCourseLen(3));
    [rcaDataALL1920, ~, ~]  = rcaRunProject(eegSrc1920, rca_path_c1920, timeCourseLen(1));
    [rcaDataALL2122, ~, ~]  = rcaRunProject(eegSrc2122, rca_path_c2122, timeCourseLen(2));
    [rcaDataALL2324, ~, ~]  = rcaRunProject(eegSrc2324, rca_path_c2324, timeCourseLen(3));    
    [rcaDataALL2526, ~, ~]  = rcaRunProject(eegSrc2526, rca_path_c2526, timeCourseLen(1));
    [rcaDataALL2728, ~, ~]  = rcaRunProject(eegSrc2728, rca_path_c2728, timeCourseLen(2));
    [rcaDataALL2930, ~, ~]  = rcaRunProject(eegSrc2930, rca_path_c2930, timeCourseLen(3));
    [rcaDataALL3132, ~, ~]  = rcaRunProject(eegSrc3132, rca_path_c3132, timeCourseLen(1));
    [rcaDataALL3334, ~, ~]  = rcaRunProject(eegSrc3334, rca_path_c3334, timeCourseLen(2));    
    [rcaDataALL3536, ~, ~]  = rcaRunProject(eegSrc3536, rca_path_c3536, timeCourseLen(3));
 
    [mu12, s12] = prepData(rcaDataALL12);
    [mu34, s34] = prepData(rcaDataALL34);
    [mu56, s56] = prepData(rcaDataALL56);
    [mu78, s78] = prepData(rcaDataALL78);
    [mu910, s910] = prepData(rcaDataALL910);
    [mu1112, s1112] = prepData(rcaDataALL1112);
    [mu1314, s1314] = prepData(rcaDataALL1314);
    [mu1516, s1516] = prepData(rcaDataALL1516);
    [mu1718, s1718] = prepData(rcaDataALL1718);
    [mu1920, s1920] = prepData(rcaDataALL1920);
    [mu2122, s2122] = prepData(rcaDataALL2122);
    [mu2324, s2324] = prepData(rcaDataALL2324);
    [mu2526, s2526] = prepData(rcaDataALL2526);
    [mu2728, s2728] = prepData(rcaDataALL2728);
    [mu2930, s2930] = prepData(rcaDataALL2930);
    [mu3132, s3132] = prepData(rcaDataALL3132);
    [mu3334, s3334] = prepData(rcaDataALL3334);
    [mu3536, s3536] = prepData(rcaDataALL3536);


    
    %% plot different conditions
    
      
    tc1 = linspace(0, timeCourseLen(1), size(rcaDataALL12{1, 1}, 1));
    tc2 = linspace(0, timeCourseLen(2), size(rcaDataALL34{1, 1}, 1));
    tc3 = linspace(0, timeCourseLen(3), size(rcaDataALL56{1, 1}, 1));
    peakSign = -1;
    
    peakRange_time = [200, 250];
    peakRangeInd1 = find(tc1>peakRange_time(1) & tc1<peakRange_time(2));
    peakRangeInd2 = find(tc2>peakRange_time(1) & tc2<peakRange_time(2));
    peakRangeInd3 = find(tc3>peakRange_time(1) & tc3<peakRange_time(2));
    
    sgn12 = flipWaveform(mu12(:, 1), peakSign, peakRangeInd1);
    sgn34= flipWaveform(mu34(:, 1), peakSign, peakRangeInd2);
    sgn56 = flipWaveform(mu56(:, 1), peakSign, peakRangeInd3);
    sgn78 = flipWaveform(mu78(:, 1), peakSign, peakRangeInd1);
    sgn910 = flipWaveform(mu910(:, 1), peakSign, peakRangeInd2);
    sgn1112 = flipWaveform(mu1112(:, 1), peakSign, peakRangeInd3);
    sgn1314 = flipWaveform(mu1314(:, 1), peakSign, peakRangeInd1);
    sgn1516= flipWaveform(mu1516(:, 1), peakSign, peakRangeInd2);
    sgn1718 = flipWaveform(mu1718(:, 1), peakSign, peakRangeInd3);
    sgn1920 = flipWaveform(mu1920(:, 1), peakSign, peakRangeInd1);
    sgn2122 = flipWaveform(mu2122(:, 1), peakSign, peakRangeInd2);
    sgn2324 = flipWaveform(mu2324(:, 1), peakSign, peakRangeInd3);
    sgn2526 = flipWaveform(mu2526(:, 1), peakSign, peakRangeInd1);
    sgn2728 = flipWaveform(mu2728(:, 1), peakSign, peakRangeInd2);
    sgn2930 = flipWaveform(mu2930(:, 1), peakSign, peakRangeInd3);
    sgn3132 = flipWaveform(mu3132(:, 1), peakSign, peakRangeInd1);
    sgn3334 = flipWaveform(mu3334(:, 1), peakSign, peakRangeInd2);
    sgn3536 = flipWaveform(mu3536(:, 1), peakSign, peakRangeInd3);
    
    
    close all;
    % 2.72Hz [1 2 7 8 13 14 19 20 25 26 31 32];
    legendList = {'F1 Lower', 'P1 Lower', 'P2 Lower', 'F1 Upper', 'P1 Upper', 'P2 Upper'};
    [legendColor_27, gcf_27] = plotMultipleComponents(tc1', sgn12*mu12, sgn12*s12, ...
        sgn78*mu78, sgn78*s78, sgn1314*mu1314, sgn1314*s1314, ...
        sgn1920*mu1920, sgn1920*s1920, sgn2526*mu2526, sgn2526*s2526,...
        sgn3132*mu3132, sgn3132*s3132);
    
    legend([legendColor_27{:}], legendList(:));
    saveas(gcf_27, fullfile(rca_path.rootFolder, 'ALL_2.7Hz.fig'));    
    close(gcf_27);
    
    % 3 Hz [3 4 9 10 15 16 21 22 27 28 33 34];
    [legendColor_3, gcf_3] = plotMultipleComponents(tc2', sgn34*mu34, sgn34*s34, ...
        sgn910*mu910, sgn910*s910, sgn1516*mu1516, sgn1516*s1516, ...
        sgn2122*mu2122, sgn2122*s2122, sgn2728*mu2728, sgn2728*s2728, ...
        sgn3334*mu3334, sgn3334*s3334);
    legend([legendColor_3{:}], legendList(:));
    saveas(gcf_3, fullfile(rca_path.rootFolder, 'ALL_3Hz.fig'));    
    close(gcf_3);

    % 3.75 Hz [5 6 11 12 17 18 23 24 29 30 35 36];    
    [legendColor_37, gcf_37] = plotMultipleComponents(tc3', sgn56*mu56, sgn56*s56, ... 
        sgn1112*mu1112, sgn1112*s1112, sgn1718*mu1718, sgn1718*s1718, ...
        sgn2324*mu2324, sgn2324*s2324, sgn2930*mu2930, sgn2930*s2930, ...
        sgn3536*mu3536, sgn3536*s3536);
    legend([legendColor_37{:}], legendList(:));
    saveas(gcf_37, fullfile(rca_path.rootFolder, 'ALL_3.7Hz.fig'));    
    close(gcf_37);
    
    
    group_upper_27_mu = [mu12(:, 1) mu78(:, 1) mu1314(:, 1)];
    group_lower_27_mu = [mu1920(:, 1) mu2526(:, 1) mu3132(:, 1)];   
    group_upper_27_s = [s12(:, 1) s78(:, 1) s1314(:, 1)];
    group_lower_27_s = [s1920(:, 1) s2526(:, 1) s3132(:, 1)];
    
    group_upper_3_mu = [mu34(:, 1) mu910(:, 1) mu1516(:, 1)];
    group_lower_3_mu = [mu2122(:, 1) mu2728(:, 1) mu3334(:, 1)];   
    group_upper_3_s = [s34(:, 1) s910(:, 1) s1516(:, 1)];
    group_lower_3_s = [s2122(:, 1) s2728(:, 1) s3334(:, 1)];
    
    [h1, h2] = subplotMultipleComponents(tc1, {group_upper_27_mu, group_upper_27_s}, ...
        {group_lower_27_mu, group_lower_27_s}, {'Lower 2.73 Hz', 'Upper 2.73Hz'});
    saveas(h1, 'Group  2.73 Hz _byGroup.fig');
    saveas(h2, 'Group  2.73 Hz _byRegion.fig');

    [h3, h4] = subplotMultipleComponents(tc2, {group_upper_3_mu, group_upper_3_s}, ...
        {group_lower_3_mu, group_lower_3_s}, {'Lower 3Hz', 'Upper 3Hz'});
    saveas(h3, 'Group  3 Hz _byGroup.fig');
    saveas(h4, 'Group  3 Hz _byRegion.fig');
    
end

function [muData, semData] = prepData(rcaDataIn)
    catData = cat(3, rcaDataIn{:});
    muData = nanmean(catData, 3);
    muData = muData - repmat(muData(1, :), [size(muData, 1) 1]);
    semData = nanstd(catData, [], 3)/(sqrt(size(catData, 3)));
end

function signFlip = flipWaveform(w, peakSign, peakRange)    
    signFlip = 1;
%     [~, l] = findpeaks(w);
%     peaksInRange = intersect(l, peakRange);
%     if (~isempty(peaksInRange)) 
%         peakVar = max(w(peaksInRange));  
%         signFlip = peakSign*sign(peakVar);
%     end
end



