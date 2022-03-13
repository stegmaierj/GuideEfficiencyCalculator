%%
% GuideEfficiencyCalculator.
% Copyright (C) 2017 Christelle Etard, Swarnima Joshi, Johannes Stegmaier, Ralf Mikut, and Uwe Strähle
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the Liceense at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Please refer to the documentation for more information about the software
% as well as for installation instructions.
%
% If you use this application for your work, please cite the repository and one
% of the following publications:
%
% Etard, C., Joshi, S., Stegmaier, J., Mikut, R., & Strähle, U. "Tracking
% of indels by decomposition is a simple and effective method to assess
% efficiency of guide RNAs in zebrafish". Zebrafish, 14(6), 586-588. 2017.
%
%%

function [efficiency, stdDev, stdErr, figureHandle] = GuideEfficiencyCalculator(inputFileName, pamPosition, radius, previewMode)

    %% set additional parameters
    debugFigures = false;
    showCompleteSequence = false;
    if (min(radius) == -1 || min(pamPosition) == -1)
        showCompleteSequence = true;
    end
    
    %% prepare inputs for batch processing of a folder or a single file
    numInputFiles = 1;
    if (isfolder(inputFileName))
        inputFiles = dir(inputFileName, '*.scf');
        numInputFiles = length(inputFiles);
    else
        inputFiles{1} = inputFileName;
    end
    
    %% process all input files
    for f=1:numInputFiles

        %% set current input file
        if (isfolder(inputFileName))
            currentInputFileName = [inputFileName filesep inputFiles(f).name];
        else
            currentInputFileName = inputFileName;
        end

        %% setup default parameters
        reconstructMissingBases = true;
        lineWidth = 2;
        fontSize = 14;
        maxPhredScore = 60;

        %% read the scf file
        [A,C,G,T,P_A,P_C,P_G,P_T,~,PEAK_INDEX,BASE] = scfread(currentInputFileName);

        %% calculate the phred score
        phredScore = 10.^(-(-10*log10(max(max(max(P_A, P_C) ,P_G), P_T)))/10)';

        %% calculate total counts and relative values
        totalCounts = A + C + G + T;
        relativeA = A ./ totalCounts;
        relativeC = C ./ totalCounts;
        relativeG = G ./ totalCounts;
        relativeT = T ./ totalCounts;
        relativeA(isnan(relativeA)) = 0;
        relativeC(isnan(relativeC)) = 0;
        relativeG(isnan(relativeG)) = 0;
        relativeT(isnan(relativeT)) = 0;
        
        %% select the entire sequence if radius or pam location are not specified
        if (showCompleteSequence == true)
            pamPosition(f) = length(PEAK_INDEX);
            radius(f) = 200;
        end

        %% reconstruct missing peaks
        if (reconstructMissingBases == true)

            %% find the remaining peaks if the original basecalling fails. Phred score for these values is zero to indicate potentially unreliable values
            uncalledSequenceRange = (max(PEAK_INDEX)+1):length(A);
            uncalledA = A(uncalledSequenceRange);
            uncalledC = C(uncalledSequenceRange);
            uncalledG = G(uncalledSequenceRange);
            uncalledT = T(uncalledSequenceRange);
            [peakValues, peakIndices] = findpeaks(max([uncalledA, uncalledC, uncalledG, uncalledT], [], 2));
            uncalledSeq = [];
            for j=1:length(peakIndices)
                if (peakValues(j) == uncalledA(peakIndices(j)))
                    uncalledSeq = [uncalledSeq, 'A']; %#ok<*AGROW> 
                elseif (peakValues(j) == uncalledC(peakIndices(j)))
                    uncalledSeq = [uncalledSeq, 'C'];
                elseif (peakValues(j) == uncalledG(peakIndices(j)))
                    uncalledSeq = [uncalledSeq, 'G'];
                elseif (peakValues(j) == uncalledT(peakIndices(j)))
                    uncalledSeq = [uncalledSeq, 'T'];
                else
                    uncalledSeq = [uncalledSeq, 'N'];
                end
            end

            if (debugFigures == true)
                figure; hold on;
                plot(uncalledA, '-r');
                plot(uncalledC, '-g');
                plot(uncalledG, '-b');
                plot(uncalledT, '-k');
                plot(peakIndices, peakValues, '*m');
            end

            PEAK_INDEX = [PEAK_INDEX; (max(PEAK_INDEX)+1)+peakIndices];
            seqData = [BASE', uncalledSeq];
        else
            seqData = BASE';
        end

        %% sort peaks in descending order to get the maximum value peaks
        %sortedPeaks = sort([A'; C'; G'; T'], 1, 'descend');
        sortedPeaksRelative = sort([relativeA'; relativeC'; relativeG'; relativeT'], 1, 'descend');

        %% identify the maximum peak and normalize everything according to this maximum peak
        %maxValue = quantile(sortedPeaks(1,:), 0.95);
        %maxValue = max(sortedPeaks(1,:));
        %sortedPeaks = (sortedPeaks - 0.25) / (maxValue - 0.25);

        %% get the values to left and right of the pam location
        leftValues = sortedPeaksRelative(1, PEAK_INDEX(max(1, pamPosition(f)-radius):pamPosition(f)));
        rightValues = sortedPeaksRelative(1, PEAK_INDEX(pamPosition(f):min(pamPosition(f)+radius, length(PEAK_INDEX))));
        %leftValuesAbsolute = sortedPeaks(1, PEAK_INDEX(max(1, pamPosition(f)-radius):pamPosition(f)));
        %rightValuesAbsolute = sortedPeaks(1, PEAK_INDEX(pamPosition(f):min(pamPosition(f)+radius, length(PEAK_INDEX))));
        minLength = min(length(leftValues), length(rightValues));

        %% calculate median values and std error for left and right part of the pam
        leftValuesMedian = median(leftValues);
        rightValuesMedian = median(rightValues);
        leftValuesStd = std(leftValues);
        rightValuesStd = std(rightValues);
        leftValuesStdErr = std(leftValues) / sqrt(length(leftValues));
        rightValuesStdErr = std(rightValues) / sqrt(length(leftValues));

        %% open a new figure and plot the results
        fh = figure(1); clf; 
        
        if (previewMode == false)
            set(gcf, 'renderer', 'painters'); 
            set(gcf, 'Visible', 'off');
        end

        for i=1:2
            if (previewMode == false)
                if (i==1)
                    subplot(3,4,1:4);
                else
                    subplot(3,4,5:8);
                end
            else
                if (i==1)
                    subplot(2,1,1);
                    set(gca, 'Position', [0.025,0.575, 0.965, 0.4], 'Units', 'normalize');
                else
                    subplot(2,1,2);
                    set(gca, 'Position', [0.025,0.075, 0.965, 0.4], 'Units', 'normalize');
                end
            end
            hold on;

            %% get the maximum peak index
            minPeakIndex = min(PEAK_INDEX(max(1, pamPosition(f)-radius):min(length(PEAK_INDEX), pamPosition(f)+radius)));
            maxPeakIndex = max(PEAK_INDEX(max(1, pamPosition(f)-radius):min(length(PEAK_INDEX), pamPosition(f)+radius)));


            %% plot the intensities for the different bases
            if (i==1)
                maxOriginalIntensity = maxPhredScore / max([A(minPeakIndex:maxPeakIndex); C(minPeakIndex:maxPeakIndex); G(minPeakIndex:maxPeakIndex); T(minPeakIndex:maxPeakIndex)]);
                plot(minPeakIndex:maxPeakIndex, A(minPeakIndex:maxPeakIndex) * maxOriginalIntensity, '-g');
                plot(minPeakIndex:maxPeakIndex, C(minPeakIndex:maxPeakIndex) * maxOriginalIntensity, '-b');
                plot(minPeakIndex:maxPeakIndex, G(minPeakIndex:maxPeakIndex) * maxOriginalIntensity, '-k');
                plot(minPeakIndex:maxPeakIndex, T(minPeakIndex:maxPeakIndex) * maxOriginalIntensity, '-r');
            else
                maxOriginalIntensity = maxPhredScore / max([relativeA(minPeakIndex:maxPeakIndex); relativeC(minPeakIndex:maxPeakIndex); relativeG(minPeakIndex:maxPeakIndex); relativeT(minPeakIndex:maxPeakIndex)]);
                plot(minPeakIndex:maxPeakIndex, relativeA(minPeakIndex:maxPeakIndex) * maxOriginalIntensity, '-g');
                plot(minPeakIndex:maxPeakIndex, relativeC(minPeakIndex:maxPeakIndex) * maxOriginalIntensity, '-b');
                plot(minPeakIndex:maxPeakIndex, relativeG(minPeakIndex:maxPeakIndex) * maxOriginalIntensity, '-k');
                plot(minPeakIndex:maxPeakIndex, relativeT(minPeakIndex:maxPeakIndex) * maxOriginalIntensity, '-r');
            end

            %% plot the phred score and the pam location used for the calculations
            baseRange = max(1, pamPosition(f)-radius):min(length(PEAK_INDEX), pamPosition(f)+radius);
            plot(PEAK_INDEX(baseRange(baseRange <= length(phredScore))), phredScore(baseRange(baseRange <= length(phredScore))), '--k', 'LineWidth', lineWidth);
            plot([PEAK_INDEX(pamPosition(f)), PEAK_INDEX(pamPosition(f))], [0, maxPhredScore], '-k', 'LineWidth', lineWidth);

            %% switch x axis labels depending on the plot
            if (i==1)
                labels = cell(1, length(baseRange));
                for j=1:length(baseRange)
                    labels{1,j} = seqData(baseRange(j));
                end
                set(gca, 'XTick', PEAK_INDEX(baseRange) ,'XTickLabels', labels);
                xlabel('Base');
            else
                labels = cell(1, length(baseRange));
                stepSize = 10;
                currentLabel = 1;
                for j=1:stepSize:length(baseRange)
                    labels{1,currentLabel} = num2str(baseRange(j));
                    currentLabel = currentLabel+1;
                end
                set(gca, 'XTick', PEAK_INDEX(baseRange(1:stepSize:length(baseRange))),'XTickLabels', labels);
                xlabel('Base Number');
            end
            legend('A', 'C', 'G', 'T', 'PhredScore', 'PAM Location');
            ylabel('Phred Score / Intensity (a.u.)');
            set(gca, 'FontSize', fontSize);
            axis tight;
            box off;
        end

        if (previewMode == false)
            %% plot scatter plot
            subplot(3,4,9);
            hold on;

            %% identify the peak values
            peakValues = sortedPeaksRelative(1,PEAK_INDEX(max(1, pamPosition(f)-radius):min(length(PEAK_INDEX), pamPosition(f)+radius)));

            validRange = (max(1, pamPosition(f)-radius):min(length(PEAK_INDEX), pamPosition(f)+radius)) - pamPosition(f);
            
            %% plot the peak positions within the given radius
            plot(validRange, peakValues, '.k');
            plot([-radius, 0], [leftValuesMedian, leftValuesMedian], '-c', 'LineWidth', lineWidth);
            plot([0, radius], [rightValuesMedian, rightValuesMedian], '-m', 'LineWidth', lineWidth);
            ylabel('Max. Relative Intensity');
            axis([min(validRange), max(validRange), 0, 1]);
            box off;
            set(gca, 'FontSize', fontSize);

            %% box plots of the distribution of values before and after the pam
            subplot(3,4,10);
            hold on;            
            boxplot([leftValues(1:minLength)', rightValues(1:minLength)']);
            set(gca, 'XTick', [1,2], 'XTickLabel', {'<- PAM', 'PAM ->'});
            set(gca, 'YTick', [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2], 'YTickLabel', {'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0', ''});
            ylabel('Relative Intensity');
            axis([0.5, 2.5, 0, 1.0]);
            box off;
            set(gca, 'FontSize', fontSize);

            %% bar plots of the median values and the std. err. bars
            subplot(3,4,11);
            hold on
            bar((1:2)',[leftValuesMedian, rightValuesMedian]);
            errorbar(1:2,[leftValuesMedian, rightValuesMedian],[leftValuesStdErr, rightValuesStdErr], '.k', 'LineWidth', lineWidth)
            set(gca, 'XTick', [1,2], 'XTickLabel', {'<- PAM', 'PAM ->'});
            set(gca, 'YTick', [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2], 'YTickLabel', {'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0', ''});
            ylabel('Median of Max. Relative Intensity');
            axis([0.5, 2.5, 0, 1.0]);
            box off;
            set(gca, 'FontSize', fontSize);

            %% estimated cutting efficiency
            subplot(3,4,12);
            normalizedLeftEfficiency = 1 - (leftValuesMedian - 0.25) / (max(leftValuesMedian, rightValuesMedian) - 0.25);
            normalizedRightEfficiency = 1 - (rightValuesMedian - 0.25) / (max(leftValuesMedian, rightValuesMedian) - 0.25);
            normalizedLeftStd = std((leftValues - 0.25) / (max(leftValuesMedian, rightValuesMedian) - 0.25));
            normalizedRightStd = std((rightValues - 0.25) / (max(leftValuesMedian, rightValuesMedian) - 0.25));
            normalizedLeftStdErr = std((leftValues - 0.25) / (max(leftValuesMedian, rightValuesMedian) - 0.25)) / sqrt(length(leftValues));
            normalizedRightStdErr = std((rightValues - 0.25) / (max(leftValuesMedian, rightValuesMedian) - 0.25)) / sqrt(length(rightValues));

            hold on;
            bar((1:2)',[normalizedLeftEfficiency, normalizedRightEfficiency]);
            errorbar(1:2,[normalizedLeftEfficiency, normalizedRightEfficiency], [normalizedLeftStdErr, normalizedRightStdErr], '.k', 'LineWidth', lineWidth);
            set(gca, 'XTick', [1,2], 'XTickLabel', {'<- PAM', 'PAM ->'});
            set(gca, 'YTick', [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2], 'YTickLabel', {'0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0', ''});
            ylabel('Estimated Cutting Efficiency');
            axis([0.5, 2.5, 0, 1.0]);
            box off;
            set(gca, 'FontSize', fontSize);
        end

        %% export the png image
        [folder, file, ~] = fileparts(currentInputFileName);
        if (previewMode == false)
            resultsFolder = [folder filesep  'results' filesep];
            if (~exist(resultsFolder, 'dir'))
                mkdir(resultsFolder);
            end
            
            %% write result image
            outputFileName = strrep([folder filesep 'results' filesep file '_Results.png'], '//', '/');
            outFileNameCSV = strrep(outputFileName, '.png', '.csv');
            outFileNameTXT = strrep(outputFileName, '.png', '.txt');
            export_png_without_border(outputFileName, fh, 1920, 1080);
            close(fh);
            previewImage = imread(outputFileName);
            imwrite(255-previewImage, strrep(outputFileName, '_Results.png', '_Preview.png'));
            
            dlmwrite(outFileNameCSV, [leftValues(1:minLength)', rightValues(1:minLength)'], ';', 0, 0); %#ok<DLMWT> 
            prepend2file('Left Values; Right Values', outFileNameCSV, true);
            
            %% write the overview text file
            [h,p,ci,stats] = ttest2(leftValues, rightValues);
            fileID = fopen(outFileNameTXT,'w');            
            fprintf(fileID,'Median Left: %f\n', leftValuesMedian);
            fprintf(fileID,'Median Right: %f\n', rightValuesMedian);
            fprintf(fileID,'Std. Dev. Left: %f\n', leftValuesStd);
            fprintf(fileID,'Std. Dev. Right: %f\n', rightValuesStd);
            fprintf(fileID,'Std. Err. Left: %f\n', leftValuesStdErr);
            fprintf(fileID,'Std. Err. Right: %f\n', rightValuesStdErr);
            fprintf(fileID,'Estimated Cutting Efficiency: %f\n', max(normalizedLeftEfficiency, normalizedRightEfficiency));
            fprintf(fileID,'ttest h (0: not rejected, 1: rejected, alpha=0.05): %f\n', h);
            fprintf(fileID,'ttest p-value: %f\n', p);
            fprintf(fileID,'ttest confidence interval: [%f, %f]\n', ci(1), ci(2));
            fprintf(fileID,'ttest stats (tstat, df, sd): %f, %f, %f\n', stats.tstat, stats.df, stats.sd);
            fclose(fileID);
        end

        %% compute stats
        if (previewMode == false)
            efficiency = [normalizedLeftEfficiency, normalizedRightEfficiency];
            stdDev = [normalizedLeftStd, normalizedRightStd];
            stdErr = [normalizedLeftStdErr, normalizedRightStdErr];
        else
            efficiency = 0; stdDev = 0; stdErr = 0;
        end
    end

    %% return figure handle
    figureHandle = fh;
end