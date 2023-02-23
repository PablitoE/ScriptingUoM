% Test to analyse different focusing scores. Using information obtained
% from refocusing old hist files.

clear variables
classdata_file = 'D:\HoloICE\ICE-D_mainz\holograms\warm\seq12\seq12_cD.mat'; % ClassificationData file with hand-labelled or predicted particles.
config_file = 'D:\HoloICE\ICE-D_mainz\holograms\warm\seq12\kotelett.cfg';
output_dir = 'corrected_hists';

[dir_data, cD_filename, cD_ext] = fileparts(classdata_file);
output_path = fullfile(dir_data,output_dir);
ValidatedByHandFilename = 'validatedFocuses';
if ~endsWith(ValidatedByHandFilename,'.mat'), ValidatedByHandFilename = [ValidatedByHandFilename '.mat']; end
loadingValidated = fullfile(output_path,ValidatedByHandFilename);
load(loadingValidated,"validatedPrtcls","focusingScores","chosenFocus","namesOfFocusingScores")

kScores = find(cellfun(@(x) ~isempty(x), focusingScores(:,1)));
nScores = length(kScores);
nTypesScore = size(focusingScores,2);
selectedPeak = zeros(nScores,nTypesScore);
for ks = 1:nScores
   kfs = kScores(ks);
   selectedPeak(ks,:) = cellfun(@(x) find(x == max(x)), focusingScores(kfs,:));
end

goodSelection = selectedPeak == chosenFocus(kScores);
Performance = mean(goodSelection) * 100;
disp("Performance: "), disp( strjoin(namesOfFocusingScores))
disp(Performance)

%%%
% Output:
% Performance: 
% scores_mean ratioscores scores difscores scores_rescaled scores_mean_rescaled difscores_rescaled ratioscores_rescaled prtclPerimeterSobel
%    95.8209    7.1642   99.1045   67.1642   95.2239    8.0597   51.6418    5.6716   83.8806
   