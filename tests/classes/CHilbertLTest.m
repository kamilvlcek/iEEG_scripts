function tests = CHilbertLTest
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.hilbert = CHilbertL('..\test-data\test_CHilb.mat');
end

function teardownOnce(testCase)
testCase.TestData.hilbert = [];
end

%% getenvelopes
function testGetanvelopesAll(testCase)
originalSize = size(testCase.TestData.hilbert.HFreqEpochs);
envelopes = testCase.TestData.hilbert.getenvelopes();
assertNotEmpty(testCase, envelopes);

verifyEqual(testCase, size(envelopes), originalSize);
verifyEqual(testCase, envelopes, testCase.TestData.hilbert.HFreqEpochs);
end

% getting only select channel
function testGetenvelopesChannel(testCase)
originalSize = size(testCase.TestData.hilbert.HFreqEpochs);
envelopes = testCase.TestData.hilbert.getenvelopes('channels', [1 3]);
assertNotEmpty(testCase, envelopes);

envSize = size(envelopes);
verifyEqual(testCase, envSize(2), 2);
verifyEqual(testCase, envSize([1 3 4]), originalSize([1 3 4]));
end

% getting only select frequencies
function testGetenvelopesFrequencies(testCase)
originalSize = size(testCase.TestData.hilbert.HFreqEpochs);
envelopes = testCase.TestData.hilbert.getenvelopes('frequencies', [1 3]);
assertNotEmpty(testCase, envelopes);

envSize = size(envelopes);
verifyEqual(testCase, envSize(3), 2);
verifyEqual(testCase, envSize([1 2 4]), originalSize([1 2 4]));

envelopes2 = testCase.TestData.hilbert.getenvelopes('frequencies',...
    testCase.TestData.hilbert.Hfmean([1 3]));
assertNotEmpty(testCase, envelopes2);
verifyEqual(testCase, envelopes, envelopes2);
end

% Getting only select categories
function testGetenvelopesCategories(testCase)
hilbert = testCase.TestData.hilbert;
% Validate passing wrong parameters
verifyError(testCase, @()hilbert.getenvelopes('category',[]),...
    'MATLAB:InputParser:ArgumentFailedValidation');

originalSize = size(testCase.TestData.hilbert.HFreqEpochs);
categoryNum = hilbert.PsyData.Categories(false);
categoryNames = hilbert.PsyData.CategoryName(categoryNum, []);

envelopes = hilbert.getenvelopes('categories', categoryNum([1 3]));
assertNotEmpty(testCase, envelopes);

envelopes2 = hilbert.getenvelopes('categories',categoryNames([1 3]));
assertNotEmpty(testCase, envelopes);
verifyEqual(testCase, envelopes, envelopes2);
end


% Getting only specific timeframes
function testGetenvelopesTime(testCase)
hilbert = testCase.TestData.hilbert;
verifyError(testCase, @()hilbert.getenvelopes('time',[]),...
    'MATLAB:InputParser:ArgumentFailedValidation');
end


% Rejecting bad epochs
function testGetenvelopesReject(testCase)
category = 1;
[~, ~, ~, iEpochs] = testCase.TestData.hilbert.CategoryData(category);
originalSize = size(testCase.TestData.hilbert.HFreqEpochs);
originalSize(4) = sum(iEpochs(:, 1));
envelopes = testCase.TestData.hilbert.getenvelopes('categories', category, 'reject', true);
verifySize(testCase, envelopes, originalSize);
end

%% Wilcox tests
function testWilcoxbaseline(testCase)
% ADD checks for size
wp = testCase.TestData.hilbert.wilcoxbaseline('baseline', [-0.1 0], ...
    'response', [0.01 0.1], 'frequencies', 1, 'categories', {'Ovoce'});
verifyNotEmpty(testCase, wp);
%verifySize(testCase, wp, size(tesCate.TestData.hilbert));
wp = testCase.TestData.hilbert.wilcoxbaseline('baseline', [-0.2 0],...
    'response', [0.01 0.1], 'frequencies', 1, 'categories', 1);
verifyNotEmpty(testCase, wp);
% test squeezing
end

function testWilcoxcategories(testCase)
originalSize = size(testCase.TestData.hilbert.HFreqEpochs);
wp = testCase.TestData.hilbert.wilcoxcategories({'Ovoce' 'Scene'},...
    'frequencies', 1);
assertNotEmpty(testCase, wp);
verifySize(testCase, wp, originalSize([1 2 3]));

wp = testCase.TestData.hilbert.wilcoxcategories({'Ovoce' 'Scene'},...
    'channels', 1, 'frequencies', 1);
assertNotEmpty(testCase, wp);
verifySize(testCase, wp, originalSize([1 2 3]));
% Test erroring in case multiple categories are passed
end

function testWilcoxaveragebaseline(testCase)
wp = testCase.TestData.hilbert.wilcoxaveragebaseline([1:10]);
end


function testWilcoxaveragecategories(testCase)

end

