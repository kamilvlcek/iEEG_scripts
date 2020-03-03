function tests = CHilbertLTest
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.hilbert = CHilbertL('..\example-data\p073_PPA CHilbert 50-150Hz -0.2-0.8 refBipo Ep2019-11 FEX_CHilb.mat');
end

function teardownOnce(testCase)
testCase.TestData.hilbert = [];
end

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
% Getting only select categories
end