function tests = CHilbertLTest
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
% create and change to temporary folder
testCase.TestData.hilbert = CHilbertL('..\..\..\example-data\p073_PPA CHilbert 50-150Hz -0.2-0.8 refBipo Ep2019-11 FEX_CHilb.mat');
end

function teardownOnce(testCase)
testCase.TestData.hilbert = [];
end

function testGetanvelopes(testCase)
env = testCase.TestData.hilbert.getenvelopes();
verifyEqual(testCase, size(env), size(testCase.TestData.hilbert.HFreqEpochs));
end