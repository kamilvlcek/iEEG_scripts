function tests = CiEEGDataTest
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.hilbert = CHilbertL('..\example-data\p073_PPA CHilbert 50-150Hz -0.2-0.8 refBipo Ep2019-11 FEX_CHilb.mat');
end

function teardownOnce(testCase)
testCase.TestData.hilbert = [];
end

function testCategoryDataRejection(testCase)
for category = testCase.TestData.hilbert.PsyData.Categories
    [~, ~, RjEpCh, iEpochy] = testCase.TestData.hilbert.CategoryData(category);
    verifyEqual(testCase, sum(sum(iEpochy)), sum(sum(~RjEpCh)));
end
end
