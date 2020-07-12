function tests = CiEEGDataTest
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.hilbert = CHilbertL('..\test-data\test_CHilb.mat');
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
