function tests = wilcox3dTest
tests = functiontests(localfunctions);
end

function testArguments(testCase)

end

function testFunctionality(testCase)
    testOnes = ones(10,10,20);
    testZeros = zeros(10,10,20);
    wp = wilcox3d(testOnes,testZeros);
    assertSize(testCase, wp, [10 10]);
    verifyTrue(testCase, all(wp(:)== wp(1)));
    verifyTrue(testCase, all(wp(:) < 0.01));
    
    wp2 = wilcox3d(testOnes, testOnes-0.01);
    verifyEqual(testCase, wp2, wp);
    
    testFlipFlop = repmat(cat(3,ones(10,10) + 0.1,ones(10,10) - 0.1), 1, 1, 10);
    wp = wilcox3d(testOnes, testFlipFlop);
    verifyTrue(testCase, all(wp(:) == 1));
end

function testFdr(testCase)

end

function testPaired(testCase)
end