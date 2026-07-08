function tests = test_outerSegmentClasses()
tests = functiontests(localfunctions);
end

function testFactoryReturnsConcreteSubclasses(testCase)
testCase.verifyClass(osCreate('linear'), 'osLinear');
testCase.verifyClass(osCreate('biophys'), 'osBioPhys');
testCase.verifyClass(osCreate('identity'), 'osIdentity');
testCase.verifyClass(osCreate('displayrgb'), 'osDisplayRGB');
end

function testSharedOuterSegmentContract(testCase)
objects = {osLinear(), osBioPhys(), osIdentity(), osDisplayRGB()};
for ii = 1:numel(objects)
    obj = objects{ii};
    testCase.verifyTrue(isa(obj, 'outerSegment'));
    testCase.verifyEqual(obj.timeStep, 1e-4);
    testCase.verifyEqual(obj.noiseFlag, 'random');
    obj.noiseFlag = 'none';
    testCase.verifyEqual(obj.noiseFlag, 'none');
end
end

function testIdentityAndDisplayDataAccessors(testCase)
data = reshape(1:24, [2 3 4]);

identity = osIdentity();
identity.set('photon rate', data);
testCase.verifyEqual(identity.get('photon rate'), data);

displayRGB = osDisplayRGB();
displayRGB.set('rgb data', data);
testCase.verifyEqual(displayRGB.get('rgb data'), data);
end
