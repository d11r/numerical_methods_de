function [ xre, xrie, xrke, re, rie, rke ] = GetAllSolutionsWithI(params, h, i)
    [xre, re]  = EulerMethod('Problem12', params(i,1), params(i,2), params(i,3), h);
    [xrie, rie] = ImprovedEulerMethod('Problem12', params(i,1), params(i,2), params(i,3), h);
    [xrke, rke] = RungeKuttaMethod('Problem12', params(i,1), params(i,2), params(i,3), h);    
end

