function [ t, t1, t2, t3, params, case_no ] = NumberOfSections( a, as1, as2, b, h, c, asymptote_exists, ya )
    t = []; t1 = []; t2 = []; t3 = []; params = []; case_no = 0;
    % determine # of sections
    if (b<=as1 || (a>=as1 && b<=as1) || a>=as2 || ~asymptote_exists)
        t = a:h:b; params = [params, [a,b,ya]];
        case_no = 1;
    elseif (a<as1 && b<=as2)
        % lower asymptote
        t1 = a:h:as1; t2 = as1+h:h:as2;
        params = [[a,as1,ya];[as1+h,b,DESolution(as1+h,c)]];
        case_no = 2;
    elseif (a>=as1 && b>as2)
        % only upper
        t1 = a:h:as2; t2 = as2+h:h:b;
        params = [[a,as2,ya];[as2+h,b,DESolution(as2+h,c)]];
        case_no = 2;
    elseif (a<as1 && b>as2)
        % both
        t1 = a:h:as1; t2 = as1+h:h:as2; t3 = as2+h:h:b;
        params = [[a,as1,ya];[as1+h,as2,DESolution(as1+h,c)];[as2+h,b,DESolution(as2+h,c)]];
        case_no = 3;
    end
end

