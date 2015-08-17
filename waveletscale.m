function [scale] = waveletscale(n,linlog_flag)
%
% Generate array of scales for wavelet analysis.
% Created: Prabu, 8/13/2015
%
% n - length of original signal
% linlog_flag - set to 0 for linear scale, 1 for log scale
%
if ~rem(n,2)==0
    disp('Error! n is not even.')
    scale = [];
    return
end

a0 = 1; na = n/2; a1 = na; scale = zeros(1,na);
if linlog_flag %log scale
    %     disp('log scale')
    da = (a1/a0)^(1/(na-1));
    %     delta = log(da);
    for i = 1:na
        a = a0*(da^(i-1));
        scale(i) = 2*pi()*a/n;
    end
else %linear scale
    %     disp('linear scale')
    da = (a1-a0+1)/na;
    for i = 1:na
        a = a0+((i-1)*da);
        scale(i) = 2*pi()*a/n;
    end
    %     delta = da./a;
end
% plot(scale,'*r')

end
