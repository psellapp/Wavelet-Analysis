function [scale] = waveletscale(n,linlog_flag,a0,a1)
%
% Generate array of scales for wavelet analysis.
% Created: Prabu, 8/13/2015
% modified: Prabu, 8/17/2015. Included option to specify scale limits
%
% n - length of original signal
% linlog_flag - set to 0 for linear scale, 1 for log scale
% a0,a1 - min,max values of scale a - defaults to 1 and n/2
%
if ~rem(n,2)==0
    disp('Error! n is not even.')
    scale = [];
    return
end

na = n/2;

if isempty(a0)
    a0 = 1;
end
if isempty(a1)
    a1 = na;
end
na = a1-a0+1;
scale = zeros(1,na);

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
