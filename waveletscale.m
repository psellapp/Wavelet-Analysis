function [scale,delta,da] = waveletscale(n,dt,ds,a0,a1,linlog_flag)
%
% Generate array of scales for wavelet analysis.
% Created: Prabu, 8/13/2015
% modified: Prabu, 8/17/2015.-Included option to specify scale limits
%           Prabu, 8/25/2015.-Mod to be more in line with matlab cwtft
%           function
%           Prabu, 9/3/2015. -output delta and da
%
% n - length of original signal
% linlog_flag - set to 0 for linear scale, 1 for log scale
% dt,ds - time and scale interval
% a0, a1 - min and max scale
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

% ===============================================
% if isempty(ds)
%     ds = 0.4875;% arbitrary value. Check Matlab documentation for function cwtft
% end
% if isempty(a0)
%     a0 = 2*dt;
% end
% if isempty(a1)
%     %     a1 = na*dt;
%     numScales = fix(log2(n)/ds)+1;
% else
%     numScales = length(a0:ds:a1);
% end
% 
% 
% if linlog_flag %log scale
%     %     disp('log scale')
%     scale = a0*2.^((0:numScales-1)*ds);
% else %linear scale
%     %     disp('linear scale')
%     scale = 2*pi*(a0+(0:numScales-1)*ds)/n;
% end
% 
% ================================================

na = a1-a0+1;
scale = zeros(1,na);

if linlog_flag %log scale
    %     disp('log scale')
    da = (a1/a0)^(1/(na-1));
        delta = log(da);
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
        delta = da./a;
end
figure(2); plot(scale,'*r')

end
