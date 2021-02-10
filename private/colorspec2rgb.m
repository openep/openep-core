function rgb = colorspec2rgb(c)
% COLORSPEC2RGB converts a colorspec to truecolor.
% Usage:
%   rgb = colorspec2rgb(c)
% Where:
%   c is a color - either rgb or any of the letters in 'rgbwcmyk'
%
% COLORSPEC2RGB converts a character colorspec into a truecolor rgb
% representation. If c is not a character then it is returned as rgb - so
% that it will not affect index or truecolor variables that are input.
%
% Author: Nick Linton (2011)
% Modifications - 


rgbspec = [1 0 0;0 1 0;0 0 1;1 1 1;0 1 1;1 0 1;1 1 0;0 0 0];
cspec = 'rgbwcmyk';

% Deal with string color specifications.
if ischar(c),
  k = find(cspec==c(1));
  if isempty(k)
      error('COLORSPEC2RGB: Unknown color string.'); 
  end
  if k~=3 || length(c)==1,
    rgb = rgbspec(k,:);
  elseif length(c)>2,
    if strcmpi(c(1:3),'bla')
      rgb = [0 0 0];
    elseif strcmpi(c(1:3),'blu')
      rgb = [0 0 1];
    else
      error('COLORSPEC2RGB: Unknown color string.'); 
    end
  end
elseif isreal(c) && size(c,2)==3
    rgb = c;
else
    error('COLORSPEC2RGB: Unknown color format.');
end














