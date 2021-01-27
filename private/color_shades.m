% color_shades
% Function to create custom colormaps. Requires rgb function.
% Returns a matrix with matlab colors (0 to 1).
% Input must be a cell string like in the example below. To see the full
% list of colors names availables execute rgb chart (from rgb.m) function.
% Example: colormap_aux = color_shades({'red','yellow','green'})
% Author: Gonzalo A. Ferrada (gonzalo-ferrada@uiowa.edu)
% March 2019
%

function color_shades_all = color_shades(input_colors)

N        =  numel(input_colors);
Tshades  =  200;
Nshades  =  (Tshades-1)/(N-1);  % shades per color

% Assign colors
for i = 1:N
    %disp([' Assigning ' input_colors{i} ' color'])
    col(i,:)              =  rgb(input_colors{i});
    color_shades{i}(1,:)  =  col(i,:); %Assign first color of colormap
end

% Interpolate other colors of shade:
for i = 1:N-1
    incr    =  (col(i+1,:) - col(i,:)) ./ Nshades;
for j = 2:Nshades
    previous             =  color_shades{i}(j-1,:);
    color_shades{i}(j,:) =  previous + incr;
end
end

% Appending to create a single matrix:
color_shades_all = color_shades{1};

for i = 2:N
    row_i    =  size(color_shades_all,1) + 1;
    row_new  =  size(color_shades{i},1) + row_i - 1;
    color_shades_all(row_i:row_new,:) = color_shades{i};
end

end
