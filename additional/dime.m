function dime(height,width)
set(gcf, 'units', 'centimeters')
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) height*0.9677 width*0.9677]);
end