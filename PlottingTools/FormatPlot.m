ax = gca;
% make all lines thicker
lines = findobj(gcf,'Type','Line');
for i = 1:numel(lines)
  lines(i).LineWidth = 2.0;
end
% Make all ticks and lines thicker
ax.TickLength = [.02,.05];
ax.LineWidth = 1;
% Make font size larger
ax.FontSize = 16;
%ax.TickLabelInterpreter = 'latex';