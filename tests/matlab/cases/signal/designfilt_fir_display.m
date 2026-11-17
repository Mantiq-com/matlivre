% designfilt (FIR) and the digitalFilter display. MATLAB lists the specification properties in a
% run-dependent order, so the display lines are sorted before they are printed.
d1 = designfilt('lowpassfir','FilterOrder',20,'PassbandFrequency',0.2,'StopbandFrequency',0.3,'DesignMethod','ls');
d2 = designfilt('lowpassfir','FilterOrder',20,'PassbandFrequency',0.2,'StopbandFrequency',0.3);
d3 = designfilt('highpassfir','FilterOrder',20,'CutoffFrequency',200,'SampleRate',1000,'DesignMethod','window','Window','hann');
d4 = designfilt('bandpassfir','FilterOrder',30,'StopbandFrequency1',0.1,'PassbandFrequency1',0.2,'PassbandFrequency2',0.5,'StopbandFrequency2',0.6);
ds = {d1, d2, d3, d4};
for k = 1:numel(ds)
  txt = evalc('disp(ds{k})');
  ln = sort(strsplit(strtrim(txt), newline));
  fprintf('%s\n', ln{:});
  fprintf('%s %d\n', class(ds{k}), numel(ds{k}.Coefficients));
end
clear d1 d2 d3 d4 ds k txt ln
