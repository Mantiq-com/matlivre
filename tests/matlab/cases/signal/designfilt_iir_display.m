% designfilt (IIR) and the digitalFilter display. MATLAB lists the specification properties in a
% run-dependent order, so the display lines are sorted before they are printed.
d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',0.3);
d2 = designfilt('highpassiir','StopbandFrequency',0.2,'PassbandFrequency',0.3,'PassbandRipple',1,'StopbandAttenuation',40,'DesignMethod','ellip');
d3 = designfilt('bandpassiir','FilterOrder',6,'PassbandFrequency1',0.2,'PassbandFrequency2',0.5,'PassbandRipple',1);
ds = {d1, d2, d3};
for k = 1:numel(ds)
  txt = evalc('disp(ds{k})');
  ln = sort(strsplit(strtrim(txt), newline));
  fprintf('%s\n', ln{:});
  fprintf('%s %d %d\n', class(ds{k}), size(ds{k}.Coefficients));
end
clear d1 d2 d3 ds k txt ln
