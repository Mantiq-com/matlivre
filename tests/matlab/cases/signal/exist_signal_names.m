% exist() sees the signal functions (MATLAB returns 2, MatLivre 5 for built-ins: only nonzero is pinned) and they are not shadowed by variables.
r = double(all([exist('window') exist('db') exist('thd') exist('sfdr') exist('toi') exist('snr') exist('sinad') exist('fillgaps') exist('stftmag2sig') exist('scaleFilterSections')] > 0))
db = 3;
q = db(1)
clear db
