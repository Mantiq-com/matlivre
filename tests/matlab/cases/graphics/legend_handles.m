% legend with a vector of plotted objects, legend hide, and a vector SizeData
f = figure('Visible', 'off'); x = 1:3;
p1 = plot(x, x); hold on; p2 = plot(x, 2*x); p3 = plot(x, 3*x);
l = legend([p2 p1], 'B', 'A');
fprintf('String: %s\n', strjoin(l.String, ','));
fprintf('DisplayName: %s | %s | [%s]\n', p1.DisplayName, p2.DisplayName, p3.DisplayName);
fprintf('children=%d\n', numel(f.Children));
l2 = legend([p3 p1]);
fprintf('no labels: %s\n', strjoin(l2.String, ','));
l3 = legend([p1 p2 p3], {'one', 'two', 'three'});
fprintf('cell labels: %s\n', strjoin(l3.String, ','));
fprintf('visible before hide: %s\n', char(l3.Visible));
legend hide
fprintf('after hide: children=%d visible=%s\n', numel(f.Children), char(l3.Visible));
legend show
fprintf('after show: children=%d visible=%s\n', numel(f.Children), char(l3.Visible));
legend('hide')
fprintf('function form: children=%d visible=%s\n', numel(f.Children), char(l3.Visible));
legend off
fprintf('after off: children=%d\n', numel(f.Children));
g = figure('Visible', 'off');
s1 = scatter(1:3, 1:3, [10 20 30]);
fprintf('SizeData=%s size=%s\n', mat2str(s1.SizeData), mat2str(size(s1.SizeData)));
s2 = scatter(1:3, 1:3, [10; 20; 30]);
fprintf('SizeData=%s size=%s\n', mat2str(s2.SizeData), mat2str(size(s2.SizeData)));
s3 = scatter(1:3, 1:3, 40);
fprintf('SizeData=%s\n', mat2str(s3.SizeData));
s3.SizeData = [1 2 3];
fprintf('set vector: SizeData=%s\n', mat2str(s3.SizeData));
s3.SizeData = 9;
fprintf('set scalar: SizeData=%s\n', mat2str(s3.SizeData));
close all
clear f x p1 p2 p3 l l2 l3 g s1 s2 s3
