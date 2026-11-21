% Axes with nothing in them have the limits [0 1]; a 2-D plot has z limits [-1 1]
ax = axes;
xl = ax.XLim;
yl = ax.YLim;
zl = ax.ZLim;
zt = ax.ZTick;
plot(1:3);
ax2 = gca;
zl2 = ax2.ZLim;
zt2 = ax2.ZTick;
clear ax ax2
