function labelAxes(ax, ff, theXLabel, theYLabel)

    ax.XLabel.String =  theXLabel;
    ax.XLabel.FontAngle = ff.axisFontAngle;

    ax.YLabel.String = theYLabel;
    ax.YLabel.FontAngle = ff.axisFontAngle;
end