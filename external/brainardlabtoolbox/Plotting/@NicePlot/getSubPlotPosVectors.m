function posVectors = getSubPlotPosVectors(varargin)

    self.rowsNum        = 2;
    self.colsNum        = 2;
    self.widthMargin    = 0.01;
    self.heightMargin   = 0.01;
    self.leftMargin     = 0.06;
    self.rightMargin    = 0.06;
    self.bottomMargin   = 0.08;
    self.topMargin      = 0.08;
    
    % parse inputs
    parser = inputParser;
    parser.addParamValue('rowsNum',         self.rowsNum);
    parser.addParamValue('colsNum',         self.colsNum);
    parser.addParamValue('widthMargin',     self.widthMargin);
    parser.addParamValue('heightMargin',    self.heightMargin);
    parser.addParamValue('leftMargin',      self.leftMargin);
    parser.addParamValue('rightMargin',     self.rightMargin);
    parser.addParamValue('bottomMargin',    self.bottomMargin); 
    parser.addParamValue('topMargin',       self.topMargin);

    
    % Execute the parser to make sure input is good
    parser.parse(varargin{:});
    % Copy the parse parameters to the ExperimentController object
    pNames = fieldnames(parser.Results);
    for k = 1:length(pNames)
       self.(pNames{k}) = parser.Results.(pNames{k}); 
    end
    
    plotWidth  = ((1.0-self.leftMargin-self.rightMargin) - self.widthMargin*(self.colsNum-1) - 0.01)/self.colsNum;
    plotHeight = ((1.0-self.bottomMargin-self.topMargin) - self.heightMargin*(self.rowsNum-1) - 0.01)/self.rowsNum;
    
    for row = 1:self.rowsNum
        yo = 0.99 - self.topMargin - (row)*(plotHeight+self.heightMargin) + self.heightMargin;
        for col = 1:self.colsNum
            xo = self.leftMargin + (col-1)*(plotWidth+self.widthMargin);
            posVectors(row,col).v = [xo yo plotWidth plotHeight];
        end
    end
end
