%
%PAL_PFLR_setupMC       Create or update structure MC specifying which
%   model comparison is to be performed by PAL_PFLR_TLR.
%
%Syntax: MC = PAL_PFLR_setupMC(params, {optional arguments})
%
%Internal function
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.6.3 (see History.m)

function MC = PAL_PFLR_setupMC(params, varargin)

NumOpts = length(varargin);

if mod(NumOpts,2) == 0
    MC.argsAfuller = PAL_Contrasts(size(params,1),'identity');
    MC.argsBfuller = PAL_Contrasts(size(params,1),'identity');
    MC.argsGfuller = [];
    MC.argsLfuller = [];
    MC.argsAlesser = ones(1,size(params,1));
    MC.argsBlesser = ones(1,size(params,1));
    MC.argsGlesser = [];
    MC.argsLlesser = [];
else
    MC = varargin{1};
end

if NumOpts > 1
    for n = [1:2:NumOpts-mod(NumOpts,2)]
        n = n+mod(NumOpts,2);
        valid = 0;
        if strncmpi(varargin{n}, 'ThresholdLesser', 10)
            valid = 1;
            if isnumeric(varargin{n+1})
                MC.argsAlesser = varargin{n+1};
            else
                if strncmpi(varargin{n+1}, 'constrained', 4)
                    MC.argsAlesser = ones(1,size(params,1));
                end
                if strncmpi(varargin{n+1}, 'unconstrained', 4)
                    MC.argsAlesser = PAL_Contrasts(size(params,1),'identity');
                end
                if strcmpi(varargin{n+1}, 'fixed')
                    MC.argsAlesser = [];
                end             
            end
        end
        if strncmpi(varargin{n}, 'SlopeLesser', 8)
            valid = 1;
            if isnumeric(varargin{n+1})
                MC.argsBlesser = varargin{n+1};
            else
                if strncmpi(varargin{n+1}, 'constrained', 4)
                    MC.argsBlesser = ones(1,size(params,1));
                end
                if strncmpi(varargin{n+1}, 'unconstrained', 4)
                    MC.argsBlesser = PAL_Contrasts(size(params,1),'identity');
                end
                if strcmpi(varargin{n+1}, 'fixed')
                    MC.argsBlesser = [];
                end             
            end
        end
        if strncmpi(varargin{n}, 'guessLesser', 8)
            valid = 1;
            if isnumeric(varargin{n+1})
                MC.argsGlesser = varargin{n+1};
            else
                if strncmpi(varargin{n+1}, 'constrained', 4)
                    MC.argsGlesser = ones(1,size(params,1));
                end
                if strncmpi(varargin{n+1}, 'unconstrained', 4)
                    MC.argsGlesser = PAL_Contrasts(size(params,1),'identity');
                end
                if strcmpi(varargin{n+1}, 'fixed')
                    MC.argsGlesser = [];
                end             
            end
        end
        if strncmpi(varargin{n}, 'lapseLesser', 8)
            valid = 1;
            if isnumeric(varargin{n+1})
                MC.argsLlesser = varargin{n+1};
            else
                if strncmpi(varargin{n+1}, 'constrained', 4)
                    MC.argsLlesser = ones(1,size(params,1));
                end
                if strncmpi(varargin{n+1}, 'unconstrained', 4)
                    MC.argsLlesser = PAL_Contrasts(size(params,1),'identity');
                end
                if strcmpi(varargin{n+1}, 'fixed')
                    MC.argsLlesser = [];
                end             
            end
        end
        if strncmpi(varargin{n}, 'ThresholdFuller', 12)
            valid = 1;
            if isnumeric(varargin{n+1})
                MC.argsAfuller = varargin{n+1};
            else
                if strncmpi(varargin{n+1}, 'constrained', 4)
                    MC.argsAfuller = ones(1,size(params,1));
                end
                if strncmpi(varargin{n+1}, 'unconstrained', 4)
                    MC.argsAfuller = PAL_Contrasts(size(params,1),'identity');
                end
                if strcmpi(varargin{n+1}, 'fixed')
                    MC.argsAfuller = [];
                end             
            end
        end
        if strncmpi(varargin{n}, 'Slopefuller', 8)
            valid = 1;
            if isnumeric(varargin{n+1})
                MC.argsBfuller = varargin{n+1};
            else
                if strncmpi(varargin{n+1}, 'constrained', 4)
                    MC.argsBfuller = ones(1,size(params,1));
                end
                if strncmpi(varargin{n+1}, 'unconstrained', 4)
                    MC.argsBfuller = PAL_Contrasts(size(params,1),'identity');
                end
                if strcmpi(varargin{n+1}, 'fixed')
                    MC.argsBfuller = [];
                end             
            end
        end
        if strncmpi(varargin{n}, 'guessfuller', 8)
            valid = 1;
            if isnumeric(varargin{n+1})
                MC.argsGfuller = varargin{n+1};
            else
                if strncmpi(varargin{n+1}, 'constrained', 4)
                    MC.argsGfuller = ones(1,size(params,1));
                end
                if strncmpi(varargin{n+1}, 'unconstrained', 4)
                    MC.argsGfuller = PAL_Contrasts(size(params,1),'identity');
                end
                if strcmpi(varargin{n+1}, 'fixed')
                    MC.argsGfuller = [];
                end             
            end
        end
        if strncmpi(varargin{n}, 'lapsefuller', 8)
            valid = 1;
            if isnumeric(varargin{n+1})
                MC.argsLfuller = varargin{n+1};
            else
                if strncmpi(varargin{n+1}, 'constrained', 4)
                    MC.argsLfuller = ones(1,size(params,1));
                end
                if strncmpi(varargin{n+1}, 'unconstrained', 4)
                    MC.argsLfuller = PAL_Contrasts(size(params,1),'identity');
                end
                if strcmpi(varargin{n+1}, 'fixed')
                    MC.argsLfuller = [];
                end             
            end
        end

        if valid == 0
            warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{n})
        end        

    end
end