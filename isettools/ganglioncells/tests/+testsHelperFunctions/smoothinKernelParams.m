function [smoothinKernelSizePixels, smoothinKernelSigmaPixels] = smoothinKernelParams()
    smoothinKernelSizePixels = [];
    smoothinKernelSigmaPixels = [];

    %smoothinKernelSizePixels = 9*2+1;
    %smoothinKernelSigmaPixels = 1.35*2;

    if (1==2)
        % Or specify a value for both to employ a user-specificed RF smoothing kernel
        smoothinKernelSizePixels = 11/3;
        smoothinKernelSigmaPixels = smoothinKernelSizePixels / 6;
        smoothinKernelSizePixels = round(smoothinKernelSizePixels);
        if (mod(smoothinKernelSizePixels,2) == 0)
            % Make it odd
            smoothinKernelSizePixels = smoothinKernelSizePixels + 1;
        end
    end
end