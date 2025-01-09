% Firfilt_raw.m: pads constants at the begging and at the end, amount: same as group delay, which is the delay by which the FIR filter is shifted
%
% Also it uses zi parameter in the MATLAB filter function to do continuous 
% filtering, it does three separate filters. First filter filters the front 
% padding and equivalent length worth of data, the third filter filters the 
% back padding and equivalent length worth of data, and the second filter 
% filters the rest. The reason the filtering process is separated to three 
% separate steps is because the filter funcion in MATLAB has a parameter 
% called zi, which can be passed on to subsequent filter operations to 
% ensure two filter function's intial conditions are continuous as oppose 
% to being reinitilaised everytime.  
% 
% 
% 
% 
% 
% 
% 
% 
% firfilt() - Pad data with DC constant, filter data with FIR filter,
%             and shift data by the filter's group delay
%
% Usage:
%   >> EEG = firfilt(EEG, b, nFrames);
%
% Inputs:
%   EEG           - EEGLAB EEG structure
%   b             - vector of filter coefficients
%
% Optional inputs:
%   nFrames       - number of frames to filter per block {default 1000}
%
% Outputs:
%   EEG           - EEGLAB EEG structure
%
% Note:
%   Higher values for nFrames increase speed and working memory
%   requirements.
%
% Author: Andreas Widmann, University of Leipzig, 2005
%
% See also:
%   filter, findboundaries

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function data = firfilt_raw(data, b, nFrames)

if nargin < 2
    error('Not enough input arguments.');
end
if nargin < 3 || isempty(nFrames)
    nFrames = 1000;
end

% Filter's group delay
if mod(length(b), 2) ~= 1
    error('Filter order is not even.');
end
groupDelay = (length(b) - 1) / 2;

dcArray = [1,length(data)];
% Initialize progress indicator
nSteps = 20;
step = 0;
fprintf(1, 'firfilt(): filtering...');
% fprintf(1, 'firfilt(): |');
% strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '|   0%%']);
% tic

for iDc = 1:(length(dcArray) - 1)

        % Pad beginning of data with DC constant and get initial conditions
        ziDataDur = min(groupDelay, dcArray(iDc + 1) - dcArray(iDc));
        [~, zi] = filter(b, 1, double([data(:, ones(1, groupDelay) * dcArray(iDc)) ...
                                  data(:, dcArray(iDc):(dcArray(iDc) + ziDataDur - 1))]), [], 2);

        blockArray = [(dcArray(iDc) + groupDelay):nFrames:(dcArray(iDc + 1) - 1) dcArray(iDc + 1)];
        %disp(blockArray)
        for iBlock = 1:(length(blockArray) - 1)

            % Filter the data
            [data(:, (blockArray(iBlock) - groupDelay):(blockArray(iBlock + 1) - groupDelay - 1)), zi] = ...
                filter(b, 1, double(data(:, blockArray(iBlock):(blockArray(iBlock + 1) - 1))), zi, 2);

            % Update progress indicator
            % [step, strLength] = mywaitbar((blockArray(iBlock + 1) - groupDelay - 1), size(EEG.data, 2), step, nSteps, strLength);
        end

        % Pad end of data with DC constant
        temp = filter(b, 1, double(data(:, ones(1, groupDelay) * (dcArray(iDc + 1) - 1))), zi, 2);
        data(:, (dcArray(iDc + 1) - ziDataDur):(dcArray(iDc + 1) - 1)) = ...
            temp(:, (end - ziDataDur + 1):end);

        % Update progress indicator
        % [step, strLength] = mywaitbar((dcArray(iDc + 1) - 1), size(EEG.data, 2), step, nSteps, strLength);

end

% Deinitialize progress indicator
fprintf(1, '\n')

end