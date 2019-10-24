% eeg_eventformat() - Convert the event information of a dataset from struct
%                 to array or vice versa.
%
% Usage: >> [eventout fields] = eeg_eventformat( event, 'format', fields );
%
% Inputs:
%   event  - event array or structure
%   format - ['struct'|'array'] see below
%   fields - [optional] cell array of strings containing the names of
%            the event struct fields. If this field is empty, it uses 
%            the following list for
%            the names of the fields { 'type' 'latency' 'var1' ...
%            'var2' ... }.
% Output:
%   eventout  - output event array or structure
%   fields    - output cell array with the name of the fields
%
% Event formats:
%   struct - Events are organised as an array of structs with at
%            least two fields ('type' and 'latency')
%            (Ex: reaction_time may be type 1).
%   array  - events are organized as an array, the first column
%            representing the type, the second the latency and the
%            other ones user-defined variables.
%
% Note: 1) The event structure is defined only for continuous data
%          or epoched data derived from continuous data.
%       2) The event 'struct' format is more comprehensible.
%          For instance, to see all the properties of event 7,
%          type >> EEG.event(7)
%          Unfortunately, structures are awkward for expert users to deal
%          with from the command line (Ex: To get an array of latencies,
%           >> [ EEG.event(:).latency ] )
%          In array format, the same information is obtained by typing
%           >> EEG.event(:,2)
%       3) This function automatically updates the 'eventfield'
%          cell array depending on the format.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 27 Jan 2002
%
% See also: eeglab(), pop_selectevent(), pop_importevent()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 27 Jan 2002, arno@salk.edu
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

% 2/06/02 modifed header - sm & ad
% 2/08/02 add field input - ad
% 2/12/02 reprogrammed function using epochformat.m - ad

function [event, eventfield] = eeg_eventformat(event, format, fields);

if nargin < 2
   help eeg_eventformat;
   return;
end;	

if exist('fields') ~= 1, fields = { 'type', 'latency' }; end;

[event eventfield] = eeg_epochformat( event, format, fields);

