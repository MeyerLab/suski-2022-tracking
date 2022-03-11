function [inds] = gate_vals(data, gates, intervals, checkField)
%GATE_VALS gate indices based on parameters
%
%   INPUTS: 
%   DATA      - Data structure containing data to be gated on
%   GATES     - Names of structure fields to be gated on, in a cell array
%   INTERVALS - Corresponding intervals for each field to be gated on (in
%               cell array, e.g. {[0 1],[10 50]}
%   OUTPUTS:
%   INDS      - Outputs indices which satisfy all gates

    if nargin < 4
        checkField = 'area';
    end
    inds = true(size(data.(checkField),1),1);
    for i = 1:length(gates)
        gateData = data.(gates{i});
        inds = inds & gateData < intervals{i}(2) & gateData > intervals{i}(1);
    end

end

