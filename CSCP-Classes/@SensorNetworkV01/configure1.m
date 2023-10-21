%{
SOFTWARE LICENSE
----------------
Copyright (c) 2023 by 
	Raghvendra V Cowlagi
	Bejamin Cooper
	Prakash Poudel

Permission is hereby granted to any person obtaining a copy of this
software and associated documentation files (the "Software"), to deal in
the Software, including the rights to use, copy, modify, merge, copies of
the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:  

* The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.
* The Software, and its copies or modifications, may not be distributed,
published, or sold for profit. 
* The Software, and any substantial portion thereof, may not be copied or
modified for commercial or for-profit use.

The software is provided "as is", without warranty of any kind, express or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose and noninfringement. In no event shall the
authors or copyright holders be liable for any claim, damages or other
liability, whether in an action of contract, tort or otherwise, arising
from, out of or in connection with the software or the use or other
dealings in the software.      


PROGRAM DESCRIPTION
-------------------
Reconfiguration function in class definition of sensor network:
	* Grid locations
	* Placement with basis ID matching
	* Works with a threat model defined by ParametricThreat class
%}

function obj = configure1(obj,threat_, grid_, optimalPath, timestep)
    sensorCombination = nchoosek(1:obj.gridWorld.nPoints,obj.nSensors);
    numberCombination = numel((sensorCombination(:,1)));
    mutualInformationHistory  = [];
    mutualInformation_array   = zeros(1,obj.nSensors + 1);                 %[index, MI]
%     allConf_MI = [];                                                       % Array to store MI for all configurations
%     j     = 1;
    pathLength = length(optimalPath);

    for i = 1:numberCombination
        possibleConfigurations	= sensorCombination(i,:);
        H = threat_.calc_rbf_value(grid_.coordinates(:, possibleConfigurations));
        pNext = obj.threatModel.pNext;
        pNextHist = obj.threatModel.pNextHist;
        size(pNext);
        size(pNextHist);

        sum  = zeros(threat_.nStates,1);
        sum1 = zeros(1,1);
        sum2 = zeros(1,1);
        obj.truepathCost = 0;
        obj.estimatedpathCost = 0;
        for k    = 1:pathLength
            phi  = threat_.calc_rbf_value(grid_.coordinates(:, optimalPath(:,k)));
            sum  = sum + (pNextHist(:,:,k)) * phi';
            sum1 = sum1 + phi *(pNextHist(:,:,k)) * phi';
            obj.truepathCost = obj.truepathCost  + phi * obj.threatModel.state;
            obj.estimatedpathCost = obj.estimatedpathCost + phi * obj.threatModel.stateEstimate;

            for ii = 1: obj.pathLength-1
                for jj = 2: obj.pathLength
                    phi1 = threat_.calc_rbf_value(grid_.coordinates(:, optimalPath(:,ii)));
                    phi2 = threat_.calc_rbf_value(grid_.coordinates(:, optimalPath(:,jj)));
                    sum2 = phi1 *((pNextHist(:,:,ii)) + (pNextHist(:,:,jj)))*phi2';
                end
            end
        end
        
        tau = H * obj.gridWorld.spacing * sum;
        tau = tau';
        obj.varpathCost  = (obj.gridWorld.spacing)^2 * (sum1 + 2* sum2);
        Xi     = H * pNext*H' + obj.noiseVariance * eye(obj.nSensors);
        obj.truepathCost =  pathLength + obj.gridWorld.spacing  *obj.truepathCost; 
        obj.estimatedpathCost =  pathLength + obj.gridWorld.spacing  *obj.estimatedpathCost;  
        obj.pathRisk = obj.estimatedpathCost + sqrt(obj.varpathCost);
        
       % Calculate mutual information between the cost and measurement.
        currentmutualInformation = 0.5 * log(det(obj.varpathCost)/(det(obj.varpathCost - tau * pinv(Xi) * tau')));
        mutualInformationHistory = cat(2, mutualInformationHistory,  currentmutualInformation);
        if currentmutualInformation > mutualInformation_array(:,end)
           mutualInformation_array(:,1:obj.nSensors) = sensorCombination(i,:);
           mutualInformation_array(:,obj.nSensors +1) = currentmutualInformation;
        end      
    end
    obj.configuration = mutualInformation_array(:, 1:obj.nSensors);
    obj.configHistory = cat(2, obj.configHistory,  obj.configuration');
    obj.truepathCostHistory = [obj.truepathCostHistory		obj.truepathCost];
    obj.estimatedpathCostHistory = [obj.estimatedpathCostHistory		obj.estimatedpathCost];
    obj.varpathCostHistory = [obj.varpathCostHistory		obj.varpathCost];
    obj.pathRiskHistory = [obj.pathRiskHistory		obj.pathRisk];
end
