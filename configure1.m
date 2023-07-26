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

function obj = configure1(obj, optimalPath)
    sensorCombination = nchoosek(1:obj.gridWorld.nPoints,obj.nSensors);
    numberCombination = numel((sensorCombination(:,1)));
    mutualInformationHistory  = [];
    mutualInformation_array   = zeros(1,obj.nSensors + 1);                 %[index, MI]
%     allConf_MI = [];                                                       % Array to store MI for all configurations
%     j     = 1;
    pathLength = length(optimalPath);

    for i = 1:numberCombination
        possibleConfigurations	= sensorCombination(i,:);
        H     = obj.calc_rbf_value(possibleConfigurations);
        pNext = obj.threatModel.pNext;
        pNextHist = obj.threatModel.pNextHist;
        size(pNext);
        size(pNextHist);

%         pNextHist(:,:,1);
%         diag(pNextHist(:,:,1))
        sum  = zeros(9,1);
        sum1 = zeros(1,1);
        sum2 = zeros(1,1);
        
        for k    = 1:pathLength
            phi  = obj.calc_rbf_value(optimalPath(:,k));
            sum  = sum + (pNextHist(:,:,k)) * phi';
            sum1 = sum1 + phi *(pNextHist(:,:,k)) * phi';
            obj.truepathCost = obj.threatModel.offset + phi * obj.threatModel.state;
            obj.estimatedpathCost = obj.threatModel.offset + phi * obj.threatModel.stateEstimate;
            for ii = 1: obj.threatModel.pathLength-1
                for jj = 2: obj.threatModel.pathLength
                    phi1 = obj.calc_rbf_value(optimalPath(:,ii));
                    phi2 = obj.calc_rbf_value(optimalPath(:,jj));
                    sum2 = phi1 *(pNextHist(:,:,ii)) * (pNextHist(:,:,jj))*phi2';
                end
            end
        end
        tau = H * sum;
        tau = tau';
%         obj.calc_rbf_value(optimalPath(:,2));
        Xi  = sum1 + 2* sum2;
        obj.pathRisk = obj.estimatedpathCost + sqrt(Xi);
        
       % tau   = pNext * H';
%        Xi     = H * pNext*H' + obj.noiseVariance * eye(obj.nSensors);
        currentmutualInformation = 0.5 * log(det(pNext)/(det(pNext - tau * pinv(Xi) * tau')));
        mutualInformationHistory = cat(2, mutualInformationHistory,  currentmutualInformation);
%         allConf_MI(j,1) = i;
%         allConf_MI(j,2) = currentmutualInformation;

        if currentmutualInformation > mutualInformation_array(:,end)
           mutualInformation_array(:,1:obj.nSensors) = sensorCombination(i,:);
           mutualInformation_array(:,obj.nSensors +1) = currentmutualInformation;
        end
%         j = j + 1;       
    end
%     allConf_MI;
    obj.configuration = mutualInformation_array(:, 1:obj.nSensors);
    obj.configHistory = cat(2, obj.configHistory,  obj.configuration');
    obj.truepathCostHistory = [obj.truepathCostHistory		obj.truepathCost];
    obj.estimatedpathCostHistory = [obj.estimatedpathCostHistory		obj.estimatedpathCost];

end
