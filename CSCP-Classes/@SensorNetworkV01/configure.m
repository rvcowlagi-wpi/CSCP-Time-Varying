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

function obj = configure(obj, optimalPath)
    numberCombination = nchoosek(obj.gridWorld.nPoints,obj.nSensors);
    sensorCombination = nchoosek(1:obj.gridWorld.nPoints,obj.nSensors);
    mutualInformationHistory = [];
    for   i     = 1:numberCombination
          obj.allConfigurations	= sensorCombination(i,:);
          H     = obj.calc_rbf_value(obj. allConfigurations);
          pNext = obj.threatModel.pNext;
          tau   = pNext * H';
          E     = H * pNext*H' + obj.noiseVariance * eye(obj.nSensors);
          mutualInformation = 0.5 * log(det(pNext)/(det(pNext - tau * pinv(E) * tau')));
          mutualInformationHistory = cat(2, mutualInformationHistory,  mutualInformation);
    end
%         MI_arr(j,1) = i;
%         MI_arr(j,2) = curr_mutualInformation;
%         threatStateHat_k = threat_.stateEstimate;
%     
%         if curr_mutualInformation>max_mutualInf(:,1)
%             max_mutualInf = [curr_mutualInformation, i];
%         end
%         j = j+1;
end

