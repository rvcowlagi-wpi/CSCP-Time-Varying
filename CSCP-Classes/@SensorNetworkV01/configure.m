PROGRAM DESCRIPTION
-------------------
Reconfiguration function in class definition of sensor network:
	* Grid locations
	* Placement with basis ID matching
	* Works with a threat model defined by ParametricThreat class
%}

function obj = configure(obj,threat_, grid_, optimalPath, timestep)
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
        sum1 = zeros(1,1);
        sum2 = zeros(1,1);
        obj.truepathCost = 0;
        obj.estimatedpathCost = 0;
        for k    = 1:pathLength
            phi  = threat_.calc_rbf_value(grid_.coordinates(:, optimalPath(:,k)));
            sum1 = sum1 + phi *(pNextHist(:,:,k)) * phi';
            obj.truepathCost = obj.truepathCost + phi * obj.threatModel.state;
            obj.estimatedpathCost = obj.estimatedpathCost + phi * obj.threatModel.stateEstimate;
            for ii = 1: obj.pathLength-1
                for jj = 2: obj.pathLength
                    phi1 = threat_.calc_rbf_value(grid_.coordinates(:, optimalPath(:,ii)));
                    phi2 = threat_.calc_rbf_value(grid_.coordinates(:, optimalPath(:,jj)));
                    sum2 = phi1 *(pNextHist(:,:,ii)) * (pNextHist(:,:,jj))*phi2';
                end
            end
        end
        
        tau   = pNext * H';
        obj.varpathCost  = (obj.gridWorld.spacing)^2 * (sum1 + 2* sum2);
        Xi     = H * pNext*H' + obj.noiseVariance * eye(obj.nSensors);
        obj.truepathCost =   pathLength + obj.gridWorld.spacing  *obj.truepathCost; 
        obj.estimatedpathCost =  pathLength + obj.gridWorld.spacing  *obj.estimatedpathCost;  
        obj.pathRisk = obj.estimatedpathCost + sqrt(obj.varpathCost);

        % Calculate mutual information between the state and measurement.
        currentmutualInformation = 0.5 * log(det(pNext)/(det(pNext - tau * pinv(Xi) * tau')));
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
