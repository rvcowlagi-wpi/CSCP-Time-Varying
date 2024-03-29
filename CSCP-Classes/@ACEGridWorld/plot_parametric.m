function plot_parametric(obj, threat_, sensor_, flags_)

if flags_.SHOW_TRUE && flags_.SHOW_ESTIMATE
	if flags_.JUXTAPOSE
		figure('Name', 'True', 'Units','normalized', 'Position', [0.3 0.1 0.3*[1.8 1.6]]);
		axisTrue= subplot(1,2,1);
		axisEst	= subplot(1,2,2);
	else
		figure('Name', 'True', 'Units','normalized', 'Position', [0.7 0.1 0.25*[0.9 1.6]]);
		axisTrue= gca;
		figure('Name', 'Estimate', 'Units','normalized', 'Position', [0.7 0.4 0.25*[0.9 1.6]]);
		axisEst	= gca;
	end
else
	figure('Name', 'True', 'Units','normalized', 'Position', [0.6 0.1 0.25*[0.9 1.6]]);
	if flags_.SHOW_TRUE
		axisTrue= gca;
	elseif flags_.SHOW_ESTIMATE
		axisEst	= gca;
	end
end

nPlotPts		= 200;
xGridPlot		= linspace(-obj.halfWorkspaceSize, obj.halfWorkspaceSize, nPlotPts);
yGridPlot		= linspace(-obj.halfWorkspaceSize, obj.halfWorkspaceSize, nPlotPts);
[xMesh, yMesh]	= meshgrid(xGridPlot, yGridPlot);
locationsMesh(:, :, 1) = xMesh;
locationsMesh(:, :, 2) = yMesh;

if flags_.SHOW_TRUE
	threatMesh	= threat_.calculate_at_locations(...
		locationsMesh, threat_.stateHistory(:, 1));
	imageMax	= max(threatMesh(:));
	imageMin	= min(threatMesh(:));
	imageClims	= [0.8*imageMin 1.5*imageMax];
	
	grHdlSurf	= surfc(axisTrue, xMesh, yMesh, threatMesh,'LineStyle','none');
	clim(imageClims); colorbar; view(2);
	axis equal; axis tight; hold on;
	set(gca, 'Color', '#E0E0E0')

	xlim(1.05*[-obj.halfWorkspaceSize, obj.halfWorkspaceSize]); 
	ylim(1.05*[-obj.halfWorkspaceSize, 1.45*obj.halfWorkspaceSize]);
	zlim(imageClims);
	
	timeText = ['$t = $ ' num2str(0) ' units'];
	grHdlTimeText	= text(axisTrue, ...
		-0.98*obj.halfWorkspaceSize, 1.3*obj.halfWorkspaceSize, 2*imageMax, timeText, ...
		'Color', 'k', 'FontName', 'Times New Roman', ...
		'FontSize', 12, 'Interpreter','latex');

	%----- Plot grid
	plot3(...
		obj.coordinates(1, :), obj.coordinates(2, :), ...
		imageMax*ones(1, size(obj.coordinates, 2)), ...
		'.', 'Color', 'w', 'MarkerSize', 20);

	%----- Plot centers of basis functions
	plot3(...
		threat_.basisCenter(1, :), threat_.basisCenter(2, :), ...
		imageMax*ones(1, size(threat_.basisCenter, 2)), ...
		'.', 'Color', 'k', 'MarkerSize', 30);
	for m2 = 1:threat_.nStates
		text(axisTrue, ...
			threat_.basisCenter(1, m2), (threat_.basisCenter(2, m2) + 0.05), ...
			2*imageMax, num2str(m2), 'Color', 'k', 'FontName', 'Times New Roman', ...
			'FontSize', 12, 'Interpreter','latex')
	end

	drawnow();

	for m1 = 2:length(threat_.timeStampState)
		delete(grHdlSurf);
		delete(grHdlTimeText);

		%----- Plot threat field as a surface
		threatMesh	= threat_.calculate_at_locations(...
			locationsMesh, threat_.stateHistory(:, m1));
		surfc(axisTrue, xMesh, yMesh, threatMesh,'LineStyle','none');
		hold on;



		%----- Indicate time step
		timeText = ['$t = $ ' num2str(threat_.timeStampState(m1)) ' units'];
		grHdlTimeText	= text(axisTrue, ...
			-0.98*obj.halfWorkspaceSize, 1.3*obj.halfWorkspaceSize, 2*imageMax, timeText, ...
			'Color', 'k', 'FontName', 'Times New Roman', ...
			'FontSize', 12, 'Interpreter','latex');

		drawnow();
	end
end