steps = [1 2 3];
modelname = 'TWG';

% Mesh
meshsize = 200;

% Friction setup
friction_coefficient_min = 20;
friction_law = 'Weertman'; % 'Budd' | 'Weertman'
% Specify exponents using the following convention:
%     taub = C_Budd * N * ub^m
%     taub = C_Weertman * ub^m
m_Budd = 1;     % 1 | 1/3   [ISSM convention is p=1/m, q=p, s=1/p (s=m), r=q/p (r=1)]
m_Weertman = 1; % 1 | 1/3   [ISSM convention is m=1/m]

% Terminus setup
terminus0_x = 26000;
magnitude = 1250;

% Cluster parameters
cluster = generic('name', oshostname(), 'np', 2);
cluster.interactive = 1;
waitonlock = 10;

% Run steps
org=organizer('repository',['./Models_' modelname],'prefix',['MISMIP_' modelname '_' num2str(meshsize) 'm_'],'steps',steps);

% Initialization
if perform(org,'Mesh'),% {{{1  STEP 1

   md = bamg(model,'domain',['Exp/Domain_' modelname '.exp'],'hmax',meshsize,'splitcorners',1);
	md.miscellaneous.name=['MISMIP_' modelname];

	savemodel(org,md);
end % }}}
if perform(org, ['Parameterization_frictionCoeffMin' num2str(friction_coefficient_min)]),% {{{1  STEP 2

	md=loadmodel(org,'Mesh');

	md=setmask(md,'','');
	md=parameterize(md, ['MISMIP_' modelname '.par']);

	savemodel(org,md);
end% }}}
if perform(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min)]),% {{{1  STEP 3

	md=loadmodel(org, ['Parameterization_frictionCoeffMin' num2str(friction_coefficient_min)]);

	md=setflowequation(md,'SSA','all');

   md.transient.requested_outputs={'default','IceVolumeAboveFloatation'};

	md.timestepping.time_step=1;
	md.timestepping.final_time=100; %50000; %200000;
	
   if meshsize == 100
      md.timestepping.time_step = 0.025;
      md.settings.output_frequency=1000;
   elseif meshsize == 200
      md.timestepping.time_step = 0.050;
      md.settings.output_frequency=500;
   elseif meshsize == 400
      md.timestepping.time_step = 0.100;
      md.settings.output_frequency=500;
   elseif meshsize == 1000
      md.timestepping.time_step = 0.025;
      md.settings.output_frequency=500;
   else
      fprintf('Enter time_step into runme!\n');
      return
   end
	md.timestepping.final_time=50;
	md.settings.output_frequency=10;
   
   % Adaptive timestepping
   %md.timestepping = timesteppingadaptive(md.timestepping);
	%md.timestepping.time_step_min = 0.05;

	md.stressbalance.abstol=NaN;
	md.stressbalance.restol=1e-4;
   md.stressbalance.maxiter=10;
	md.verbose=verbose('convergence',false,'solution',true);
	md.cluster=cluster;
   md.settings.waitonlock=waitonlock;
	md=solve(md,'tr');
   if md.settings.waitonlock == 0
      fprintf('\n \033[103;30m Load results with: md = loadresultsfromcluster(md,''runtimename'',''%s''); \033[0m \n', md.private.runtimename);
      fprintf(' \033[103;30m Save results with: save(''%s.mat'', ''md'', ''-v7.3''); \033[0m \n\n', [org.repository filesep org.prefix org.steps(end).string]);
      return
   end

	savemodel(org,md);
end% }}}
if perform(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_noFloating']),% {{{1  STEP 4

   md = loadmodel(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min)]);

   % Transfer resuts to model fields
   md = transientrestart(md);
   md.timestepping.final_time = md.timestepping.start_time + 100;
   if meshsize == 100
      md.timestepping.time_step = 0.005;
      md.settings.output_frequency = 100;
   elseif meshsize == 200
      md.timestepping.time_step = 0.010;
      md.settings.output_frequency = 50;
   elseif meshsize == 400
      md.timestepping.time_step = 0.020;
      md.settings.output_frequency = 25;
   elseif meshsize == 1000
      md.timestepping.time_step = 0.050;
      md.settings.output_frequency = 10;
   else
      fprintf('Enter time_step into runme!\n');
      return
   end

   % Activate moving boundary
   md.transient.ismovingfront = 1;

   % Start spclevelset with the original ice_levelset
   md.levelset.spclevelset          = [md.mask.ice_levelset; md.timestepping.start_time];

   % Remove all floating ice over the first year
   levelset = md.mask.groundedice_levelset;
   levelset(md.mask.groundedice_levelset<0) = +1;
   levelset(md.mask.groundedice_levelset>0) = -1;
   md.levelset.spclevelset(:,end+1) = [levelset; md.timestepping.start_time+1];
   md.calving.calvingrate = zeros(md.mesh.numberofvertices,1);
   md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);

   % Convert to signed distance fields
   if size(md.levelset.spclevelset,2)>1,
      disp('Converting levelsets to signed distance fields');
      for i=1:size(md.levelset.spclevelset,2)
         levelset = md.levelset.spclevelset(1:end-1,i);
         pos      = find(levelset<0);

         if exist('TEMP.exp','file'), delete('TEMP.exp'); end
         expcontourlevelzero(md,levelset,0,'TEMP.exp');
         levelset = abs(ExpToLevelSet(md.mesh.x,md.mesh.y,'TEMP.exp'));
         delete('TEMP.exp');
         levelset(pos) = -levelset(pos);
         md.levelset.spclevelset(1:end-1,i) = levelset;
      end
   end

   % Set the requested outputs
   md.transient.requested_outputs={'default','IceVolumeAboveFloatation'};

   % Go solve
   md.verbose.solution=1;
   %md.verbose = verbose('all'); md.verbose.solver = true;
   md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   md = solve(md,'tr');
   if md.settings.waitonlock == 0
      fprintf('\n \033[103;30m Load results with: md = loadresultsfromcluster(md,''runtimename'',''%s''); \033[0m \n', md.private.runtimename);
      fprintf(' \033[103;30m Save results with: save(''%s.mat'', ''md'', ''-v7.3''); \033[0m \n\n', [org.repository filesep org.prefix org.steps(end).string]);
      return
   end

   % Save
   savemodel(org,md);

end % }}}
if perform(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_noFloating_extra']),% {{{1  STEP 5

   md = loadmodel(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_noFloating']);

   % Transfer resuts to model fields
   md = transientrestart(md);
   md.timestepping.final_time = md.timestepping.start_time + 100;

	if meshsize == 100 & contains(cluster.name, 'discover'), cluster.time = 10 * 3600; end

   % Go solve
   md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   md = solve(md,'tr');
   if md.settings.waitonlock == 0
      fprintf('\n \033[103;30m Load results with: md = loadresultsfromcluster(md,''runtimename'',''%s''); \033[0m \n', md.private.runtimename);
      fprintf(' \033[103;30m Run the following:\n');
      fprintf(' \033[103;30m  tmpResults = [md.results.TransientSolution3 md.results.TransientSolution];\n');
      fprintf(' \033[103;30m  md.results.TransientSolution = tmpResults;\n');
      fprintf(' \033[103;30m  md.results = rmfield(md.results, ''TransientSolution3'');\n');
      fprintf(' \033[103;30m  clear tmpResults\n');
      fprintf(' \033[103;30m Save results with: save(''%s.mat'', ''md'', ''-v7.3''); \033[0m \n\n', [org.repository filesep org.prefix ...
         strrep(org.steps(end).string, '_extra', '')]);
      return
   end

   % Combine with Transient_Steadystate_noFloating results
   tmpResults = [md.results.TransientSolution3 md.results.TransientSolution];
   md.results.TransientSolution = tmpResults;
   md.results = rmfield(md.results, 'TransientSolution3');
   clear tmpResults

   % Save - cheat and save it back to Transient_Steadystate_noFloating
   save([org.repository filesep org.prefix strrep(org.steps(end).string, '_extra', '') '.mat'], 'md', '-v7.3');
   %savemodel(org,md);

end % }}}

%  Steady-state
if perform(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm']),%    STEP 6{{{

   % Start with Transient_Steadystate_noFloating because this gets us close to the terminus position
   % that we want to start with
   md = loadmodel(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_noFloating']);
   %md = loadmodel(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min)]);

   % Transfer results to model fields
   md = transientrestart(md);
   md.timestepping.final_time = md.timestepping.start_time + 30;
   if meshsize == 100
      md.timestepping.time_step = 0.005;
      md.settings.output_frequency = 100;
   elseif meshsize == 200
      md.timestepping.time_step = 0.010;
      md.settings.output_frequency = 50;
   elseif meshsize == 400
      md.timestepping.time_step = 0.020;
      md.settings.output_frequency = 25;
   elseif meshsize == 1000
      md.timestepping.time_step = 0.050;
      md.settings.output_frequency = 10;
   else
      fprintf('Enter time_step into runme!\n');
      return
   end

   % Set the terminus position
   md.transient.ismovingfront = 0;
   distance_x = md.mesh.x - terminus0_x;
   if exist('./TEMP.exp','file'), delete('./TEMP.exp'); end
   expcontourlevelzero(md,distance_x,0,'./TEMP.exp');
   levelset = abs(ExpToLevelSet(md.mesh.x,md.mesh.y,'./TEMP.exp'));
   delete('./TEMP.exp');
   md.mask.ice_levelset(md.mesh.x >terminus0_x) = +abs(levelset(md.mesh.x >terminus0_x));
   md.mask.ice_levelset(md.mesh.x<=terminus0_x) = -abs(levelset(md.mesh.x<=terminus0_x));
   %md.mask.ice_levelset(md.mesh.x >terminus0_x) = +1;
   %md.mask.ice_levelset(md.mesh.x<=terminus0_x) = -1;

   % Convert to signed distance fields -- NOT NEEDED ANYMORE BECAUSE THE CODE ABOVE NOW SETS UP THE LEVELSET AS A SIGNED DISTANCE
   %disp('Converting levelsets to signed distance fields');
   %pos      = find(md.mask.ice_levelset<0);
   %if exist('./TEMP.exp','file'), delete('./TEMP.exp'); end
   %expcontourlevelzero(md,md.mask.ice_levelset,0,'./TEMP.exp');
   %levelset = abs(ExpToLevelSet(md.mesh.x,md.mesh.y,'./TEMP.exp'));
   %delete('./TEMP.exp');
   %levelset(pos) = -levelset(pos);
   %md.mask.ice_levelset = levelset;

   % Go solve
   md.verbose.solution=1;
   md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   md = solve(md,'tr');
   if md.settings.waitonlock == 0
      fprintf('\n \033[103;30m Load results with: md = loadresultsfromcluster(md,''runtimename'',''%s''); \033[0m \n', md.private.runtimename);
      fprintf(' \033[103;30m Save results with: save(''%s.mat'', ''md'', ''-v7.3''); \033[0m \n\n', [org.repository filesep org.prefix org.steps(end).string]);
      return
   end

   % Save model
   savemodel(org,md);

end % }}}
if perform(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_extra']),%    STEP 7{{{

   md = loadmodel([org.repository '/' org.prefix 'Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm.mat']);

   % Transfer results to model fields
   md = transientrestart(md);
   md.timestepping.final_time = md.timestepping.start_time + 50;
   %if meshsize == 100
   %   md.timestepping.time_step = 0.005;
   %   md.settings.output_frequency = 100;
   %elseif meshsize == 200
   %   md.timestepping.time_step = 0.010;
   %   md.settings.output_frequency = 50;
   %elseif meshsize == 400
   %   md.timestepping.time_step = 0.020;
   %   md.settings.output_frequency = 25;
   %elseif meshsize == 1000
   %   md.timestepping.time_step = 0.250;
   %   md.settings.output_frequency = 10;
   %else
   %   fprintf('Enter time_step into runme!\n');
   %   return
   %end

   % Stress-balance parameters
   md.stressbalance.restol = 1e-4;
   md.stressbalance.maxiter = 10;
   %md.verbose.convergence = 1;

   % Go solve
   md.verbose.solution=1;
   md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   md = solve(md,'tr');
   if md.settings.waitonlock == 0
      fprintf('\n \033[103;30m Load results with: md = loadresultsfromcluster(md,''runtimename'',''%s''); \033[0m \n', md.private.runtimename);
      fprintf(' \033[103;30m Run the following:\n');
      fprintf(' \033[103;30m  tmpResults = [md.results.TransientSolution4 md.results.TransientSolution];\n');
      fprintf(' \033[103;30m  md.results.TransientSolution = tmpResults;\n');
      fprintf(' \033[103;30m  md.results = rmfield(md.results, ''TransientSolution4'');\n');
      fprintf(' \033[103;30m  clear tmpResults \033[0m\n');
      fprintf(' \033[103;30m Save results with: save(''%s.mat'', ''md'', ''-v7.3''); \033[0m \n\n', [org.repository filesep org.prefix ...
         strrep(org.steps(end).string, '_extra', '')]);
      return
   end

   % Combine results
   tmpResults = [md.results.TransientSolution4 md.results.TransientSolution];
   md.results.TransientSolution = tmpResults;
   md.results = rmfield(md.results, 'TransientSolution4');
   clear tmpResults

   % Save - cheat and save it back to Transient_Steadystate_noFloating
   save([org.repository filesep org.prefix strrep(org.steps(end).string, '_extra', '') '.mat'], 'md', '-v7.3');

end % }}}

% Switch sliding law
switch friction_law %%{{{
   case 'Budd'
      m = m_Budd;
   case 'Weertman'
      m = m_Weertman;
end %%}}}
if perform(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m)]),%    STEP 8{{{

   md = loadmodel(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm']);
   md = transientrestart(md);

	md.cluster = cluster; generic('name', oshostname(), 'np', 1);
	md.cluster.interactive = 1;
	md.settings.waitonlock = 10;
   
   % Calculate velocity before sliding law swap
   %[~, ~, taub_orig] = basalstress(md);
   md.timestepping.final_time = md.timestepping.start_time + md.timestepping.time_step;
   md = solve(md, 'sb');
   md.initialization.vx  = md.results.StressbalanceSolution.Vx;
   md.initialization.vy  = md.results.StressbalanceSolution.Vy;
   md.initialization.vel = md.results.StressbalanceSolution.Vel;
   %[~, ~, taub1] = basalstress(md);
   ub1 = md.results.StressbalanceSolution.Vel;

   
   % Switch the friction coefficient
   if strcmp(friction_law, 'Budd') && m == 1
      % skip
      fprintf('\n');
      fprintf('no change to sliding\n');
      fprintf('\n');
   else
      switch friction_law
         case 'Budd'
            fprintf('\n');
            fprintf(['Switching to Budd sliding (m=' num2str(m_Budd) ')\n']);
            fprintf('\n');
            md = friction_coefficient_conversion(md, 'budd', 'budd', 'p', 1/m_Budd, 'q', 1/m_Budd);
         case 'Weertman'
            fprintf('\n');
            fprintf(['Switching to Weertman sliding (m=' num2str(m_Weertman) ')\n']);
            fprintf('\n');
            md = friction_coefficient_conversion(md, 'budd', 'weertman', 'm', 1/m_Weertman);
         case 'Schoof'
            fprintf('\n');
            fprintf(['Switching to Schoof sliding (m=' num2str(m_Schoof) ')\n']);
            fprintf('\n');
            %md = friction_coefficient_conversion(md, 'budd', 'schoof', 'm', 1/m_Schoof, 'C_max', 0.4);
            
            % TODO DEBUG
            [~, ~, b_budd] = basalstress(md);
            md_schoof = friction_coefficient_conversion(md, 'budd', 'schoof', 'm', 1, 'C_max', 0.4);
            [~, ~, b_schoof] = basalstress(md_schoof);
            return
            % TODO DEBUG
      end
   end

   % Calculate velocity after sliding law swap
   md = solve(md, 'sb');
   md.initialization.vx  = md.results.StressbalanceSolution.Vx;
   md.initialization.vy  = md.results.StressbalanceSolution.Vy;
   md.initialization.vel = md.results.StressbalanceSolution.Vel;
   %[~, ~, taub2] = basalstress(md);
   ub2 = md.results.StressbalanceSolution.Vel;

   ub_diff_rel = (ub2-ub1) ./ ub1;
   fprintf('\n');
   fprintf('Mean absolute relative difference in velocity = %10.8f m/yr\n', mean(abs(ub_diff_rel), 'omitnan'));
   fprintf('\n');

   % Save
   savemodel(org,md);

end % }}}

if perform(org, ['Transient_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m) '_noMotion']),% {{{

   md = loadmodel(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m)]);

   % Transfer results to model fields
   md.timestepping.final_time = md.timestepping.start_time + 30;
   if meshsize == 100
      md.timestepping.time_step = 0.005;
      md.settings.output_frequency = 100;
   elseif meshsize == 200
      md.timestepping.time_step = 0.010;
      md.settings.output_frequency = 50;
   elseif meshsize == 400
      md.timestepping.time_step = 0.020;
      md.settings.output_frequency = 25;
   else
      fprintf('Enter time_step into runme!\n');
      return
   end

   % Stress-balance parameters
   md.stressbalance.restol = 1e-4;
   md.stressbalance.maxiter = 10;
   %md.verbose.convergence = 1;

   % Set the requested outputs
   md.stressbalance.requested_outputs={'default','StrainRatexx','StrainRateyy','StrainRatexy','StrainRateeffective'};
   md.stressbalance.requested_outputs={'default','StrainRatexx','StrainRateyy','StrainRatexy','StrainRateeffective', ...
      'DeviatoricStressxx','DeviatoricStressxy','DeviatoricStressyy'};
   md.transient.requested_outputs={'default','IceVolumeAboveFloatation'};

   % Go solve
   md.verbose.solution=1;
   md.cluster = cluster;
   md.settings.waitonlock = waitonlock;
   md = solve(md,'transient');
   if md.settings.waitonlock == 0
      fprintf('\n \033[103;30m Load results with: md = loadresultsfromcluster(md,''runtimename'',''%s''); \033[0m \n', md.private.runtimename);
      fprintf(' \033[103;30m Save results with: save(''%s.mat'', ''md'', ''-v7.3''); \033[0m \n\n', [org.repository filesep org.prefix org.steps(end).string]);
      return
   end

   % Save
   savemodel(org,md);

end % }}}
if perform(org, ['Transient_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m) '_1yearly' num2str(magnitude) 'mChannel']),% {{{

   md = loadmodel(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m)]);

   % Transfer results to model fields
   md.timestepping.final_time = md.timestepping.start_time + 30;
   if meshsize == 100
      md.timestepping.time_step = 0.005;
      md.settings.output_frequency = 10;
   elseif meshsize == 200
      md.timestepping.time_step = 0.010;
      md.settings.output_frequency = 5;
   elseif meshsize == 400
      md.timestepping.time_step = 0.020;
      md.settings.output_frequency = 2;
   elseif meshsize == 1000
      md.timestepping.time_step = 0.050;
      md.settings.output_frequency = 1;
   else
      fprintf('Enter time_step into runme!\n');
      return
   end
   
   % Stress-balance parameters
   md.stressbalance.restol = 1e-4;
   md.stressbalance.maxiter = 20; 10;

	% % TODO DEBUG
	% md.verbose.convergence = true;
	% %md.stressbalance.abstol = 10;
   % %md.stressbalance.restol = 1;
	% %md.timestepping.time_step = 0.005;
	% md.timestepping.final_time = md.timestepping.start_time + 5;
   % md.settings.output_frequency = 1;
   % md.stressbalance.maxiter = 20;

   % Activate moving boundary
   md.transient.ismovingfront = 1;
   % md.calving.calvingrate = zeros(md.mesh.numberofvertices,1);
   % md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);

   % Set up spclevelset
   % NOTE: Terminus retreats/re-advances by [magnitude] over 1 year
   levelset0 = md.mask.ice_levelset;
   md.levelset.spclevelset = [];
   md.levelset.spclevelset(:,end+1) = [levelset0; md.timestepping.start_time];
   for time = md.timestepping.start_time : 1 : md.timestepping.final_time
      signeddistance = move_terminus_levelset(md, levelset0, magnitude, -1, true);
     
      signeddistance(md.geometry.bed>0 & levelset0<0) = -1;
      pos      = find(signeddistance<0);

      if exist('TEMP.exp','file'), delete('TEMP.exp'); end
      expcontourlevelzero(md, signeddistance, 0, 'TEMP.exp');
      signeddistance = abs(ExpToLevelSet(md.mesh.x, md.mesh.y, 'TEMP.exp'));
      delete('TEMP.exp');
      signeddistance(pos) = -signeddistance(pos);
      
      md.levelset.spclevelset(:,end+1) = [signeddistance; time + 0.5];
      %signeddistance = move_terminus_levelset(md, levelset, magnitude, +1, true);
      md.levelset.spclevelset(:,end+1) = [levelset0; time + 1.0];
   end
   
   % Set the requested outputs
   md.stressbalance.requested_outputs={'default','StrainRatexx','StrainRateyy','StrainRatexy','StrainRateeffective'};
   md.transient.requested_outputs={'default','IceVolumeAboveFloatation'};

   % Go solve
   md.verbose.solution=1;
   md.cluster = cluster; if meshsize == 100 & contains(cluster.name, 'discover'), md.cluster.time = 12*3600; end
   md.settings.waitonlock = waitonlock;
   md = solve(md,'transient');
   if md.settings.waitonlock == 0
      fprintf('\n \033[103;30m Load results with: md = loadresultsfromcluster(md,''runtimename'',''%s''); \033[0m \n', md.private.runtimename);
      fprintf(' \033[103;30m Save results with: save(''%s.mat'', ''md'', ''-v7.3''); \033[0m \n\n', [org.repository filesep org.prefix org.steps(end).string]);
      return
   end

   % Save
   savemodel(org,md);

end % }}}
if perform(org, ['Transient_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m) ...
      '_1yearly' num2str(magnitude) 'mChannel_Seasonality1']),% {{{

   md = loadmodel(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m)]);

   % Transfer results to model fields
   md.timestepping.final_time = md.timestepping.start_time + 30;
   if meshsize == 100
      md.timestepping.time_step = 0.005;
      md.settings.output_frequency = 10;
   elseif meshsize == 200
      md.timestepping.time_step = 0.010;
      md.settings.output_frequency = 5; %50;
   elseif meshsize == 400
      md.timestepping.time_step = 0.020;
      md.settings.output_frequency = 25;
   elseif meshsize == 1000
      md.timestepping.time_step = 0.050;
      md.settings.output_frequency = 1;
   else
      fprintf('Enter time_step into runme!\n');
      return
   end

   % Stress-balance parameters
   md.stressbalance.restol = 1e-4;
   md.stressbalance.maxiter = 20;
   %md.verbose.convergence = 1;

   % Activate moving boundary
   md.transient.ismovingfront = 1;

   % Set up spclevelset
   % NOTE: Terminus retreats/re-advances by [magnitude] over 1 year
   levelset0 = md.mask.ice_levelset;
   md.levelset.spclevelset = [];
   md.levelset.spclevelset(:,end+1) = [levelset0; md.timestepping.start_time];
   for time = md.timestepping.start_time : 1 : md.timestepping.final_time
      signeddistance = move_terminus_levelset(md, levelset0, magnitude, -1, true);
     
      signeddistance(md.geometry.bed>0 & levelset0<0) = -1;
      pos      = find(signeddistance<0);

      if exist('TEMP.exp','file'), delete('TEMP.exp'); end
      expcontourlevelzero(md, signeddistance, 0, 'TEMP.exp');
      signeddistance = abs(ExpToLevelSet(md.mesh.x, md.mesh.y, 'TEMP.exp'));
      delete('TEMP.exp');
      signeddistance(pos) = -signeddistance(pos);
      
      md.levelset.spclevelset(:,end+1) = [signeddistance; time + 0.25];
      %signeddistance = move_terminus_levelset(md, levelset, magnitude, +1, true);
      md.levelset.spclevelset(:,end+1) = [levelset0; time + 1.0];
   end
   
   % Set the requested outputs
   md.stressbalance.requested_outputs={'default','StrainRatexx','StrainRateyy','StrainRatexy','StrainRateeffective'};
   md.transient.requested_outputs={'default','IceVolumeAboveFloatation'};

   % Go solve
   md.verbose.solution=1;
   md.cluster = cluster; if meshsize == 100 & contains(cluster.name, 'discover'), md.cluster.time = 12*3600; end
   md.settings.waitonlock = waitonlock;
   md = solve(md,'transient');
   if md.settings.waitonlock == 0
      fprintf('\n \033[103;30m Load results with: md = loadresultsfromcluster(md,''runtimename'',''%s''); \033[0m \n', md.private.runtimename);
      fprintf(' \033[103;30m Save results with: save(''%s.mat'', ''md'', ''-v7.3''); \033[0m \n\n', [org.repository filesep org.prefix org.steps(end).string]);
      return
   end

   % Save
   savemodel(org,md);

end % }}}
if perform(org, ['Transient_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m) ...
      '_1yearly' num2str(magnitude) 'mChannel_Seasonality2']),% {{{

   md = loadmodel(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m)]);

   % Transfer results to model fields
   md.timestepping.final_time = md.timestepping.start_time + 30;
   if meshsize == 100
      md.timestepping.time_step = 0.005;
      md.settings.output_frequency = 10;
   elseif meshsize == 200
      md.timestepping.time_step = 0.010;
      md.settings.output_frequency = 5; %50;
   elseif meshsize == 400
      md.timestepping.time_step = 0.020;
      md.settings.output_frequency = 25;
   elseif meshsize == 1000
      md.timestepping.time_step = 0.050;
      md.settings.output_frequency = 1;
   else
      fprintf('Enter time_step into runme!\n');
      return
   end

   % Stress-balance parameters
   md.stressbalance.restol = 1e-4;
   md.stressbalance.maxiter = 20;
   %md.verbose.convergence = 1;

   % Activate moving boundary
   md.transient.ismovingfront = 1;

   % Set up spclevelset
   % NOTE: Terminus retreats/re-advances by [magnitude] over 1 year
   levelset0 = md.mask.ice_levelset;
   md.levelset.spclevelset = [];
   md.levelset.spclevelset(:,end+1) = [levelset0; md.timestepping.start_time];
   for time = md.timestepping.start_time : 1 : md.timestepping.final_time
      signeddistance = move_terminus_levelset(md, levelset0, magnitude, -1, true);
     
      signeddistance(md.geometry.bed>0 & levelset0<0) = -1;
      pos      = find(signeddistance<0);

      if exist('TEMP.exp','file'), delete('TEMP.exp'); end
      expcontourlevelzero(md, signeddistance, 0, 'TEMP.exp');
      signeddistance = abs(ExpToLevelSet(md.mesh.x, md.mesh.y, 'TEMP.exp'));
      delete('TEMP.exp');
      signeddistance(pos) = -signeddistance(pos);
      
      md.levelset.spclevelset(:,end+1) = [signeddistance; time + 0.75];
      %signeddistance = move_terminus_levelset(md, levelset, magnitude, +1, true);
      md.levelset.spclevelset(:,end+1) = [levelset0; time + 1.0];
   end
   
   % Set the requested outputs
   md.stressbalance.requested_outputs={'default','StrainRatexx','StrainRateyy','StrainRatexy','StrainRateeffective'};
   md.transient.requested_outputs={'default','IceVolumeAboveFloatation'};

   % Go solve
   md.verbose.solution=1;
   md.cluster = cluster; if meshsize == 100 & contains(cluster.name, 'discover'), md.cluster.time = 12*3600; end
   md.settings.waitonlock = waitonlock;
   md = solve(md,'transient');
   if md.settings.waitonlock == 0
      fprintf('\n \033[103;30m Load results with: md = loadresultsfromcluster(md,''runtimename'',''%s''); \033[0m \n', md.private.runtimename);
      fprintf(' \033[103;30m Save results with: save(''%s.mat'', ''md'', ''-v7.3''); \033[0m \n\n', [org.repository filesep org.prefix org.steps(end).string]);
      return
   end

   % Save
   savemodel(org,md);

end % }}}
if perform(org, ['Transient_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m) '_meanTerminus' num2str(magnitude) 'mChannel']),% {{{

   md = loadmodel(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m)]);

   % Transfer results to model fields
   md.timestepping.final_time = md.timestepping.start_time + 50;
   if meshsize == 100
      md.timestepping.time_step = 0.005;
      md.settings.output_frequency = 10;
   elseif meshsize == 200
      md.timestepping.time_step = 0.010;
      md.settings.output_frequency = 50;
      %md.settings.output_frequency = 1;
   elseif meshsize == 400
      md.timestepping.time_step = 0.020;
      md.settings.output_frequency = 25;
   elseif meshsize == 1000
      md.timestepping.time_step = 0.050;
      md.settings.output_frequency = 10;
   else
      fprintf('Enter time_step into runme!\n');
      return
   end

	% % TODO DEBUG - high res run
	% md.timestepping.final_time = md.timestepping.start_time + 10;
	% md.settings.output_frequency = 1;

   % Stress-balance parameters
   md.stressbalance.restol = 1e-4;
   md.stressbalance.maxiter = 10;
   %md.verbose.convergence = 1;

   % Activate moving boundary
   md.transient.ismovingfront = 1;

   % Set the terminus position
   levelset0 = md.mask.ice_levelset;
   
   signeddistance = move_terminus_levelset(md, levelset0, magnitude/2, -1, true);
   signeddistance(md.geometry.bed>0 & levelset0<0) = -1;
   pos      = find(signeddistance<0);

   if exist('TEMP.exp','file'), delete('TEMP.exp'); end
   expcontourlevelzero(md, signeddistance, 0, 'TEMP.exp');
   signeddistance = abs(ExpToLevelSet(md.mesh.x, md.mesh.y, 'TEMP.exp'));
   delete('TEMP.exp');
   signeddistance(pos) = -signeddistance(pos);
   
   if false
   levelset = zeros(md.mesh.numberofvertices,1);
   levelset(md.mesh.x >  terminus0_x-magnitude/2) = +1;
   levelset(md.mesh.x <= terminus0_x-magnitude/2) = -1;

   levelset(md.geometry.bed>0 & levelset0<0) = -1;
   
   % Create signed distance field
   disp('Converting levelsets to signed distance fields');
   pos      = find(levelset<0);

   if exist('TEMP.exp','file'), delete('TEMP.exp'); end
   expcontourlevelzero(md, levelset, 0, 'TEMP.exp');
   % (TODO possibly) Check that there is only one zero-level contour
   signeddistance = abs(ExpToLevelSet(md.mesh.x, md.mesh.y, 'TEMP.exp'));
   delete('TEMP.exp');
   signeddistance(pos) = -signeddistance(pos);
   end
   
   md.levelset.spclevelset = zeros(length(levelset0)+1,3);
   md.levelset.spclevelset(1:end-1,1) = levelset0;
   md.levelset.spclevelset(end,1) = md.timestepping.start_time;
   md.levelset.spclevelset(1:end-1,2) = signeddistance;
   md.levelset.spclevelset(end,2) = md.timestepping.start_time + 1.5;
   md.levelset.spclevelset(1:end-1,3) = signeddistance;
   md.levelset.spclevelset(end,3) = md.timestepping.final_time;

   % Set the requested outputs
   md.stressbalance.requested_outputs={'default','StrainRatexx','StrainRateyy','StrainRatexy','StrainRateeffective'};
   md.stressbalance.requested_outputs={'default','StrainRatexx','StrainRateyy','StrainRatexy','StrainRateeffective', ...
      'DeviatoricStressxx','DeviatoricStressxy','DeviatoricStressyy'};
   md.transient.requested_outputs={'default','IceVolumeAboveFloatation'};

   % Go solve
   md.verbose.solution=1;
   md.cluster = cluster; if meshsize == 100 & contains(cluster.name, 'discover'), md.cluster.time = 12*3600; end
   md.settings.waitonlock = waitonlock;
   md = solve(md,'transient');
   if md.settings.waitonlock == 0
      fprintf('\n \033[103;30m Load results with: md = loadresultsfromcluster(md,''runtimename'',''%s''); \033[0m \n', md.private.runtimename);
      fprintf(' \033[103;30m Save results with: save(''%s.mat'', ''md'', ''-v7.3''); \033[0m \n\n', [org.repository filesep org.prefix org.steps(end).string]);
      return
   end

   % Save
   savemodel(org,md);

end % }}}

if perform(org, ['Transient_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m) '_mostRetreatedTerminus' num2str(magnitude) 'mChannel']),% {{{

   md = loadmodel(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m)]);

   % Transfer results to model fields
   md.timestepping.final_time = md.timestepping.start_time + 50;
   if meshsize == 100
      md.timestepping.time_step = 0.005;
      md.settings.output_frequency = 10;
   elseif meshsize == 200
      md.timestepping.time_step = 0.010;
      md.settings.output_frequency = 50;
      %md.settings.output_frequency = 1;
   elseif meshsize == 400
      md.timestepping.time_step = 0.020;
      md.settings.output_frequency = 25;
   elseif meshsize == 1000
      md.timestepping.time_step = 0.050;
      md.settings.output_frequency = 10;
   else
      fprintf('Enter time_step into runme!\n');
      return
   end

	% % TODO DEBUG - high res run
	% md.timestepping.final_time = md.timestepping.start_time + 10;
	% md.settings.output_frequency = 1;

   % Stress-balance parameters
   md.stressbalance.restol = 1e-4;
   md.stressbalance.maxiter = 10;
   %md.verbose.convergence = 1;

   % Activate moving boundary
   md.transient.ismovingfront = 1;

   % Set the terminus position
   levelset0 = md.mask.ice_levelset;
   
   signeddistance = move_terminus_levelset(md, levelset0, magnitude, -1, true);
   signeddistance(md.geometry.bed>0 & levelset0<0) = -1;
   pos      = find(signeddistance<0);

   if exist('TEMP.exp','file'), delete('TEMP.exp'); end
   expcontourlevelzero(md, signeddistance, 0, 'TEMP.exp');
   signeddistance = abs(ExpToLevelSet(md.mesh.x, md.mesh.y, 'TEMP.exp'));
   delete('TEMP.exp');
   signeddistance(pos) = -signeddistance(pos);
   
   if false
   levelset = zeros(md.mesh.numberofvertices,1);
   levelset(md.mesh.x >  terminus0_x-magnitude/2) = +1;
   levelset(md.mesh.x <= terminus0_x-magnitude/2) = -1;

   levelset(md.geometry.bed>0 & levelset0<0) = -1;
   
   % Create signed distance field
   disp('Converting levelsets to signed distance fields');
   pos      = find(levelset<0);

   if exist('TEMP.exp','file'), delete('TEMP.exp'); end
   expcontourlevelzero(md, levelset, 0, 'TEMP.exp');
   % (TODO possibly) Check that there is only one zero-level contour
   signeddistance = abs(ExpToLevelSet(md.mesh.x, md.mesh.y, 'TEMP.exp'));
   delete('TEMP.exp');
   signeddistance(pos) = -signeddistance(pos);
   end
   
   md.levelset.spclevelset = zeros(length(levelset0)+1,3);
   md.levelset.spclevelset(1:end-1,1) = levelset0;
   md.levelset.spclevelset(end,1) = md.timestepping.start_time;
   md.levelset.spclevelset(1:end-1,2) = signeddistance;
   md.levelset.spclevelset(end,2) = md.timestepping.start_time + 1.5;
   md.levelset.spclevelset(1:end-1,3) = signeddistance;
   md.levelset.spclevelset(end,3) = md.timestepping.final_time;

   % Set the requested outputs
   md.stressbalance.requested_outputs={'default','StrainRatexx','StrainRateyy','StrainRatexy','StrainRateeffective'};
   md.stressbalance.requested_outputs={'default','StrainRatexx','StrainRateyy','StrainRatexy','StrainRateeffective', ...
      'DeviatoricStressxx','DeviatoricStressxy','DeviatoricStressyy'};
   md.transient.requested_outputs={'default','IceVolumeAboveFloatation'};

   % Go solve
   md.verbose.solution=1;
   md.cluster = cluster; if meshsize == 100 & contains(cluster.name, 'discover'), md.cluster.time = 12*3600; end
   md.settings.waitonlock = waitonlock;
   md = solve(md,'transient');
   if md.settings.waitonlock == 0
      fprintf('\n \033[103;30m Load results with: md = loadresultsfromcluster(md,''runtimename'',''%s''); \033[0m \n', md.private.runtimename);
      fprintf(' \033[103;30m Save results with: save(''%s.mat'', ''md'', ''-v7.3''); \033[0m \n\n', [org.repository filesep org.prefix org.steps(end).string]);
      return
   end

   % Save
   savemodel(org,md);

end % }}}

