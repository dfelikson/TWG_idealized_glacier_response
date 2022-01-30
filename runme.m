steps = [8];
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
terminus0_x = 26500;
magnitude = 1000;

% Cluster parameters
cluster = generic('name', oshostname(), 'np', 2);
cluster.interactive = 1;
waitonlock = 10;

% NOTE
cluster=discover;
cluster.name='discover.nccs.nasa.gov';
cluster.login='dfelikso';
cluster.project='s2133';
cluster.numnodes=4;
cluster.cpuspernode=nan; %16;
cluster.time=5.0*60*60;
cluster.processor='sand';
cluster.queue='allnccs';
cluster.codepath='/discover/nobackup/dfelikso/Software/ISSM/trunk-jpl/bin';
cluster.executionpath='/discover/nobackup/dfelikso/Software/ISSM/trunk-jpl/execution';
cluster.email='denis.felikson@nasa.gov';

cluster.interactive = 0;
waitonlock = 0;
% NOTE

% % NOTE
% cluster = generic('name', 'gs615-oibserve.ndc.nasa.gov', 'np', 28, ...
%    'login', 'dfelikso', ...
%    'codepath', '/home/dfelikso/Software/ISSM/trunk-jpl/bin', ...
%    'etcpath', '/home/dfelikso/Software/ISSM/trunk-jpl/etc', ...
%    'executionpath', '/home/dfelikso/Projects/GrIS_Calibrated_SLR/ISSM/execution');
% cluster.interactive = 0;
% waitonlock = 0;
% % NOTE

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

   md.timestepping.time_step = 0.050;
   md.timestepping.final_time=50;
	md.settings.output_frequency=10;
   
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
   md.timestepping.final_time = md.timestepping.start_time + 200;
   md.timestepping.time_step = 0.010;
   md.settings.output_frequency = 50;

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

%  Steady-state
if perform(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm']),%    STEP 5{{{

   % Start with Transient_Steadystate_noFloating because this gets us close to the terminus position
   % that we want to start with
   md = loadmodel(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_noFloating']);

   % Transfer results to model fields
   md = transientrestart(md);
   md.timestepping.final_time = md.timestepping.start_time + 80;
   md.timestepping.time_step = 0.010;
   md.settings.output_frequency = 50;

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

% Switch sliding law
switch friction_law %%{{{
   case 'Budd'
      m = m_Budd;
   case 'Weertman'
      m = m_Weertman;
end %%}}}
if perform(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m)]),%    STEP 6{{{

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
   md.timestepping.time_step = 0.010;
   md.settings.output_frequency = 50;

   % Stress-balance parameters
   md.stressbalance.restol = 1e-4;
   md.stressbalance.maxiter = 10;

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
if perform(org, ['Transient_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m) '_1yearly' num2str(magnitude) 'm']),% {{{

   md = loadmodel(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m)]);

   % Transfer results to model fields
   md.timestepping.final_time = md.timestepping.start_time + 30;
   md.timestepping.time_step = 0.010;
   md.settings.output_frequency = 5;
   
   % Stress-balance parameters
   md.stressbalance.restol = 1e-4;
   md.stressbalance.maxiter = 20; 10;

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
if perform(org, ['Transient_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m) '_meanTerminus' num2str(magnitude) 'm']),% {{{

   md = loadmodel(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m)]);

   % Transfer results to model fields
   md.timestepping.final_time = md.timestepping.start_time + 30;
   md.timestepping.time_step = 0.010;
   md.settings.output_frequency = 50;

   % Stress-balance parameters
   md.stressbalance.restol = 1e-4;
   md.stressbalance.maxiter = 10;

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
if perform(org, ['Transient_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m) '_mostRetreatedTerminus' num2str(magnitude) 'm']),% {{{

   md = loadmodel(org, ['Transient_Steadystate_frictionCoeffMin' num2str(friction_coefficient_min) '_terminus' num2str(terminus0_x) 'm_friction' friction_law sprintf('%3.1f',m)]);

   % Transfer results to model fields
   md.timestepping.final_time = md.timestepping.start_time + 30;
   md.timestepping.time_step = 0.010;
   md.settings.output_frequency = 50;

   % Stress-balance parameters
   md.stressbalance.restol = 1e-4;
   md.stressbalance.maxiter = 10;

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

