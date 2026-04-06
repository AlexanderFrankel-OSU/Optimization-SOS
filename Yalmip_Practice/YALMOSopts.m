%% Solve SOS Problem
%strPrnt=[dirRslts,strID(1:end-4),'SolverOutput.txt'];
%diary(strPrnt)
% disp(['Solving uncertain system with Re=',num2str(Re),...
% ' and alpha=',num2str(alpha)]);
options = sdpsettings('solver','mosek','sos.congruence',1, 'sos.reuse', 0); %,'sos.newton',0,'sos.congruence',0);
% options = sdpsettings('solver', 'mosek', ...
% 'mosek.MSK_DPAR_OPTIMIZER_MAX_TIME', 100.0, ...
% 'mosek.MSK_IPAR_INTPNT_SOLVE_FORM', 'MSK_SOLVE_DUAL', ...
% 'verbose', 1, ...
% 'savedebug', 1, ...
% 'mosektaskfile', 'dump.task.gz', ...
% 'mosek.MSK_IPAR_NUM_THREADS', 16);
% options.sdpt3.maxit=46;
% options.mosek.MSK_IPAR_NUM_THREADS = 16;
dtol = 0.001;
% options.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS = dtol*options.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS; % Dual feasibility tolerance used by the conic interior-point optimizer
% options.mosek.MSK_DPAR_INTPNT_CO_TOL_MU_RED = dtol*options.mosek.MSK_DPAR_INTPNT_CO_TOL_MU_RED; % Optimality tolerance for the conic solver
% options.mosek.MSK_DPAR_INTPNT_CO_TOL_NEAR_REL = 100/dtol; % Optimality tolerance for the conic solver
% options.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS = dtol*options.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS; % Primal feasibility tolerance used by the conic interior-point optimizer
% options.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = dtol*options.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP; % Relative gap termination tolerance used by the conic interior-point optimizer
% options.mosek.MSK_DPAR_INTPNT_CO_TOL_INFEAS = 1e-15; % Smaller means less likely to declare ill-posed or infeasible
% options.mosek.MSK_IPAR_INTPNT_ORDER_METHOD = 'MSK_ORDER_METHOD_FREE';
% options.mosek.MSK_IPAR_INTPNT_ORDER_METHOD = 'MSK_ORDER_METHOD_APPMINLOC';
% time3=tic;%preprocessing
% [cnstr,obj]=compilesos(Fcnstr,optgoal,options,optvar);
% options.mosek
%options.savesolverinput = 1;
%profile clear;
%profile on -memory;
[sol,v,Q]=solvesos(cnstr,[],options,optvar);
%profile report;