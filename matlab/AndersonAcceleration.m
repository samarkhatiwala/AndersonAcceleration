function [xsol,iter,aa] = AndersonAcceleration(g,x,AAparams,histParams,restartFile,fileSuff,aa)

% This function implements Anderson Acceleration for fixed point iteration. It is essentially 
% based on AndAcc.m originally written by Homer F. Walker (listed in H. Walker, Anderson acceleration: 
% Algorithms and implementations, Worcester Polytechnic Institute Mathematical Sciences Department 
% Research Report MS-6-15-50, June, 2011) with extensive modifications as per Khatiwala (2022) 
% (https://doi.org/xxx) to apply it to the computation of equilibrium solutions of periodically-forced 
% ocean models on batch HPC systems. If you use this code please cite BOTH sources.

% Input arguments:
%  Required:
%   g: function handle for the fixed point map (see below for arguments)
%   x: initial guess vector
%  Optional:
%   AAparams: struct with fields defining algorithm parameters (see defaults below)
%   histParams: struct with fields defining history parameters (see defaults below)
%   restartFile: name of file to which to save restart data (or [] if none required)
%   fileSuff: filename suffix to tag checkpoint files
%   aa: struct containing state of algorithm if restarting from a previous iteration 
%       (or [] if this is the first iteration)

% Output:
%   xsol: final solution
%   iter: final iteration number
%   aa: struct containing algorithm state

% Definition of the g function:
%  [gx,gv,vnorms,externalconvergence]=g(x,fetchOutput)
%  Inputs:
%    x: iterate vector at which to evaluate g
%    fetchOutput: flag (0, 1 or []) to indicate 'mode' in which to evaluate g
%  Outputs:
%    gx: output vector from evaluating g
%    gv: dimensional (rescaled) output (not used)
%    vnorms: dimensional residual norm (a vector if this is a multi-tracer problem, otherwise a scalar)
%    externalconvergence: struct of external convergence data

if nargin < 2
  error('AndAcc requires at least two arguments.')
end

% Figure out whether to instantiate the aa structure
if nargin<7 || isempty(aa)
  iniAA=1;
  disp('Initializing AA')
else
  iniAA=0;
end

if nargin>6 && ~isempty(aa)
  reuseAA=1;
  disp('AA is being reused')
else
  reuseAA=0;
end

if nargin<5
  disp('No restart')
  doRestart=0;
  doRestartWithFile=0;
  doRestartWithStruct=0;
else
  if isempty(restartFile) && nargin<7
    doRestart=0;
    disp('No restart')
  else
    if ~isempty(restartFile) && reuseAA
      error('ERROR: Cannot specify both restart file and aa struct')
    end
    if nargin>6  % AA takes priority
      doRestart=1;
      doRestartWithFile=0;
      doRestartWithStruct=1;
      disp('Restart with AA')
    elseif ~isempty(restartFile)
      doRestart=1;
      doRestartWithFile=1;
      doRestartWithStruct=0;
      disp('Restart with file')
	  if nargin<6
		fileSuff=[];
	  end
    else
      error('ERROR: must specify one of restart file or aa struct')
    end
  end    
end

if ~doRestartWithFile || (doRestartWithFile && ~exist(restartFile,'file'))
% Instantiate the aa structure
% Set algorithm parameter defaults and overwrite with values from AAparams: 
% mMax: maximum number of stored residuals (non-negative integer).
%       NOTE: mMax = 0 => no acceleration.
% itmax: maximum allowable number of iterations.
% atol: absolute error tolerance.
% rtol: relative error tolerance.
% droptol: tolerance for dropping stored residual vectors to improve conditioning
%          NOTE: If droptol > 0, drop residuals if the condition number exceeds droptol; 
%          if droptol <= 0, do not drop residuals.
% beta: damping factor: If beta > 0 (and beta ~= 1), then the step is damped by beta; 
%       otherwise, the step is not damped.
%       NOTE: beta can be a function handle; form beta(iter), where iter is the iteration 
%       number and 0 < beta(iter) <= 1.
% AAstart: acceleration delay factor: If AAstart > 0, start acceleration when iter = AAstart.
% restartAANormStagnation: flag to indicate whether to monitor for stagnation of the method
% restartAANormDiff: not currently used.
% restartAASuccessiveIters: maximum number of successive iteration for which the norm increases 
% before previous iterates are discarded (only used if restartAANormStagnation=1).
% restartAAPeriodic: If restartAAPeriodic>0, the previous iterates are discarded every 
%                    restartAAPeriodic iterations. (If restartAAPeriodic>mMax, an error is 
%                    thrown.)
  if iniAA
	aa.x = x;  
	aa.mMax = min(50, size(aa.x,1));
	aa.itmax = 100;
	aa.atol = 1.e-10;
	aa.rtol = 1.e-10;
	aa.droptol = 1.e10;
	aa.beta = 1;
	aa.AAstart = 0;
	aa.restartAANormStagnation = 0;
	aa.restartAANormDiff = 0.01;
	aa.restartAASuccessiveIters = 2;
	aa.restartAAPeriodic = 0;

	if nargin>2 && ~isempty(AAparams)
	  if isfield(AAparams,'mMax') && ~isempty(AAparams.mMax)
		aa.mMax=AAparams.mMax;
	  end
	  if isfield(AAparams,'itmax') && ~isempty(AAparams.itmax)
		aa.itmax=AAparams.itmax;
	  end
	  if isfield(AAparams,'atol') && ~isempty(AAparams.atol)
		aa.atol=AAparams.atol;
	  end
	  if isfield(AAparams,'rtol') && ~isempty(AAparams.rtol)
		aa.rtol=AAparams.rtol;
	  end
	  if isfield(AAparams,'droptol') && ~isempty(AAparams.droptol)
		aa.droptol=AAparams.droptol;
	  end
	  if isfield(AAparams,'beta') && ~isempty(AAparams.beta)
		aa.beta=AAparams.beta;
	  end
	  if isfield(AAparams,'AAstart') && ~isempty(AAparams.AAstart)
		aa.AAstart=AAparams.AAstart;
	  end
	  if isfield(AAparams,'restartAANormStagnation') && ~isempty(AAparams.restartAANormStagnation)
		aa.restartAANormStagnation=AAparams.restartAANormStagnation;
	  end  
	  if isfield(AAparams,'restartAANormDiff') && ~isempty(AAparams.restartAANormDiff)
		aa.restartAANormDiff=AAparams.restartAANormDiff;
	  end
	  if isfield(AAparams,'restartAASuccessiveIters') && ~isempty(AAparams.restartAASuccessiveIters)
		aa.restartAASuccessiveIters=AAparams.restartAASuccessiveIters;
	  end
	  if isfield(AAparams,'restartAAPeriodic') && ~isempty(AAparams.restartAAPeriodic)
		aa.restartAAPeriodic=AAparams.restartAAPeriodic;
	  end
	end

	if aa.mMax == 0
	  disp('Warning: mMax is 0; no acceleration will be applied!')
	end  
	if aa.mMax < 0
	  error('ERROR: mMax must be >=0!')
	end

    if aa.restartAAPeriodic>aa.mMax
       error(['ERROR: restartAAPeriodic exceeds mMax!'])
    end

	aa.DG = []; % Storage of g-value differences.
	aa.Q =[];
	aa.R =[];
	aa.restart_count = 0;
	aa.flushAA = 0;
	aa.mAA = 0;
	aa.nfeval=0;
	aa.iter0=0;
	aa.fetchOutput=0;
	aa.res_hist = []; % Storage of residual history (a matrix with columns: 
%                   	[iter nfeval vnorms res_norm]). This will be overwriiten below.
    
%   Set history defaults and overwrite with values from histParams
%   nhistfreq = frequency with which to save history
%   ncheckpointfreq = frequency with which to write checkpoint (if < 0 
%                     then no checkpoints are written; this is the default)
	aa.nhistfreq=20;
	aa.ncheckpointfreq=-1;
	if nargin>3 && ~isempty(histParams)
	  if isfield(histParams,'nhistfreq') && ~isempty(histParams.nhistfreq)
		aa.nhistfreq=histParams.nhistfreq;
	  end
	  if isfield(histParams,'ncheckpointfreq') && ~isempty(histParams.ncheckpointfreq)
		aa.ncheckpointfreq=histParams.ncheckpointfreq;
	  end
	end
    aa.nhistmax=ceil(aa.itmax/aa.nhistfreq)+5; % space to store history
	
%   initialize history
	aa.aahist.iterhist=zeros(aa.nhistmax,1);
	aa.aahist.nfevalhist=zeros(aa.nhistmax,1);
	aa.aahist.xhist=cell(aa.nhistmax,1);
	aa.aahist.vnormshist=cell(aa.nhistmax,1);
	aa.aahist.convdatahist=cell(aa.nhistmax,1);
	aa.ithist=0;	
  end
end

if (doRestartWithFile && exist(restartFile,'file'))
  disp(['Loading restart file ...'])
  load(restartFile,'aa')
  disp(['iter0=' num2str(aa.iter0)])
end

aa.endOfRun=0;

% Top of the iteration loop.
for iter = aa.iter0:aa.itmax
  if exist('killaa','file')
	sprintf('Terminating AndAcc because of external signal at iter %d\n',iter)
	xsol=aa.x;
	return
  end
  % Apply g and compute the current residual norm.
  if ~doRestart
    [gval,~,vnorms,externalconvergence] = g(aa.x,[]);
  else
    if aa.fetchOutput
      [gval,~,vnorms,externalconvergence] = g(aa.x,1);
      aa.fetchOutput=0;
    else  
%     we save the restart first before submitting the next job
%     fetchOutput is set to 1 for the subsequent (restarted) iteration but note that 
%     the function g is called below with the fetchOutput argument set to 0
      aa.fetchOutput=1;
	  aa.iter0=iter;
	  if doRestartWithFile
%       Note: we exclude g from the restart file (and the aa struct) as it will overwrite 
%       the argument to AndAcc, which may change from one (restarted) iteration to the next 
%       depending on the problem being solved.
		save(restartFile,'aa','-v7.3')
	  end
	  g(aa.x,0); % submit job
	  break
	end  
  end  
  aa.nfeval=aa.nfeval+1;
  fval = gval - aa.x;
  res_norm = norm(fval);
  fprintf(' %d %e \n', iter, res_norm);
  if iter==0 % initialize storage of residual history
    aa.res_hist=zeros(aa.itmax+1,3+length(vnorms));
  end  
  aa.res_hist(iter+1,:) = [iter aa.nfeval vnorms res_norm];
  % Set the residual tolerance on the initial iteration.
  if iter == 0
    aa.tol = max(aa.atol,aa.rtol*res_norm);
  end
  if iter==0
    aa.ithist=aa.ithist+1;
    aa.aahist.iterhist(aa.ithist)=iter;
    aa.aahist.nfevalhist(aa.ithist)=aa.nfeval;
    aa.aahist.xhist{aa.ithist}=aa.x;
    aa.aahist.vnormshist{aa.ithist}=[vnorms res_norm];
    if (~isempty(externalconvergence)) && (isfield(externalconvergence,'convdata') && ~isempty(externalconvergence.convdata))
      aa.aahist.convdatahist{aa.ithist}=externalconvergence.convdata;
    end
    if (aa.ncheckpointfreq>0)
	  fn=['aa_spinup_checkpoint_' num2str(iter) fileSuff];
	  writeCheckpoint(fn,iter,aa.nfeval,aa.x,gval,vnorms,res_norm)
	end  
  end
  % Test for stopping.
  if (~isempty(externalconvergence)) && (isfield(externalconvergence,'converged') && ~isempty(externalconvergence.converged))
%   Check external convergence criteria  
	if externalconvergence.converged==0
	  if res_norm <= aa.tol
		fprintf('Overwriting default convergence reason as external convergence criteria have not been met: %e,  %d\n', res_norm, externalconvergence.converged); 
	  end
	  aa.endOfRun=0;  
	elseif externalconvergence.converged>0
	  fprintf('External convergence criteria have been met: %e,  %d\n', res_norm, externalconvergence.converged);
	  aa.endOfRun=1;
	  break;
	else
	  fprintf('Divergence indicated by external convergence check: %e,  %d\n', res_norm, externalconvergence.converged); 
	  aa.endOfRun=1;    
	  break;
	end
  else	
%   Use default convergence criterion  
	if res_norm <= aa.tol
	  fprintf('Terminate with residual norm = %e \n\n', res_norm);
	  aa.endOfRun=1;
	  break;
	end
  end	
  if aa.mMax == 0 || iter < aa.AAstart
    % Without acceleration, update x <- g(x) to obtain the next
    % approximate solution.
    aa.x = gval;
  else
    % Apply Anderson acceleration.
    if iter>0 && aa.mAA>0
	  if exist('flushaa','file')
		fprintf('Restart of AA is being triggered at iter %d, mAA=%d due to external file signal\n',iter,aa.mAA);
		aa.flushAA = 1;
		evalExternalCommand('rm -f flushaa');	  
	  elseif (~isempty(externalconvergence)) && (isfield(externalconvergence,'flushAA') && ~isempty(externalconvergence.flushAA))
        aa.flushAA=externalconvergence.flushAA;
        if aa.flushAA
		  fprintf('Restart of AA is being triggered at iter %d, mAA=%d due to external convergence signal\n',iter,aa.mAA);
        end 
	  elseif (aa.restartAAPeriodic>0) && mod(iter,aa.restartAAPeriodic)==0
		aa.flushAA=1;		  
		fprintf('Restart of AA is being triggered at iter %d due to periodic restart condition\n',iter);
      elseif (aa.restartAANormStagnation==1)
		itc=iter+1;
		itp=max(1,itc-1);
		rres = aa.res_hist(itp:itc,end);
		rd=rres(2)-rres(1); % change in residual norm
		if (rd>0) %|| (rd<0 && abs(rd)<aa.restartAANormDiff*res_norm)
		  aa.restart_count=aa.restart_count+1;
          if (aa.restart_count>=aa.restartAASuccessiveIters)
            aa.flushAA=1;
		    fprintf('Restart of AA is being triggered at iter %d, mAA = %d due to stagnation: %g, %g, %g\n',iter,aa.mAA,rres(1), rres(2), rd);
          end		  
		end
      end
    end
          
	if aa.flushAA
	  aa.AAstart = iter;
	  aa.mAA = 0;
	  aa.DG = [];
	  aa.R = [];
	  aa.Q = [];
      aa.restart_count=0;
	  aa.flushAA = 0;
	end  

    % Update the df vector and the DG array.
    if iter > aa.AAstart
      df = fval-aa.f_old;
      if aa.mAA < aa.mMax
        aa.DG = [aa.DG gval-aa.g_old];
      else
        aa.DG = [aa.DG(:,2:aa.mAA) gval-aa.g_old];
      end
      aa.mAA = aa.mAA + 1;
    end
    aa.f_old = fval;
    aa.g_old = gval;
    if aa.mAA == 0
      % If mAA == 0, update x <- g(x) to obtain the next approximate solution.
      aa.x = gval;
    else
      % If mAA > 0, solve the least-squares problem and update the solution
      if aa.mAA == 1
        % If mAA == 1, form the initial QR decomposition.
        aa.R(1,1) = norm(df);
        aa.Q = aa.R(1,1)\df;
      else
        % If mAA > 1, update the QR decomposition.
        if aa.mAA > aa.mMax
          % If the column dimension of Q is mMax, delete the first column and update the decomposition
          [aa.Q,aa.R] = qrdelete(aa.Q,aa.R,1);
          aa.mAA = aa.mAA - 1;
          % The following treats the qrdelete quirk described below.
          if size(aa.R,1) ~= size(aa.R,2)
            aa.Q = aa.Q(:,1:aa.mAA-1); aa.R = aa.R(1:aa.mAA-1,:);
          end
          % Explanation: If Q is not square, then qrdelete(Q,R,1) reduces the
          % column dimension of Q by 1 and the column and row
          % dimensions of R by 1. But if Q *is* square, then the
          % column dimension of Q is not reduced and only the column
          % dimension of R is reduced by one. This is to allow for
          % MATLAB's default "thick" QR decomposition, which always
          % produces a square Q.
        end
        % Now update the QR decomposition to incorporate the new
        % column.
        for j = 1:aa.mAA - 1
          aa.R(j,aa.mAA) = aa.Q(:,j)'*df;
          df = df - aa.R(j,aa.mAA)*aa.Q(:,j);
        end
        aa.R(aa.mAA,aa.mAA) = norm(df);
        aa.Q = [aa.Q,aa.R(aa.mAA,aa.mAA)\df];
      end
      if aa.droptol > 0
        % Drop residuals to improve conditioning if necessary.
        condDF = cond(aa.R);
        while condDF > aa.droptol && aa.mAA > 1
          fprintf(' cond(D) = %e, reducing mAA to %d \n', condDF, aa.mAA-1);
          [aa.Q,aa.R] = qrdelete(aa.Q,aa.R,1);
          aa.DG = aa.DG(:,2:aa.mAA);
          aa.mAA = aa.mAA - 1;
          % The following treats the qrdelete quirk described above.
          if size(aa.R,1) ~= size(aa.R,2)
            aa.Q = aa.Q(:,1:aa.mAA); aa.R = aa.R(1:aa.mAA,:);
          end
          condDF = cond(aa.R);
        end
      end
      % Solve the least-squares problem.
      gamma = aa.R\(aa.Q'*fval);
      % Update the approximate solution.
      aa.x = gval - aa.DG*gamma;
      % Apply damping if beta is a function handle or if beta > 0
      % (and beta ~= 1).
      if isa(aa.beta,'function_handle')
        aa.x = aa.x - (1-aa.beta(iter))*(fval - aa.Q*aa.R*gamma);
      else
        if aa.beta > 0 && aa.beta ~= 1
          aa.x = aa.x - (1-aa.beta)*(fval - aa.Q*aa.R*gamma);
        end
      end
    end
  end
  if iter>0 && mod(iter,aa.nhistfreq)==0
    aa.ithist=aa.ithist+1;
    aa.aahist.iterhist(aa.ithist)=iter;
    aa.aahist.nfevalhist(aa.ithist)=aa.nfeval;
    aa.aahist.xhist{aa.ithist}=aa.x;
    aa.aahist.vnormshist{aa.ithist}=[vnorms res_norm];
    if (~isempty(externalconvergence)) && (isfield(externalconvergence,'convdata') && ~isempty(externalconvergence.convdata))
	  aa.aahist.convdatahist{aa.ithist}=externalconvergence.convdata;
	end  
  end
  if (aa.ncheckpointfreq>0) && (iter>0 && mod(iter,aa.ncheckpointfreq)==0)
    fn=['aa_spinup_checkpoint_' num2str(iter) fileSuff];
    writeCheckpoint(fn,iter,aa.nfeval,aa.x,gval,vnorms,res_norm)    
%     save(['aahist_checkpoint' fileSuff],'aahist','-v7.3')
  end  
  if (iter==aa.itmax)
    aa.endOfRun=1;
  end
end

if aa.endOfRun
% Bottom of the iteration loop.
  if res_norm > aa.tol && iter == aa.itmax
	fprintf('\n Terminate after itmax = %d iterations. \n', aa.itmax);
	fprintf(' Residual norm = %e \n\n', res_norm);
  end

  aa.res_hist=aa.res_hist(1:(iter+1),:);

  % store final iteration
  aa.ithist=aa.ithist+1;
  aa.aahist.iterhist(aa.ithist)=iter;
  aa.aahist.nfevalhist(aa.ithist)=aa.nfeval;
  aa.aahist.xhist{aa.ithist}=aa.x;
  aa.aahist.vnormshist{aa.ithist}=[vnorms res_norm];
  if (~isempty(externalconvergence)) && (isfield(externalconvergence,'convdata') && ~isempty(externalconvergence.convdata))
	aa.aahist.convdatahist{aa.ithist}=externalconvergence.convdata;
  end

  if (aa.ncheckpointfreq>0)
	fn=['aa_spinup_checkpoint_final' fileSuff];
	writeCheckpoint(fn,iter,aa.nfeval,aa.x,gval,vnorms,res_norm)  
  end

% Trim off unused space  
  aa.aahist.iterhist=aa.aahist.iterhist(1:aa.ithist);
  aa.aahist.nfevalhist=aa.aahist.nfevalhist(1:aa.ithist);
  aa.aahist.xhist=aa.aahist.xhist(1:aa.ithist);
  aa.aahist.vnormshist=aa.aahist.vnormshist(1:aa.ithist);
  aa.aahist.convdatahist=aa.aahist.convdatahist(1:aa.ithist);

  xsol=aa.x;
else
  xsol=[];
end

end % main function

function writeCheckpoint(fn,iter,nfeval,x,gval,vnorms,res_norm)

save(fn,'iter','nfeval','x','gval','vnorms','res_norm')

end
