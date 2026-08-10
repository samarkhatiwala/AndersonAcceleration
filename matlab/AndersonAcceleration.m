% This file implements Anderson Acceleration (AA). It is specifically designed for the 
% computation of equilibrium solutions of periodically-forced ocean and land biogeochemical 
% models (see Khatiwala (2023; https://doi.org/10.1029/2022MS003447) and Khatiwala (2024; 
% https://doi.org/10.1126/sciadv.adn2839)) but can be used for any fixed point iteration 
% problem. It is based on AndAcc.m originally written by Homer F. Walker (listed in H. Walker, 
% Anderson acceleration: Algorithms and implementations, Worcester Polytechnic Institute 
% Mathematical Sciences Department Research Report MS-6-15-50, June, 2011). The original 
% code has been reimplemented as an object oriented Matlab handle class with extensive 
% modifications as per Khatiwala (2023, 2024). 
% To use, first instantiate an AndersonAcceleration object before calling the solve method:
% aa=AndersonAcceleration;
% ...=aa.solve(...)

classdef AndersonAcceleration < handle   
  properties
    x
    y = []
    x_prev
    mMax = 10
    itmax = 100
    atol = 1.e-10
    rtol = 1.e-10
    droptol = 1.e10
    beta = 1.0
    AAstart = 0
    restartAANormStagnation = 0
    restartAANormDiff = 0.01
    restartAASuccessiveIters = 2
    restartAAPeriodic = 0
    DG = [] % Storage of g-value differences.
    Q  = []
    R = []
    f_old = 0
    g_old = 0
    DGy = [] % Storage of g-value differences.
    fy_old = 0
    gy_old = 0
    f_old_prev
    g_old_prev
    restart_count = 0
    flushAA = 0
    mAA = 0
    nfeval = 0
    iter0 = 0
    tol = 0
    fetchOutput = 0
    endOfRun = 0
    res_hist = [] % Storage of residual history (a matrix with columns: 
%                   	[iter nfeval vnorms res_norm]). This will be overwriiten below.
    nhistfreq = -1
    ncheckpointfreq = -1
    nhistmax
    aahist
    ithist	
  end    
  methods
    function [xsol,iter,ysol,converged] = solve(aa,g,x,AAparams,histParams,doRestart,fileSuff,y)
%                                                  1 2    3         4         5        6      7
    % This is the main solver function.
    % Input arguments:
    %  Required:
    %   g: function handle for the fixed point map (see below for arguments)
    %   x: initial guess vector
    %  Optional:
    %   AAparams: struct with fields defining algorithm parameters (see defaults below)
    %   histParams: struct with fields defining history parameters (see defaults below)
    %   doRestart: flag to indicate whether to restart or not
    %   fileSuff: string to use for tagging checkpointing files
    %   y: initial guess if using an "x/y split"
    % Output:
    %   xsol: final solution
    %   iter: final iteration number
    %   ysol: final y solution or empty vector
    %   converged: flag indicating convergence

    % Definition of the g function:
    %  [gx,gy,vnorms,externalconvergence,v,gv]=g(x,y,fetchOutput)
    %  Inputs:
    %    x: iterate vector at which to evaluate g
    %    fetchOutput: flag (0, 1 or []) to indicate 'mode' in which to evaluate g
    %  Outputs:
    %    gx: output vector from evaluating g
    %    gv: dimensional (rescaled) output (not used)
    %    vnorms: dimensional residual norm (a ROW vector if this is a multi-tracer problem, otherwise a scalar)
    %    externalconvergence: struct of external convergence data

    % Note: when this function is invokes as obj.solve(...), the first argument passed is obj 
    %       and is included in nargin
    numInputargs = nargin - 1;
    if numInputargs<2
      error('AndersonAcceleration.solve requires at least two arguments.')
    end

    % Figure out whether we are doing restarts
    if numInputargs<5 || isempty(doRestart)
      doRestart=0;
    end
    if ~ismember(doRestart,[0 1])
      error('ERROR: doRestart must be 0 or 1!')
    end  

    % Instantiate a default histParams struct if one is not provided
    if numInputargs<4 || isempty(histParams)
      histParams=struct();
    end

    trimHistory=1;
    if numInputargs>3 && ~isempty(histParams)
      if isfield(histParams,'trimHistory') && ~isempty(histParams.trimHistory)
        trimHistory=histParams.trimHistory;
      end
    end

    if numInputargs<6 || isempty(fileSuff)
      fileSuff='';
    end

    % Instantiate a default AAparams struct if one is not provided
    if numInputargs<3 || isempty(AAparams)
      AAparams=struct();
      AAparams.killFile='killaa';
      AAparams.wallTime=[];
    end

    if isempty(aa.x)
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
      % Note: we exclude g from the aa object as it may change from one (restarted) 
      % iteration to the next depending on the problem being solved.

      aa.x = x;
      if numInputargs>6 && ~isempty(y)
        aa.y = y;
      end
      aa.mMax = min(50, size(aa.x,1));

      if numInputargs>2 && ~isempty(AAparams)
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
        if isfield(AAparams,'killFile') && isempty(AAparams.killFile)
          AAparams.killFile='killaa';
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

    %   Set history defaults and overwrite with values from histParams
    %   nhistfreq = frequency with which to save history
    %   ncheckpointfreq = frequency with which to write checkpoint (if < 0 
    %                     then no checkpoints are written; this is the default)
      if numInputargs>3 && ~isempty(histParams)
        if isfield(histParams,'nhistfreq') && ~isempty(histParams.nhistfreq)
          aa.nhistfreq=histParams.nhistfreq;
        end
        if isfield(histParams,'ncheckpointfreq') && ~isempty(histParams.ncheckpointfreq)
          aa.ncheckpointfreq=histParams.ncheckpointfreq;
        end
        if isfield(histParams,'nhistmax') && ~isempty(histParams.nhistmax)
          if (histParams.nhistmax<(ceil(aa.itmax/aa.nhistfreq)+5))
            error('ERROR: histParams.nhistmax is too small!')
          end  
          aa.nhistmax=histParams.nhistmax;
        else
          aa.nhistmax=ceil(aa.itmax/aa.nhistfreq)+5; % space to store history	  
        end
      end
      % SPK this needs to be fixed for when I don't want to save history
    % initialize history
      if aa.nhistfreq>0
        aa.aahist.iterhist=zeros(aa.nhistmax,1);
        aa.aahist.nfevalhist=zeros(aa.nhistmax,1);
        aa.aahist.xhist=cell(aa.nhistmax,1);
        aa.aahist.yhist=cell(aa.nhistmax,1);
        aa.aahist.vnormshist=cell(aa.nhistmax,1);
        aa.aahist.convdatahist=cell(aa.nhistmax,1);
        aa.ithist=0;
      end
    else
      % This is a restart so check if any parameters have changed
      if aa.iter0>0
        aaHasChanged=0;
        if numInputargs>2 && ~isempty(AAparams)
          if isfield(AAparams,'itmax') && ~isempty(AAparams.itmax)
            if aa.itmax~=AAparams.itmax
              disp('WARNING: itmax has changed!')
              if aa.itmax>AAparams.itmax
                error('ERROR: new itmax < old itmax!')
              else
                aaHasChanged=1;
                disp('Resizing res_hist')
                tmp=aa.res_hist;   
                aa.res_hist=zeros(AAparams.itmax+1,size(tmp,2));
                aa.res_hist(1:aa.itmax+1,:)=tmp;
                aa.itmax=AAparams.itmax;
                clear tmp
              end
            end  
          end
        end
        if size(aa.res_hist,1)~=(aa.itmax+1)
          aaHasChanged=1;
          disp('Resizing res_hist')
          tmp=aa.res_hist;
          aa.res_hist=zeros(aa.itmax+1,size(tmp,2));
          aa.res_hist(1:size(tmp,1),:)=tmp;
          clear tmp
        end    
        if aaHasChanged
          if aa.nhistfreq>0
            if (aa.nhistmax<(ceil(aa.itmax/aa.nhistfreq)+5))
              error('ERROR: aa.nhistmax is too small!')
            end
          end
        end
      end	
    end

    % Finished initialization

    aa.endOfRun=0;
    converged=0;

    % Top of the iteration loop.
    for iter = aa.iter0:aa.itmax
      if isfield(AAparams,'wallTime') && ~isempty(AAparams.wallTime)
        if (toc(AAparams.startTime)>AAparams.wallTime)
          sprintf('Terminating AndAcc because of wall time limit at iter %d\n',iter)
          aa.iter0=iter;
          xsol=aa.x;
          ysol=aa.y;
          return
        end
      end
      if exist(AAparams.killFile,'file')
       sprintf('Terminating AndAcc because of external signal at iter %d\n',iter)
       aa.iter0=iter;
       xsol=aa.x;
       ysol=aa.y;
       return
      end
      % Apply g and compute the current residual norm.
      if ~doRestart
        [gval,gyval,vnorms,externalconvergence,v,gv] = g(aa.x,aa.y,[]);
      else
        if aa.fetchOutput
          [gval,gyval,vnorms,externalconvergence,v,gv] = g(aa.x,aa.y,1);
          aa.fetchOutput=0;
        else  
    %     we save the restart information first before submitting the next job
    %     fetchOutput is set to 1 for the subsequent (restarted) iteration but note that 
    %     the function g is called below with the fetchOutput argument set to 0
          aa.fetchOutput=1;
          aa.iter0=iter;
          g(aa.x,aa.y,0); % submit job
          break
        end  
      end  
      aa.nfeval=aa.nfeval+1;
      fval = gval - aa.x;
      if ~isempty(aa.y)
        fyval = gyval - aa.y;
      end
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
        if aa.nhistfreq>0
          aa.ithist=aa.ithist+1;
          aa.aahist.iterhist(aa.ithist)=iter;
          aa.aahist.nfevalhist(aa.ithist)=aa.nfeval;
          aa.aahist.xhist{aa.ithist}=aa.x;
          aa.aahist.yhist{aa.ithist}=aa.y;
          aa.aahist.vnormshist{aa.ithist}=[vnorms res_norm];
          if (~isempty(externalconvergence)) && (isfield(externalconvergence,'convdata') && ~isempty(externalconvergence.convdata))
            aa.aahist.convdatahist{aa.ithist}=externalconvergence.convdata;
          end
        end
        if (aa.ncheckpointfreq>0)
          fn=['aa_spinup_checkpoint_' num2str(iter) fileSuff];
          writeCheckpoint(fn,iter,aa.nfeval,aa.x,v,gval,gv,vnorms,res_norm,aa.y)
        end  
      end
      % Test for stopping.
      if (~isempty(externalconvergence)) && (isfield(externalconvergence,'converged') && ~isempty(externalconvergence.converged))
        %   Check external convergence criteria  
        if externalconvergence.converged==0
          if res_norm <= aa.tol
         fprintf('Overwriting default convergence reason as external convergence criteria have not been met: %e,  %d\n', res_norm, externalconvergence.converged); 
          end
        elseif externalconvergence.converged>0
          fprintf('External convergence criteria have been met: %e,  %d\n', res_norm, externalconvergence.converged);
          aa.endOfRun=1;
          converged=externalconvergence.converged;
          break;
        else
          fprintf('Divergence indicated by external convergence check: %e,  %d\n', res_norm, externalconvergence.converged); 
          aa.endOfRun=1;
          converged=externalconvergence.converged;
          break;
        end
      else	
       %   Use default convergence criterion  
        if res_norm <= aa.tol
          fprintf('Terminate with residual norm = %e \n\n', res_norm);
          aa.endOfRun=1;
          converged=1;
          break;
        end
      end	
      if aa.mMax == 0 || iter < aa.AAstart
        % Without acceleration, update x <- g(x) to obtain the next
        % approximate solution.
        aa.x = gval;
        if ~isempty(aa.y)
          aa.y = gyval;
        end
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
            else
              aa.restart_count=0;
            end  
          end
        end
          
        if aa.flushAA
          aa.AAstart = iter;
          aa.mAA = 0;
          aa.DG = [];
          aa.R = [];
          aa.Q = [];
          if ~isempty(aa.y)
            aa.DGy = [];
          end	  
          aa.restart_count=0;
          aa.flushAA = 0;
        end  

        % Update the df vector and the DG array.
        if iter > aa.AAstart
          df = fval-aa.f_old;
          if aa.mAA < aa.mMax
            aa.DG = [aa.DG gval-aa.g_old];
            if ~isempty(aa.y)
              aa.DGy = [aa.DGy gyval-aa.gy_old];
            end
          else
            aa.DG = [aa.DG(:,2:aa.mAA) gval-aa.g_old];
            if ~isempty(aa.y)
              aa.DGy = [aa.DGy(:,2:aa.mAA) gyval-aa.gy_old];
            end      
          end
          aa.mAA = aa.mAA + 1;
        end
        aa.f_old = fval;
        aa.g_old = gval;
        if ~isempty(aa.y)
          aa.fy_old = fyval;
          aa.gy_old = gyval;
        end	
        if aa.mAA == 0
          % If mAA == 0, update x <- g(x) to obtain the next approximate solution.
          if ~aa.flushAA
            aa.x = gval;
            if ~isempty(aa.y)
              aa.y = gyval;
            end
          else
    %         fprintf('Backtracking to previous solution\n');
    %   	    aa.x=aa.x_prev;
    %         aa.f_old=aa.f_old_prev;
    %         aa.g_old=aa.g_old_prev;
    %   	    aa.flushAA=0;
          end  
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
              if ~isempty(aa.y)
                aa.DGy = aa.DGy(:,2:aa.mAA);
              end          
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
          if ~isempty(aa.y)
            aa.y = gyval - aa.DGy*gamma;
          end  
          % Apply damping if beta is a function handle or if beta > 0
          % (and beta ~= 1).
          if isa(aa.beta,'function_handle')
            aa.x = aa.x - (1-aa.beta(iter))*(fval - aa.Q*aa.R*gamma);
            if ~isempty(aa.y)
              error('This is not supported!')
            end  
          else
            if aa.beta > 0 && aa.beta ~= 1
              aa.x = aa.x - (1-aa.beta)*(fval - aa.Q*aa.R*gamma);
              if ~isempty(aa.y)
                error('This is not supported!')
              end
            end
          end
        end
        aa.x_prev=aa.x;
        aa.f_old_prev=aa.f_old;
        aa.g_old_prev=aa.g_old;
      end % end of apply AA block
      if aa.nhistfreq>0
        if iter>0 && mod(iter,aa.nhistfreq)==0
          aa.ithist=aa.ithist+1;
          aa.aahist.iterhist(aa.ithist)=iter;
          aa.aahist.nfevalhist(aa.ithist)=aa.nfeval;
          aa.aahist.xhist{aa.ithist}=aa.x;
          aa.aahist.yhist{aa.ithist}=aa.y;
          aa.aahist.vnormshist{aa.ithist}=[vnorms res_norm];
          if (~isempty(externalconvergence)) && (isfield(externalconvergence,'convdata') && ~isempty(externalconvergence.convdata))
            aa.aahist.convdatahist{aa.ithist}=externalconvergence.convdata;
          end  
        end
      end
      if (aa.ncheckpointfreq>0) && (iter>0 && mod(iter,aa.ncheckpointfreq)==0)
        fn=['aa_spinup_checkpoint_' num2str(iter) fileSuff];
        writeCheckpoint(fn,iter,aa.nfeval,aa.x,v,gval,gv,vnorms,res_norm,aa.y)  
      end
      if isfield(histParams,'checkpointingFunction') && isa(histParams.checkpointingFunction,'function_handle')
        histParams.checkpointingFunction(iter,aa)
      end    
      if (iter==aa.itmax)
        aa.endOfRun=1;
      end
    end % end of iteration loop

    if aa.endOfRun
    % Bottom of the iteration loop.
      if res_norm > aa.tol && iter == aa.itmax
        fprintf('Terminated after itmax = %d iterations. \n', aa.itmax);
       fprintf(' Residual norm = %e \n', res_norm);
      end

    %   aa.res_hist=aa.res_hist(1:(iter+1),:);

      % store final iteration
      if aa.nhistfreq>0
        aa.ithist=aa.ithist+1;
        aa.aahist.iterhist(aa.ithist)=iter;
        aa.aahist.nfevalhist(aa.ithist)=aa.nfeval;
        aa.aahist.xhist{aa.ithist}=aa.x;
        aa.aahist.yhist{aa.ithist}=aa.y;
        aa.aahist.vnormshist{aa.ithist}=[vnorms res_norm];
        if (~isempty(externalconvergence)) && (isfield(externalconvergence,'convdata') && ~isempty(externalconvergence.convdata))
          aa.aahist.convdatahist{aa.ithist}=externalconvergence.convdata;
        end
      end
     
    % This is to allow continuation of a previous solve
      aa.iter0=iter+1;
  
      if (aa.ncheckpointfreq>0)
        fn=['aa_spinup_checkpoint_final' fileSuff];
        writeCheckpoint(fn,iter,aa.nfeval,aa.x,v,gval,gv,vnorms,res_norm,aa.y)
      end

    % Trim off unused space  
      if trimHistory
        aa.res_hist=aa.res_hist(1:(iter+1),:);
        if aa.nhistfreq>0
          aa.aahist.iterhist=aa.aahist.iterhist(1:aa.ithist);
          aa.aahist.nfevalhist=aa.aahist.nfevalhist(1:aa.ithist);
          aa.aahist.xhist=aa.aahist.xhist(1:aa.ithist);
          aa.aahist.yhist=aa.aahist.yhist(1:aa.ithist);
          aa.aahist.vnormshist=aa.aahist.vnormshist(1:aa.ithist);
          aa.aahist.convdatahist=aa.aahist.convdatahist(1:aa.ithist);
        end
      end

      xsol=aa.x;
      ysol=aa.y;
    else
      xsol=[];
      ysol=[];
    end

    end % main function

    function writeCheckpoint(fn,iter,nfeval,x,v,gval,gv,vnorms,res_norm,y)
      save(fn,'iter','nfeval','x','v','gval','gv','vnorms','res_norm','y')
    end

  end % end methods
end % end classdef
