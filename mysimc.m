function [cumRV,fluidLevel] = mysimc(n,mycont,myObj,lam)
%mysimc simulates state transitions from an arbitrary T(t) matrix.
% the simulation approach is based on the rejection method.  See 
% Introduction to Probability Models, Ross, Ch. 11.2.2 ed 14
% amended from mysimc to add entries for time when zero is hit.
%  inputs
%   n is number of jumps
%   mycont - current state of system [phase, fluid, time]
%   myObj - CyclicFluidQueue object
%   lam value greater than or equal to any of the transition rates
%  outputs
%   cumRV - array of transition times and phases of background process
%   fluidLevel - counts times when fluid level is zero

%start_time = cputime;

T=@(t)myObj.Toft(t);
% goofed up n; making n input again
%n = myObj.Mesh();

states = max(size(T(0)));

% scale based on flow rates
cvec = myObj.flowRates();
c=max(abs(cvec));

% get lam: want constant larger than largest transition rate times
% fluid flow rate.  lam is input as larger than largest transition 
% rate of background process

lam = c*lam*1.01;
stateTrans = zeros(states,states-1);
copyarray = 1:states;
% stateTrans is a statesx(states-1) array with diagonal stripped out.  Rows
% consist of natural numbers\0 without diagonal entry
for j=1:states
    stateTrans(j,:) =copyarray([1:j-1, j+1:end]);
end

% initialize random vector
cumRV = zeros(n,2);
% we start in state given by mycont
currState = mycont(1);
% cumRV = [time state flowrate]
cumRV(1,2)=currState;
mytimer = mycont(3);
% fprintf('mysimb starts with %7.2f  days\n for model %2.0f\n',...
%     mytimer,myObj.Ver());

%***? start at j=2? *** was 1
for j=2:n
    mytimer=mytimer-log(rand())/lam;
    newrand = rand();
    discard = T(mytimer);
    newtest = abs(discard(currState,currState))/lam;
    % we test whether U>lambda(t)/lambda.  If it is we update time
    while newrand>newtest
        newrand = rand();
        mytimer=mytimer-log(rand())/lam;
        discard = T(mytimer);
        newtest = abs(discard(currState,currState))/lam;       
    end
    % we create a probability vector from the rate matrix using currState
    % as the row and all states excluding diag for the probabilities.
    %probarray = discard(currState,[1:currState-1, currState+1:end])/lam;
    % divide by |T(currState,currState)| to get discrete probability mass
    % function.
    denom = abs(discard(currState,currState));
    probarray = discard(currState,[1:currState-1, currState+1:end])/denom;
    %sample from the probability mass function to get the new state
    currState = datasample(stateTrans(currState,:),1,'Weights',probarray);
    % update the currRV vector with the new state
    cumRV(j,:) = [mytimer currState];
end

% cycle through the currRV vector to get fluid levels
fluidLevel = zeros(n,1);
% initial fluid level given by mycont
fluidLevel(1) = mycont(2);
% *** maybe a while loop since entries are added
% or put number of zeros in the count ***
steps = 2;
j=2;
while steps<n
    newFluidLevel = fluidLevel(j-1)...
        +(cumRV(j,1)-cumRV(j-1,1))*cvec(cumRV(j-1,2));
    fluidLevel(j) = max(newFluidLevel,0);
    steps = steps+1;
    if newFluidLevel<0
        timeinc = -fluidLevel(j-1)/cvec(cumRV(j-1,2));
        timewhenhit = timeinc+cumRV(j-1,1);
        cumRV = [cumRV(1:j-1,:); ...
            [timewhenhit cumRV(j-1,2)]; ...
            cumRV(j:end,:)];
        fluidLevel = [fluidLevel(1:j-1); 0; fluidLevel(j:end)];
        j=j+1;
    end
    j=j+1;
end

maxdays = floor(cumRV(end,1));

% end_time = cputime;
% sim_time = end_time - start_time;
% fprintf('mysimc required %7.2f minutes to simulate %7.2f days\n for model %2.0f\n',...
%     sim_time/60.,maxdays,myObj.Ver());

end