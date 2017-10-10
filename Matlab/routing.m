function [qRout] = routing(qSim, maxbas)
%  ROUTING performs routing using a triangular function
%  
%  qSim = Resulting dishcarge [m³/s]
%  maxbas = (integer) Transformation function parameter
%  qRout
%  
%  
%  
%  
  maxbas = max(int8(maxbas), 2);
	nSim = length(qSim);
	c = mytrimf(1:maxbas, [0, (maxbas+1)/2, maxbas+1]); %(Seibert,1997)
	c = c/sum(c) % vector of normalized coefficients - (1,MAXBAS)
	
  qRout = qSim;
	for t = maxbas:nSim
		qRout(t) = c * qSim(t-maxbas+1:t); % (Seibert,1997)
	end
  
  function f = mytrimf(x,param)
    % implements triangular-shaped membership function
    % (available in Matlab Fuzzy Logic Toolbox as 'trimf')

    f = zeros(size(x)) ; 
    idx = (x>param(1)) & (x<=param(2)) ;
    f(idx) = (x(idx)-param(1))/(param(2)-param(1)) ;
    idx = (x>param(2)) & (x<=param(3)) ;
    f(idx) = (param(3)-x(idx))/(param(3)-param(2)) ;
  end
end