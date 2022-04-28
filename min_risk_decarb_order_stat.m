function [str,abs_vol,rel_vol,RC,RC_rel,red] = min_risk_decarb_order_stat(sigs,corr,bench,CEs,Y,k,varargin)
   % Inputs:
   % sigs = vector of volatilities
   % corr = matrix of correlation
   % bench = benchmark structure of portfolio
   % CEs = carbon emissions/intensities
   % Y = vector of (Net) revenues
   
   % ---------------------
   % Outputs
   % str = structure of the portfolio
   % abs_vol = absolute volatility 
   % rel_vol = relative volatility
   % RC = risk contributions
   % RC_rel = relative risk contributions
   
   p = inputParser;
   options = {'relative','absolute'};
   %Above it is referred as relative (tracking error vol) or absolute vol
   %(x*Sigma*x')^(1/2)
   carbon_options = {'WACI','emissions'};
   addOptional(p,'type','relative',@(x) any(validatestring(x,options)));
   addOptional(p,'carbon','WACI',@(x) any(validatestring(x,carbon_options)));
   parse(p,varargin{:});
   opt = p.Results.type;
   %opt_carb = p.Results.carbon;
   cov_mat = sigs'.*sigs.*corr;
   carb_intensities = CEs./Y;
   if strcmpi(opt,'relative')==1
      func = @(x) 1/2*(x-bench)*cov_mat*(x-bench)';
   else
      func = @(x) 1/2*x*cov_mat*x';
   end
   n = length(sigs);
   pos = arrayfun(@(i) find(carb_intensities(i)==sort(carb_intensities,'desc')),1:k);
   z = zeros(1,n);
   z(pos) = ones(1,k);
   Aeq = [z;ones(1,n)];beq = [0;1];
   x0 = ones(1,n)/n;
   str = fmincon(func,x0,[],[],Aeq,beq,zeros(1,n),ones(1,n));
   abs_vol = sqrt(str*cov_mat*str');
   if isempty(bench)==0
      rel_vol = sqrt((str-bench)*cov_mat*(str-bench)');
      RC = (str-bench).*(cov_mat * (str-bench)')'/rel_vol;
      RC_rel = RC/sum(RC);
   else
      rel_vol = [];
      RC = str.*(cov_mat*str')'/abs_vol;
      RC_rel = RC/sum(RC);
   end
   red = 1-str*carb_intensities'/(bench*carb_intensities');
end