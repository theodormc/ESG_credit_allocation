function [str,abs_vol,rel_vol,RC,RC_rel,ret] = min_risk_decarb(sigs,corr,bench,CEs,Y,R,HCIS,varargin)
   % Inputs:
   % sigs = vector of volatilities
   % corr = matrix of correlation
   % bench = benchmark structure of portfolio
   % CEs = carbon emissions/intensities
   % Y = vector of (Net) revenues
   % R = target reduction rate
   % HCIS (High Climate Impact Sector): vector of binary values 0/1 where 0
   % represents that the stock (entity) i does not belong to HCIS and 1 if
   % entity i does belong to HCIS. 
   % varargin : contains information on the:
   % 1) type of risk that has to be minimized:   
   % It can be absolute (absolute vol) and relative (tracking error
   % volatility)
   % the type of constrain on carbon Emissions: Weighted Average Carbon Intensity
   % or Absolute Carbon Emission ('WACI'/'emission'). 
   
   % Outputs:
   % str = structure of the portfolio
   % abs_vol = absolute volatility
   % rel_vol = tracking error volatility 
   % RC = risk contributions 
   % RC_rel = relative (normalized) risk contributions
   
   p = inputParser;
   options = {'relative','absolute'};
   %Above it is referred as relative (tracking error vol) or absolute vol
   %(x*Sigma*x')^(1/2)
   carbon_options = {'WACI','emissions'};
   addOptional(p,'type','relative',@(x) any(validatestring(x,options)));
   addOptional(p,'carbon','WACI',@(x) any(validatestring(x,carbon_options)));
   addOptional(p,'returns',[],@(x) isvector(x));
   parse(p,varargin{:});
   opt = p.Results.type;
   opt_carb = p.Results.carbon;
   rets = p.Results.returns;
   cov_mat = sigs'.*sigs.*corr;
   carb_intensities = CEs./Y;
   if strcmpi(opt,'relative')==1
      func = @(x) (x-bench)*cov_mat*(x-bench)';
   else
      func = @(x) x*cov_mat*x';
   end
   n = length(sigs);
   % C and D represent restrictions. Work with carbon intensities
   if strcmpi(opt_carb,'WACI')==1
      C = carb_intensities;D = (1-R)*bench*carb_intensities';
   else
      C = CEs./bench ;
      D = (1-R) * sum(CEs);
   end
   if isempty(HCIS)==1
      str = fmincon(func,ones(1,n)/n,C,D,ones(1,n),1,zeros(1,n),ones(1,n));
   else
      HCIS_bench_exposure = bench*HCIS';
      C = [C;-HCIS];
      D = [D;-HCIS_bench_exposure];
      str = fmincon(func,ones(1,n)/n,C,D,ones(1,n),1,zeros(1,n),ones(1,n));
   end
   
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
   
   if isempty(rets)==1
       ret = [];
   else
       ret = rets*str';
   end
end

