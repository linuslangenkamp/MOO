model include_test
  Real x1(start = 1, fixed=true) "kerogen";
  Real x2(start = 0, fixed=true) "pyrolytic bitumen";
  Real x3(start = 0, fixed=true) "oil & gas";
  Real x4(start = 0, fixed=true) "organic carbon";
  Real k1;
  Real k2;
  Real k3;
  Real k4;
  Real k5;
  Real final_c(max=0.1) = x2 annotation(isFinalConstraint=true);
  input Real u(start=1, min=650/700, max = 748.15/700);
  output Real cost1 = -x2 annotation(isMayer=true);
equation
  k1 = exp(8.86 - 20300 / 1.9872 / (700 * u));
  k2 = exp(24.25 - 37400 / 1.9872 / (700 * u));
  k3 = exp(23.67 - 33800 / 1.9872 / (700 * u));
  k4 = exp(18.75 - 28200 / 1.9872 / (700 * u));
  k5 = exp(20.70 - 31000 / 1.9872 / (700 * u));
  der(x1) = (-k1 * x1) - (k3 + k4 + k5) * x1 * x2;
  der(x2) = k1 * x1 - k2 * x2 + k3 * x1 * x2;
  der(x3) = k2 * x2 + k4 * x1 * x2;
  der(x4) = k5 * x1 * x2;
  annotation(
    experiment(StartTime = 0, StopTime = 8, Tolerance = 1e-10, Interval = 0.05),
    __OpenModelica_simulationFlags(s = "optimization", optimizerNP = "3", ipopt_init="CONST",  lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS,LOG_IPOPT"),
    __OpenModelica_commandLineOptions = "+g=Optimica --matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
end include_test;
/*
model include_test
  Real x1(start = 1, fixed=true);
  Real x2(start = 0, fixed=true);
  input Real u(min=0, max=5);
  output Real cost_m = -x2 annotation(isMayer=true);
equation
  der(x1) = -(u + u^2 / 2) * x1;
  der(x2) = u * x1;
  annotation(
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-12, Interval = 0.05),
    __OpenModelica_simulationFlags(s = "optimization", optimizerNP = "3", lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS"),
    __OpenModelica_commandLineOptions = "+g=Optimica --matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
end include_test;
*/