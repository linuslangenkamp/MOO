/*model A
  Real x1(start = 1, fixed=true) "kerogen";
  Real x2(start = 0, fixed=true) "pyrolytic bitumen";
  Real k1;
  Real k2;
  Real k3;
  Real k4;
  Real k5;
  // Real z;
  Real CONSTR(min=-1.25, max=1000) = x1 + u annotation(isConstraint=true);
  output Real FINALCONSTR(min=-1.25, max=1.25) = x1 * x2 annotation(isFinalConstraint=true);
  input Real u(start=700, min=650, max = 748.15);
  output Real cost_m = -x1 annotation(isMayer=true);
  output Real cost_l = -x2 + u annotation(isLagrange=true);
equation
  // x1 + u^2 = z;
  k1 = exp( 8.86 - 20300 / 1.9872 / u);
  k2 = exp(24.25 - 37400 / 1.9872 / u);
  k3 = exp(23.67 - 33800 / 1.9872 / u);
  k4 = exp(18.75 - 28200 / 1.9872 / u);
  k5 = exp(20.70 - 31000 / 1.9872 / u);
  der(x1) = (-k1 * x1) - (k3 + k4 + k5) * x1 * x2;
  der(x2) = k1 * x1 - k2 * x2 + k3 * x1 * x2 ;
  annotation(
    experiment(StartTime = 0, StopTime = 16, Tolerance = 1e-12, Interval = 0.5),
    __OpenModelica_simulationFlags(solver = "optimization", optimizerNP = "3", lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS"),
    __OpenModelica_commandLineOptions = "+g=Optimica --matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
end A;*/

model include_test
  Real x1(start = 1, fixed=true);
  Real x2(start = 0, fixed=true);
  Real g(min=0, max=1.5) = x1 + u annotation(isConstraint=true);
  input Real u(min=0, max=5);
  output Real cost_l = -u * x1 annotation(isLagrange=true);
equation
  der(x1) = -(u + u^2 / 2) * x1;
  der(x2) = u * x1;
  annotation(
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-12, Interval = 0.5),
    __OpenModelica_simulationFlags(solver = "optimization", optimizerNP = "3", lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS"),
    __OpenModelica_commandLineOptions = "+g=Optimica --matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian");
end include_test;