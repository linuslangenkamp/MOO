model include_test
  Real x(start=1, fixed=true);
equation
  der(x) = x;
end include_test;
