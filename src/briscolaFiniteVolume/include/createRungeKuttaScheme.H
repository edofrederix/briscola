autoPtr<RungeKuttaScheme> rkPtr(RungeKuttaScheme::New(fvMsh));

RungeKuttaScheme& rk = rkPtr();
