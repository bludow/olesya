double integral_gauss_phi0(double r_max, double b) {

	double a = 0;
	double b = ...;
	double sum  = 0.;
	double V = potential();
	double E = ...;
	//3 - порядок метода 
	for (double i = a; i < 2; ++i) {
			double r1 = ((a + b) / 2.) - (b - a) / (2. * sqrt(3.));
			double dr1 = ... ;
			double fun1 = ((b / r1 * r1) / sqrt((1 - b * b) / r1 * r1 - (V / E))) * dr1;


			double r2 = ((a + b) / 2.) + (b - a) / (2. * sqrt(3.));
			double dr2 = ... ;
			double fun2 = ((b / r2 * r2) / sqrt((1 - b * b) / r2 * r2 - (V / E))) * dr2;

			double r = fun1 + fun2;
			sum += ((b - a) / 2.) * (((b / r1 * r1) / sqrt((1 - b * b) / r1 * r1 - (V / E))) * dr1);
	}

}
