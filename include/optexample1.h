#include "cfcmINLP.h"
#include <math.h>

class myExample1 : public OptInterface{

public:
	myExample1();
	~myExample1() {}

	void eval_f(uint n, double *x, double& obj_value);

	void eval_g(uint n, double* x, uint m, double* g);   

};