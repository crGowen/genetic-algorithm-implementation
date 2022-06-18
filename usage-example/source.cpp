// UsesLib.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <genetic-alg.h>

// NOTE: the algorithm is not necessarily going to produce the most 'elegant' solution, it will simply go for a solution that works well, for instance: in this example program,
// the obvious solution to this fitness function (produce 15 from x * y) is x=3, y=5, however a result of something like
// x = 2.6660446, y = 5.6263123 is just as likely as they both produce 15 within the floating point margin of error

// fitness function must always be defined by the user as a double with a single (unsigned __int32*) parameter
double EvalSolut(unsigned __int32* g) {
	/*
	suitable for lower level of precision
	float x = *reinterpret_cast<float*>(&g[0]);
	float y = *reinterpret_cast<float*>(&g[1]);

	BUT FOR DOUBLE PRECISION WE CAN DO THIS INSTEAD:
	*/

	double x = *reinterpret_cast<double*>(&g[0]);
	double y = *reinterpret_cast<double*>(&g[2]);


	// reinterpret the 8 bytes of 2  __int32 types as the 8 bytes of a double


	double res = 15.0 - (x * y); // use that to find the difference from 15

	/* force the floating point to focus on a narrow range with high precision
	(if we didn't do this, the floating point would likely just reduce its own precision and
	"win" by cancelling out large exponents, which is lazy and useless)
	*/

	double tempA = (abs(x) + 1) / 10;
	double tempB = 1 / abs(x);

	// if either x or y genes move to being less than 1 or more than 10, punish this solution heavily by magnifying the res value:
	// we want to keep both x and y within the range 1-10, otherwise it will
	// "cheat" by reducing its own precision through the use of high or low exponents
	if (abs(x) < 1 || abs(x) > 10) {
		tempA = (abs(x) + 1) / 10;
		tempB = 1 / abs(x);
		if (tempA > tempB)
			res = res * tempA;
		else
			res = res * tempB;
	}

	// do the same for y
	if (abs(y) < 1 || abs(y) > 10) {
		tempA = (abs(y) + 1) / 10;
		tempB = 1 / abs(y);
		if (tempA > tempB)
			res = res * tempA;
		else
			res = res * tempB;
	}

	// return that difference as negative (because the gen alg logic always wants to optimise the fitness value to be as great as possible, which means
	// for negative numbers the closer to 0, the greater the value
	return (-1) * abs(res);
}

int main()
{
	/*
	create algorithm with:

	a population of 300,000
	80 generations
	4 genes (which are x and y in our fitness evaluation that we defined above)
	we point to our defined fitness evaluation "EvalSolut",
	tournament selection group size of 6
	70,000/100,000 (70%) mutation rate
	50,000/100,000 (50%) crossover rate
	*/

	GeneticAlgorithm::Initialise(300000, 80, 4, &EvalSolut, 6, 70000, 50000);

	//run the algorithm with output as true (output means it will print the result of each generation, eg: 'Generation 42 completed. Best fitness = -4.21e-19')
	GeneticAlgorithm::RunGeneticAlgorithm(true);

	// now that alg is finished, grab best solution and print the value of its genes (use the same translation as in the fitness eval for consistency)
	for (int i = 0; i < 4; i++) {
		if (i % 2 == 0)
			std::cout << "Gene " << i << " = " << *reinterpret_cast<double*>(&(GeneticAlgorithm::bestSolution.genes[i])) << "\n";
	}

	// clear dynamic memory used
	GeneticAlgorithm::ClearObject();
}
