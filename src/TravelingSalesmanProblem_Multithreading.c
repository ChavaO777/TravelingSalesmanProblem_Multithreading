// TravelingSalesmanProblem_Multithreading.c

struct individual{

	int size;
	int* permutation;
	double value;
}

struct coordinates{

	double latitude;
	double longitude;
}

int* generateRandomPopulation(int amountOfCities);
int populationAlreadyExists(int* newPermutation, int** previousPermutations);
void saveValidPopulation(int* newPermutation, int** previousPermutations);

void readInput();
double computeDistance();

int main(){

	readInput();

	return 0;
}