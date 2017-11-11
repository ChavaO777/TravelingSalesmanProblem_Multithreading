#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define EARTH_RADIUS_KM 6371 // authalic radius based on/extracted from surface area;

/**
 * Structure of a city. It includes an integer id, its latitude and longitude
 */
struct city{

    int id; //id of the city
    double latitude; //latitude of the cities
    double longitude; //longitude  of the cities

};//end of the city struct

/** Structure of a chromosome. It includes an array of integers representing
*   the permutation of the n cities (0 to n - 1) that form that solution and
*   the total distance traveled to visit all n cities in that order and come
*   back to the first city.
*/
struct chromosome {

    int citiesAmount;
    int* citiesPermutation; // a pointer to an integer corresponding to the permutation of the n cities
    double totalDistance; // the total distance traveled
    int generation;
};//end of the chromosome struct

double computeChromosomeScore(struct chromosome chromosome){

    return chromosome.totalDistance;
}

/**
 * Function to compare chromosomes according to their distance and generation
 * 
 */ 
int cmpfunc(const void *ptrChromo1, const void *ptrChromo2){

    const struct chromosome* ptrChromosome1 = (struct chromosome*) ptrChromo1;
    const struct chromosome* ptrChromosome2 = (struct chromosome*) ptrChromo2;
    return computeChromosomeScore(*ptrChromosome1) - computeChromosomeScore(*ptrChromosome2);
}

void displayChromosome(struct chromosome chromosome){
    
    printf("\n");
    printf("citiesAmount = %d\n", chromosome.citiesAmount);
    printf("citiesPermutation = ");

    if(chromosome.citiesPermutation == NULL){

        printf("The permutation is NULL.\n");
    }

    else{

        for(int i = 0; i < chromosome.citiesAmount; i++)
            printf("%d ", chromosome.citiesPermutation[i]);

        printf("\n");
    }

    printf("totalDistance = %lf\n", chromosome.totalDistance);
    printf("generation = %d\n", chromosome.generation);
    printf("\n");
}

//Function that reads the input
void readInput(struct city citiesArray[], int numberCities){

    for(int i = 0; i < numberCities; i++){

        citiesArray[i].id = i; //Set the id of the city
        scanf("%lf %lf", &citiesArray[i].latitude, &citiesArray[i].longitude); //Read the city's latitude and longitude
    }
}//end of the readInput function

//function that reads the number of cities
int readCities(){

    int totalCities;
    scanf("%d", &totalCities);
    return totalCities;
}

//Function to display the read input
void displayCities(struct city citiesArray[], int numberCities){

    printf("\n");
    printf("The input is the following:\n\n"); //Break line for aesthetics

    for(int i = 0; i < numberCities; i++){

        //Print all the information of a given city
        printf("City id: %d; Latitude: %lf; Longitude: %f\n", citiesArray[i].id, citiesArray[i].latitude, citiesArray[i].longitude);
    }

    printf("\n"); //Break line for aesthetics
}//End of the displayCities function

// This function converts decimal degrees to radians
double convertDegreesToRadians(double degrees) {

    return (degrees * M_PI / 180);
}

// This function converts radians to decimal degrees
double convertRadiansToDegrees(double radians) {

    return (radians * 180.0 / M_PI);
}

/**
*   Computes the distance between two cities given their indices. It is an implementation
*   of the Haversine formula (https://en.wikipedia.org/wiki/Haversine_formula)
*/
double computeDistanceBetweenCities(int indexCoord1, int indexCoord2, struct city citiesArray[]){

    double lat1deg, lng1deg, lat2deg, lng2deg;

    lat1deg = citiesArray[indexCoord1].latitude; //Read the latitude of the first city in degrees
    lng1deg = citiesArray[indexCoord1].longitude; //Read the longitude of the first city in degrees
    lat2deg = citiesArray[indexCoord2].latitude; //Read the latitude of the second city in degrees
    lng2deg = citiesArray[indexCoord2].longitude; //Read the longitude of the second city in degrees

    double lat1rad, lng1rad, lat2rad, lng2rad, u, v;
    lat1rad = convertDegreesToRadians(lat1deg); //Convert the latitude of the first city to degrees
    lng1rad = convertDegreesToRadians(lng1deg); //Convert the longitude of the first city to degrees
    lat2rad = convertDegreesToRadians(lat2deg); //Convert the latitude of the second city to degrees
    lng2rad = convertDegreesToRadians(lng2deg); //Convert the longitude of the second city to degrees
    u = sin((lat2rad - lat1rad)/2.0); //Compute the sine of the latitude differential
    v = sin((lng2rad - lng1rad)/2.0); //Compute the sine of the longitude diferential

    //Haversine formula
    return 2.0 * EARTH_RADIUS_KM * asin(sqrt(u * u + cos(lat1rad) * cos(lat2rad) * v * v));
}//end of the computeDistance function

/**
 * Function that sets the total distance traveled using the cities permutation of a given chromosome.
 *
 * @param ptrChromosome a pointer to a chromosome struct whose total distance traveled we want to calculte
 * @param citiesArray the array of cities given by the input
 * @param numberCities an integer number representing the total amount of cities in the input
 */
void setChromosomeTotalDistance(struct chromosome* ptrChromosome, struct city citiesArray[], int numberCities){

    double totalDistance = 0.0;
    int* citiesPermutation = (*ptrChromosome).citiesPermutation;//Pointer where chromosome are

    for(int citiesIndex = 1; citiesIndex < numberCities; citiesIndex++){

        //call the computeDistanceBetweenCities function for these two cities and add to the totalDistance variable
        totalDistance += computeDistanceBetweenCities(citiesPermutation[citiesIndex - 1], citiesPermutation[citiesIndex], citiesArray);
    }

    //The last trip is from the last city to the first one
    totalDistance += computeDistanceBetweenCities(citiesPermutation[numberCities - 1], citiesPermutation[0], citiesArray);
    (*ptrChromosome).totalDistance = totalDistance;//store the total distance in the chromosome structure
} //End of the setChromosomeTotalDistance() function

/**
*   Function to calculate the total distance traveled by all chromosomes in the array.
*   The function receives, the structure array where the chromosomes are, the
*   structure array where the cities are and the total number of chromosomes.
*/
void setAllChromosomesTotalDistance(struct chromosome chromosomeArray[], struct city citiesArray[], int numberCities, int numberOfChromosomes) {

    for(int chromosomeIndex = 0; chromosomeIndex < numberOfChromosomes; chromosomeIndex++){

        setChromosomeTotalDistance(&(chromosomeArray[chromosomeIndex]), citiesArray, numberCities);
    }
}//end of computeOverallDistance

//Function to swap to values in an array
void swap(int* array, int index1, int index2){

    int tmp = array[index2]; //Initialize a temporary variable
    array[index2] = array[index1];
    array[index1] = tmp;
}

//Implementation of the Fisher-Yates algorithm
void shuffle(int *array, int amountOfCities, int seedParameter){

    time_t t;
    /* initialize random seed. Make i a parameter. */
    srand((unsigned) time(&t) + seedParameter*seedParameter);
    int i, j, tmp;
    for(i = amountOfCities - 1; i > 0; i--){

        j = rand()%(i + 1); //Pick a random number between 0 and i
        swap(array, i, j); //Swap array[i] and array[j]
    }
}

void displayArray(int* arr, int size){

    if(arr == NULL){

        printf("Null array!\n");
        return;
    }

    for(int i = 0; i < size; i++)
        printf("%d ", *(arr + i));

    printf("\n");
}

int areEqualPermutations(int *permutation1, int* permutation2, int size){   

    if(permutation1 == NULL || permutation2 == NULL)
        return 0;

    for(int i = 0; i < size; i++){ //Iterate through the whole permutations

        if(permutation1[i] != permutation2[i]) //If there is a mismatch, they're not equal{
            return 0;
    }

    return 1; //If there was no mismatch, then they're equal
}

/**
 * Function that returns a 1 if a permutation had already been generated before or a
 * zero if it is a new permutation
 *
 * @param citiesPermutation the permutation to be checked
 * @param size the size of the permutation from 0 to (size - 1)
 *
 * @return 1 if the permutation already exists. Else, 0.
 */
 int permutationAlreadyExists(int* citiesPermutation, int size, int totalChromosomes, struct chromosome chromosomesArray[]){

    for(int i = 0; i < totalChromosomes; i++){

        if(chromosomesArray[i].citiesPermutation == NULL){

            continue;
        }

        if(areEqualPermutations(citiesPermutation, chromosomesArray[i].citiesPermutation, size))
            return 1;
    }
   
   return 0;
 }

/**
 * Function to generate a random permutation of a given size by creating an array
 * with the integers from 0 to size - 1 and shuffling that array.
 * 
 * @param size the size of the desired permutation
 * @return a pointer to an integer representing the array of the shuffled numbers
 * from 0 to size - 1.
 */ 
int* generateRandomPermutation(int size, int seedParameter, int totalChromosomes, struct chromosome chromosomesArray[]){
    
    int* citiesPermutation = malloc(sizeof(int)*size);
    
    for(int i = 0; i < size; i++)
        citiesPermutation[i] = i; //Assign the values from 0 to n - 1

    int counter = 0;

    do{
        shuffle(citiesPermutation, size, seedParameter); //Shuffle the array
    }while(permutationAlreadyExists(citiesPermutation, size, totalChromosomes, chromosomesArray));
    
    return citiesPermutation; //Return the array
}

int valueExistsInArray(int size, int value, int* array){

    for(int i = 0; i < size; i++)
        if(array[i] == value)
            return 1;

    return 0;
}

int indexWasFilledBefore(int index, int leftSubsetIndex, int rightSubsetIndex){

    return leftSubsetIndex <= index && index <= rightSubsetIndex;
}

// Idea taken from: http://www.theprojectspot.com/tutorial-post/applying-a-genetic-algorithm-to-the-travelling-salesman-problem/5
int* orderedCrossover(int size, int* arr1, int* arr2, int seedParameter){

    srand(time(NULL));
    seedParameter += rand()%1000;

    time_t t;
    /* initialize random seed. Make i a parameter. */
    srand((unsigned) time(&t) + seedParameter*seedParameter);

    int arr1SubsetLeftIndex = rand()%(size/2);
    int arr1SubsetRightIndex = arr1SubsetLeftIndex + rand()%(size/2);

    int* newArr = malloc(sizeof(int)*size);

    for(int i = 0; i < size; i++){

        //Just pass this value to the new array
        if(arr1SubsetLeftIndex <= i && i <= arr1SubsetRightIndex)
            newArr[i] = arr1[i];

        else //Initialize with a fake negative value
            newArr[i] = -1; 
    }

    int arr2Index = 0;
    for(int i = 0; i < size; i++){

        if(indexWasFilledBefore(i, arr1SubsetLeftIndex, arr1SubsetRightIndex))
            continue;

        while(arr2Index < size && valueExistsInArray(size, arr2[arr2Index], newArr))
            arr2Index++;

        newArr[i] = arr2[arr2Index];
    }

    return newArr;
}

void performMutation(int size, int* arr, int seedParameter){

    srand(time(NULL));
    seedParameter += rand()%1000;

    time_t t;
    /* initialize random seed. Make i a parameter. */
    srand((unsigned) time(&t) + seedParameter*seedParameter);

    int index1 = -1;
    int index2 = -1;

    do{
        index1 = rand()%size;
        index2 = rand()%size;
    }while(index1 == index2);

    swap(arr, index1, index2);
}

int* createPermutationFromParents(struct chromosome c1, struct chromosome c2){

    /*
        INSERT LOGIC ABOUT MERGING CHROMOSOMES HERE
    */

    int primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};

    int* childPermutation = orderedCrossover(c1.citiesAmount, c1.citiesPermutation, c2.citiesPermutation, c1.generation*c2.generation + primes[rand()%15]);
    performMutation(c1.citiesAmount, childPermutation, c1.generation*c2.generation + primes[rand()%15]);
    
    return childPermutation;
}

int* setToNegativesArray(int size){

    int* arr = malloc(sizeof(int)*size);

    for(int i = 0; i < size; i++){

        arr[i] = -1;
    }

    return arr;
}

/**
 * Function to create a chromosome given the amount of cities its permutation requires, the cities array
 * to calculate the total distance traveled and the generation of the chromosome.
 *
 * @param ptrChromosome a pointer to the chrosome whose attributes are to be initialized
 * @param amountOfCities the amount of cities that the chromosome's permutation requires
 * @param citiesArray the array of cities in which the coordinates and the id of each city
 * are stored.
 * @generation the generation of the chromosome since the start of the program's execution
 */ 
void createRandomChromosome(struct chromosome* ptrChromosome, int seedParameter, int amountOfCities, struct city citiesArray[], int generation, int totalChromosomes, struct chromosome chromosomesArray[]){

    (*ptrChromosome).citiesPermutation = generateRandomPermutation(amountOfCities, seedParameter, totalChromosomes, chromosomesArray); //Set the population
    (*ptrChromosome).citiesAmount = amountOfCities; //Set the amout of cities
    (*ptrChromosome).generation = generation; //Set the generation
    setChromosomeTotalDistance(ptrChromosome, citiesArray, amountOfCities); //Set the total distance traveled
}

struct chromosome createChildChromosome(struct chromosome c1, struct chromosome c2, int generation, int amountOfCities, struct city citiesArray[]){
    
    struct chromosome childChromosome;
    childChromosome.citiesPermutation = createPermutationFromParents(c1, c2);
    childChromosome.citiesAmount = c1.citiesAmount;
    childChromosome.generation = generation;
    setChromosomeTotalDistance(&childChromosome, citiesArray, amountOfCities); //Set the total distance traveled

    return childChromosome;  
}

void displayChromosomesArray(int totalChromosomes, struct chromosome chromosomesArray[]){

    for(int i = 0; i < totalChromosomes; i++){

        printf("Chromosome #%d:", i);
        displayChromosome(chromosomesArray[i]);
    }

    printf("\n");
}

/**
 * Function in which a genetic algorithm will be implemented to solve the TSP.
 *
 * @param amountOfCities the total amount of cities the traveler must travel
 * @param citiesArray the array of cities in which the coordinates and the id of each city
 * are stored.
 *
 * @return bestChromosome the chromosome with the best route for the traveler
 */
struct chromosome solve(int amountOfCities, struct city citiesArray[]){

    int totalChromosomes = 100;
    struct chromosome chromosomesArray[totalChromosomes];

    for(int i = 0; i < totalChromosomes; i++) //Initialize all the permutations with fake values
        chromosomesArray[i].citiesPermutation = setToNegativesArray(amountOfCities);

    for(int i = 0; i < totalChromosomes; i++){

        createRandomChromosome(&(chromosomesArray[i]), i + 1, amountOfCities, citiesArray, 0, totalChromosomes, chromosomesArray);
    }

    /*
        INSERT LOGIC OF GA HERE
    */

    int totalGenerations = 100;
    int bestToChromosomesToBeTaken = totalChromosomes/4;

    for(int i = 0; i < totalGenerations; i++){

        //Sort the array of chromosomes
        qsort(chromosomesArray, totalChromosomes, sizeof(struct chromosome), cmpfunc);
        
        /*  
            Pick the best bestToChromosomesToBeTaken chromosomes and then make (totalChromosomes - bestToChromosomesToBeTaken) new chromosomes with the following pairing schemas:

                1. The i-th chromosome with the (i + 1)-th chromosome
                2. The i-th chromosome with the (i + 2)-th chromosome
                3. The i-th chromosome with the (i + 3)-th chromosome
        */

        for(int k = 1; k <= 3; k++){

            for(int j = 0; j < 25; j++){

                // Create child chromosomes for the indices between (bestToChromosomesToBeTaken) and (totalChromosomes - 1), inclusive
                printf("chromosomesArray[%d] -> ", j);
                displayArray(chromosomesArray[j].citiesPermutation, chromosomesArray[j].citiesAmount);

                printf("chromosomesArray[%d] -> ", (j + k)%bestToChromosomesToBeTaken);
                displayArray(chromosomesArray[(j + k)%totalChromosomes].citiesPermutation, chromosomesArray[(j + k)%bestToChromosomesToBeTaken].citiesAmount);
                chromosomesArray[k*bestToChromosomesToBeTaken + j] = createChildChromosome(chromosomesArray[j], chromosomesArray[(j + k)%bestToChromosomesToBeTaken], i, amountOfCities, citiesArray);
            }
        }
    }

    //Sort the array one last time
    qsort(chromosomesArray, totalChromosomes, sizeof(struct chromosome), cmpfunc);

    displayChromosomesArray(totalChromosomes, chromosomesArray);

    //Return the first element of the array after being sorted, i.e. the chromosome with the least distance
    return chromosomesArray[0];
}

//main function
int main(){

    int numberCities = readCities();//store the number of cities in the file
    struct city citiesArray[numberCities];//array of structures of city
    readInput(citiesArray, numberCities);//read the values of latitude and longitude in the file
    displayCities(citiesArray, numberCities);//Display the cities information.
    struct chromosome shortestPathChromosome = solve(numberCities, citiesArray);
    
    printf("Winning chromosome:");
    displayChromosome(shortestPathChromosome);

    return 0;
}//end of the main function
