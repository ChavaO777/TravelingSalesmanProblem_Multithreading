#define _GNU_SOURCE
#include <pthread.h>
#include <syscall.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <assert.h>
#include <math.h>
#include <time.h>
int coreAssign(int core_id); //declararion of coreAssign
#define EARTH_RADIUS_KM 6371 // authalic radius based on/extracted from surface area;
/**
 * Structure of a argumenst. these structure will use to pass parameters to the threads 
 */
struct arguments{

    int numberCities;
    struct city* citiesList;
    int core;
};//end of the arguments struct
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

/**
 * Function to compute the score of a chromosome.
 * 
 * @param chromosome the chromosome whose score is to be computed
 * @return the score of the chromosome
 */ 


double computeChromosomeScore(struct chromosome chromosome){

    // The score of the chromosome is only its total distance; the lower, the better
    return chromosome.totalDistance;
}

/**
 * Function to compare chromosomes according to their distance and generation
 * 
 */ 
int cmpfunc(const void *ptrChromo1, const void *ptrChromo2){

    //Get the first chromosome
    const struct chromosome* ptrChromosome1 = (struct chromosome*) ptrChromo1;
    //Get the second chromosome
    const struct chromosome* ptrChromosome2 = (struct chromosome*) ptrChromo2;
    // Compute the subtraction of their scores
    return computeChromosomeScore(*ptrChromosome1) - computeChromosomeScore(*ptrChromosome2);
}

/**
 *  Function to display the information of a given chromosome
 * 
 *  @param chromosome the chromosome whose information is to be displayed
 */
void displayChromosome(struct chromosome chromosome){
    
    printf("\n");
    printf("citiesAmount = %d\n", chromosome.citiesAmount);
    printf("citiesPermutation = ");

    //In case the permutation is NULL
    if(chromosome.citiesPermutation == NULL){

        printf("The permutation is NULL.\n");
    }

    else{

        //Print the whole permutation of the cities
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

        //Set the id of the city
        citiesArray[i].id = i; 
        //Read the city's latitude and longitude
        scanf("%lf %lf", &citiesArray[i].latitude, &citiesArray[i].longitude);
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

//Function to swap to values in an array
void swap(int* array, int index1, int index2){

    int tmp = array[index2]; //Initialize a temporary variable
    array[index2] = array[index1];
    array[index1] = tmp;
}

/**Implementation of the Fisher-Yates algorithm to shuffle an array
 * 
 * @param array the array that is to be shuffled
 * @param size the size of the array
 * @param seedParameter the parameter for initializing the random seed
 */ 
void shuffle(int *array, int size, int seedParameter){

    time_t t;
    /* initialize random seed. Make i a parameter. */
    srand((unsigned) time(&t) + seedParameter*seedParameter);

    int i, j, tmp;
    for(i = size - 1; i > 0; i--){

        j = rand()%(i + 1); //Pick a random number between 0 and i
        swap(array, i, j); //Swap array[i] and array[j]
    }
}

/**
 *  Function to display an array
 *  
 *  @param arr the array to be displayed
 *  @size the size of the array
 */ 
void displayArray(int* arr, int size){

    //Check if the array is NULL
    if(arr == NULL){

        printf("Null array!\n");
        return;
    }

    //Print the whole array
    for(int i = 0; i < size; i++)
        printf("%d ", *(arr + i));

    printf("\n");
}

/**
 *  Function to determine if to permutations are equal
 *  
 *  @param permutation1 the first permutation to be compared
 *  @param permutation2 the second permutation to be compared
 *  @param the size of the permutations
 * 
 *  @return 1 if the permutations are equal. 0, otherwise.
 */ 
int areEqualPermutations(int *permutation1, int* permutation2, int size){   

    if(permutation1 == NULL || permutation2 == NULL)
        return 0;

    for(int i = 0; i < size; i++){ //Iterate through the whole permutations

        if(permutation1[i] != permutation2[i]) //If there is a mismatch, then they're not equal
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
 * @param totalChromosomes the amount of chromosomes that are being stored
 * @param chromosomesArray the array of chromosomes we currently have
 *
 * @return 1 if the permutation already exists. Else, 0.
 */
 int permutationAlreadyExists(int* citiesPermutation, int size, int totalChromosomes, struct chromosome chromosomesArray[]){

    for(int i = 0; i < totalChromosomes; i++){

        //If the current permutation is NULL, just skip it
        if(chromosomesArray[i].citiesPermutation == NULL){

            continue;
        }

        //If they're equal, then the permutation already exists
        if(areEqualPermutations(citiesPermutation, chromosomesArray[i].citiesPermutation, size))
            return 1;
    }
    
    //In case the permutation didn't exist yet
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

/**
 *  Function to determine if a given value exists in an array
 *  
 *  @param size the size of the array
 *  @param value the value we're looking for
 *  @param array the array in which we're going to look for the value
 * 
 *  @return 1 in case the value exists. Else, 0.
 */ 
int valueExistsInArray(int size, int value, int* array){

    for(int i = 0; i < size; i++)
        if(array[i] == value)
            return 1;

    return 0;
}

/**
 * Function to determine if a given value is in a given range
 * 
 * @param value the value to be checked
 * @param loValue the minimum value of the range, inclusive
 * @param hiValue the maximum value of the range, inclusive
 * 
 * @return 1 in case the value is in range. Else, 0.
 */ 
int isInRange(int value, int loValue, int hiValue){

    return loValue <= value && value <= hiValue;
}

/**
 * Function to pass a subset of the values of a given array directly
 * 
 * @param size the size of the array
 * @param newArr the new array in which the values are going to be assigned
 * @param leftIndex the first index for the copy procedure
 * @param rightIndex the last index for the copy procedure
 * @param parentArr the array from which the values are being passed to the new array
 */ 
void passSomeValuesDirectly(int size, int* newArr, int leftIndex, int rightIndex, int* parentArr){

    // Initialize the new array with arr1's values in the picked range and -1's elswhere
    for(int i = 0; i < size; i++){

        //If the index is in range, just pass this value form arr1 to the new array
        if(isInRange(i, leftIndex, rightIndex))
            newArr[i] = parentArr[i];

        else //Initialize with a fake negative value
            newArr[i] = -1; 
    }
}

/**
 * Function to pass the remaining required values from a parent array to a new array
 * 
 * @param size the size of the array
 * @param newArr the new array in which the values are going to be assigned
 * @param leftIndex the first index of the range that has to be avoided
 * @param rightIndex the last index of the range that has to be avoided
 * @param parentArr the array from which the values are being passed to the new array
 */ 
void passTheRestOfTheValues(int size, int* newArr, int leftIndex, int rightIndex, int* parentArr){

    //Index of parentArr
    int parentArrIndex = 0;
    for(int i = 0; i < size; i++){

        //If the value was already set, just skip this index
        if(isInRange(i, leftIndex, rightIndex))
            continue;

        //While the current value of parentArr exists in the new array, move to the right 
        while(parentArrIndex < size && valueExistsInArray(size, parentArr[parentArrIndex], newArr))
            parentArrIndex++;

        //Assign this value to the i-th index of the new array
        newArr[i] = parentArr[parentArrIndex];
    }
}

/**
 * Function to pick a random index range. Both indices could be equal,
 * thus generating a range of size 1
 * 
 * @param ptrInd1 a pointer to the first index
 * @param ptrInd2 a pointer to the second index
 * @param size the maximum size of the range
 * @param seedParameter the parameter for initializing the random seed
 */ 
void pickRandomIndexRange(int* ptrInd1, int* ptrInd2, int size, int seedParameter){

    /* initialize random seed.*/
    srand(time(NULL));
    seedParameter += rand()%1000;

    time_t t;
    /* initialize random seed. Use the seedParameter. */
    srand((unsigned) time(&t) + seedParameter*seedParameter);

    //Pick an index in arr1
    *ptrInd1 = rand()%(size/2);
    //Pick another index in arr1 which is to the right of the first picked index
    *ptrInd2 = (*ptrInd1) + rand()%(size/2);
}

/**
 * Function to pick two different indices in a given range
 * 
 * @param ptrInd1 a pointer to the first index
 * @param ptrInd2 a pointer to the second index
 * @param size the maximum size of the range
 * @param seedParameter the parameter for initializing the random seed
 */ 
void pickRandomDifferentIndicesInRange(int* ptrInd1, int* ptrInd2, int size, int seedParameter){

    /* initialize random seed.*/
    srand(time(NULL));
    //Add to the seedParameter a random value
    seedParameter += rand()%1000;

    time_t t;
    /* initialize random seed. Use the parameter. */
    srand((unsigned) time(&t) + seedParameter*seedParameter);

    do{
        //Pick the first index
        *ptrInd1 = rand()%size;
        //Pick the second index
        *ptrInd2 = rand()%size;
    }while(*ptrInd1 == *ptrInd2); //Cycle while they're equal
}

// Idea taken from: http://www.theprojectspot.com/tutorial-post/applying-a-genetic-algorithm-to-the-travelling-salesman-problem/5
/**
 * Function to perform an ordered crossover between two arrays
 * 
 * @param size the size of the arrays
 * @param arr1 the first array 
 * @param arr2 the second array
 * @param seedParameter the parameter for initializing the random seed
 * 
 * @return a new array which is a child of arr1 and arr2 after the crossover has been made
 */ 
int* performOrderedCrossover(int size, int* arr1, int* arr2, int seedParameter){

    //Declare variables for the indices
    int arr1SubsetLeftIndex, arr1SubsetRightIndex;
    //Assign values to those indices
    pickRandomIndexRange(&arr1SubsetLeftIndex, &arr1SubsetRightIndex, size, seedParameter);

    //Create a new array of size "size".
    int* newArr = malloc(sizeof(int)*size);

    //Pass some values directly from arr1
    passSomeValuesDirectly(size, newArr, arr1SubsetLeftIndex, arr1SubsetRightIndex, arr1);
    //Pass the remaining values from arr2
    passTheRestOfTheValues(size, newArr, arr1SubsetLeftIndex, arr1SubsetRightIndex, arr2);

    return newArr;
}

/**
 *  Function to perform a mutation in an array
 * 
 *  @param size the size of the array
 *  @param arr the array in which a mutation is to be performed
 *  @param seedParameter the parameter for initializing the random seed
 */
void performMutation(int size, int* arr, int seedParameter){

    //Declare to indices
    int index1 = -1;
    int index2 = -1;
    //Assign values to those indices
    pickRandomDifferentIndicesInRange(&index1, &index2, size, seedParameter);
    //Swap the values in those indices
    swap(arr, index1, index2);
}

/**
 * Function to create a child permutation from two given parent chromosomes 
 * 
 * @param c1 the first parent chromosome
 * @param c2 the second parent chromosome
 * 
 * @return a child permutation created by crossing over the parent's permutations
 */
int* createPermutationFromParents(struct chromosome c1, struct chromosome c2){

    //An array of primes to add randomness to the crossover
    int primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};

    //Create the cross over
    int* childPermutation = performOrderedCrossover(c1.citiesAmount, c1.citiesPermutation, c2.citiesPermutation, c1.generation*c2.generation + primes[rand()%15]);
    
    //Total mutations to be performed
    const int totalMutations = 1 + rand()%3;

    //Perform mutations in the array returned by the crossover
    for(int i = 0; i < totalMutations; i++){

        performMutation(c1.citiesAmount, childPermutation, c1.generation*c2.generation + primes[rand()%15]);
    }
    
    return childPermutation;
}

/**
 * Function to generate an array filled with -1's
 *  
 * @param size the size of the desired array
 * @return an array of size "size" filled with -1's
 */
int* createArrayWithNegativeValues(int size){

    //Create the array
    int* arr = malloc(sizeof(int)*size);

    //Initialize with -1's
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

/**
 * Function to create a child chromosome form two parent chromosomes
 * 
 * @param c1 the first parent chromosome
 * @param c2 the second parent chromosome
 * @param generation the generation of the child chromosome
 * @param amountOfCities the amount of cities for the chromosomes
 * @param citiesArray the array of cities for calculating the total distance of a chromosome
 * @return the child chromosome
 */ 
struct chromosome createChildChromosome(struct chromosome c1, struct chromosome c2, int generation, int amountOfCities, struct city citiesArray[]){

    // Declare a new chromosome
    struct chromosome childChromosome;
    // Create a child permutation for the child chromosome
    childChromosome.citiesPermutation = createPermutationFromParents(c1, c2);
    // Assign the amount of cities
    childChromosome.citiesAmount = c1.citiesAmount;
    // Assign the generation of the child chromosome
    childChromosome.generation = generation;
    setChromosomeTotalDistance(&childChromosome, citiesArray, amountOfCities); //Set the total distance traveled

    return childChromosome;  
}

/**
 * Function to display an array of chromosomes
 * @param size the size of the array
 * @param chromosomesArray the array of chromosomes
 */
void displayChromosomesArray(int size, struct chromosome chromosomesArray[]){

    for(int i = 0; i < size; i++){

        printf("Chromosome #%d:", i);
        //Display the information of this chromosome
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
   
    // Amount of chromosomes we will be working with at any given generation
    const int totalChromosomes = 100;

    //Array of chromosomes we will be working with
    struct chromosome chromosomesArray[totalChromosomes];
 
    for(int i = 0; i < totalChromosomes; i++){ 
        
        //Initialize all the permutations with fake values
        chromosomesArray[i].citiesPermutation = createArrayWithNegativeValues(amountOfCities);
    }
    
    for(int i = 0; i < totalChromosomes; i++){

        //Create a random chromosome
        createRandomChromosome(&(chromosomesArray[i]), i + 1, amountOfCities, citiesArray, 0, totalChromosomes, chromosomesArray);
    }

    // Amount of chromosome generations we will be working with
    const int totalGenerations = 1000;
    // Amount of "elite" chromosomes we will be taking; the elite ones are the ones with the least total distance
    const int bestToChromosomesToBeTaken = totalChromosomes/4;

    // Iterate through the generations
    for(int i = 0; i < totalGenerations; i++){

        //Sort the array of chromosomes
        qsort(chromosomesArray, totalChromosomes, sizeof(struct chromosome), cmpfunc);
        
        /*  
            Pick the best bestToChromosomesToBeTaken chromosomes and then make (totalChromosomes - bestToChromosomesToBeTaken) new chromosomes with the following pairing schemas:

                1. The i-th chromosome with the (i + 1)-th chromosome
                2. The i-th chromosome with the (i + 2)-th chromosome
                3. The i-th chromosome with the (i + 3)-th chromosome
        */

        // Create child chromosomes for the indices between (bestToChromosomesToBeTaken) and (totalChromosomes - 1), inclusive
        for(int k = 1; k <= 3; k++){

            for(int j = 0; j < 25; j++){

                //Declare the index of the first parent
                int parent1Index = j;
                //Declare the index of the second parent
                int parent2Index = (j + k)%bestToChromosomesToBeTaken;
                //Declare the index of the child chromosome
                int childChromosomeIndex = k*bestToChromosomesToBeTaken + j;
                //Create the child chromosome
                chromosomesArray[childChromosomeIndex] = createChildChromosome(chromosomesArray[parent1Index], chromosomesArray[parent2Index], i, amountOfCities, citiesArray);
            }
        }
    }

    //Sort the array one last time
    qsort(chromosomesArray, totalChromosomes, sizeof(struct chromosome), cmpfunc);

    //Display the final chromosome array
    //displayChromosomesArray(totalChromosomes, chromosomesArray);

    //Return the first element of the array after being sorted, i.e. the chromosome with the least distance
    return chromosomesArray[0];
}

 void* threadSolution(void *arg)
{
    struct arguments* id = (struct arguments*) arg;
    int error = coreAssign((*id).core);//call to coreAssign function and check if are aviable core
    if(err == 0)
    {
        struct chromosome shortestPathChromosome = solve( (*id).numberCities, (*id).citiesList);
        int sid = syscall(SYS_gettid);//get the number of the current thread
        printf("In thread %d  Winning chromosome:",sid);
        displayChromosome(shortestPathChromosome );
        printf("\n");
    }
    else{
        int sid = syscall(SYS_gettid);//get the number of the current thread
        printf("No core Avaibale for thread %d \n",sid);
    }
}

int coreAssign(int core_id) {
   int numberOfCores = sysconf(_SC_NPROCESSORS_ONLN);//get number of processors which are currently online 
   if (core_id < 0 || core_id >= numberOfCores)
      return 1;//retunr 1 if its an error in the core_id
   cpu_set_t processor;//declare a set of CPUs
   CPU_ZERO(&processor);// Clears set, so that it contains no CPUs.
   CPU_SET(core_id, &processor);// Add CPU cpu to set.
   pthread_t actualThread = pthread_self(); //call the actual thread id
   return pthread_setaffinity_np(actualThread, sizeof(cpu_set_t), &processor);//sets the thread in the cpu assigned
}
//main function
int main(){
    int numbersProcessors = 4;//set numer of processor
    int numberCities = readCities();//store the number of cities in the file
    struct city citiesArray[numberCities];//array of structures of city
    readInput(citiesArray, numberCities);//read the values of latitude and longitude in the file
    displayCities(citiesArray, numberCities);//Display the cities information.
    struct arguments args; //create a new srtucture for the threads
    args.numberCities = numberCities;
    args.citiesList = citiesArray;   
    pthread_t tid[numbersProcessors];//declaration for threads
    for(int threadFor = 0; threadFor < numbersProcessors;threadFor++ ){
        args.core =threadFor;//set the procesor number
        pthread_create(&(tid[threadFor]), NULL, threadSolution, &args);
   }   
    for(int threadFor = 0; threadFor < numbersProcessors;threadFor++ ){
         pthread_join(tid[threadFor], NULL);
   }   
  
    return 0;
}//end of the main function
