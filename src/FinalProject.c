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

int cmpfunc(const void *ptrChromo1, const void *ptrChromo2){

    const struct chromosome* ptrChromosome1 = (struct chromosome*) ptrChromo1;
    const struct chromosome* ptrChromosome2 = (struct chromosome*) ptrChromo2;

    return (*ptrChromosome1).totalDistance - (*ptrChromosome2).totalDistance;
}

void displayChromosome(struct chromosome chromosome){

    printf("\n");
    printf("citiesAmount = %d\n", chromosome.citiesAmount);
    printf("citiesPermutation = ");

    for(int i = 0; i < chromosome.citiesAmount; i++)
        printf("%d ", chromosome.citiesPermutation[i]);

    printf("\n");

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

    double totalDistance = 0;
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
static int rand_int(int n)
{
  int limit = RAND_MAX - RAND_MAX % n;
  int rnd;

  do {
    rnd = rand();
  }
  while(rnd >= limit);
  return rnd %n;
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
 int permutationAlreadyExists(int *newPopulation,int ** previousPopulations, int numberOfCities, int filas,struct chromosome chromosomeArray[]){

   int i = 0;
   int j = 0;
   int checksum = 0;
   int checktrue = 50;
   //Si el arreglo esta vacio, lo llena con la primera permutacion

   if((*(*previousPopulations+0)+0) == NULL){

     for(i = 0; i < numberOfCities; i++){

       checktrue = newPopulation[i];
       previousPopulations[0][i] = checktrue;

     }
     chromosomeArray[0].citiesPermutation = previousPopulations[0];
     return filas;
   }

   //Reasigna el tamaño del arreglo y le agrega filas


   filas++;

   for(i = 0; i < numberOfCities; i++){

     for(j = 0; j < filas - 1; j++){

       if(previousPopulations[j][i] != newPopulation[i]){

         checksum++;
       }
     }
   }

   if(checksum != 0){

     for(i = 0; i < numberOfCities; i++){

       previousPopulations[filas - 1][i] = newPopulation[i];

     }
     chromosomeArray[filas-1].citiesPermutation = previousPopulations[filas-1];
   }

   return filas;
 }

/**
 * Function to generate a random permutation of a given size by creating an array
 * with the integers from 0 to size - 1 and shuffling that array.
 *
 * @param size the size of the desired permutation
 * @return a pointer to an integer representing the array of the shuffled numbers
 * from 0 to size - 1.
 */
 void generateRandomPermutation(int *array, int amountOfCities)
 {
   int i, j, tmp;
   for (i=amountOfCities-1;i>0;i--)
   {
     j = rand_int(i+1);
     tmp = array[j];
     array[j] = array[i];
     array[i] = tmp;
   }
 }

struct chromosome createChildChromosome(struct chromosome c1, struct chromosome c2){

    /*
        INSERT LOGIC ABOUT MERGING CHROMOSOMES HERE
    */

    return c1;
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
void createChromosome(struct chromosome* ptrChromosome, int amountOfCities, struct city citiesArray[], int generation){

  //  (*ptrChromosome).citiesPermutation = generateRandomPermutation(amountOfCities); //Set the population
    (*ptrChromosome).citiesAmount = amountOfCities; //Set the amout of cities
    (*ptrChromosome).generation = generation; //Set the generation
    setChromosomeTotalDistance(ptrChromosome, citiesArray, amountOfCities); //Set the total distance traveled
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

    time_t t;

    for(int i = 0; i < totalChromosomes; i++){

        /* initialize random seed. Make i a parameter. */
        srand((unsigned) time(&t) + i*i);
        createChromosome(&(chromosomesArray[i]), amountOfCities, citiesArray, 0);
        // displayChromosome(chromosomesArray[i]);
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
            Pick the best 25 chromosomes and then make 75 new chromosomes with the following pairing schemas:

                1. The i-th chromosome with the (i + 1)-th chromosome
                2. The i-th chromosome with the (i + 2)-th chromosome
                3. The i-th chromosome with the (i + 3)-th chromosome
        */

        for(int k = 1; k <= 3; k++){

            for(int j = 0; j < 25; j++){

                //Create child chromosomes for the indices between (bestToChromosomesToBeTaken) and (totalChromosomes - 1), inclusive
                chromosomesArray[k*bestToChromosomesToBeTaken + j] = createChildChromosome(chromosomesArray[j], chromosomesArray[(j + k)%totalChromosomes]);
            }
        }
    }

    //Sort the array once last time
    qsort(chromosomesArray, totalChromosomes, sizeof(struct chromosome), cmpfunc);

    //Return the first element of the array after being sorted, i.e. the chromosome with the least distance
    return chromosomesArray[0];
}

//main function
int main(){
    int filas = 1;
    int i;
    int numberCities = readCities();//store the number of cities in the file
    struct city citiesArray[numberCities];//array of structures of city
    struct chromosome chromosomeArray[10000];
    int *numbers; //Es la permutacion
    int **previousPopulations; //arreglo permutaciones
    previousPopulations = (int **)malloc(10000 * sizeof(int*));
    for(i=0;i<10000;i++)
      previousPopulations[i] = (int*)malloc(numberCities * sizeof(int));

    numbers = malloc(numberCities * sizeof(int));
    for(int i = 0; i < numberCities; i++)
      numbers[i] = i; //Fillas array with numbers in order
    for(i=0;i<10000;i++)
    {
      generateRandomPermutation(numbers,numberCities); //makes random array
      filas = permutationAlreadyExists(numbers,previousPopulations,numberCities,filas,chromosomeArray); //stores array
    }
    readInput(citiesArray, numberCities);//read the values of latitude and longitude in the file
    displayCities(citiesArray, numberCities);//Display the cities information.
    struct chromosome shortestPathChromosome = solve(numberCities, citiesArray);
    displayChromosome(shortestPathChromosome);

    return 0;
}//end of the main function
