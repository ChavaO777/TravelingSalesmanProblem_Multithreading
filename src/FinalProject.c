#include <stdio.h>
#include <math.h>

int totalCities = 0; // total amount of city coordinates in the input file 
double latitude[100]; // Array for storing the latitude of the cities 
double longitude[100]; // Array for storing the longitude  of the cities

/** Structure of an individual. It includes an array of integers representing
*   the permutation of the n cities (0 to n - 1) that form that solution and 
*   the total distance traveled to visit all n cities in that order and come
*   back to the first city.
*/ 
struct individual {
    
    int* citiesPermutation; // a pointer to an integer corresponding to the permutation of the n cities
    double totalDistance; // the total distance traveled
};//end of the individual struct

//Function that reads the input
void readInput(){

    scanf("%d", &totalCities);

    for(int i = 0; i < totalCities; i++){

        scanf("%lf %lf", &latitude[i], &longitude[i]);
    }
}//end of the readInput function 

/**
*   Computes the distance between two cities given their indices. It is an implementation
*   of the Haversine formula (https://en.wikipedia.org/wiki/Haversine_formula)
*/
double computeDistance(int indexCoord1, int indexCoord2){
    
   
}//end of the computeDistance function

struct individual* solution computeOverallDistance(int* population) {


}   

//main function
int main(){

    readInput();
    return 0;
}//end of the main function