#include <stdio.h>
#include <math.h>
#include <stdlib.h>

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
    
    int* citiesPermutation; // a pointer to an integer corresponding to the permutation of the n cities
    double totalDistance; // the total distance traveled
};//end of the chromosome struct

//Function that reads the input
void readInput(struct city* citiesArray[], int numberCities){
  
    for(int i = 0; i < numberCities; i++){

        citiesArray[i] = (struct city*) malloc(sizeof(struct city));
        citiesArray[i]->id = i;
        scanf("%lf %lf", &citiesArray[i]->latitude, &citiesArray[i]->longitude);
    }
}//end of the readInput function 

//function that reads the number of cities
int readCities(){
    
    int totalCities;
    scanf("%d", &totalCities);
    return totalCities;
}

//Function to display the read input
void displayCities(struct city* citiesArray[], int numberCities){

    printf("\n"); //Break line for aesthetics

    for(int i = 0; i < numberCities; i++){
        
        //Print all the information of a given city
        printf("City id: %d; Latitude: %lf; Longitude: %f\n", citiesArray[i]->id, citiesArray[i]->latitude, citiesArray[i]->longitude);
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
double computeDistance(int indexCoord1, int indexCoord2, struct city* citiesArray[]){

    double lat1deg, lng1deg, lat2deg, lng2deg; 
    
    lat1deg = citiesArray[indexCoord1]->latitude; //Read the latitude of the first city in degrees
    lng1deg = citiesArray[indexCoord1]->longitude; //Read the longitude of the first city in degrees    
    lat2deg = citiesArray[indexCoord2]->latitude; //Read the latitude of the second city in degrees
    lng2deg = citiesArray[indexCoord2]->longitude; //Read the longitude of the second city in degrees

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
*   Function to calculate the total distance traveled by all chromosomes in the array.
*   The function receives, the structure array where the chromosomes are, the 
*   structure array where the cities are and the total number of chromosomes.
*/
void computeOverallDistance(struct chromosome* chromosomeArray[], struct city* citiesArray[], int numberCities, int numberOfChromosomes) {
    
    for(int chromosomeIndex = 0; chromosomeIndex < numberOfChromosomes; chromosomeIndex++){
        
        int* chromosome = chromosomeArray[chromosomeIndex]->citiesPermutation;//Pointer where chromosome are
        double totalDistance = 0;

        for(int citiesIndex = 1; citiesIndex < numberCities; citiesIndex++){
            
            //call the computeDistance function for these two cities and add to the totalDistance variable
            totalDistance += computeDistance(chromosome[citiesIndex - 1], chromosome[citiesIndex], citiesArray);
        }

        chromosomeArray[chromosomeIndex]->totalDistance = totalDistance;//store the total distance in the chromosome structure
    }
}//end of computeOverallDistance

//main function
int main(){
    
    int numberCities = readCities();//store the number of cities in the file
    struct city* citiesArray[numberCities];//array of structures of city 
    readInput(citiesArray, numberCities);//read the values of latitude and longitude in the file
    displayCities(citiesArray, numberCities);//Display the cities information.
    return 0;
}//end of the main function