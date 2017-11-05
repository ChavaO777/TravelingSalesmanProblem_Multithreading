#include <stdio.h>
#include <math.h>
#define EARTH_RADIUS_KM 6371 // authalic radius based on/extracted from surface area;
#define radianes (M_PI / 180.0) //approximation of the radius of the average circumference

/** Structure of an individual. It includes an array of integers representing
*   the permutation of the n cities (0 to n - 1) that form that solution and 
*   the total distance traveled to visit all n cities in that order and come
*   back to the first city.
*/ 
struct city{
    int id;//id of the city
    double latitude;//latitude of the cities
    double longitude;//longitude  of the cities
};//end of the city struct

struct individual {
    
    int* citiesPermutation; // a pointer to an integer corresponding to the permutation of the n cities
    double totalDistance; // the total distance traveled
};//end of the individual struct

//Function that reads the input
void readInput(struct city citiesArray[], int numberCities){
  
    for(int i = 0; i < numberCities; i++){
        
        citiesArray[i].id = i + 1;
        scanf("%lf %lf", &citiesArray[i].latitude, &citiesArray[i].longitude);
    }
}//end of the readInput function 

//function that reads the number of cities
int readCities(){
    
    int totalCities;
    scanf("%d", &totalCities);
    return totalCities;
}

void displayCities(struct city citiesArray[], int numberCities){

    printf("\n");

    for(int i = 0; i < numberCities; i++){

        printf("City id: %d; Latitude: %lf; Longitude: %f\n", citiesArray[i].id, citiesArray[i].latitude, citiesArray[i].longitude);
    }

    printf("\n");
}

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
double computeDistance(int indexCoord1, int indexCoord2, struct city citiesArray[]){

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

/*struct individual* solution computeOverallDistance(int* population) {


}  */ 

//main function
int main(){
    
    int numberCities = readCities();//store the number of cities in the file
    struct city citiesArray[numberCities];//array of structures of city 
    readInput(citiesArray, numberCities);//read the values of latitude and longitude in the file
    displayCities(citiesArray, numberCities);
    
    return 0;
}//end of the main function