#include <stdio.h>
#include <math.h>
#define radius 6371 // authalic radius based on/extracted from surface area;
#define radianes (3.1415926536 / 180) //approximation of the radius of the average circumference
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
void readInput(struct city citys[], int numberCities){
  //  int length = sizeof(citys)/sizeof(struct);
    for(int i = 0; i < numberCities; i++){
        citys[i].id = i+1;
        scanf("%lf %lf", &citys[i].latitude, &citys[i].longitude);
    }
}//end of the readInput function 
//fnction that reads the number of cities
int readCitys(){
    int totalCities;
    scanf("%d", &totalCities);
    return totalCities;
}
/**
*   Computes the distance between two cities given their indices. It is an implementation
*   of the Haversine formula (https://en.wikipedia.org/wiki/Haversine_formula)
*/
double computeDistance(int indexCoord1, int indexCoord2, struct city cities[]){

    double part1, part2, part3, longitude1,longitude2 ,latitude1,latitude2;
    longitude1 = cities[indexCoord1].longitude; 
    longitude2 = cities[indexCoord2].longitude;
    latitude1 = cities[indexCoord1].latitude;
    latitude2 = cities[indexCoord2].latitude;
    longitude1 *= radianes;//transform to radians
    longitude2 *= radianes;//transform to radians
    latitude1*= radianes; //transform to radians
    latitude2 *= radianes;//transform to radians
    double latitude = latitude2 -latitude1;//diferential betwen latitudes
    double longitude = longitude2 -longitude1;//diferential betwen longitudes
    //Haversine formula
	return asin(sqrt(pow(sin(latitude/2),2)+pow(sin(longitude/2),2)*cos(latitude1)*cos(latitude2))) * 2 * radius;
   
}//end of the computeDistance function

/*struct individual* solution computeOverallDistance(int* population) {


}  */ 

//main function
int main(){
    int numberCities = readCitys();//store the number of cities in the file
    struct city cities[numberCities];//array of structures of city 
    readInput(cities,numberCities);//read the values of latitude and longitude in the file
    return 0;
}//end of the main function