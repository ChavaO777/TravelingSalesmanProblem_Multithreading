#include <stdio.h>
#include <math.h>
int totalCitys = 0; //total de la ciudades que ahy en el archivo
 double latitude[100];//array de las pociciones en x
 double longitude[100];//array de as pociciones en y
//sructura de los chromosomas
struct individual {
    double chromosoma[100];//se guarda aqui un chromosama aqui que es la solucion
    double distance;//aqui se guarda la distancia total de la solucion
 };//fin strutura
//funcion que le el archivo
void readInput(){
   // printf("dame la ciudades");//esqeuletopara formar el esqueleto
    scanf("%d",&totalCitys);
    int scan = 0;//contador para el for
   while(scan < totalCitys){
        scanf("%lf",&latitude[scan]);
       // printf("dame la ciudades");//prntprueba
        scanf("%lf",&longitude[scan]);
        scan++; 
    }
};//fin del readInput
//computa la distancia entre los dos puntos tomando la pocicion de los dos array longitud y latitud
double computeDistance(int indexCoord1, int indexCoord2){
    double result = sqrt(pow(latitude[indexCoord2]-latitude[indexCoord1],2)+pow(longitude[indexCoord2]-longitude[indexCoord1],2));
    return result;
}//fin computeDistance
struct individual* solution computeOverallDistance(int* population) {

}
//main
int main()
{
    readInput();
//printf("%Lf", latitude[0]);//print de pruebas
}//fin Main