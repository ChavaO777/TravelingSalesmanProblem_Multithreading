#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


//Esta función define el límite en cuanto al random del numero de ciudades
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
//Esta funcion genera la permutacion
void generateRandomPopulation(int *array, int amountOfCities)
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

//Esta funcion verifica si ya existe una permutacion existente en el arreglo
int populationAlreadyExists(int *newPopulation,int ** previousPopulations, int numberOfCities, int filas)
{
  int i =0;
  int j=0;
  int checksum = 0;
  int checktrue = 50;
  //Si el arreglo esta vacio, lo llena con la primera permutacion
  if((*(*previousPopulations+0)+0)==NULL)
  {
    for(i=0;i<numberOfCities;i++)
    {
      checktrue = newPopulation[i];
      previousPopulations[0][i] = checktrue;
    }
    return filas;
  }
  else
  {
    //Reasigna el tamaño del arreglo y le agrega filas
    int **tmp = realloc(previousPopulations, sizeof *previousPopulations * 2);
    if(tmp)
    {
      previousPopulations = tmp;
      for(size_t i = 0; i<1;i++)
      {
        previousPopulations[filas + i] = malloc(sizeof *previousPopulations[filas + i] * numberOfCities);
      }
      filas++;

    }
      for(i=0;i<numberOfCities;i++)
      {
        for(j=0;j<filas-1;j++)
        {
        if(previousPopulations[j][i]!=newPopulation[i])
        {
          checksum++;
        }
      }
      }
      if(checksum!=0)
      {
        for(i=0;i<numberOfCities;i++)
        {
          previousPopulations[filas-1][i] = newPopulation[i];
        }
      }
      return filas;
  }

}


int main(void)
{
  int i = 0;
  int filas = 1;
  int numberOfCities = 50; //ASigna 50 ciudades como ejemplo
  int *numbers; //Es la permutacion
  int **previousPopulations; //Arreglo de permutaciones
  previousPopulations = (int **)malloc(filas *sizeof(int *));
  for (i=0; i<filas; i++)
         previousPopulations[i] = (int *)malloc(numberOfCities * sizeof(int));
  numbers = malloc(numberOfCities * sizeof(int));
  for(i=0;i<numberOfCities;i++)
    numbers[i] = i;
  generateRandomPopulation(numbers,numberOfCities); //Las siguientes son pruebas de que funciona
  filas = populationAlreadyExists(numbers,previousPopulations,numberOfCities,filas);
  printf("\nArray after shuffling is: \n");
  for(i=0;i<numberOfCities;i++)
    printf("%d ",previousPopulations[0][i]);
  generateRandomPopulation(numbers,numberOfCities);
  filas = populationAlreadyExists(numbers,previousPopulations,numberOfCities,filas);
  generateRandomPopulation(numbers,numberOfCities);
  filas = populationAlreadyExists(numbers,previousPopulations,numberOfCities,filas);
  printf("\nArray after shuffling is:%d \n",filas);
  for(i=0;i<numberOfCities;i++)
    printf("%d ",previousPopulations[2][i]);
  return 0;
}
