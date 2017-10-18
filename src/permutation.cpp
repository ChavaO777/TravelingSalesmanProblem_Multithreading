// permutation.cpp

#include <iostream>
using namespace std;

int main(){

	int i, j, k, t;
	unsigned long n = 3;

	int* p = new int[n + 1];

	// starting permutation
	// identity 1, 2, 3, ... , n -> 1, 2, ... , n

	for(int i = 1; i <= n; i++){

		p[i] = i;
		// cout<<"p["<<i<<"] = "<<p[i]<<" ";
		cout<<p[i]<<" ";
	}

	cout<<endl;

	int test = 1;
	do{

		i = n - 1;

		while(p[i] > p[i + 1])
			i = i - 1;

		if(i > 0)
			test = 1;

		else
			test = 0;

		j = n;

		while(p[j] <= p[i])
			j = j - 1;

		t = p[i];
		p[i] = p[j];
		p[j] = t;
		i = i + 1;
		j = n;

		while(i < j){

			t = p[i];
			p[i] = p[j];
			p[j] = t;
			i = i + 1;
			j = j - 1;
		}

		//display result

		for(int tau = 1; tau <= n; tau++){

			if(p[tau] == 0){

				test = 0;
				break;
			}

			// cout<<"p["<<tau<<"] = "<<p[tau]<<" ";
			cout<<p[tau]<<" ";
		}

		if(test)
			cout<<endl;

	}while(test == 1);

	delete[] p;

	return 0;
}







