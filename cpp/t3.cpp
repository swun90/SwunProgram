#include<iostream>
using namespace std;

int addition_1(int a, int b);
int factorial_1(int n);
void dupli_1(int &a);

int main()
{
	int z,x=2;
	long long f1;
	f1=factorial_1(10);
	z=addition_1(4,3);
	cout<<"z is "<<z<<endl;
	cout<<"factorial(10) is  "<<f1<<endl;
	cout<<"original x is "<<x<<endl;
	dupli_1(x);
	cout<<"after a function having variable & before, x is "
		<<x<<endl;
	return 0;
}

int addition_1(int a, int b)
{
	return a+b;
}

int factorial_1(int n)
{
	if (n==0 || n==1)
		return 1;
	else
		return n*factorial_1(n-1);
}
void dupli_1(int &a)
{
	a *= 2;
}
