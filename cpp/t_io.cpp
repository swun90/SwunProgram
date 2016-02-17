#include<iostream>
#include<string>
#include<math.h>
// #include<fstream>
using namespace std;

#define PI 3.14159265
int main(int argc, char const *argv[])
{
	// string mystring = 
	// 	"This is a string. \\\"include<string>\"";
	// cout<<mystring<<endl;
	// int a=015;
	// int b=0x12;
	// double c=6.020e-23;
	// cout<<"a(015) = "<<a<<endl;
	// cout<<"b(0x12) = "<<b<<endl;
	// cout<<"c = "<<c<<endl;
	// double r=1.23;
	// cout<<"circle length is "<<2*PI*r<<endl;
	// int g=12%5;
	// int h=10;
	// h %= 3;
	// cout<<"12\%5 = "<<g<<endl;
	// cout<<" (10)h \%= 3 makes h = "<<h<<endl;
	// int a=2,b,c;
	// cout<<"a = "<<a<<endl;
	// b=a++;
	// cout<<"b=a++ makes b = "<<b<<endl;
	// a=2;
	// c=++a;
	// cout<<"c=++a makes c = "<<c<<endl;
	// a=2; b=3;
	// c=(a>b)?a:b;
	// cout<<c<<endl;

	// int a,b;
	// a=(b=3,b+2);
	// cout<<"b = "<<b<<"\n"<<"a = "<<a<<endl;
	
	// int a,b,c;
	// float p=3.14;
	// a=(int) p;
	// b=int (p);
	// cout<<"(int) p = "<<a<<endl;
	// cout<<"int (p) = "<<b<<endl;
	// c=sizeof(char);
	// cout<<"sizeof(char) = "<<c<<endl;
	// c=sizeof(int);
	// cout<<"sizeof(int) = "<<c<<endl;
	// c=sizeof(long);
	// cout<<"sizeof(long) = "<<c<<endl;
	// c=sizeof(long long);
	// cout<<"sizeof(long long) = "<<c<<endl;
	// c=sizeof(double);
	// cout<<"sizeof(double) = "<<c<<endl;
	// c=sizeof(float);
	// cout<<"sizeof(float) = "<<c<<endl;
	// c=sizeof(string);
	// cout<<"sizeof(string) = "<<c<<endl;

	string mystr ("1024");
	int myint;
	stringstream(mystr)>>myint;
	return 0;
}
// int main(int argc, char const *argv[])
// {
// 	int a[10], max, i, order;
// 	ofstream outfile("d1.dat",ios::out);
// 	if(!outfile){
// 		cerr<<"open error!"<<endl;
// 		exit(1);
// 	}
// 	cout<<"enter 10 integer numbers: "<<endl;
// 	for (int i = 0; i < 10; ++i)
// 	{
// 		cin>>a[10];
// 		outfile<<a[i]<<" ";
// 	}
// 	outfile.close();
// 	return 0;
// }

// int main(int argc, char const *argv[])
// {
// 	int a;
// 	cout<<"input number a: ";
// 	cin>>a;
// 	cout<<"dec: "<<dec<<a<<endl;
// 	cout<<"hex: "<<hex<<a<<endl;
// 	cout<<"oct: "<<setbase(8)<<a<<endl;
// 	const char *pt="China";
// 	cout<<setw(10)<<pt<<endl;
// 	cout<<setfill('*')<<setw(10)<<pt<<endl;
// 	double pi=22.0/7.0;
// 	cout<<setiosflags(ios::scientific)<<setprecision(8);
// 	cout<<"pi = "<<pi<<endl;
// 	cout<<"pi = "<<setprecision(4)<<pi<<endl;
// 	cout<<"pi = "<<setiosflags(ios::fixed)<<pi<<endl;

// 	return 0;
// }
// {
// 	const char *str="BASIC";
// 	for (int i = 4; i >= 0; --i)
// 	{
// 		cout.put(*(str+i));
// 	}
// 	cout<<endl;
// 	return 0;
// }
// {
// 	float grade;
// 	cout<<"enter grade: ";
// 	while(cin>>grade){
// 		if(grade>=85) cout<<grade<<"\tGood!"<<endl;
// 		if(grade<60) cout<<grade<<"\tFail!"<<endl;
// 		cout<<"enter next grade: ";
// 	}
// 	cout<<"The end!"<<endl;
// 	return 0;
// }
// {
// 	int c;
// 	cout<<"enter a sentence: "<<endl;
// 	while((c=cin.get())!=EOF) cout.put(c);
// 	return 0;
// }
// {
// 	char c;
// 	while(!cin.eof())
// 		if((c=cin.get())!=' ')
// 			cout.put(c);
// 	return 0;
// }