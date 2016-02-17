#include<iostream>
using namespace std;

class Student{
private:
	const char *name;
	int age;
	float score;

public:
	Student(const char *,int,float);
	~Student();
	void say();
};

// Student::Student(const char *nam1, int ag1, float scor1){
// 	name=nam1;
// 	age=ag1;
// 	score=scor1;
// }
Student::Student(const char *nam1, int ag1, float scor1):
name(nam1),age(ag1),score(scor1){}
Student::~Student(){
	cout<<name<<" goodbye"<<endl;
}
void Student::say(){
	cout<<name<<"\tage: "<<age<<"\tscore: "<<score<<endl;
}

int main(int argc, char const *argv[])
{
	Student stu1("Kate",16,87.5);
	stu1.say();

	Student stu2("Bella",15,92.5);
	stu2.say();
	return 0;
}
// int main(){
// 	int sum = 0;
// 	int val = 0;
// 	cout<<"(input non-integer to stop input!)"<<endl;
// 	cout<<"input a number:\t";
// 	while(cin>>val){
// 		sum += val;
// 		cout<<"input next number:\t";
// 	}
// 	cout<<"the sum of all input number is "<<sum<<endl;

// 	return 0;
// }

// int valplus(int &a);

// int main(int argc, char const *argv[])
// {
// 	int n1=10;
// 	int n2;
// 	n2=valplus(n1);
// 	cout<<"n1 = "<<n1<<"\nn2 = "<<n2<<endl;
// 	return 0;
// }

// int valplus(int &a){
// 	a = a+5;
// 	return a;
// }

// class Student{
// private:
// 	const char *name;
// 	int age;
// 	float score;

// public:
// 	void setname(const char *name1){
// 		name=name1;
// 	}
// 	void setage(int age1){
// 		age=age1;
// 	}
// 	void setscore(float score1){
// 		score=score1;
// 	}
// 	void say(){
// 		cout<<name<<"\tage: "<<age<<"\tscore: "<<score<<endl;
// 	}
// };

// int main(int argc, char const *argv[])
// {
// 	Student stu1;
// 	stu1.setname("Bella");
// 	stu1.setage(15);
// 	stu1.setscore(92.5f);
// 	stu1.say();
// 	return 0;
// }