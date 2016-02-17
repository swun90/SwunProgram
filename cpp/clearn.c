#include<stdio.h>

int main(){
	struct Student{
		char *name;
		int age;
		float score;
	};

	struct Student stu1;
	stu1.name="Abel";
	stu1.age=15;
	stu1.score=92.5;

	printf("%s\tage: %d, score: %f\n",stu1.name,stu1.age,stu1.score);
	return 0;
}