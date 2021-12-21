#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#define PI 3.14159265
#define N 2



using namespace std;

struct float12 
{
   float f01;
   float f02;
   float f03; 
   float f1;
   float f2; 
   float f3;
   float f4;
   float f5;
   float f6; 
   float f7;
   float f8;
   float f9;     
};

struct uint0 
{
   short int u0;
};

struct uint1 
{
   unsigned int u;
};

   uint0 pack00;
   uint1 pack0;
   float12 pack;


  struct S {
    double a[N];
    double b[N];
    double c[N];
  };

vector<S> vec;

int main(int argc, char* argv[]) {
std::ifstream infile(argv[1]);

//Vector example
float x1, x2, y1, y2, z1, z2;
 S s1;
 while (infile >> x1 >> x2 >> y1 >> y2 >> z1 >> z2)
{
      s1.a[0] ={x1};
  s1.a[1] ={x2};
  s1.b[0] ={y1};
  s1.b[1] ={y2};
  s1.c[0] ={z1};
  s1.c[1] ={z2};
  vec.push_back(s1);
}
//  s1.a[0] ={0};
//  s1.a[1] ={10};
//  s1.b[0] ={0};
//  s1.b[1] ={10};
//  s1.c[0] ={0};
//  s1.c[1] ={10};
//  vec.push_back(s1);
//  s1.a[0] ={20};
//  s1.a[1] ={30};
//  vec.push_back(s1);
//  s1.a[0] ={40};
//  s1.a[1] ={50};
//  vec.push_back(s1);


int boxes_amount=vec.size()*12;



//Вставьте свой путь
string fullPath = argv[2];


FILE *j;
string b;
j = fopen(fullPath.c_str(),"w");
ofstream dout (fullPath);

      pack00.u0=0; 
      pack0.u=0;    
      for (int i=0; i<20; i++)
   {
   fwrite(&pack0, sizeof(uint1), 1, j);  
   } 
   
   pack0.u=boxes_amount;   
   fwrite(&pack0, sizeof(uint1), 1, j);  

for (int i = 0; i < vec.size(); i++)
	{ 

float flx1=vec[i].a[0];
float flx2=vec[i].a[1];
float fly1=vec[i].b[0];
float fly2=vec[i].b[1];
float flz1=vec[i].c[0];
float flz2=vec[i].c[1];
   
   pack.f01=1;
   pack.f02=0;
   pack.f03=0;     
   pack.f1=flx2;
   pack.f2=fly2;
   pack.f3=flz2;
   pack.f4=flx2;
   pack.f5=fly1;
   pack.f6=flz2;
   pack.f7=flx2;
   pack.f8=fly2;   
   pack.f9=flz1; 
   fwrite(&pack, sizeof(float12), 1, j);
   fwrite(&pack00, sizeof(uint0), 1, j); 

   pack.f1=flx2;
   pack.f2=fly1;
   pack.f3=flz1;
   fwrite(&pack, sizeof(float12), 1, j);    
   fwrite(&pack00, sizeof(uint0), 1, j);   

      
   pack.f01=-1;
   pack.f02=0;
   pack.f03=0;     
   pack.f1=flx1;
   pack.f2=fly2;
   pack.f3=flz1;
   pack.f4=flx1;
   pack.f5=fly2;
   pack.f6=flz2;
   pack.f7=flx1;
   pack.f8=fly1;   
   pack.f9=flz1; 
   fwrite(&pack, sizeof(float12), 1, j);
   fwrite(&pack00, sizeof(uint0), 1, j);  
   
   pack.f1=flx1;
   pack.f2=fly1;
   pack.f3=flz2;
   fwrite(&pack, sizeof(float12), 1, j);
   fwrite(&pack00, sizeof(uint0), 1, j);  
      
   pack.f01=0;
   pack.f02=-1;
   pack.f03=0; 
   pack.f1=flx1;
   pack.f2=fly1;
   pack.f3=flz1; 
   pack.f4=flx1;
   pack.f5=fly1;
   pack.f6=flz2;
   pack.f7=flx2;
   pack.f8=fly1;   
   pack.f9=flz1; 
   fwrite(&pack, sizeof(float12), 1, j);
   fwrite(&pack00, sizeof(uint0), 1, j);
   
   pack.f1=flx2;
   pack.f2=fly1;
   pack.f3=flz2;
   fwrite(&pack, sizeof(float12), 1, j);
   fwrite(&pack00, sizeof(uint0), 1, j);   
      
   pack.f01=0;
   pack.f02=1;
   pack.f03=0;     
   pack.f1=flx2;
   pack.f2=fly2;
   pack.f3=flz1;
   pack.f4=flx2;
   pack.f5=fly2;
   pack.f6=flz2;
   pack.f7=flx1;
   pack.f8=fly2;   
   pack.f9=flz1; 
   fwrite(&pack, sizeof(float12), 1, j);
   fwrite(&pack00, sizeof(uint0), 1, j);  
      
   pack.f1=flx1;
   pack.f2=fly2;
   pack.f3=flz2;
   fwrite(&pack, sizeof(float12), 1, j);
   fwrite(&pack00, sizeof(uint0), 1, j);  
      
   pack.f01=0;
   pack.f02=0;
   pack.f03=-1;
   pack.f1=flx1;
   pack.f2=fly2;
   pack.f3=flz1;
   pack.f4=flx1;
   pack.f5=fly1;
   pack.f6=flz1;
   pack.f7=flx2;
   pack.f8=fly2;   
   pack.f9=flz1; 
   fwrite(&pack, sizeof(float12), 1, j);
   fwrite(&pack00, sizeof(uint0), 1, j);     
   
   pack.f1=flx2;
   pack.f2=fly1;
   pack.f3=flz1;
   fwrite(&pack, sizeof(float12), 1, j);
   fwrite(&pack00, sizeof(uint0), 1, j);
      
   pack.f01=0;
   pack.f02=0;
   pack.f03=1; 
   pack.f1=flx1;
   pack.f2=fly2;
   pack.f3=flz2;
   pack.f4=flx2;
   pack.f5=fly2;
   pack.f6=flz2;
   pack.f7=flx1;
   pack.f8=fly1;   
   pack.f9=flz2; 
   fwrite(&pack, sizeof(float12), 1, j);
   fwrite(&pack00, sizeof(uint0), 1, j);  
   
   pack.f1=flx2;
   pack.f2=fly1;
   pack.f3=flz2;
   fwrite(&pack, sizeof(float12), 1, j);
   fwrite(&pack00, sizeof(uint0), 1, j);
   
					}
fclose(j);

  return 0;
}
