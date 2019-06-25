#ifndef HW3
#define HW3

#include <iostream>
#include <cmath>
#include <vector>
#include <GL/glut.h>
#include <stdarg.h>
#include <numeric>
#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <fstream>

#define LEN 8192  //  Maximum length of text string

#define Cos(th) cos(M_PI/180*(th))
#define Sin(th) sin(M_PI/180*(th))

extern "C" {
unsigned int LoadTexBMP(const char* file);
void Fatal(const char* format , ...);
void ErrCheck(const char* where);
}

#endif
