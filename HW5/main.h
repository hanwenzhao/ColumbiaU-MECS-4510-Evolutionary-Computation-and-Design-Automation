//
// Created by Hanwen Zhao on 12/8/18.
// Copyright (c) 2018 Hanwen Zhao. All rights reserved.
//

#ifndef HW5_MAIN_H
#define HW5_MAIN_H

#include <iostream>
#include <GL/freeglut.h>
#include <GL/glut.h>
#include <vector>
#include <stdarg.h>
#include <fstream>
#include "omp.h"
#include <algorithm>
#include <math.h>
#include <cmath>
#include <numeric>
#include <random>
#include <queue>
#include <ctime>
#include <chrono>
#include <string>

PFNGLWINDOWPOS2IPROC glWindowPos2i;

#define LEN 8192  //  Maximum length of text string

#define Cos(th) cos(M_PI/180*(th))
#define Sin(th) sin(M_PI/180*(th))

struct MASS
{
    double m;       // mass
    double p[3];    // 3D position
    double v[3];    // 3D velocity
    double a[3];    // 3D acceleration
    int springFactor;
};

struct SPRING
{
    double k;       // spring constant
    double L_0;     // rest length
    double L;
    int m1;         // first mass connected
    int m2;         // second mass connected
    int type;
};

struct GENE
{
    bool exist;
    int springFactor;
    std::vector<double> p;
};

#endif //HW5_MAIN_H
