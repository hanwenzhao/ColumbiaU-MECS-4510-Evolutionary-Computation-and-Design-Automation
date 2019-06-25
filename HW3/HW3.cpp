#include "HW3.h"

// OPENGL Variables
int th = 0;            //  Azimuth of view angle
int ph = 0;           //  Elevation of view angle
int axes = 0;         //  Display axes
int light = 1;
double asp = 1;     // aspect ratio
int fov = 45;         //  Field of view (for perspective)
double dim = 1.0;  // size of the workd
double skyBoxScale = 0.5;

int emission  = 100;  // Emission intensity (%)
int ambient   = 100;  // Ambient intensity (%)
int diffuse   = 100;  // Diffuse intensity (%)
int specular  = 100;  // Specular intensity (%)
int shininess = 128;  // Shininess (power of two)
float shiny   =   1;    // Shininess (value)
float white[] = {1,1,1,1};
float black[] = {0,0,0,1};

unsigned int grassTexture;
unsigned int slimeTexture;
unsigned int skyBoxTexture[10]; // Texture for Sky Box

// Physics Simluator Variables
double mass = 0.1;
double length = 0.1;
double gravity = 9.81;
double T = 0;
double springPotentialEnergy = 0;
double gravityPotentialEnergy = 0;
double totalPotentialEnergy = 0;
double kineticEnergy = 0;
double totalEnergy = 0;
double groundEnergy = 0;

double timeStep = 0.001;
double restoreConstant = 10000;
double springConstant = 1000;
double dampingConstant = 0.999;
double frictionCoefficient = 0;//0.8;

GLfloat worldRotation[16] = {1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1};
std::clock_t start = std::clock();
double duration;

std::ofstream myFile("breathing.txt");

struct MASS
{
    double m;       // mass
    double p[3];    // 3D position
    double v[3];    // 3D velocity
    double a[3];    // 3D acceleration
};

struct SPRING
{
    double k;       // spring constant
    double L_0;     // rest length
    int m1;         // first mass connected
    int m2;         // second mass connected
};

std::vector<MASS> generateCubeMasses(double mass, double length, double initialX, double initialY, double initialZ)
{
    std::vector<MASS> cubeMasses(8);
    cubeMasses[0] = {mass,{initialX+length/2,initialY+length/2,initialZ},{0,0,0},{0,0,0}};
    cubeMasses[1] = {mass,{initialX+length/2,initialY-length/2,initialZ},{0,0,0},{0,0,0}};
    cubeMasses[2] = {mass,{initialX-length/2,initialY-length/2,initialZ},{0,0,0},{0,0,0}};
    cubeMasses[3] = {mass,{initialX-length/2,initialY+length/2,initialZ},{0,0,0},{0,0,0}};
    cubeMasses[4] = {mass,{initialX+length/2,initialY+length/2,initialZ+length},{0,0,0},{0,0,0}};
    cubeMasses[5] = {mass,{initialX+length/2,initialY-length/2,initialZ+length},{0,0,0},{0,0,0}};
    cubeMasses[6] = {mass,{initialX-length/2,initialY-length/2,initialZ+length},{0,0,0},{0,0,0}};
    cubeMasses[7] = {mass,{initialX-length/2,initialY+length/2,initialZ+length},{0,0,0},{0,0,0}};
    return cubeMasses;
}

std::vector<SPRING> generateCubeSprings(double springConstant)
{
    double length = 0.1;
    double diagonal2 = sqrt(pow(length,2)+pow(length,2));
    double diagonal3 = sqrt(pow(length,2)+pow(length,2)+pow(length,2));
    std::vector<SPRING> cubeSprings(28);
    cubeSprings[0] = {springConstant,length,   0,1};
    cubeSprings[1] = {springConstant,diagonal2,0,2};
    cubeSprings[2] = {springConstant,length,   0,3};
    cubeSprings[3] = {springConstant,length,   0,4};
    cubeSprings[4] = {springConstant,diagonal2,0,5};
    cubeSprings[5] = {springConstant,diagonal3,0,6};
    cubeSprings[6] = {springConstant,diagonal2,0,7};

    cubeSprings[7] = {springConstant,length,    1,2};
    cubeSprings[8] = {springConstant,diagonal2, 1,3};
    cubeSprings[9] = {springConstant,diagonal2, 1,4};
    cubeSprings[10] = {springConstant,length,   1,5};
    cubeSprings[11] = {springConstant,diagonal2,1,6};
    cubeSprings[12] = {springConstant,diagonal3,1,7};

    cubeSprings[13] = {springConstant,length,   2,3};
    cubeSprings[14] = {springConstant,diagonal3,2,4};
    cubeSprings[15] = {springConstant,diagonal2,2,5};
    cubeSprings[16] = {springConstant,length,   2,6};
    cubeSprings[17] = {springConstant,diagonal2,2,7};

    cubeSprings[18] = {springConstant,diagonal2,3,4};
    cubeSprings[19] = {springConstant,diagonal3,3,5};
    cubeSprings[20] = {springConstant,diagonal2,3,6};
    cubeSprings[21] = {springConstant,length,   3,7};

    cubeSprings[22] = {springConstant,length,   4,5};
    cubeSprings[23] = {springConstant,diagonal2,4,6};
    cubeSprings[24] = {springConstant,length,   4,7};

    cubeSprings[25] = {springConstant,length,   5,6};
    cubeSprings[26] = {springConstant,diagonal2,5,7};

    cubeSprings[27] = {springConstant,length,6,7};
    return cubeSprings;
}

// create a list for masses
std::vector<MASS> cubeMasses = generateCubeMasses(mass, length, 0.0, 0.0, 0.5);
// assign values as a cube
std::vector<SPRING> cubeSprings = generateCubeSprings(springConstant);

// calcualte norm for vector
double norm( double x[], std::size_t sz )
{
    return std::sqrt( std::inner_product( x, x+sz, x, 0.0 ) ) ;
}

void drawCube(std::vector<MASS>& cubeMasses, std::vector<SPRING>& cubeSprings)
{
    glColor3f(1,0,0);
    /*
    GLUquadric *quad;
    quad = gluNewQuadric();
    for (int i = 0; i < (int)cubeMasses.size(); i++)
    {
        glPushMatrix();
        glMultMatrixf(worldRotation);
        glTranslated(cubeMasses[i].p[0],cubeMasses[i].p[1],cubeMasses[i].p[2]);
        gluSphere(quad,length/10,10,10);
        glPopMatrix();
    }

    for (int i = 0; i < (int)cubeSprings.size(); i++)
    {
        glPushMatrix();
        glMultMatrixf(worldRotation);
        glBegin(GL_LINES);
        glVertex3f(cubeMasses[cubeSprings[i].m1].p[0],cubeMasses[cubeSprings[i].m1].p[1],cubeMasses[cubeSprings[i].m1].p[2]);
        glVertex3f(cubeMasses[cubeSprings[i].m2].p[0],cubeMasses[cubeSprings[i].m2].p[1],cubeMasses[cubeSprings[i].m2].p[2]);
        glEnd();
        glPopMatrix();
    }
    */
    glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
    glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);
    //  Save transformation
    glPushMatrix();
    glMultMatrixf(worldRotation);
    // Enable Texture
    glEnable(GL_TEXTURE_2D);
    glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
    glColor3f(1,1,1);
    glBindTexture(GL_TEXTURE_2D,slimeTexture);

    //  Front
    glColor3f(1,1,1);
    glBindTexture(GL_TEXTURE_2D,slimeTexture);
    glBegin(GL_QUADS);
    glNormal3f( 0, 0, 1);
    glTexCoord2f(0.0f,0.0f);    glVertex3f(cubeMasses[2].p[0],cubeMasses[2].p[1],cubeMasses[2].p[2]);
    glTexCoord2f(1.0f,0.0f);    glVertex3f(cubeMasses[3].p[0],cubeMasses[3].p[1],cubeMasses[3].p[2]);
    glTexCoord2f(1.0f,1.0f);    glVertex3f(cubeMasses[7].p[0],cubeMasses[7].p[1],cubeMasses[7].p[2]);
    glTexCoord2f(0.0f,1.0f);    glVertex3f(cubeMasses[6].p[0],cubeMasses[6].p[1],cubeMasses[6].p[2]);
    glEnd();

    //  Back
    glColor3f(1,1,1);
    glBindTexture(GL_TEXTURE_2D,slimeTexture);
    glBegin(GL_QUADS);
    glNormal3f( 0, 0, -1);
    glTexCoord2f(0.0f,0.0f);    glVertex3f(cubeMasses[0].p[0],cubeMasses[0].p[1],cubeMasses[0].p[2]);
    glTexCoord2f(1.0f,0.0f);    glVertex3f(cubeMasses[1].p[0],cubeMasses[1].p[1],cubeMasses[1].p[2]);
    glTexCoord2f(1.0f,1.0f);    glVertex3f(cubeMasses[5].p[0],cubeMasses[5].p[1],cubeMasses[5].p[2]);
    glTexCoord2f(0.0f,1.0f);    glVertex3f(cubeMasses[4].p[0],cubeMasses[4].p[1],cubeMasses[4].p[2]);
    glEnd();

    //  Right
    glColor3f(1,1,1);
    glBindTexture(GL_TEXTURE_2D,slimeTexture);
    glBegin(GL_QUADS);
    glNormal3f( 1, 0, 0);
    glTexCoord2f(0.0f,0.0f);    glVertex3f(cubeMasses[0].p[0],cubeMasses[0].p[1],cubeMasses[0].p[2]);
    glTexCoord2f(1.0f,0.0f);    glVertex3f(cubeMasses[3].p[0],cubeMasses[3].p[1],cubeMasses[3].p[2]);
    glTexCoord2f(1.0f,1.0f);    glVertex3f(cubeMasses[7].p[0],cubeMasses[7].p[1],cubeMasses[7].p[2]);
    glTexCoord2f(0.0f,1.0f);    glVertex3f(cubeMasses[4].p[0],cubeMasses[4].p[1],cubeMasses[4].p[2]);
    glEnd();


    //  Left
    glColor3f(1,1,1);
    glBindTexture(GL_TEXTURE_2D,slimeTexture);
    glBegin(GL_QUADS);
    glNormal3f( -1, 0, 0);
    glTexCoord2f(0.0f,0.0f);    glVertex3f(cubeMasses[2].p[0],cubeMasses[2].p[1],cubeMasses[2].p[2]);
    glTexCoord2f(1.0f,0.0f);    glVertex3f(cubeMasses[1].p[0],cubeMasses[1].p[1],cubeMasses[1].p[2]);
    glTexCoord2f(1.0f,1.0f);    glVertex3f(cubeMasses[5].p[0],cubeMasses[5].p[1],cubeMasses[5].p[2]);
    glTexCoord2f(0.0f,1.0f);    glVertex3f(cubeMasses[6].p[0],cubeMasses[6].p[1],cubeMasses[6].p[2]);
    glEnd();


    //  Top
    glColor3f(1,1,1);
    glBindTexture(GL_TEXTURE_2D,slimeTexture);
    glBegin(GL_QUADS);
    glNormal3f( 0, 1, 0);
    glTexCoord2f(0.0f,0.0f);    glVertex3f(cubeMasses[4].p[0],cubeMasses[4].p[1],cubeMasses[4].p[2]);
    glTexCoord2f(1.0f,0.0f);    glVertex3f(cubeMasses[5].p[0],cubeMasses[5].p[1],cubeMasses[5].p[2]);
    glTexCoord2f(1.0f,1.0f);    glVertex3f(cubeMasses[6].p[0],cubeMasses[6].p[1],cubeMasses[6].p[2]);
    glTexCoord2f(0.0f,1.0f);    glVertex3f(cubeMasses[7].p[0],cubeMasses[7].p[1],cubeMasses[7].p[2]);
    glEnd();

    //  Bottom
    glColor3f(1,1,1);
    glBindTexture(GL_TEXTURE_2D,slimeTexture);
    glBegin(GL_QUADS);
    glNormal3f( 0, -1, 0);
    glTexCoord2f(0.0f,0.0f);    glVertex3f(cubeMasses[0].p[0],cubeMasses[0].p[1],cubeMasses[0].p[2]);
    glTexCoord2f(1.0f,0.0f);    glVertex3f(cubeMasses[1].p[0],cubeMasses[1].p[1],cubeMasses[1].p[2]);
    glTexCoord2f(1.0f,1.0f);    glVertex3f(cubeMasses[2].p[0],cubeMasses[2].p[1],cubeMasses[2].p[2]);
    glTexCoord2f(0.0f,1.0f);    glVertex3f(cubeMasses[3].p[0],cubeMasses[3].p[1],cubeMasses[3].p[2]);
    glEnd();

    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
    
}

void cubeUpdate(std::vector<MASS>& cubeMasses, std::vector<SPRING>& cubeSprings, int option)
{
    // initialize the force vector with value 0
    std::vector<std::vector<double>> cubeForces((int)cubeMasses.size(),std::vector<double>(3));
    // loop through all springs to calculate spring forces
    for (int i = 0; i < (int)cubeSprings.size(); i++)
    {
        // option 1 to have breathing cube
        if (option == 1)
        {
            if (T > 0.2)
            {
              // 5 12 14 19
                cubeSprings[0].L_0 = 1.0 * length + 0.8 * length * sin(10*T);
                //cubeSprings[13].L_0 = 1.0 * length + 0.5 * length * cos(100*T);
                //cubeSprings[9].L_0 = 1.0 * length + 0.5 * length * sin(T);
                //cubeSprings[22].L_0 = 1.0 * length + 0.5 * length * sin(T);
            }
        }
        MASS mass1 = cubeMasses[cubeSprings[i].m1];
        MASS mass2 = cubeMasses[cubeSprings[i].m2];
        double positionDiff[3] = {mass2.p[0] - mass1.p[0],mass2.p[1] - mass1.p[1],mass2.p[2] - mass1.p[2]};
        double L = norm(positionDiff,3);
        double force = cubeSprings[i].k * fabs(cubeSprings[i].L_0 - L);
        double direstion[3] = {positionDiff[0]/L,positionDiff[1]/L,positionDiff[2]/L};
        // contraction case
        if (L > cubeSprings[i].L_0)
        {
            cubeForces[cubeSprings[i].m1][0] = cubeForces[cubeSprings[i].m1][0] + direstion[0]*force;
            cubeForces[cubeSprings[i].m1][1] = cubeForces[cubeSprings[i].m1][1] + direstion[1]*force;
            cubeForces[cubeSprings[i].m1][2] = cubeForces[cubeSprings[i].m1][2] + direstion[2]*force;
            cubeForces[cubeSprings[i].m2][0] = cubeForces[cubeSprings[i].m2][0] - direstion[0]*force;
            cubeForces[cubeSprings[i].m2][1] = cubeForces[cubeSprings[i].m2][1] - direstion[1]*force;
            cubeForces[cubeSprings[i].m2][2] = cubeForces[cubeSprings[i].m2][2] - direstion[2]*force;
        }
        // expansion case
        else if (L < cubeSprings[i].L_0)
        {
            cubeForces[cubeSprings[i].m1][0] = cubeForces[cubeSprings[i].m1][0] - direstion[0]*force;
            cubeForces[cubeSprings[i].m1][1] = cubeForces[cubeSprings[i].m1][1] - direstion[1]*force;
            cubeForces[cubeSprings[i].m1][2] = cubeForces[cubeSprings[i].m1][2] - direstion[2]*force;
            cubeForces[cubeSprings[i].m2][0] = cubeForces[cubeSprings[i].m2][0] + direstion[0]*force;
            cubeForces[cubeSprings[i].m2][1] = cubeForces[cubeSprings[i].m2][1] + direstion[1]*force;
            cubeForces[cubeSprings[i].m2][2] = cubeForces[cubeSprings[i].m2][2] + direstion[2]*force;
        }
        // calcualte spring potential energy
        springPotentialEnergy = springPotentialEnergy + cubeSprings[i].k * pow((L - cubeSprings[i].L_0),2)/2;
    }
    // loop through all masses
    for (int i = 0; i < (int)cubeMasses.size(); i++)
    {
        // add gravity
        cubeForces[i][2] = cubeForces[i][2] - cubeMasses[i].m * gravity;
        // if the mass is below ground, add restroration force and calculate friction
        if (cubeMasses[i].p[2] <= 0)
        {
            cubeForces[i][2] = cubeForces[i][2] + restoreConstant * fabs(cubeMasses[i].p[2]);
            groundEnergy = groundEnergy + restoreConstant * pow(cubeMasses[i].p[2],2)/2;
            // calculate horizontal force and vertical force
            double F_h = sqrt(pow(cubeForces[i][0],2)+pow(cubeForces[i][1],2));
            double F_v = cubeForces[i][2];
            if (F_h < F_v * frictionCoefficient)
            {
                cubeForces[i][0] = 0;
                cubeForces[i][1] = 0;
                cubeMasses[i].v[0] = 0;
                cubeMasses[i].v[1] = 0;
            }
            else
            {
                for (int j = 0; j < 2; j ++)
                {
                    if (cubeForces[i][j] < 0)
                    {
                        cubeForces[i][j] = cubeForces[i][j] + F_v * frictionCoefficient * cubeForces[i][j]/F_h;
                        if (cubeForces[i][j] > 0){cubeForces[i][j] = 0;}
                    }
                    else
                    {
                        cubeForces[i][j] = cubeForces[i][j] - F_v * frictionCoefficient * cubeForces[i][j]/F_h;
                        if (cubeForces[i][j] < 0){cubeForces[i][j] = 0;}
                    }
                }
            }
        }
        // update acceleration
        cubeMasses[i].a[0] = cubeForces[i][0]/cubeMasses[i].m;
        cubeMasses[i].a[1] = cubeForces[i][1]/cubeMasses[i].m;
        cubeMasses[i].a[2] = cubeForces[i][2]/cubeMasses[i].m;
        // update velocity
        cubeMasses[i].v[0] = dampingConstant * (cubeMasses[i].v[0] + cubeMasses[i].a[0] * timeStep);
        cubeMasses[i].v[1] = dampingConstant * (cubeMasses[i].v[1] + cubeMasses[i].a[1] * timeStep);
        cubeMasses[i].v[2] = dampingConstant * (cubeMasses[i].v[2] + cubeMasses[i].a[2] * timeStep);
        // update position
        cubeMasses[i].p[0] = cubeMasses[i].p[0] + cubeMasses[i].v[0] * timeStep;
        cubeMasses[i].p[1] = cubeMasses[i].p[1] + cubeMasses[i].v[1] * timeStep;
        cubeMasses[i].p[2] = cubeMasses[i].p[2] + cubeMasses[i].v[2] * timeStep;
        // calculate gravity potential energy
        gravityPotentialEnergy = gravityPotentialEnergy + cubeMasses[i].m * gravity * cubeMasses[i].p[2];
        // calculate kinetic energy
        kineticEnergy = kineticEnergy + cubeMasses[i].m * pow(norm(cubeMasses[i].v,3),2) / 2;
    }
    totalPotentialEnergy = springPotentialEnergy + gravityPotentialEnergy + groundEnergy;
    totalEnergy = totalPotentialEnergy + kineticEnergy;
    std::cout << "Time: " << T << std::endl;
    std::cout << "Total Kinetic Energy: " << kineticEnergy << std::endl;
    std::cout << "Total Potential Energy: " << totalPotentialEnergy << std::endl;
    std::cout << "Total Energy: " << totalEnergy << std::endl;
    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    myFile << duration << " " << T << " " << totalPotentialEnergy << " " << kineticEnergy << " " << totalEnergy << std::endl;
    drawCube(cubeMasses, cubeSprings);
    // update time
    T = T + timeStep;
    // re-initilize energy variable
    kineticEnergy = 0;
    springPotentialEnergy = 0;
    gravityPotentialEnergy = 0;
    totalPotentialEnergy = 0;
    groundEnergy = 0;
    totalEnergy = 0;
}

void drawGrass()
{
    glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
    glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);
    //  Save transformation
    glPushMatrix();
    //  Offset
    //glTranslated(x,y,z);
    //glRotated(th,0,1,0);
    //glScaled(dx,dy,dz);
    // Enable Texture
    glEnable(GL_TEXTURE_2D);
    glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
    glColor3f(1,1,1);
    glBindTexture(GL_TEXTURE_2D,grassTexture);

    //  Front
    glColor3f(1,1,1);
    glBindTexture(GL_TEXTURE_2D,grassTexture);
    glBegin(GL_QUADS);
    glNormal3f( 0, 1, 0);
    glTexCoord2f(0.0,0.0);  glVertex3f(-1.0,+0.0,-1.0);
    glTexCoord2f(0.0,1.0);  glVertex3f(+1.0,+0.0,-1.0);
    glTexCoord2f(1.0,1.0);  glVertex3f(+1.0,+0.0,+1.0);
    glTexCoord2f(1.0,0.0);  glVertex3f(-1.0,+0.0,+1.0);
    glEnd();
    //  Undo transformations
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
}

static void skyBox(double scale)
{
    glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
    glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);
    //  Save transformation
    glPushMatrix();
    //  Offset
    glScaled(scale,scale,scale);
    // Enable Texture
    glEnable(GL_TEXTURE_2D);
    glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
    glColor3f(1,1,1);
    glBindTexture(GL_TEXTURE_2D,skyBoxTexture[0]);
    // Front
    glBindTexture(GL_TEXTURE_2D,skyBoxTexture[0]);
    glBegin(GL_QUADS);
    glNormal3f( 0, 0,-1);
    glTexCoord2f(0.0f,0.0f);    glVertex3f( +5,-3,+5);    // Lower Left
    glTexCoord2f(1.0f,0.0f);    glVertex3f( -5,-3,+5);    // Lower Right
    glTexCoord2f(1.0f,1.0f);    glVertex3f( -5,+7,+5);    // Upper Right
    glTexCoord2f(0.0f,1.0f);    glVertex3f( +5,+7,+5);    // Upper Left
    glEnd();
    // Back
    glBindTexture(GL_TEXTURE_2D,skyBoxTexture[1]);
    glBegin(GL_QUADS);
    glNormal3f( 0, 0,+1);
    glTexCoord2f(0.0f,0.0f);    glVertex3f( -5,-3,-5);    // Lower Left
    glTexCoord2f(1.0f,0.0f);    glVertex3f( +5,-3,-5);    // Lower Right
    glTexCoord2f(1.0f,1.0f);    glVertex3f( +5,+7,-5);    // Upper Right
    glTexCoord2f(0.0f,1.0f);    glVertex3f( -5,+7,-5);    // Upper Left
    glEnd();
    // Left
    glBindTexture(GL_TEXTURE_2D,skyBoxTexture[2]);
    glBegin(GL_QUADS);
    glNormal3f(+1, 0, 0);
    glTexCoord2f(0.0f,0.0f);    glVertex3f( -5, -3, +5);    // Lower Left
    glTexCoord2f(1.0f,0.0f);    glVertex3f( -5, -3, -5);    // Lower Right
    glTexCoord2f(1.0f,1.0f);    glVertex3f( -5, +7,-5);    // Upper Right
    glTexCoord2f(0.0f,1.0f);    glVertex3f( -5, +7,+5);    // Upper Left
    glEnd();
    // Right
    glBindTexture(GL_TEXTURE_2D,skyBoxTexture[3]);
    glColor3f(1,1,1);
    glBegin(GL_QUADS);
    glNormal3f(-1, 0, 0);
    glTexCoord2f(0.0f,0.0f);    glVertex3f( +5, -3, -5);    // Lower Left
    glTexCoord2f(1.0f,0.0f);    glVertex3f( +5, -3, +5);    // Lower Right
    glTexCoord2f(1.0f,1.0f);    glVertex3f( +5, +7,+5);    // Upper Right
    glTexCoord2f(0.0f,1.0f);    glVertex3f( +5, +7,-5);    // Upper Left
    glEnd();
    // Top
    glBindTexture(GL_TEXTURE_2D,skyBoxTexture[4]);
    glColor3f(1,1,1);
    glBegin(GL_QUADS);
    glNormal3f( 0,-1, 0);
    glTexCoord2f(0.0f,0.0f);    glVertex3f( +5,+7,-5);
    glTexCoord2f(1.0f,0.0f);    glVertex3f( +5,+7,+5);
    glTexCoord2f(1.0f,1.0f);    glVertex3f( -5,+7,+5);
    glTexCoord2f(0.0f,1.0f);    glVertex3f( -5,+7,-5);
    glEnd();
    // Bottom
    glBindTexture(GL_TEXTURE_2D,skyBoxTexture[5]);
    glColor3f(1,1,1);
    glBegin(GL_QUADS);
    glNormal3f( 0,+1, 0);
    glTexCoord2f(0.0f,0.0f);    glVertex3f( +5,+0,-5);
    glTexCoord2f(1.0f,0.0f);    glVertex3f( +5,+0,+5);
    glTexCoord2f(1.0f,1.0f);    glVertex3f( -5,+0,+5);
    glTexCoord2f(0.0f,1.0f);    glVertex3f( -5,+0,-5);
    glEnd();
    //  Undo transformations
    glPopMatrix();
}

static void ball(double x,double y,double z,double r)
{
    //  Save transformation
    glPushMatrix();
    //  Offset, scale and rotate
    glTranslated(x,y,z);
    glScaled(r,r,r);
    //  White ball
    glColor3f(1,1,1);
    glutSolidSphere(1.0,16,16);
    //  Undo transofrmations
    glPopMatrix();
}

void Print(const char* format , ...)
{
    char    buf[LEN];
    char*   ch=buf;
    va_list args;
    //  Turn the parameters into a character string
    va_start(args,format);
    vsnprintf(buf,LEN,format,args);
    va_end(args);
    //  Display the characters one at a time at the current raster position
    while (*ch)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,*ch++);
}

/*
 *  OpenGL (GLUT) calls this routine to display the scene
 */
void display()
{
    const double len=1.5;  //  Length of axes
    //  Erase the window and the depth buffer
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    //  Enable Z-buffering in OpenGL
    glEnable(GL_DEPTH_TEST);
    //  Undo previous transformations
    glLoadIdentity();
    //  Eye position
    double Ex = -2*dim*Sin(th)*Cos(ph);
    double Ey = +2*dim        *Sin(ph);
    double Ez = +2*dim*Cos(th)*Cos(ph);
    gluLookAt(Ex,Ey,Ez, 0,0,0, 0,Cos(ph),0);;
    // draw light
    if (light)
    {
        //  Translate intensity to color vectors
        float Ambient[]   = {0.01*ambient ,0.01*ambient ,0.01*ambient ,2.0};
        float Diffuse[]   = {0.01*diffuse ,0.01*diffuse ,0.01*diffuse ,2.0};
        float Specular[]  = {0.01*specular,0.01*specular,0.01*specular,2.0};
        //  Light direction
        //float Position[]  = {0,10,0};
        //float Position[]  = {5*Cos(zh),ylight,5*Sin(zh),1};
        float Position[]  = {0,0.5,0,1};
        //  Draw light position as ball (still no lighting here)
        glColor3f(16,16,16);
        ball(Position[0],Position[1],Position[2], 0.001);
        //  OpenGL should normalize normal vectors
        glEnable(GL_NORMALIZE);
        //  Enable lighting
        glEnable(GL_LIGHTING);
        //  glColor sets ambient and diffuse color materials
        glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
        glEnable(GL_COLOR_MATERIAL);
        //  Enable light 0
        glEnable(GL_LIGHT0);
        //  Set ambient, diffuse, specular components and position of light 0
        glLightfv(GL_LIGHT0,GL_AMBIENT ,Ambient);
        glLightfv(GL_LIGHT0,GL_DIFFUSE ,Diffuse);
        glLightfv(GL_LIGHT0,GL_SPECULAR,Specular);
        glLightfv(GL_LIGHT0,GL_POSITION,Position);
    }
    // draw scene
    cubeUpdate(cubeMasses, cubeSprings, 0);
    //drawGrass();
    skyBox(skyBoxScale);
    //  Draw axes
    glColor3f(1,1,1);
    if (axes)
    {
        glBegin(GL_LINES);
        glVertex3d(0.0,0.0,0.0);
        glVertex3d(len,0.0,0.0);
        glVertex3d(0.0,0.0,0.0);
        glVertex3d(0.0,len,0.0);
        glVertex3d(0.0,0.0,0.0);
        glVertex3d(0.0,0.0,len);
        glEnd();
        //  Label axes
        glRasterPos3d(len,0.0,0.0);
        Print("X");
        glRasterPos3d(0.0,len,0.0);
        Print("Y");
        glRasterPos3d(0.0,0.0,len);
        Print("Z");
    }
    //  Render the scene
    glFlush();
    //  Make the rendered scene visible
    glutSwapBuffers();
}

/*
 *  GLUT calls this routine when an arrow key is pressed
 */
void special(int key,int x,int y)
{
    //  Right arrow key - increase angle by 5 degrees
    if (key == GLUT_KEY_RIGHT)
        th += 5;
    //  Left arrow key - decrease angle by 5 degrees
    else if (key == GLUT_KEY_LEFT)
        th -= 5;
    //  Up arrow key - increase elevation by 5 degrees
    else if (key == GLUT_KEY_UP)
    {
        if (ph +5 < 90)
            {
                ph += 5;
            }
    }
    //  Down arrow key - decrease elevation by 5 degrees
    else if (key == GLUT_KEY_DOWN)
    {
        if (ph-5>0)
            {
                ph -= 5;
            }
    }
    //  Keep angles to +/-360 degrees
    th %= 360;
    ph %= 360;
    //  Tell GLUT it is necessary to redisplay the scene
    glutPostRedisplay();
}

/*
 *  Set projection
 */
void Project(double fov,double asp,double dim)
{
   //  Tell OpenGL we want to manipulate the projection matrix
   glMatrixMode(GL_PROJECTION);
   //  Undo previous transformations
   glLoadIdentity();
   //  Perspective transformation
   if (fov)
      gluPerspective(fov,asp,dim/16,16*dim);
   //  Orthogonal transformation
   else
      glOrtho(-asp*dim,asp*dim,-dim,+dim,-dim,+dim);
   //  Switch to manipulating the model matrix
   glMatrixMode(GL_MODELVIEW);
   //  Undo previous transformations
   glLoadIdentity();
}

/*
 *  GLUT calls this routine when a key is pressed
 */
void key(unsigned char ch,int x,int y)
{
    //  Exit on ESC
    if (ch == 27)
        exit(0);
    //  Reset view angle
    else if (ch == '0')
        th = ph = 0;
    //  Toggle axes
    else if (ch == 'a' || ch == 'A')
        axes = 1-axes;
    //  Change field of view angle
    else if (ch == '-' && ch>1)
        fov++;
    else if (ch == '=' && ch<179)
        fov--;
    //  PageUp key - increase dim
    else if (ch == GLUT_KEY_PAGE_DOWN){
        dim += 0.1;
    }
    //  PageDown key - decrease dim
    else if (ch == GLUT_KEY_PAGE_UP && dim>1){
        dim -= 0.1;
    }
    else if (ch == 'w'){
        cubeMasses[0].v[2] = 1;
    }
    //  Keep angles to +/-360 degrees
    th %= 360;
    ph %= 360;
    //  Reproject
    Project(fov,asp,dim);
    //  Tell GLUT it is necessary to redisplay the scene
    glutPostRedisplay();
}

/*
 *  GLUT calls this routine when the window is resized
 */
void reshape(int width,int height)
{
    //  Ratio of the width to the height of the window
    asp = (height>0) ? (double)width/height : 1;
    //  Set the viewport to the entire window
    glViewport(0,0, width,height);
    //  Set projection
    Project(fov,asp,dim);
}

/*
 *  GLUT calls this toutine when there is nothing else to do
 */
void idle()
{
    glutPostRedisplay();
}

int main(int argc,char* argv[])
{
    // Initialize GLUT and process user parameters
    glutInit(&argc, argv);
    // double buffered, true color 600*600
    glutInitWindowSize(1000,800);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
    // create the window
    glutCreateWindow("hz2547");
    //  Tell GLUT to call "idle" when there is nothing else to do
    glutIdleFunc(idle);
    //  Tell GLUT to call "display" when the scene should be drawn
    glutDisplayFunc(display);
    //  Tell GLUT to call "reshape" when the window is resized
    glutReshapeFunc(reshape);
    //  Tell GLUT to call "special" when an arrow key is pressed
    glutSpecialFunc(special);
    //  Tell GLUT to call "key" when a key is pressed
    glutKeyboardFunc(key);
    // load texture
    grassTexture = LoadTexBMP("./texture/grass.bmp");
    slimeTexture = LoadTexBMP("./texture/slime.bmp");
    // Skybox texture
    skyBoxTexture[0] = LoadTexBMP("./texture/skyboxes/front.bmp");
    skyBoxTexture[1] = LoadTexBMP("./texture/skyboxes/back.bmp");
    skyBoxTexture[2] = LoadTexBMP("./texture/skyboxes/left.bmp");
    skyBoxTexture[3] = LoadTexBMP("./texture/skyboxes/right.bmp");
    skyBoxTexture[4] = LoadTexBMP("./texture/skyboxes/up.bmp");
    skyBoxTexture[5] = LoadTexBMP("./texture/skyboxes/down.bmp");
    //  Pass control to GLUT so it can interact with the user
    glutMainLoop();
    return 0;
};
