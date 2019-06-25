#include "main.h"

#define SPACE_DIM 3
#define TREE_LEVEL 4
#define W 4
#define NUM_THREADS 16
#define pMutation 0.9

#define GRAPHICS
#define UNIQUEPOSITION

// OPENGL Variables
int th = 0;            //  Azimuth of view angle
int ph = 0;           //  Elevation of view angle
int axes = 0;         //  Display axes
int light = 0;
double asp = 1;     // aspect ratio
int fov = 40;         //  Field of view (for perspective)
double dim = 30.0;  // size of the workd
double skyBoxScale = 1.0;
double cx=0;            //  Camera Location
double cy=5;
double cz=4;
double view=1000;
double viewlr=90;
int mode = 1;

int generationNumber = 1;
int robotNumber = 4;
int simulationTime = 10;

float emission  = 60;  // Emission intensity (%)
float ambient   = 60;  // Ambient intensity (%)
float diffuse   = 60;  // Diffuse intensity (%)
float specular  = 60;  // Specular intensity (%)
float shininess = 64;  // Shininess (power of two)
float shiny   =   1;    // Shininess (value)
float white[] = {1,1,1,1};
float black[] = {0,0,0,1};

unsigned int grassTexture;
unsigned int slimeTexture;
unsigned int skyBoxTexture[10]; // Texture for Sky Box

// Physics Simluator Variables
double mass = 0.5;
double length = 1.0;
double gravity = 9.8;
double T = 0;

double timeStep = 0.001;
double restoreConstant = 10000;
double springConstant = 5000;
double dampingConstant = 0.99;
double frictionCoefficient = 0.8;

static GLint Frames = 0;
static GLfloat fps = -1;
static GLint T0 = 0;

double start_time;

GLfloat worldRotation[16] = {1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1};


double norm( double x[], std::size_t sz )
{
    return std::sqrt( std::inner_product( x, x+sz, x, 0.0 ) ) ;
}

std::vector<int> sort_indexes(const std::vector<double> &v) {

    // initialize original index locations
    std::vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](int i1, int i2) {return v[i1] > v[i2];});
    return idx;
}

std::vector<int> generateChildrenIndex(int crossOverPoint){
    std::vector<int> childrenIndex;
    childrenIndex.push_back(crossOverPoint);
    std::queue<int> myQueue;
    myQueue.push(crossOverPoint);
    while (!myQueue.empty()){
        int index = myQueue.front();
        if (3*index+3 <= (pow(3, TREE_LEVEL)-1)/2) {
            int indexArray[3] = {3*index+1, 3*index+2, 3*index+3};
            for (int i = 0; i < 3; i++){
                childrenIndex.push_back(indexArray[i]);
                myQueue.push(indexArray[i]);
            }
        }
        myQueue.pop();
    }
    return childrenIndex;
}

std::vector<int> generateCrossOverRange(int crossOverLevel){
    std::vector<int> crossOverRange;
    int start = (pow(3, crossOverLevel+1)-1)/2;
    int end = (pow(3, crossOverLevel+2)-1)/2;
    for (int i = start; i < end; i++){
        crossOverRange.push_back(i);
    }
    return crossOverRange;
}

int randomNumber(int min, int max){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(min, max);
    return dis(gen);
}


int randomSpringFactor()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 3);
    return dis(gen);
}

class ROBOT
{
public:
    std::vector<GENE> robotGene;
    std::vector<MASS> robotMasses;
    std::vector<SPRING> robotSprings;
    double initialLocation[3] = {0,0,0};

    ROBOT(double initialX, double initialY, double initialZ, std::vector<GENE> robotGenes){
        robotGene = robotGenes;
        generateRobot(initialX, initialY, initialZ);
        generateSpring();
        calculateInitialLocation();
    }

    void calculateInitialLocation(){
        double x = 0; double y = 0;
        for (int j = 0; j < robotMasses.size(); j++){
            x = x + robotMasses[j].p[0];
            y = y + robotMasses[j].p[1];
        }
        x = x / robotMasses.size();
        y = y / robotMasses.size();
        initialLocation[0] = x; initialLocation[1] = y;
    }


    void generateSpring(){

        for (int i = 0; i < robotMasses.size() - 1; i++){
            for (int j = i+1; j < robotMasses.size(); j++){
                double positionDiff[3] = {robotMasses[j].p[0] - robotMasses[i].p[0],robotMasses[j].p[1] - robotMasses[i].p[1],robotMasses[j].p[2] - robotMasses[i].p[2]};
                if (norm(positionDiff,3) < 2*length){
                    robotSprings.push_back({springConstant,norm(positionDiff,3),norm(positionDiff,3),i,j,(robotMasses[i].springFactor+robotMasses[j].springFactor)/2});
                }
            }
        }
        /*
        std::queue<int> myQueue;
        myQueue.push(0);
        while (!myQueue.empty()){
            int index = myQueue.front();
            //std::cout << index << std::endl;
            //printf("X: %f, Y: %f\n", robotGene[index].p[0], robotGene[index].p[1]);
            if (3*index+3 <= (pow(3, TREE_LEVEL)-1)/2) {

                int indexArray[4] = {index, 3*index+1, 3*index+2, 3*index+3};
                for (int i = 0; i < 3; i++){
                    for (int j = i+1; j < 4; j++){
                        double positionDiff1[3] = {robotGene[indexArray[i]].p[0] - robotGene[indexArray[j]].p[0],
                                                   robotGene[indexArray[i]].p[1] - robotGene[indexArray[j]].p[1],
                                                   robotGene[indexArray[i]].p[2] - robotGene[indexArray[j]].p[2]};
                        double L1 = norm(positionDiff1, 3);
                        //std::cout << L1 << std::endl;
                        //std::cout << (robotGene[indexArray[i]].springFactor+robotGene[indexArray[j]].springFactor)/2 << std::endl;
                        robotSprings.push_back({springConstant, L1, L1, indexArray[i], indexArray[j], (robotGene[indexArray[i]].springFactor+robotGene[indexArray[j]].springFactor)/2});
                    }
                }

                for (int i = 1; i < 4; i++){
                    std::random_device rd;
                    std::mt19937 gen(rd());
                    std::uniform_int_distribution<> dis(1, 10);
                    if (dis(gen) > 5){
                        double positionDiff1[3] = {robotGene[indexArray[i]].p[0] - robotGene[0].p[0],
                                                   robotGene[indexArray[i]].p[1] - robotGene[0].p[1],
                                                   robotGene[indexArray[i]].p[2] - robotGene[0].p[2]};
                        double L1 = norm(positionDiff1, 3);
                        robotSprings.push_back({springConstant, L1, L1, indexArray[i], 0, (robotGene[indexArray[i]].springFactor+robotGene[0].springFactor)/2});
                        //std::cout <<  "Add Spring to Root" << std::endl;
                    }
                }

                int indexArray2[3] = {3*index+1, 3*index+2, 3*index+3};
                for (int i = 0; i < 3; i++){
                    if (robotGene[indexArray2[i]].exist == true){
                        myQueue.push(indexArray2[i]);
                    }
                }
            }
            myQueue.pop();
        }*/
    }

    void generateRobot(double iX, double iY, double iZ){
        std::queue<int> myQueue;
        myQueue.push(0);

        while (!myQueue.empty()){
            int index = myQueue.front();
            //std::cout << index << std::endl;
            //printf("X: %f, Y: %f\n", robotGene[index].p[0], robotGene[index].p[1]);
            robotMasses.push_back(
                    {0.5, {iX+robotGene[index].p[0], iY+robotGene[index].p[1], iZ+robotGene[index].p[2]}, {0, 0, 0}, {0, 0, 0}, robotGene[index].springFactor});

            if (3*index+3 <= (pow(3, TREE_LEVEL)-1)/2) {
                int indexArray2[3] = {3*index+1, 3*index+2, 3*index+3};
                for (int i = 0; i < 3; i++){
                    if (robotGene[indexArray2[i]].exist == true){
                        myQueue.push(indexArray2[i]);
                    }
                }
            }
            myQueue.pop();
        }
    }


    void robotDraw() {
        glColor3f(1, 0, 0);

        GLUquadric *quad;
        quad = gluNewQuadric();
        for (int i = 0; i < (int) robotMasses.size(); i++) {
            glPushMatrix();
            glMultMatrixf(worldRotation);
            glTranslated(robotMasses[i].p[0], robotMasses[i].p[1], robotMasses[i].p[2]);
            gluSphere(quad, length / 20, 10, 10);
            glPopMatrix();
        }

        for (int i = 0; i < (int) robotSprings.size(); i++) {
            if (robotSprings[i].type == 0) {glColor3f(0.0,0.4,0.8);}
            if (robotSprings[i].type == 1) {glColor3f(1.0,0.4,0.4);}
            if (robotSprings[i].type == 2) {glColor3f(0.4,0.6,0.0);}
            if (robotSprings[i].type == 3) {glColor3f(0.4,0.7,1.0);}
            glPushMatrix();
            glMultMatrixf(worldRotation);
            glBegin(GL_LINES);
            glLineWidth(1000);
            glVertex3f(GLfloat(robotMasses[robotSprings[i].m1].p[0]), GLfloat(robotMasses[robotSprings[i].m1].p[1]),
                       GLfloat(robotMasses[robotSprings[i].m1].p[2]));
            glVertex3f(GLfloat(robotMasses[robotSprings[i].m2].p[0]), GLfloat(robotMasses[robotSprings[i].m2].p[1]),
                       GLfloat(robotMasses[robotSprings[i].m2].p[2]));
            //printf("Drawing Line between {%f, %f, %f} and {%f, %f, %f}\n", robotMasses[robotSprings[i].m1].p[0], robotMasses[robotSprings[i].m1].p[1], robotMasses[robotSprings[i].m1].p[2], robotMasses[robotSprings[i].m2].p[0], robotMasses[robotSprings[i].m2].p[1], robotMasses[robotSprings[i].m2].p[2]);
            glEnd();
            glPopMatrix();
        }
        // draw line between middle point and initial position
        double x = 0; double y = 0;
        for (int j = 0; j < robotMasses.size(); j++){
            x = x + robotMasses[j].p[0];
            y = y + robotMasses[j].p[1];
        }
        x = x / robotMasses.size();
        y = y / robotMasses.size();
        glPushMatrix();
        glMultMatrixf(worldRotation);
        glBegin(GL_LINES);
        glColor3f(0.0, 1.0, 0.0);
        glVertex3f(GLfloat(x), GLfloat(y), GLfloat(0.0));
        glVertex3f(GLfloat(initialLocation[0]), GLfloat(initialLocation[1]), GLfloat(0.0));
        glEnd();
        glPopMatrix();

    }

    void robotUpdate()
    {
        // initialize the force vector with value 0
        std::vector<std::vector<double>> robotForces((int)robotMasses.size(),std::vector<double>(3));
        // loop through all springs to calculate spring forces
        for (int i = 0; i < (int)robotSprings.size(); i++) {
            if (robotSprings[i].type == 0){
                // hard
                robotSprings[i].k = 2000;
            }else if(robotSprings[i].type == 1){
                // cos
                robotSprings[i].k = 1000;
                robotSprings[i].L = robotSprings[i].L_0 - 0.3 * length * sin(W*T);
            }else if(robotSprings[i].type == 2){
                // sin
                robotSprings[i].k = 1000;
                robotSprings[i].L = robotSprings[i].L_0 + 0.3 * length * sin(W*T);
            }else if(robotSprings[i].type == 3){
                // soft
                robotSprings[i].k = 1000;
            }
            MASS mass1 = robotMasses[robotSprings[i].m1];
            MASS mass2 = robotMasses[robotSprings[i].m2];
            double positionDiff[3] = {mass2.p[0] - mass1.p[0], mass2.p[1] - mass1.p[1], mass2.p[2] - mass1.p[2]};
            double L = norm(positionDiff, 3);
            double force = robotSprings[i].k * fabs(robotSprings[i].L - L);
            double direction[3];
            if (L == 0){
                direction[0] = 0;direction[1] = 0;direction[2] = 0;
            }else{
                direction[0] = positionDiff[0] / L;direction[1] = positionDiff[1] / L;direction[2] = positionDiff[2] / L;
            }
            // contraction case
            if (L > robotSprings[i].L) {
                robotForces[robotSprings[i].m1][0] = robotForces[robotSprings[i].m1][0] + direction[0] * force;
                robotForces[robotSprings[i].m1][1] = robotForces[robotSprings[i].m1][1] + direction[1] * force;
                robotForces[robotSprings[i].m1][2] = robotForces[robotSprings[i].m1][2] + direction[2] * force;
                robotForces[robotSprings[i].m2][0] = robotForces[robotSprings[i].m2][0] - direction[0] * force;
                robotForces[robotSprings[i].m2][1] = robotForces[robotSprings[i].m2][1] - direction[1] * force;
                robotForces[robotSprings[i].m2][2] = robotForces[robotSprings[i].m2][2] - direction[2] * force;
            }
                // expansion case
            else if (L < robotSprings[i].L) {
                robotForces[robotSprings[i].m1][0] = robotForces[robotSprings[i].m1][0] - direction[0] * force;
                robotForces[robotSprings[i].m1][1] = robotForces[robotSprings[i].m1][1] - direction[1] * force;
                robotForces[robotSprings[i].m1][2] = robotForces[robotSprings[i].m1][2] - direction[2] * force;
                robotForces[robotSprings[i].m2][0] = robotForces[robotSprings[i].m2][0] + direction[0] * force;
                robotForces[robotSprings[i].m2][1] = robotForces[robotSprings[i].m2][1] + direction[1] * force;
                robotForces[robotSprings[i].m2][2] = robotForces[robotSprings[i].m2][2] + direction[2] * force;
            }
        }

        for (int i = 0; i < (int)robotMasses.size(); i++) {
            // add gravity
            robotForces[i][2] = robotForces[i][2] - robotMasses[i].m * gravity;
            // if the mass is below ground, add restroration force and calculate friction
            if (robotMasses[i].p[2] <= 0) {
                robotForces[i][2] = robotForces[i][2] + restoreConstant * fabs(robotMasses[i].p[2]);
                // calculate horizontal force and vertical force
                double F_h = sqrt(pow(robotForces[i][0], 2) + pow(robotForces[i][1], 2));
                double F_v = robotForces[i][2];
                if (F_h < F_v * frictionCoefficient) {
                    robotForces[i][0] = 0;
                    robotForces[i][1] = 0;
                    robotMasses[i].v[0] = 0;
                    robotMasses[i].v[1] = 0;
                }
            }
            // update acceleration
            robotMasses[i].a[0] = robotForces[i][0] / robotMasses[i].m;
            robotMasses[i].a[1] = robotForces[i][1] / robotMasses[i].m;
            robotMasses[i].a[2] = robotForces[i][2] / robotMasses[i].m;

            // update velocity
            robotMasses[i].v[0] = dampingConstant * (robotMasses[i].v[0] + robotMasses[i].a[0] * timeStep);
            robotMasses[i].v[1] = dampingConstant * (robotMasses[i].v[1] + robotMasses[i].a[1] * timeStep);
            robotMasses[i].v[2] = dampingConstant * (robotMasses[i].v[2] + robotMasses[i].a[2] * timeStep);

            // update position
            robotMasses[i].p[0] = robotMasses[i].p[0] + robotMasses[i].v[0] * timeStep;
            robotMasses[i].p[1] = robotMasses[i].p[1] + robotMasses[i].v[1] * timeStep;
            robotMasses[i].p[2] = robotMasses[i].p[2] + robotMasses[i].v[2] * timeStep;
        }
    }
};

class Simulation{
private:
    int populationSize;
    std::vector<double> populationDistance;
    std::vector<std::vector<GENE>> populationGene;
    std::vector<std::vector<GENE>> newPopulationGene;
    std::vector<ROBOT> robots;
    std::vector<std::vector<double>> massPositions;
    std::ofstream bestGene;
    std::ofstream popDis;
public:
    double averageDistance;
    double maxDistance;

    Simulation(int popSize) {
        populationSize = popSize;
        generateGenes();
        //generateBestGene();
        generateRobots();
        std::time_t t = std::time(0);
        std::tm* now = std::localtime(&t);
        std::string hour =  std::to_string(now->tm_hour);
        std::string min =  std::to_string(now->tm_min);
        std::string popDisFileName = "populationDistance_" + hour + "_" + min + ".txt";
        std::cout << popDisFileName << std::endl;
        std::string bestGeneFileName = "bestGene_" + hour + "_" + min + ".txt";
        std::cout << bestGeneFileName << std::endl;
        popDis.open (popDisFileName); popDis.close();
        bestGene.open(bestGeneFileName); bestGene.close();
        popDis.open (popDisFileName);
        bestGene.open(bestGeneFileName);
    }
    void startSim(double time) {
        if (T < time) {
            simUpdate();
        #ifdef  GRAPHICS
            simDraw();
        #endif
        } else {
            printf("### Generation %d ###", generationNumber);
            double time = omp_get_wtime() - start_time;
            printf(" Time: %f ###\n", time);
            calculatePopulationDistance();
            selection();
            crossOver();
            populationGene.clear();
            populationGene.shrink_to_fit();
            populationGene = newPopulationGene;
            generationNumber++;
            robots.clear();
            robots.shrink_to_fit();
            generateRobots();
            T = 0;
            start_time = omp_get_wtime();
        }
    }

    void generateBestGene(){
        int gene98[52] = {0,1,1,0,1,1,1,1,0,2,1,0,2,2,2,0,3,1,0,2,1,1,0,0,0,1,2,0,1,1,0,0,3,0,1,2,3,1,1,2,2,1,0,2,1,0,0,2,2,0,0,1};
        int gene190[52] = {0,1,1,0,1,1,1,1,1,1,0,0,1,1,0,0,3,1,0,2,1,2,2,2,3,1,0,2,3,1,0,2,1,1,0,0,1,1,0,0,3,1,0,2,1,1,0,0,3,1,0,2};
        int gene299[52] = {0,1,1,0,1,1,1,1,1,1,1,1,1,1,2,2,3,1,0,2,1,1,0,0,1,1,0,0,3,1,0,2,1,1,0,0,3,1,0,2,3,1,2,1,1,1,0,0,1,1,0,0};
        int gene417[52] = {0,1,1,0,1,1,1,1,1,1,2,2,1,1,2,2,3,1,0,2,1,1,0,0,3,1,0,2,3,1,2,1,1,1,1,1,1,1,0,1,3,1,2,1,1,1,0,0,1,1,0,0};
        int gene500[52] = {0,1,1,0,1,1,0,1,1,1,2,2,1,1,2,2,3,1,0,2,1,1,0,0,1,1,2,2,1,1,2,2,1,1,2,2,1,1,0,0,3,1,2,1,1,1,2,2,1,1,1,1};
        int gene1708[160] = {1,2,2,1,2,3,2,0,2,3,2,0,3,0,0,3,3,0,3,0,0,2,0,1,0,3,1,1,1,1,3,1,0,2,0,1,1,2,0,0,3,0,3,0,2,0,3,1,3,2,1,0,3,1,3,1,0,2,1,2,2,1,2,2,3,1,3,1,0,1,2,1,2,1,0,0,3,0,3,0,0,3,1,1,2,1,0,0,0,3,1,1,0,2,0,1,0,2,0,1,3,0,3,0,0,1,3,2,3,0,3,0,2,1,0,1,3,3,1,1,2,1,3,0,3,1,3,1,0,1,3,0,2,1,2,2,1,1,1,1,0,0,3,0,0,1,2,1,1,3,2,1,0,3,3,3,0,0,0,1};
        int gene1600[484] = {0,4,4,1,2,4,3,3,3,4,0,1,2,2,3,3,3,4,0,1,3,1,0,1,1,1,4,3,3,0,2,4,3,2,0,1,3,3,4,2,3,0,1,4,1,1,4,3,1,2,1,1,0,3,1,0,2,3,4,2,2,4,2,0,3,3,4,2,3,4,2,1,2,0,3,3,2,1,3,2,2,0,2,0,0,4,0,0,0,1,0,3,2,3,2,0,2,1,0,2,2,3,1,2,3,4,0,1,1,4,1,4,0,1,0,3,2,3,3,0,2,1,1,4,0,2,0,1,2,3,2,0,2,2,0,2,2,4,1,1,0,4,0,0,0,0,1,3,1,2,2,2,0,1,1,1,0,2,1,3,0,0,4,3,0,0,3,4,3,1,0,3,2,3,3,0,1,0,3,2,2,3,1,0,1,1,0,1,1,3,3,4,0,2,0,2,3,1,4,3,1,3,3,4,3,0,4,2,1,0,0,4,1,1,4,4,1,2,4,2,2,3,3,0,1,4,2,1,2,1,1,2,2,0,3,1,2,3,1,0,2,4,2,1,2,3,3,3,1,3,2,1,0,0,0,4,2,2,3,2,0,1,3,4,2,0,1,1,0,0,0,0,0,4,3,4,3,3,3,0,2,4,1,2,1,4,3,0,2,2,2,3,2,4,3,0,1,3,2,1,0,1,2,3,2,3,3,0,1,4,2,1,2,1,0,2,1,0,4,3,1,0,0,4,1,2,3,3,2,0,3,1,3,3,3,0,2,3,1,0,0,0,4,3,0,0,2,1,3,3,0,1,2,3,3,0,2,0,4,4,3,4,3,3,0,4,0,1,1,1,1,4,0,3,4,3,0,4,2,2,0,0,1,3,0,1,3,2,2,3,3,0,1,0,0,4,2,3,1,1,2,4,2,2,1,2,3,4,0,1,3,0,2,3,3,0,1,4,1,1,2,0,3,4,1,4,1,3,2,1,2,4,3,3,0,0,2,2,3,2,0,1,3,4,2,1,0,2,2,1,0,1,0,3,3,2,0,1,1,1,1,1,4,2,2,0,1,1,3,4,2,2,2,0,2,0,0,4,2,4,0,1,2,3};
        std::vector<GENE> temp;
        for (int i = 0; i < 484; i = i+4){
            std::vector<double> p{(double)gene1600[i+1], (double)gene1600[i+2], (double)gene1600[i+3]};
            temp.push_back({true, gene1600[i], p});
        }
        populationGene.push_back(temp);
    }

    void selection(){
        std::vector<int> index = sort_indexes(populationDistance);
        newPopulationGene.clear();
        newPopulationGene.shrink_to_fit();
        for (int i = 0; i < index.size()/2; i++) {
            newPopulationGene.push_back(populationGene[index[i]]);
        }
        for (int i = 0; i < newPopulationGene[0].size(); i++){
            bestGene << newPopulationGene[0][i].springFactor << " " << newPopulationGene[0][i].p[0] << " " << newPopulationGene[0][i].p[1] << " " << newPopulationGene[0][i].p[2] << " ";
        }
        bestGene << "\n";
        bestGene.flush();
    }

    void crossOver(){
        for (int n = 0; n < populationSize/4; n++){
            int parentIndex1 = randomNumber(0, newPopulationGene.size()-1);
            int parentIndex2 = randomNumber(0, newPopulationGene.size()-1);
            //printf("New population size: %d \n", newPopulationGene.size());
            //printf("Parent Index: %d %d \n", parentIndex1, parentIndex2);
            std::vector<GENE> parent1 = newPopulationGene[parentIndex1];
            std::vector<GENE> parent2 = newPopulationGene[parentIndex2];
            int maxLevel = std::max(log(parent1.size()*2-1)/log(3), log(parent2.size()*2-1)/log(3));
            int crossOverLevel = randomNumber(1, maxLevel);
            //printf("Crossover Level: %d \n", crossOverLevel);
            std::vector<int> crossOverRange = generateCrossOverRange(crossOverLevel-1);
            int crossOverPoint1 = crossOverRange[randomNumber(0, crossOverRange.size()-1)];
            int crossOverPoint2 = crossOverRange[randomNumber(0, crossOverRange.size()-1)];
            //printf("Crossover Point: %d, %d\n", crossOverPoint1, crossOverPoint2);
            std::vector<int> crossOverLocations1 = generateChildrenIndex(crossOverPoint1);
            std::vector<int> crossOverLocations2 = generateChildrenIndex(crossOverPoint2);
            std::vector<GENE> subHeap1; std::vector<GENE> subHeap2;
            //printf("Subheap size: %d, %d\n", subHeap1.size(), subHeap2.size());
            for (int i = 0; i < crossOverLocations1.size(); i++){
                //printf("Cross Point: %d\n", crossOverLocations1[i]);
                //printf("Parent 1 Size: %d\n", parent1.size());
                subHeap1.push_back(parent1[crossOverLocations1[i]]);
            }
            for (int i = 0; i < crossOverLocations2.size(); i++){
                subHeap2.push_back(parent2[crossOverLocations2[i]]);
            }
            //printf("Subheap size: %d, %d\n", subHeap1.size(), subHeap2.size());
            std::vector<GENE> offSpring1(parent1);
            std::vector<GENE> offSpring2(parent2);
            for (int i = 0; i < crossOverLocations1.size(); i++){
                offSpring1[crossOverLocations1[i]] = subHeap2[i];
                offSpring2[crossOverLocations2[i]] = subHeap1[i];
            }
            offSpring1 = mutation(offSpring1);
            offSpring2 = mutation(offSpring2);
            newPopulationGene.push_back(offSpring1);
            newPopulationGene.push_back(offSpring2);
        }
    }

    std::vector<GENE> mutation(std::vector<GENE> offSpring){
        for (int i = 0; i < offSpring.size(); i++) {
            float r = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/1.0));
            if (r > 0.9){
                offSpring[i].p[0] = randomNumber(0, SPACE_DIM);
                offSpring[i].p[1] = randomNumber(0, SPACE_DIM);
                offSpring[i].p[2] = randomNumber(0, SPACE_DIM);
            }
        }
        return offSpring;
    }

    std::vector<double> generatePosition(int MAX){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, MAX);
        std::vector<double> position;
        #ifndef UNIQUEPOSITION
        for (int i = 0; i < 3; i++)
        {
            int randN = dis(gen);
            position.push_back(randN);
        }
        #endif
        #ifdef UNIQUEPOSITION
        while (1){
            position.clear();
            for (int i = 0; i < 3; i++)
            {
                int randN = dis(gen);
                position.push_back(randN);
            }
            if (std::find(massPositions.begin(), massPositions.end(), position) == massPositions.end()){
                break;
            }
        }
        massPositions.push_back(position);
        #endif
        return position;
    }

    void generateGenes(){
        for (int i = 0; i < populationSize; i++) {
            massPositions.clear();
            std::vector<GENE> robotGene;
            for (int i = 0; i < (pow(3, TREE_LEVEL) - 1) / 2; i++) {
                int level = (int) (log(i * 2 - 1) / log(3));
                bool existence = false;
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_int_distribution<> dis(1, 20);
                if (level < dis(gen)) {
                    existence = true;
                }
                robotGene.push_back({existence, randomSpringFactor(), generatePosition(SPACE_DIM)});
            }
            populationGene.push_back(robotGene);
        }
    }

    void generateRobots(){
        for (int i = 0; i < populationSize; i++){
            double X = -20 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/40));
            double Y = -20 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/40));
            robots.push_back(ROBOT(X, Y, 0.1, populationGene[i]));
        }
    }

    void simUpdate(){
        #ifndef GRAPHICS
        #pragma omp parallel for num_threads(NUM_THREADS)
        #endif
        for (int i = 0; i < populationSize; i++){
            robots[i].robotUpdate();
        }
    }

    void simDraw(){
        for (int i = 0; i < populationSize; i++){
            robots[i].robotDraw();
        }
    }

    void calculatePopulationDistance(){
        populationDistance.clear();
        populationDistance.shrink_to_fit();
        for (int i = 0; i < populationSize; i++){
            double x = 0; double y = 0;
            for (int j = 0; j < robots[i].robotMasses.size(); j++){
                //printf("X: %f, Y: %f \n", robots[i].robotMasses[j].p[0], robots[i].robotMasses[j].p[1]);
                x = x + robots[i].robotMasses[j].p[0];
                y = y + robots[i].robotMasses[j].p[1];
            }
            x = x / robots[i].robotMasses.size();
            y = y / robots[i].robotMasses.size();
            double distance[2] = {fabs(x - robots[i].initialLocation[0]), fabs(y - robots[i].initialLocation[1])};
            double distanceNorm = norm(distance, 2);
            //printf("Distance %f\n", distanceNorm);
            populationDistance.push_back(distanceNorm);
        }
        averageDistance = 0;
        maxDistance = 0;
        int counter = 0;
        for (int i = 0; i < populationSize; i++){
            averageDistance = averageDistance + populationDistance[i];
            maxDistance = std::max(maxDistance, populationDistance[i]);
            popDis << populationDistance[i] << " ";
        }
        popDis << "\n"; popDis.flush();
        averageDistance = averageDistance/populationSize;
        std::cout << "Maximum Distance: " << maxDistance << std::endl;
        std::cout << "Average Distance: " << averageDistance << std::endl;
    }

};

#ifdef GRAPHICS
Simulation sim1(robotNumber);
#endif

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
        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24,*ch++);
}

void drawGrid(){
    for (int i = -dim/2; i < dim/2 + 1; i++) {
        for (int j = -dim / 2; j < dim / 2 + 1; j++) {
            float white[] = {1,1,1,1};
            float black[] = {0,0,0,1};
            glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
            glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
            glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,black);

            glPushMatrix();
            glTranslatef(i*2, 0, j*2);
            glBegin(GL_QUADS);
            //
            glNormal3f(0, 1, 0);
            glColor3f(0.1, 0.1, 0.1);
            glVertex3f(+0, -0.01, +0);
            glVertex3f(+1, -0.01, +0);
            glVertex3f(+1, -0.01, +1);
            glVertex3f(+0, -0.01, +1);
            //
            glNormal3f(0, 1, 0);
            glColor3f(0.15, 0.15, 0.15);
            glVertex3f(-1, -0.01, +0);
            glVertex3f(+0, -0.01, +0);
            glVertex3f(+0, -0.01, +1);
            glVertex3f(-1, -0.01, +1);
            //
            glNormal3f(0, 1, 0);
            glColor3f(0.1, 0.1, 0.1);
            glVertex3f(-1, -0.01, -1);
            glVertex3f(+0, -0.01, -1);
            glVertex3f(+0, -0.01, +0);
            glVertex3f(-1, -0.01, +0);
            //
            glNormal3f(0, 1, 0);
            glColor3f(0.15, 0.15, 0.15);
            glVertex3f(-0, -0.01, -1);
            glVertex3f(+1, -0.01, -1);
            glVertex3f(+1, -0.01, +0);
            glVertex3f(-0, -0.01, +0);
            glEnd();
            glPopMatrix();
        }
    }

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

void display()
{
    const double len=2.0;  //  Length of axes
    //  Erase the window and the depth buffer
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    //  Enable Z-buffering in OpenGL
    glEnable(GL_DEPTH_TEST);
    //  Undo previous transformations
    glLoadIdentity();
    //  Eye position
    if (mode == 1) {
        //  Eye position
        double Ex = -2*dim*Sin(th)*Cos(ph);
        double Ey = +2*dim        *Sin(ph);
        double Ez = +2*dim*Cos(th)*Cos(ph);
        gluLookAt(Ex,Ey,Ez, 0,0,0, 0,Cos(ph),0);
    }
    if (mode == 2) {
        gluLookAt(cx,cy,cz,cx+view*Cos(viewlr),cy,cz+view*Sin(viewlr),0,1,0);
    }

    if (light)
    {
        //  Translate intensity to color vectors
        float Ambient[]   = {(float)0.01*ambient,(float)0.01*ambient,(float)0.01*ambient ,2.0};
        float Diffuse[]   = {(float)0.01*diffuse,(float)0.01*diffuse,(float)0.01*diffuse ,2.0};
        float Specular[]  = {(float)0.01*specular,(float)0.01*specular,(float)0.01*specular,2.0};
        //  Light direction
        //float Position[]  = {0,10,0};
        //float Position[]  = {5*Cos(zh),ylight,5*Sin(zh),1};
        float Position[]  = {95,90,60,1};
        //  Draw light position as ball (still no lighting here)
        glColor3f(16,16,16);
        ball(Position[0],Position[1],Position[2], 5);
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
    else {
        glDisable(GL_LIGHTING);
    }

    drawGrid();



#ifdef GRAPHICS
    glColor3f(1,1,1);
    glWindowPos2i(0,0);
    Print("Generation: %d", generationNumber);
    glWindowPos2i(850, 0);
    Print("Max: %4.2f", sim1.maxDistance);
    glWindowPos2i(850, 50);
    Print("Average %4.2f", sim1.averageDistance);
    sim1.startSim(simulationTime);
    T = T + timeStep;
#endif



    Frames++;
    GLint t = glutGet(GLUT_ELAPSED_TIME);
    if (t - T0 >= 1000) {
        GLfloat seconds = (t - T0) / 1000.0;
        fps = Frames / seconds;
        //printf("%d frames in %6.3f seconds = %6.3f FPS\n", Frames, seconds, fps);
        T0 = t;
        Frames = 0;
    }

    glRasterPos3d(0.0,2,0.0);
    //if (fps>0) Print("FPS %.3f", fps);
    //drawGrass();
    //skyBox(skyBoxScale);
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
        //axes = 1-axes;
        int x;
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
    else if (ch == 'w' || ch == 'W'){
        cx += 0.3*Cos(viewlr);
        cz += 0.3*Sin(viewlr);
    }
    else if (ch == 'a' || ch == 'A'){
        cx += 0.3*Sin(viewlr);
        cz -= 0.3*Cos(viewlr);
    }
    else if (ch == 's' || ch == 'S'){
        cx -= 0.3*Cos(viewlr);
        cz -= 0.3*Sin(viewlr);
    }
    else if (ch == 'd' || ch == 'D'){
        cx -= 0.3*Sin(viewlr);
        cz += 0.3*Cos(viewlr);
    }
    else if (ch == 'r' || ch == 'R'){
        if (cy+0.3 < skyBoxScale*10){
            cy += 0.3;
        }
    }
    else if (ch == 'f' || ch == 'F'){
        if (cy-0.3 > 0){
            cy -= 0.3;
        }
    }
    else if (ch == 'q' || ch=='Q'){
        viewlr-=15;
    }
    else if (ch == 'e' || ch=='E'){
        viewlr+=15;
    }

    if (ch == '1')
    {
        mode = 1;
    }
    else if (ch == '2')
    {
        mode = 2;
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

int main(int argc,char* argv[]) {
    #ifdef GRAPHICS
    // Initialize GLUT and process user parameters
    glutInit(&argc, argv);
    glWindowPos2i =  (PFNGLWINDOWPOS2IPROC) glutGetProcAddress("glWindowPos2i");
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
    //  Pass control to GLUT so it can interact with the user
    glutMainLoop();
    #endif

    #ifndef GRAPHICS
    Simulation sim1(robotNumber);
    while (1) {
        double start_time = omp_get_wtime();
        sim1.startSim(simulationTime);
        T = T + timeStep;
        //printf("%f\n", T);
        double time = omp_get_wtime() - start_time;
    }
    #endif

    return 0;
}
