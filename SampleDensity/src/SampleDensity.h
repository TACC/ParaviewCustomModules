#include <stdio.h>
#include <stdlib.h>

struct Point3
{
  Point3() : x(0), y(0), z(0) {}
  Point3(float x, float y, float z) : x(x), y(y), z(z) {}
  Point3 operator-(const Point3 p) { return Point3(this->x - p.x, this->y - p.y, this->z - p.z); }
  float dot(const Point3 p) { return (x * p.x) + (y * p.y) + (z * p.z); }
  Point3 cross(const Point3 p) { return Point3(y*p.z - z*p.y, z*p.x - x*p.z, x*p.y - y*p.x); }
  float x, y, z;
};
  
// Never freed; owned by caller

static float *pdf = NULL;
static Point3 *points = NULL;
static int *tets = NULL;

// always freed

static float *cum_volume = NULL;
static Point3 *samples = NULL;
static int *sampleIndices = NULL;

static int nTets;
static int nPoints;
static int nSamples;
static int allocSamples = 0;


