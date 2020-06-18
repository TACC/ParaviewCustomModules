#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _WIN32
#define EXTERN  __declspec(dllexport)
#else
#define EXTERN  extern 
#endif

void PyInit_Sample(){}

struct Point3
{
  float x, y, z;
};

struct Point3 cross(struct Point3 a, struct Point3 b)
{
  struct Point3 r;
  r.x = a.y*b.z - a.z*b.y;
  r.y = a.z*b.x - a.x*b.z;
  r.z = a.x*b.y - a.y*b.x;
  return r;
}

float dot(struct Point3 a, struct Point3 b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
}
  
// Never freed; owned by caller

static float *pdf = NULL;
static struct Point3 *points = NULL;
static int *tets = NULL;
static int nTets;
static int nPoints;

struct Point3 *samples = NULL;
static int nSamples;

#define RAND (((float)rand())/RAND_MAX)

float *
ComputeCumulativeWeightedTetVolumes(int pdep)
{
  float *cum_volume = (float *)malloc(nTets*sizeof(float));

  for (int i = 0; i < nTets; i++)
  {
    int p = tets[5*i + 1];
    int q = tets[5*i + 2];
    int r = tets[5*i + 3];
    int s = tets[5*i + 4];

    struct Point3 pp = points[p];
    struct Point3 pq = points[q];
    struct Point3 pr = points[r];
    struct Point3 ps = points[s];

    struct Point3 a;
    a.x = pq.x - pp.x;
    a.y = pq.y - pp.y;
    a.z = pq.z - pp.z;

    struct Point3 b;
    b.x = pr.x - pp.x;
    b.y = pr.y - pp.y;
    b.z = pr.z - pp.z;
    
    struct Point3 c;
    c.x = ps.x - pp.x;
    c.y = ps.y - pp.y;
    c.z = ps.z - pp.z;

    float cval = (pdep) ? (pdf[p] + pdf[q] + pdf[r] + pdf[s]) / 4.0 : pdf[i];
    float v = cval * dot(a, cross(b, c));

    cum_volume[i] = (i == 0) ? v : v + cum_volume[i-1];
  }

  return cum_volume;
}

static void
subdivide(int start, int end, int n, float *cum_volume, int *sampleIndices)
{ 
  float pvol = cum_volume[end] - ((start == 0) ? 0 : cum_volume[start-1]);
  if (pvol == 0)
    return;

  float nExpected = (pvol / cum_volume[nTets-1]) * n;

  if (nExpected == 0)
    return;
  
  // If we expect less than one sample in this subset, then decide whether we have one.
  if (nExpected < 1)
  {
    if (RAND < nExpected)
    {
      // then we will put one sample somewhere in this set of cells.

      for (int i = start, done = 0; i < end && !done; i++)
      {
        float v = cum_volume[i] - ((i == 0) ? 0 : cum_volume[i-1]);
        float p = v / pvol;
        if (RAND < p)
        {
          sampleIndices[nSamples++] = i;
          done = 1;
        }
      }
    }
  }
  else if (start == end)
  {
    for (int i = 0; i < nExpected; i++)
    {
      sampleIndices[nSamples++] = start;
    }

    if (RAND < (nExpected - floor(nExpected)))
    {
      sampleIndices[nSamples++] = start;
    }
  }
  else if (start + 1 == end)
  {
    subdivide(start, start, n, cum_volume, sampleIndices);
    subdivide(end, end, n, cum_volume, sampleIndices);
  }
  else
  {
    int mid = (start + end) / 2;
    subdivide(start, mid-1, n, cum_volume, sampleIndices);
    subdivide(mid, end, n, cum_volume, sampleIndices);
  }
  return;
}

EXTERN void
SampleTetrahedra(int pdep, int n, int np, int nt, float *p, float *d, int *t)
{
  nPoints = np;
  nTets = nt;
  points = (struct Point3*)p;
  tets   = (int*)t;

  float min_pdf = d[0];
  float max_pdf = d[0];

  int npdf = pdep ? nPoints : nTets;

  for (int i = 1; i < npdf; i++)
  {
    if (min_pdf > d[i]) min_pdf = d[i];
    if (max_pdf < d[i]) max_pdf = d[i];
  }

  fprintf(stderr, "minmax done... %f %f\n", min_pdf, max_pdf);

  pdf = (float *)malloc(nTets*sizeof(float));
  if (min_pdf == max_pdf)
    for (int i = 0; i < npdf; i++)
      pdf[i] = 1;
  else
    for (int i = 0; i < npdf; i++)
      pdf[i] = (d[i] - min_pdf) / (max_pdf - min_pdf);

  fprintf(stderr, "normalization done\n");

  float *cum_volume = ComputeCumulativeWeightedTetVolumes(pdep);

  fprintf(stderr, "weighted volumes  done\n");

  int *sampleIndices = (int *)malloc(2 * n * sizeof(int));

  nSamples = 0;
  subdivide(0, nTets-1, n, cum_volume, sampleIndices);

  samples = (struct Point3*)malloc(nSamples*sizeof(struct Point3));

  for (int i = 0; i < nSamples; i++)
  {
    int *tet = tets + (sampleIndices[i]*5);

    struct Point3 v0 = points[tet[1]];
    struct Point3 v1 = points[tet[2]];
    struct Point3 v2 = points[tet[3]];
    struct Point3 v3 = points[tet[4]];

    double s = RAND;
    double t = RAND;
    double u = RAND;

    if (s+t > 1.0)
    {
      s = 1.0 - s;
      t = 1.0 - t;
    }

    if (t+u > 1.0)
    {
      double tmp = u;
      u = 1.0 - s - t;
      t = 1.0 - tmp;
    }
    else if (s+t+u > 1.0)
    {
      double tmp = u;
      u = s + t + u - 1.0;
      s = 1 - t - tmp;
    }

    double a = 1 - (s+t+u);

    float x = v0.x*a + v1.x*s + v2.x*t + v3.x*u;
    float y = v0.y*a + v1.y*s + v2.y*t + v3.y*u;
    float z = v0.z*a + v1.z*s + v2.z*t + v3.z*u;

    samples[i].x = x;
    samples[i].y = y;
    samples[i].z = z;
  }

  free((void*)pdf);
  free((void*)cum_volume);
  free((void*)sampleIndices);
}


EXTERN void
SampleRectilinear(int pdep, int n, int *dim, float *x_array, float *y_array, float *z_array, float *pdf)
{
  int nCells = (dim[0] - 1) * (dim[1] - 1) * (dim[2] - 1);

  int istep, jstep, kstep;

  if (pdep)
  {
    istep = 1;
    jstep = dim[0];
    kstep = dim[0] * dim[1];
  }
  else
  {
    istep = 1;
    jstep = (dim[0] - 1);
    kstep = (dim[0] - 1) * (dim[1] - 1);
  }

  float *cell_pdf = (float *)malloc(nCells*sizeof(float));

  float tot_cval = 0;
  float min_cval = 0;
  float max_cval = 0;
  int first = 1;

  float *p = cell_pdf;

  for (int i = 0; i < dim[0] - 1; i++)
  {
    float dx = x_array[i+1] - x_array[i];

    for (int j = 0; j < dim[1] - 1; j++)
    {
      float dy = y_array[j+1] - y_array[j];

      for (int k = 0; k < dim[2] - 1; k++)
      {
        float dz = z_array[k+1] - z_array[k];

        float cval;

        if (pdep)
        {
          int corner = i*istep + j*jstep + k*kstep;
 
          cval = pdf[corner + 0*istep + 0*jstep + 0*kstep] +
                 pdf[corner + 1*istep + 0*jstep + 0*kstep] +
                 pdf[corner + 0*istep + 1*jstep + 0*kstep] +
                 pdf[corner + 1*istep + 1*jstep + 0*kstep] +
                 pdf[corner + 0*istep + 0*jstep + 1*kstep] +
                 pdf[corner + 1*istep + 0*jstep + 1*kstep] +
                 pdf[corner + 0*istep + 1*jstep + 1*kstep] +
                 pdf[corner + 1*istep + 1*jstep + 1*kstep];
        }
        else
        {
          int cell_index = i*istep + j*jstep + k*kstep;
          cval = pdf[cell_index];
        }

        cval *= (dx * dy * dz);

        *p++ = cval;
        tot_cval += cval;

        if (first)
        {     
          first = 0;
          min_cval = max_cval = cval;
        }
        else
        {
          if (min_cval > cval) min_cval = cval;
          if (max_cval < cval) max_cval = cval;
        }
      }
    }
  }

  min_cval = max_cval = cell_pdf[0];
  for (int i = 1; i < nCells; i++)
  {
    
    if (min_cval > cell_pdf[i]) min_cval = cell_pdf[i];
    if (max_cval < cell_pdf[i]) max_cval = cell_pdf[i];
  }

  float tot_pdf = 0;
  if (min_cval == max_cval)
  {
    for (int i = 0; i < nCells; i++)
      cell_pdf[i] = 1;
    tot_pdf = nCells;
  }
  else
  {
    for (int i = 0; i < nCells; i++)
    {
      cell_pdf[i] = (cell_pdf[i] - min_cval) / (max_cval - min_cval);
      tot_pdf += cell_pdf[i];
    }
  }

  int *cell_samples = (int *)malloc(nCells*sizeof(int));

  nSamples = 0;
  for (int i = 0; i < nCells; i++)
  {
    float f = (cell_pdf[i] / tot_pdf) * n;
    cell_samples[i] = ((int)f) + ((RAND < (f - (int)f)) ? 1 : 0);
    nSamples += cell_samples[i];
  }

  free(cell_pdf);

  if (samples) free(samples);
  samples = (struct Point3*)malloc(nSamples*sizeof(struct Point3));

  int s = 0, c = 0;
  for (int i = 0; i < dim[0] - 1; i++)
  {
    float x = x_array[i];
    float dx = x_array[i+1] - x_array[i];

    for (int j = 0; j < dim[1] - 1; j++)
    {
      float y = y_array[j];
      float dy = y_array[j+1] - y_array[j];

      for (int k = 0; k < dim[2] - 1; k++)
      {
        float z = z_array[k];
        float dz = z_array[k+1] - z_array[k];

        for (int l = 0; l < cell_samples[c]; l++, s++)
        {
          samples[s].x = x + RAND*dx;
          samples[s].y = y + RAND*dy;
          samples[s].z = z + RAND*dz;
        }
      }
    }
  }

  free(cell_samples);
}
            
EXTERN void
SampleVTI(int pdep, int n, int *dim, float *orig, float *spacing, float *pdf)
{
  int nCells = (dim[0] - 1) * (dim[1] - 1) * (dim[2] - 1);

  int istep = 1;
  int jstep = dim[0];
  int kstep = dim[0] * dim[1];

  float *cell_pdf = (float *)malloc(nCells*sizeof(float));

  float tot_cval = 0;
  float min_cval = 0;
  float max_cval = 0;
  int first = 1;

  float *p = cell_pdf;

  int cell_id = 0;
  for (int i = 0; i < dim[0] - 1; i++)
  {
    for (int j = 0; j < dim[1] - 1; j++)
    {
      for (int k = 0; k < dim[2] - 1; k++, cell_id++)
      {
        float cval;

        if (pdep)
        {
          int corner = i*istep + j*jstep + k*kstep;
          cval = pdf[corner + 0*istep + 0*jstep + 0*kstep] +
                 pdf[corner + 1*istep + 0*jstep + 0*kstep] +
                 pdf[corner + 0*istep + 1*jstep + 0*kstep] +
                 pdf[corner + 1*istep + 1*jstep + 0*kstep] +
                 pdf[corner + 0*istep + 0*jstep + 1*kstep] +
                 pdf[corner + 1*istep + 0*jstep + 1*kstep] +
                 pdf[corner + 0*istep + 1*jstep + 1*kstep] +
                 pdf[corner + 1*istep + 1*jstep + 1*kstep];
        }
        else
          cval = pdf[cell_id];


        *p++ = cval;
        tot_cval += cval;

        if (first)
        {     
          first = 0;
          min_cval = max_cval = cval;
        }
        else
        {
          if (min_cval > cval) min_cval = cval;
          if (max_cval < cval) max_cval = cval;
        }
      }
    }
  }

  float tot_pdf = 0;
  if (min_cval == max_cval)
  {
    for (int i = 0; i < nCells; i++)
      cell_pdf[i] = 1;
    tot_pdf = nCells;
  }
  else
  {
    for (int i = 0; i < nCells; i++)
    {
      cell_pdf[i] = (cell_pdf[i] - min_cval) / (max_cval - min_cval);
      tot_pdf += cell_pdf[i];
    }
  }

  int *cell_samples = (int *)malloc(nCells*sizeof(int));

  nSamples = 0;
  for (int i = 0; i < nCells; i++)
  {
    float f = (cell_pdf[i] / tot_pdf) * n;
    cell_samples[i] = ((int)f) + ((RAND < (f - (int)f)) ? 1 : 0);
    nSamples += cell_samples[i];
  }

  free(cell_pdf);

  if (samples) free(samples);
  samples = (struct Point3*)malloc(nSamples*sizeof(struct Point3));

  int s = 0, c = 0;
  for (int i = 0; i < dim[0] - 1; i++)
    for (int j = 0; j < dim[1] - 1; j++)
      for (int k = 0; k < dim[2] - 1; k++, c++)
      {
        float x = orig[0] + i*spacing[0]; 
        float y = orig[1] + j*spacing[1]; 
        float z = orig[2] + k*spacing[2]; 

        for (int l = 0; l < cell_samples[c]; l++, s++)
        {
          samples[s].x = x + RAND*(spacing[0]);
          samples[s].y = y + RAND*(spacing[1]);
          samples[s].z = z + RAND*(spacing[2]);
        }
      }

  free(cell_samples);
}
            
EXTERN void *
GetSamples() { return (void *)samples; }

EXTERN int
GetNumberOfSamples() { return nSamples; }

EXTERN void
Cleanup()
{
  if (samples)
  {
    free((void *)samples);
    samples = NULL;
  }
}
