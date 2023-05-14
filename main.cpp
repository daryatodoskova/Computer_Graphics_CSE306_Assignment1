#define _CRT_SECURE_NO_WARNINGS 1
#include <iostream>
#include <vector>
#include <chrono>
#include <string>
#include <cstring>
#include <stdio.h>
#include <iostream> 
#include <algorithm>
#include <list>
#include <random>
#include <random>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "/usr/local/opt/libomp/include/omp.h"

static std::default_random_engine engine(10); // random seed = 10
static std::uniform_real_distribution<double> uniform(0., 1.);

//VECTOR 

class Vector {
  public:
    explicit Vector(double x = 0., double y = 0., double z = 0.);
    Vector& operator+=(const Vector& b);
    Vector &operator-=(const Vector& b);
    Vector& operator/=(double t);
    const double& operator[](int i) const;
    double& operator[](int i);

  private:
    double coords[3];
};

Vector::Vector(double x, double y, double z) {
  coords[0] = x;
  coords[1] = y;
  coords[2] = z;
}

Vector& Vector::operator+=(const Vector& b) {
  coords[0] += b[0];
  coords[1] += b[1];
  coords[2] += b[2];
  return *this;
}

Vector& Vector::operator-=(const Vector& b) {
  coords[0] -= b[0];
  coords[1] -= b[1];
  coords[2] -= b[2];
  return *this;
}

Vector& Vector::operator/=(double t) {
  coords[0] /= t;
  coords[1] /= t;
  coords[2] /= t;
  return *this;
}

const double& Vector::operator[](int i) const { return coords[i]; }
double& Vector::operator[](int i) { return coords[i]; }


Vector operator+(const Vector& a, const Vector& b) {
  return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector& a, const Vector& b) {
  return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator*(const Vector& a, double t) {
  return Vector(a[0]*t, a[1]*t, a[2]*t);
}

Vector operator*(double t, const Vector& a) {
  return Vector(a[0]*t, a[1]*t, a[2]*t);
}

Vector operator*(const Vector &a, const Vector &b) {
  return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}

Vector operator/(const Vector& a, double t) {
  return Vector(a[0]/t, a[1]/t, a[2]/t);
}

double dot(const Vector& a, const Vector& b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

Vector cross(const Vector& a, const Vector& b) {
  return Vector(a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0]);
}

double norm(const Vector& a) {
  return sqrt(dot(a, a));
}

Vector normalize(const Vector& a) {
  return a/sqrt(dot(a, a));
}

// INTERSECTION ----------------------------------------------------------------

struct Intersection {
  bool is_inter = false;
  bool mirror = false;
  bool is_light = false;
  double refractive_index = 1.;
  double t;
  Vector sphere_center;
  Vector P;
  Vector N;
  Vector albedo;
};

// BOUNDING BOX ----------------------------------------------------------------

struct BoundingBox {
  Vector m, M; // the two extremas of our bounding box
};

// RAY -------------------------------------------------------------------------

class Ray {
    public:
        Ray(const Vector& O, const Vector& u) : O(O), u(u) {};
        Vector O; //origin
        Vector u; //unit direction
};

// GEOMETRY --------------------------------------------------------------------

class Geometry {
  public:
    virtual Intersection intersect(const Ray &ray) = 0;
    bool mirror;
};

// SPHERE ----------------------------------------------------------------------

class Sphere : public Geometry {
  public:
    Sphere(
      const Vector& C, //center
      const double& R, //radius
      const Vector& albedo, //color
      bool mirror = false,
      double refractive_index = 1.,
      bool invert_normal = false,
      bool is_light_source = false
    ) {
        this->C = C;
        this->R = R;
        this->albedo = albedo;
        this->mirror = mirror;
        this->refractive_index = refractive_index;
        this->invert_normal = invert_normal;
        this->is_light_source = is_light_source;
    }
    //RAY-SPHERE INTERSECTION
    //similar to Scene's intersection function
    Intersection intersect(const Ray &ray) override {
      Intersection intersection;
      //reduced discriminant formula (p.20 of pdf)
      double delta = pow(dot(ray.u, (ray.O - C)), 2) - (pow(norm((ray.O - C)), 2) - pow(R, 2)); 
      intersection.is_inter = delta >= 0;
      intersection.t = 0; 
      if (intersection.is_inter) {
          // two possible intersection parameters
          double t_1 = dot(ray.u, (ray.O - C)*(-1.)) - sqrt(delta);
          double t_2 = dot(ray.u, (ray.O - C)*(-1.)) + sqrt(delta);
          if (t_2 < 0) { // the ray does not intersect the sphere
            intersection.is_inter = false;
          }
          else {
            if (t_1 >= 0) {
              intersection.t = t_1;
            } 
            else {
              intersection.t = t_2;
            }
          }
      }
      intersection.P = ray.O + ray.u*intersection.t; //location of intersection point P
      //for further lighting computation
      intersection.N = normalize(intersection.P - C);  //unit normal N at P
      intersection.albedo = albedo;
      intersection.sphere_center = this->C;
      intersection.refractive_index = refractive_index;
      if (this->mirror) {
        intersection.mirror = true;
      }
      if (this->invert_normal) {
        intersection.N = (-1.)*intersection.N;
      }
      if (this->is_light_source){
        intersection.is_light = true;
      }
      return intersection;
    }
    Vector albedo;
  private:
    Vector C;
    double R;
    double refractive_index;
    bool invert_normal;
    bool is_light_source;
};

//BOX MULLER
//lecture notes p.31
void boxMuller(double stdev, double &x, double &y) {
  double r1 = uniform(engine);
  double r2 = uniform(engine);
  // double r1 = ((double) rand() / (RAND_MAX));
  // double r2 = ((double) rand() / (RAND_MAX));
  x = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2) * stdev;
  y = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2) * stdev;
}

Vector random_cos(const Vector &N) {
    // int current_thread = omp_get_thread_num();  //add include to use this  //change in vsc settings too to make faster
    // double r1 = uniform(engine);
    // double r2 = uniform(engine);

    double r1 = ((double) rand() / (RAND_MAX));
    double r2 = ((double) rand() / (RAND_MAX));

    //lecture pdf p.32
    double x = sqrt(1-r1) * cos(2. * M_PI * r2);
    double y = sqrt(1-r1) * sin(2. * M_PI * r2);
    double z = sqrt(r2);

    Vector T1;
    if (std::abs(N[0]) <= std::abs(N[1]) && std::abs(N[0]) <= std::abs(N[2])) {
        T1 = Vector(0,N[2], -N[1]);
    }
    else {
        if ((std::abs(N[1]) <= std::abs(N[0])) && (std::abs(N[1]) <= std::abs(N[2]))) {
            T1 = Vector(N[2], 0, -N[0]);
        }
        else{
            T1 = Vector(N[1], -N[0], 0);
        }
    }
    T1 = normalize(T1);
    Vector T2 = cross(T1, N);
    return x*T1+y*cross(N,T1)+z*N;
}

// MESH ------------------------------------------------------------------------

//from https://pastebin.com/CAgp9r15

class TriangleIndices { // CLASS or STRUCT?
public:
	TriangleIndices(
	  	int vtxi = -1, int vtxj = -1, int vtxk = -1,
	  	int ni = -1, int nj = -1, int nk = -1,
	  	int uvi = -1, int uvj = -1, int uvk = -1,
	  	int group = -1,
      bool added = false
	) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk),
		uvi(uvi), uvj(uvj), uvk(uvk),
		ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

struct Node {
  Node *child_left;
  Node *child_right;
  BoundingBox bounding_box;
  int starting_triangle;
  int ending_triangle;
};

class TriangleMesh : public Geometry {
  public:
    ~TriangleMesh() {}
    TriangleMesh(
      double scaling_factor,
      const Vector &translation,
      const Vector &albedo,
      bool mirror = false
    ) {
      this->mirror = mirror;
      this->scaling_factor = scaling_factor;
      this->translation = translation;
      this->albedo = albedo;
      this->bvh_root = new Node;
      // double scaling_factor;
        // Vector translation;
        // Vector albedo;
        // Node* bvh_root;
        // void build_bvh(Node *node, int starting_triangle, int ending_triangle);
        // bool intersect_bounding_box(const Ray &ray, BBox bounding_box, double &t);
        // std::vector<TriangleIndices> indices;
        // std::vector<Vector> vertices;
        // std::vector<Vector> normals;
        // std::vector<Vector> uvs;
        // std::vector<Vector> vertexcolors;
    };

    void readOBJ(const char* obj) {

      char matfile[255];
      char grp[255];

      FILE* f;
      f = fopen(obj, "r");
      int curGroup = -1;
      while (!feof(f)) {
        char line[255];
        if (!fgets(line, 255, f)) break;

        std::string linetrim(line);
        linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
        strcpy(line, linetrim.c_str());

        if (line[0] == 'u' && line[1] == 's') {
          sscanf(line, "usemtl %[^\n]\n", grp);
          curGroup++;
        }

        if (line[0] == 'v' && line[1] == ' ') {
          Vector vec;

          Vector col;
          if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
            col[0] = std::min(1., std::max(0., col[0]));
            col[1] = std::min(1., std::max(0., col[1]));
            col[2] = std::min(1., std::max(0., col[2]));

            vertices.push_back(vec);
            vertexcolors.push_back(col);

          } else {
            sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
            vertices.push_back(vec);
          }
        }
        if (line[0] == 'v' && line[1] == 'n') {
          Vector vec;
          sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
          normals.push_back(vec);
        }
        if (line[0] == 'v' && line[1] == 't') {
          Vector vec;
          sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
          uvs.push_back(vec);
        }
        if (line[0] == 'f') {
          TriangleIndices t;
          int i0, i1, i2, i3;
          int j0, j1, j2, j3;
          int k0, k1, k2, k3;
          int nn;
          t.group = curGroup;

          char* consumedline = line + 1;
          int offset;

          nn = sscanf(consumedline, "%i/%i/%i %i/%i/%i %i/%i/%i%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
          if (nn == 9) {
            if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
            if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
            if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
            if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
            if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
            if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
            if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
            if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
            if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
            indices.push_back(t);
          } else {
            nn = sscanf(consumedline, "%i/%i %i/%i %i/%i%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
            if (nn == 6) {
              if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
              if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
              if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
              if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
              if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
              if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
              indices.push_back(t);
            } else {
              nn = sscanf(consumedline, "%i %i %i%n", &i0, &i1, &i2, &offset);
              if (nn == 3) {
                if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
                if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
                if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
                indices.push_back(t);
              } else {
                nn = sscanf(consumedline, "%i//%i %i//%i %i//%i%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
                if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
                if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
                if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
                if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
                if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
                indices.push_back(t);
              }
            }
          }

          consumedline = consumedline + offset;

          while (true) {
            if (consumedline[0] == '\n') break;
            if (consumedline[0] == '\0') break;
            nn = sscanf(consumedline, "%i/%i/%i%n", &i3, &j3, &k3, &offset);
            TriangleIndices t2;
            t2.group = curGroup;
            if (nn == 3) {
              if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
              if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
              if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
              if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
              if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
              if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
              if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
              if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
              if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
              indices.push_back(t2);
              consumedline = consumedline + offset;
              i2 = i3;
              j2 = j3;
              k2 = k3;
            } else {
              nn = sscanf(consumedline, "%i/%i%n", &i3, &j3, &offset);
              if (nn == 2) {
                if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
                if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
                if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
                if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
                if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
                if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
                consumedline = consumedline + offset;
                i2 = i3;
                j2 = j3;
                indices.push_back(t2);
              } else {
                nn = sscanf(consumedline, "%i//%i%n", &i3, &k3, &offset);
                if (nn == 2) {
                  if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
                  if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
                  if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
                  if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
                  if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
                  if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
                  consumedline = consumedline + offset;
                  i2 = i3;
                  k2 = k3;
                  indices.push_back(t2);
                } else {
                  nn = sscanf(consumedline, "%i%n", &i3, &offset);
                  if (nn == 1) {
                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
                    consumedline = consumedline + offset;
                    i2 = i3;
                    indices.push_back(t2);
                  } else {
                    consumedline = consumedline + 1;
                  }
                }
              }
            }
          }

        }

      }
      fclose(f);
      build_bvh(this->bvh_root, 0, indices.size());
    }

    Intersection intersect(const Ray &ray) override {
      Intersection intersection;
      intersection.is_inter = false;
      double t;
      double min_t = MAXFLOAT;

      if (not intersect_bounding_box(ray, this->bvh_root->bounding_box, t)){
        return intersection;
      }
      
      // NO BVH
      // Vector A, B, C, N, e1, e2;
      // for (int i = 0; i < indices.size(); i++) {
      //   TriangleIndices triangle = this->indices[i];
      //   A = scaling_factor*vertices[triangle.vtxi] + translation;
      //   B = scaling_factor*vertices[triangle.vtxj] + translation;
      //   C = scaling_factor*vertices[triangle.vtxk] + translation;
      //   e1 = B - A;
      //   e2 = C - A;
      //   N = cross(e1, e2);

      //   double beta  =  dot(e2, cross(A - ray.O, ray.u)) / dot(ray.u, N);
      //   double gamma = -dot(e1, cross(A - ray.O, ray.u)) / dot(ray.u, N);
      //   double alpha = 1. - beta - gamma;

      //   if (alpha > 0. && beta > 0. && gamma > 0.) {
      //     double t = dot(A - ray.O, N) / dot(ray.u, N);
      //     if (0 < t && t < min_t) {
      //       min_t = t;
      //       intersection.is_inter = true;
      //       intersection.t = t;
      //       intersection.P = A + beta*e1 + gamma*e2;
      //       intersection.N = N;
      //       if (this->mirror) intersection.mirror = true;
      //       intersection.albedo = this->albedo;
      //     }
      //   }
      // }
      // return intersection;


      // LAB 4
      // BVH
      // pdf p.49
      std::list<Node*> nodes_to_visit;
      nodes_to_visit.push_front(this->bvh_root);
      while (not nodes_to_visit.empty()) {
        Node* curNode = nodes_to_visit.back();
        nodes_to_visit.pop_back();
        //if there is one child, then it is not a leaf, so test the bounding box
        if (curNode->child_left) {
          if (intersect_bounding_box(ray, curNode->child_left->bounding_box, t)) {
            if (t < min_t) nodes_to_visit.push_back(curNode->child_left);
          }
          if (intersect_bounding_box(ray, curNode->child_right->bounding_box, t)) {
            if (t < min_t) nodes_to_visit.push_back(curNode->child_right);
          }
        } else {
          // test all triangles between curNode->starting triangle 
          // and curNode->ending triangle as before
          // if an intersection is found, update best inter distance if needed
          Vector A, B, C, N, e1, e2;  //Normals
          for (int i = curNode->starting_triangle; i < curNode->ending_triangle; i++) {
            TriangleIndices triangle = this->indices[i];
            A = scaling_factor*vertices[triangle.vtxi] + translation;
            B = scaling_factor*vertices[triangle.vtxj] + translation;
            C = scaling_factor*vertices[triangle.vtxk] + translation;
            e1 = B - A;
            e2 = C - A;
            N = cross(e1, e2);

            double beta  =  dot(e2, cross(A - ray.O, ray.u)) / dot(ray.u, N);
            double gamma = -dot(e1, cross(A - ray.O, ray.u)) / dot(ray.u, N);
            double alpha = 1. - beta - gamma;

            if (alpha > 0. && beta > 0. && gamma > 0.) {
              double t = dot(A - ray.O, N) / dot(ray.u, N);
              if (0 < t && t < min_t) {
                min_t = t;
                intersection.is_inter = true;
                intersection.t = t;
                intersection.P = A + beta*e1 + gamma*e2;
                intersection.N = N;
                if (this->mirror) intersection.mirror = true;
                intersection.albedo = this->albedo;
              }
            }
          }
        }
      }
      return intersection;
    }

  private:

    void build_bvh(Node *node, int starting_triangle, int ending_triangle) {
      //BBox from starting triangle included to ending triangle excluded
      node->bounding_box =  compute_bounding_box(starting_triangle, ending_triangle);
      node->starting_triangle = starting_triangle;
      node->ending_triangle = ending_triangle;

      Vector diag = node->bounding_box.M - node->bounding_box.m;
      Vector middle_diag = node->bounding_box.m + diag*0.5;

      int longest_axis = 0;
      double max = abs(diag[0]);
      for (int i = 1; i < 3; i++) {
        if (abs(diag[i]) > max) {
          max = abs(diag[i]);
          longest_axis = i;
        }
      }

      int pivot_index = starting_triangle;
      for (int i = starting_triangle; i < ending_triangle; i++) {
        Vector vertex1 = scaling_factor*this->vertices[this->indices[i].vtxi] + translation;
        Vector vertex2 = scaling_factor*this->vertices[this->indices[i].vtxj] + translation;
        Vector vertex3 = scaling_factor*this->vertices[this->indices[i].vtxk] + translation;
        Vector barycenter = (vertex1 + vertex2 + vertex3)/3.;
        //the swap below guarantees triangles whose barycenter are smaller than middle diag are before ”pivot index”
        if (barycenter[longest_axis] < middle_diag[longest_axis]) {
          std::swap(indices[i], indices[pivot_index]);
          pivot_index++;
        }
      }
      // stopping criterion
      if (pivot_index <= starting_triangle || pivot_index >= ending_triangle - 5 || ending_triangle - starting_triangle < 5) {
        return;
      }
      node->child_left = new Node();
      node->child_right = new Node();
      build_bvh(node->child_left, starting_triangle, pivot_index);
      build_bvh(node->child_right, pivot_index, ending_triangle);
    }

    BoundingBox compute_bounding_box(int starting_triangle, int ending_triangle) {
      double min_x = MAXFLOAT, min_y = MAXFLOAT, min_z = MAXFLOAT;
      double max_x = - MAXFLOAT, max_y = - MAXFLOAT, max_z = - MAXFLOAT;
      for (int i = starting_triangle; i < ending_triangle; i++) {
        auto original_triangle_vertices = {
          this->vertices[this->indices[i].vtxi],
          this->vertices[this->indices[i].vtxj],
          this->vertices[this->indices[i].vtxk]
        };
        for (auto const& vertex: original_triangle_vertices) {
          Vector V = scaling_factor*vertex + translation;
          if (V[0] < min_x) {
            min_x = V[0];
          }
          else if (V[0] > max_x) {
            max_x = V[0];
          }
          if (V[1] < min_y) {
            min_y = V[1];
          }
          else if (V[1] > max_y) {
            max_y = V[1];
          }
          if (V[2] < min_z) {
            min_z = V[2];
          }
          else if (V[2] > max_z) {
            max_z = V[2];
          }
        }
      }
      BoundingBox bounding_box;
      bounding_box.m = Vector(min_x, min_y, min_z);
      bounding_box.M = Vector(max_x, max_y, max_z);
      return bounding_box;
    }

    bool intersect_bounding_box(const Ray &ray, BoundingBox bounding_box, double &t) {

      double tx1 = (bounding_box.m[0] - ray.O[0]) / ray.u[0];
      double tx2 = (bounding_box.M[0] - ray.O[0]) / ray.u[0];
      double txmin = std::min(tx1, tx2);
      double txmax = std::max(tx1, tx2);

      double ty1 = (bounding_box.m[1] - ray.O[1]) / ray.u[1];
      double ty2 = (bounding_box.M[1] - ray.O[1]) / ray.u[1];
      double tymin = std::min(ty1, ty2);
      double tymax = std::max(ty1, ty2);

      double tz1 = (bounding_box.m[2] - ray.O[2]) / ray.u[2];
      double tz2 = (bounding_box.M[2] - ray.O[2]) / ray.u[2];
      double tzmin = std::min(tz1, tz2);
      double tzmax = std::max(tz1, tz2);

      double actualT = std::max(txmin, std::max(tymin,tzmin));
      double furthestT = std::min(txmax, std::min(tymax,tzmax));

      if (furthestT > actualT > 0) {  
        t = actualT;
        return true;
      }
      return false;
    }
    
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    Vector albedo;
    Vector translation;
    double scaling_factor;
    BoundingBox box;
    Node* bvh_root;
};

// SCENE  ----------------------------------------------------------------------

class Scene {
  public:
    explicit Scene(
      const Vector& light_center,
      double light_radius,
      double light_intensity
    ){
      //Sphere S(Vector(0,0,0), 10., Vector (0., 0.5, 1.));
        // Sphere* left_wall = new Sphere(Vector(-1000, 0, 0), 940, Vector (0.5, 0.8, 0.1));
        // //Sphere left_wall(Vector(-1000, 0, 0), 940, Vector (0.5, 0.8, 0.1));
        // Sphere* right_wall = new Sphere(Vector(1000, 0, 0), 940, Vector (0.9, 0.2, 0.3));
        // Sphere* ceiling = new Sphere(Vector(0, 1000, 0), 940, Vector (0.3, 0.5, 0.3));
        // Sphere* floor = new Sphere(Vector(0, -1000, 0), 990, Vector (0.6, 0.5, 0.7));
        // Sphere* front_wall = new Sphere(Vector(0, 0, -1000), 940, Vector (0.1, 0.6, 0.7));
        // Sphere* behind_wall = new Sphere(Vector(0, 0, 1000), 940, Vector (0.8, 0.2, 0.9));
        // geometries.push_back(left_wall);
        // geometries.push_back(right_wall);
        // geometries.push_back(ceiling);
        // geometries.push_back(floor);
        // geometries.push_back(front_wall);
        // geometries.push_back(behind_wall);
        // this->light_center = light_center;
        // this->light_radius = light_radius;
        // this->light_intensity = light_intensity;

      Sphere* left_wall = new Sphere(Vector(0, 1000, 0), 940, Vector(1, 0, 0));
      Sphere* right_wall = new Sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1));
      Sphere* ceiling = new Sphere(Vector(1000, 0, 0), 940, Vector(1, 1, 0));
      Sphere* floor = new Sphere(Vector(-1000, 0, 0), 940, Vector(0, 1, 1));
      Sphere* front_wall = new Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0));
      Sphere* behind_wall = new Sphere(Vector(0, 0, 1000), 940, Vector(1, 0, 1));
      geometries.push_back(left_wall);
      geometries.push_back(right_wall);
      geometries.push_back(ceiling);
      geometries.push_back(floor);
      geometries.push_back(front_wall);
      geometries.push_back(behind_wall);
      this->light_center = light_center;
      this->light_radius = light_radius;
      this->light_intensity = light_intensity;
    }

    void addGeometry(Geometry* geometry) {
      geometries.push_back(geometry);
    }
    //RAY-SCENE INTERSECTION
    //computes the point of intersection between a Ray and the sphere, if any
    //returns an Intersection structure that contains all the relevant information including a bool flag
    Intersection intersect(const Ray& ray) {
      Intersection best_inters;
      best_inters.is_inter = false;
      double min_t = MAXFLOAT;
      for (auto &geometry : geometries) {
        Intersection intersection = geometry->intersect(ray);
        if (intersection.is_inter && intersection.t < min_t) {
          min_t = intersection.t;
          best_inters = intersection;
        }
      }
      //ray-sphere intersection that is closest to the ray origin among all 
      return best_inters;
    }

    Vector getColor(const Ray& ray, int ray_depth, bool last_bounce_diffuse = false) {
      if (ray_depth < 0) return Vector(0., 0., 0.);  // //terminates recursion at some point
      Intersection intersection = intersect(ray);
      if (not intersection.is_inter)
      {
        return Vector(0, 0, 0);
      }
      Vector Lo(0., 0., 0.);
      if (intersection.is_inter) {
        double eps = 1e-10;
        Vector N = intersection.N;
        Vector P = intersection.P + N*eps;

        if (intersection.is_light) {
          double R = light_radius;
          if (last_bounce_diffuse) {  //if this is an indirect diffuse bounce
            //if we hit a light source by chance via an indirect diffuse bounce, return 0 to avoid counting it twice
            return Vector(0., 0. ,0.);
          }
          else {
            return Vector(1., 1. ,1.)*light_intensity/(4*M_PI*M_PI*R*R);
          }
        }

        if (intersection.mirror) {
          // Reflection
          Ray reflection_ray = Ray(P, ray.u - (2*dot(N,ray.u)) * N);
          return getColor(reflection_ray, ray_depth - 1);
        }
        //diffuse surfaces
        else if (intersection.refractive_index != 1.) {
          // Refraction
          double n1, n2; //refraction indices
          if (dot(N,ray.u) > 0) {
            N = (-1.)*N;
            n1 = intersection.refractive_index;
            n2 = 1.;
          } else {
            n1 = 1.;
            n2 = intersection.refractive_index;
          }
          //FRESNEL LAW
          double k0 = pow(n1 - n2, 2.)/pow(n1 + n2, 2.); //reflection coefficient at normal incidence
          P = intersection.P - N*eps;
          if (1. - pow((n1/n2), 2.) * (1 - pow(dot(N,ray.u), 2.)) > 0) {
            // Standard Refraction
            Vector w_T = (n1/n2) * (ray.u - dot(N,ray.u)*N);
            Vector w_N = (-1.)*N * sqrt(1 - pow((n1/n2),2.)*(1 - pow(dot(N,ray.u),2.)));
            Vector w = w_T + w_N;
            double uni = uniform(engine); // (uniform) random number u between 0 and 1
            double R = k0 + (1-k0)*pow(1 - abs(dot(N,w)),5.); //reflection coefficient for incidence w
            if (uni < R) { 
              // launch a reflection ray
              Ray reflection_ray = Ray(P, ray.u - (2*dot(intersection.N,ray.u)) * intersection.N);
              return getColor(reflection_ray, ray_depth - 1);
            } else {
              //launch a refraction ray
              Ray refraction_ray = Ray(P, w);
              return getColor(refraction_ray, ray_depth - 1);
            }
          } else {
            // Total Internal Reflection
            Ray internal_reflection_ray = Ray(P, ray.u - (2*dot(intersection.N,ray.u)) * intersection.N);
            return getColor(internal_reflection_ray, ray_depth - 1);
          }
        }
        else {
          double I = light_intensity;
          double R = light_radius;
          Vector x = intersection.sphere_center;
          double squared_distance_light = R*R;
          //Shading and shadows computation
          //Direct Lighting
          Vector xprime = R*random_cos(normalize(x - light_center)) + light_center;
          Vector Nprime = normalize(xprime - light_center);
          double d = norm(xprime - P);
          Vector omega = normalize(xprime - P);
          Ray lightRay = Ray(P, omega);
          Intersection lightIntersection = intersect(lightRay);
          //computes the visibility term by launching a ray towards the light source
          double visibility = 1; 
          if (lightIntersection.is_inter && lightIntersection.t <= d){
            double visibility = 0;
          }
          Vector rho = intersection.albedo;
          //formula p.20 of pdf
          Lo = I/(4* M_PI* squared_distance_light)* rho/M_PI * visibility*std::max(dot(N,omega),0.)*std::max(dot(Nprime,(-1.)*omega),0.)/(pow(norm(xprime-P),2.)*(dot(Nprime, normalize(x-light_center))/(M_PI*R*R)));
          // Indirect Lighting
          Ray randomRay = Ray(P, random_cos(N)); //randomly sample ray using random_cos
          Lo += rho * getColor(randomRay, ray_depth - 1, true);
        }
      }
      return Lo;
    }

  private:
    std::vector<Geometry*> geometries;
    Vector light_center;
    double light_radius;
    double light_intensity;
};

//void addSphere(const Sphere& s){objects.push_back(s);}

    // bool intersect(const Ray& r, Vector& P, Vector& N, double& t, int& sphere_id){
    //     bool has_inter = false;
    //     t = std::numeric_limits<double>::max();

    //     for (int i=0; i<objects.size();i++) {

    //         Vector localP, localN;
    //         double localt;

    //         if (objects[i].intersect(r, localP, localN, localt)) {
    //             if (localt<t) {
    //                 t = localt;
    //                 P = localP;
    //                 N = localN;
    //                 has_inter = true;
    //                 sphere_id = i;
    //             }
    //         }
    //     }
    //     return has_inter;

    // }
    // private:
    //     std::vector<Geometry*> geometries;
    //     Vector light_center;
    //     double light_radius;
    //     double light_intensity;



// MAIN ------------------------------------------------------------------------

int main() {

  auto start = std::chrono::high_resolution_clock::now();

  //VARIABLES

    Vector light_center(10, 20, 5);   //good
    double light_radius = 3; // > 0   //good
    double light_intensity = 1e5;  //good

    Scene scene = Scene(light_center, light_radius, light_intensity);

  // Light Source
  // Sphere* light = new Sphere(light_center, light_radius, Vector(1., 0., 0.), false, 1., false, true);
  // scene.addGeometry(light);

  // Cat
  TriangleMesh* cat = new TriangleMesh(0.6, Vector(-10, -10, 0), Vector(0.3, 0.2, 0.25), false);
  cat->readOBJ("cat.obj");
  scene.addGeometry(cat);

  //Dog
  TriangleMesh* dog = new TriangleMesh(4, Vector(15, -10, 0), Vector(0.3, 0.2, 0.25), false);
  dog->readOBJ("dog2.obj");
  scene.addGeometry(dog);

  // Soccer ball
  // TriangleMesh* ball = new TriangleMesh(15, Vector(30, 5, 0), Vector(0.1, 0.1, 0.1), false);
  // ball->readOBJ("soccer_ball.obj");
  // scene.addGeometry(ball);

  // White Sphere
  Sphere* white_sphere = new Sphere(Vector(-5, -7, 5), 3, Vector(1., 1., 1.));
  scene.addGeometry(white_sphere);

  // Big Reflective Sphere
  Sphere* reflective_sphere = new Sphere(Vector(8, 20, -15), 13, Vector(1., 1., 1.), true);
  scene.addGeometry(reflective_sphere);

  // Solid Refractive Sphere
  Sphere* solid_refractive_sphere = new Sphere(Vector(-15, -7, 15), 3, Vector(1., 1., 1.), false, 1.5);
  scene.addGeometry(solid_refractive_sphere);

  // Hollow Refractive Sphere
  Sphere* hollow_refractive_sphere_outer = new Sphere(Vector(7, 0, 20), 4, Vector(1., 1., 1.), false, 1.5);
  Sphere* hollow_refractive_sphere_inner = new Sphere(Vector(7, 0, 20), 3.75, Vector(1., 1., 1.), false, 1.5, true);
  scene.addGeometry(hollow_refractive_sphere_outer);
  scene.addGeometry(hollow_refractive_sphere_inner);
    
    //512x512 images
    int W = 512; 
    int H = 512;

    //double I = 2*10000000000;
    double alpha = M_PI / 3.;

    int max_path_length = 5;  //4?  //good
    int n_paths = 64; //change to 64
    double fov = 60. * M_PI / 180; // field of view 60 deg 
    double gamma = 2.2;  //for gamma correction //color^(1/2.2)
    Vector camera_center = Vector(0, 0, 55);
    std::vector<unsigned char> image(W*H * 3, 0);  //good

    //PARALLELIZATION
  #pragma omp parallel for schedule(dynamic,1)
  //computing ray direction
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector pixelColor;
            for (int k = 0; k < n_paths; k++) {
                Vector rand_dir;
                rand_dir[0] = (camera_center[0] + j + 0.5 - W / 2);
                rand_dir[1] = (camera_center[1] + (H-i-1) + 0.5 - H / 2);
                rand_dir[2] = camera_center[2] - W / (2 * tan(fov / 2));
                if (random){
                  double x, y;
                  boxMuller(0.5, x, y);
                  rand_dir[0] += x;
                  rand_dir[1] += y;
                }

                Ray ray = Ray(camera_center, normalize(rand_dir - camera_center)); //cast a ray from the camera center to pixel i,j
                pixelColor += scene.getColor(ray, max_path_length);
        }

        pixelColor /= n_paths;

        image[(i*W + j) * 3 + 0] = std::max(std::min(255., pow(pixelColor[0], 1./gamma)*255),0.); //red
        image[(i*W + j) * 3 + 1] = std::max(std::min(255., pow(pixelColor[1], 1./gamma)*255),0.); //green
        image[(i*W + j) * 3 + 2] = std::max(std::min(255., pow(pixelColor[2], 1./gamma)*255),0.); //blue
        }
    }

    //save image and return 0

    stbi_write_png("final.png", W, H, 3, &image[0], 0);

    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  
    std::cout << "Function duration: " << duration.count() << " ms" << std::endl;
    return 0;
}