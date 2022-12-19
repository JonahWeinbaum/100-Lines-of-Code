#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <float.h>
#include <string.h>
#define BAR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
typedef enum { LAMBERTIAN, METAL, DIELECTRIC } SurfaceType_t;
static constexpr float epsilon = 0.01f;
typedef struct Vector { float x, y, z; } Vector;
Vector operator+(const Vector& v1, const Vector& v2) { return Vector({v1.x + v2.x, v1.y + v2.y, v1.z + v2.z}); }
Vector operator-(const Vector& v1, const Vector& v2) { return Vector({v1.x - v2.x, v1.y - v2.y, v1.z - v2.z}); }
Vector operator*(const Vector& v1, const Vector& v2) { return Vector({v1.x * v2.x, v1.y * v2.y, v1.z * v2.z}); }
Vector operator*(const float &c, const Vector& v) { return Vector({c * v.x, c * v.y, c * v.z}); }
float dot(const Vector& l, const Vector& r) { return l.x * r.x + l.y * r.y + l.z * r.z; }
Vector normalize(const Vector& v) { return (1 / std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z)) * v; }
typedef struct Ray { Vector o, d; float maxt; } Ray;
Vector rand_sph(Vector &d) {
    Vector p;
    do { p = 2.0f * Vector({rand() / (float)RAND_MAX, rand() / (float)RAND_MAX, rand() / (float)RAND_MAX}) - Vector({1, 1, 1}); } 
    while (dot(p, p) >= 1.0f || dot(p, d) < 0.f);
    return p;
}
typedef struct Surface {
    Surface(float _r, Vector _c, Vector _atten, float _emit, SurfaceType_t _type) : r(_r), c(_c), atten(_atten), emit(_emit), type(_type) {}
    float intersect(const Ray& wi) {
        float disc = dot(wi.o - c, wi.d)*dot(wi.o - c, wi.d) - dot(wi.d, wi.d)*(dot(wi.o - c, wi.o - c) - r*r); 
        float t1 = (-dot(wi.o - c, wi.d) - std::sqrt(disc)) * (1 / dot(wi.d, wi.d));
        float t2 = (-dot(wi.o - c, wi.d) + std::sqrt(disc)) * (1 / dot(wi.d, wi.d));
        return disc < 0.f ? FLT_MAX : ((t1 < epsilon || wi.maxt < t1) ? ((t2 < epsilon || wi.maxt < t2) ? FLT_MAX : t2) :  t1);
    }
    Vector scatter(const Ray& wi, Vector hit_point) {
        Vector n = dot(wi.d, hit_point - c) < 0.f ? normalize(hit_point - c) : -1.f * normalize(hit_point - c); 
        switch (type) {
            case LAMBERTIAN : {
                return n + normalize(rand_sph(n));
            } case METAL : {
                return normalize((normalize(wi.d) - 2 * dot(normalize(wi.d), n) * n) + 0.1*normalize(rand_sph(n))); 
            } case DIELECTRIC : {         
                float refr_ratio = dot(wi.d, hit_point - c) < 0.f ? (1.0f/1.5f) : 1.5f;
                float discrim = 1.0f - refr_ratio * refr_ratio * (1.0f - dot(normalize(wi.d), n) * dot(normalize(wi.d), n));
                float cos_theta = dot(-1.f * normalize(wi.d), n) < 1.f ? dot(-1.f * normalize(wi.d), n) : 1.f;
                if (std::pow((1-refr_ratio)/(1+refr_ratio),2) + (1.f-std::pow((1-refr_ratio)/(1+refr_ratio),2))*(std::pow(1.f - cos_theta, 5.f)) > (rand() / (float)RAND_MAX)) { 
                    return (normalize(wi.d) - 2 * dot(normalize(wi.d), n) * n); 
                } else { 
                    return refr_ratio * (normalize(wi.d) - dot(normalize(wi.d), n) * n) - std::sqrt(discrim) * n;
                }
            } default : { return Vector({0.f, 0.f, 0.f}); }
        }
    }
    float r, emit;
    SurfaceType_t type;
    Vector c, atten;
} Surface;
Vector trace(Ray& ray, Surface** surfaces, int length, int depth) {
    float hit_ind = 0, t = FLT_MAX;
    for(int i = length; i--;) if((t = surfaces[i]->intersect(ray)) && t < ray.maxt) {ray.maxt = t; hit_ind = i;} 
    if (ray.maxt == FLT_MAX) { return Vector({0, 0, 0}); }
    if (depth >= 64) { return surfaces[(int)hit_ind]->emit*Vector({7.5, 7.5, 7.5}); }
    Ray scattered = Ray({ray.o + ray.maxt*ray.d, surfaces[(int)hit_ind]->scatter(ray, ray.o + ray.maxt*ray.d), FLT_MAX});
    Vector out = surfaces[(int)hit_ind]->emit*Vector({7.5, 7.5, 7.5}) + surfaces[(int)hit_ind]->atten*trace(scattered, surfaces, length, depth + 1);
    out.x = out.x < 0.f ? 0.f : (out.x > 1.f) ? 1.f : out.x; 
    out.y = out.y < 0.f ? 0.f : (out.y > 1.f) ? 1.f : out.y; 
    out.z = out.z < 0.f ? 0.f : (out.z > 1.f) ? 1.f : out.z; 
    return out;
}
int main(int argc, char** argv) {
    int width = 640, height = 480, num_samp = (argc == 2 ? atoi(argv[1]) : 20), num_surf = 7, count = 1; 
    Surface** surfaces = new Surface*[num_surf];
    surfaces[0] = new Surface(100.f, Vector({0.f, -100.202, -1.f}), Vector({0.8f, 0.8f, 0.8f}), 0.f, LAMBERTIAN);  //Floor    
    surfaces[1] = new Surface(100.f, Vector({-100.202, 0.f, -1.f}), Vector({0.8f, 0.28f, 0.28f}), 0.f, LAMBERTIAN);  //Left Wall 
    surfaces[2] = new Surface(100.f, Vector({100.202, 0.f, -1.f}), Vector({0.28f, 0.28f, 0.8f}), 0.f, LAMBERTIAN);  //Right Wall  
    surfaces[3] = new Surface(100.f, Vector({0.f, 100.202, -1.f}), Vector({1.f, 1.f, 1.f}), 1.f, LAMBERTIAN);  //Ceiling (Light)   
    surfaces[4] = new Surface(100.f, Vector({0.f, 0.f, -100.8}), Vector({0.8f, 0.8f, 0.8f}), 0.f, LAMBERTIAN);  //Back Wall  
    surfaces[5] = new Surface(0.168f / 2, Vector({-.10f, -.118, -0.6f}), Vector({0.9f, 0.9f, 0.9f}), 0.f, METAL);  //Metal   
    surfaces[6] = new Surface(0.168f / 2, Vector({.10, -.118, -0.5f}), Vector({1.f, 1.f, 1.f}), 0.f, DIELECTRIC);  //Glass   
    FILE *f = fopen("cornellbox.ppm", "w");  
    fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
    Vector *image = new Vector[width*height];
    for (int pix = 0; pix < width*height; pix++){
        if (pix % width == 0) { printf("\r%6.2f%% [%.*s%*s]", (count / (float)height) * 100.f, (int)((count / (float)height) * strlen(BAR)), BAR, 
				(int)(strlen(BAR) - (int)((count / (float)height) * strlen(BAR))), ""); count++; fflush(stdout); }
        Vector color = {0, 0, 0};
        for (int i = 0; i < num_samp; i++) { 
            float t_x = ((pix % width) + (rand() / (float)RAND_MAX)) / ((float)width - 1.f); //Parameter for LERP
            float t_y = ((pix / width) + (rand() / (float)RAND_MAX)) / ((float)height - 0.f); // Parameter for LERP
            float ray_dir_x = (-0.7145315005)*(1-t_x) + (0.7145315005)*t_x; //LERP Equation
            float ray_dir_y = (0.53589862537)*(1-t_y) + (-0.53589862537)*t_y; //LERP Equation
            Ray ray = Ray({Vector({0.f, 0.f, 0.f}), normalize(Vector({ray_dir_x, ray_dir_y, -1.f})), FLT_MAX});
            color = trace(ray, surfaces, num_surf, 0) + color;
        }
        image[pix] = Vector({(float)(pow(color.x / (float)num_samp,1/2.2f)*255+.5f), (float)(pow(color.y / (float)num_samp,1/2.2f)*255+.5f), (float)(pow(color.z / (float)num_samp,1/2.2f)*255+.5f)});
   } 
   for (int pix = 0; pix < width*height; pix++) { fprintf(f, "%d %d %d ", (int)image[pix].x, (int)image[pix].y, (int)image[pix].z); }
}

//References
// - Progress bar adapted from https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf
// - Optimizations adapted from https://www.kevinbeason.com/smallpt/ 
// - Based on code from personal github @JonahWeinbaum (cannot be made public as per class policy), 
//   although base code for the project is publicly available on https://github.com/cs87-dartmouth/darts-2022 
// - Optimizations in quadratic equation root solving based on methods from Ray Tracing in One Weekend by Peter Shirley