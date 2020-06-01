
#include "Ray.h"
#include "hitable.h"
#include "Raytracing.h"

#include <stdlib.h>  // Needed for drand48()
#include <random>
#include <cmath>

namespace rt {

glm::vec3 unit_vector(glm::vec3 v) {
    float v_length= sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return v / v_length;
}

float dot(const glm::vec3 &u, const glm::vec3 &v) {
    return u[0] * v[0]
         + u[1] * v[1]
         + u[2] * v[2];
}
glm::vec3 reflect(const glm::vec3& v, const glm::vec3& n) {
     return v - 2*dot(v,n)*n;
}

glm::vec3 random_in_unit_sphere_mat(){
    glm::vec3 p;
    do{
        p = 2.0f*glm::vec3(rand(),rand(),rand()),glm::vec3(1,1,1);
    }while((p[0]*p[0]+p[1]*p[1],p[2]*p[2]) >= 1,0);
    return p;
}



float random_double_aux() {
    return rand() / (RAND_MAX + 1.0);
}


float random_double(double min, double max) {
    // Returns a random real in [min,max).
    return min + (max-min)*random_double_aux();
}

float length_squared(glm::vec3 e)  {
    return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
}
float schlick(float cosine, float ref_idx) {
    auto r0 = (1-ref_idx) / (1+ref_idx);
    r0 = r0*r0;
    return r0 + (1-r0)*pow((1 - cosine),5);
}


/*

glm::vec3 random_unit_vector() {
    auto a = random_double_aux();
    auto z = random_double_aux();
    auto r = sqrt(1 - z*z);
    return glm::vec3(r*cos(a), r*sin(a), z);
}
*/
glm::vec3 random_unit_vector() {
    auto a = random_double(0, 2*3.14);
    auto z = random_double(-1, 1);
    auto r = sqrt(1 - z*z);
    return glm::vec3(r*cos(a), r*sin(a), z);
}
bool refract(const glm::vec3 &v, const glm::vec3 &n, float ni_over_nt, glm::vec3 &refracted){
    glm::vec3 uv = glm::normalize(v);
    float dt = glm::dot(uv,n);
    float discriminat = 1.0 - ni_over_nt*ni_over_nt*(1-dt*dt);
    if(discriminat>0){
        refracted = ni_over_nt*(uv-n*dt)-n*sqrt(discriminat);
        return true;
    }
    else {
        return false;
    }
}

class material {
    public:
        virtual bool scatter(const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered) const = 0;
};

class metal : public material {
    public:
        metal(const glm::vec3& a,float f) : albedo(a) ,fuzz(f < 1 ? f : 1)  {}

        virtual bool scatter(const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered) const {
            glm::vec3 reflected = glm::reflect(glm::normalize(r_in.direction()),rec.normal);

            scattered = Ray(rec.p, (reflected+  fuzz*random_unit_vector()));//<- funkar inte

            attenuation = albedo;
            return (glm::dot(scattered.direction(), rec.normal) > 0.0f);
        }

    public:
        glm::vec3 albedo;
        float fuzz;
};


class lambertian  : public material {
    public:
        lambertian(const glm::vec3& a) : albedo(a) {}

        virtual bool scatter(const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered) const {
            glm::vec3 reflected = reflect(glm::normalize(r_in.direction()),rec.normal);
            glm::vec3 scatter_direction = rec.normal +  random_unit_vector();
            scattered = Ray(rec.p,scatter_direction);

            attenuation = albedo;
            return true;
        }

    public:
        glm::vec3 albedo;
};
/*
class dielectric : public material {
    public:
        dielectric(float ri) : ref_idx(ri) {}
    virtual bool scatter(const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered) const {

        
        attenuation = glm::vec3(1.0, 1.0, 1.0);
        auto etai_over_etat = (rec.front_face) ? (1.0 / ref_idx) : (ref_idx);

        glm::vec3 unit_direction = glm::normalize(r_in.direction());
        auto cos_theta = float(fmin(dot(-unit_direction, rec.normal), 1.0));
        auto sin_theta = float(sqrt(1.0 - cos_theta*cos_theta));
        if (etai_over_etat * sin_theta > 1.0 ) {
            glm::vec3 reflected = glm::reflect(glm::normalize(r_in.direction()), rec.normal);

            scattered = Ray(rec.p, reflected);
            return true;
        }
    
        float reflect_prob = schlick(cos_theta, etai_over_etat);

        if (random_double_aux() < reflect_prob)
        {
            glm::vec3 reflected = glm::reflect(glm::normalize(r_in.direction()), rec.normal);

            scattered = Ray(rec.p, reflected);
            return true;
        }
            
        glm::vec3 refracted = refract(glm::normalize(r_in.direction()), rec.normal, etai_over_etat);
        scattered = Ray(rec.p, refracted);

        return true;
        }
    public:
        float ref_idx;
};
*/

class dielectric : public material {
    public:
        dielectric(float ri) : ref_idx(ri) {}
        virtual bool scatter(const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered) const {
        glm::vec3  outward_normal;
        glm::vec3 reflected = reflect(glm::normalize(r_in.direction()),rec.normal);
        float ni_over_nt;
        
        
        attenuation = glm::vec3(1.0, 1.0, 1.0);
        glm::vec3 refracted;
        float reflect_prob;
        float cosi;
        if(glm::dot(r_in.direction(),rec.normal)>0){
            outward_normal = -rec.normal;
            ni_over_nt = ref_idx;
            cosi = ref_idx*glm::dot(r_in.direction(),rec.normal)/r_in.direction().length();
        }
        else {
            outward_normal = rec.normal;
            ni_over_nt = 1.0f / ref_idx;
            cosi = -glm::dot(r_in.direction(),rec.normal)/r_in.direction().length();

        }
        if(refract(r_in.direction(),outward_normal,ni_over_nt,refracted)){
            reflect_prob = schlick(cosi,ref_idx);
        }
        else{
            reflect_prob = 1.0f;
        }
        if(random_double_aux()<reflect_prob){
            scattered = Ray(rec.p,reflected);
        }
        else {
            scattered = Ray(rec.p,refracted);

        }
        return true;
        }
    public:
        float ref_idx;
};

class ground_sphere  : public material {
    public:
        ground_sphere(const glm::vec3& a) : albedo(a) {}

        virtual bool scatter(const Ray &r_in, const HitRecord &rec, glm::vec3 &attenuation, Ray &scattered) const {
            attenuation = albedo;

            return true;
        }

    public:
        glm::vec3 albedo;

};

} // namespace rt