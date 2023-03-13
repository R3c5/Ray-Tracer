#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>
#include <iostream>
#include <random>

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    // TODO: implement the Phong shading model.
    if (features.enableShading) {
        glm::vec3 vertexPos = ray.origin + ray.direction * ray.t;

        glm::vec3 diffuseComponent;
        float dotProductDiffuse = glm::dot(glm::normalize(hitInfo.normal), glm::normalize(lightPosition - vertexPos));
        if (dotProductDiffuse < 0)
            diffuseComponent = glm::vec3 { 0.0f };
        else
            diffuseComponent = lightColor * hitInfo.material.kd * dotProductDiffuse;

        glm::vec3 specularComponent;
        Ray reflectionRay = computeReflectionRay(ray, hitInfo);
        glm::vec3 reflectionRayVector = reflectionRay.origin + reflectionRay.t * reflectionRay.direction;
        float dotProductSpecular = glm::dot(glm::normalize(ray.origin - vertexPos), glm::reflect(glm::normalize(vertexPos - lightPosition), glm::normalize(hitInfo.normal)));  
        if (dotProductSpecular < 0)
            specularComponent = glm::vec3 { 0.0f };
        else {
            
            specularComponent = lightColor * hitInfo.material.ks * pow(dotProductSpecular, hitInfo.material.shininess);
        }
            

        return diffuseComponent + specularComponent;
    }
    return hitInfo.material.kd;
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    Ray reflectionRay {};
    glm::vec3 givenRay = ray.origin + ray.direction * ray.t;
    glm::vec3 reflection = glm::normalize(ray.direction) - 2 * glm::dot(glm::normalize(ray.direction), hitInfo.normal) * hitInfo.normal;
    float eps = 1e-5;
    reflectionRay.origin = givenRay + eps * reflection;
    reflectionRay.direction = glm::normalize(reflection);
    reflectionRay.t = std::numeric_limits<float>::max();
    // TODO: implement the reflection ray computation.
    return reflectionRay;
}

std::vector<Ray> getGlossyRays(Ray reflectionRay, HitInfo hitInfo) {
    int amountOfRays = 500;
    float a = 1 / hitInfo.material.shininess; //degrees of blur

    std::vector<Ray> rays;
    for (int i = 0; i < amountOfRays; i++) {
        rays.push_back(getGlossyRay(reflectionRay, a));
    }
    return rays;
}

Ray getGlossyRay(Ray reflectionRay, float a) {
    std::srand(time(NULL));
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<float> distribution(0.0f, std::nextafter(1.0f, 2.0f));
    float xi1 = distribution(generator);
    float xi2 = distribution(generator);
    // std::cout << xi1 << " " << xi2 << "\n"; //xi is a greek letter
    float uCoefficient = -(1.0f / 2.0f) * a + xi1 * a;
    float vCoefficient = -(1.0f / 2.0f) * a + xi2 * a;
    glm::vec3 w = glm::normalize(reflectionRay.direction);
    glm::vec3 t = w;
    if (t.x <= t.y && t.x <= t.z) {
        t.x = 1.0f;
    } else if (t.y <= t.x && t.y <= t.z) {
        t.y = 1.0f;
    } else {
        t.z = 1.0f;
    }
    t = glm::normalize(t);
    glm::vec3 u = glm::normalize(glm::cross(w, t));
    glm::vec3 v = glm::normalize(glm::cross(w, u));
    Ray glossyRay = reflectionRay;
    glossyRay.direction = reflectionRay.direction + uCoefficient * u + vCoefficient * v;
    return glossyRay;
}

//const Ray computeReflectionRay(Ray ray, HitInfo hitInfo)
//{
//    // Do NOT use glm::reflect!! write your own code.
//    Ray reflectionRay {};
//    glm::vec3 givenRay = ray.origin + ray.direction * ray.t;
//    glm::vec3 reflection = glm::normalize(ray.direction) - 2 * glm::dot(glm::normalize(ray.direction), hitInfo.normal) * glm::normalize(hitInfo.normal);
//    reflection = glm::normalize(reflection);
//    reflectionRay.origin = ray.origin + FLT_EPSILON * ray.direction + ray.t * ray.direction;
//    reflectionRay.t = FLT_MAX;
//    reflectionRay.direction = reflection;
//    // TODO: implement the reflection ray computation.
//    return reflectionRay;
//}