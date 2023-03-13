#include "intersect.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <limits>


bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p) {
    float areaABC = glm::length(glm::cross(v1 - v0, v2 - v0)) / 2.0f;
    float areaPBC = glm::length(glm::cross(v1 - p, v2 - p)) / 2.0f;
    float areaPAC = glm::length(glm::cross(v0 - p, v2 - p)) / 2.0f;
    float areaPAB = glm::length(glm::cross(v0 - p, v1 - p)) / 2.0f;
    if (areaABC - areaPAB - areaPAC - areaPBC <= FLT_EPSILON && areaABC - areaPAB - areaPAC - areaPBC >= -FLT_EPSILON)
        return true;
    return false;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    if (glm::dot(plane.normal, glm::normalize(ray.origin + ray.direction)) <= FLT_EPSILON && glm::dot(plane.normal, glm::normalize(ray.origin + ray.direction)) >= -FLT_EPSILON) {
        return false;
    }
    float t = (plane.D - glm::dot(ray.origin, plane.normal)) / glm::dot(ray.direction, plane.normal);
    if (t <= FLT_EPSILON) {
        // std::cout << "false ";
        return false;
    }
    if (t > ray.t + FLT_EPSILON)
        return false;
    ray.t = t;
    // std::cout << "true ";
    return true;
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    if (glm::cross(v0 - v2, v1 - v2) == glm::vec3 { 0, 0, 0 })
        plane.normal = glm::vec3(1, 1, 1);
    else
        plane.normal = glm::normalize(glm::cross(v0 - v2, v1 - v2));
    plane.D = glm::dot(plane.normal, v0);
    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    // TODO: implement this function.
    float t = ray.t;
    Plane p = trianglePlane(v0, v1, v2);
    if (intersectRayWithPlane(p, ray) == false) {
        // std::cout << "plane intersect false ";
        return false;
    }
    if (pointInTriangle(v0, v1, v2, p.normal, ray.origin + ray.t * ray.direction) == false) {
        ray.t = t;
        // std::cout << "point in triangle false ";
        return false;
    }
    // std::cout << "true ";
    return true;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    // TODO: implement this function.
    return false;
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    return false;
}
