#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>

// samples a segment light source
// you should fill in the vectors position and color with the sampled position and color
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{
    float r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX); // random between 0 and 1
    position = segmentLight.endpoint1 * (1.0f - r) + segmentLight.endpoint0 * r;
    color = (1.0f - r) * segmentLight.color1 + r * segmentLight.color0;
}

// samples a parallelogram light source
// you should fill in the vectors position and color with the sampled position and color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{
    float a = static_cast<float>(rand()) / static_cast<float>(RAND_MAX); // random between 0 and 1
    float b = static_cast<float>(rand()) / static_cast<float>(RAND_MAX); // random between 0 and 1

    glm::vec3 col1 = a * parallelogramLight.edge01;
    glm::vec3 col2 = b * parallelogramLight.edge02;

    float total = glm::length(glm::cross(parallelogramLight.edge01, parallelogramLight.edge02));
    float d1 = glm::length(glm::cross(col1, col2));
    float d2 = glm::length(glm::cross(parallelogramLight.edge01 - col1, col2));
    float d3 = glm::length(glm::cross(parallelogramLight.edge01 - col1, parallelogramLight.edge02 - col2));
    float d4 = glm::length(glm::cross(col1, parallelogramLight.edge02 - col2));

    position = parallelogramLight.v0 + col1 + col2;
    color = (d1 / total) * parallelogramLight.color0 + (d2 / total) * parallelogramLight.color1 + (d3 / total) * parallelogramLight.color2 + (d4 / total) * parallelogramLight.color3;
}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableHardShadow || features.enableSoftShadow) {
        // compute the point
        glm::vec3 point = ray.origin + ray.t * ray.direction;

        // direction of point to light(samplePos)
        glm::vec3 vec = samplePos - point;
        float length = glm::length(vec);
        vec = glm::normalize(vec); // normalize it

        ray.direction = vec;
        ray.t = length;
        ray.origin = point + 0.0001f * vec;  //adding an offset

        if (bvh.intersect(ray, hitInfo, features)) {
            drawRay(ray, { 1, 0, 0 });
            return 0.0f;
        }
    }
    drawRay(ray, debugColor);
    return 1.0f;
}

// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)

// Lights are stored in a single array (scene.lights) where each item can be either a PointLight, SegmentLight or ParallelogramLight.
// You can check whether a light at index i is a PointLight using std::holds_alternative:
// std::holds_alternative<PointLight>(scene.lights[i])
//
// If it is indeed a point light, you can "convert" it to the correct type using std::get:
// PointLight pointLight = std::get<PointLight>(scene.lights[i]);
//
//
// The code to iterate over the lights thus looks like this:
// for (const auto& light : scene.lights) {
//     if (std::holds_alternative<PointLight>(light)) {
//         const PointLight pointLight = std::get<PointLight>(light);
//         // Perform your calculations for a point light.
//     } else if (std::holds_alternative<SegmentLight>(light)) {
//         const SegmentLight segmentLight = std::get<SegmentLight>(light);
//         // Perform your calculations for a segment light.
//     } else if (std::holds_alternative<ParallelogramLight>(light)) {
//         const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
//         // Perform your calculations for a parallelogram light.
//     }
// }
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading) {
        // If shading is enabled, compute the contribution from all lights.
        glm::vec3 pointColor { 0, 0, 0 };
        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);
                pointColor += testVisibilityLightSample(pointLight.position, pointLight.color, bvh, features, ray, hitInfo) * computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);
                // Perform your calculations for a point light.
            } else if (std::holds_alternative<SegmentLight>(light)) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);
                glm::vec3 lightPos = { 0.0f, 0.0f, 0.0f };
                glm::vec3 lightColor { 0, 0, 0 };
                for (int i = 1; i <= features.numberOfSamples; i++) {
                    sampleSegmentLight(segmentLight, lightPos, lightColor);
                    pointColor += testVisibilityLightSample(lightPos, lightColor, bvh, features, ray, hitInfo) * computeShading(lightPos, lightColor, features, ray, hitInfo);
                }
                pointColor = pointColor / (float)features.numberOfSamples;

                // Perform your calculations for a segment light.
            } else if (std::holds_alternative<ParallelogramLight>(light)) {
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                glm::vec3 lightPos { 0, 0, 0 };
                glm::vec3 lightColor { 0, 0, 0 };
                for (int i = 1; i <= features.numberOfSamples; i++) {

                    sampleParallelogramLight(parallelogramLight, lightPos, lightColor);
                    pointColor += testVisibilityLightSample(lightPos, lightColor, bvh, features, ray, hitInfo) * computeShading(lightPos, lightColor, features, ray, hitInfo);
                }
                pointColor = pointColor / (float)features.numberOfSamples; // averaging the color

                // Perform your calculations for a parallelogram light.
            }
        }
        // TODO: replace this by your own implementation of shading
        return glm::vec3 { std::min(pointColor.x, 1.0f), std::min(pointColor.y, 1.0f), std::min(pointColor.z, 1.0f) };
    } else if (features.enableHardShadow) {
        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);
                testVisibilityLightSample(pointLight.position, pointLight.color, bvh, features, ray, hitInfo);
            }
        }

    } else if (features.enableSoftShadow) {
        for (const auto& light : scene.lights) {
            if (std::holds_alternative<SegmentLight>(light)) {
                const SegmentLight segmentLight = std::get<SegmentLight>(light);
                glm::vec3 lightPos = segmentLight.endpoint0;
                glm::vec3 lightColor { 0, 0, 0 };
                for (int i = 1; i <= features.numberOfSamples; i++) {
                    float r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                    sampleSegmentLight(segmentLight, lightPos, lightColor);
                    testVisibilityLightSample(lightPos, lightColor, bvh, features, ray, hitInfo);
                }
            }
            if (std::holds_alternative<ParallelogramLight>(light)) {
                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                glm::vec3 lightPos { 0, 0, 0 };
                glm::vec3 lightColor { 0, 0, 0 };
                for (int i = 1; i <= features.numberOfSamples; i++) {
                    sampleParallelogramLight(parallelogramLight, lightPos, lightColor);
                    testVisibilityLightSample(lightPos, lightColor, bvh, features, ray, hitInfo);
                }
            }
        }
    } else {
        // If shading is disabled, return the albedo of the material.
        return hitInfo.material.kd;
    }
}
