#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif
#include <nativefiledialog/nfd.h>

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);
    
        if (features.enableRecursive) {
            std::vector<Ray> rays;
            if (hitInfo.material.ks != glm::vec3{ 0, 0, 0 } && rayDepth != 0) {
                Ray reflection = computeReflectionRay(ray, hitInfo);
                if (features.extra.enableGlossyReflection) {
                    rays = getGlossyRays(reflection, hitInfo);
                    for (int i = 0; i < rays.size(); i++) {
                        if (glm::dot(hitInfo.normal, rays[i].direction) > 0.0f && std::sqrt(1.0f - glm::dot(hitInfo.normal, rays[i].direction) * glm::dot(hitInfo.normal, rays[i].direction)) > 0.0f)
                        Lo += getFinalColor(scene, bvh, rays[i], features, rayDepth - 1) * (1.0f / static_cast<float>(rays.size())) * hitInfo.material.ks;
                    }
                } else {
                    Lo = Lo + getFinalColor(scene, bvh, reflection, features, rayDepth - 1) * hitInfo.material.ks;
                }
            }  

            if (features.extra.enableTransparency) {
                if (hitInfo.material.transparency != 1.0f && rayDepth != 0) {
                    Ray extension;
                    float eps = 1e-5f;
                    extension.origin = ray.origin + ray.direction * ray.t + ray.direction * eps;
                    extension.direction = ray.direction;
                    Lo += getFinalColor(scene, bvh, extension, features, rayDepth - 1) * (1 - hitInfo.material.transparency);
                }
            }
            // TODO: put your own implementation of recursive ray tracing here.
        }

        // Draw a white debug ray if the ray hits.
        drawRay(ray, Lo);
        //drawRay(computeReflectionRay(ray,hitInfo), Lo);
        
        // Set the color of the pixel to white if the ray hits.
        return Lo;
    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}


#define numSamples 20

float randomOffset(float aperture)
{
    return ((2 * aperture) * ((float)rand() / RAND_MAX)) - aperture;
}
glm::vec3 executeDepthOfField(const Scene& scene, Ray ray, const BvhInterface& bvh, const Features& features, float focal_length, float aperture)
{
    glm::vec3 focal_point = ray.origin + focal_length * ray.direction;

    glm::vec3 sum_color = getFinalColor(scene, bvh, ray, features);
    for (int i = 0; i < numSamples; i++) {
        Ray sec;
        sec.origin = glm::vec3(ray.origin.x + randomOffset(aperture), ray.origin.y + randomOffset(aperture), ray.origin.z + randomOffset(aperture));
        sec.direction = glm::normalize(focal_point - sec.origin);
        sec.t = std::numeric_limits<float>::max();

        sum_color += getFinalColor(scene, bvh, sec, features);
    }
    glm::vec3 final_color = (1.0f / numSamples) * sum_color;
    final_color = glm::min(final_color, glm::vec3(1));
    return final_color;
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
     

    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            glm::vec3 color = getFinalColor(scene, bvh, cameraRay, features);
            if (features.extra.enableDepthOfField) {
                color = executeDepthOfField(scene, cameraRay, bvh, features, features.dof.focal_length, features.dof.aperture);
            }
            screen.setPixel(x, y, color);
        }
    }
    
    if (features.extra.enableBloomEffect) {
        bloomFilter(screen, features);
    }
}

void bloomFilter(Screen& screen, Features features) {
    glm::ivec2 windowResolution = screen.resolution();
    Screen screenThreshold = thresholdScreen(screen);
    Screen boxFilteredScreen = applyBoxFilter(screenThreshold);
    Screen scaledScreen = scaleScreen(boxFilteredScreen);

    if (features.extra.enableThreshold) {
        screen = screenThreshold;
    } else if (features.extra.enableBoxFilter) {
        screen = applyBoxFilter(screen);
    } else if (features.extra.enableScale) {
        screen = scaleScreen(screen);
    } else {
        for (int y = 0; y < windowResolution.y; y++) {
            for (int x = 0; x != windowResolution.x; x++) {
                screen.setPixel(x, y, screen.pixels()[screen.indexAt(x, y)] + scaledScreen.pixels()[screen.indexAt(x, y)]);
            }
        }
    }
}


Screen thresholdScreen(Screen& screen) {
    Screen screenThreshold = screen;
    glm::ivec2 windowResolution = screen.resolution();
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            glm::vec3 color = screen.pixels()[screen.indexAt(x, y)];
            if ((color.r + color.b + color.g) / 3.0f > threshold) {
                screenThreshold.setPixel(x, y, color);
            } else {
                screenThreshold.setPixel(x, y, glm::vec3 { 0.0f });
            }
        }
    }
    return screenThreshold;
}


Screen applyBoxFilter(Screen& source) {
    glm::ivec2 windowResolution = source.resolution();
    Screen result = source;
    for (int i = 0; i < windowResolution.x; i++) {
        for (int j = 0; j < windowResolution.y; j++) {
                result.setPixel(i,j, boxFilter(source, i, j));
        }
    }
    return result;
}

glm::vec3 boxFilter(Screen& source, int i, int j) {
    filtersize = std::max(1, filtersize);
    glm::vec3 sum { 0.0f };
    for (int x = -filtersize; x <= filtersize; x++) {
        for (int y = -filtersize; y <= filtersize; y++) {
            if (i + x >= 0 && j + y >= 0 && i + x < source.resolution().x && j + y < source.resolution().y) {
                sum += source.pixels()[source.indexAt(i + x, j + y)];
            } else {
                sum += glm::vec3 { 0.0f };
            }
            
        }
    }
    sum /= (2*filtersize + 1) * (2*filtersize + 1);
    return sum;
}

Screen scaleScreen(Screen& screen) {
    glm::ivec2 windowResolution = screen.resolution();
    Screen scaledScreen = screen;
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++){
            scaledScreen.pixels()[scaledScreen.indexAt(x,y)] = screen.pixels()[screen.indexAt(x, y)] * scale;
        }
    }
    return scaledScreen;
}