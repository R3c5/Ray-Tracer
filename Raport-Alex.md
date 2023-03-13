# Report

# 3.1 Shading

The shading part is implemented over 4 functions: 

```computeShading(...)``` in shading.cpp  
```computeLightContribution(...)``` in light.cpp  
```getFinalColor(...)``` in render.cpp  
`intersect` in boundingVolumeHierarcy.cpp

## computeShading
This function applies the Phong model to the point of intersection between the ray and surface by adding the diffuse and specular components together:

$` Diffuse = I \cdot K_d \cdot \cos(\theta) `$  
$` Specular = I \cdot K_s \cdot (\cos\phi)^n `$  
 
Where:  
- $\theta$ is the angle formed by light direction and the surface normal  
- $\phi$ is the angle between direction to the camera and the reflection of the light direction over the surface normal  
- n is the shininess of the surface

The $\cos$ is calculated using the dot product over normalized vectors. It is also used to check whether the light is coming from behind the surface.

## computeLightContribution

If shading is enabled, the function will calculate the influence of each light source over the point on the surface by invoking ```computeShading``` for each light and adding the results. Moreover, by using ```testVisibilityLightSample```,  it checks if a light source actually has influence over that point (is visible from that point or it is in shadow).  Otherwise it will return the albedo of the material.


## getFinalColor 
This function calls ```computeLightContribution``` for each ray that intersects a surface.

## intersect

Here we compute the normal for each intersected surface. This is done by doing the cross product of 2 intersecting sides of the intersected triangle. Afterwards, this vector is, of course, normalised.

## Additional info

For the visual debug, a ray of the computed color will be drawn(see pictures).
Sources used: Lecture 04 - Shading. (n.d.). Brightspace. https://brightspace.tudelft.nl/d2l/le/content/499418/viewContent/2765037/View

<img src="Report_images\Shading_ray_tracing_cornell_box_with_mirror.PNG" alt="Shading_1" width="600"/>

### Shading, ray trace mode, Cornell Box (with mirror)

<img src="Report_images\Shading_ray_tracing_cornell_box_with_mirror_and_parallelogram_light.PNG" alt="Shading_2" width="600"/>

### Shading, 257 light samples, ray-trace mode, Cornell Box (parallelogram and mirror)

<img src="Report_images\Shading_rasterization_cornell_box_with_mirror_red.PNG" alt="Shading_3" width="600"/>

### Shading, rasterization mode, Cornell Box (with mirror)

<img src="Report_images\Shading_rasterization_cornell_box_with_mirror_green.PNG" alt="Shading_4" width="600"/>

### Shading, rasterization mode, Cornell Box (with mirror)


# 3.2 Recusive ray-tracer

The recursive ray-tracer is split over 3 functions:  
```getFinalColor(...)``` in render.cpp  
```computeReflectionRay(...)``` in shading.cpp  
`intersect` in boundingVolumeHierarcy.cpp
## getFinalColor

If the recursive ray-tracer is enabled, the function will recursively calculate the reflection of a ray(while the materials hit are reflective, meaning $K_s \ne 0$  and `rayDepth` $\ne 0$ ) using `computeReflectionRay`. For each ray, it will compute the color(accounting for the effect of all the reflected rays by multiplying the color computed for the reflected ray with $K_s$), draw it, and return the computed color to be used for the color computation of the previous rays. The color is computed using `computeLightContribution`.   

## computeReflectionRay

In this function the reflected ray is computed by calculating each of the 3 components:  
- t - it is set to infinity(or max float) because the ray didn't intersect any surface yet(this will happen in get final color using the `bvh.intersect`)  
- direction - it is computed by reflecting the incoming ray over the normal using the following formula: $r = d - 2(d\cdot n)n$, where $r$ is the reflected ray direction, $d$ is the incoming ray direction and $n$ is the surface normal.
- origin - it is the point of intersection of the incoming ray with the surface displaced by an epsilon(really small value) in the direction of the reflection. The displacement is done to prevent auto-intersection(reflected ray intersects with the incoming ray resulting in the reflected ray's t to be equal to 0).

## intersect
!Same as shading(This is used for both features)
Here we compute the normal for each intersected surface. This is done by getting the direction of the perpendicular line to the surface by doing the cross product of 2 intersecting sides of the intersected triangle. This is, of course, normalised.
## Additional info
A ray starts with a certain, ```rayDepth```, assigned in render.h and every time it is reflected, the ```rayDepth``` is decreased. When it reaches 0, it will no longer be reflected. In the submited solution, the default `rayDepth` is 5.  
Visual debug shows a ray and its reflection and in ray-trace mode a working mirror(see pictures)


<img src="Report_images\Recursive_ray-tracer_reflected_ray_rasterization.PNG" alt="Recursive-1" width="600"/>

### Recursive ray-tracing, rasterization mode, Cornell Box (with mirror)

<img src="Report_images\Recursive_ray-tracer_reflected_ray_rasterization_to_inf.PNG" alt="Recursive-2" width="600"/>

### Recursive ray-tracing, rasterization mode, Cornell Box (with mirror)

<img src="Report_images\Recursive_ray-tracer_reflected_ray_rasterization_in_cube.PNG" alt="Recursive-3" width="600"/>

### Recursive ray-tracing, rasterization mode, Cornell Box (with mirror)

<img src="Report_images\Recursive_ray_tracer_working_mirror_1.PNG" alt="Recursive-4" width="600"/>

### Recursive ray-tracing, shading, ray-trace mode, Cornell Box (with mirror)

<img src="Report_images\Recursive_ray_tracer_working_mirror_2.PNG" alt="Recursive-5" width="600"/>

### Recursive ray-tracing, shading, ray-trace mode, Cornell Box (with mirror)

# 3.6 Barycentric coordinates for normal interpolation

Split over 3 methods:

`computeBarycentricCoord(...)` in interpolate.cpp  
`interpolateNormal(...)` in interpolate.cpp  
`intersect` in boundingVolumeHierarcy.cpp

## computeBarycentricCoord 

Calculates the barycentric coordinates dividing the areas of the 3 triangles formed by a
point P inside a triangle ABC with the area of the triangle ABC

## interpolateNormal

Calculates the interpolated normal by multipling each vertex normal with the coresponding barycentric coordinate

## intersect

Here we add the visual debug in 2 ways:
- In rasterization mode: draw the 3 vertice normals and the interpolated normal computed using the vertices from the intersected triangle and the point in which the triangle was intersected.  
- In ray-trace mode: use the interpolated normal instead of the surface normal

<img src="Report_images\Normal_interpolation_triangle_monkey_rasterization.PNG" alt="Normal_interp_1" width="600"/>

### Normal interpolation, Rasterizaton mode, Monkey

<img src="Report_images\Normal_interpolation_triangle_1.PNG" alt="Normal_interp_2" width="600"/>

### Normal interpolation, Rasterizaton mode, Single triangle

<img src="Report_images\Normal_interpolation_triangle_2.PNG" alt="Normal_interp_3" width="600"/>

### Normal interpolation, Rasterizaton mode, Single triangle

<img src="Report_images\No-normal_interp_ray_tracing.PNG" alt="Normal_interp_4" width="600"/>

### NO Normal interpolation, Shading, Ray-trace mode, Monkey

<img src="Report_images\Normal_interpolation_monkey_ray_tracing.PNG" alt="Normal_interp_5" width="600"/>

### Normal interpolation, Shading, Ray-trace mode, Monkey

# 3.7 Textures 
`interpolateTexCoords(...)` in interpolate.cpp  
`aquireTexel(...)` in texture.cpp  
`intersect` in boundingVolumeHierarcy.cpp
## interpolateTexCoords
Computes the texture coordinates of a certain point in a triangle by multipling each vertex's texture coord with its respective barycentric coordinate

## aquireTexel

I begin by computing the image's pixel coordinates which is done by multiplying the texture coordinates with the image's width in pixels, and respectively with the height. After finding the pixel coordinates, I retrive and return the pixel color from the image's vector by using the following formula $width \cdot j + i$

## intersect

If textures are enabled, I map the retrived color from the texture to its respective point

Sources used: Lecture 2 - Images and Algebra. (n.d.). Brightspace. https://brightspace.tudelft.nl/d2l/le/content/499418/viewContent/2764991/View

<img src="Report_images\Textures_quad.PNG" alt="Texture_1" width="600"/>\
<img src="Report_images\Textures_cube.PNG" alt="Texture_2" width="600"/>

# Extra Features

# Bloom Filter
This part is implemented in the following functions:  
`renderRayTracing` in render.cpp  
`bloomFilter`in render.cpp  
`thresholdScreen`in render.cpp  
`applyBoxFilter`in render.cpp  
`boxFilter`in render.cpp  
`scaleScreen`in render.cpp  
main.cpp
common.h  
render.h
## Main
In main, I added 3 sliders that control the filter size, the scaling factor, and the threshold for the visual debug and 3 buttons(they must be used individually, if more than one is activated, only the higher one in the list will work) that allow you to see each component of the bloom filter(threshold, boxFilter, scaling) individually.

## common.h
I added some flags in order to use the buttons in main.
## render.h
I added some extern variables in order to pass the slider's values to the function.

## renderRayTracing

Passes a reference of the screen to `bloomFIlter` if the bloom effects are enabled

## bloomFilter

Takes the screen and applies 3 steps(treshold, boxFilter, scale) in order to create the bloom filter. Each step is done in a separate function(`thresholdScreen`, `applyBoxFilter`, `scaleScreen`). If any of the 3 visual debug buttons are pressed, the function passes the screen to the respective function and assigns the result to the screen. If none of them are pressed, it adds the bloom filter to the screen's values.

## thresholdScreen

Recevies a `Screen` and goes through its pixels assigning the ones whose average of the RGB components is above the threshold to a black `Screen` . Returns the created `Screen`. 

## applyBoxFilter
Goes through the received `Screen`'s pixels and applies boxFilter to each one of them by using the `boxFilter` function, then it assigns them to a new screen which, in the end, is returned.

## scaleScreen

Goes through the received `Screen`'s pixels and multiplies them with the scaling factor, assigning them to a new screen, which, at the end, is returned.

## Additional info
If none of the buttons are pressed, `bloomFilter` applies the box filter and the scaling on the thresholded screen and then adds it back to the original screen, but if any of the debug buttons is pressed, it will applied the specified function directly to the screen. If more than one button is pressed, the highest one in the list will be displayed.  
Source: Lecture 2 - Images and Algebra. (n.d.). Brightspace. https://brightspace.tudelft.nl/d2l/le/content/499418/viewContent/2764991/View

<img src="Report_images\Bloom_high_threshold.PNG" alt="Bloom_1" width="600"/>\
<img src="Report_images\Bloom_low_threshold_filter1.PNG" alt="Bloom_2" width="600"/>\
<img src="Report_images\Bloom_low_threshold_filter7.PNG" alt="Bloom_3" width="600"/>

### !!!The pictures above show the bloom filter with different parameters. The first one has a higher threshold than the 2 bellow. The last 2 pictures have the same threshold but different filter size(the middle one has 1, the last has 7).  
### !!! The pictures bellow show each part of the filter individually(can be accesed using the designated buttons). First one shows only the threshold function applied directly to the screen, the second one shows the box filter applied directly to the screen and the third one shows a low scaling factor applied directly on the screen.

<img src="Report_images\Bloom_only_threshold_low.PNG" alt="Bloom_4" width="600"/>\
<img src="Report_images\Bloom_only_boxfilter.PNG" alt="Bloom_5" width="600"/>\
<img src="Report_images\Bloom_only_scale_low_scale_factor.PNG" alt="Bloom_6" width="600"/>


# Glossy Reflections
`getFinalColor` in render.cpp
`getGlossyRays` in shading.cpp
`getGlossyRay` in shading.cpp

## getGlossyRay
This function generates a orthonormal basis that defines a square perpendicular to the reflected ray. Afterwards, using a uniform real distribution and the degrees of blur(1/material.shininess) generates 2 coefficients `uCoefficient` and `vCoefficient`. These 2 coefficients are then multiplied with the orthonormal basis and added to the direction of the reflected ray. This operation perturbs the direction of the ray by no more than `a`/2 in each direction. This ray is returned to `getGlossyRays`

## getGlossyRays

This function samples a vector of 500 glossyRays that will be used in `getFinalColor` to average the color of a ray. By averaging the color of the square around the ray, we create the effect of a non-ideal mirror, such as a brushed metal.
  
  Source used: Marschner, S., & Shirley, P. (2016). Fundamentals of Computer Graphics, Fourth Edition. CRC Press. , pages 333-334
<img src="Report_images\Glossy_reflections_sampling.PNG" alt="Glossy_1" width="600"/>

### This picture highlights how a vector of 500 rays averges the color of a square of width `a`(degrees of blur)

<img src="Report_images\Glossy_reflections.png" alt="Glossy_2" width="600"/>

### This picture shows how a glossy surface is rendered using ray trace mode.
 
# Transparency

`getFinalColor` in render.cpp
! Needs recursive ray tracer in order to work.  
If transparency is enabled, it will check if the surface is transparent and that the ray hasn't reached its maximum depth. If the conditions are met, it will shoot a ray through the transparent panel and it will calculate its color using alpha blending. When the ray is extended, the origin of the extension is displaced by a factor of epsilon(a really small factor) in its direction, in order to avoid auto intersection.

<img src="Report_images\Transparency_raster_box.PNG" alt="Transparency_1" width="600"/>\
<img src="Report_images\Transparency_raster_box_inside.PNG" alt="Transparency_2" width="600"/>

### In this 2 pictures, the same ray can be seen from outside the cube(first picture) and from inside the cube(second picture). The ray that enters the triangle has slightly lighter color than the one inside(better seen if looking at the 2 dots in the second picture) because it alpha blends both cube faces that it passes and the background, while the one inside only alpha blends the green face with the background. The ray that exists is red beacuse it does not intersect any surface and therfore it has an infinite length.

<img src="Report_images\Transparency_cube.PNG" alt="Transparency_3" width="600"/>

### This last picture shows a transparent cube in ray trace mode, highlighting the alpha blending between the faces.


