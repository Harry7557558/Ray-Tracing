# Ray-Tracing

![ ](Cover.jpg)

Ray-tracing, without using any third-party libraries.

An image of this large (1920x1080) requires about 20 seconds and 20M memories. (4 threads, Intel(R) Core(TM) i7-8550U CPU @ 1.80GHz, Windows 10 home)

All vector objects and intersection algorithms are in [Object.h](Object.h).

Support intersections of vector planes, triangles, parallelograms, spheres, circles, and cylinders; refraction of planes and spheres; 

Comparison of non-rendered and rendered images: 

![ ](compare.png)



# Triangles

Creating shapes with millions of triangles, although it may be low efficiency. 

![_](ring1.jpg)
![_](ring2.jpg)

The rings goes through the most wonderful optimization I have ever made. <br/>
<pre>CPU:     44min => 5.4s
Memory:  1.26GB => 56MB</pre>

![_](Γ.jpg)

This picture shows <a href="https://en.wikipedia.org/wiki/Gamma_function" target="_blank">Γ function</a> on complex plane. It consists of 1,600,000 triangles. It requires about 20 seconds and about 200M memories to render. 

![_](beads.jpg)
![_](pyramid.jpg)

This two pictures are vector spheres and cylindars, each requires less than 1 seconds and less than 1M memories. 



