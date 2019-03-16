# Ray-Tracing
<img src="https://raw.githubusercontent.com/Harry7557558/Ray-Tracing/master/RT005.png" alt="" />

Ray-tracing, without using any third-party libraries. <br/>
A image of this large requires about 20 seconds and 20M memories. (4 threads, Intel(R) Core(TM) i7-8550U CPU @ 1.80GHz, Windows 10 home)

Core sources are in "Object.h", other source files are debugging. <br/>
Support intersection of vector planes, triangles, parallelograms, spheres, circles, and cylinders. <br/>
Temporary don't support diffuse reflection and refraction. 

Compare of non-rendered and rendered images: 
<img src="https://raw.githubusercontent.com/Harry7557558/Ray-Tracing/master/compare.png" alt="" />

Some result images: <br/>
<img src="https://raw.githubusercontent.com/Harry7557558/Ray-Tracing/Image/RayTracing_3.f.05%2B.png" width="420px" height="280px" alt="" />
<img src="https://raw.githubusercontent.com/Harry7557558/Ray-Tracing/Image/RayTracing_3.f.02%2B.png" width="420px" height="280px" alt="" />
<img src="https://raw.githubusercontent.com/Harry7557558/Ray-Tracing/Image/RayTracing_4.e.png" width="420px" height="280px" alt="" />
<img src="https://raw.githubusercontent.com/Harry7557558/Ray-Tracing/Image/RayTracing_5.c.png" width="420px" height="280px" alt="" />
<img src="https://raw.githubusercontent.com/Harry7557558/Ray-Tracing/Image/%CE%93.png" width="420px" height="420px" alt="" />

The first and the second consists of over 400,000 triangles. Each needs 16-20s but more than 1.6G memories; (Since mutex lock don't work well in this kind of situation, I need to copy data of all 400,000 objects for 4 threads)

The third and fourth image are vector. Each requires less than 2 seconds and less than 2M memories; 

The last one is complex Γ function. Some details were not very well, because my computer don't have that much memories for that much triangles.  :-(

I had done a lot of optimization in CPU usage. I used to spent more than 20 minutes to render an object consists of 6000 triangles. Now, I can do it in 5 seconds. But it requires 4 times of memory. 

More images are in "Image/comment.zip". 
