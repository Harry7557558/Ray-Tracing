# My AVI1O CPT: Creating an Animation

Create a PSA of any kind of animation, at least 30 seconds / 12 fps.

Our animation type is 3D animation / CGI.

**My Task: Rendering animation with my own ray-tracer**

**Work Due: June 18th, 2019**

**Final Product: [YouTube Video](https://www.youtube.com/watch?v=ZQ4x4_ZViqo&feature=youtu.be); [My Part](Final Product/CGI.mp4)**


# Creation Process

My art teacher took this CPT on May 27th. Its due date was said to be June 14th. According to the document, our goal is to make a 30-second PSA using stop motion, but my art teacher said it can be any type of animation (2D, 3D, flash, etc.) He gave a list of various topics we could choose, and drew us into 2-people groups to complete this CPT. 

At first, what I can do is barely rendering some simple geometrical shapes. (Ex. regular spheres, planes, triangles, etc.) But I decided to make what I like, which is CGI animation. My partner Carrie thought we can do Bullying. But I didn’t choose a topic from the list, but a topic of my own, which was Temptation. I think bullying is currently not a very serious problem in our school, but there are all kinds of temptations around us - video games, cigarettes and drugs, money, etc. I wrote a short story about a glass man tempting by 3 temptations - a golden coin, video game, and an image of a sexy woman, finally falls to the enormous depth and darkness beneath it. My decision is submitted the following day and can’t be changed. 

But soon I regret. I couldn’t construct that complex shapes. I didn’t know how to calculate the intersection of rays and spline surfaces. And in a console, without a real-time visualizer, it would be very hard to design objects. Carrie designed a very complex glass man and a golden coin with complicated patterns on paper. I didn’t use his golden coin, and discarded the man’s face with an excuse of “A man without a face can represent all kinds of people”. 

At that evening I learned tracing objects defined by SDF. Within one day, I completed SDF core of various geometrical shapes - spheres, standard and truncated cylinders / cones, toruses, standard boxes… I also implemented lots of operations of SDF - union / intersect / subtraction, rounding, scaling… I decided to use an affine-transformed box with rounding to construct components of the glass man. The SDF must be exact in order to round it. I got stuck on it until the following Tuesday, which we only have a week to complete it. 

By Thursday I made a glass man class (See [GlassMan.h](GlassMan.h), although it may contain a lot of bugs) with fixed coefficients of all body components. Since there’s no formula about human walking attitudes, I found a video of human walking on YouTube, split it into frames, and sent it to Carrie with a [tool](GetCoordinate.html) that measures the coordinates of a pixel, told her to help me measuring the coordinates of all body components of the man in each frame. But it didn’t work on her Mac computer. I tried to explain, but she got furious. She said I when I select this type of animation, I didn’t think of her at all. I should select an animation type that everyone of us can do. Our art teacher said he’ll make this project as an individual mark, he told us he needs our work dividing and schedule information. If I do all the things and she have done nothing, she won’t get any marks. She said some groups have finished all of this CPT and theirs look very nice, but I haven’t rendered even a single frame. I felt very upset, but I could do nothing. 

At that night, I didn’t sleep until 4:30 a.m. to measure all the coordinates, convert them into unit vectors, and interpolated them using Fourier series. I fell asleep on math class the following day. My math teacher asked if I stayed up late that night. I just nodded my head and said nothing. 

A warning comes that my religion CPT is due next Monday. It was supposed to start 3 weeks ago, but what I had done was nothing. I spent all class times doing those math works for my art CPT. The whole weekend I was working on my religion CPT and the final part of ESL CPT and paused this art CPT. 

On Tuesday, I finished rendering the first scene of all people walking on the path. It took me 7 hours and half to render. (In commercial software, it may take just half an hour. I need more optimizations.) The next task is to construct the golden coin. I implemented the SDF of 2D Bezier curves to construct the pattern of the coin. It was my first time to solve cubic equations, and it took me the whole night. (My work is due Friday!)

I fell asleep at 6:30 on Wednesday morning. As my mother woke me up, it was 8 o’clock, and I realized that I have been late for the EQAO exam. I rushed to school, and saw my math teacher was waiting for me anxiously. It was a very important math exam. It was a one-hour exam, and I had been late for 35 minutes. (I knew I did two multiple choice problems wrong out of 14 problems after this CPT due. I supposed to be full mark.) 

I apologized to Carrie online on Thursday and said I couldn’t finish it. She said the art teacher have changed the due date to the next Tuesday. My art teacher received a phone call. He said if I couldn’t finish it, Carrie can do the remaining part using stop motion. He said he knows CG wastes a lot of time. 

I rendered a title page on Saturday. My mother said she’ll help me measure the coordinate of a human running video. She finished it on Sunday evening, although it contains lots of mistakes. (Ex. 741→141, 53.5→23.5) I made a 2-second part of the glass man rushing toward the golden coin out of 7 hours of working. I fell asleep again. As I woke up, it was one o’clock in the afternoon. I have been absent for school. 

I downloaded and installed Adobe Premiere from a Chinese pirated software manufacturer. I decided to do post production at school, but I definitely don’t have that much time. On Monday night, I rendered my work last night, as well as 6 seconds of new works. Carrie sent me the background music and a short walking steps sound effect. He also sent me some photos she took for stop motion of the latter part. The final product was encoded on 8:30. I uploaded it and rushed to school. Carrie have wrote a full paragraph of analysis. We presented it just 20 minutes after it was released. 

It was a very successful presentation. The music have linked the beginning and latter part very well. We both got a mark of 49 out of 50. I have been sleeping for the next two classes, and my teachers didn’t mark me. 

My origin art teacher left two months after the semester began (some said he was hit by a car), and a new art teacher took his class. This art teacher loves animations, especially hand-painted animation. A rumor (I’m not sure if it’s true) says he dreamed to work for Walt Disney animation studio when he was young and achieved it in 2006. But as Disney purchased Pixar, Disney stopped making 2D animations and switched to 3D animations. My art teacher couldn’t find anything to do, and walked the plank in 2018. 


# Thoughts and Feelings

I didn’t finish this CPT, but I had a very great process. In the beginning, I thought I could never finish a 3D animation of my own. I have heard the difficulties of making animations. This time, although I didn’t finish it, but it was my first time using splines and rendering objects defined by SDF. In these two weeks, I used three full packages of scratch papers, and typed more than 300K sources. If there’s no pressure behind me, I could never make it. 

I think sometimes in corporations, we need to apologize to others. As I argued with Carrie that evening, for nearly a week we didn’t say anything. We just avoid eye contact when meet. If I didn’t apologize initiatively, she couldn’t help me with the stop motion and background music, I may stubbornly make my CGI and wouldn’t finish it by time, and finally both of us get no mark. 

Third, I think when starting a new project, we need to think of our abilities, and care about our partners. A group made a stop motion about climate change. They didn’t take much time, but their final effects were very shocking. Some of hand-painted animations last century were much better than 3D animations or CGI movies today. Modern technologies are not always better than old methods. 

Finally, I think I should care about my health, no matter how busy I am. I stayed up late, couldn’t wake up, and missed lots of courses and exams, and caught a cold during these days. 


# Future Plans

I didn’t finish this project, and I won’t continue working on it. I decided to continue my love in computer graphics, make a real-time designer with GUI, write clear sources, do more optimizations, and make more creative pictures or videos. It doesn't need to be my career, and I don’t need to compare myself with those CGI movies. I just keep it as a hobby, and be proud of it. 

I have just started learning coding last April, and now I have done great works. In the beginning, I didn’t know clearly what I can do. I have tried making small tools, computer algebra system, convolutional neural network, or even hacking. None of them are successful, and some of them are really boring. I have spent two months on drawing the Mandelbrot set, and now I didn’t use it at all. Three months ago, my friend in grade 12 helped me discover the interest in computer graphics. This hobby haven’t lasted long, but I think it wouldn’t grow bored that I can continue it all my lifetime. 

I know almost half of students in grade 11 or 12 don’t have a goal. They just spend all their free time playing video games, and their lives are chaos. But I have a goal, or a specific hobby. This goal will promote my working in math, physics, and English. I can learn to manage my time and interpersonal communications. 

From this CPT, I learned a lot. 

