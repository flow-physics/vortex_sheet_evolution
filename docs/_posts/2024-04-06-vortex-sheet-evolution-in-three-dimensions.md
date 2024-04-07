---
author: Rajesh Venkatesan
comments: true
date: 2024-04-06 00:00:00+00:00
layout: post
link: http://flowphysics.com/vortex_sheet_evolution/vortex-sheet-evolution-in-three-dimensions
slug: vortex-sheet-evolution-in-three-dimensions
title: Vortex Sheet Evolution in Three Dimensions
tags: vorticity 3D
---

# Vortex Sheet Evolution in three dimensions

## Outline of the problem

We start with,

$$
\begin{equation} \tag{1} % \label{eq:continuity} 
\nabla \cdot \vec{u} = 0
\end{equation} $$ 

$$ \begin{equation} \tag{2} % \label{eq:vort_def}
\nabla \times \vec{u} =\vec{\omega}
\end{equation} $$

Taking the curl of the Equation (2), we get,

$$ \begin{equation*}
\nabla \times \nabla \times \vec{u} = \nabla \times \vec{\omega}
\end{equation*}$$

$$ \begin{equation*}
\nabla (\nabla \cdot \vec{u}) -(\nabla \cdot \nabla) \vec{u} = \nabla \times \vec{\omega}
\end{equation*}$$

Using Equation (1), we have,

$$\begin{equation} \tag{3} % \label{eq:vel_vort}
\nabla ^2 \vec{u} = - \nabla \times \vec{\omega} 
\end{equation}$$

By solving this Poisson equation (3), we get the velocity field corresponding to a given vorticity field. Note that this is only a kinematic relationship and no time evolution is involved. The time evolution of vorticity in three dimensions is given by,

$$ \begin{equation} \tag{4} % \label{eq:vort_3d}
\frac{D\vec{\omega}}{Dt}=(\vec{\omega}.\nabla)\vec{u}+\nu\nabla^2\vec{\omega}
\end{equation}$$

If we want to see how a vortex sheet evolves in 3D, we need to solve Equation (4). The corresponding velocity field at each instant of time can be obtained from Equation (3). The last term in Equation (4) is the viscous diffusion of vorticity and will be neglected in our analysis.

## Two dimensional case

Some simplifications are possible in two dimensions. In this case, the dynamical equation (4) becomes,

$$ \begin{equation} \tag{5} % \label{eq:vort_2d}
\frac{D\omega_z}{Dt}=0
\end{equation}$$

The stretching and tilting terms have vanished and only one component of vorticity can exist in two-dimensions. This tells us that vorticity moves with the material element in 2D. So, the staretegy for solving the vortex sheet evolution in 2D is as follows:

- **Step (1):** Given the initial vorticity distribution $\omega_z$ at time $t_0$, we solve the equation (3) in the following form,

$$\begin{equation}
\nabla^2 u = -\frac{\partial \omega_z}{\partial y} 
\end{equation}$$

$$ \begin{equation}
\nabla^2 v = \frac{\partial \omega_z}{\partial x} 
\end{equation}$$

- **Step (2):** We find the displacement of material elements using the above velocity field at the time $t_0 + \Delta t$. i.e.,

$$\begin{equation}
\Delta x = u \Delta t \quad ; \quad \Delta y = v \Delta t
\end{equation}$$
- **Step (3):** We displace all the vortex sheet elements by the above displacement amount. Since we know from Equation (5) that vorticity moves with the fluid, we assume that the displaced vortex sheet element will have the same vorticity as in time $t_0$. The only difference is that the vortex elements are at different positions at time $t_0 +\Delta t$. We now go back to step (1) and proceed further in time.

It is important to note that the assumption of **Step (3)** amounts to solving the dynamical evolution equation for vorticity in 2D. But in three dimensions, there are vortex stretching and tilting terms in the governing equation (4). Vorticity does _not_ move with the fluid in three dimensions and we need to solve the equation (4) explicitly.

This also explains why taking each vorticity component separately will not resolve the issue. Even though vorticity components are scalars, they do not still move with the material element in 3D. 

## Introduction of a passive scalar
If we introduce a passive scalar in the flow, the governing equation of the scalar concentration ($c$) is, 

$$\begin{equation}
\frac{Dc}{Dt}=\kappa \nabla ^2 c
\end{equation}$$

where $\kappa$ is the molecular diffusivity of the scalar. 
Or in tensor notation,

$$ \begin{equation}
\frac{\partial c}{\partial t}+u_j\frac{\partial c}{\partial x_j}=\kappa \frac{\partial^2 c }{\partial x_j \partial x_j}
\end{equation}$$

Taking the gradient of the above equation,

$$\begin{equation}
\frac{\partial}{\partial x_i}\frac{\partial c}{\partial t}+\frac{\partial}{\partial x_i}\left( u_j\frac{\partial c}{\partial x_j} \right) = \kappa \frac{\partial}{\partial x_i} \left( \frac{\partial^2 c }{\partial x_j \partial x_j} \right)
\end{equation}$$

$$\begin{align*}
\frac{\partial}{\partial t} \frac{\partial c}{\partial x_i}+\frac{\partial u_j}{\partial x_i}\frac{\partial c}{\partial x_j} +u_j \frac{\partial^2 c }{\partial x_i \partial x_j}=&\kappa \frac{\partial^2 }{\partial x_j \partial x_j} \left(  \frac{\partial c}{\partial x_i} \right) \\
\frac{\partial}{\partial t} \frac{\partial c}{\partial x_i}+\frac{\partial u_j}{\partial x_i}\frac{\partial c}{\partial x_j} +u_j \frac{\partial }{\partial x_j}\left(\frac{\partial c}{\partial x_i}\right)=&\kappa \frac{\partial^2 }{\partial x_j \partial x_j} \left(  \frac{\partial c}{\partial x_i} \right)
\end{align*}$$

$$\begin{equation} \tag{6} % \label{eq:scalar_grad}
\frac{\partial}{\partial t} \left(\frac{\partial c}{\partial x_i}\right) +u_j \frac{\partial }{\partial x_j}\left(\frac{\partial c}{\partial x_i}\right)= -\frac{\partial u_j}{\partial x_i}\frac{\partial c}{\partial x_j} +\kappa \frac{\partial^2 }{\partial x_j \partial x_j} \left(  \frac{\partial c}{\partial x_i} \right) 
\end{equation}$$

Or in vector notation,

$$ \begin{equation} \tag{7} % \label{eq:scalar_grad1}
\frac{D}{Dt}\nabla c=-\nabla \vec{u} \cdot \nabla c + \kappa \nabla^2 (\nabla c)
\end{equation}$$

This is the governing equation for the gradient of the passive scalar. Note that there are similarities between this equation and the three-dimensional vorticity transport equation (4). The first term on the right hand side of this equation (7) is _like_ the vorticity tilting and stretching term, although with a negative sign. To see this more clearly, let us write the equation for the $x$-component scalar gradient in expanded form as,

$$ \begin{equation}
\frac{D}{Dt} \left( \frac{\partial c}{\partial x} \right) = -\left( \frac{\partial u}{\partial x} \frac{\partial c}{\partial x} +\frac{\partial v}{\partial x} \frac{\partial c}{\partial y} + \frac{\partial w}{\partial x} \frac{\partial c}{\partial z}\right) + \kappa \nabla ^2 \left( \frac{\partial c}{\partial x} \right)
\end{equation}$$

It is obvious that the first term on the right hand side of this equation is the strectching of scalar concentration in the $x$-direction and a decrease in $\frac{\partial c}{\partial x}$ results from this term. For the second term, we can imagine the presence of $\frac{\partial c}{\partial y}$ as the iso-contours of the scalar concentration lying on/parallel to the $x-z$ plane. The strain rate $\frac{\partial v}{\partial x}$ tilts these  surfaces like a rotation about the $z$- axis, thereby producing a component $-\frac{\partial c}{\partial x}$.

There is an interesting connection between scalar transport and vorticity transport. In order to see this, we first take the dot product of vorticity and Equation (6) to get,

$$ \begin{equation}
\omega_i \frac{\partial}{\partial t} \left(\frac{\partial c}{\partial x_i}\right) + \omega_i u_j \frac{\partial }{\partial x_j}\left(\frac{\partial c}{\partial x_i}\right)= -\omega_i \frac{\partial u_j}{\partial x_i}\frac{\partial c}{\partial x_j} +\kappa \omega_i \frac{\partial^2 }{\partial x_j \partial x_j} \left(  \frac{\partial c}{\partial x_i} \right) 
\end{equation}$$


$$\begin{equation} \tag{8} % \label{eq:temp1}
\vec{\omega} \cdot \frac{D}{Dt}\nabla c=- \vec{\omega} \cdot (\nabla \vec{u} \cdot \nabla c) + \kappa \vec{\omega} \cdot \nabla^2 (\nabla c)
\end{equation}$$

Now, we take the vorticity equation in the form,

$$ \begin{equation} \tag{9} % \label{eq:vort_3d_tensor}
\frac{\partial \omega_i}{\partial t}+u_j \frac{\partial \omega_i}{\partial x_j} = \omega_j \frac{\partial u_i}{\partial x_j} + \nu \frac{\partial ^2 \omega_i}{\partial x_j \partial x_j}
\end{equation}$$

Taking the dot product of this Equation (9) with the gradient of scalar concentration, we have,

$$ \begin{equation}
\frac{\partial c}{\partial x_i} \frac{\partial \omega_i}{\partial t}+\frac{\partial c}{\partial x_i} u_j \frac{\partial \omega_i}{\partial x_j} = \frac{\partial c}{\partial x_i} \omega_j \frac{\partial u_i}{\partial x_j} + \nu \frac{\partial c}{\partial x_i} \frac{\partial^2 \omega_i }{\partial x_j \partial x_j}
\end{equation}$$

$$ \begin{equation} \tag{10} % \label{eq:temp2}
\nabla c \cdot \frac{D \vec{\omega}}{Dt}= \left(\vec{\omega}.\nabla \vec{u}\right)\cdot \nabla c + \nu \nabla c \cdot \nabla ^2\vec{\omega}
\end{equation}$$

Adding Equation (8) \& Equation (10), after neglecting the effects of viscosity and molecular diffusion of the passive scalar, we have,

$$ \begin{equation} \tag{11} %\label{eq:scalar_vort}
\frac{D}{Dt}(\vec{\omega} \cdot \nabla c) = 0
\end{equation}$$

The conclusion from equation (11) is that the scalar quantity $\vec{\omega}\cdot \nabla c$ moves with the material element in 3D. If we start with a given distribution of the vorticity and scalar concentration, we can find out what is the scalar field $\vec{\omega}\cdot \nabla c$ at the initial time. As the flow evolves in time, each material element will preserve its initial value of $\vec{\omega}\cdot \nabla c$. But, this does not tell us how the material element moves. Without solving the Euler equation for velocity, or the equation for vorticity, I do not see any way in which we can get the evolution of the velocity field. Even if we know $\vec{\omega}\cdot \nabla c$ at a given time, it is not possible to reconstruct $\vec{\omega}$ from this scalar field.

A physical interpretation I have in mind for the Equation (11) is that the velocity gradient tensor affects the vorticity field in the exact opposite way in which it affects the scalar gradient field. Without the molecular diffusion, the scalar becomes a material property and gets advected with the flow. So we may assume it to have the characteristic of the velocity field $\vec{u}$. The gradient of this is affected in a particular way by the velocity gradient tensor. But vorticity being the anti-symmetric part of the velocity gradient tensor, it is affected by velocity gradient tensor in the opposite way (i.e. with a negative sign).
