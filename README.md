# KernelEstimation

Supplementary material to the articles [''Online estimation of Hilbert-Schmidt operators and application to kernel reconstruction of neural fields''](https://www.google.com/url?q=https%3A%2F%2Fdoi-org.ezproxy.universite-paris-saclay.fr%2F10.1109%2FCDC51059.2022.9992414&sa=D) and
''Adaptive observer and control of spatiotemporal delayed neural fields''
by L. Brivadis, A. Chaillet and J. Auriol.

This page contains code for the numerical experiments of these papers in the branches main and No-delay, respectively.

## How to reproduce the experiments of the paper

From inside the main folder, run
	```
	matlab KernelEstimation.m
	```
	or
	```
	matlab KernelEstimation_delay.m
	```
to obtain the figures of the paper.

<br/><br/>

<p align="center">
	<img src="https://github.com/brivadis/KernelEstimation/blob/No-delays/fig2.gif" title="Evolution of the estimated kernel">
</p>
<figure>

<p align="center">
	<img src="https://github.com/brivadis/KernelEstimation/blob/main/obs.jpg" title="Convergence of the adaptive observer">
</p>
<figure>
