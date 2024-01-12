# KernelEstimation

Supplementary material to the articles ''[Online estimation of Hilbert-Schmidt operators and application to kernel reconstruction of neural fields](https://ieeexplore-ieee-org.ezproxy.universite-paris-saclay.fr/document/9992414))'' and
''[Adaptive observer and control of spatiotemporal delayed neural fields](https://hal.science/hal-04106785)''
by L. Brivadis, A. Chaillet and J. Auriol.

This page contains code for the numerical experiments of these papers in the branches main and No-delay, respectively. Additional pertubations can be added in the branch With-delays for the second paper.

## How to reproduce the experiments of the paper

From inside the main folder, run
	```
	matlab KernelEstimation.m
	```
	or
	```
	matlab KernelEstimation_delay.m
	```
	or
	```
	matlab KernelEstimation_delay_pertubations.m
	```
to obtain the figures of the paper.

<br/><br/>

<p align="center">
	<img src="https://github.com/brivadis/KernelEstimation/blob/main/obs_gif.gif" title="Evolution of the estimated kernel">
</p>
<figure>

<p align="center">
	<img src="https://github.com/brivadis/KernelEstimation/blob/main/obs.jpg" title="Convergence of the adaptive observer">
</p>
<figure>
