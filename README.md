# Gyunity
This project's origin name is SmallSPPM since I have answered someones' question about sppm. But now I decide to extend it to a graphics library, so I rename this project to "Gyunity".

## Features
1. PathTracing
2. Stochastic progressive photon mapping (SPPM)
3. LightTracing
4. Bidirectional PathTracing (BDPT)
5. Multiple importance sampling (MIS)
6. Microfacet model with GGX distribution (GGX)
7. Ambient Occlusion (AO)
8. Low Discrepancy Sequence (Halton, Sobol...)
9. Kdtree, BVH
10. Fast vector and matrix calculation (SIMD-SSE)

## Some Results
<div align="center">
  <img width="400px" height="600px" src="https://github.com/Microsoft-Virtual-Ads/assets/blob/main/gallery/bunny2.png?raw=true"/>
  <img width="400px" height="600px" src="https://github.com/Microsoft-Virtual-Ads/assets/blob/main/gallery/bunny.png?raw=true"/>
</div>

<br/>
<div align="center">
  <img width="400px" height="600px" src="https://github.com/Microsoft-Virtual-Ads/assets/blob/main/gallery/coca2.png?raw=true"/>
  <img width="400px" height="600px" src="https://raw.githubusercontent.com/Microsoft-Virtual-Ads/assets/main/gallery/coca1.png"/>
</div>

<br/>
<div align="center">
  <img width="400px" height="600px" src="https://raw.githubusercontent.com/Microsoft-Virtual-Ads/assets/main/gallery/EiffelTower_GGX_RoughGlass_EnvMap.png"/>
  <img width="400px" height="600px" src="https://raw.githubusercontent.com/Microsoft-Virtual-Ads/assets/main/gallery/EiffelTower_Glass_EnvMap.png"/>
</div>

<br/>
<div align="center">
  <img width="400px" height="600px" src="https://raw.githubusercontent.com/Microsoft-Virtual-Ads/assets/main/gallery/Lucy_Bunny.png"/>
  <img width="400px" height="600px" src="https://raw.githubusercontent.com/Microsoft-Virtual-Ads/assets/main/gallery/water_cornell_box.png"/>
</div>

## Build
1. git clone --recursive https://github.com/xiaoyxue/Gyunity.git
2. In Developer Command Prompt for VS to run init.cmd
3. Open Gyunity.sln and build the solution