# Gyunity
This project's origin name is SmallSPPM since I have answered someones' question about sppm. But now I decide to extend it to a graphics library, so I rename this project to "Gyunity".

## Features
- PathTracing
- Stochastic progressive photon mapping (SPPM)
- LightTracing
- Bidirectional PathTracing (BDPT)
- Multiple importance sampling (MIS)
- Microfacet model with GGX distribution (GGX)
- Ambient Occlusion (AO)
- Low Discrepancy Sequence (Halton, Sobol...)
- Kdtree, BVH
- Fast vector and matrix calculation (SIMD-SSE)

## Some Results
<div align="center">
  <img width="400px" height="400px" src="https://github.com/xiaoyxue/GyunityAssets/blob/main/GyunityAssets/Gallery/bunny2.png?raw=true"/>
  <img width="400px" height="400px" src="https://github.com/xiaoyxue/GyunityAssets/blob/main/GyunityAssets/Gallery/bunny.png?raw=true"/>
</div>

<br/>
<div align="center">
  <img width="400px" height="400px" src="https://github.com/xiaoyxue/GyunityAssets/blob/main/GyunityAssets/Gallery/coca1.png?raw=true"/>
  <img width="400px" height="400px" src="https://github.com/xiaoyxue/GyunityAssets/blob/main/GyunityAssets/Gallery/coca2.png?raw=true"/>
</div>

<br/>
<div align="center">
  <img width="400px" height="400px" src="https://github.com/xiaoyxue/GyunityAssets/blob/main/GyunityAssets/Gallery/EiffelTower_GGX_RoughGlass_EnvMap.png?raw=true"/>
  <img width="400px" height="400px" src="https://github.com/xiaoyxue/GyunityAssets/blob/main/GyunityAssets/Gallery/EiffelTower_Glass_EnvMap.png?raw=true"/>
</div>

<br/>
<div align="center">
  <img width="400px" height="400px" src="https://github.com/xiaoyxue/GyunityAssets/blob/main/GyunityAssets/Gallery/Lucy_Bunny.png?raw=true"/>
  <img width="400px" height="400px" src="https://github.com/xiaoyxue/GyunityAssets/blob/main/GyunityAssets/Gallery/water_cornell_box.png?raw=true"/>
</div>

<br/>
<div align="center">
  <img width="400px" height="400px" src="https://github.com/xiaoyxue/GyunityAssets/blob/main/GyunityAssets/Gallery/bunny_envmap.gif?raw=true">
  <img width="400px" height="400px" src="https://github.com/xiaoyxue/GyunityAssets/blob/main/GyunityAssets/Gallery/sppm_cornell_box_area_light.gif?raw=true">
</div>

## Build
```bash
1. git clone --recursive https://github.com/xiaoyxue/Gyunity.git
2. In Developer Command Prompt for VS to run init.cmd
3. Open Gyunity.sln and build the solution
```
