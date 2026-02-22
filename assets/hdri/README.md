# HDRI assets

`studio_small_09_1k.hdr` is [Studio Small 09](https://polyhaven.com/a/studio_small_09)
by Sergej Majboroda, distributed under CC0 by Poly Haven.

`studio_small_09_1k_soft.hdr` is derived from it by an equirectangular-aware
gaussian blur (sigma 2.5 px, horizontal sigma scaled by 1/sin(theta)). The
softboxes in the original carry a honeycomb diffuser grid, and a glass sphere
images the environment like a fisheye, so that grid shows up as a moire inside
every refractive object. Blurring removes it while preserving total irradiance
to within 0.3%. Use the original for scenes with no clear-glass hero object.
