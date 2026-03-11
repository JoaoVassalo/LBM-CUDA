nvcc -G -rdc=true main.cu functions/*.cu presets/stencil.cu -o main
 ./main