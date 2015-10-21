################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/matrixMulCUBLAS.cpp 

OBJS += \
./src/matrixMulCUBLAS.o 

CPP_DEPS += \
./src/matrixMulCUBLAS.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-6.5/bin/nvcc -I"/usr/local/cuda-6.5/samples/0_Simple" -I"/usr/local/cuda-6.5/samples/common/inc" -I"/home/tchen44/~cd" -G -g -O0 -gencode arch=compute_30,code=sm_30  -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-6.5/bin/nvcc -I"/usr/local/cuda-6.5/samples/0_Simple" -I"/usr/local/cuda-6.5/samples/common/inc" -I"/home/tchen44/~cd" -G -g -O0 --compile  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


