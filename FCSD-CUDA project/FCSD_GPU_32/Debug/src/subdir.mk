################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/FCSD_detection_GPU.cpp \
../src/FCSD_ordering_CPU.cpp \
../src/communication_line.cpp 

CU_SRCS += \
../src/FCSD_decoding.cu \
../src/FCSD_ordering.cu \
../src/MATRIX_INVERSE_GJE.cu \
../src/chol.cu \
../src/chol_without_onchipmem.cu \
../src/row_column_transformation.cu 

CU_DEPS += \
./src/FCSD_decoding.d \
./src/FCSD_ordering.d \
./src/MATRIX_INVERSE_GJE.d \
./src/chol.d \
./src/chol_without_onchipmem.d \
./src/row_column_transformation.d 

OBJS += \
./src/FCSD_decoding.o \
./src/FCSD_detection_GPU.o \
./src/FCSD_ordering.o \
./src/FCSD_ordering_CPU.o \
./src/MATRIX_INVERSE_GJE.o \
./src/chol.o \
./src/chol_without_onchipmem.o \
./src/communication_line.o \
./src/row_column_transformation.o 

CPP_DEPS += \
./src/FCSD_detection_GPU.d \
./src/FCSD_ordering_CPU.d \
./src/communication_line.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-6.5/bin/nvcc -G -g -O0 -prec-div false -prec-sqrt false -gencode arch=compute_30,code=sm_30  -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-6.5/bin/nvcc -G -g -O0 -prec-div false -prec-sqrt false --compile --relocatable-device-code=false -gencode arch=compute_30,code=compute_30 -gencode arch=compute_30,code=sm_30  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-6.5/bin/nvcc -G -g -O0 -prec-div false -prec-sqrt false -gencode arch=compute_30,code=sm_30  -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-6.5/bin/nvcc -G -g -O0 -prec-div false -prec-sqrt false --compile  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


