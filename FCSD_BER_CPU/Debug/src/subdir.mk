################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../src/FCSD_detection_CPU.o \
../src/fullfactorial.o 

CPP_SRCS += \
../src/FCSD_detection_CPU.cpp \
../src/communication_line.cpp \
../src/fullfactorial.cpp 

OBJS += \
./src/FCSD_detection_CPU.o \
./src/communication_line.o \
./src/fullfactorial.o 

CPP_DEPS += \
./src/FCSD_detection_CPU.d \
./src/communication_line.d \
./src/fullfactorial.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


