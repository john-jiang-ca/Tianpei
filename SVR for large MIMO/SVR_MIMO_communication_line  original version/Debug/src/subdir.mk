################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/SVR_DETECTOR.cpp \
../src/SVR_MIMO_communication_line.cpp 

OBJS += \
./src/SVR_DETECTOR.o \
./src/SVR_MIMO_communication_line.o 

CPP_DEPS += \
./src/SVR_DETECTOR.d \
./src/SVR_MIMO_communication_line.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


