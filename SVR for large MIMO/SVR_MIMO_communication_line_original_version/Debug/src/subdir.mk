################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/SVR_DETECTOR.c \
../src/SVR_MIMO_communication_line.c 

OBJS += \
./src/SVR_DETECTOR.o \
./src/SVR_MIMO_communication_line.o 

C_DEPS += \
./src/SVR_DETECTOR.d \
./src/SVR_MIMO_communication_line.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


